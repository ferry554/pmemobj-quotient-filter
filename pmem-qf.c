/*
 * qf.c
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

#include <stdlib.h>
#include <string.h>

#include "pmem-qf.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define LOW_MASK(n) ((1ULL << (n)) - 1ULL)
//ULL is used for Unsigned Long Long which is defined using 64 bits which can store large values.
//用于取出一个long long的低n位的掩码

TOID(uint64_t) qf_table_toid;
PMEMoid qf_table_oid;

//需要写入，是根API
bool qf_init(PMEMobjpool *pop,TOID(struct quotient_filter) qf, uint32_t q, uint32_t r)
{
	if (q == 0 || r == 0 || q + r > 64) {
		//如果商或余数长度为0，或者指纹总长度超过64，则是无效初始化
		return false;
	}

	
	bool ret;

    TX_BEGIN(pop){
		//要写入整个QF结构体的内容，故ADD整个qf，这样使用的元数据最少
		//新分配的内存不需要添加
        TX_ADD(qf);

        D_RW(qf)->qf_qbits = q;//商的长度
        D_RW(qf)->qf_rbits = r;//余数长度
        D_RW(qf)->qf_elem_bits = D_RO(qf)->qf_rbits + 3;//一个slot中存储长度，r+3
        D_RW(qf)->qf_index_mask = LOW_MASK(q);
        D_RW(qf)->qf_rmask = LOW_MASK(r);
        D_RW(qf)->qf_elem_mask = LOW_MASK(D_RO(qf)->qf_elem_bits);

        D_RW(qf)->qf_entries = 0; 
        D_RW(qf)->qf_max_size = 1 << q;//当前已有0个元素，最多有2^q个元素

		//如果分配失败，事务会自动abort
		qf_table_toid=TX_ZALLOC(uint64_t,qf_table_size(q, r));
		qf_table_oid=qf_table_toid.oid;
        D_RW(qf)->qf_table=(uint64_t*)pmemobj_direct(qf_table_oid);

	}TX_ONABORT{
		ret=false;
	}TX_ONCOMMIT{
		ret=true;
	}TX_END;

    return ret;
}

/* Return QF[idx] in the lower bits. */
//根据在QF中的id来索引到桶，不需要写入
static uint64_t get_elem(TOID(struct quotient_filter) qf, uint64_t idx)
{
	uint64_t elt = 0;
	//bit position
	size_t bitpos = D_RO(qf)->qf_elem_bits * idx;

	//tab position，在第几个uint64里
	size_t tabpos = bitpos / 64;

	//在所在的unint64里的偏移是多少
	size_t slotpos = bitpos % 64;
	int spillbits = (slotpos + D_RO(qf)->qf_elem_bits) - 64;

	//结果是根据tab找到相应uint64，将其读出
	elt = D_RO(qf)->qf_table [tabpos] >> slotpos & D_RO(qf)->qf_elem_mask;
	if (spillbits > 0) {
		//大于0说明是跨tab存储
		++tabpos;
		uint64_t x = D_RO(qf)->qf_table [tabpos] & LOW_MASK(spillbits);
		elt |= x << (D_RO(qf)->qf_elem_bits - spillbits);
	}
	return elt;
}

/* Store the lower bits of elt into QF[idx]. */
//根据在QF中的id来索引到桶，并设置r+3 bit的数据，需要写入，不是根API
static void set_elem(TOID(struct quotient_filter) qf, uint64_t idx, uint64_t elt)
{
    size_t bitpos = D_RO(qf)->qf_elem_bits * idx;
    size_t tabpos = bitpos / 64;
    size_t slotpos = bitpos % 64;
    int spillbits = (slotpos + D_RO(qf)->qf_elem_bits) - 64;
    elt &= D_RO(qf)->qf_elem_mask;
    D_RO(qf)->qf_table [tabpos] &= ~(D_RO(qf)->qf_elem_mask << slotpos);
    D_RO(qf)->qf_table [tabpos] |= elt << slotpos;
    if (spillbits > 0) {
        ++tabpos;
        D_RO(qf)->qf_table [tabpos] &= ~LOW_MASK(spillbits);
        D_RO(qf)->qf_table [tabpos] |= elt >> (D_RO(qf)->qf_elem_bits - spillbits);
    }

}

//针对特定QF的游标增减，实现了wrap
//不需写入
static inline uint64_t incr(TOID(struct quotient_filter) qf, uint64_t idx)
{
	return (idx + 1) & D_RO(qf)->qf_index_mask;
}

static inline uint64_t decr(TOID(struct quotient_filter) qf, uint64_t idx)
{
	return (idx - 1) & D_RO(qf)->qf_index_mask;
}

//以下为处理elt中的三个标志位的函数，对一个bit有is,set,clr操作
static inline int is_occupied(uint64_t elt)
{
	return elt & 1;
}

static inline uint64_t set_occupied(uint64_t elt)
{
	return elt | 1;
}

static inline uint64_t clr_occupied(uint64_t elt)
{
	return elt & ~1;
}

static inline int is_continuation(uint64_t elt)
{
	return elt & 2;
}

static inline uint64_t set_continuation(uint64_t elt)
{
	return elt | 2;
}

static inline uint64_t clr_continuation(uint64_t elt)
{
	return elt & ~2;
}

static inline int is_shifted(uint64_t elt)
{
	return elt & 4;
}

static inline uint64_t set_shifted(uint64_t elt)
{
	return elt | 4;
}

static inline uint64_t clr_shifted(uint64_t elt)
{
	return elt & ~4;
}

//从一个elt中提取余数
static inline uint64_t get_remainder(uint64_t elt)
{
	return elt >> 3;
}

//三个标志位均为0时说明是空槽
static inline bool is_empty_element(uint64_t elt)
{
	return (elt & 7) == 0;
}

//100说明是cluster的开始
static inline bool is_cluster_start(uint64_t elt)
{
	return is_occupied(elt) && !is_continuation(elt) && !is_shifted(elt);
}

//100和x01说明是run的开始
static inline bool is_run_start(uint64_t elt)
{
	return !is_continuation(elt) && (is_occupied(elt) || is_shifted(elt));
}

//根据hash生成对这个QF的商和余数
//不需写入
static inline uint64_t hash_to_quotient(TOID(struct quotient_filter) qf,
		uint64_t hash)
{
	//先把余数部分去了，再取商的掩码
	return (hash >> D_RO(qf)->qf_rbits) & D_RO(qf)->qf_index_mask;
}

static inline uint64_t hash_to_remainder(TOID(struct quotient_filter) qf,
		uint64_t hash)
{
	//余数在最低位，直接取掩码
	return hash & D_RO(qf)->qf_rmask;
}

//定位一个商所属的run的实际位置
/* Find the start index of the run for fq (given that the run exists). */
//不需写入
static uint64_t find_run_index(TOID(struct quotient_filter) qf, uint64_t fq)
{
	/* Find the start of the cluster. */
	//从本位开始向左扫描到cluster的开始
	uint64_t b = fq;
	while (is_shifted(get_elem(qf, b))) {
		b = decr(qf, b);
	}

	/* Find the start of the run for fq. */
	//有几个occupied，就说明有几种不同的商q的指纹被插入，那就是有多少个run
	uint64_t s = b;
	while (b != fq) {
		do {
			s = incr(qf, s);
		} while (is_continuation(get_elem(qf, s)));

		do {
			b = incr(qf, b);
		} while (!is_occupied(get_elem(qf, b)));
	}//向右扫描到商所属run的开始
	return s;
}

/* Insert elt into QF[s], shifting over elements as necessary. */
//需要写入，不是根API
static void insert_into(TOID(struct quotient_filter) qf, uint64_t s, uint64_t elt)
{
	uint64_t prev;
	uint64_t curr = elt;
	bool empty;

	//在s处插入elt，然后把原有的数据挤到下一个桶中，直到挤到一个空桶里
	//使用O(1)的空间来存储temp数据，
	do {
		prev = get_elem(qf, s);//prev是当前s处的元素
		empty = is_empty_element(prev);//prev是否是000的空位
		if (!empty) {
			/* Fix up `is_shifted' and `is_occupied'. */
			prev = set_shifted(prev);//要移位了
			if (is_occupied(prev)) {//isO和桶对应，而不和桶中存储的余数对应
				curr = set_occupied(curr);
				prev = clr_occupied(prev);
			}
		}
		set_elem(qf, s, curr);
		curr = prev;
		s = incr(qf, s);
	} while (!empty);
}

//需要写入，是根API
bool qf_insert(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint64_t hash)
{
	if (D_RO(qf)->qf_entries >= D_RO(qf)->qf_max_size) {
		//QF已满
		return false;
	}

	//根据hash得到商和余数，根据商得到本位elt，并将余数和000结合准备插入
	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);
	uint64_t entry = (fr << 3) & ~7;

    bool ret;
	uint64_t start;
	uint64_t s;

    TX_BEGIN(pop) {
        //要修改qf_table中的内容
		/*pmemobj_tx_add_range_direct(D_RO(qf)->qf_table,
            qf_table_size(D_RO(qf)->qf_qbits, D_RO(qf)->qf_rbits));
        //要修改qf的qf_entries字段
        TX_ADD_FIELD(qf,qf_entries);*/
		TX_ADD(qf);
        

        /* Special-case filling canonical slots to simplify insert_into(). */
        if (is_empty_element(T_fq)) {
            //如果本位是000，说明本位是空位，设置为100直接插入本位
            set_elem(qf, fq, set_occupied(entry));
            ++D_RW(qf)->qf_entries;
            goto end;
        }

        if (!is_occupied(T_fq)) {
            //如果本位是001或011，说明该run还不存在但是本位已经被其他run侵占。
            //先设置本位的isO，使得该商有run。
            set_elem(qf, fq, set_occupied(T_fq));
        }
        //之后遵从和查询类似的方式，先找到该商的run
        start = find_run_index(qf, fq);
        s = start;

        if (is_occupied(T_fq)) {
            //本位的isO为1，说明在执行本次插入操作之前，该商的run已经存在了。
            /* Move the cursor to the insert position in the fq run. */
            //从该run的起始处开始，查询每个桶中存储的余数
            //应该是升序存储。如果等于，说明发生了硬冲突，可以中止插入，返回
            //如果大于，说明找到了应插入的位置
            do {
                uint64_t rem = get_remainder(get_elem(qf, s));
                if (rem == fr) {
                    goto end;
                } else if (rem > fr) {
                    break;
                }
                s = incr(qf, s);
            } while (is_continuation(get_elem(qf, s)));

            //s此时是要插入的位置
            if (s == start) {
                //应该插到该run的起始处。之前的起始应该后移。
                //找到起始的elt，将其isC设为1
                /* The old start-of-run becomes a continuation. */
                uint64_t old_head = get_elem(qf, start);
                set_elem(qf, start, set_continuation(old_head));
            } else {
                /* The new element becomes a continuation. */
                //不需要插入起始处，设置要插入的isC为1即可
                entry = set_continuation(entry);
            }
        }

        /* Set the shifted bit if we can't use the canonical slot. */
        //设置isS,插入
        if (s != fq) {
            entry = set_shifted(entry);
        }

        insert_into(qf, s, entry);
        ++D_RW(qf)->qf_entries;
        //pmemobj_tx_process();
		end:
		uint8_t tee=1;

    }TX_END;

    return true;
}

//不需写入
bool qf_may_contain(TOID(struct quotient_filter) qf, uint64_t hash)
{
	//得到hash的商和余数
	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);

	//根据商得到本位中存储的数据
	uint64_t T_fq = get_elem(qf, fq);

	/* If this quotient has no run, give up. */
	if (!is_occupied(T_fq)) {
		//如果isO为0，说明该商的run不存在，元素也一定不存在
		return false;
	}

	/* Scan the sorted run for the target remainder. */
	//否则run存在，定位这个商的run的起始位置
	uint64_t s = find_run_index(qf, fq);
	do {
		//根据位置，先得到elt，再得到余数
		uint64_t rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			return true;//存在这个余数，可能存在
		} else if (rem > fr) {
			return false;//按序查找已经直接超过了，说明一定不存在
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));//直到该run结束也未找到
	return false;//一定不存在
}

/* Remove the entry in QF[s] and slide the rest of the cluster forward. */
//需要写入，不是根API
static void delete_entry(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint64_t s, uint64_t quot)
{
	uint64_t next;
	uint64_t curr = get_elem(qf, s);
	uint64_t sp = incr(qf, s);
	uint64_t orig = s;

	/*
	 * FIXME(vsk): This loop looks ugly. Rewrite.
	 */
	while (true) {
		next = get_elem(qf, sp);
		bool curr_occupied = is_occupied(curr);

		if (is_empty_element(next) || is_cluster_start(next) || sp == orig) {
			set_elem(qf, s, 0);
			return;
		} else {
			/* Fix entries which slide into canonical slots. */
			uint64_t updated_next = next;
			if (is_run_start(next)) {
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));

				if (curr_occupied && quot == s) {
					updated_next = clr_shifted(next);
				}
			}

			set_elem(qf, s, curr_occupied ?
					set_occupied(updated_next) :
					clr_occupied(updated_next));
			s = sp;
			sp = incr(qf, sp);
			curr = next;
		}
	}
}

//需要写入，是根API
bool qf_remove(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint64_t hash)
{
	uint64_t highbits = hash >> (D_RO(qf)->qf_qbits + D_RO(qf)->qf_rbits);
	if (highbits) {
		return false;
	}

	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);

	if (!is_occupied(T_fq) || !D_RO(qf)->qf_entries) {
		return true;
	}

	uint64_t start = find_run_index(qf, fq);
	uint64_t s = start;
	uint64_t rem;

	/* Find the offending table index (or give up). */
	do {
		rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			break;
		} else if (rem > fr) {
			return true;
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));
	if (rem != fr) {
		return true;
	}

	uint64_t kill = (s == fq) ? T_fq : get_elem(qf, s);
	bool replace_run_start = is_run_start(kill);

    bool ret;

    TX_BEGIN(pop) {
        //要修改qf_table中的内容
		/*pmemobj_tx_add_range_direct(D_RO(qf)->qf_table,
            qf_table_size(D_RO(qf)->qf_qbits, D_RO(qf)->qf_rbits));
        //要修改qf的qf_entries字段
        TX_ADD_FIELD(qf,qf_entries);*/
		TX_ADD(qf);
        

        /* If we're deleting the last entry in a run, clear `is_occupied'. */
        if (is_run_start(kill)) {
            uint64_t next = get_elem(qf, incr(qf, s));
            if (!is_continuation(next)) {
                T_fq = clr_occupied(T_fq);
                set_elem(qf, fq, T_fq);
            }
        }

        delete_entry(pop,qf, s, fq);

        if (replace_run_start) {
            uint64_t next = get_elem(qf, s);
            uint64_t updated_next = next;
            if (is_continuation(next)) {
                /* The new start-of-run is no longer a continuation. */
                updated_next = clr_continuation(next);
            }
            if (s == fq && is_run_start(updated_next)) {
                /* The new start-of-run is in the canonical slot. */
                updated_next = clr_shifted(updated_next);
            }
            if (updated_next != next) {
                set_elem(qf, s, updated_next);
            }
        }

        --D_RW(qf)->qf_entries;

    }TX_END;

    return true;
}

//清空QF的存储空间，是根API
void qf_clear(PMEMobjpool *pop, TOID(struct quotient_filter) qf)
{
    TX_BEGIN(pop) {
        //要修改qf的qf_entries字段
        //TX_ADD_FIELD(qf,qf_entries);
		TX_ADD(qf);
        D_RW(qf)->qf_entries = 0;

        //要修改qf_table中的内容
        //pmemobj_tx_add_range_direct(D_RO(qf)->qf_table,
            //qf_table_size(D_RO(qf)->qf_qbits, D_RO(qf)->qf_rbits));
        TX_MEMSET(D_RO(qf)->qf_table, 0, 
            qf_table_size(D_RO(qf)->qf_qbits, D_RO(qf)->qf_rbits));
        
    } TX_END;

}

//销毁QF，是根API
void qf_destroy(PMEMobjpool *pop, TOID(struct quotient_filter) qf)
{
    TX_BEGIN(pop) {
        //分配和释放内存都要添加整个qf
        TX_ADD(qf);
        TX_FREE(qf_table_toid);
        D_RW(qf)->qf_table=NULL;
    } TX_END; 
}

//QF的存储空间，即2^q*(r+3)，返回向上取整的字节大小
size_t qf_table_size(uint32_t q, uint32_t r)
{
	size_t bits = (1 << q) * (r + 3);
	size_t bytes = bits / 8;
	return (bits % 8) ? (bytes + 1) : bytes;
}

//根API，需要写入
bool qf_merge(PMEMobjpool *pop, TOID(struct quotient_filter) qf1, 
	TOID(struct quotient_filter) qf2, TOID(struct quotient_filter) qfout)
{
	uint32_t q = 1 + MAX(D_RO(qf1)->qf_qbits, D_RO(qf2)->qf_qbits);
	uint32_t r = MAX(D_RO(qf1)->qf_rbits, D_RO(qf2)->qf_rbits);

	bool ret;
	struct qf_iterator qfi;

	TX_BEGIN(pop) {
		if (!qf_init(pop,qfout, q, r)) {
			pmemobj_tx_abort(-1);
		}
		qfi_start(qf1, &qfi);
		while (!qfi_done(qf1, &qfi)) {
			qf_insert(pop,qfout, qfi_next(qf1, &qfi));
		}
		qfi_start(qf2, &qfi);
		while (!qfi_done(qf2, &qfi)) {
			qf_insert(pop,qfout, qfi_next(qf2, &qfi));
		}
	}TX_ONABORT{
		ret=false;
	}TX_ONCOMMIT{
		ret=true;
	}TX_END;

	return ret;
}


//迭代部分，只在QF的merge中用到
void qfi_start(TOID(struct quotient_filter) qf, struct qf_iterator *i)
{
	/* Mark the iterator as done. */
	i->qfi_visited = D_RO(qf)->qf_entries;

	if (D_RO(qf)->qf_entries == 0) {
		return;
	}

	/* Find the start of a cluster. */
	uint64_t start;
	for (start = 0; start < D_RO(qf)->qf_max_size; ++start) {
		if (is_cluster_start(get_elem(qf, start))) {
			break;
		}
	}

	i->qfi_visited = 0;
	i->qfi_index = start;
}

bool qfi_done(TOID(struct quotient_filter) qf, struct qf_iterator *i)
{
	return D_RO(qf)->qf_entries == i->qfi_visited;
}

uint64_t qfi_next(TOID(struct quotient_filter) qf, struct qf_iterator *i)
{
	while (!qfi_done(qf, i)) {
		uint64_t elt = get_elem(qf, i->qfi_index);

		/* Keep track of the current run. */
		if (is_cluster_start(elt)) {
			i->qfi_quotient = i->qfi_index;
		} else {
			if (is_run_start(elt)) {
				uint64_t quot = i->qfi_quotient;
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));
				i->qfi_quotient = quot;
			}
		}

		i->qfi_index = incr(qf, i->qfi_index);

		if (!is_empty_element(elt)) {
			uint64_t quot = i->qfi_quotient;
			uint64_t rem = get_remainder(elt);
			uint64_t hash = (quot << D_RO(qf)->qf_rbits) | rem;
			++i->qfi_visited;
			return hash;
		}
	}

	abort();
}

