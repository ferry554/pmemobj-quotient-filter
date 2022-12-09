/*
 * qf.h
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <libpmemobj.h>

//#define LAYOUT_NAME "pmem_qf"

POBJ_LAYOUT_BEGIN(pmem_qf);
POBJ_LAYOUT_ROOT(pmem_qf,struct my_root);
POBJ_LAYOUT_TOID(pmem_qf,struct quotient_filter);
POBJ_LAYOUT_TOID(pmem_qf,uint64_t);
POBJ_LAYOUT_END(pmem_qf);

struct my_root {
	TOID(struct quotient_filter) qf1_bench;
	TOID(struct quotient_filter) qf1_test;
	TOID(struct quotient_filter) qf2_test;
	TOID(struct quotient_filter) qf21_test;
	TOID(struct quotient_filter) qf22_test;
};


struct quotient_filter {
    //元数据
	uint8_t qf_qbits;//商长度
	uint8_t qf_rbits;//余数长度
	uint8_t qf_elem_bits;//整个elt长度，r+3
	uint64_t qf_index_mask;//取出一个elt
	uint64_t qf_rmask;//取出余数
	uint64_t qf_elem_mask;//取出商

    uint32_t qf_entries;//已有元素个数n？
	uint64_t qf_max_size;//最多元素个数m=2^q
    uint64_t* qf_table;

    //实现是以64bit为单位，但概念上是r+3 bit为单位
};

struct qf_iterator {
	uint64_t qfi_index;
	uint64_t qfi_quotient;
	uint64_t qfi_visited;
};

/*
 * Initializes a quotient filter with capacity 2^q.
 * Increasing r improves the filter's accuracy but uses more space.
 * 
 * Returns false if q == 0, r == 0, q+r > 64, or on ENOMEM.
 */
//需要写入，分配内存
bool qf_init(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint32_t q, uint32_t r);

/*
 * Inserts a hash into the QF.
 * Only the lowest q+r bits are actually inserted into the QF table.
 *
 * Returns false only if the QF is full.
 */
//需要写入
bool qf_insert(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint64_t hash);

/*
 * Returns true if the QF may contain the hash. Returns false otherwise.
 */
//查询是只读的
bool qf_may_contain(TOID(struct quotient_filter) qf, uint64_t hash);

/*
 * Removes a hash from the QF.
 *
 * Caution: If you plan on using this function, make sure that your hash
 * function emits no more than q+r bits. Consider the following scenario;
 *
 *	insert(qf, A:X)   # X is in the lowest q+r bits.
 *	insert(qf, B:X)   # This is a no-op, since X is already in the table.
 *	remove(qf, A:X)   # X is removed from the table.
 *
 * Now, may-contain(qf, B:X) == false, which is a ruinous false negative.
 *
 * Returns false if the hash uses more than q+r bits.
 */
//需要写入
bool qf_remove(PMEMobjpool *pop, TOID(struct quotient_filter) qf, uint64_t hash);


/*
 * Resets the QF table. This function does not deallocate any memory.
 */
//需要写入
void qf_clear(PMEMobjpool *pop, TOID(struct quotient_filter) qf);

/*
 * Deallocates the QF table.
 */
//需要写入，释放内存
void qf_destroy(PMEMobjpool *pop, TOID(struct quotient_filter) qf);

/*
 * Finds the size (in bytes) of a QF table.
 *
 * Caution: sizeof(struct quotient_filter) is not included.
 */
size_t qf_table_size(uint32_t q, uint32_t r);


/*
 * Initializes qfout and copies over all elements from qf1 and qf2.
 * Caution: qfout holds twice as many entries as either qf1 or qf2.
 *
 * Returns false on ENOMEM.
 */
bool qf_merge(PMEMobjpool *pop, TOID(struct quotient_filter) qf1, 
	TOID(struct quotient_filter) qf2, TOID(struct quotient_filter) qfout);



/*
 * Initialize an iterator for the QF.
 */
void qfi_start(TOID(struct quotient_filter) qf, struct qf_iterator *i);

/*
 * Returns true if there are no elements left to visit.
 */
bool qfi_done(TOID(struct quotient_filter) qf, struct qf_iterator *i);

/*
 * Returns the next (q+r)-bit fingerprint in the QF.
 *
 * Caution: Do not call this routine if qfi_done() == true.
 */
uint64_t qfi_next(TOID(struct quotient_filter) qf, struct qf_iterator *i);

