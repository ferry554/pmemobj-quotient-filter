/*
 * test.cc
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

extern "C"
{
#include "pmem-qf.c"
}

#include <set>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <unistd.h>

#define POOL_SIZE	(1024 * 1024 * 1024) /* 1GB */

using namespace std;

/* I need a more powerful machine to increase these parameters... */
const uint32_t Q_MAX = 4;
const uint32_t R_MAX = 2;
const uint32_t ROUNDS_MAX = 20;

static void fail(TOID(struct quotient_filter) qf, const char *s)
{
	printf("qf(q=%u, r=%u): %s\n", D_RO(qf)->qf_qbits, D_RO(qf)->qf_rbits, s);
	abort();
}

static uint64_t rand64()
{
	return (((uint64_t)rand()) << 32) | ((uint64_t)rand());
}

static void qf_print(TOID(struct quotient_filter) qf)
{
	char buf[32];
	uint32_t pad = uint32_t(ceil(float(D_RO(qf)->qf_qbits) / logf(10.f))) + 1;

	for (uint32_t i = 0; i < pad; ++i)
	{
		printf(" ");
	}
	printf("| is_shifted | is_continuation | is_occupied | remainder"
		   " nel=%u\n",
		   D_RO(qf)->qf_entries);

	for (uint64_t idx = 0; idx < D_RO(qf)->qf_max_size; ++idx)
	{
		snprintf(buf, sizeof(buf), "%lu", idx);
		printf("%s", buf);

		int fillspace = pad - strlen(buf);
		for (int i = 0; i < fillspace; ++i)
		{
			printf(" ");
		}
		printf("| ");

		uint64_t elt = get_elem(qf, idx);
		printf("%d          | ", !!is_shifted(elt));
		printf("%d               | ", !!is_continuation(elt));
		printf("%d           | ", !!is_occupied(elt));
		printf("%lu\n", get_remainder(elt));
	}
}

/* Check QF structural invariants. */
static void qf_consistent(TOID(struct quotient_filter) qf)
{
	assert(D_RO(qf)->qf_qbits);
	assert(D_RO(qf)->qf_rbits);
	assert(D_RO(qf)->qf_qbits + D_RO(qf)->qf_rbits <= 64);
	assert(D_RO(qf)->qf_elem_bits == (D_RO(qf)->qf_rbits + 3));
	assert(D_RO(qf)->qf_table);

	uint64_t idx;
	uint64_t start;
	uint64_t size = D_RO(qf)->qf_max_size;
	assert(D_RO(qf)->qf_entries <= size);
	uint64_t last_run_elt;
	uint64_t visited = 0;

	if (D_RO(qf)->qf_entries == 0)
	{
		for (start = 0; start < size; ++start)
		{
			assert(get_elem(qf, start) == 0);
		}
		return;
	}

	for (start = 0; start < size; ++start)
	{
		if (is_cluster_start(get_elem(qf, start)))
		{
			break;
		}
	}

	assert(start < size);

	idx = start;
	do
	{
		uint64_t elt = get_elem(qf, idx);

		/* Make sure there are no dirty entries. */
		if (is_empty_element(elt))
		{
			assert(get_remainder(elt) == 0);
		}

		/* Check for invalid metadata bits. */
		if (is_continuation(elt))
		{
			assert(is_shifted(elt));

			/* Check that this is actually a continuation. */
			uint64_t prev = get_elem(qf, decr(qf, idx));
			assert(!is_empty_element(prev));
		}

		/* Check that remainders within runs are sorted. */
		if (!is_empty_element(elt))
		{
			uint64_t rem = get_remainder(elt);
			if (is_continuation(elt))
			{
				assert(rem > last_run_elt);
			}
			last_run_elt = rem;
			++visited;
		}

		idx = incr(qf, idx);
	} while (idx != start);

	assert(D_RO(qf)->qf_entries == visited);
}

/* Generate a random 64-bit hash. If @clrhigh, clear the high (64-p) bits. */
static uint64_t genhash(TOID(struct quotient_filter) qf, bool clrhigh,
						set<uint64_t> &keys)
{
	uint64_t hash;
	uint64_t mask = clrhigh ? LOW_MASK(D_RO(qf)->qf_qbits + D_RO(qf)->qf_rbits) : ~0ULL;
	uint64_t size = D_RO(qf)->qf_max_size;

	/* If the QF is overloaded, use a linear scan to find an unused hash. */
	if (keys.size() > (3 * (size / 4)))
	{
		uint64_t probe;
		uint64_t start = rand64() & D_RO(qf)->qf_index_mask;
		for (probe = incr(qf, start); probe != start; probe = incr(qf, probe))
		{
			if (is_empty_element(get_elem(qf, probe)))
			{
				uint64_t hi = clrhigh ? 0 : (rand64() & ~mask);
				hash = hi | (probe << D_RO(qf)->qf_rbits) | (rand64() & D_RO(qf)->qf_rmask);
				if (!keys.count(hash))
				{
					return hash;
				}
			}
		}
	}

	/* Find a random unused hash. */
	do
	{
		hash = rand64() & mask;
	} while (keys.count(hash));
	return hash;
}

/* Insert a random p-bit hash into the QF. */
static void ht_put(PMEMobjpool *pop, TOID(struct quotient_filter) qf, set<uint64_t> &keys)
{
	uint64_t hash = genhash(qf, true, keys);
	assert(qf_insert(pop,qf, hash));
	keys.insert(hash);
}

/* Remove a hash from the filter. */
static void ht_del(PMEMobjpool *pop, TOID(struct quotient_filter) qf, set<uint64_t> &keys)
{
	set<uint64_t>::iterator it;
	uint64_t idx = rand64() % keys.size();
	for (it = keys.begin(); it != keys.end() && idx; ++it, --idx)
		;
	uint64_t hash = *it;
	assert(qf_remove(pop,qf, hash));
	assert(!qf_may_contain(qf, hash));
	keys.erase(hash);
}

/* Check that a set of keys are in the QF. */
static void ht_check(TOID(struct quotient_filter) qf, set<uint64_t> &keys)
{
	qf_consistent(qf);
	set<uint64_t>::iterator it;
	for (it = keys.begin(); it != keys.end(); ++it)
	{
		uint64_t hash = *it;
		assert(qf_may_contain(qf, hash));
	}
}

static void qf_test_basic(PMEMobjpool *pop, TOID(struct quotient_filter) qf)
{
	// 此时传入的QF是刚初始化的
	// 基本插入和查询
	/* Basic get/set tests. */
	uint64_t idx;
	uint64_t size = D_RO(qf)->qf_max_size;
	for (idx = 0; idx < size; ++idx)
	{
		assert(get_elem(qf, idx) == 0);
		set_elem(qf, idx, idx & D_RO(qf)->qf_elem_mask);
	}
	for (idx = 0; idx < size; ++idx)
	{
		assert(get_elem(qf, idx) == (idx & D_RO(qf)->qf_elem_mask));
	}
	qf_clear(pop,qf);

	/* Random get/set tests. */
	// 随机插入和查询
	vector<uint64_t> elements(size, 0);
	for (idx = 0; idx < size; ++idx)
	{
		uint64_t slot = rand64() % size;
		uint64_t hash = rand64();
		set_elem(qf, slot, hash & D_RO(qf)->qf_elem_mask);
		elements[slot] = hash & D_RO(qf)->qf_elem_mask;
	}
	for (idx = 0; idx < elements.size(); ++idx)
	{
		assert(get_elem(qf, idx) == elements[idx]);
	}
	qf_clear(pop,qf);

	/* Check: forall x, insert(x) => may-contain(x). */
	// 如果插入了就必定存在，即不存在假阴性。
	set<uint64_t> keys;
	for (idx = 0; idx < size; ++idx)
	{
		uint64_t elt = genhash(qf, false, keys);
		assert(qf_insert(pop,qf, elt));
		keys.insert(elt);
	}
	ht_check(qf, keys);
	keys.clear();
	qf_clear(pop,qf);

	/* Check that the QF works like a hash set when all keys are p-bit values. */
	for (idx = 0; idx < ROUNDS_MAX; ++idx)
	{
		while (D_RO(qf)->qf_entries < size)
		{
			ht_put(pop,qf, keys);
		}

		while (D_RO(qf)->qf_entries > (size / 2))
		{
			ht_del(pop,qf, keys);
		}
		ht_check(qf, keys);

		struct qf_iterator qfi;
		qfi_start(qf, &qfi);
		while (!qfi_done(qf, &qfi))
		{
			uint64_t hash = qfi_next(qf, &qfi);
			assert(keys.count(hash));
		}
	}
}

/* Fill up the QF (at least partially). */
static void random_fill(PMEMobjpool *pop, TOID(struct quotient_filter) qf)
{
	set<uint64_t> keys;
	uint64_t elts = ((uint64_t)rand()) % D_RO(qf)->qf_max_size;
	while (elts)
	{
		ht_put(pop,qf, keys);
		--elts;
	}
	qf_consistent(qf);
}

/* Check if @lhs is a subset of @rhs. */
static void subsetof(TOID(struct quotient_filter) lhs, TOID(struct quotient_filter) rhs)
{
	struct qf_iterator qfi;
	qfi_start(lhs, &qfi);
	while (!qfi_done(lhs, &qfi))
	{
		uint64_t hash = qfi_next(lhs, &qfi);
		assert(qf_may_contain(rhs, hash));
	}
}

/* Check if @qf contains both @qf1 and @qf2. */
static void supersetof(TOID(struct quotient_filter) qf, TOID(struct quotient_filter) qf1,
					   TOID(struct quotient_filter) qf2)
{
	struct qf_iterator qfi;
	qfi_start(qf, &qfi);
	while (!qfi_done(qf, &qfi))
	{
		uint64_t hash = qfi_next(qf, &qfi);
		assert(qf_may_contain(qf1, hash) || qf_may_contain(qf2, hash));
	}
}

static void qf_bench(PMEMobjpool *pop,TOID(struct quotient_filter) qf1_bench)
{
	//struct quotient_filter qf;
	const uint32_t q_large = 10;
	const uint32_t q_small = 5;
	const uint32_t nlookups = 10000;
	struct timeval tv1, tv2;
	uint64_t sec;

	/* Test random inserts + lookups. */
	uint32_t ninserts = (3 * (1 << q_large) / 4); // 随机插入75%的元素
	printf("Testing %u random inserts and %u lookups\n", ninserts, nlookups);
	fflush(stdout);
	//在初始化之前，qf1_bench还全都为0
	qf_init(pop,qf1_bench, q_large, 1); // 初始化QF
	gettimeofday(&tv1, NULL);
	while (D_RO(qf1_bench)->qf_entries < ninserts)
	{
		assert(qf_insert(pop,qf1_bench, (uint64_t)rand())); // 还没满时每次插入都应该成功
		if (D_RO(qf1_bench)->qf_entries % 1000 == 0)
		{
			printf(".");
			fflush(stdout);
		}
	}
	for (uint32_t i = 0; i < nlookups; ++i)
	{
		qf_may_contain(qf1_bench, (uint64_t)rand()); // 随机查询
	}
	gettimeofday(&tv2, NULL);
	sec = tv2.tv_sec - tv1.tv_sec;
	printf(" done (%lu seconds).\n", sec);
	fflush(stdout);
	qf_destroy(pop,qf1_bench);

	/* Create a large cluster. Test random lookups. */
	qf_init(pop,qf1_bench, q_small, 1);
	printf("Testing %u contiguous inserts and %u lookups", 1 << q_small,
		   nlookups); // 更少的哈希桶
	fflush(stdout);
	gettimeofday(&tv1, NULL);
	for (uint64_t quot = 0; quot < (1 << (q_small - 1)); ++quot)
	{
		// 这一个cluster就是整个QF
		uint64_t hash = quot << 1;
		assert(qf_insert(pop,qf1_bench, hash));
		assert(qf_insert(pop,qf1_bench, hash | 1));
		if (quot % 2000 == 0)
		{
			printf(".");
			fflush(stdout);
		}
	}
	for (uint32_t i = 0; i < nlookups; ++i)
	{
		qf_may_contain(qf1_bench, (uint64_t)rand()); // 这再查询，就慢很多了
		if (i % 50000 == 0)
		{
			printf(".");
			fflush(stdout);
		}
	}
	gettimeofday(&tv2, NULL);
	sec = tv2.tv_sec - tv1.tv_sec;
	printf(" done (%lu seconds).\n", sec);
	fflush(stdout);
	qf_destroy(pop,qf1_bench);
}

static void qf_test(PMEMobjpool *pop,TOID(struct quotient_filter) qf1_test,
	TOID(struct quotient_filter) qf2_test,TOID(struct quotient_filter) qf21_test,TOID(struct quotient_filter) qf22_test)
{
	
	for (uint32_t q = 1; q <= Q_MAX; ++q)
	{
		printf("Starting rounds for qf_test::q=%u\n", q);
#pragma omp parallel for
		for (uint32_t r = 1; r <= R_MAX; ++r)
		{
			// 使用不同的q和r初始化QF，然后test
			//struct quotient_filter qf;
			if (!qf_init(pop,qf1_test, q, r))
			{
				fail(qf1_test, "init-1");
			}
			qf_test_basic(pop,qf1_test);
			qf_destroy(pop,qf1_test);
		}
	}

	for (uint32_t q1 = 1; q1 <= Q_MAX; ++q1)
	{
		for (uint32_t r1 = 1; r1 <= R_MAX; ++r1)
		{
			for (uint32_t q2 = 1; q2 <= Q_MAX; ++q2)
			{
				printf("Starting rounds for qf_merge::q1=%u,q2=%u\n", q1, q2);

#pragma omp parallel for
				for (uint32_t r2 = 1; r2 <= R_MAX; ++r2)
				{
					//struct quotient_filter qf;
					//struct quotient_filter qf1, qf2;
					if (!qf_init(pop,qf21_test, q1, r1) || !qf_init(pop,qf22_test, q2, r2))
					{
						fail(qf21_test, "init-2");
					}

					random_fill(pop,qf21_test);
					random_fill(pop,qf22_test);
					assert(qf_merge(pop,qf21_test, qf22_test, qf2_test));
					qf_consistent(qf2_test);
					subsetof(qf21_test, qf2_test);
					subsetof(qf22_test, qf2_test);
					supersetof(qf2_test, qf21_test, qf22_test);
					qf_destroy(pop,qf21_test);
					qf_destroy(pop,qf22_test);
					qf_destroy(pop,qf2_test);
				}
			}
		}
	}
}

int main(int argc, char *argv[])
{
	static PMEMobjpool *pop;

	srand(0);

	if(argc!=3||(strcmp(argv[2],"bench")&&strcmp(argv[2],"test")))
	{
		printf("usage: ./test <pmemfile> <mode>\nmode : bench or test\n");
		exit(1);
	}

	if (access(argv[1], F_OK))
	{
		// NVM文件还不存在，就创建
		if ((pop = pmemobj_create(argv[1],POBJ_LAYOUT_NAME(pmem_qf),POOL_SIZE, 0666)) == NULL)
		{
			fprintf(stderr, "%s", pmemobj_errormsg());
			exit(1);
		}
	}
	else
	{
		//NVM文件已经存在，直接打开
		if ((pop = pmemobj_open(argv[1],POBJ_LAYOUT_NAME(pmem_qf))) == NULL)
		{
			fprintf(stderr, "%s", pmemobj_errormsg());
			exit(1);
		}
	}
	TOID(struct my_root) root=POBJ_ROOT(pop,struct my_root);
	TX_BEGIN(pop) {
		TX_ADD(root);
		D_RW(root)->qf1_bench=TX_NEW(struct quotient_filter);
		D_RW(root)->qf1_test=TX_NEW(struct quotient_filter);
		D_RW(root)->qf2_test=TX_NEW(struct quotient_filter);
		D_RW(root)->qf21_test=TX_NEW(struct quotient_filter);
		D_RW(root)->qf22_test=TX_NEW(struct quotient_filter);
	}TX_END;

	if(!strcmp(argv[2],"bench"))
	{
		qf_bench(pop,D_RW(root)->qf1_bench);
	}
	else if(!strcmp(argv[2],"test"))
	{
		qf_test(pop,D_RW(root)->qf1_test,D_RW(root)->qf2_test,
			D_RW(root)->qf21_test,D_RW(root)->qf22_test);
	}
	else
	{
		printf("usage: ./test <pmemfile> <mode>\nmode : bench or test\n");
		exit(1);
	}
	puts("[PASSED] qf tests");
	pmemobj_close(pop);
	return 0;
}
