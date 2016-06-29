/*
 * juliet.h
 *
 *  Created on: Jun 29, 2016
 *      Author: saliya
 */

#ifndef JULIET_H_
#define JULIET_H_

#include <sched.h>

#endif /* JULIET_H_ */

void set_bit_mask(int rank, int thread_id, int tpp, cpu_set_t *mask) {
	/* Hard coded values for Juliet*/
	int cps = 12; // cores per socket
	int spn = 2; // sockets per node
	int htpc = 2; // hyper threads per core
	int cpn = cps * spn; // cores per node

	int ppn = cpn / tpp; // process per node
	int cpp = cpn / ppn; // cores per process

	// Assuming continuous ranking within a node
	int node_local_rank = rank % ppn;

	int j;
	for (j = 0; j < htpc; ++j) {
		CPU_SET(node_local_rank * cpp + thread_id + (cpn * j), mask);
	}
}
