#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "gmp.h"

// #define DEBUG

#ifdef DEBUG
#define DEBUG_CMD(cmd) cmd;
#else
#define DEBUG_CMD(CMD) ;
#endif

const mp_bitcnt_t precision = 200;
const unsigned int n = 1000;

void fact_calc(mpf_t* res, unsigned long num) {
	mpf_set_d(*res, 1);
	for (unsigned long i = 1; i <= num; i++) {
		mpf_mul_ui(*res, *res, i);
	}
}

int main(int argc, char** argv) {
	int root_proc_id = 0;
	int proc_id, num_procs;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	mpf_set_default_prec(precision);
	long block_size = n / num_procs;
	mpf_t fact;
	mpf_init(fact);
	fact_calc(&fact, block_size * proc_id + 1);
	mpf_t e;
	mpf_init(e);
	if (proc_id == root_proc_id) {
		mpf_set_d(e, 1);
	} else {
		mpf_set_d(e, 0);
	}
	unsigned long limit = proc_id != num_procs - 1 ? block_size * (proc_id + 1) : n;
	for (unsigned long i = block_size * proc_id + 1; i <= limit; i++) {
		mpf_t tmp;
		mpf_t one;
		mpf_init(tmp);
		mpf_init(one);
		mpf_set_d(one, 1);
		mpf_div(tmp, one, fact);
		mpf_add(e, e, tmp);
		mpf_mul_ui(fact, fact, i + 1);
	}
	if (proc_id != root_proc_id) {
		char* buff = calloc(2*precision, sizeof(char));
		MPI_Recv(buff, 2*precision, MPI_CHAR, proc_id - 1, 1, MPI_COMM_WORLD, &status);
		mpf_t recv;
		mpf_init(recv);
		gmp_sscanf(buff, "%Fe", recv);
		DEBUG_CMD(
		printf("proc id: %d, recv ", proc_id);
		mpf_out_str(stdout, 10, 0, recv);
		printf("\n");)
		mpf_add(e, e, recv);
		free(buff);
	}
	if (proc_id != num_procs - 1) {
		DEBUG_CMD(
		printf("proc id: %d, send ", proc_id);
		mpf_out_str(stdout, 10, 0, e);
		printf("\n");)
		char* buff = calloc(2*precision, sizeof(char));
		gmp_snprintf(buff, 2*precision, "%.Fe", e);
		MPI_Send(buff, 2*precision, MPI_CHAR, proc_id + 1, 1, MPI_COMM_WORLD);
		free(buff);
	} else {
		printf("e = ");
		mpf_out_str(stdout, 10, 0, e);
		printf("\n");
	}
	MPI_Finalize();
	return 0;
}
