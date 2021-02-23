#include <stdio.h>
#include <mpi.h>

// #define DEBUG

#ifdef DEBUG
#define DEBUG_CMD(cmd) cmd;
#else
#define DEBUG_CMD(CMD) ;
#endif

double fact_calc(long num) {
	double res = 1;
	for (long i = 1; i <= num; i++) {
		res *= i;
	}
	return res;
}

int main(int argc, char** argv) {
	long n = 1000;
	int root_proc_id = 0;
	int proc_id, num_procs;
	double recv = 0;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	long block_size = n / num_procs;
	double fact = fact_calc(block_size * proc_id + 1);
	double e = proc_id == root_proc_id ? 1 : 0; 
	long limit = proc_id != num_procs - 1 ? block_size * (proc_id + 1) : n;
	for (long i = block_size * proc_id + 1; i <= limit; i++) {
		e += 1 / fact;
		fact *= (i + 1);
	}
	if (proc_id != num_procs - 1) {
		MPI_Recv(&recv, 1, MPI_DOUBLE, proc_id + 1, 1, MPI_COMM_WORLD, &status);
		e += recv;
		DEBUG_CMD(printf("proc id: %d, recv %.20lf\n", proc_id, recv);)
	}
	if (proc_id != root_proc_id) {
		DEBUG_CMD(printf("proc id: %d, send %.20lf\n", proc_id, e);)
		MPI_Send(&e, 1, MPI_DOUBLE, proc_id - 1, 1, MPI_COMM_WORLD);
	} else {
		printf("e = %.20lf\n", e);
	}
	MPI_Finalize();
	return 0;
}
