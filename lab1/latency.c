#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int proc_id, proc_num;
    int n = 100000;
    int count = 1024;
    int buf_size = count*sizeof(int);
    int* buf = malloc(buf_size);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    if (proc_num == 2) {
        if (proc_id == 0) {
            for (int i = 0; i < n; i++) {
                MPI_Recv(buf, count, MPI_INT, 1, 0, MPI_COMM_WORLD, NULL);
            }
        }
        else {
            clock_t t;
            t = clock();
            for (int i = 0; i < n; i++) {
                MPI_Send(buf, count, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            t = clock() - t;
            double time = (double)t/(CLOCKS_PER_SEC/1000);
            printf("Latency %lf ms\n", time/n);
            FILE* fout = fopen("latency.txt", "w");
            fprintf(fout, "%d;%lf", count, time/n);
            fclose(fout);
        }
    }
    free(buf);
    MPI_Finalize();
}
