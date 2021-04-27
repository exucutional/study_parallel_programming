#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#define DEBUG

#ifdef DEBUG
#define DEBUG_CMD(cmd) cmd;
#else
#define DEBUG_CMD(CMD) ;
#endif

const double ranget = 1;
const double rangex = 1;
const double initt = 0;
const double initx = 0;
const double h = 0.2; // 2t >= h
const double t = 0.1;
const double a = 2;

typedef struct {
    double value;
} Cell;

typedef struct {
    Cell** cell;
    long dimt;
} Grid;

Grid* grid_init(long dimt, long dimx, int block_size, bool isRootGrid) {
    Cell** cell = calloc(dimt, sizeof(Cell*));
    if (!isRootGrid) {
        Cell* cell_rows = calloc(dimt*block_size, sizeof(Cell));
        for (int i = 0; i < dimt; i++) {
            cell[i] = &(cell_rows[i*block_size]);
        }
    }
    else {
        Cell* cell_rows = calloc(dimt*(dimx+dimt-1), sizeof(Cell));
        for (int i = 0; i < dimt; i++) {
            cell[i] = &(cell_rows[i*(dimx+dimt-1)]);
        }
    }
    Grid* grid = calloc(1, sizeof(Grid));
    grid->cell = cell;
    grid->dimt = dimt;
    return grid;
}

void grid_free(Grid* grid) {
    free(grid->cell[0]);
    free(grid->cell);
    free(grid);
}

double f_init_x(double x) {
    return cos(3.1415*x);
}

double f_init_t(double t) {
    return exp(-t);
}

double f(double t, double x) {
    return x + t;
}

void calc_cell(Grid* grid, long ti, long xi, double t, double h, long block_size, long offset, int proc_id, int proc_num, long tag_recv, long tag_send)
{
    Cell** cell = grid->cell;
    double u;

    if (ti == 0) {
        u = f_init_x((xi+offset)*h);
    }
    else if (xi == 0 && offset == 0) {
        u = f_init_t(ti*t);
    }
    else {
        double ul, ur;
        if (xi-1 < 0) {
            MPI_Recv(&ul, 1, MPI_DOUBLE, proc_id-1, tag_recv-1, MPI_COMM_WORLD, NULL);
            //DEBUG_CMD(printf("[rank: %d] [tag %ld] Recv %.2lf\n", proc_id, tag_recv-1, ul);)
        }
        else {
            ul = cell[ti-1][xi-1].value;
        }
        if (xi+1 > block_size ) {
            MPI_Recv(&ur, 1, MPI_DOUBLE, proc_id+1, tag_recv+block_size, MPI_COMM_WORLD, NULL);
            //DEBUG_CMD(printf("[rank: %d] [tag %ld] Recv %.2lf\n", proc_id, tag_recv+block_size, ur);)
        }
        else {
            ur = cell[ti-1][xi+1].value;
        }

        double tmp1 = 0.5*(ur+ul);
        double tmp2 = -a*0.5*t/h*(ur-ul);
        u = t*f(ti*t-t, (xi+offset)*h)+tmp1+tmp2;
    }

    cell[ti][xi].value = u;

    MPI_Request request;

    if (xi-1 < 0 && offset != 0) {
        //DEBUG_CMD(printf("[rank: %d] [tag %ld] Sent %.2lf\n", proc_id, tag_send, u);)
        MPI_Isend(&cell[ti][xi].value, 1, MPI_DOUBLE, proc_id-1, tag_send, MPI_COMM_WORLD, &request);
    }
    if (xi+1 > block_size-1 && proc_id+1 != proc_num) {
        //DEBUG_CMD(printf("[rank: %d] [tag %ld] Sent %.2lf\n", proc_id, tag_send, u);)
        MPI_Isend(&cell[ti][xi].value, 1, MPI_DOUBLE, proc_id+1, tag_send, MPI_COMM_WORLD, &request);
    }
}

void calc_grid(Grid* grid, long dimt, long dimx, long block_size, long offset, int proc_id, int proc_num) {
    for (int ti = 0; ti < dimt; ti++) {
        long t_offset = block_size+offset-dimt;
        for (int xi = 0; xi < block_size; xi++) {
            if (dimt+dimx-1-(xi+offset)-ti > 0) {
                long tag_recv = (ti-1)*dimt+offset+xi;
                long tag_send = tag_recv+dimt;
                //DEBUG_CMD(printf("[rank: %d] [ti: %d] [xi: %d] Start calculation\n", proc_id, ti, xi));
                calc_cell(grid, ti, xi, t, h, block_size, offset, proc_id, proc_num, tag_recv, tag_send);
            }
        }
    }
}

void root_grid_copy(Grid* root_grid, Grid* grid, long dimt, long dimx, long block_size, long offset) {
    for (int ti = 0; ti < dimt; ti++) {
        for (int xi = 0; xi < block_size; xi++) {
            root_grid->cell[ti][xi+offset].value = grid->cell[ti][xi].value;
        }
    }
}

long get_block_size(long dimt, long dimx, int proc_id, int proc_num, long* offset) {
    long block_size = (dimt+dimx-1)/proc_num;
    if (block_size*proc_id > dimx) {
        block_size += block_size*proc_id / dimx;
    }
    block_size--;
    if (offset) {
        *offset = 0;
        int rank = proc_id;
        while (rank > 0) {
            rank--;
            *offset += get_block_size(dimt, dimx, rank, proc_num, NULL);
        }
    }
    if (proc_num == proc_id+1) {
        block_size = dimx+dimt-1-*offset;
    }
    return block_size;
}

void file_dump(Grid* grid, long dimt, long dimx)
{
    FILE* fout = fopen("grid_par.csv", "w");
    for (int ti = 0; ti < dimt; ti++) {
        fprintf(fout, "%lf;", ti*t);
    }
    fprintf(fout, "\n");
    for (int xi = 0; xi < dimx; xi++) {
        fprintf(fout, "%lf;", xi*h);
    }
    fprintf(fout, "\n");
    for (int ti = 0; ti < dimt; ti++) {
        for (int xi = 0; xi < dimx; xi++) {
            fprintf(fout, "%lf;", grid->cell[ti][xi].value);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}

void stdout_dump(Grid* grid_root, long dimt, long dimx, long block_size, long offset, int proc_num)
{
    printf("Grid\n");
    for (int ti = dimt-1; ti > -1; ti--) {
        for (int xi = 0; xi < dimx+dimt-1-ti; xi++) {
            printf("|%5.2lf ", grid_root->cell[ti][xi].value);
        }
        printf("|\n");
    }
    int rank = 0;
    for (int tx = 0; tx < dimt+dimx-1; tx++) {
        if (offset == 0 && tx == 0) {
            printf("|rank %d", rank);
        }
        else if (tx > block_size+offset-1) {
            rank++;
            printf("|rank %d", rank);
            block_size = get_block_size(dimt, dimx, rank, proc_num, &offset);
        }
        else {
            printf("_______");
        }
    }
    printf("|\n");
}

int main(int argc, char** argv) {
    int root_proc_id = 0;
    int proc_id, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    long dimt = (long)(ranget-initt)/t;
    long dimx = (long)(rangex-initx)/h;
    long offset = 0;
    long block_size = get_block_size(dimt, dimx, proc_id, proc_num, &offset);
    //DEBUG_CMD(printf("[rank: %d] [block_size: %ld] [offset: %ld]\n", proc_id, block_size, offset);)
    Grid* grid = grid_init(dimt, dimx, block_size, false);
    Cell** cell = grid->cell;
    calc_grid(grid, dimt, dimx, block_size, offset, proc_id, proc_num);
    if (proc_id != root_proc_id) {
        MPI_Send(&(cell[0][0]), dimt*block_size, MPI_DOUBLE, root_proc_id, 0, MPI_COMM_WORLD);
        //DEBUG_CMD(printf("[rank: %d] Sent grid with size %ld\n", proc_id, dimt*block_size);)
        grid_free(grid); 
    }
    if (proc_id == root_proc_id) {
        Grid* grid_root = grid_init(dimt, dimx, block_size, true);
        // Копирование для того чтобы собрать результат и наглядно отобразить в консоль
        // Однако сам результат можно вывести например в файл без соединения всех сеток в одну
        // Поэтому при подсчете времени не будут учитываться эта функция копирования с пересылкой, а также дамп в консоль и дамп в файл.
        root_grid_copy(grid_root, grid, dimt, dimx, block_size, offset);
        for (int proc_iter = root_proc_id+1; proc_iter < proc_num; proc_iter++) {
            long offset_copy = 0;
            long block_size_copy = get_block_size(dimt, dimx, proc_iter, proc_num, &offset_copy);
            Grid* grid_tmp = grid_init(dimt, dimx, block_size_copy, false);
            MPI_Recv(&(grid_tmp->cell[0][0]), dimt*block_size_copy, MPI_DOUBLE, proc_iter, 0, MPI_COMM_WORLD, NULL);
            //DEBUG_CMD(printf("[rank: %d] Recv grid with size %ld\n", proc_id, dimt*block_size_copy);)
            root_grid_copy(grid_root, grid_tmp, dimt, dimx, block_size_copy, offset_copy);
            grid_free(grid_tmp);
        }
        DEBUG_CMD(stdout_dump(grid_root, dimt, dimx, block_size, offset, proc_num);)
        //file_dump(grid_root, dimt, dimx);
        grid_free(grid);
    }
    MPI_Finalize();
    return 0;
}
