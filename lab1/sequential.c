#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
const double h = 0.2; //2t >= h
const double t = 0.1;
const double a = 2;

typedef struct {
    double value;
} Cell;

typedef struct {
    Cell** cell;
    long dimx; 
    long dimt;
} Grid;

Grid* grid_init(long dimt, long dimx) {
    Cell** cell = calloc(dimt, sizeof(Cell*));
    for (int i = 0; i < dimt; i++) {
        cell[i] = calloc(dimx + dimt - 1 - i, sizeof(Cell));
    }
    Grid* grid = calloc(1, sizeof(Grid));
    grid->cell = cell;
    grid->dimx = dimx;
    grid->dimt = dimt; 
    return grid;
}

void grid_free(Grid* grid) {
    for (int i = 0; i < grid->dimt; i++) {
        free(grid->cell[i]);
    }
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

void calc_cell(Grid* grid, long ti, long xi, double t, double h)
{
    Cell** cell = grid->cell;
    if (ti == 0) {
        cell[ti][xi].value = f_init_x(xi*h);
        return;
    }
    if (xi == 0) {
        cell[ti][xi].value = f_init_t(ti*t);
        return;
    }
    double tmp1 = 0.5*(cell[ti-1][xi+1].value+cell[ti-1][xi-1].value);
    double tmp2 = -a*0.5*t/h*(cell[ti-1][xi+1].value-cell[ti-1][xi-1].value);
    cell[ti][xi].value = t*f(ti*t-t, xi*h)+tmp1+tmp2;
}

void file_dump(Grid* grid, long dimt, long dimx) {
    FILE* fout = fopen("grid_seq.csv", "w");
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

void stdout_dump(Grid* grid, long dimt, long dimx) {
    printf("Grid\n");
    for (int ti = dimt-1; ti > -1; ti--) {
        for (int xi = 0; xi < dimx + dimt - 1 - ti; xi++) {
            printf("|%5.2lf ", grid->cell[ti][xi].value);
        }
        printf("|\n");
    }
}

int main() {
    long dimt = (long)(ranget - initt)/t;
    long dimx = (long)(rangex - initx)/h;
    Grid* grid = grid_init(dimt, dimx);
    Cell** cell = grid->cell;
    for (int ti = 0; ti < dimt; ti++) {
        for (int xi = 0; xi < dimx + dimt - 1 - ti; xi++) {
            calc_cell(grid, ti, xi, t, h);
        }
    }
    DEBUG_CMD(stdout_dump(grid, dimt, dimx);)
    //file_dump(grid, dimt, dimx);
    grid_free(grid); 
    return 0;
}
