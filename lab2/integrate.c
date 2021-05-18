#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>

#define NDEBUG

#ifdef NDEBUG
#define DEBUG_CMD(cmd);
#else
#define DEBUG_CMD(cmd) cmd;
#endif

typedef struct {
    long double A;
    long double B;
    long double fA;
    long double fB;
    long double s;
} task_t;

typedef struct _thread_data_t {
    int tid;
} thread_data_t;

sem_t sem_res;
sem_t sem_task;
sem_t sem_task_present;

typedef void* thr_sub(void*);
pthread_t *tid_procs;

task_t* tasks;
long task_n;
long task_max_n = 100;
long task_max_n_per_thread = 10;
long active_thread_n;
long double res;

const long double x1 = 0.001;
const long double x2 = 1;
const long double acc = 1e-5;
const long num_threads = 4;

long double f(long double x) {
    return sinl(1.0 / x);
}

int break_cond(long double sACB, long double sAB) {
    return fabs(sAB - sACB) < acc * fabs(sACB);
}

task_t* tasks_init(long* _task_n, long _task_max_n) {
    task_t* _tasks = calloc(_task_max_n, sizeof(task_t));
    *_task_n = 0;
    return _tasks;
}

void get_task(task_t* _tasks, long* _task_n, long double* _A, long double* _B, long double* _fA, long double* _fB, long double* _s) {
    assert(*_task_n > 0);
    (*_task_n)--;
    *_A = _tasks[*_task_n].A;
    *_B = _tasks[*_task_n].B;
    *_fA = _tasks[*_task_n].fA;
    *_fB = _tasks[*_task_n].fB;
    *_s = _tasks[*_task_n].s;
}

void put_task(task_t* _tasks, long* _task_n, long* _task_max_n, long double _A, long double _B, long double _fA, long double _fB, long double _s) {
    assert(*_task_n <= *_task_max_n);
    if (*_task_n == *_task_max_n) {
        *_task_max_n *= 2;
        _tasks = realloc(_tasks, *_task_max_n * sizeof(task_t));
    }
    _tasks[*_task_n].A = _A;
    _tasks[*_task_n].B = _B;
    _tasks[*_task_n].fA = _fA;
    _tasks[*_task_n].fB = _fB;
    _tasks[*_task_n].s = _s;
    (*_task_n)++;
}

void tasks_free(task_t* _tasks, long* _task_n) {
    free(_tasks);
    *_task_n = 0;
}

void* thread_job(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    int tid = data->tid;
    long task_n_local;
    long task_max_n_local = 100;
    long double s_local = 0;
    task_t* tasks_local = tasks_init(&task_n_local, task_max_n_local);
    while(1) {
        sem_wait(&sem_task_present);
        sem_wait(&sem_task);
        long double A, B, fA, fB, sAB;
        DEBUG_CMD(fprintf(stderr, "[thread %d] Get global task: ", tid));
        get_task(tasks, &task_n, &A, &B, &fA, &fB, &sAB);
        DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", A, B, fA, fB, sAB));
        if (task_n)
            sem_post(&sem_task_present);
        
        if (A <= B)
            active_thread_n++;

        sem_post(&sem_task);

        if (A > B)
            break;
        
        while (1) {
            long double C = (A + B) / 2;
            long double fC = f(C);
            long double sAC = (fA + fC) * (C - A) / 2;
            long double sCB = (fC + fB) * (B - C) / 2;
            long double sACB = sAC + sCB;
            if (break_cond(sACB, sAB)) {
                s_local += sACB;
                if (!task_n_local)
                    break;

                DEBUG_CMD(fprintf(stderr, "[thread %d] Get local task: ", tid));
                get_task(tasks_local, &task_n_local, &A, &B, &fA, &fB, &sAB);
                DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", A, B, fA, fB, sAB));
            }
            else {
                DEBUG_CMD(fprintf(stderr, "[thread %d] Put local task: ", tid));
                put_task(tasks_local, &task_n_local, &task_max_n_local, A, C, fA, fC, sAC);
                DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", A, C, fA, fC, sAC));
                A = C;
                fA = fC;
                sAB = sCB;
            }
            if ((task_n_local > task_max_n_per_thread) && (!task_n)) {
                sem_wait(&sem_task);
                if (!task_n)
                    sem_post(&sem_task_present);

                while ((task_n_local > 1) && (task_n < task_max_n)) {
                    DEBUG_CMD(fprintf(stderr, "[thread %d] Get local task: ", tid));
                    get_task(tasks_local, &task_n_local, &A, &B, &fA, &fB, &sAB);
                    DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", A, B, fA, fB, sAB));
                    DEBUG_CMD(fprintf(stderr, "[thread %d] Put global task: ", tid));
                    put_task(tasks, &task_n, &task_max_n, A, B, fA, fB, sAB);
                    DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", A, B, fA, fB, sAB));
                }

                sem_post(&sem_task);
            }
        }

        sem_wait(&sem_task);
        active_thread_n--;

        if (!active_thread_n && !task_n) {
            for (int i = 0; i < num_threads; i++) {
                DEBUG_CMD(fprintf(stderr, "[thread %d] Put global task: ", tid));
                put_task(tasks, &task_n, &task_max_n, 2, 1, 0, 0, 0);
                DEBUG_CMD(fprintf(stderr, "A = 2 B = 1 fA = 0 fB = 0 s = 0\n"));
            }
            sem_post(&sem_task_present);
        }
        sem_post(&sem_task);
    }
    sem_wait(&sem_res);
    res += s_local;
    sem_post(&sem_res);
    tasks_free(tasks_local, &task_n_local);
}

void start_threads(thr_sub* tsub) {
    tid_procs = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&(tid_procs[i]), NULL, tsub, &i);
    }
}

void wait_threads() {
    for (int i = 0; i < num_threads; i++) {
        pthread_join(tid_procs[i], NULL);
    }
}

int main() {
    tasks = tasks_init(&task_n, task_max_n);
    active_thread_n = 0;
    sem_init(&sem_res, 0, 1);
    sem_init(&sem_task, 0, 0);
    sem_init(&sem_task_present, 0, 0);
    res = 0;
    DEBUG_CMD(fprintf(stderr, "[thread main] Put global task: "));
    put_task(tasks, &task_n, &task_max_n, x1, x2, f(x1), f(x2), (f(x1) + f(x2)) * (x2 - x1) / 2);
    DEBUG_CMD(fprintf(stderr, "A = %Lf B = %Lf fA = %Lf fB = %Lf s = %Lf\n", x1, x2, f(x1), f(x2), (f(x1) + f(x2)) * (x2 - x1) / 2));
    sem_post(&sem_task);
    sem_post(&sem_task_present);
    start_threads(thread_job);
    wait_threads();
    sem_destroy(&sem_res);
    sem_destroy(&sem_task);
    sem_destroy(&sem_task_present);
    printf("Result: %.20Lf\n", res);
    return 0;
}
