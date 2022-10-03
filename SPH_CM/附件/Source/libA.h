#pragma once



//mpich instance of mpi, not openmpis
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#ifdef _MSC_VER
#include<Windows.h>
#include"mpi.h"
//#pragma comment(lib,"C:\\Program Files (x86)\\MPI\\Lib\\x64\\msmpi.lib")
#endif

#ifdef __GNUC__
#include <unistd.h>
#include<mpich/mpi.h>
#endif


struct WORKER;
typedef struct WORKER* Worker;


int* WorkerRank(struct WORKER* worker);
int* WorkerNproc(struct WORKER* worker);
long long* WorkerNdata(struct WORKER* worker);
double** WorkerData(struct WORKER* worker);
int WorkerStart(struct WORKER ** result);
int WorkerEnd(struct WORKER * result);



void mpiwait(struct WORKER* wkr, int argc, char **argv);





