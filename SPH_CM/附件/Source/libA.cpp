#include "libA.h"




struct WORKER
{
    int rank;
    int nproc;
    double* data;
    long long ndata;
};

int* WorkerRank(struct WORKER* worker)
{
    return &(worker->rank);
}
int* WorkerNproc(struct WORKER* worker)
{
    return &(worker->nproc);
}
long long* WorkerNdata(struct WORKER* worker)
{
    return &(worker->ndata);
}
double** WorkerData(struct WORKER* worker)
{
    return &(worker->data);
}

int WorkerStart(struct WORKER ** result)
{   
    *result = (struct WORKER*)(malloc(sizeof(struct WORKER)));
    MPI_Comm_rank(MPI_COMM_WORLD,WorkerRank(*result));
    MPI_Comm_size(MPI_COMM_WORLD,WorkerNproc(*result));
    if(result)
        return 0;
    else 
        return 1;

}

int WorkerEnd(struct WORKER * result)
{
    free(result->data);
    free(result);
    return 0;
}


/*
void mpiwait(struct WORKER* wkr, int argc, char **argv)
{
    if(*WorkerRank(wkr)==0)
    {
        int ifdbg = 0;
        
        for(int i=0;i<argc;i++)
        {
                //printf("DEBUGGING/n");
                //fflush(stdout);
            if(! strcmp(argv[i],"debug") )
            {
                printf("DEBUGGING\n");
                ifdbg = 1;
                fflush(stdout);
            }
        }
        
        int i;
        for(i = 1; i<(*WorkerNproc(wkr)); i++)
        {
            MPI_Send(&ifdbg,1,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        if(ifdbg)
        {
            printf("Worker %d started waiting.\n",*WorkerRank(wkr));
        }
        fflush(stdout);
        MPI_Status stat;
        while(ifdbg)
        {
            /// SET BREAK POINT AT MPISEND
            MPI_Send(&ifdbg,1,MPI_INT,(*WorkerNproc(wkr))-1,0,MPI_COMM_WORLD);
            /// SET BREAK POINT AT MPISEND


            MPI_Recv(&ifdbg,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
            if(!ifdbg)
            {
                printf("Worker %d quit waiting.\n",*WorkerRank(wkr));
                fflush(stdout);

            }
        };
    }
    else
    {
        MPI_Status stat;
        int ifdbg;
        MPI_Recv(&ifdbg,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        //printf("Worker %d received ifdbg %llx.\n",*WorkerRank(wkr),(long long)ifdbg);
        //fflush(stdout);
        if(ifdbg)
        {
            printf("Worker %d started waiting.\n",*WorkerRank(wkr));
        }
        fflush(stdout);
        while(ifdbg)
        {
            MPI_Recv(&ifdbg,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

            /// SET BREAK POINT AT MPISEND
            MPI_Send(&ifdbg,1,MPI_INT,(*WorkerRank(wkr))-1,0,MPI_COMM_WORLD);
            /// SET BREAK POINT AT MPISEND

            if(!ifdbg)
            {
                printf("Worker %d quit waiting.\n",*WorkerRank(wkr));
                fflush(stdout);
            }
        }
    }
}
*/
void mpiwait(struct WORKER* wkr, int argc, char** argv)
{
    int ifdbg = 0;
    int i;

    for (int i = 0; i < argc; i++)
    {
        //printf("DEBUGGING/n");
        //fflush(stdout);
        if (!strcmp(argv[i], "debug"))
        {
            printf("DEBUGGING\n");
            ifdbg = 1;
            fflush(stdout);
        }
    }

    if (ifdbg)
    {
        printf("Worker %d started waiting.\n", *WorkerRank(wkr));
    }
    fflush(stdout);
    
    while (ifdbg)
    {
        for (i = 0; i < wkr->nproc; i++)
        {
            if (i == wkr->rank) 
            {
#ifdef _MSC_VER
                Sleep(100);
#endif

#ifdef __GNUC__
                sleep(1);
#endif
            }
            MPI_Bcast(&ifdbg,1,MPI_INT,i,MPI_COMM_WORLD);
            

        }
    }
}
//#include"libA.h" 
//
//
//#ifndef libA_DBG
//#define libA_DBG
//void mpiwait(struct WORKER* wkr, int argc, char** argv)
//{
//    if (*WorkerRank(wkr) == 0)
//    {
//        int ifdbg = 0;
//
//        for (int i = 0; i < argc; i++)
//        {
//            //printf("DEBUGGING/n");
//            //fflush(stdout);
//            if (!strcmp(argv[i], "debug"))
//            {
//                printf("DEBUGGING\n");
//                ifdbg = 1;
//                fflush(stdout);
//            }
//        }
//
//        int i;
//        for (i = 1; i < (*WorkerNproc(wkr)); i++)
//        {
//            MPI_Send(&ifdbg, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//        }
//        if (ifdbg)
//        {
//            printf("Worker %d started waiting.\n", *WorkerRank(wkr));
//        }
//        fflush(stdout);
//        MPI_Status stat;
//        while (ifdbg)
//        {
//            /// SET BREAK POINT AT MPISEND
//            MPI_Ssend(&ifdbg, 1, MPI_INT, (*WorkerNproc(wkr)) - 1, 0, MPI_COMM_WORLD);
//            /// SET BREAK POINT AT MPISEND
//
//
//            MPI_Recv(&ifdbg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
//            if (!ifdbg)
//            {
//                printf("Worker %d quit waiting.\n", *WorkerRank(wkr));
//                fflush(stdout);
//
//            }
//#ifdef _MSC_VER
//            Sleep(1000);
//#endif
//
//#ifdef __GNUC__
//            sleep(1);
//#endif
//            
//        };
//    }
//    else
//    {
//        MPI_Status stat;
//        int ifdbg;
//        MPI_Recv(&ifdbg, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
//        //printf("Worker %d received ifdbg %llx.\n",*WorkerRank(wkr),(long long)ifdbg);
//        //fflush(stdout);
//        if (ifdbg)
//        {
//            printf("Worker %d started waiting.\n", *WorkerRank(wkr));
//        }
//        fflush(stdout);
//        while (ifdbg)
//        {
//            MPI_Recv(&ifdbg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
//
//            /// SET BREAK POINT AT MPISEND
//            MPI_Ssend(&ifdbg, 1, MPI_INT, (*WorkerRank(wkr)) - 1, 0, MPI_COMM_WORLD);
//            /// SET BREAK POINT AT MPISEND
//
//            if (!ifdbg)
//            {
//                printf("Worker %d quit waiting.\n", *WorkerRank(wkr));
//                fflush(stdout);
//            }
//#ifdef _MSC_VER
//            Sleep(1000);
//#endif
//
//#ifdef __GNUC__
//            sleep(1);
//#endif
//        }
//    }
//}


