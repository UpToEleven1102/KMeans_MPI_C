#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    //init mpi
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // declare input variables
    srand(rank);
    const int DIM = 4;
    const int N_DATA = 100;

    double *data;

    data = (double*) malloc(DIM*N_DATA*sizeof(double));

    for (int i = 0; i < N_DATA * DIM; ++i) {
        data[i] = (double)(rand()%10)/10;
    }

    printf("%d \n", rank);
    for (int i = 0; i < 10; ++i) {
        printf("%f\n",data[i]);
    }




    printf("Process %i - num process: %i\n", rank, num_procs);

    MPI_Finalize();
}