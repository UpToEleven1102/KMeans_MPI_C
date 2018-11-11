#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "lib/Util.h"

#define dim 4
#define n_data 100

int initializeCentroids(int k, int rank, double *data, double **cluster_centroids) {
    if (rank == 0) {
        mpi_sendDoublePointer(k, data, rank, 1);
    } else if(rank == 1) {
        double *result = (double*)malloc(k*sizeof(double));
        mpi_recvDoublePointer(k, result, rank, 0);
        for(int i = 0; i<k; i++) {
            printf("%f", result[i]);
        }
    } else if(rank == 2) {

    } else if(rank == 3) {

    }
}

int main(int argc, char **argv) {
    //init mpi
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Status status;

    // declare input variables
    srand(rank);

    int k = 8;


    double *data;
    double **cluster_centroids = (double**)malloc(k*sizeof(double*));
    for (int i = 0; i < k; ++i) {
        cluster_centroids[i] = (double*)malloc(dim*sizeof(double));
    }

    data = (double*) malloc(dim*n_data*sizeof(double));

    for (int i = 0; i < n_data * dim; ++i) {
        data[i] = (double)(rand()%10)/10;
    }

    initializeCentroids(k, rank, data, cluster_centroids);



    /*for (int i = 0; i < dim; ++i) {
        MPI_Send(&data[i], 1, MPI_DOUBLE, 1, i, MPI_COMM_WORLD);
    }
         MPI_Send( void *data, int send_count, MPI_Datatype send_type, int destination_ID, int tag, MPI_Comm comm);snip

        MPI_Recv(void *received_data, int receive_count, MPI_Datatype receive_type, int sender_ID, int tag, MPI_Comm comm, MPI_Status *status);
    */


    //free pointers
    for (int i = 0; i < k; ++i) {
        free(cluster_centroids[i]);
    }
    free(cluster_centroids);
    free(data);
    MPI_Finalize();
}