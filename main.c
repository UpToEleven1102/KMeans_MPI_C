#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "lib/Util.h"
#include "lib/KMeans_MPI.h"

#define dim 4
#define n_data 25

int main(int argc, char **argv) {
    //init mpi
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Status status;

    // declare input variables
    srand(rank);

    //number of clusters
    int k = 8;


    double *data;

    double **cluster_centroids = (double**)malloc(k*sizeof(double*));
    for (int i = 0; i < k; ++i) {
        cluster_centroids[i] = (double*)malloc(dim*sizeof(double));
    }

    int *cluster_start, *cluster_size;
    cluster_start = (int*)malloc(k * sizeof(int));
    cluster_size = (int*)malloc(k*sizeof(int));

    //initialize data set
    data = (double*) malloc(dim*n_data*sizeof(double));

    for (int i = 0; i < n_data * dim; ++i) {
        data[i] = (double)(rand()%10)/10;
    }

    KMeans(dim, n_data, k, rank, num_procs, data, cluster_start, cluster_size, cluster_centroids);




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