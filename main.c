#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "lib/Util.h"
#include "lib/KMeans_MPI.h"
#include "lib/KMeans_Search_MPI.h"

#define dim 4
#define n_data 100

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

    double *cluster_radius = (double*)malloc(k*sizeof(double));

    //initialize data set
    data = (double*) malloc(dim*n_data*sizeof(double));

    for (int i = 0; i < n_data * dim; ++i) {
        data[i] = (double)(rand()%10)/10;
    }

    //start K-means
    k = KMeans(dim, n_data, k, rank, num_procs, data, cluster_start, cluster_size, cluster_centroids, cluster_radius);

    double *result_pt = (double*)malloc(dim*sizeof(double));
    double *query = (double*) malloc(dim * sizeof(double));
    if (rank == 0) {
        for (int i = 0; i < dim; ++i) {
            query[i] = (double)(rand()%10)/10;
        }
    }

    MPI_Bcast(query, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < n_data; ++i) {
            printf("distance %d - %f\n", i, distanceOf2Points(dim, query, getElement(i, dim, data )));
        }
    }

    int counter = searchKMeans(dim, n_data, k, rank, num_procs, data, cluster_start, cluster_size,
            cluster_radius,
            cluster_centroids,
            query,
            result_pt);

    if(rank ==0) {
        printf("Closest point: counter %d -  distance %f\n",counter, distanceOf2Points(dim, query, result_pt));
        printArray(dim, result_pt);
    }

    for (int i = 0; i < k; ++i) {
        free(cluster_centroids[i]);
    }
    free(cluster_centroids);
    free(data);
    MPI_Finalize();
}