//
// Created by h on 11/11/18.
//
#include <mpi.h>
#ifndef KMEANS_PARALLEL_KMEANS_MPI_H
#define KMEANS_PARALLEL_KMEANS_MPI_H

int KMeans(int dim, int n_data, int k, int rank, int num_procs, double *data, int *cluster_start, int *cluster_size, double **cluster_centroids);

#endif //KMEANS_PARALLEL_KMEANS_MPI_H
