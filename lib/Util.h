//
// Created by Huyen V on 11/10/18.
//

#ifndef KMEANS_PARALLEL_UTIL_H
#define KMEANS_PARALLEL_UTIL_H
void mpi_sendDoublePointer(int size, double *data, int rank);
MPI_Status mpi_recvDoublePointer(int size, double *result, int rank, int source){

#endif //KMEANS_PARALLEL_UTIL_H