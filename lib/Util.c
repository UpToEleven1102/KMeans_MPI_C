//
// Created by h on 11/10/18.
//

#include <arpa/nameser.h>
#include "Util.h"

void mpi_sendDoublePointer(int size, double *data, int rank, int target){
    for (int i = 0; i < size; ++i) {
        MPI_Send(&data[i], 1, MPI_DOUBLE, target, rank, MPI_COMM_WORLD);
    }
}

MPI_Status mpi_recvDoublePointer(int size, double *result, int rank, int source){
    MPI_Status status;
    for (int i = 0; i < size; ++i) {
        MPI_Recv(&result[i], 1, MPI_DOUBLE, source, source, MPI_COMM_WORLD, &status);
    }
    return status;
}