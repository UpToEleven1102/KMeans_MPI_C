//
// Created by h on 11/10/18.
//

#include <mpi.h>
#include "Util.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>

double square(double n){
    return n*n;
}

int compareDataPoint(int dim, const double *point1, const double *point2) {
    for (int i = 0; i < dim; ++i) {
        if (point1[i] != point2[i]) {
            return -1;
        }
    }
    return 0;
}

int compareArrayToPointer(int dim, const double point1[], const double *point2){
    for (int i = 0; i < dim; ++i) {
        if (point1[i] != point2[i]) {
            return -1;
        }
    }
    return 0;
}

double distanceOf2Points(int dim, double *point1, double *point2) {
    double distance = 0;
    for (int i = 0; i < dim; ++i) {
        distance += square(point1[i] - point2[i]);
    }
    return sqrt(distance);
}

int printArray(int size, double *array){
    for (int i = 0; i < size; ++i) {
        printf("%f \n",array[i]);
    }
}

double *getElement(int index, int dim, double *data){
    double *element = (double*)malloc(dim*sizeof(double));
    for (int i = index * dim; i < index*dim + dim ; ++i) {
        element[i - index*dim] = data[i];
    }
    return element;
}

void printDataSet(int dim, int ndata, double *data) {
    for (int i = 0; i < ndata; ++i) {
        printf("--%d--\n", i);
        printArray(dim, getElement(i, dim, data));
    }
}

void mpi_sendDoublePointer(int size, double *data, int rank, int target){
    for (int i = 0; i < size; ++i) {
        MPI_Send(&data[i], 1, MPI_DOUBLE, target, rank, MPI_COMM_WORLD);
    }
}

MPI_Status mpi_recvDoublePointer(int size, double *result, int rank, int source){
    MPI_Status status;
    for(int i = 0; i < size; ++i) {
        MPI_Recv(&result[i], 1, MPI_DOUBLE, source, source, MPI_COMM_WORLD, &status);
    }
    return status;
}

int mpi_allGatherDoublePointer(int size, double *data, double *result) {
    double array[size];
    for (int i = 0; i < size; ++i) {
        array[i] = data[i];
    }
    MPI_Allgather(&array, size, MPI_DOUBLE, result, size, MPI_DOUBLE, MPI_COMM_WORLD);
}