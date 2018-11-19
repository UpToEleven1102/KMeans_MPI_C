//
// Created by Huyen V on 11/10/18.
//

#include <mpi.h>

#ifndef K_MEANS_UTILITIES_H
#define K_MEANS_UTILITIES_H

typedef enum {
    false, true
} bool;

typedef struct dataPoint dataPoint;
struct dataPoint {
    double *data;
};

typedef struct metaData metaData;
struct metaData {
    int idx;
    metaData *next;
};

int compareDataPoint(int dim, const double *point1, const double *point2);

int compareArrayToPointer(int dim, const double point1[], const double *point2);

double distanceOf2Points(int dim, double *point1, double *point2);

int printArray(int size, double *array);

void printDataSet(int dim, int ndata, double *data);

double square(double n);

double *getElement(int index, int dim, double *data);

void mpi_sendDoublePointer(int size, double *data, int rank, int target);

MPI_Status mpi_recvDoublePointer(int size, double *result, int rank, int source);

int mpi_bCastDoublePointer(int size, double *data, int rank);

int mpi_allGatherDoublePointer(int size, double *data, double *result);

#endif //KMEANS_PARALLEL_UTIL_H
