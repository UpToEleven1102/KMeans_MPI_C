//
// Created by h on 11/11/18.
//

#include <mpi.h>
#include <malloc.h>
#include "KMeans_Search_MPI.h"
#include "Util.h"

double distanceToCluster(int dim, double *query, double *cluster_centroid, double cluster_radius){
    double distance = distanceOf2Points(dim, query,cluster_centroid);
    return distance - cluster_radius;
}

double *searchCluster(int dim, int *counter, double *query_pt, int cluster_start, int cluster_size, double *data) {
    double *result_pt = NULL;
    double *temp = NULL;
    double min_distance = 999;
    double distance;
    for (int i = cluster_start; i < cluster_start+cluster_size; ++i) {
        *counter = *counter+1;
        dataPoint *element = getElement(i, dim, data);
        temp = element->data;
        free(element);
        distance = distanceOf2Points(dim, query_pt, temp);
        if (distance < min_distance) {
            min_distance = distance;
            free(result_pt);
            result_pt = temp;
        }
    }
    free(temp);
    return result_pt;
}

int searchKMeans(int dim, int n_data, int k, int rank, int num_procs,
        double *data,
        int *cluster_start,
        int *cluster_size,
        double *cluster_radius,
        double **cluster_centroids,
        double *query_pt,
        double *result_pt){

    double *buffer = NULL;

    double *distances = (double *) malloc(k* sizeof(double));
    int cluster_idx = 0;
    distances[0] = distanceToCluster(dim, query_pt, cluster_centroids[0], cluster_radius[0]);
    double shortest_distance = distances[0];

    for (int i = 1; i < k; ++i) {
        distances[i] = distanceToCluster(dim, query_pt, cluster_centroids[i], cluster_radius[i]);
        if(distances[i] < shortest_distance) {
            shortest_distance = distances[i];
            cluster_idx = i;
        }
    }

    int *counter = (int*)malloc(sizeof(int));
    *counter = 0;

    double *closest_point = searchCluster(dim, counter, query_pt, cluster_start[cluster_idx],
            cluster_size[cluster_idx], data);

    shortest_distance = distanceOf2Points(dim, closest_point, query_pt);

    double distance = 0;
    for (int i = 0; i < k; ++i) {
        if (cluster_idx != i) {
            if(shortest_distance > distances[i]) {
                buffer = searchCluster(dim, counter, query_pt, cluster_start[i], cluster_size[i], data);
                distance = distanceOf2Points(dim, query_pt, buffer);
                if(distance < shortest_distance) {
                    shortest_distance = distance;
                    free(closest_point);
                    closest_point = buffer;
                }
                free(buffer);
            }
        }
    }

    buffer = (double*)malloc(dim*num_procs*sizeof(double));

    mpi_allGatherDoublePointer(dim, closest_point, buffer);

    double *temp = NULL;
    for (int i = 0; i < num_procs; ++i) {
        temp = (double*)malloc(dim* sizeof(double));
        for (int j = 0; j < dim; ++j) {
            temp[j] = buffer[i*dim + j];
            distance = distanceOf2Points(dim, temp, query_pt);
            if(distance < shortest_distance){
                shortest_distance = distance;
                free(closest_point);
                closest_point = temp;
            }
        }
    }

    for (int i = 0; i < dim; ++i) {
        result_pt[i] = closest_point[i];
    }

    free(temp);
    free(closest_point);
    free(distances);
    return *counter;
}