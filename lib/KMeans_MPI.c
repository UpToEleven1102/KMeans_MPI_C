//
// Created by h on 11/11/18.
//

#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <stdlib.h>
#include "KMeans_MPI.h"
#include "Util.h"

bool isCentroid(int dim, int k, double *datum, double **cluster_centroids) {
    for (int i = 0; i < k; ++i) {
        if (compareDataPoint(dim, datum, cluster_centroids[i]) == 0) {
            return true;
        }
    }
    return false;
}

int nextCentroid(int k, int dim, int n_data, int rank, int num_procs, double *data, double **cluster_centroids) {
    //array save min distance of each data to all clusters
    double *minDistances = (double *) malloc(n_data * sizeof(double));
    for (int i = 0; i < n_data; ++i) {
        minDistances[i] = RAND_MAX;
    }
    double *datum = NULL;
    double distance;
    double maxDistance = 0;
    int max_idx = 0;

    //find max of mins
    for (int i = 0; i < n_data; ++i) {
        datum = getElement(i, dim, data);
        if (isCentroid(dim, k, datum, cluster_centroids)) {
            continue;
        }
        //k is current number of clusters
        for (int j = 0; j < k; ++j) {
            distance = distanceOf2Points(dim, datum, cluster_centroids[j]);
            if (minDistances[i] > distance)
                minDistances[i] = distance;
        }

        if (minDistances[i] > maxDistance) {
            maxDistance = minDistances[i];
            max_idx = i;
        }
        free(datum);
    }

    //buff saves info of local max of mins
    double buff[dim + 1];
    datum = getElement(max_idx, dim, data);
    for (int i = 0; i < dim; ++i) {
        buff[i] = datum[i];
    }
    buff[dim] = maxDistance;

    //buffer save info of all max of mins of all processes
    double *buffer = (double *) malloc(num_procs * (dim + 1) * sizeof(double));
    MPI_Allgather(buff, dim + 1, MPI_DOUBLE, buffer, dim + 1, MPI_DOUBLE, MPI_COMM_WORLD);

    maxDistance = buffer[dim];
    for (int i = 0; i < dim; ++i) {
        datum[i] = buffer[i];
    }

    //find max of all maxs
    dim = dim + 1;
    for (int i = 1; i < num_procs; ++i) {
        if (maxDistance < buffer[i * dim + dim - 1]) {
            maxDistance = buffer[i * dim + dim - 1];
            for (int j = 0; j < dim - 1; ++j) {
                datum[j] = buffer[i * dim + j];
            }
        }
    }

    cluster_centroids[k] = datum;
    free(buffer);
    free(minDistances);
}

int initializeCentroids(int k, int dim, int n_data, int rank, int num_procs, double *data, double **cluster_centroids) {
    int counter = 1;

    //generate first centroid
    if (rank == 0) {
        cluster_centroids[0] = getElement(rand() % n_data, dim, data);
    }

    MPI_Bcast(cluster_centroids[0], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    //find next centroids
    while (counter < k) {
        nextCentroid(counter, dim, n_data, rank, num_procs, data, cluster_centroids);
        counter++;
    }
}

double *sumElements(int dim, int cluster_size, metaData *cluster_info, double *data) {
    double *centroid = (double *) malloc(dim * sizeof(double));
    for (int i = 0; i < dim; ++i) {
        centroid[i] = 0;
    }

    while (cluster_info != NULL) {
        double *element = getElement(cluster_info->idx, dim, data);
        for (int i = 0; i < dim; ++i) {
            centroid[i] = centroid[i] + element[i];
        }
        free(element);
        cluster_info = cluster_info->next;
    }

    return centroid;
}

int removeCluster(int n_cluster, int idx, metaData **cluster_info, int *cluster_size) {
    for (int i = idx; i < n_cluster; ++i) {
        cluster_info[i] = cluster_info[i + 1];
        cluster_size[i] = cluster_size[i + 1];
    }
}

double findClusterRadius(int dim, int cluster_start, int cluster_size, double *cluster_centroid, double *data) {
    double result = 0;
    double *datum;
    for (int i = cluster_start; i < cluster_start + cluster_size; ++i) {
        datum = getElement(i, dim, data);
        double distance = distanceOf2Points(dim, datum, cluster_centroid);
        if (distance > result) {
            result = distance;
        }
        free(datum);
    }
    return result;
}

int rearrangeData(int dim, int n_data, int n_cluster, double **cluster_centroids, double *data, metaData **cluster_info,
                  int *cluster_size, int *cluster_start, double *cluster_radius) {
    double *cluster_assign = (double *) malloc(dim * n_data * sizeof(double));
    int idx = 0;
    cluster_start[0] = 0;

    double *datum;

    for (int i = 0; i < n_cluster; ++i) {
        if (i < n_cluster - 1)
            cluster_start[i + 1] = cluster_start[i] + cluster_size[i];
        while (cluster_info[i] != NULL) {
            datum = getElement(cluster_info[i]->idx, dim, data);
            for (int j = 0; j < dim; ++j) {
                cluster_assign[idx * dim + j] = datum[j];
            }
            free(datum);
            idx++;
            cluster_info[i] = cluster_info[i]->next;
        }
    }

    for (int i = 0; i < dim * n_data; ++i) {
        data[i] = cluster_assign[i];
    }

    for (int i = 0; i < n_cluster; ++i) {
        cluster_radius[i] = findClusterRadius(dim, cluster_start[i], cluster_size[i], cluster_centroids[i], data);

    }
    free(cluster_assign);
}


int KMeans(int dim, int n_data, int k, int rank, int num_procs,
           double *data,
           int *cluster_start,
           int *cluster_size,
           double **cluster_centroids, double *cluster_radius) {

    initializeCentroids(k, dim, n_data, rank, num_procs, data, cluster_centroids);

    metaData **cluster_info = (metaData **) malloc(k * sizeof(metaData *));

    int n_cluster = k;
    int cluster_idx = 0;
    double distance;
    bool noChange = false;
    int isEmpty = 0;
    double *datum = NULL;
    double localNewCentroid[dim + 1];
    double *globalNewCentroids = NULL;
    int *buffer_int;

    while (!noChange) {
        noChange = true;
        for (int i = 0; i < n_cluster; ++i) {
            cluster_info[i] = NULL;
            cluster_size[i] = 0;
        }

        for (int i = 0; i < n_data; ++i) {
            datum = getElement(i, dim, data);
            distance = RAND_MAX;

            for (int j = 0; j < n_cluster; ++j) {
                double tempDistance = distanceOf2Points(dim, datum, cluster_centroids[j]);
                if (distance > tempDistance) {
                    distance = tempDistance;
                    cluster_idx = j;
                }
            }

            //found cluster idx for data at index i, register to clusterInfo
            cluster_size[cluster_idx]++;
            metaData *dataInfo = (metaData *) malloc(sizeof(metaData));
            dataInfo->idx = i;

            if (cluster_info[cluster_idx] == NULL) {
                cluster_info[cluster_idx] = dataInfo;
                dataInfo->next = NULL;
            } else {
                dataInfo->next = cluster_info[cluster_idx];
                cluster_info[cluster_idx] = dataInfo;
            }

            free(datum);
        }

        for (int i = 0; i < n_cluster; ++i) {
            isEmpty = 0;
            if (cluster_info[i] != NULL) {
                isEmpty = 1;
                //inform other processes
            }

            MPI_Allreduce(&isEmpty, &isEmpty, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (isEmpty == 0) {
                n_cluster--;
                removeCluster(n_cluster, i, cluster_info, cluster_size);
                i--;
            } else {
                if (cluster_info[i] != NULL) {
                    datum = sumElements(dim, cluster_size[i], cluster_info[i], data);

                    for (int j = 0; j < dim; ++j) {
                        localNewCentroid[j] = datum[j];
                    }
                    free(datum);
                } else {
                    for (int j = 0; j < dim; ++j) {
                        localNewCentroid[j] = 0;
                    }
                }
                localNewCentroid[dim] = cluster_size[i];

                //broadcast new centroid to global variable pool
                globalNewCentroids = (double *) malloc((dim + 1) * num_procs * sizeof(double));
                MPI_Allgather(localNewCentroid, dim + 1, MPI_DOUBLE, globalNewCentroids, (dim + 1), MPI_DOUBLE,
                              MPI_COMM_WORLD);
                //calculate new centroid
                for (int j = 0; j < dim; ++j) {
                    localNewCentroid[j] = 0;
                }
                int globalClusterSize = 0;
                for (int j = 0; j < num_procs; ++j) {
                    for (int l = 0; l < dim; ++l) {
                        localNewCentroid[l] += globalNewCentroids[j * (dim + 1) + l];
                    }
                    globalClusterSize += (int) (globalNewCentroids[j * (dim + 1) + dim]);
                }

                for (int j = 0; j < dim; ++j) {
                    localNewCentroid[j] /= globalClusterSize;
                }
                if (compareArrayToPointer(dim, localNewCentroid, cluster_centroids[i]) != 0) {
                    for (int j = 0; j < dim; ++j) {
                        cluster_centroids[i][j] = localNewCentroid[j];
                        noChange = false;
                    }
                }
            }
        }
    }

    rearrangeData(dim, n_data, n_cluster, cluster_centroids, data, cluster_info, cluster_size, cluster_start,
                  cluster_radius);

    return n_cluster;

}