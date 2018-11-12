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
    double *minDistances = (double *) malloc(n_data * sizeof(double));
    for (int i = 0; i < n_data; ++i) {
        minDistances[i] = 999;
    }
    double *datum = NULL;
    double distance;
    double maxDistance = 0;
    int max_idx = 0;

//    printf("Debug %i", k);

    for (int i = 0; i < n_data; ++i) {
        dataPoint *element = getElement(i, dim, data);
        datum = element->data;
        if (isCentroid(dim, k, datum, cluster_centroids)) {
            continue;
        }
        for (int j = 0; j < k; ++j) {
            distance = distanceOf2Points(dim, datum, cluster_centroids[j]);
            if (minDistances[i] > distance)
                minDistances[i] = distance;
        }
        if (minDistances[i] > maxDistance) {
            maxDistance = minDistances[i];
            max_idx = i;
        }
        free(element);
        free(datum);
    }


    //buff globally store data being shared by all processes
    double buff[dim + 1];
    dataPoint *element = getElement(max_idx, dim, data);
    datum = element->data;
    free(element);
    for (int i = 0; i < dim; ++i) {
        buff[i] = datum[i];
    }
    buff[dim] = maxDistance;

    double *buffer = (double *) malloc(num_procs * (dim + 1) * sizeof(double));
    MPI_Allgather(&buff, dim + 1, MPI_DOUBLE, buffer, dim + 1, MPI_DOUBLE, MPI_COMM_WORLD);


    maxDistance = buffer[dim];
    for (int i = 0; i < dim; ++i) {
        datum[i] = buffer[i];
    }

    dim = dim + 1;
    for (int i = 1; i < num_procs; ++i) {
        if (maxDistance < buffer[i * dim + dim - 1]) {
            maxDistance = buffer[i * dim + dim - 1];
            for (int j = 0; j < dim - 1; ++j) {
                datum[j] = buffer[i * dim + j];
            }
        }
    }

    for (int i = 0; i < dim; ++i) {
        cluster_centroids[k][i] = datum[i];
    }
    free(buffer);
    free(minDistances);
    free(datum);
}

int initializeCentroids(int k, int dim, int n_data, int rank, int num_procs, double *data, double **cluster_centroids) {
    int counter = 1;

    //generate first centroid
    if (rank == 0) {
        int idx = rand() % n_data;
        dataPoint *ele = getElement(idx, dim, data);
        double *datum = ele->data;
        for (int i = 0; i < dim; ++i) {
            cluster_centroids[0][i] = datum[i];
        }
        free(datum);
        free(ele);
    }

    mpi_bCastDoublePointer(dim, cluster_centroids[0], 0);

    //find next centroids
    while (counter < k) {
        nextCentroid(counter, dim, n_data, rank, num_procs, data, cluster_centroids);
        counter++;
    }

//    if (rank == 1) {
//        for (int i = 0; i < k; ++i) {
//            for (int j = 0; j < dim; ++j) {
//                printf("%f\n", cluster_centroids[i][j]);
//            }
//        }
//    }
}

double *newCentroid(int dim, int cluster_size, metaData *cluster_info, double *data) {
    double *centroid = (double *) malloc(dim * sizeof(double));
    for (int i = 0; i < dim; ++i) {
        centroid[i] = 0;
    }

    while (cluster_info != NULL) {
        dataPoint *element = getElement(cluster_info->idx, dim, data);
        for (int i = 0; i < dim; ++i) {
            centroid[i] = centroid[i] + element->data[i];
        }
        free(element->data);
        free(element);
        cluster_info = cluster_info->next;
    }

    for (int i = 0; i < dim; ++i) {
        centroid[i] = centroid[i] / cluster_size;
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
         datum = getElement(i, dim, data) -> data;
        double distance = distanceOf2Points(dim, datum, cluster_centroid);
        if (distance > result) {
            result = distance;
        }
        free(datum);
    }
    return result;
}

int rearrangeData(int dim,int n_data, int n_cluster, double **cluster_centroids,double *data, metaData **cluster_info, int* cluster_size, int *cluster_start, double *cluster_radius){
    double *cluster_assign = (double*) malloc(dim * n_data * sizeof(double));
    int idx = 0;
    cluster_start[0] = 0;
    double *datum;
    for (int i = 0; i < n_cluster; ++i) {
        if( i < n_cluster -1)
            cluster_start[i+1] = cluster_start[i] + cluster_size[i];
        while(cluster_info[i] !=NULL) {
            dataPoint *element = getElement(cluster_info[i]->idx, dim, data);
            datum= element->data;
            for (int j = 0; j < dim; ++j) {
                cluster_assign[idx * dim +j] = datum[j];
            }
            free(element); free(datum);
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
    int notEmpty = 0;
    double *datum = NULL;
    double localNewCentroid[dim];
    double *globalNewCentroids = NULL;
    int *buffer_int;

    while (!noChange) {
        noChange = true;
        for (int i = 0; i < n_cluster; ++i) {
            cluster_info[i] = NULL;
            cluster_size[i] = 0;
        }

        for (int i = 0; i < n_data; ++i) {
            dataPoint *element = getElement(i, dim, data);
            datum = element->data;
            distance = 9999;

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
            free(element);
        }

//        if (rank == 1)
//        for (int i = 0; i < n_cluster; ++i) {
//            printf("cluster - %d\n", i);
//            while(cluster_info[i] != NULL) {
//                printf("idx - %d\n", cluster_info[i]->idx);
//                cluster_info[i] = cluster_info[i]->next;
//            }
//        }

        for (int i = 0; i < n_cluster; ++i) {
            notEmpty = 0;
            if (cluster_info[i] != NULL) {
                notEmpty = 1;
                //inform other processes
            }

            buffer_int = (int *) malloc(num_procs * sizeof(int));

            MPI_Allgather(&notEmpty, 1, MPI_INT, buffer_int, 1, MPI_INT, MPI_COMM_WORLD);
//
//            if (rank == 1){
//                for (int i = 0; i < num_procs; ++i) {
//                    printf("buffer - %d", buffer_int[i]);
//                }
//                printf("\n");
//            }

            notEmpty = 0;
            for (int j = 0; j < num_procs; ++j) {
                if (buffer_int[j] == 1) {
                    notEmpty = 1;
                    break;
                }
            }

            free(buffer_int);

            if (notEmpty == 0) {
                n_cluster--;
                removeCluster(n_cluster, i, cluster_info, cluster_size);
                i--;
            } else {
                if (cluster_info[i] != NULL) {
                    //dont' want to allocate a temp variable, use this to hold temp value
                    globalNewCentroids = newCentroid(dim, cluster_size[i], cluster_info[i], data);

                    for (int j = 0; j < dim; ++j) {
                        localNewCentroid[j] = globalNewCentroids[j];
                    }
                    free(globalNewCentroids);
                } else {
                    for (int j = 0; j < dim; ++j) {
                        localNewCentroid[j] = 0;
                    }
                }

                //broadcast new centroid to global variable pool
                globalNewCentroids = (double *) malloc(dim * num_procs * sizeof(double));
                MPI_Allgather(&localNewCentroid, dim, MPI_DOUBLE, globalNewCentroids, dim, MPI_DOUBLE, MPI_COMM_WORLD);

                //calculate new centroid
                for (int j = 0; j < dim; ++j) {
                    localNewCentroid[j] = 0;
                }

                for (int j = 0; j < num_procs; ++j) {
                    for (int l = 0; l < dim; ++l) {
                        localNewCentroid[l] += globalNewCentroids[j * dim + l];
                    }
                }

                for (int j = 0; j < dim; ++j) {
                    localNewCentroid[j] /= num_procs;
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

    rearrangeData(dim, n_data, n_cluster, cluster_centroids, data, cluster_info, cluster_size, cluster_start, cluster_radius);
//    if (rank == 1) {
//        printf("debug cluster: %d - new centroids: \n", n_cluster);
//        for (int i = 0; i < n_cluster; ++i) {
//            for (int j = 0; j < dim; ++j) {
//                printf("%f\n", cluster_centroids[i][j]);
//            }
//            printf("cluster %d size %d\n", i, cluster_size[i]);
//        }
//
//    }

    return n_cluster;

}