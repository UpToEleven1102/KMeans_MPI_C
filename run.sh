#!/bin/bash
mpicc main.c lib/Util.h lib/Util.c  lib/KMeans_MPI.h lib/KMeans_MPI.c -o kmeans -lm
mpirun -np 4 kmeans
