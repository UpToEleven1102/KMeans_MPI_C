#!/bin/bash
mpicc main.c -o kmeans
mpirun -np 4 kmeans

