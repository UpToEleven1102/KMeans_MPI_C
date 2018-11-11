#!/bin/bash
mpicc main.c lib/Util.h  lib/Util.c -o kmeans
mpirun -np 4 kmeans
