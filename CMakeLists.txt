cmake_minimum_required(VERSION 3.12)
project(KMeans_parallel C)

set(CMAKE_C_STANDARD 99)

add_executable(KMeans_parallel main.c lib/Util.c lib/Util.h lib/KMeans_MPI.c lib/KMeans_MPI.h lib/KMeans_Search_MPI.c lib/KMeans_Search_MPI.h)