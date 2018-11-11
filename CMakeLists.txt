cmake_minimum_required(VERSION 3.12)
project(KMeans_parallel C)

set(CMAKE_C_STANDARD 99)

add_executable(KMeans_parallel main.c lib/Util.c lib/Util.h)