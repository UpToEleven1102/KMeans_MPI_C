#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int rank, num_procs;

    printf("before \n");
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    printf("Process %i - num process: %i\n", rank, num_procs);

    MPI_Finalize();

    printf("after \n");
}