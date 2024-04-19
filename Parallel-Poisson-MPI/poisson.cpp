/* Student: Aleksandart Mitic */

#include <stdio.h>
#include <iostream>
#include <cstdio> 
#include "methods.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
    /* Start the MPI environment*/
    MPI_Init(&argc, &argv);
    int my_rank,p;
    // float *local_nodes = (float*)malloc(N * sizeof(float));
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm cart_comm;
    int coords[2], dims[2];
    create_grid(&cart_comm, coords, dims,p,my_rank);
    MPI_Finalize();
    return 0;
}