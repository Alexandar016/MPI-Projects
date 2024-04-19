/* Student: Aleksandart Mitic */

#include "methods.h"
#include <iostream>
#include <cstdio>
#include <math.h>
#include <cmath>

void read_distribute_input(int p, int my_rank,float& N,float& M) {
    float *tmp_arr;
    if (my_rank == 0) {
        FILE* inputFile = fopen("input.txt", "r");
        if (inputFile == nullptr) {
            std::cerr << "Error: Unable to open input file." << std::endl;
            return; 
        }
        fscanf(inputFile, "%f", &N);
        fscanf(inputFile, "%f", &M);
        std::cout << "Size N is: " << N << " and size M: " << M << std::endl;
        fclose(inputFile);
    }
}

void create_grid(MPI_Comm *cart_comm, int *coords, int *dims,int p,int my_rank) {
    // Calculate the dimensions of the grid
    MPI_Dims_create(p, 2, dims); // Assuming a 2D grid

    // Create a Cartesian communicator
    int periods[2] = {1, 1}; // Wrap around in both dimensions
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, cart_comm);

    // Get the coordinates of this process in the grid
    MPI_Cart_coords(*cart_comm, my_rank, 2, coords);

    // Print out the coordinates of this process
    std::cout << "Process " << my_rank << " is at position (" << coords[0] << ", " << coords[1] << ")" << std::endl;
}
