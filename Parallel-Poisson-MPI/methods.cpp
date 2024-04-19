/* Student: Aleksandart Mitic */

#include "methods.h"
#include <iostream>
#include <cstdio>
#include <math.h>
#include <cmath>

double f(double x, double y) {
    return 2 * ((1 + x) * sin(x + y) - cos(x + y));
}

double g(double x, double y) {
    return (1 + x) * sin(x + y);
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
    MPI_Barrier(MPI_COMM_WORLD);
}
