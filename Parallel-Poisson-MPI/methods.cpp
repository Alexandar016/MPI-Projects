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

void constructCoefficientMatrix(double* A, int n) {
    // Fill the array with zeros
    for (int i = 0; i < n * n; ++i) {
        A[i] = 0.0;
    }

    // Fill the diagonal block (T)
    for (int i = 0; i < n; ++i) {
        A[i * n + i] = 4.0; // Diagonal element
        if (i > 0) A[i * n + i - 1] = -1.0; // Lower diagonal
        if (i < n - 1) A[i * n + i + 1] = -1.0; // Upper diagonal
    }

    // Fill the off-diagonal blocks (-I)
    for (int i = 0; i < n - 1; ++i) {
        A[i * n + i + 1] = -1.0;
        A[(i + 1) * n + i] = -1.0;
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
    MPI_Barrier(MPI_COMM_WORLD);
}
