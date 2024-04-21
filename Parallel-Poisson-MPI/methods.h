/* Student: Aleksandart Mitic */
#ifndef METHODS_H
#define METHODS_H

#include "mpi.h"

double f(double x, double y);
double g(double x, double y);
void create_grid(MPI_Comm *cart_comm, int *coords, int *dims,int p,int my_rank);
void constructCoefficientMatrix(double* A, int n);

#endif