/* Student: Aleksandart Mitic */
#ifndef METHODS_H
#define METHODS_H

#include "mpi.h"

void create_grid(MPI_Comm *cart_comm, int *coords, int *dims,int p,int my_rank);
void read_distribute_input(int p, int my_rank, float& N, float& M);

#endif