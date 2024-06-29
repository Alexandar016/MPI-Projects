#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --exclusive
#SBATCH --partition=compute2008
# exec <input.d
mpirun -np 4 ./a.out
