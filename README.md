# Parallel MPI Projects

This project contains three separate implementations focusing on parallel numerical methods using MPI:

1. **Parallel Implementation of Poisson Equation using Jacobi Iterative Solver and MPI in Fortran**
2. **Parallel Implementation of Poisson Equation using Jacobi Iterative Solver and MPI in C++**
3. **Parallel Implementation of a Random Walk of a Point in a 3D Lattice using MPI (Parallelized in 2D)**

## Table of Contents
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Installation and Compilation](#installation-and-compilation)
  - [Fortran Implementation](#fortran-implementation)
  - [C++ Implementation](#c-implementation)
- [Running the Programs](#running-the-programs)
  - [Fortran Implementation](#running-fortran-implementation)
  - [C++ Implementation](#running-c-implementation)
  - [Random Walk Simulation](#running-random-walk-simulation)
- [Description of the Algorithms](#description-of-the-algorithms)
  - [Jacobi Iterative Solver](#jacobi-iterative-solver)
  - [Random Walk Simulation](#random-walk-simulation)
- [Authors](#authors)

## Project Structure

```plaintext
.
├── fortran_poisson_jacobi/
│   ├── Makefile
│   ├── poisson_jacobi.f90
│   └── README.md
├── cpp_poisson_jacobi/
│   ├── Makefile
│   ├── poisson_jacobi.cpp
│   └── README.md
└── random_walk/
    ├── random_walk_3d_mpi.cpp
    └── README.md
