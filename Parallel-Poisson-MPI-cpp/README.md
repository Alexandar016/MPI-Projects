# Parallel-Poisson-MPI
The Jacobi method is an iterative technique for solving linear equations. In the context of solving the Poisson equation, it involves discretizing the domain into a grid of points and updating the value of each grid point based on the values of its neighboring points until the solution converges to a desired accuracy.

## Usage
To compile and run the project you will need a mpic++ (gcc) compiler
```bash
# Example of configuration and installation on wished number of processes "p"
make
make run
```
After running one can also delete compiled ouptuts with 
```bash
make clean
```
