!> @file poisson_jacobi.f90
!! @brief Parallel solution of the Poisson equation using Jacobi method with MPI.
!! 
!! This program solves the 2D Poisson equation on a Cartesian grid of MPI processes.
!! It uses the Jacobi iterative method for solving the equation in parallel.
!! Boundary conditions are not periodic and are implemented using MPI sendrecv operations.
!!
!! The solution approach involves:
!! - Decomposing the domain across MPI processes in a Cartesian grid.
!! - Initializing the grid and source terms.
!! - Iteratively applying the Jacobi method and exchanging boundaries with neighboring processes.
!! - Checking for convergence based on a tolerance value.
!!
!! @author Aleksandar Mitic
!! @June 2024
!!
program poisson_jacobi
  use mpi
  implicit none
  
  integer :: ierr, my_rank, num_procs
  integer :: cart_comm
  integer :: coords(2)
  integer :: ndims
  integer :: dims(2)
  logical :: periods(2)
  logical :: reorder
  integer :: north, south, east, west

  ! Define problem parameters
  integer, parameter :: N = 100
  integer, parameter :: max_iter = 10000
  double precision, parameter :: tol = 1.0e-6
  double precision, parameter :: dx = 1.0d0 / (N - 1)
  double precision :: delta, global_delta

  ! MPI status and request arrays
  integer :: status(MPI_STATUS_SIZE)

  ! Arrays for the solution and source terms
  double precision, allocatable :: u(:,:), u_new(:,:), f_array(:,:)
  double precision :: x, y
  integer :: i, j, iter

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  ! Set up dimensions and periods wraparound for the Cartesian grid
  ndims = 2   ! 2D grid
  dims(1) = int(sqrt(real(num_procs)))
  dims(2) = num_procs / dims(1)
  periods(1) = .true.  
  periods(2) = .true.
  reorder = .false.

  ! Create the Cartesian communicator
  call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, cart_comm, ierr)
  ! Get the coordinates of this process in the grid
  call MPI_Cart_coords(cart_comm, my_rank, ndims, coords, ierr)
  ! Get ranks of neighboring processes
  call MPI_Cart_shift(cart_comm, 0, 1, north, south, ierr)
  call MPI_Cart_shift(cart_comm, 1, 1, west, east, ierr)

  ! Allocate arrays
  allocate(u(0:N+1, 0:N+1))
  allocate(u_new(0:N+1, 0:N+1))
  allocate(f_array(0:N+1, 0:N+1))

  ! Initialize the source term f_array and the initial guess for u
  do i = 1, N
    do j = 1, N
      x = (coords(1) * (N-2) + i - 1) * dx
      y = (coords(2) * (N-2) + j - 1) * dx
      call compute_f(x, y, f_array(i,j))
      u(i,j) = 0.0d0
    end do
  end do

  ! Main iteration loop
  do iter = 1, max_iter
    ! Jacobi update
    do i = 1, N
      do j = 1, N
        u_new(i,j) = 0.25d0 * (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - dx*dx*f_array(i,j))
      end do
    end do

    ! Boundary exchange
    ! Send and receive boundaries with neighboring processes

    ! Exchange with north and south
    call MPI_Sendrecv(u(1,:), N+2, MPI_DOUBLE_PRECISION, north, 1, &
                      u(N+1,:), N+2, MPI_DOUBLE_PRECISION, south, 1, &
                      cart_comm, status, ierr)

    call MPI_Sendrecv(u(N,:), N+2, MPI_DOUBLE_PRECISION, south, 2, &
                      u(0,:), N+2, MPI_DOUBLE_PRECISION, north, 2, &
                      cart_comm, status, ierr)

    ! Exchange with west and east
    call MPI_Sendrecv(u(:,1), N+2, MPI_DOUBLE_PRECISION, west, 3, &
                      u(:,N+1), N+2, MPI_DOUBLE_PRECISION, east, 3, &
                      cart_comm, status, ierr)

    call MPI_Sendrecv(u(:,N), N+2, MPI_DOUBLE_PRECISION, east, 4, &
                      u(:,0), N+2, MPI_DOUBLE_PRECISION, west, 4, &
                      cart_comm, status, ierr)

    ! Compute local maximum difference for convergence check
    delta = 0.0d0
    do i = 1, N
      do j = 1, N
        delta = max(delta, abs(u_new(i,j) - u(i,j)))
      end do
    end do

    ! Global convergence check
    call MPI_Allreduce(delta, global_delta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (global_delta < tol) exit

    ! Update solution
    u = u_new
  end do

  ! Print final solution (example for process 0)
  if (my_rank == 0) then
    print *, 'Converged after ', iter, ' iterations with delta = ', global_delta
  end if

  ! Deallocate arrays
  deallocate(u)
  deallocate(u_new)
  deallocate(f_array)

  ! Finalize MPI
  call MPI_FINALIZE(ierr)
  
contains
  subroutine compute_f(x, y, result)
    implicit none
    double precision, intent(in)  :: x, y
    double precision, intent(out) :: result
    result = 2 * ((1 + x) * sin(x + y) - cos(x + y))
  end subroutine compute_f

  subroutine g(x, y, result)
    implicit none
    double precision, intent(in)  :: x, y
    double precision, intent(out) :: result
    result = (1 + x) * sin(x + y)
  end subroutine g
end program poisson_jacobi