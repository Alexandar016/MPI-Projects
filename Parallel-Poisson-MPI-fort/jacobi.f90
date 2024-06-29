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

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  ! Set up dimensions and periods for the Cartesian grid
  ndims = 2   ! 2D grid
  dims(1) = num_procs/2 ! Size of dimension 1 (example: 2 processes in dimension 1)
  dims(2) = num_procs/2 ! Size of dimension 2 (example: 2 processes in dimension 2)
  periods(1) = .true. ! Wrap around in dimension 1
  periods(2) = .true. ! Wrap around in dimension 2
  reorder = .false. ! Do not reorder ranks

  ! Create the Cartesian communicator
  call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, cart_comm, ierr)
  ! Get the coordinates of this process in the grid
  call MPI_Cart_coords(cart_comm, my_rank, ndims, coords, ierr)

  ! Print out the coordinates of this process
  print *, 'Process', my_rank, 'is at position (', coords(1), ',', coords(2), ')'

  ! Finalize MPI
  call MPI_FINALIZE(ierr)
  
contains
  subroutine f(x, y, result)
    implicit none
    double precision, intent(in)  :: x, y
    double precision, intent(out) :: result
    result = 2 * ((1 + x) * sin(x + y) - cos(x + y))
  end subroutine f

  subroutine g(x, y, result)
    implicit none
    double precision, intent(in)  :: x, y
    double precision, intent(out) :: result
    result = (1 + x) * sin(x + y)
  end subroutine g
end program poisson_jacobi
