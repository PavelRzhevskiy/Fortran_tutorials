program integral
  use IntegralModule
  use mpi
  implicit none

  
  integer :: ierr, myrank, mysize
  integer(kind=8) :: N1dim, q, r, Ntotal, start_idx, stop_idx, pos, i, j, k
  real(kind=8) :: h, sum, global_sum
  real(kind=8) :: t_start, t_stop, t_summed, t_summed_sq


  call MPI_INIT(ierr)

  !define # of pts in each dimension and grid step 
  N1dim = 2000
  h=1.D0/(N1dim-1)
    
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr) 

  !divide points and calculate start and stop indices for a process
  Ntotal = N1dim**3
  q = Ntotal/mysize
  r = Ntotal - q*mysize

  call divide_points(q, r, myrank, start_idx, stop_idx)
  !print*, "myrank = ", myrank, "start_idx = ", start_idx, "stop_idx = ", stop_idx 
  
  sum = 0
  t_start = MPI_WTIME()
  do pos = start_idx, stop_idx
     !find i,j,k:
     k = (pos-1)/(N1dim**2) + 1
     j = (pos - (k-1)*N1dim**2 - 1)/(N1dim) + 1
     i = (pos - (k-1)*N1dim**2 - (j-1)*N1dim)

     sum = sum + w(i,j,k,N1dim)*dsin((i+j+k-3)*h)*dexp(-(i+j+k-3)*h)
  end do
  t_stop = MPI_WTIME()
 
  sum = sum*h**3
  print*, "myrank = ", myrank, "sum = ", sum

  !pass sum, time, and time^2 to process 0
  call MPI_REDUCE(sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(t_stop-t_start, t_summed, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE((t_stop-t_start)**2, t_summed_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  if (myrank == 0) then
     t_stop = MPI_WTIME()
     print*, "Integral = ", global_sum, "error ~ ", h**2
     print*, "Average elapsed time = ", t_summed/mysize, " +- ", sqrt(t_summed_sq/mysize - (t_summed/mysize)**2)
  end if
    
  
  call MPI_FINALIZE(ierr)

end program integral
