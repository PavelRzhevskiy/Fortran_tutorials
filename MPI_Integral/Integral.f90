program integral
  use mpi
  implicit none

  
  integer :: ierr, myrank, mysize
  integer :: q, r, x_size, N
  integer :: i,j,k, start_index
  real(kind=8) :: h, sum, global_sum, x0, x, y
  real(kind=8) :: t_start, t_stop, t_sum, t_sum_sq

  
  call MPI_INIT(ierr)

  ! define # of pts and grid step 
  N = 2000
  h=1.D0/(N-1)
    
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

       
  ! Get the size of x-sliced chunk for the process
  q = N/mysize
  r = N - q*mysize
  if (myrank + 1 <= r) then
     x_size = q + 1
     start_index = (q + 1)*myrank
  else
     x_size = q
     start_index = q*myrank + r
  end if
  
  !Calculate local sum
  sum = 0
  x0 = start_index*h
  t_start = MPI_WTIME()
  !print*, "rank = ", myrank, "start_index = ", start_index, "arr size = ", arr_size
  do i = 0, x_size-1
     x = x0 + i*h

     do j = 0, N-1
        y = j*h

        do k = 0, N-1
           sum = sum + dsin(x + y + k*h)*dexp(-x - y - k*h)*w(i+1, j+1, k+1)
           
        end do
        
     end do
     
  end do
  t_stop = MPI_WTIME()
  sum = sum*h**3
  print*, "myrank = ", myrank,  "sum =", sum
  
  
  
  call MPI_REDUCE(sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  call MPI_REDUCE(t_stop - t_start, t_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  call MPI_REDUCE((t_stop - t_start)**2, t_sum_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
  
  if (myrank == 0) then
     print*, "Integral = ", global_sum, "error ~ ", h**2
     print*, "Average elapsed time = ", t_sum/mysize, " +- ", sqrt(t_sum_sq/mysize - (t_sum/mysize)**2)
  end if
  
  call MPI_FINALIZE(ierr)

  
contains
  function w(i,j,k)
    integer, intent(in) :: i,j,k
    real(kind=8) :: w
    logical :: x_plane, y_plane, z_plane, edge

    y_plane = (j == 1 .or. j == N)
    z_plane = (k == 1 .or. k == N)
    edge = (myrank == 0 .or. myrank == mysize - 1)

    
    if (.not. edge) then
       ! if not rightmost or leftmost chunk

       !conditions for surfaces and edge lines
       if (y_plane .and. z_plane) then
          w = 1.D0/4
       else if (y_plane .or. z_plane) then
          w = 1.D0/2
       else
          w = 1.D0
       end if
       
    else
       ! if rightmost or leftmost chunk
       if (myrank == 0) then 
          x_plane = (i == 1)
       else
          x_plane = (i == x_size)
       end if

       ! if only one process is involved
       if (mysize == 1) x_plane = ((i == 1) .or. (i == N))

       ! conditions for surfaces, edge lines and corners
       if (x_plane .and. y_plane .and. z_plane) then
          w = 1.D0/8
       else if ((x_plane .and. y_plane) .or. (x_plane .and. z_plane) .or. (y_plane .and. z_plane)) then 
          w = 1.D0/4
       else if (x_plane .or. y_plane .or. z_plane) then
          w = 1.D0/2
       else
          w = 1.D0
       end if


    end if
    
  end function w
  
end program integral
