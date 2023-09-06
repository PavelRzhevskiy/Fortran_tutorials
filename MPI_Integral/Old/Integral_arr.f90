program integral
  use mpi

  implicit none
  integer :: ierr, myrank, mysize
  integer :: quotient, remainder, arr_size, N
  integer :: i,j,k, start_index
  real(kind=8) :: h, sum, x0, x, y
  real(kind=8), allocatable :: w(:,:,:), sum_arr(:)
  call MPI_INIT(ierr)

  ! define # of pts and grid step 
  N = 300
  h=1.D0/(N-1)
  print*, h

  
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

  quotient = N/mysize
  remainder = N - quotient*mysize
  if (myrank + 1 <= remainder) then
     arr_size = quotient + 1
     start_index = (quotient + 1)*myrank
  else
     arr_size = quotient
     start_index = quotient*myrank + remainder
  end if
  
  print*, "rank = ", myrank, "size = ", mysize, "arr_size = ", arr_size


  allocate(w(arr_size, N, N))

  ! initialize local array for weights at inner points, edge surfaces, edge bounds, and (possibly) corners and x boundaries
  w = 1.D0 !inner points
  w(:, :, [1,N]) = 1.D0/2 !xoy surfaces
  w(:, [1,N], :) = 1.D0/2 !xoz surfaces
  w(:, [1,N], [1,N]) = 1.D0/4 !x edges
  !leftmost chunk
  if (myrank == 0) then
     w(1, :, :) = 1.D0/2 !yoz surface
     w(1, [1,N], :) = 1.D0/4 !z-edges
     w(1, :, [1,N]) = 1.D0/4 !y-edges
     w(1, [1,N], [1,N]) = 1.D0/8 !corners
  end if
  !rightmost chunk
  if (myrank == mysize - 1) then
     w(arr_size, :, :) = 1.D0/2 !yoz surface
     w(arr_size, [1,N], :) = 1.D0/4 !z-edges
     w(arr_size, :, [1,N]) = 1.D0/4 !y-edges
     w(arr_size, [1,N], [1,N]) = 1.D0/8 !corners
  end if
  
!  print*, w(1,3,5)
 ! print*, w(N,3,3)
  !print*, w(3,4,2)
  



  !Calculate local sum
  sum = 0
  x0 = start_index*h
  do i = 0, arr_size-1
     x = x0 + i*h
     if (myrank == 2) then
        !print*, x
     end if
     do j = 0, N-1
        y = j*h
        
        do k = 0, N-1
           sum = sum + dsin(x + y + k*h)*w(i+1, j+1, k+1)
           
        end do
        
     end do
     
  end do

  deallocate(w)
  sum = sum*h*h*h

  print*, "myrank = ", myrank, "sum = ", sum

  
  call MPI_FINALIZE(ierr)

  
end program integral
