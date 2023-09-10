program Mandelbrot
  use mpi
  use MandelbrotModule
  use GlobalVars
  
  implicit none
  integer, allocatable :: grid_c(:,:)
  integer :: ierr, myrank, mysize
  integer :: i, j, i0, j0, i1, j1
  integer :: Nx, Ny, Nx0, Ny0, N_total, temp, start_loc, stop_loc, counter
  integer, allocatable :: stop_idx(:), data_out(:), data_gathered_out(:), sendcounts(:), displs(:)
  real :: hx, hy, hx0, hy0
  real(kind=8) :: t_start, t_stop, t_summed, t_summed_sq
  integer :: job_size
  complex :: c
  integer(kind=8) :: op_counter, sum, sum_div
  
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

  !Fill coarse grid by process 0  
  Nx0 = 300
  Ny0 = 300
  hx0 = (xN - x1)/(Nx0-1)
  hy0 = (yN - y1)/(Ny0-1)

  allocate(grid_c(Nx0, Ny0))
  if (myrank == 0) call scan_grid(Nx0, Ny0, grid_c)

  
  !define fine grid
  Nx = 3000
  Ny = 3000
  hx = (xN - x1)/(Nx-1)
  hy = (yN - y1)/(Ny-1)
  N_total = Nx*Ny
  allocate(stop_idx(mysize))
  
  if (myrank == 0) sum = calc_weight(Nx, Ny, Nx0, Ny0, grid_c)
  sum_div = sum/mysize
   
  if (myrank == 0) then
     allocate(data_gathered_out(N_total))
     !find domains for each process
     sum = 0
     temp = 1 !rank holder
     do j=0,Ny-1
        do i=0,Nx-1
           sum = sum + grid_c(int(i*hx/hx0 + 0.5) + 1, int(j*hy/hy0 + 0.5) + 1)
           if (sum > sum_div) then
              stop_idx(temp) = j*Nx + i + 1
              temp = temp + 1
              sum = 0
           end if
        end do
     end do
  
     stop_idx(mysize) = N_total
  end if

  !Broadcast domain limits 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_BCAST(stop_idx, mysize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   

  !Start parallel execution
  t_start = MPI_WTIME()
  if (myrank == 0) then
     job_size = stop_idx(1)
     start_loc = 1
  else
     job_size = stop_idx(myrank+1) - stop_idx(myrank)
     start_loc = stop_idx(myrank) + 1
  end if

  stop_loc = stop_idx(myrank + 1)

  call index2point(start_loc, Nx, i0, j0)

  !Some tests!
  !call index2point(stop_loc, Nx, i1, j1)
  !print*, "myrank = ", myrank, "start_loc = ", start_loc, "stop_loc = ", stop_loc
  !print*, "i0 = ", i0, "j0 = ", j0, "i1 = ", i1, "j1 = ", j1
  !print*, ""
  !End tests!

  counter = 0
  allocate(data_out(job_size))
  i=i0
  j=j0
  !"Loop" from (i0,j0) to the end point
  op_counter = 0
  do while (counter < job_size)
     c = cmplx(x1 + (i-1)*hx, y1 + (j-1)*hy)
     data_out(counter + 1) = test_point(c)
     op_counter = op_counter + data_out(counter+1)
     if (i+1 > Nx) then
        j = j+1
        i = 1
     else
        i = i+1
     end if
     counter = counter + 1
  end do
  t_stop = MPI_WTIME()
  print*, "myrank = ", myrank, "total number of iterations = ", op_counter


  
  !Merge all arrays
  allocate(sendcounts(mysize)) !send(and receive) counts
  allocate(displs(mysize))  !displacements
  
  call MPI_GATHER(job_size, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
     displs(1) = 0
     do i=2,mysize
        displs(i) = sendcounts(i-1) + displs(i-1)
     end do
     !print*, "displs = ", displs
     !print*, "sendcounts = ", sendcounts
  end if
  
  call MPI_GATHERV(data_out, size(data_out), MPI_INTEGER, data_gathered_out, sendcounts, displs, &
       MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  if (myrank == 0) call output_to_file("out.dat", Nx, Ny, data_gathered_out)



  !Calculate avarage times and deviations
  call MPI_REDUCE(t_stop-t_start, t_summed, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE((t_stop-t_start)**2, t_summed_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
     print*, "Average elapsed time = ", t_summed/mysize, " +- ", sqrt(t_summed_sq/mysize - (t_summed/mysize)**2)
  end if
  
  
  call MPI_FINALIZE(ierr)
  
    
  
end program Mandelbrot

