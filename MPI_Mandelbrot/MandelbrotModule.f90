module MandelbrotModule
  use GlobalVars
  implicit none

contains

  subroutine scan_grid(Nx, Ny, grid)
    integer, intent(inout) :: grid(:,:)
    integer, intent(in) :: Nx, Ny
    real :: x, hx, hy
    integer :: i,j
    complex :: c

    hx = (xN - x1)/(Nx - 1)
    hy = (yN - y1)/(Ny - 1)
    do i=0,Nx-1
        x = x1 + i*hx
        do j=0,Ny-1
           c = cmplx(x, y1 + j*hy)
           grid(i+1,j+1) = test_point(c)
        end do
     end do
   end subroutine scan_grid


   function calc_weight(Nx, Ny, Nx0, Ny0, grid_c) result(sum)
     integer, intent(in) :: Nx, Ny, Nx0, Ny0
     integer, intent(in) :: grid_c(:,:)
     real :: hx, hy, hx0, hy0
     integer(kind=8) :: sum
     integer :: i,j

     hx0 = (xN - x1)/(Nx0-1)
     hy0 = (yN - y1)/(Ny0-1)
     hx = (xN - x1)/(Nx-1)
     hy = (yN - y1)/(Ny-1)

     sum=0
     do i=0,Nx-1
        do j=0,Ny-1
           sum = sum + grid_c(int(i*hx/hx0 + 0.5) + 1, int(j*hy/hy0 + 0.5) + 1)
        end do
     end do
   end function calc_weight
     

     
    function test_point(c) result(niter)
    complex, intent(in) :: c
    complex :: z
    integer :: i,niter

    z = c
    do i=1,N_MAX_ITER
       z = z*z + c
       if (abs(z) > REF_VAL) then
          niter = i
          return
        end if
    end do

    niter = N_MAX_ITER
  end function test_point
 

  subroutine index2point(idx, Nx, i0, j0)

    integer, intent(in) :: idx, Nx
    integer :: i0,j0
    j0 = (idx-1)/Nx + 1 !y-coordinate
    i0 = idx - (j0-1)*Nx !x-coordinate

  end subroutine index2point


  subroutine output_to_file(filename, Nx, Ny, data)
    integer(kind=4), intent(in) :: data(:)
    character(len=20), intent(in) :: filename
    integer, intent(in) :: Nx, Ny
    real :: hx, hy
    real :: x,y
    integer(kind=4) :: i,j,k

    hx = (xN - x1)/(Nx-1)
    hy = (yN - y1)/(Ny-1)

    open(unit = 1, file = filename, status = 'replace', action = 'write')

    do k=1,size(data)
       call index2point(k,Nx,i,j)
       x = x1 + (i-1)*hx
       y = y1 + (j-1)*hy
       write(1,*) x, y, data(k)
    end do

  end subroutine output_to_file

  
end module MandelbrotModule
  
