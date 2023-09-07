module IntegralModule
  implicit none

contains

  function w(i,j,k,N)
    integer(kind=8), intent(in) :: i,j,k,N
    real(kind=8) :: w
    logical :: x_plane, y_plane,  z_plane

    x_plane = (i == 1 .or. i == N)
    y_plane = (j == 1 .or. j == N)
    z_plane = (k == 1 .or. k == N)

    if (x_plane .and. y_plane .and. z_plane) then
       w = 1.D0/8
    else if ((x_plane .and. y_plane) .or. (x_plane .and. z_plane) .or. (y_plane .and. z_plane)) then 
       w = 1.D0/4
    else if (x_plane .or. y_plane .or. z_plane) then
       w = 1.D0/2
    else
       w = 1.D0
    end if
 
  end function w
  
  subroutine divide_points(q, r, myrank, start_idx, stop_idx)
    integer, intent(in) :: myrank
    integer(kind=8), intent(in) :: q,r
    integer(kind=8), intent(out) :: start_idx, stop_idx
       
    if (myrank + 1 <= r) then
       start_idx = (q + 1)*myrank + 1
       stop_idx = start_idx + (q+1) -1
    else
       start_idx = q*myrank + r + 1
       stop_idx = start_idx + q - 1
    end if
     
  end subroutine divide_points


end module IntegralModule
