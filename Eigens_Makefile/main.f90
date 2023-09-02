program eigen
  ! Reads matrix from file and solves eigenvalue problem !
  use workmodule
  implicit none
  
  character(len=20) :: filename
  real(kind=8), allocatable :: m_arr(:,:), wr(:), wi(:)
  
  filename = 'input.dat'
  
  call read_matrix_from_file(filename, m_arr)
  call find_eigenval(m_arr, wr, wi)
  call print_eigenval(wr, wi)
  
end program eigen

