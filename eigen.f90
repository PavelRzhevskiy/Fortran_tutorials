program eigen
  ! Reads matrix from file and solves eigenvalue problem !

  implicit none

  character(len=20) :: filename
  real(kind=8), allocatable :: m_arr(:,:), vr(:,:), wr(:), wi(:), work_arr(:)
  real(kind=8) :: empty(1,1), work_size(1)
  real(kind=8) :: x
  integer :: i, stat, dim, lda, ldvr, lwork, info, opt_size
  
  filename = 'input.dat'
  open(unit = 1, file = filename, status = 'old', action = 'read')
  
  !Read # of rows from file to find dimension, rewind file and read matrix
  stat=0
  dim=-1
  do while(stat == 0)
     dim = dim + 1
     read(1, *, iostat=stat) x
  enddo
  print '("Matrix dimension is",I3)', dim
  rewind(1)
  allocate(m_arr(dim, dim), vr(dim, dim), wr(dim), wi(dim))
  read(1,*) m_arr
  m_arr = transpose(m_arr)
  print*, "Matrix: "
  do i=1,dim
     print *, m_arr(i,:)
  enddo
  
  !Find optimal WORK size
  ldvr = dim
  lda = dim
  lwork = -1
  call dgeev('N', 'V', dim, m_arr, lda, wr, wi, empty, 1, vr, ldvr, work_size, lwork, info)
  opt_size = int(work_size(1))
  print *, "Optimal WORK size = ", opt_size
  allocate(work_arr(opt_size))
  
  !Find eigenvalues (and right eigenvactors)
  call dgeev('N', 'V', dim, m_arr, lda, wr, wi, empty, 1, vr, ldvr, work_arr, opt_size, info)
  print*, "dgeev status: ", info
  print*, "Eigenvalues:"
  do i=1,dim
     print ' ( "(", es10.3, "," e10.3, ")" ) ', wr(i), wi(i)
  enddo
  


end program eigen

