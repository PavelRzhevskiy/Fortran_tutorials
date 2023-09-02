  ! module containing functions to read matrix from file, find eigenvalues and print eigenvalues
module workmodule
  implicit none

contains


  subroutine read_matrix_from_file(filename, m_arr)
    character(len=20), intent(in) :: filename
    real(kind=8), allocatable, intent(out) :: m_arr(:,:)
    
    integer :: stat, i, dim
    real :: x
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
    allocate(m_arr(dim, dim))
    read(1,*) m_arr
    m_arr = transpose(m_arr)
    print*, "Matrix: "
    do i=1,dim
       print *, m_arr(i,:)
    enddo

  end subroutine read_matrix_from_file


  subroutine find_eigenval(m_arr, wr, wi)
    real(kind=8), allocatable, intent(in) :: m_arr(:,:)
    real(kind=8), allocatable, intent(out) :: wr(:), wi(:)
    
    real(kind=8), allocatable :: vr(:,:)
    real(kind=8), allocatable :: work_arr(:)
    real(kind=8) :: empty(1,1), work_size(1)
    integer :: dim, lda, ldvr, lwork, info, opt_size
    integer :: arr_shape(2)

    arr_shape = shape(m_arr)
    dim = arr_shape(1)
    allocate(vr(dim, dim), wr(dim), wi(dim))
  
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
    
  end subroutine find_eigenval


  subroutine print_eigenval(wr, wi)
    
    real(kind=8), allocatable, intent(in) :: wr(:), wi(:)
    integer :: i, dim
    
    dim = size(wr)
    print*, "Eigenvalues:"
    do i=1,dim
       print ' ( "(", es10.3, "," e10.3, ")" ) ', wr(i), wi(i)
    enddo

  end subroutine print_eigenval

end module workmodule

