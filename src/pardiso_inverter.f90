module pardiso_mod
use real_kind_mod
use csr_mod
use MKL_PARDISO
implicit none

contains

subroutine invert(A,Ainv,N,res)
real(REAL_KIND) :: A(N,N), Ainv(N,N)
integer :: N, res
TYPE(MKL_PARDISO_HANDLE), allocatable :: PT(:)
integer, allocatable :: perm(:), ivect(:), jvect(:)
real(REAL_KIND), allocatable :: b(:), x(:), AS(:)
INTEGER :: maxfct, mnum, mtype, phase, nrhs, msglvl, i, j, k
INTEGER :: iparm(64)
integer :: NN, NE, err
maxfct = 1
mnum = 1
nrhs = N    ! matrix inversion
! need to choose NN
NN = N*N

allocate( PT(64), perm(N), b(nrhs*N), x(nrhs*N), AS(NN), ivect(NN), jvect(NN) )
! ivect(:) is row_index(:)
! jvect(:) is col_index(:)
  
do i = 1, 64
    iparm(i) = 0
    PT(i)%DUMMY = 0
end do
iparm(1) = 0 ! solver default

err = 0    ! initialize error flag
msglvl = 0 ! 1 = print statistical information
mtype = 11 ! real unsymmetric 
phase = 13 ! solve

call to_csr3_1(A,N,NN,AS,ivect,jvect,NE,err)  ! 1-based
if (err /= 0) then
    deallocate(PT, perm, b, x, AS, ivect, jvect)
    res = 1
    return
endif

!write(*,'(10f6.1)') AS(1:NE)
!write(*,'(10i6)') ivect(1:N+1)
!write(*,'(10i6)') jvect(1:NN)

b = 0
do i = 1,N
    b(i+(i-1)*N) = 1
enddo
err = 0
CALL pardiso (PT, maxfct, mnum, mtype, phase, N, AS, ivect, jvect, perm, nrhs, iparm, msglvl, b, x, err)
if (err /= 0) then
    deallocate(PT, perm, b, x, AS, ivect, jvect)
    res = 2
    return
endif

k = 0
do j = 1,N
    do i = 1,N
        k = k+1
        Ainv(i,j) = x(k)
    enddo
enddo
phase = -1
CALL pardiso (PT, maxfct, mnum, mtype, phase, N, AS, ivect, jvect, perm, nrhs, iparm, msglvl, b, x, err)
deallocate(PT, perm, b, x, AS, ivect, jvect)
res = 0
end subroutine

end module

!PROGRAM PardisoTest
!use csr_mod
!use pardiso_mod
!IMPLICIT NONE
!  
!integer, parameter :: N = 5
!real(REAL_KIND) :: A(N,N), Ainv(N,N)
!integer :: error
!
!A = reshape([ &
!      1, -1,  0, -3,  0,  &
!     -2,  5,  0,  0,  0,  &
!	  0,  0,  4,  6,  4,  &
!	 -3,  0,  6,  7,  0,  &
!	  0,  8,  0,  0, -5   &
!    ],shape(A))
!A = transpose(A)
!
!call invert(A,N,Ainv,error)
!write(*,*)
!if (error /= 0) then
!    write(*,*) 'error: ',error
!else
!    write(*,'(5f10.4)') Ainv
!endif
!
!END PROGRAM PardisoTest
