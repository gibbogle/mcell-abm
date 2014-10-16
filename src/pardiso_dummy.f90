SUBROUTINE PARDISO_D( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA, PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
USE MKL_PARDISO_PRIVATE
use real_kind_mod
TYPE(MKL_PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
integer :: MAXFCT, MNUM, MTYPE, PHASE, N, IA(*), JA(*), PERM(*), NRHS, IPARM(*), MSGLVL, ERROR
real(REAL_KIND) :: A(*), B(*), X(*)

write(*,*) 'No pardiso!'
end subroutine

