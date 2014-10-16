module csr_mod
use real_kind_mod

implicit none

contains

!------------------------------------------------------------------------------------------------------------
! To generate a CSR3 representation of a matrix A(N,N)
!    Output, real AS(NN), the stored entries of the sparse matrix A.
!    NE is the number of nonzeros 
!    Output, integer ( kind = 4 ) row_index(NN), pointers to specify rows for the stored nonzero entries in AS.
!    Output, integer ( kind = 4 ) col_index(NN), pointers to specify columns for the stored nonzero entries in AS.
!    Return error code: err = 0 if OK, = 1 if number of non-zero entries exceeds NN
! This version is 0-based
!------------------------------------------------------------------------------------------------------------
subroutine to_csr3_0(A,N,NN,AS,row_index,col_index,NE,err)
real(REAL_KIND) :: A(N,N)
integer :: N, NN, err
real(REAL_KIND) :: AS(NN)
integer :: row_index(NN), col_index(NN)
integer :: NE
integer :: row, col, ir, cnt, nrow

cnt = 0
ir = 1
do row = 1,N
	nrow = 0
	do col = 1,N
		if (A(row,col) /= 0) then
			if (nrow == 0) then		! first entry in the row
				row_index(ir) = cnt
				ir = ir+1
			endif
			nrow = nrow + 1
			cnt = cnt + 1
            if (cnt > NN) then
                err = 1
                return
            endif
			col_index(cnt) = col-1
			AS(cnt) = A(row,col)
		endif
	enddo
	if (nrow == 0) then
		if (row == 1) then
			row_index(ir) = 0
			ir = ir+1
		else
			row_index(ir) = cnt
			ir = ir+1
		endif
	endif
	row_index(ir) = cnt
enddo
NE = cnt
err = 0
end subroutine

!------------------------------------------------------------------------------------------------------------
! This version is 1-based
!------------------------------------------------------------------------------------------------------------
subroutine to_csr3_1(A,N,NN,AS,row_index,col_index,NE,err)
real(REAL_KIND) :: A(N,N)
integer :: N, NN, err
real(REAL_KIND) :: AS(NN)
integer :: row_index(NN), col_index(NN)
integer :: NE
integer :: row, col, ir, cnt, nrow

cnt = 0
ir = 1
do row = 1,N
	nrow = 0
	do col = 1,N
		if (A(row,col) /= 0) then
			if (nrow == 0) then		! first entry in the row
				row_index(ir) = cnt + 1
				ir = ir+1
			endif
			nrow = nrow + 1
			cnt = cnt + 1
            if (cnt > NN) then
                err = 1
                return
            endif
			col_index(cnt) = col-1 + 1
			AS(cnt) = A(row,col)
		endif
	enddo
	if (nrow == 0) then
		if (row == 1) then
			row_index(ir) = 0 + 1
			ir = ir+1
		else
			row_index(ir) = cnt + 1
			ir = ir+1
		endif
	endif
	row_index(ir) = cnt + 1
enddo
NE = cnt
err = 0
end subroutine

end module

!program main
!use csr_mod
!integer, parameter :: N = 5, NN = 20
!real(REAL_KIND) :: A(N,N), AS(NN)
!integer :: row_index(NN), col_index(NN), NE, err
!
!A = reshape([ &
!      1, -1,  0, -3,  0,  &
!     -2,  5,  0,  0,  0,  &
!	  0,  0,  4,  6,  4,  &
!	 -3,  0,  6,  7,  0,  &
!      0,  0,  0,  0,  0   &
!!	  0,  8,  0,  0, -5   &
!    ],shape(A))
!A = transpose(A)
!
!call to_csr3(A,N,NN,AS,row_index,col_index,NE,err)
!
!write(*,*) 'NE: ',NE
!write(*,'(10f6.1)') AS(1:NE)
!write(*,'(10i6)') col_index(1:NE)
!write(*,'(10i6)') row_index(1:NE)
!
!end