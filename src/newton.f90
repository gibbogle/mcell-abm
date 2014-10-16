! Newton's method
module newton
use global

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Newtonsolve(v,res)
use :: iso_c_binding
integer :: n, nd, k, kx, ky, kz, i, j, kd, it
real(REAL_KIND) :: v(:), res, F(3), sum
real(REAL_KIND), allocatable :: fv(:), del(:)
integer, parameter :: nit = 3
real(REAL_KIND), parameter :: beta = 1.0	! under-relaxation

nd = Ndim
n = nd*Ncirc*Nlong
allocate(fv(n))
allocate(del(n))
!call nlopt_solve2D(Ncirc,Nlong,v,res)
!call nlopt_solve3D(nd,Ncirc,Nlong,v,res)

do j = 1,Nlong
	do i = 1,Ncirc
		write(*,*)
		call getForce(nd,v,i,j,F)
		write(*,'(a,2i4,3e12.3)') 'i,j,F: ',i,j,F
	enddo
enddo
!stop
do it = 1,nit
	write(*,*) 'makeJacobian'
	call makeJacobian
!	do i = 1,3*Ncirc*Nlong
!		do j = 1,3*Ncirc*Nlong
!			if (abs(Jac(i,j)) > 1) write(*,'(2i4,e12.3)') i,j,Jac(i,j)
!		enddo
!	enddo
	write(*,*) 'invertJacobian'
	call invertJacobian
!	do i = 1,3*Ncirc*Nlong
!		do j = 1,3*Ncirc*Nlong
!			if (abs(Jac(i,j)) > 1.0e4) write(*,'(2i4,e12.3)') i,j,Jac(i,j)
!		enddo
!	enddo
	! Note: Jac(:,:) now holds Jac^-1
	! Set up fv(x)
	do i = 1,Ncirc
		do j = 1,Nlong
			call getForce(nd,v,i,j,F)
			do kd = 1,nd
				k = (j-1)*nd*Ncirc + (i-1)*nd + kd 
				fv(k) = F(kd)
			enddo
		enddo
	enddo
	del = -beta*matmul(Jac,fv)
	sum = 0
	do i = 1,n
		sum = sum + del(i)**2
		v(i) = v(i) + del(i)
	enddo
	sum = sqrt(sum)
	write(*,*) 'del'
	write(*,'(3f8.3)') del
	write(*,'(a,e12.3)') 'norm: ',sum/n
enddo
write(*,*) 'f'
write(*,'(3e12.3)') fv
write(*,*) 'v'
write(*,'(3f8.3)') v
deallocate(fv)
deallocate(del)
stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine invertJacobian
external dlange
real(REAL_KIND) :: dlange
!integer n, nrhs
!parameter (n=10000,lda=2,nrhs=1,ldb=2,nwork = n*1000)
!real(REAL_KIND) :: A(n,n)
!real(REAL_KIND) :: b(n,nrhs)
integer, parameter :: nw = 100
real(REAL_KIND), allocatable :: work(:)
real(REAL_KIND) :: anorm, rcond, eps, det, A(2,2)
integer, allocatable :: ipiv(:)
integer :: n, nwork, info, i, j
logical, parameter :: GET_CONDITION = .false.

n = Ndim*Ncirc*Nlong
allocate(ipiv(n))
nwork = nw*n
allocate(work(nwork))

!A(1,1) = 1
!A(1,2) = 2
!A(2,1) = 3
!A(2,2) = 4
!n = 2
!anorm = DLANGE( '1', n, n, A, n, work )
!call dgetrf(n,n,A,n,ipiv,info)
!call DGECON( '1', n, A, n, anorm, rcond, work, nwork, info )
!write(*,'(a,e12.3)') 'anorm: ',anorm
!write(*,'(a,e12.3)') 'Condition number: ',rcond

!write(*,'(6e12.3)') Jac
eps = 1.0e-10
det = dmgt(eps,n,Jac)
write(*,*)
write(*,'(a,e12.3)') 'determinant of Jac: ',det
anorm = DLANGE( '1', n, n, Jac, n, work )
write(*,'(a,e12.3)') 'anorm: ',anorm
call dgetrf(n,n,Jac,n,ipiv,info)
if (info /= 0) then
	write(*,*) 'dgetrf failed: ',info
	stop
else
	write(*,*) 'dgetrf successful'
endif

if (GET_CONDITION) then
	call DGECON( '1', n, Jac, n, anorm, rcond, work, nwork, info )
	write(*,'(a,e12.3)') 'Condition number: ',rcond
	stop
endif

call dgetri(n,Jac,n,ipiv,work,nwork,info)
deallocate(ipiv)
deallocate(work)
if (info /= 0) then
	write(*,*) 'dgetri failed: ',info
	stop
else
	write(*,*) 'dgetri successful'
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The force F(:) on cell (i,j) depends in general on positions of:
!   (i,j)			
!			kx = (j-1)*nd*Ncirc + (i-1)*nd + 1 
!			ky = (j-1)*nd*Ncirc + (i-1)*nd + 2 
!			kz = (j-1)*nd*Ncirc + (i-1)*nd + 3 
!   (i-1,j)
!			kx = (j-1)*nd*Ncirc + (i-2)*nd + 1 
!			ky = (j-1)*nd*Ncirc + (i-2)*nd + 2 
!			kz = (j-1)*nd*Ncirc + (i-2)*nd + 3 
!   (i+1,j)
!			kx = (j-1)*nd*Ncirc + i*nd + 1 
!			ky = (j-1)*nd*Ncirc + i*nd + 2 
!			kz = (j-1)*nd*Ncirc + i*nd + 3 
!   (i,j-1)
!			kx = (j-2)*nd*Ncirc + (i-1)*nd + 1 
!			ky = (j-2)*nd*Ncirc + (i-1)*nd + 2 
!			kz = (j-2)*nd*Ncirc + (i-1)*nd + 3 
!   (i,j+1)
!			kx = j*nd*Ncirc + (i-1)*nd + 1 
!			ky = j*nd*Ncirc + (i-1)*nd + 2 
!			kz = j*nd*Ncirc + (i-1)*nd + 3 
!-----------------------------------------------------------------------------------------
subroutine makeJacobian
integer :: i, j, istep, in, jn, kx, kxn, nd
integer :: step(0:4,2)
real(REAL_KIND) :: c0(3), cn(3)
real(REAL_KIND), pointer :: pv(:)

step(0,:) = [0,0]	! centre
step(1,:) = [-1,0]	! Left
step(2,:) = [1,0]	! Right
step(3,:) = [0,-1]	! Down
step(4,:) = [0,1]	! Up

nd = Ndim
if (nd == 2) then
	pv => v2D
else
	pv => v3D
endif
Jac = 0
do i = 1,Ncirc
	do j = 1,Nlong
		c0 = mcell(i,j)%centre
		kx = (j-1)*nd*Ncirc + (i-1)*nd + 1 
!			write(*,*) 'i,j: ',i,j
		do istep = 0,4
			in = i + step(istep,1)
			jn = j + step(istep,2)
!				write(*,'(3i4)') istep,in,jn
			if (in >= 1 .and. jn >= 1 .and. in <= Ncirc .and. jn <= Nlong) then
				cn = mcell(in,jn)%centre
				if (istep == 0) then
					kxn = (j-1)*nd*Ncirc + (i-1)*nd + 1
				elseif (istep == 1) then
					kxn = (j-1)*nd*Ncirc + (i-2)*nd + 1 
				elseif (istep == 2) then
					kxn = (j-1)*nd*Ncirc + i*nd + 1 
				elseif (istep == 3) then
					kxn = (j-2)*nd*Ncirc + (i-1)*nd + 1 
				elseif (istep == 4) then
					kxn = j*nd*Ncirc + (i-1)*nd + 1 
				endif
				call derivs(i,j,kx,kxn,c0,cn,pv)
			endif
		enddo
	enddo
enddo
write(*,*) 'Jacobian: ',nd*Ncirc*Nlong,' x ',nd*Ncirc*Nlong
end subroutine

!-----------------------------------------------------------------------------------------
! k0 = index of mcell(i,j)%centre(1)		cell that force derivatives are computed for
! kn0 = index of mcell(in,jn)%centre(1)		neighbour cell
! Don't need c0, cn?
!-----------------------------------------------------------------------------------------
subroutine derivs(i,j,k0,kn0,c0,cn,v)
integer :: i, j, k0, kn0
real(REAL_KIND) :: v(*), c0(3), cn(3)
real(REAL_KIND) :: F0(3), F(3)
integer :: dk, dkn, nd, k, kn
real(REAL_KIND) :: delta = 0.001

nd = Ndim
call getForce(nd,v,i,j,F0)
!write(*,'(a,3f8.3)') 'F0: ',F0
do 	dkn = 0, nd-1
	kn = kn0 + dkn
	v(kn) = v(kn) + delta
	call getForce(nd,v,i,j,F)
	v(kn) = v(kn) - delta
	do dk = 0,nd-1
		k = k0 + dk
		Jac(k,kn) = (F(dk+1) - F0(dk+1))/delta
!		write(*,'(4i4,f10.3)') i,j,k,kn,Jac(k,kn)
	enddo
enddo
end subroutine
	
end module