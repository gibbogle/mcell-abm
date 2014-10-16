! To solve cell motion using RK (r8_rkf45)

module solve
use ellipsoid
use rkf45
use pardiso_mod

implicit none

real(REAL_KIND), allocatable :: F(:,:,:)    ! cell-cell forces
real(REAL_KIND), allocatable :: M(:,:,:)    ! cell-cell moments
!real(REAL_KIND), allocatable :: cp(:,:,:)  ! offset of the contact points from centres
!integer, allocatable :: nbrs(:)             ! number of cell neighbours
!integer, allocatable :: nbrlist(:,:)        ! cell neighbours

integer :: np
logical :: fderiv_ok

    contains

!------------------------------------------------------------------------------------------
! To construct the Jacobian matrix:
!
!         col 1        col 2      col 3      ...  col n
! row 1: df(1)/dx(1) df(1)/dx(2) df(1)/dx(3) ... df(1)/dx(n)
! row 2: df(2)/dx(1) df(2)/dx(2) df(2)/dx(3) ... df(2)/dx(n)
! ...
! row n: df(n)/dx(1) df(n)/dx(2) df(n)/dx(3) ... df(n)/dx(n)
!
! where dx(i)/dt = f(i)
! Jac(i,j) = df(i)/dx(j) = Jac(row,col), 
!------------------------------------------------------------------------------------------
subroutine JacobianSlow(Jac,f,x,n,ok)
integer :: n
real(REAL_KIND) :: x(n), Jac(n,n), f(n)
logical :: ok
real(REAL_KIND), allocatable :: xx(:), ff(:)
real(REAL_KIND) :: t
integer :: i, j, icell, jcell
real(REAL_KIND), parameter :: tol = 1.0e-8
real(REAL_KIND), parameter :: dx = 0.001

allocate(xx(n))
!allocate(f(n))
allocate(ff(n))

Jac = 0
j_deriv = 0
call fderiv(t,x,f)
ok = fderiv_ok
if (.not.ok) return
!write(*,'(10e12.3)') f
xx = x
do j = 1,n
    j_deriv = j
	xx(j) = x(j) + dx
	call fderiv(t,xx,ff)
    ok = fderiv_ok
    if (.not.ok) return
    !if (j == 11) then
    !    write(*,*) 'f and ff for j = ',j
    !    do i = 1,n
    !        write(*,'(i4,3e12.3)') i,f(i),ff(i),(ff(i) - f(i))/dx
    !    enddo
    !endif
	do i = 1,n
    	Jac(i,j) = (ff(i) - f(i))/dx
        if (abs(Jac(i,j)) < tol) then
            Jac(i,j) = 0
        endif
        icell = (i-1)/3 + 1
        jcell = (j-1)/3 + 1
        !if (j == 11) then
        !    write(*,'(a,4i4,3f10.3)') 'Slow Jac: ',i,j,icell,jcell,ff(i),f(i),Jac(i,j)
        !endif
        if (abs(Jac(i,j)) > 2) then
            write(*,'(a,4i4,3f10.4)') 'Big Jac: ',i,j,icell,jcell,ff(i),f(i),Jac(i,j)
        endif
	enddo
	xx(j) = x(j)
!    write(*,*) 'j, Jac(8,11): ',j,Jac(8,11)
enddo
deallocate(xx)
!deallocate(f)
deallocate(ff)
ok = .true.
end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
logical function inlist(i,list,n)
integer :: i, n, list(:)
integer :: k

do k = 1,n
	if (list(k) == i) then
		inlist = .true.
		return
	endif
enddo
inlist = .false.
end function

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine remove(i,list,n)
integer :: i, n, list(:)
integer :: k

do k = 1,n
	if (list(k) == i) then
		list(k) = -i
		return
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
logical function SameCentre(c1,c2)
real(REAL_KIND) :: c1(3), c2(3)
integer :: i
real(REAL_KIND) :: tol = 1.0e-3

SameCentre = .false.
if (abs(c1(1)-c2(1)) > tol) return
if (abs(c1(2)-c2(2)) > tol) return
if (abs(c1(3)-c2(3)) > tol) return
write(*,'(a,6f8.4)') 'Same centres: ',c1,c2
SameCentre = .true.
end function

!------------------------------------------------------------------------------------------
! Note: only those forces (F, M) on cells that are neighbours of the cell associated with
! variable j are changed when xx(j) is varied.  Therefore it is necessary to compute
! only a small subset of the cell-cell forces.  The other F, M are zero.
! df(i)/dx(j) is non-zero (variable i belonging to cell i1) only in the case that cell i1
! owns variable j or has a neighbour cell i2 that owns variable j.  This means that when
! looking at the variation of x(j), we need only to compute interactions for cell i1 and
! neighbours cells of i1 that own variable j. 
! We can precompute the list of cells that need to be considered for variable j.  This 
! could perhaps be done in the j=0 pass, but for clarity maybe not.
!------------------------------------------------------------------------------------------
subroutine JacobianFast(Jac, f0, x, n, ok)
real(REAL_KIND) :: Jac(n,n), f0(n), x(n)
integer :: n
logical :: ok
real(REAL_KIND), allocatable :: xx(:), f1(:), Fbase(:,:,:), Mbase(:,:,:)
real(REAL_KIND) :: t
integer :: i, j, jj
integer :: i1, i2, jn, k1, k2, k, nzero, ntzero, nfeval, kpar=0
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3), a2, b2, centre2(3), orient2(3), s1, s2, delta
real(REAL_KIND) :: FF(3), MM1(3), MM2(3), Fsum(3), Msum(3), R
type(cell_type), pointer :: p1, p2
integer, allocatable :: connected(:,:), nconnected(:)
logical :: hit
logical :: clean_list = .true.
real(REAL_KIND), parameter :: tol = 1.0e-8
real(REAL_KIND), parameter :: dx = 0.001

allocate(xx(n))
!allocate(f0(n))
allocate(f1(n))
allocate(connected(n,MAX_NBRS+1))
allocate(nconnected(n))
allocate(Fbase(ncells,ncells,3))
if (np == 6) then
    allocate(Mbase(ncells,ncells,3))
endif

Jac(1:n,1:n) = 0

! This can be improved by doing it in the j=0 pass, since that will enable dropping of cells that are 
! neighbours but too far away to interact with cell i1.
nconnected = 0
do i1 = 1,ncells
	p1 => cell_list(i1)
	k1 = (i1-1)*np
	do k = 1,np
		j = k1 + k
		if (.not.inlist(i1,connected(j,:),nconnected(j))) then
			nconnected(j) = nconnected(j) + 1
			connected(j,nconnected(j)) = i1
		endif
	enddo
	do jn = 1,p1%nbrs
		i2 = p1%nbrlist(jn)
		k2 = (i2-1)*np
		do k = 1,np
			j = k1 + k
			if (.not.inlist(i2,connected(j,:),nconnected(j))) then
				nconnected(j) = nconnected(j) + 1
				connected(j,nconnected(j)) = i2
			endif
		enddo
	enddo
enddo

!do j = 1,n
!	write(nflog,'(19i4)') j,nconnected(j),connected(j,	1:nconnected(j))
!enddo
write(*,*)
xx = x
nfeval = 0
ntzero = 0
F = 0
if (np == 6) M = 0

do j = 0,n		! column index, except for j=0 which evaluates the current f(:)
	if (j /= 0) then
		xx(j) = x(j) + dx
		F = Fbase
		if (np == 6) M = Mbase
    endif
	nzero = 0
	do i1 = 1,ncells
		p1 => cell_list(i1)
		a1 = p1%a
		b1 = p1%b
		k1 = (i1-1)*np
		centre1 = xx(k1+1:k1+3)
		if (np == 6) then
			orient1 = xx(k1+4:k1+6)
		else
			orient1 = p1%orient
		endif
		!if (np == 6) then
		!	M(i1,1:ncells,:) = 0
		!endif
		do jn = 1,p1%nbrs
			i2 = p1%nbrlist(jn)
            if (i2 == i1) cycle
			if (j > 0) then
				if (.not.inlist(i2,connected(j,:),nconnected(j))) cycle
			endif
			p2 => cell_list(i2)
			a2 = p2%a
			b2 = p2%b
			k2 = (i2-1)*np
			centre2 = xx(k2+1:k2+3)
            if (SameCentre(centre1,centre2)) stop
			if (np == 6) then
				orient2 = xx(k2+4:k2+6)
			else
				orient2 = p2%orient
			endif
			s1 = s1s2(i1,i2,1)
			s2 = s1s2(i1,i2,2)
			call CellInteraction(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,delta,FF,MM1,MM2,ok)
            if (.not.ok) then
                write(*,*) 'istep, i1, i2: ',istep,i1,i2
                return
            endif
			s1s2(i1,i2,1) = s1
			s1s2(i1,i2,2) = s2
			s1s2(i2,i1,1) = s1
			s1s2(i2,i1,2) = s2
			nfeval = nfeval + 1
			if (FF(1)==0 .and. FF(2)==0 .and. FF(3)==0) then
				nzero = nzero + 1
				ntzero = ntzero + 1
				
				if (clean_list) then
					do k = 1,np
						jj = k1 + k
						call remove(i2,connected(jj,:),nconnected(jj))
					enddo
					do k = 1,np
						jj = k2 + k
						call remove(i1,connected(jj,:),nconnected(jj))
					enddo
				endif
	
			endif
			if (.not.ok) then
				write(*,*) 'istep, i1, i2: ',istep,i1,i2
				return
			endif
			F(i1,i2,:) = FF
			F(i2,i1,:) = -F(i1,i2,:)
			if (np == 6) then
		        M(i1,i2,:) = MM1
		        M(i2,i1,:) = MM2
		    endif
		enddo
	enddo
	if (j == 0) then
		Fbase = F
		if (np == 6) Mbase = M
	endif
	
	f1 = 0
	do i1 = 1,ncells
		Fsum = 0
		Msum = 0
		p1 => cell_list(i1)
		do jn = 1,p1%nbrs
			i2 = p1%nbrlist(jn)
			Fsum = Fsum + F(i1,i2,:)
			if (np == 6) then
		        Msum = Msum + M(i1,i2,:)
            endif
		enddo
		k1 = (i1-1)*np
		f1(k1+1:k1+3) = Fsum/Fdrag
		if (np == 6) then
    		orient1 = xx(k1+4:k1+6)
			call cross_product(Msum/Mdrag,orient1,f1(k1+4:k1+6))
		endif
	enddo
	if (j == 0) then
		f0 = f1
	else
		do i = 1,n
			Jac(i,j) = (f1(i) - f0(i))/dx
            if (abs(Jac(i,j)) < tol) then
                Jac(i,j) = 0
            endif
		enddo
		xx(j) = x(j)	
	endif
enddo

deallocate(xx)
!deallocate(f0)
deallocate(f1)
deallocate(connected)
deallocate(nconnected)
deallocate(Fbase)
if (np == 6) deallocate(Mbase)
ok = .true.
end subroutine

!------------------------------------------------------------------------------------------
! Note: Jac is replaced by I - dt.Jac
!------------------------------------------------------------------------------------------
subroutine DoInversion(Jac,Ainv,n,dt,ok)
real(REAL_KIND) :: Jac(n,n), Ainv(n,n), dt
integer :: n
logical :: ok
integer :: i, res

Jac = -dt*Jac
do i = 1,n
    Jac(i,i) = 1 + Jac(i,i)
enddo
call invert(Jac,Ainv,n,res)
ok = (res==0)
end subroutine

!------------------------------------------------------------------------------------------
! x(:) holds the position (centre) and orientation (orient) of each cell.
! The neighbour list nbrlist(:,:) holds each cell's neighbours, assumed not to change
! over the duration of this solving step.
! The force between two cells is stored in F(:,:,:), and the contact points on the two cells
! is in cp(:,:,:).  If the force between cells i1 and i2 has already been computed: F(i1,i2,:),
! then F(i2,i1,:) is set = -F(i1,i2,:) and is not recomputed later when cell i2 is treated.
!
!Note: could take account of growth of a, b with t!
!------------------------------------------------------------------------------------------
subroutine fderiv(t,x,xp)
real(REAL_KIND) :: t, x(*), xp(*)
integer :: i1, i2, j, k1, k2, kpar=0
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3), a2, b2, centre2(3), orient2(3), s1, s2, delta
real(REAL_KIND) :: FF(3), MM1(3), MM2(3), Fsum(3), Msum(3), R, mamp, vm(3), dangle, bigsum1(3), bigsum2(3)
type(cell_type), pointer :: p1, p2
logical :: ok

overlap_average = 0
overlap_max = 0
F = 0
if (np == 6) then
    M = 0
endif
do i1 = 1,ncells
    p1 => cell_list(i1)
    a1 = p1%a
    b1 = p1%b
    k1 = (i1-1)*np
    centre1 = x(k1+1:k1+3)
    if (np == 6) then
		orient1 = x(k1+4:k1+6)
	else
		orient1 = p1%orient
	endif
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        if (i2 == i1) cycle
        if (F(i1,i2,1) /= 0) cycle		!??????????  !  try removing this, repeat computation with i2,i1
        p2 => cell_list(i2)
        a2 = p2%a
        b2 = p2%b
        k2 = (i2-1)*np
        centre2 = x(k2+1:k2+3)
        if (SameCentre(centre1,centre2)) stop
	    if (np == 6) then
			orient2 = x(k2+4:k2+6)
		else
			orient2 = p2%orient
        endif
!		s1 = 0.5	! initial guesses
!		s2 = 0.5
		s1 = s1s2(i1,i2,1)
		s2 = s1s2(i1,i2,2)
        call CellInteraction(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,delta,FF,MM1,MM2,ok)
        if (.not.ok) then
            write(*,*) 'istep, i1, i2: ',istep,i1,i2
            fderiv_ok = .false.
            return
        endif
		s1s2(i1,i2,1) = s1
		s1s2(i1,i2,2) = s2
		s1s2(i2,i1,1) = s1
		s1s2(i2,i1,2) = s2
        F(i1,i2,:) = FF
        F(i2,i1,:) = -F(i1,i2,:) + (/1.0e-12,0.0,0.0/)      !  try removing this, repeat computation with i2,i1
        if (np == 6) then
			M(i1,i2,:) = MM1
			M(i2,i1,:) = MM2
		endif
		if (delta < 0) then
			overlap_average = overlap_average - delta
			overlap_max = max(overlap_max,-delta)
		endif
    enddo
enddo
overlap_average = overlap_average/ncells
bigsum1 = 0
bigsum2 = 0
do i1 = 1,ncells
    p1 => cell_list(i1)
    Fsum = 0
    Msum = 0
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        Fsum = Fsum + F(i1,i2,:)
        if (np == 6) then
	        Msum = Msum + M(i1,i2,:)
        endif
        bigsum1 = bigsum1 + F(i1,i2,:)
        bigsum2 = bigsum2 + F(i2,i1,:)
    enddo
    k1 = (i1-1)*np
    xp(k1+1:k1+3) = Fsum/Fdrag
!!    mamp = sqrt(dot_product(Msum,Msum))
!!    vm = Msum/mamp     ! unit vector of moment axis
!!    dangle = mamp/Mdrag
!!    orient1 = x(k1+4:k1+6)
!!    call rotate(orient1,vm,dangle,xp(k1+4:k1+6))
	if (np == 6) then
		orient1 = x(k1+4:k1+6)
		call cross_product(Msum/Mdrag,orient1,xp(k1+4:k1+6))
    endif
enddo
!write(*,*) 'xp:'
!write(*,'(3f7.3,4x,3f7.3)') xp(1:6*ncells)
if (abs(bigsum1(1)) > 1.0) then
    write(*,'(a,6f10.4)') 'bigsum1,bigsum2: ', bigsum1,bigsum2
    stop
endif
fderiv_ok = .true.
ok = .true.
end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine solver(dt, nt,ok)
real(REAL_KIND) :: dt
integer :: nt
logical :: ok
integer :: kcell, k, i, j, kk, icell, jcell, nvars, flag, kmax, res
type(cell_type), pointer :: p
real(REAL_KIND), allocatable :: x(:)        ! the cell position and orientation variables
real(REAL_KIND), allocatable :: xp(:)       ! derivs of cell position and orientation variables
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND) :: amp, sum, dmax
real(REAL_KIND), allocatable :: Jac(:,:), Ainv(:,:), dx(:), JacF(:,:)
real(REAL_KIND) :: dxlimit = 1

if (simulate_rotation) then
	np = 6
else
	np = 3
endif
nvars = np*ncells
allocate(x(nvars),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: x'
	stop
endif
allocate(xp(nvars),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: xp'
	stop
endif
allocate(F(ncells,ncells,3),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: F'
	stop
endif
allocate(M(ncells,ncells,3),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: M'
	stop
endif
!write(*,*) 'Memory required: ',8*(2*nvars + 2*3*ncells*ncells)
!allocate(cp(ncells,ncells,3))
!allocate(nbrs(ncells))
!allocate(nbrlist(ncells,MAX_NBRS))

do kcell = 1,ncells
    p =>cell_list(kcell)
    k = (kcell-1)*np
    x(k+1:k+3) = p%centre
    if (np == 6) then
		x(k+4:k+6) = p%orient
	endif
enddo

if (isolver == SLOW_EULER_SOLVER .or. isolver == FAST_EULER_SOLVER) then
    allocate(Jac(nvars,nvars),stat=res)
    if (res /= 0) then
	    write(*,*) 'Allocate failed for: Jac'
	    stop
    endif
    allocate(Ainv(nvars,nvars),stat=res)
    if (res /= 0) then
	    write(*,*) 'Allocate failed for: Ainv'
	    stop
    endif
    allocate(dx(nvars),stat=res)
    if (res /= 0) then
	    write(*,*) 'Allocate failed for: dx'
	    stop
    endif
    !allocate(JacF(nvars,nvars),stat=res)
    !if (res /= 0) then
	   ! write(*,*) 'Allocate failed for: JacF'
	   ! stop
    !endif
    if (isolver == FAST_EULER_SOLVER) then
	    call JacobianFast(Jac, xp, x, nvars, ok)    ! xp(:) = f(x)
        if (.not.ok) return
    !	write(*,*) 'got Jac fast:'
    !   write(*,'(9f8.3)') Jac
    else
        call JacobianSlow(Jac, xp, x, nvars,ok)   ! xp(:) = f(x)
        if (.not.ok) return
    endif
    !do i = 1,nvars
    !    icell = (i-1)/np + 1
    !    do j = 1,nvars
    !        jcell = (j-1)/np + 1
    !        if (abs(Jac(i,j) - JacF(i,j)) > 1.0e-6) then
    !            write(*,'(4i6,3e12.3)') i,j,icell,jcell,Jac(i,j),JacF(i,j),Jac(i,j)-JacF(i,j)
    !        endif
    !    enddo
    !enddo  
    
    call DoInversion(Jac,Ainv,nvars,dt,ok)
    do i = 1,nvars
        sum = 0
        do j = 1,nvars
            sum = sum + Ainv(i,j)*xp(j)
        enddo
        dx(i) = dt*sum
        !if (istep == 450) then
        !    write(*,'(a,i4,3f10.4)') 'sum: ',i,dt,sum,dx(i)
        !endif
    enddo
    dmax = 0
    do kcell = 1,ncells
        k = (kcell-1)*np
        sum = 0
        do i = 1,3
            sum = sum + dx(k+i)**2
        enddo
        sum = sqrt(sum)
        if (sum > dmax) then
            dmax = sum
            kmax = kcell
        endif
    enddo
    if (dmax > 1) then
        k = (kmax-1)*np
        do kk = 1,3
            i = k + kk
            write(*,*) 'dmax > 1: ',kmax,i,dmax
            do j = 1,nvars
                write(*,'(i6,3f10.3)') j,Ainv(i,j),xp(j),Ainv(i,j)*xp(j)
            enddo
        enddo
    endif
                  
    do i = 1,nvars
        x(i) = x(i) + dx(i)
    enddo
    if (allocated(JacF)) deallocate(JacF)
    deallocate(Jac)
    deallocate(Ainv)
    deallocate(dx)
    
    if (dmax > run_dmax) then
        run_dmax = dmax
        run_kmax = kmax
        run_maxstep = istep
    endif
    write(*,'(a,i4,f10.3,2x,2i6,f10.3)') 'max displacement: ',kmax,dmax,run_maxstep,run_kmax,run_dmax 
    if (dmax > dxlimit) stop
else
	tstart = 0
	xp = 0
	F = 0
	M = 0
	call fderiv(tstart,x,xp)
    ok = fderiv_ok
    if (.not.ok) return
	abserr = sqrt ( epsilon ( abserr ) )
	relerr = sqrt ( epsilon ( relerr ) )
	flag = 1
	do j = 1,nt
		tstart = (j-1)*dt
		tend = tstart + dt
		call r8_rkf45 ( fderiv, nvars, x, xp, tstart, tend, relerr, abserr, flag )
		if (flag == 4) then
            write(*,*) 'flag: ',flag
			call r8_rkf45 ( fderiv, nvars, x, xp, tstart, tend, relerr, abserr, flag )
		endif
		if (flag /= 2) then
			write(logmsg,*) 'Bad flag: ',flag
			call logger(logmsg)
			deallocate(x)
			deallocate(xp)
			deallocate(F)
			deallocate(M)
			ok = .false.
			return
        endif
        write(*,'(10f7.3)') xp(1:10)
		flag = 2
	enddo
endif
do kcell = 1,ncells
    p => cell_list(kcell)
    k = (kcell-1)*np
    p%centre = x(k+1:k+3)
    if (np == 6) then
		p%orient = x(k+4:k+6)
		amp = sqrt(dot_product(p%orient,p%orient))
		p%orient = p%orient/amp
	endif
enddo
deallocate(x)
deallocate(xp)
deallocate(F)
deallocate(M)
ok = .true.
end subroutine

end module

