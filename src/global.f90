module global
use real_kind_mod
use omp_lib
use par_zig_mod
use winsock
use, intrinsic :: iso_c_binding

implicit none

!#include 'min_dist_interface.f90'

interface
	subroutine nlopt_test(res) bind(C,name="nlopt_test") 
	use :: iso_c_binding
	real(c_double) :: res
	end subroutine
end interface

interface
	subroutine nlopt_solve2D(nc,nl,v,res) bind(C,name="nlopt_solve2D") 
	use :: iso_c_binding
	integer(c_int), value :: nc, nl
	real(c_double) :: v(*)
	real(c_double) :: res
	end subroutine
end interface

interface
	subroutine nlopt_solve3D(nd,nc,nl,v,res) bind(C,name="nlopt_solve3D") 
	use :: iso_c_binding
	integer(c_int), value :: nd, nc, nl
	real(c_double) :: v(*)
	real(c_double) :: res
	end subroutine
end interface

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)

real(REAL_KIND), parameter :: PI = 4*atan(1.0d0)

type cell_type
	real(REAL_KIND) :: a, b         ! current axes (um)
	real(REAL_KIND) :: aspect_n     ! "normal" aspect ratio (at g = 0.25)
	real(REAL_KIND) :: a_n, b_n     ! "normal" axes
	real(REAL_KIND) :: V_n          ! "normal" volume
	real(REAL_KIND) :: beta         ! b-constancy parameter (0 - 1)
	real(REAL_KIND) :: centre(3)    ! ellipsoid centre position
	real(REAL_KIND) :: orient(3)    ! ellipsoid main axis unit vector
	real(REAL_KIND) :: F(3)         ! current total force
	real(REAL_KIND) :: M(3)         ! current total moment
	real(REAL_KIND) :: Fprev(3)     ! previous total force
	real(REAL_KIND) :: Mprev(3)     ! previous total moment
	real(REAL_KIND) :: birthtime
	real(REAL_KIND) :: cycletime
	real(REAL_KIND) :: growthrate
	integer :: nbrs
	integer, allocatable :: nbrlist(:)
end type

type mcell_type
	real(REAL_KIND) :: centre(3)
	real(REAL_KIND) :: volume
	real(REAL_KIND) :: vtarget
	real(REAL_KIND) :: width(3)
	real(REAL_KIND) :: F(3)         ! current total force
	real(REAL_KIND) :: vert(8,3)	! vertices
end type

type XYZ_type
    real(REAL_KIND) :: x, y, z
end type

type, bind(C) :: point_type
    real(REAL_KIND) :: x(3)
end type

type Fparam_type
    real(REAL_KIND) :: a, b, c, e, g
end type

integer, parameter :: MAX_CHEMO = 4		! just a dummy for now
type, bind(C) :: field_data
	integer(c_int) :: site(3)
	integer(c_int) :: state
	real(c_double) :: dVdt
	real(c_double) :: volume
	real(c_double) :: conc(MAX_CHEMO+1)	! This needs to agree with the definition in field.h
end type

type, bind(C) :: hexahedron
	type(point_type) :: vertex(8)
end type

integer, parameter :: nflog=10, nfin=11, nfout=12, nfres=13
integer, parameter :: MAX_CELLS = 1000
integer, parameter :: MAX_NBRS = 30
real(REAL_KIND), parameter :: CYCLETIME0 = 12*60	! 12 hours -> minutes
logical, parameter :: POLARITY = .false.
integer, parameter :: EXPLICIT_SOLVER = 0
integer, parameter :: SLOW_EULER_SOLVER = 1
integer, parameter :: FAST_EULER_SOLVER = 2
integer, parameter :: RKF45_SOLVER = 3

integer, parameter :: Ndim = 3

character*(128) :: inputfile
character*(128) :: outputfile

integer :: Ncirc, Nlong, isolver
real(REAL_KIND) :: dx_init

integer :: NX, NY, NZ
integer :: Mnodes, ncpu_input, ncells, nsteps, istep
integer :: seed(2)
real(REAL_KIND) :: DELTA_T
type(Fparam_type) :: FP1, FP2
TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
!logical :: simulate_rotation
!logical :: simulate_growth
real(REAL_KIND) :: pressure, tension, Falpha
real(REAL_KIND) :: Rinitial

character*(128) :: logfile
character*(2048) :: logmsg
character*(128) :: cellmlfile

type(mcell_type), allocatable, target :: mcell(:,:)
real(REAL_KIND), allocatable, target :: v2D(:), v3D(:), v0(:,:), v1(:,:)
real(REAL_KIND), allocatable :: Jac(:,:)
real(REAL_KIND),allocatable :: cellml_state0(:)

type(cell_type), allocatable, target :: cell_list(:)
type(cell_type), allocatable, target :: cell_list0(:)
integer, allocatable :: s1s2(:,:,:)

real(REAL_KIND) :: run_dmax, overlap_average, overlap_max
integer :: run_kmax, run_maxstep
integer :: j_deriv

logical :: dbug = .false.

!dec$ attributes dllexport :: nsteps, DELTA_T

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
if (R == -2147483648) R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!----------------------------------------------------------------------------------------
! R = A x B 
!----------------------------------------------------------------------------------------
subroutine cross_product(A, B, R)
real(REAL_KIND) :: A(3), B(3), R(3)

R(1) = A(2)*B(3) - A(3)*B(2)
R(2) = A(3)*B(1) - A(1)*B(3)
R(3) = A(1)*B(2) - A(2)*B(1)
end subroutine

!----------------------------------------------------------------------------
!The subroutine TSRGT applies to input real square matrix A(n,n) the upper
!triangularization algorithm of Gauss method with full pivoting and keeps
!trace of successive transformations done in integer vectors KP and LP.
!-----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix (real*8)
!  Output parameters:
!  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
!  C          pointer to table storing main diagonal elements and supra-
!             diagonal elements of upper triangular matrix and the multi-
!             plying coefficients used during triangularization process
!  KP         table storing informations concerning the column exchanges
!             during process (integer)
!  LP         table storing informations concerning the line exchanges
!             during process (integer)
!-----------------------------------------------------------------------------
!The table C is first initialized to A matrix, then receives at each step k
!of the triangularization process, usefull elements of A matrix at step k for
!k=1,2,...n.
!The variables po(real*8), lo and ko(integer) store respectively pivot at step k,
!its line number and its column number.
!------------------------------------------------------------------------------
Subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
  real(REAL_KIND) :: eps
  integer n,it
  real(REAL_KIND) :: A(n,n), C(n,n)
  integer Kp(n),Lp(n) 
  real(REAL_KIND) ::  po,t0
  integer :: k, i, j, lo, ko
  
  C=A
  it=1
  k=1
  do while (it==1.and.k<n)
    po=C(k,k); lo=k; ko=k
    do i=k, n
      do j=k, n
        if (dabs(C(i,j))>dabs(po)) then
          po=C(i,j); lo=i; ko=j
        end if
      end do
    end do
    Lp(k)=lo; Kp(k)=ko
    if (dabs(po)<eps) then
      it=0
    else
      if (lo.ne.k) then
        do j=k, n
          t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
        end do
      end if
      if (ko.ne.k) then
        do i=1, n
          t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
        end do
      end if 
      do i=k+1, n
        C(i,k)=C(i,k)/po
        do j=k+1, n
          C(i,j)=C(i,j)-C(i,k)*C(k,j)
        end do 
      end do
      k=k+1
    end if
  end do
  if (it==1.and.dabs(C(n,n))<eps)  it=0
End subroutine

!----------------------------------------------------------------------------
!The function DMGT returns the determinant of a real square matrix
!A(n,n) by Gauss method with full pivoting.
!----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix
!  Output parameters:
!  None
!-----------------------------------------------------------------------------
!The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
!Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
!If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
!at location i,j (j>=i) the corresponding element of the upper triangular matrix.
!Tables Lp and Kp contain informations relative to exchanges of line or column
!that occured during the process. For instance, the element number k of Lp is
!an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
!The number of exchanges of lines and columns is stored in integer L. the
!determinant of A matrix is stored in d0 (real*8).
!-----------------------------------------------------------------------------
real(REAL_KIND) Function DMGT(eps, n, A)
implicit none
integer :: n
real(REAL_KIND) :: eps, A(n,n)
real(REAL_KIND) :: d0

real(REAL_KIND), pointer :: C(:,:)
integer,pointer :: Kp(:), Lp(:)
integer :: k, it, l, ialloc

!allocate local matrix C and vectors Kp, Lp
  allocate(C(n,n),STAT=ialloc)
  allocate(Kp(n),STAT=ialloc)
  allocate(Lp(n),STAT=ialloc)

  call TSRGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine
  if (it==0) then
    d0=0.d0  !matrix singular, det=0
  else       !matrix regular, det<>0
    d0=1.d0
    do k=1, n
	  d0=d0*C(k,k)
    end do
    l=0
    do k=1, n-1
      if (Lp(k).ne.k)  l=l+1
      if (Kp(k).ne.k)  l=l+1
    end do
    if (MOD(l,2).ne.0) d0=-d0  !l is odd
  end if
  DMGT=d0   !return determinant
  return
End function

!-----------------------------------------------------------------------------------------------------
! Rotate v0 about unit vector vN through angle to get v
!-----------------------------------------------------------------------------------------------------
subroutine Rotate(v0, vN, angle, v)
real(REAL_KIND) :: v0(3), vN(3), angle, v(3)
real(REAL_KIND) :: cosa, sina, d
type(XYZ_type) :: q1, q2, u

cosa = cos(angle)
sina = sin(angle)

q1%x = v0(1)
q1%y = v0(2)
q1%z = v0(3)
u%x = vN(1)
u%y = vN(2)
u%z = vN(3)

d = sqrt(u%y*u%y + u%z*u%z)

! Step 2 
if (d /= 0) then
    q2%x = q1%x
    q2%y = q1%y * u%z / d - q1%z * u%y / d
    q2%z = q1%y * u%y / d + q1%z * u%z / d
else
  q2 = q1
endif

! Step 3
q1%x = q2%x * d - q2%z * u%x
q1%y = q2%y
q1%z = q2%x * u%x + q2%z * d

! Step 4
q2%x = q1%x * cosa - q1%y * sina
q2%y = q1%x * sina + q1%y * cosa
q2%z = q1%z

! Inverse of step 3
q1%x =   q2%x * d + q2%z * u%x
q1%y =   q2%y
q1%z = - q2%x * u%x + q2%z * d

! Inverse of step 2 
if (d /= 0) then
    v(1) =   q1%x
    v(2) =   q1%y * u%z / d + q1%z * u%y / d
    v(3) = - q1%y * u%y / d + q1%z * u%z / d
else
    v(1) = q1%x
    v(2) = q1%y
    v(3) = q1%z
endif
end subroutine

!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
subroutine test_rotate
real(REAL_KIND) :: v0(3), vN(3), angle, v(3), r(3), F(3), M(3), mamp

v0(1) = 1   ! cell axis
v0(2) = 0
v0(3) = 0
F(1) = 0    ! force
F(2) = 0
F(3) = -1
r(1) = 3  ! offset of point of contact from centre
r(2) = 0
r(3) = 0
call cross_product(r,F,M)
mamp = sqrt(dot_product(M,M))   ! moment amplitude
vN = M/mamp     ! unit vector of moment axis
angle = DELTA_T*0.01*mamp
call rotate(v0,vN,angle,v)
write(nflog,'(a,3f8.4)') 'rotate: ',v0
write(nflog,'(a,3f8.4)') 'about:  ',vN
write(nflog,'(a,3f8.4)') 'angle:  ',angle
write(nflog,'(a,3f8.4)') 'yields: ',v

end subroutine

!-----------------------------------------------------------------------------------------
! Falpha > 0, therefore x > 0 => fcell > 0
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function fcell(x)
real(REAL_KIND) :: x
logical, parameter :: linear = .false.

if (linear) then
	fcell = Falpha*x
else
	if (x >= 0) then
		fcell = Falpha*x**2
	else
		fcell = -Falpha*x**2
	endif
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function distance(p1,p2)
real(REAL_KIND) :: p1(:), p2(:)

distance = sqrt(dot_product(p1-p2,p1-p2))
end function


!-----------------------------------------------------------------------------------------
! Get total force on cell (i,j)
! 
!           o (i,j+1)
!           |
!           | F_U
!           ^
! (i-1,j)   |       (i+1,j)
!  o---->---O--->----o
!     F_L   |   F_R
!           ^
!           | F_D
!           |
!           o (i,j-1)
!
!
! nd = 2 (2D) or 3 (3D)
! 2D case: the variable numbering is:
!	kx = (j-1)*2*Ncirc + (i-1)*2 + 1 for the x value
!	ky = (j-1)*2*Ncirc + (i-1)*2 + 2 for the y value
!   The variables are v2D(k): k = 1,..,2*Ncirc*Nlong
! 3D case: the circ direction wraps around, and the
! left neighbour of icirc=1 is icirc=Ncirc
! Note: i = icirc, j = ilong
!	kx = (j-1)*3*Ncirc + (i-1)*3 + 1 for the x value
!	ky = (j-1)*3*Ncirc + (i-1)*3 + 2 for the y value
!	kz = (j-1)*3*Ncirc + (i-1)*3 + 3 for the z value
! Need to add an internal pressure force
!-----------------------------------------------------------------------------------------
subroutine getForce(nd,v,i,j,Fsum) bind(C,name='getForce')
use, intrinsic :: iso_c_binding
integer(c_int) :: i, j, nd
real(c_double) :: v(*), F_L(3), F_R(3), F_D(3), F_U(3), F_P(3), Fsum(3)
real(REAL_KIND) :: d, factor, R, c(3), cn(3), fx, fy, fz, v_F(3)
real(REAL_KIND) :: v_L(3), v_R(3), v_P(3), area, sum
integer :: kx, ky, kz		! cell (i,j)
integer :: kxn, kyn, kzn	! neighbour cell (Left, Right, Down or Up)
integer :: in, ic
logical :: dbug = .false.

!R = Rinitial
sum = 0
do ic = 1,Ncirc
	kx = (j-1)*nd*Ncirc + (ic-1)*nd + 1 
	kz = kx + 2
	sum = sum + v(kx)**2 + v(kz)**2
enddo
R = sqrt(sum/(nd*Ncirc*Nlong))
	
area = dx_init*dx_init
if (nd == 2) then
	kx = (j-1)*nd*Ncirc + (i-1)*nd + 1 
	ky = kx + 1
	c = [v(kx),v(ky),0]
else
	kx = (j-1)*nd*Ncirc + (i-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c = [v(kx),v(ky),v(kz)]
endif
Fsum = 0
if (nd == 2) then
	F_L = 0
	if (i == 1) then			! Left boundary constraint
		d = v(kx)
		fx = fcell(-d)
		F_L = [fx, 0.d0, 0.d0]
	elseif (i > 1) then			! F_L (i-1,j) -> (i,j)
		kxn = kx - nd
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = fcell((mcell(i,j)%width(1)+mcell(i-1,j)%width(1))/2 - d)
		if (d > 0) then
			F_L = (c - cn)*factor/d
		else
			write(*,*) 'Error: F_L: getForce: d = 0'
			stop
		endif
		if (nancomponent(F_L)) then
			write(*,*) 'bad F_L'
			stop
		endif
	endif
	
	F_R = 0
	if (i < Ncirc) then	! F_R (i+1,j) -> (i,j)
		kxn = kx + nd
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = fcell((mcell(i,j)%width(1)+mcell(i+1,j)%width(1))/2 - d)
		if (d > 0) then
			F_R = (c - cn)*factor/d
		else
			write(*,*) 'Error: F_R: getForce: d = 0'
			stop
		endif
		if (nancomponent(F_R)) then
			write(*,*) 'bad F_R'
			write(*,*) i,j,kx,ky,kxn,kyn
			write(*,*) 'c, cn: ',c,cn
			write(*,*) 'd: ',d
			stop
		endif
	endif
	
	F_D = 0
	if (j == 1) then			! Bottom boundary constraint
		d = v(ky)
		fy = fcell(-d)
		F_D = [0.d0, fy, 0.d0]
	elseif (j > 1) then		! F_D (i,j-1) -> (i,j)
		kxn = kx - nd*Ncirc
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = fcell((mcell(i,j)%width(2)+mcell(i,j-1)%width(2))/2 - d)
		if (d > 0) then
			F_D = (c - cn)*factor/d
		else
			write(*,*) 'Error: F_D: getForce: d = 0'
			stop
		endif
		if (nancomponent(F_D)) then
			write(*,*) 'bad F_D'
			stop
		endif
	endif
	
	F_U = 0
	if (j < Nlong) then	! F_U (i,j+1) -> (i,j)
		kxn = kx + nd*Ncirc
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = fcell((mcell(i,j)%width(2)+mcell(i,j+1)%width(2))/2 - d)
		if (d > 0) then
			F_U = (c - cn)*factor/d
		else
			write(*,*) 'Error: F_U: getForce: d = 0'
			stop
		endif
		if (nancomponent(F_U)) then
			write(*,*) 'bad F_U'
			stop
		endif
	endif
else
	F_L = 0
	if (i == 1) then			! Left wrap: (Ncirc,j) -> (1,j)
		in = Ncirc
		kxn = (j-1)*nd*Ncirc + (in-1)*nd + 1 
		kyn = kxn + 1
		kzn = kxn + 2
	elseif (i > 1) then			!            (i-1,j) -> (i,j)
		in = i - 1
		kxn = kx - nd
		kyn = kxn + 1
		kzn = kxn + 2
	endif
	cn = [v(kxn),v(kyn),v(kzn)]
	v_F = c - cn
	d = sqrt(dot_product(v_F,v_F))
	v_F = v_F/d
	v_L = v_F
	factor = fcell((mcell(i,j)%width(1)+mcell(in,j)%width(1))/2 - d)
	F_L = factor*v_F
	if (nancomponent(F_L)) then
		write(*,*) 'bad F_L'
		stop
	endif
	if (dbug) then
		write(*,*) 'F_L: i,j,kx,kxn: ',i,j,kx,kxn
		write(*,'(a,3e12.3)') 'c: ',c 
		write(*,'(a,3e12.3)') 'cn: ',cn 
		write(*,'(a,3e12.3)') 'dtarget: ',mcell(i,j)%width(1), mcell(in,j)%width(1), (mcell(i,j)%width(1)+mcell(in,j)%width(1))/2
		write(*,'(a,3e12.3)') 'd: ',d,factor,(mcell(i,j)%width(1)+mcell(in,j)%width(1))/2 - d
	endif

	F_R = 0
	if (i < Ncirc) then	!             (i,j) <- (i+1,j)
		in = i+1
		kxn = kx + nd
		kyn = kxn + 1
		kzn = kxn + 2
	else				! right wrap: (Ncirc,j) <- (1,j)
		in = 1
		kxn = (j-1)*nd*Ncirc + (in-1)*nd + 1 
		kyn = kxn + 1
		kzn = kxn + 2
	endif	
	cn = [v(kxn),v(kyn),v(kzn)]
	v_F = c - cn
	d = sqrt(dot_product(v_F,v_F))
	v_F = v_F/d
	v_R = v_F
	factor = fcell((mcell(i,j)%width(1)+mcell(in,j)%width(1))/2 - d)
	F_R = factor*V_F
	if (nancomponent(F_R)) then
		write(*,*) 'bad F_R'
		write(*,*) i,j,kx,ky,kxn,kyn
		write(*,*) 'c, cn: ',c,cn
		write(*,*) 'd: ',d
		stop
	endif
	if (dbug) then
		write(*,*) 'F_R: i,j,kx,kxn: ',i,j,kx,kxn
		write(*,'(a,3e12.3)') 'c: ',c 
		write(*,'(a,3e12.3)') 'cn: ',cn 
		write(*,'(a,3e12.3)') 'dtarget: ',mcell(i,j)%width(1), mcell(in,j)%width(1), (mcell(i,j)%width(1)+mcell(in,j)%width(1))/2
		write(*,'(a,3e12.3)') 'd: ',d,factor
	endif
	
	F_D = 0
	if (j == 1) then			! Bottom boundary constraint
		if (i == 1) then
			fx = 10*fcell(v0(i,1) - v(kx))
			fy = 10*fcell(v0(i,2) - v(ky))
			fz = 10*fcell(v0(i,3) - v(kz))
		else
			fx = 10*fcell(v0(i,1) - v(kx))
			fy = 10*fcell(v0(i,2) - v(ky))
			fz = 10*fcell(v0(i,3) - v(kz))
		endif			
		F_D = [fx, fy, fz]
	else		! F_D (i,j-1) -> (i,j)
		kxn = kx - nd*Ncirc
		kyn = kxn + 1
		kzn = kxn + 2
		cn = [v(kxn),v(kyn),v(kzn)]
		v_F = c - cn
		d = sqrt(dot_product(v_F,v_F))
		v_F = v_F/d
		factor = fcell((mcell(i,j)%width(2)+mcell(i,j-1)%width(2))/2 - d)
		F_D = factor*v_F
		if (nancomponent(F_D)) then
			write(*,*) 'bad F_D'
			stop
		endif
	endif
	
	F_U = 0
	if (j < Nlong) then	! (i,j+1) -> (i,j)
		kxn = kx + nd*Ncirc
		kyn = kxn + 1
		kzn = kxn + 2
		cn = [v(kxn),v(kyn),v(kzn)]
		v_F = c - cn
		d = sqrt(dot_product(v_F,v_F))
		v_F = v_F/d
		factor = fcell((mcell(i,j)%width(2)+mcell(i,j+1)%width(2))/2 - d)
		F_U = factor*v_F
		if (nancomponent(F_U)) then
			write(*,*) 'bad F_U'
			stop
		endif
		if (dbug) then
			write(*,*) 'F_U: i,j,kx,kxn: ',i,j,kx,kxn
			write(*,'(a,3e12.3)') 'c: ',c 
			write(*,'(a,3e12.3)') 'cn: ',cn 
			write(*,'(a,3e12.3)') 'dtarget: ',mcell(i,j)%width(1), mcell(in,j)%width(1), (mcell(i,j)%width(1)+mcell(in,j)%width(1))/2
			write(*,'(a,2e12.3)') 'd: ',d,factor
		endif
	else
		F_U = [0.0d0,tension,0.0d0]		! try adding a tension force
!		fx = 10*fcell(v1(i,1) - v(kx))
!		fy = 10*fcell(v1(i,2) - v(ky))
!		fz = 10*fcell(v1(i,3) - v(kz))
!		F_U = [fx, fy, fz]
	endif
endif
! Pressure force
v_P = c
v_P(2) = 0
d = sqrt(dot_product(v_P,v_P))
v_P = v_P/d
!factor = (Rinitial - d)
factor = 1
F_P = Pressure*area*factor*v_P
!write(*,'(a,3f8.3)') 'F_P: ',F_P
Fsum = F_L + F_R + F_D + F_U + F_P
if (dbug) then
	write(*,'(7i6)') nd,i,j,kx,ky,kxn,kyn
	write(*,'(a,3e12.3)') 'F_L: ',F_L(1:nd)
	write(*,'(a,3e12.3)') 'F_R: ',F_R(1:nd)
	write(*,'(a,3e12.3)') 'F_D: ',F_D(1:nd)
	write(*,'(a,3e12.3)') 'F_U: ',F_U(1:nd)
endif
!if (Fsum(1) == 0) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
logical function nancomponent(v)
real(REAL_KIND) :: v(3)

nancomponent =  (isnan(v(1)) .or. isnan(v(2)) .or. isnan(v(3)))
end function

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL(REAL_KIND), INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL(REAL_KIND)        :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

end module