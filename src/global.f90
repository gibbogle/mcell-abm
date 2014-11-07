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
	real(REAL_KIND) :: area
	real(REAL_KIND) :: width(3)
	real(REAL_KIND) :: vert(8,3)	! vertices
    real(REAL_KIND) :: vx(3)		! X, Y, Z axes of the block
    real(REAL_KIND) :: vy(3)
    real(REAL_KIND) :: vz(3)
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
	real(c_double) :: centre(3)
	type(point_type) :: vertex(8)
    real(c_double) :: width(3)
    real(c_double) :: vx(3)
    real(c_double) :: vy(3)
    real(c_double) :: vz(3)
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
integer, parameter :: nramp = 10

character*(128) :: inputfile
character*(128) :: outputfile

integer :: Ncirc, Nlong, isolver
real(REAL_KIND) :: dx_init

integer :: NX, NY, NZ
integer :: Mnodes, ncpu_input, ncells, nsteps, istep, nitsolver
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
real(REAL_KIND) :: pressure, tension, Falpha_axial, Falpha_shear, Falpha_bend
real(REAL_KIND) :: Rinitial

character*(128) :: logfile
character*(2048) :: logmsg
character*(128) :: cellmlfile

type(mcell_type), allocatable, target :: mcell(:,:)
type(mcell_type), allocatable, target :: mcell0(:,:)
real(REAL_KIND), allocatable, target :: v2D(:), v3D(:), vbase(:,:)	!, v1(:,:)
real(REAL_KIND), allocatable :: Jac(:,:)
real(REAL_KIND),allocatable :: cellml_state0(:)
real(REAL_KIND), allocatable :: rwork(:)
real(REAL_KIND), allocatable  :: fb(:)

!type(cell_type), allocatable, target :: cell_list(:)
!type(cell_type), allocatable, target :: cell_list0(:)
!integer, allocatable :: s1s2(:,:,:)

real(REAL_KIND) :: run_dmax, overlap_average, overlap_max
integer :: run_kmax, run_maxstep
integer :: j_deriv
integer :: ic1_min, jl1_min, ic2_min, jl2_min
logical :: may_collide, first_collide

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
  deallocate(C, Kp,Lp)
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
! Falpha_axial > 0, therefore overlap = compression = x > 0 => faxial > 0
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function faxial(x)
real(REAL_KIND) :: x
logical, parameter :: linear = .false.

if (linear) then
	faxial = Falpha_axial*x
else
	if (x >= 0) then
		faxial = Falpha_axial*x**2
	else
		faxial = -Falpha_axial*x**2
	endif
endif
end function

!-----------------------------------------------------------------------------------------
! Falpha_shear > 0, therefore compression = x > 0 => fshear > 0
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function fshear(x)
real(REAL_KIND) :: x
logical, parameter :: linear = .true.

if (linear) then
	fshear = Falpha_shear*x
else
	if (x >= 0) then
		fshear = Falpha_shear*x**2
	else
		fshear = -Falpha_shear*x**2
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
! Note: i = icirc, j = jlong
!	kx = (j-1)*3*Ncirc + (i-1)*3 + 1 for the x value
!	ky = (j-1)*3*Ncirc + (i-1)*3 + 2 for the y value
!	kz = (j-1)*3*Ncirc + (i-1)*3 + 3 for the z value
! Optionally add shear force and an internal pressure force
! Note: bending force is not added here, because it is computed as moments, and the
! moment forces associated with a given cell are also applied to its four neighbours.
!-----------------------------------------------------------------------------------------
subroutine getForce(nd,v,icirc,jlong,Fsum) !bind(C,name='getForce')
!use, intrinsic :: iso_c_binding
!integer(c_int) :: icirc, jlong, nd
!real(c_double) :: v(*), Fsum(3)
integer :: icirc, jlong, nd
real(REAL_KIND) :: v(*), Fsum(3)
real(REAL_KIND) :: F_L(3), F_R(3), F_D(3), F_U(3), F_S(3), F_P(3)
real(REAL_KIND) :: d, factor, R, c(3), cn(3), fx, fy, fz, v_F(3)
real(REAL_KIND) :: v_P(3), area, sum, press, tens, rampfactor
integer :: kx, ky, kz		! cell (icirc,j)
integer :: kxn, kyn, kzn	! neighbour cell (Left, Right, Down or Up)
integer :: in, ic

integer :: jl_D, jl_U, ic_L, ic_R, kd, k
real(REAL_KIND) :: v_L(3), v_R(3), v_D(3), v_U(3), v_N(3)
real(REAL_KIND) :: c_LD(3), c_LU(3), c_RD(3), c_RU(3)
real(REAL_KIND) :: c_L(3), c_R(3), c_D(3), c_U(3), v_LR(3), v_DU(3)
real(REAL_KIND) :: d_L, d_R, d_LR, d_DU, w_0, w_L, w_R, del, dF(3)
logical :: is_U, is_D

logical :: compute_shear = .true.
logical :: dbug = .false.

rampfactor = min(1.0,real(istep)/nramp)

!R = Rinitial
sum = 0
do ic = 1,Ncirc
	kx = (jlong-1)*nd*Ncirc + (ic-1)*nd + 1 
	kz = kx + 2
	sum = sum + v(kx)**2 + v(kz)**2
enddo
R = sqrt(sum/(nd*Ncirc*Nlong))
	
area = dx_init*dx_init
if (nd == 2) then
	kx = (jlong-1)*nd*Ncirc + (icirc-1)*nd + 1 
	ky = kx + 1
	c = [v(kx),v(ky),0]
else
	kx = (jlong-1)*nd*Ncirc + (icirc-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c = [v(kx),v(ky),v(kz)]
endif
Fsum = 0
if (nd == 2) then
	F_L = 0
	if (icirc == 1) then			! Left boundary constraint
		d = v(kx)
		fx = faxial(-d)
		F_L = [fx, 0.d0, 0.d0]
	elseif (icirc > 1) then			! F_L (i-1,j) -> (i,j)
		kxn = kx - nd
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = faxial((mcell(icirc,jlong)%width(1)+mcell(icirc-1,jlong)%width(1))/2 - d)
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
	if (icirc < Ncirc) then	! F_R (i+1,j) -> (i,j)
		kxn = kx + nd
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = faxial((mcell(icirc,jlong)%width(1)+mcell(icirc+1,jlong)%width(1))/2 - d)
		if (d > 0) then
			F_R = (c - cn)*factor/d
		else
			write(*,*) 'Error: F_R: getForce: d = 0'
			stop
		endif
		if (nancomponent(F_R)) then
			write(*,*) 'bad F_R'
			write(*,*) icirc,jlong,kx,ky,kxn,kyn
			write(*,*) 'c, cn: ',c,cn
			write(*,*) 'd: ',d
			stop
		endif
	endif
	
	F_D = 0
	if (jlong == 1) then			! Bottom boundary constraint
		d = v(ky)
		fy = faxial(-d)
		F_D = [0.d0, fy, 0.d0]
	elseif (jlong > 1) then		! F_D (i,j-1) -> (i,j)
		kxn = kx - nd*Ncirc
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = faxial((mcell(icirc,jlong)%width(2)+mcell(icirc,jlong-1)%width(2))/2 - d)
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
	if (jlong < Nlong) then	! F_U (icirc,jlong+1) -> (icirc,jlong)
		kxn = kx + nd*Ncirc
		kyn = kxn + 1
		cn = [v(kxn),v(kyn),0]
		d = distance(c,cn)
		factor = faxial((mcell(icirc,jlong)%width(2)+mcell(icirc,jlong+1)%width(2))/2 - d)
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
	if (icirc == 1) then			! Left wrap: (Ncirc,jlong) -> (1,jlong)
		in = Ncirc
	elseif (icirc > 1) then			!            (icirc-1,jlong) -> (icirc,jlong)
		in = icirc - 1
	endif
	kxn = (jlong-1)*nd*Ncirc + (in-1)*nd + 1 
	kyn = kxn + 1
	kzn = kxn + 2
	cn = [v(kxn),v(kyn),v(kzn)]
	v_F = c - cn
	d = sqrt(dot_product(v_F,v_F))
	v_F = v_F/d
	v_L = v_F
	factor = faxial((mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2 - d)
	F_L = factor*v_F
	if (nancomponent(F_L)) then
		write(*,'(a,3i4,6f8.3)') 'bad F_L: istep, icirc, jlong, c, cn: ',istep, icirc,jlong,c,cn
		stop
	endif
	if (dbug) then
		write(*,*) 'F_L: istep,icirc,jlong,kx,kxn: ',istep,icirc,jlong,kx,kxn
		write(*,'(a,3e12.3)') 'c: ',c 
		write(*,'(a,3e12.3)') 'cn: ',cn 
		write(*,'(a,3e12.3)') 'dtarget: ',mcell(icirc,jlong)%width(1), mcell(in,jlong)%width(1), (mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2
		write(*,'(a,3e12.3)') 'd: ',d,factor,(mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2 - d
	endif

	F_R = 0
	if (icirc < Ncirc) then	!             (icirc,jlong) <- (icirc+1,jlong)
		in = icirc+1
		kxn = kx + nd
		kyn = kxn + 1
		kzn = kxn + 2
	else				! right wrap: (Ncirc,j) <- (1,j)
		in = 1
		kxn = (jlong-1)*nd*Ncirc + (in-1)*nd + 1 
		kyn = kxn + 1
		kzn = kxn + 2
	endif	
	cn = [v(kxn),v(kyn),v(kzn)]
	v_F = c - cn
	d = sqrt(dot_product(v_F,v_F))
	v_F = v_F/d
	v_R = v_F
	factor = faxial((mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2 - d)
	F_R = factor*V_F
	if (nancomponent(F_R)) then
		write(*,*) 'bad F_R'
		write(*,*) icirc,jlong,kx,ky,kxn,kyn
		write(*,*) 'c, cn: ',c,cn
		write(*,*) 'd: ',d
		stop
	endif
	if (dbug) then
		write(*,*) 'F_R: icirc,jlong,kx,kxn: ',icirc,jlong,kx,kxn
		write(*,'(a,3e12.3)') 'c: ',c 
		write(*,'(a,3e12.3)') 'cn: ',cn 
		write(*,'(a,3e12.3)') 'dtarget: ',mcell(icirc,jlong)%width(1), mcell(in,jlong)%width(1), (mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2
		write(*,'(a,3e12.3)') 'd: ',d,factor
	endif
	
	F_D = 0
	if (jlong == 1) then			! Bottom boundary constraint
		fx = 0
		fz = 0
		fx = 0.1*faxial(vbase(icirc,1) - v(kx))
		fy = faxial(vbase(icirc,2) - v(ky))
		fz = 0.1*faxial(vbase(icirc,3) - v(kz))	
		F_D = [fx, fy, fz]
	else		! F_D (i,j-1) -> (i,j)
		kxn = kx - nd*Ncirc
		kyn = kxn + 1
		kzn = kxn + 2
		cn = [v(kxn),v(kyn),v(kzn)]
		v_F = c - cn
		d = sqrt(dot_product(v_F,v_F))
		v_F = v_F/d
		factor = faxial((mcell(icirc,jlong)%width(2)+mcell(icirc,jlong-1)%width(2))/2 - d)
		F_D = factor*v_F
		if (nancomponent(F_D)) then
			write(*,*) 'bad F_D'
			stop
		endif
	endif
	
	F_U = 0
	if (jlong < Nlong) then	! (icirc,jlong+1) -> (icirc,jlong)
		kxn = kx + nd*Ncirc
		kyn = kxn + 1
		kzn = kxn + 2
		cn = [v(kxn),v(kyn),v(kzn)]
		v_F = c - cn
		d = sqrt(dot_product(v_F,v_F))
		v_F = v_F/d
		factor = faxial((mcell(icirc,jlong)%width(2)+mcell(icirc,jlong+1)%width(2))/2 - d)
		F_U = factor*v_F
		if (nancomponent(F_U)) then
			write(*,*) 'bad F_U'
			stop
		endif
		if (dbug) then
			write(*,*) 'F_U: icirc,jlong,kx,kxn: ',icirc,jlong,kx,kxn
			write(*,'(a,3e12.3)') 'c: ',c 
			write(*,'(a,3e12.3)') 'cn: ',cn 
			write(*,'(a,3e12.3)') 'dtarget: ',mcell(icirc,jlong)%width(1), mcell(in,jlong)%width(1), (mcell(icirc,jlong)%width(1)+mcell(in,jlong)%width(1))/2
			write(*,'(a,2e12.3)') 'd: ',d,factor
		endif
	else
		F_U = [0.0d0,rampfactor*tension,0.0d0]		! try adding a tension force
!		fx = 10*faxial(v1(i,1) - v(kx))
!		fy = 10*faxial(v1(i,2) - v(ky))
!		fz = 10*faxial(v1(i,3) - v(kz))
!		F_U = [fx, fy, fz]
	endif
endif

ic_L = icirc - 1
if (icirc == 1) ic_L = Ncirc

ic_R = icirc + 1
if (icirc == Ncirc) ic_R = 1

if (jlong == 1) then
	jl_U = jlong + 1
	is_D = .false.
	is_U = .true.
elseif (jlong == Nlong) then
	jl_D = jlong - 1
	is_D = .true.
	is_U = .false.
else
	jl_D = jlong - 1
	jl_U = jlong + 1
	is_D = .true.
	is_U = .true.
endif

F_S = 0
if (compute_shear) then
! Shear force
w_0 = sqrt( mcell(icirc,jlong)%width(1)**2 + mcell(icirc,jlong)%width(2)**2)
if (is_D) then
	kx = (jl_D-1)*nd*Ncirc + (ic_L-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_LD = [v(kx),v(ky),v(kz)]
!	v_L = mcell(icirc,jlong)%centre - mcell(ic_L,jl_D)%centre
	v_L = c - c_LD
	d_L = sqrt(dot_product(v_L,v_L))
	v_L = v_L/d_L
	w_L = sqrt( mcell(ic_L,jl_D)%width(1)**2 + mcell(ic_L,jl_D)%width(2)**2)
	del = (w_0 + w_L)/2 - d_L
	F_S = F_S + fshear(del)*v_L
	
	kx = (jl_D-1)*nd*Ncirc + (ic_R-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_RD = [v(kx),v(ky),v(kz)]
!	v_R = mcell(icirc,jlong)%centre - mcell(ic_R,jl_D)%centre
	v_R = c - c_RD
	d_R = sqrt(dot_product(v_R,v_R))
	v_R = v_R/d_R
	w_R = sqrt( mcell(ic_R,jl_D)%width(1)**2 + mcell(ic_R,jl_D)%width(2)**2)
	del = (w_0 + w_R)/2 - d_R
	F_S = F_S + fshear(del)*v_R
endif
if (is_U) then
	kx = (jl_U-1)*nd*Ncirc + (ic_R-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_RU = [v(kx),v(ky),v(kz)]
!	v_R = mcell(icirc,jlong)%centre - mcell(ic_R,jl_U)%centre
	v_R = c - c_RU
	d_R = sqrt(dot_product(v_R,v_R))
	v_R = v_R/d_R
	w_R = sqrt( mcell(ic_R,jl_U)%width(1)**2 + mcell(ic_R,jl_U)%width(2)**2)
	del = (w_0 + w_R)/2 - d_R
	F_S = F_S + fshear(del)*v_R
	
	kx = (jl_U-1)*nd*Ncirc + (ic_L-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_LU = [v(kx),v(ky),v(kz)]
!	v_L = mcell(icirc,jlong)%centre - mcell(ic_L,jl_U)%centre
	v_L = c - c_LU
	d_L = sqrt(dot_product(v_L,v_L))
	v_L = v_L/d_L
	w_L = sqrt( mcell(ic_L,jl_U)%width(1)**2 + mcell(ic_L,jl_U)%width(2)**2)
	del = (w_0 + w_L)/2 - d_L
	F_S = F_S + fshear(del)*v_L
endif
endif

! Pressure force - needs fixing: v_P should be outward normal vector, need area
if (Pressure > 0) then
	if (is_D) then
		kx = (jl_D-1)*nd*Ncirc + (icirc-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c_D = [v(kx),v(ky),v(kz)]
	else
		c_D = c
	endif
	if (is_U) then
		kx = (jl_U-1)*nd*Ncirc + (icirc-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c_U = [v(kx),v(ky),v(kz)]
	else
		c_U = c
	endif
	v_DU = c_U - c_D
	d_DU = sqrt(dot_product(v_DU,v_DU))
	v_DU = v_DU/d_DU
	kx = (jlong-1)*nd*Ncirc + (ic_L-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_L = [v(kx),v(ky),v(kz)]
	kx = (jlong-1)*nd*Ncirc + (ic_R-1)*nd + 1 
	ky = kx + 1 
	kz = kx + 2
	c_R = [v(kx),v(ky),v(kz)]
	v_LR = c_R - c_L
	d_LR = sqrt(dot_product(v_LR,v_LR))
	v_LR = v_LR/d_LR
	! The pressure force is normal to v_LR and v_DU, and proportional to the area,
	! which is (roughly) proportional to d_LU*d_LR
	! The outward vector is v_LR x v_DU
	call cross_product(v_LR, v_DU, v_P)
	d = sqrt(dot_product(v_P,v_P))
	v_P = v_P/d
!	area = d_LR*d_DU
!	if (jlong == 1 .or. jlong == Nlong) area = 2*area
	area = mcell(icirc,jlong)%area
	F_P = Pressure*area*rampfactor*v_P
!	if (jlong == Nlong .or. jlong == Nlong - 2) then
!	if (icirc == 1) then
!		write(nflog,'(2i4,6f8.4)') icirc,jlong,v_P,F_P
!	endif
else
	F_P = 0
endif		

!write(*,'(a,3f8.3)') 'F_P: ',F_P
Fsum = F_L + F_R + F_D + F_U + F_S + F_P
if (dbug) then
	write(*,'(7i6)') nd,icirc,jlong,kx,ky,kxn,kyn
	write(*,'(a,3e12.3)') 'F_L: ',F_L(1:nd)
	write(*,'(a,3e12.3)') 'F_R: ',F_R(1:nd)
	write(*,'(a,3e12.3)') 'F_D: ',F_D(1:nd)
	write(*,'(a,3e12.3)') 'F_U: ',F_U(1:nd)
endif
!if (Fsum(1) == 0) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getShearForces(fsh)
real(REAL_KIND) :: fsh(:)
integer :: icirc, ilong, il_D, il_U, ic_L, ic_R, kd, k
real(REAL_KIND) :: v_L(3), v_R(3), d_L, d_R, w_0, w_L, w_R, del, dF(3)
logical :: is_U, is_D

do ilong = 1,Nlong
	if (ilong == 1) then
		il_U = ilong + 1
		is_D = .false.
		is_U = .true.
	elseif (ilong == Nlong) then
		il_D = ilong - 1
		is_D = .true.
		is_U = .false.
	else
		il_D = ilong - 1
		il_U = ilong + 1
		is_D = .true.
		is_U = .true.
	endif
	do icirc = 1,Ncirc
		dF = 0
		w_0 = sqrt( mcell(icirc,ilong)%width(1)**2 + mcell(icirc,ilong)%width(2)**2)
		ic_L = icirc - 1
		if (icirc == 1) ic_L = Ncirc
		ic_R = icirc + 1
		if (icirc == Ncirc) ic_R = 1
		if (is_D) then
			v_L = mcell(icirc,ilong)%centre - mcell(ic_L,il_D)%centre
			d_L = sqrt(dot_product(v_L,v_L))
			v_L = v_L/d_L
			! compression = expected separation - d
			w_L = sqrt( mcell(ic_L,il_D)%width(1)**2 + mcell(ic_L,il_D)%width(2)**2)
			del = (w_0 + w_L)/2 - d_L
			dF = dF + fshear(del)*v_L
			v_R = mcell(icirc,ilong)%centre - mcell(ic_R,il_D)%centre
			d_R = sqrt(dot_product(v_R,v_R))
			v_R = v_R/d_R
			! compression = expected separation - d
			w_R = sqrt( mcell(ic_R,il_D)%width(1)**2 + mcell(ic_R,il_D)%width(2)**2)
			del = (w_0 + w_R)/2 - d_R
			dF = dF + fshear(del)*v_R
		endif
		if (is_U) then
			v_R = mcell(icirc,ilong)%centre - mcell(ic_R,il_U)%centre
			d_R = sqrt(dot_product(v_R,v_R))
			v_R = v_R/d_R
			! compression = expected separation - d
			w_R = sqrt( mcell(ic_R,il_U)%width(1)**2 + mcell(ic_R,il_U)%width(2)**2)
			del = (w_0 + w_R)/2 - d_R
			dF = dF + fshear(del)*v_R
			v_L = mcell(icirc,ilong)%centre - mcell(ic_L,il_U)%centre
			d_L = sqrt(dot_product(v_L,v_L))
			v_L = v_L/d_L
			! compression = expected separation - d
			w_L = sqrt( mcell(ic_L,il_U)%width(1)**2 + mcell(ic_L,il_U)%width(2)**2)
			del = (w_0 + w_L)/2 - d_L
			dF = dF + fshear(del)*v_L
		endif
		do kd = 1,Ndim
			k = (ilong-1)*Ndim*Ncirc + (icirc-1)*Ndim + kd 
			fsh(k) = dF(kd)
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! We need to look at pairs of cells that are close enough to experience repulsive force.
!-----------------------------------------------------------------------------------------
subroutine getCollisionForces(v,fc)
real(REAL_KIND) :: v(*), fc(:)
integer :: nd
real(REAL_KIND) :: c0(3), c1(3), c2(3), v1(3), v2(3)
real(REAL_KIND) :: d, w1, w2, dthresh, dlimit, theta, fun, f(3), dmin
integer :: ic1, ic2, icirc1, icirc2, jlong1, jlong2, kx, ky, kz, kd, k
integer :: ic1_lo, ic1_hi, ic2_lo, ic2_hi, jl1_lo, jl1_hi, jl2_lo, jl2_hi
integer :: Nthresh = 30
real(REAL_KIND) :: alpha_collide = 0.0001

nd = Ndim
if (.not.may_collide) then
	c0 = 0.5*(mcell(1,Nlong/2)%centre + mcell(Ncirc/2,Nlong/2)%centre)
	c1 = 0.5*(mcell(1,1)%centre + mcell(Ncirc/2,1)%centre)
	c2 = 0.5*(mcell(1,Nlong)%centre + mcell(Ncirc/2,Nlong)%centre)
!	v1 = c1 - c0
!	d = sqrt(dot_product(v1,v1))
!	v1 = v1/d
!	v2 = c2 - c0
!	d = sqrt(dot_product(v2,v2))
!	v2 = v2/d
!	may_collide = (sqrt(dot_product(v1,v2)) > 0.9)	! a very crude test
!	write(nflog,'(a,6f8.3)') 'v1,v2: ',v1,v2
	v2 = c2 - c1
	may_collide = (sqrt(dot_product(v2,v2)) < 5*Rinitial)
	if (.not.may_collide) return
	first_collide = .true.
endif
!write(nflog,*) 'getCollisionForces'	
if (first_collide) then
	jl1_lo = 1
	jl1_hi = Nlong - Nthresh
	jl2_lo = 2
	jl2_hi = Nlong
	ic1_lo = 1
	ic1_hi = Ncirc
	ic2_lo = 1
	ic2_hi = Ncirc
	first_collide = .false.
else
	jl1_lo = max(1,jl1_min-3)
	jl1_hi = min(Nlong,jl1_min+3)
	ic1_lo = ic1_min - 3
	ic1_hi = ic1_min + 3
	jl2_lo = max(1,jl2_min-3)
	jl2_hi = min(Nlong,jl2_min+3)
	ic2_lo = ic2_min - 3
	ic2_hi = ic2_min + 3
endif
dmin = 1.0e10
!do jlong1 = 1,Nlong - Nthresh
!	do icirc1 = 1,Ncirc
do jlong1 = jl1_lo,jl1_hi
	do ic1 = ic1_lo,ic1_hi
		icirc1 = ic1
		if (icirc1 < 1) icirc1 = icirc1 + Ncirc
		if (icirc1 > Ncirc) icirc1 = icirc1 - Ncirc
!		c1 = mcell(icirc1,jlong1)%centre
		kx = (jlong1-1)*nd*Ncirc + (icirc1-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c1 = [v(kx),v(ky),v(kz)]
		w1 = mcell(icirc1,jlong1)%width(3)
!		do jlong2 = jlong1 + Nthresh, Nlong
!			do icirc2 = 1,Ncirc
		do jlong2 = jl2_lo,jl2_hi
			if (abs(jlong1-jlong2) < Nthresh) cycle
			do ic2 = ic2_lo,ic2_hi
				icirc2 = ic2
				if (icirc2 < 1) icirc2 = icirc2 + Ncirc
				if (icirc2 > Ncirc) icirc2 = icirc2 - Ncirc
!				c2 = mcell(icirc2,jlong2)%centre
				kx = (jlong2-1)*nd*Ncirc + (icirc2-1)*nd + 1 
				ky = kx + 1 
				kz = kx + 2
				c2 = [v(kx),v(ky),v(kz)]
				w2 = mcell(icirc2,jlong2)%width(3)
				v2 = c2 - c1
				d = sqrt(dot_product(v2,v2))
				if (d < dmin) then
					dmin = d
					ic1_min = icirc1
					jl1_min = jlong1
					ic2_min = icirc2
					jl2_min = jlong2
!					write(nflog,*) 'dmin: ',ic1_min,jl1_min,ic2_min,jl2_min,dmin
				endif
				dthresh = 1*(w1 + w2)
				dlimit = (w1 + w2)/2
				if (d < dthresh) then
					theta = (PI/2)*(dthresh - d)/(dthresh - dlimit)
					if (theta >= PI/2) then
						write(nflog,*) 'Bad theta: ',theta,d,dthresh,dlimit
						stop
					endif
					fun = tan(theta)
					f = alpha_collide*fun*v2/d
!					write(nflog,'(4i3,6f8.4)') jlong1,icirc1,jlong2,icirc2,d,theta,fun,f
					! f(:) is the force on c2, -f(:) is the force on c1
					do kd = 1,Ndim
						k = (jlong1-1)*nd*Ncirc + (icirc1-1)*nd + kd 
						fc(k) = fc(k) - f(kd)
						k = (jlong2-1)*nd*Ncirc + (icirc2-1)*nd + kd 
						fc(k) = fc(k) + f(kd)
					enddo
				endif
			enddo
		enddo
	enddo
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getBendForces(v,fb)
real(REAL_KIND) :: v(*), fb(:)
integer :: icirc, jlong, jl_D, jl_U, ic_L, ic_R, nd, kd, k, kx, ky, kz
real(REAL_KIND) :: c(3), c_L(3), c_R(3), c_D(3), c_U(3)
real(REAL_KIND) :: v_L(3), v_R(3), v_D(3), v_U(3), v_N(3), d, F_L(3), F_R(3), F_D(3), F_U(3)
logical :: is_U, is_D

nd = Ndim
fb = 0
do jlong = 1,Nlong
	if (jlong == 1) then
		jl_U = jlong + 1
		is_D = .false.
		is_U = .true.
	elseif (jlong == Nlong) then
		jl_D = jlong - 1
		is_D = .true.
		is_U = .false.
	else
		jl_D = jlong - 1
		jl_U = jlong + 1
		is_D = .true.
		is_U = .true.
	endif
	do icirc = 1,Ncirc
		kx = (jlong-1)*nd*Ncirc + (icirc-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c = [v(kx),v(ky),v(kz)]
		
		ic_L = icirc - 1
		if (icirc == 1) ic_L = Ncirc		
		kx = (jlong-1)*nd*Ncirc + (ic_L-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c_L = [v(kx),v(ky),v(kz)]
		v_L = c_L - c
		d = sqrt(dot_product(v_L,v_L))
		v_L = v_L/d
		
		ic_R = icirc + 1
		if (icirc == Ncirc) ic_R = 1
		kx = (jlong-1)*nd*Ncirc + (ic_R-1)*nd + 1 
		ky = kx + 1 
		kz = kx + 2
		c_R = [v(kx),v(ky),v(kz)]
		v_R = c_R - c
		d = sqrt(dot_product(v_R,v_R))
		v_R = v_R/d
		
!		cosa = dot_product(v_L,v_R)
!		theta = acos(cosa)
!		call cross_product(v_R,v_L,v_N)	! v_N is normal to v_R and v_L, in the up (long) direction if theta > 0
!		d_N = sqrt(dot_product(v_N,v_N))
!		v_N = v_N/d_N	! unit vector
!		! At the left point, L, the left force contribution F_L is ~ v_N x V_L, and this gives -F_L at C
!		! At the right point, R, the right force contribution F_R is ~ V_R x v_N, and this gives -F_R at C
!		! Force magnitude is proportional to theta
!		call cross_product(v_L,v_N,F_L)
!		dF = dF + alpha_bend*theta*F_L
!		call cross_product(v_R,v_N,F_R)
!		dF = dF - alpha_bend*theta*F_R
		
		call bending_force(v_L, v_R, F_L, F_R)
		
		if (is_D .and. is_U) then
			kx = (jl_D-1)*nd*Ncirc + (icirc-1)*nd + 1 
			ky = kx + 1 
			kz = kx + 2
			c_D = [v(kx),v(ky),v(kz)]
			v_D = c_D - c
			d = sqrt(dot_product(v_D,v_D))
			v_D = v_D/d
			
			kx = (jl_U-1)*nd*Ncirc + (icirc-1)*nd + 1 
			ky = kx + 1 
			kz = kx + 2
			c_U = [v(kx),v(ky),v(kz)]
			v_U = c_U - c
			d = sqrt(dot_product(v_U,v_U))
			v_U = v_U/d
			
			call bending_force(v_D, v_U, F_D, F_U)
		else
			F_D = 0
			F_U = 0
		endif
			
		
		do kd = 1,Ndim
			k = (jlong-1)*nd*Ncirc + (icirc-1)*nd + kd 
			fb(k) = fb(k) + F_L(kd) + F_R(kd)
			fb(k) = fb(k) + F_D(kd) + F_U(kd)
			k = (jlong-1)*nd*Ncirc + (ic_L-1)*nd + kd 
			fb(k) = fb(k) - F_L(kd)
			k = (jlong-1)*nd*Ncirc + (ic_R-1)*nd + kd 
			fb(k) = fb(k) - F_R(kd)
			if (is_D .and. is_U) then
				k = (jl_D-1)*nd*Ncirc + (icirc-1)*nd + kd 
				fb(k) = fb(k) - F_D(kd)
				k = (jl_U-1)*nd*Ncirc + (icirc-1)*nd + kd 
				fb(k) = fb(k) - F_U(kd)
			endif
		enddo
	enddo
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! v1 and v2 are unit vectors from centre point C to adjacent centre points (L R or D U)
! The force F1 is applied at C, -F1 at point 1
! The force F2 is applied at C, -F2 at point 2
!-----------------------------------------------------------------------------------------
subroutine bending_force(v1, v2, F1, F2)
real(REAL_KIND) :: v1(3), v2(3), F1(3), F2(3)
real(REAL_KIND) :: cosa, theta, v_N(3), d_N

cosa = dot_product(v1,v2)
theta = acos(abs(cosa))
call cross_product(v2,v1,v_N)	! v_N is normal to v2 and v1
d_N = sqrt(dot_product(v_N,v_N))
if (d_N > 1.0d-6) then
	v_N = v_N/d_N	! unit normal vector
else
	F1 = 0
	F2 = 0
	return
endif
!write(*,'(a,3f10.6)') 'v1: ',v1
!write(*,'(a,3f10.6)') 'v2: ',v2
!write(*,'(a,3f10.6)') 'vN: ',v_N
call cross_product(v1,v_N,F1)
!write(*,'(a,3f8.4)') 'F1: ',F1
F1 = Falpha_bend*theta*F1
call cross_product(v_N,v2,F2)
!write(*,'(a,3f8.4)') 'F2: ',F2
!write(*,'(a,2e12.3)') 'cosa, theta: ',cosa,theta
F2 = Falpha_bend*theta*F2
!stop
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