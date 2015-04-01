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
real(REAL_KIND) :: growth_rate(2,3)

character*(128) :: logfile
character*(2048) :: logmsg
character*(128) :: cellmlfile

type(mcell_type), allocatable, target :: mcell(:,:)
type(mcell_type), allocatable, target :: mcell0(:,:)
real(REAL_KIND), allocatable, target :: v2D(:), v3D(:), vbase(:,:)	!, v1(:,:)
!real(REAL_KIND), allocatable :: Jac(:,:)
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

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function distance(p1,p2)
real(REAL_KIND) :: p1(:), p2(:)

distance = sqrt(dot_product(p1-p2,p1-p2))
end function

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
real(REAL_KIND) function det3(A)
real(REAL_KIND) :: A(3,3)

det3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
      -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
      +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function

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