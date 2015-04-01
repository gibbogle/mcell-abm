! Solver code using nitsol
module solve

use global
use geometry

implicit none

contains

!-----------------------------------------------------------------------------------------
! xcur is the array containing the current x value, fcur 
! is f(xcur) on output, rpar and ipar are, respectively, real 
! and integer parameter/work arrays for use by the subroutine,
! and itrmf is an integer termination flag.  The meaning of
! itrmf is as follows:
!   0 => normal termination; desired function value calculated.
!   1 => failure to produce f(xcur).
! I can't any reason why this cannot be made parralel with OMP
!-----------------------------------------------------------------------------------------
subroutine fnit(n, xcur, fcur, rpar, ipar, itrmf)
integer :: n, ipar(*), itrmf
real(REAL_KIND) :: xcur(*), fcur(*), rpar(*)
integer :: icirc, ilong, nd, k, kd, nvars, i, kc
real(REAL_KIND) :: F(3), fmax(3)

nd = Ndim
nvars = nd*Ncirc*Nlong
if (n /= nvars) then
	write(*,*) 'ERROR: fnit: n != nvars: ',n,nvars
	stop
endif
nd = Ndim

!$omp parallel do private(ilong, icirc, kd, k, F)
do kc = 1,Ncirc*Nlong
ilong = mod(kc-1,Nlong) + 1
icirc = (kc - ilong)/Nlong + 1
!do ilong = 1,Nlong
!	do icirc = 1,Ncirc
		call getForce(nd,xcur,icirc,ilong,F)
		do kd = 1,nd
			k = (ilong-1)*nd*Ncirc + (icirc-1)*nd + kd 
			fcur(k) = F(kd)
		enddo
!	enddo
!enddo
enddo
!$omp end parallel do

if (istep > nramp) then
	fb = 0
	call getBendForces(xcur,fb)
	fcur(1:n) = fcur(1:n) + fb
	fb = 0
	call getCollisionForces(xcur,fb)
	fcur(1:n) = fcur(1:n) + fb
endif
itrmf = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine nitsolve(v,res)
implicit none
external ddot
external dnrm2
real(REAL_KIND) :: ddot, dnrm2
real(REAL_KIND) :: v(*)
integer :: n, k, i, j, kd
real(REAL_KIND) :: res, F(3)
integer :: input(10), info(6), iterm, ipar(1000)
real(REAL_KIND) :: ftol, stptol, rpar(1000)

n = Ndim*Ncirc*Nlong
!kdmax = 20
!nwork = n*(kdmax+5)+kdmax*(kdmax+3) + 1000	! just for good luck
!allocate(rwork(nwork))
rwork = 0
input = 0
info = 0
ipar = 0
rpar = 0
input(3) = nitsolver
ftol = 2.0d-5
stptol = 2.0d-5
call nitsol(n, v, fnit, jacv, ftol, stptol, input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)
end subroutine

!-----------------------------------------------------------------------------------------
! Dummy, not used
!-----------------------------------------------------------------------------------------
subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
integer :: n, ijob, ipar(*), itrmjv
real(REAL_KIND) :: xcur(*), fcur(*), v, z, rpar(*)

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

end module
