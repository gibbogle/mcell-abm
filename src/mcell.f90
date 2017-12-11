! ABM for mouse heart-tube cells, off-lattice 
module mcell_mod
use global
use geometry
use solve

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
real(REAL_KIND) :: hours
integer :: nb0, nt_anim, iuse_cellml
character*(256) :: cellmlfile

! Set up growth data arrays
if (allocated(growth_file)) deallocate(growth_file)
if (allocated(growth_time)) deallocate(growth_time)
if (allocated(sector)) deallocate(sector)
if (allocated(section)) deallocate(section)
if (allocated(growth_rate)) deallocate(growth_rate)
ngrowthtimes = 5
allocate(growth_file(ngrowthtimes))
allocate(growth_time(ngrowthtimes))

open(nfin,file=inputfile)
!read(nfin,*) nb0
!read(nfin,*) igrow
!read(nfin,*) t_div_median
!read(nfin,*) t_div_shape
read(nfin,*) Ncirc
read(nfin,*) Nlong
read(nfin,*) hours
read(nfin,*) DELTA_T		! mins
read(nfin,*) dx_init
read(nfin,*) rate_factor
read(nfin,*) pressure
read(nfin,*) tension
read(nfin,*) Falpha_axial
read(nfin,*) Falpha_shear
read(nfin,*) Falpha_bend
read(nfin,*) nitsolver
read(nfin,*) seed(1)
read(nfin,*) seed(2)
read(nfin,*) ncpu_input
read(nfin,*) nt_anim
!read(nfin,*) Fdrag
!read(nfin,*) Mdrag
!read(nfin,*) Malpha
!read(nfin,*) Fjigglefactor
!read(nfin,*) Mjigglefactor
read(nfin,*) iuse_cellml
read(nfin,*) cellmlfile
!read(nfin,*) growth_rate(1,1,1)	! dorsal-bottom
!read(nfin,*) growth_rate(1,2,1)	! dorsal-middle
!read(nfin,*) growth_rate(1,3,1)	! dorsal-top
!read(nfin,*) growth_rate(2,1,1)	! ventral-bottom
!read(nfin,*) growth_rate(2,2,1)	! ventral-middle
!read(nfin,*) growth_rate(2,3,1)	! ventral-top
read(nfin,*) growth_file(1)
read(nfin,*) growth_file(2)
read(nfin,*) growth_file(3)
read(nfin,*) growth_file(4)
read(nfin,*) growth_file(5)
close(nfin)

use_cellml = (iuse_cellml == 1)
if (use_cellml) then
	write(logmsg,'(a)') 'CellML model for growth is not yet implemented'
	call logger(logmsg)
	ok = .false.
	return
endif
!write(*,'(a)') growth_file(:)

!NX = nb0
!NY = nb0
!NZ = nb0
!simulate_rotation = (Mdrag > 0)
!simulate_growth = (igrow == 1)
!DELTA_T = deltat/60 ! sec -> min
nsteps = hours*60./DELTA_T
!nsteps = 10000
write(logmsg,*) 'nitsolver: ',nitsolver
call logger(logmsg)
write(logmsg,*) 'hours,nsteps: ',hours,nsteps 
call logger(logmsg)
ok = .true.
end subroutine 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetupGrowthData
integer :: itime, i, isection
real(REAL_KIND) :: t0

do itime = 1,ngrowthtimes
	open(nfgrowth,file=growth_file(itime),status='old')
	read(nfgrowth,*) growth_time(itime)
!	write(*,*) 'growth_time: ',itime,growth_time(itime)
	if (itime == 1) then
		t0 = growth_time(1)
		read(nfgrowth,*) NR
		allocate(sector(NR,2))
		do i = 1,NR
			read(nfgrowth,*) sector(i,:)
		enddo
		read(nfgrowth,*) NL
		allocate(section(NL,2))
		do i = 1,NL
			read(nfgrowth,*) section(i,:)
		enddo
		allocate(growth_rate(NR,NL,ngrowthtimes))
	else
		do i = 1,NR+NL+2
			read(nfgrowth,*)
		enddo
	endif
	do isection = 1,NL
		read(nfgrowth,*) growth_rate(:,isection,itime)
!		write(*,'(i3,8f8.4)') isection,growth_rate(:,isection,itime)
	enddo
	close(nfgrowth)
	! Convert growth_time to time from start of simulation (mins)
	growth_time(itime) = (growth_time(itime) - t0)*24*60	! days --> mins
!	write(*,*) 'itime: ',itime,growth_time(itime)
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: kcell, i, j, error

ok = .true.
initialized = .false.
par_zig_init = .false.
istep = 0

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

call SetupGrowthData
call logger('did SetupGrowthData')

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call RngInitialisation

call PlaceCells
!do kcell = 1,ncells
!	call SetNeighbours(kcell)
!    write(*,*) 'cell: ',kcell, cell_list(kcell)%nbrs
!    write(*,'(20i4)') (cell_list(kcell)%nbrlist(j),j=1,cell_list(kcell)%nbrs)
!enddo
!cell_list0 = cell_list
write(nflog,*) 'ncells: ',ncells

do i = 1,Ncirc
	do j = 1,Nlong
		write(nflog,'(2i4,3f8.3)') i,j,mcell(i,j)%centre
	enddo
enddo
!do j = 1,Nlong
!	do i = 1,Ncirc
!		call makeVertices(i,j)
!	enddo
!enddo

may_collide = .false.
first_collide = .true.

end subroutine

!-----------------------------------------------------------------------------------------
! Base centre = (0,0,0), initial radius = Rinitial
! Base sits on XZ plane, tube initially oriented parallel to the Y axis
!-----------------------------------------------------------------------------------------
subroutine PlaceCells
integer :: icirc, ilong, kpar=0
integer :: n, kdmax, nwork
real(REAL_KIND) :: x, y, z
real(REAL_KIND) :: dtheta, theta, w, area

if (allocated(mcell)) then
	deallocate(mcell)
	deallocate(mcell0)
	deallocate(vbase)
!	deallocate(v1)
endif
if (allocated(v2D)) then
	deallocate(v2D)
endif
if (allocated(v3D)) then
	deallocate(v3D)
endif
!if (allocated(Jac)) then
!	deallocate(Jac)
!endif
allocate(mcell(Ncirc,Nlong))
allocate(mcell0(Ncirc,Nlong))
allocate(vbase(Ncirc,3))
!allocate(v1(Ncirc,3))
if (Ndim == 2) then
	allocate(v2D(2*Ncirc*Nlong))
!	allocate(Jac(2*Ncirc*Nlong,2*Ncirc*Nlong))
else
	allocate(v3D(3*Ncirc*Nlong))
!	allocate(Jac(3*Ncirc*Nlong,3*Ncirc*Nlong))
endif
! nitsol work space
if (allocated(rwork)) then
	deallocate(rwork)
endif
n = Ndim*Ncirc*Nlong
kdmax = 30	! was 20
nwork = n*(kdmax+5)+kdmax*(kdmax+3) + 1000	! just for good luck
allocate(rwork(nwork))

if (allocated(fb)) then
	deallocate(fb)
endif
allocate(fb(n))

z = 0
if (Ndim == 2) then
	! 2D
	w = dx_init
	do icirc = 1,Ncirc
		x = (icirc-1)*w
		do ilong = 1,Nlong
			y = (ilong-1)*w
			mcell(icirc,ilong)%centre = [x,y,z]
			mcell(icirc,ilong)%width = [w, w, w]
		enddo
	enddo
else
	! 3D
	Rinitial = Ncirc*dx_init/(2*PI)
	dtheta = 2*PI/Ncirc
	w = 2*Rinitial*sin(dtheta/2)
	area = w*2*PI*(Rinitial-w/2)/Ncirc	! internal radius = Rinitial - w/2
	do icirc = 1,Ncirc
		theta = (icirc-0.5)*dtheta
		x = Rinitial*cos(theta)
		z = -Rinitial*sin(theta)	! for right-handed x-y-z axes (with conventional x-y, z is out of the screen)
		do ilong = 1,Nlong
			y = (ilong-1)*w
			mcell(icirc,ilong)%centre = [x,y,z]
			mcell(icirc,ilong)%width = [w, w, w]
			mcell(icirc,ilong)%area = area
			mcell(icirc,ilong)%theta = theta*180/PI
		enddo
	enddo
endif
mcell%volume = w**3
Ncells = Ncirc*Nlong
do icirc = 1,Ncirc
	vbase(icirc,:) = mcell(icirc,1)%centre
enddo
do icirc = 1,Ncirc
	do ilong = 1,Nlong
		mcell0(icirc,ilong)%centre = mcell(icirc,ilong)%centre
		mcell0(icirc,ilong)%volume = mcell(icirc,ilong)%volume
		mcell0(icirc,ilong)%width = mcell(icirc,ilong)%width
		mcell0(icirc,ilong)%theta = mcell(icirc,ilong)%theta
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Programmed test case for cell shape change.
!--------------------------------------------------------------------------------
subroutine updateWidths
integer :: ilong, icirc, ivar
real(REAL_KIND) :: x, y, z, ymid, tfactor, xfactor, yfactorx, yfactorz, zfactor, deltay, deltaz, cyz
real(REAL_KIND) :: deltay_max = 0.6	!0.6
real(REAL_KIND) :: deltaz_max = 0.35	!0.35
real(REAL_KIND) :: t_start = 60		! time until bending start (min)
real(REAL_KIND) :: t_var = 30*60	! duration of bending (min)
integer :: nstep_start
integer :: nstep_var

nstep_start = t_start/DELTA_T
nstep_var = t_var/DELTA_T
if (istep <= nstep_start) return
!nstep_var = (deltay_max/0.1)*40.		! sets a valid rate
ivar = istep - nstep_start
ymid = (1 + Nlong)/2.
if (ivar <= nstep_var) then
	tfactor = (ivar-1.0)/(nstep_var - 1.0)
else
	tfactor = 1.0
endif
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		x = mcell0(icirc,ilong)%centre(1)
		xfactor = x/Rinitial
		z = mcell0(icirc,ilong)%centre(3)
		zfactor = z/Rinitial
		yfactorx = 1 - abs(ilong - ymid)/(ymid-1)
		deltay = tfactor*xfactor*yfactorx*deltay_max
		yfactorz = (ilong - ymid)/Nlong
		deltaz = tfactor*zfactor*yfactorz*deltaz_max
		cyz = (1 + deltay)*(1 + deltaz)
!		write(*,'(a,4f8.4)') 'tfactor,xfactor,yfactorx,deltay_max: ',tfactor,xfactor,yfactorx,deltay_max
!		write(*,'(a,3f8.4)') 'deltay,deltaz,cyz: ',deltay,deltaz,cyz
		mcell(icirc,ilong)%width(2) = cyz*mcell0(icirc,ilong)%width(2)
		if (cyz < 1) then		! adjust width(1) and width(3) to preserve volume
			mcell(icirc,ilong)%width(1) = mcell0(icirc,ilong)%width(1)/sqrt(cyz)
			mcell(icirc,ilong)%width(3) = mcell0(icirc,ilong)%width(3)/sqrt(cyz)
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine updateWidths6
integer :: ilong, icirc, ivar
real(REAL_KIND) :: rate, cyz
real(REAL_KIND) :: t_start = 10		! time until bending start (min)
integer :: nstep_start

nstep_start = t_start/DELTA_T
if (istep <= nstep_start) return
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		rate = get_growthrate(ilong,icirc)
		cyz = (1 + DELTA_T*rate)
!		if (ilong == int(ymid)) write(nflog,*) 'cyz: ',cyz
		mcell(icirc,ilong)%width(2) = cyz*mcell(icirc,ilong)%width(2)
		if (cyz < 1) then		! adjust width(1) and width(3) to preserve volume (do not allow block shrinkage)
			mcell(icirc,ilong)%width(1) = mcell(icirc,ilong)%width(1)/sqrt(cyz)
			mcell(icirc,ilong)%width(3) = mcell(icirc,ilong)%width(3)/sqrt(cyz)
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
real(REAL_KIND) function get_growthrate(ilong,icirc) result(rate)
integer :: ilong, icirc
real(REAL_KIND) :: tnow, lo, hi, theta, yfraction, tfraction1, tfraction2
integer :: i, isector, isection, itime1, itime2
real(REAL_KIND) :: theta0, theta1, theta2, alfa
real(REAL_KIND) :: yfraction0, yfraction1, yfraction2, beta
integer ::  isector1, isector2, isection1, isection2

tnow = istep*DELTA_T
itime1 = 0
do i = 2,ngrowthtimes
	if (tnow < growth_time(i)) then
		itime1 = i-1
		itime2 = i
		tfraction1 = (tnow - growth_time(itime1))/(growth_time(itime2) - growth_time(itime1))
		tfraction2 = 1 - tfraction1
		exit
	endif
enddo
if (itime1 == 0) then
	itime1 = ngrowthtimes
	itime2 = ngrowthtimes
	tfraction1 = 1
	tfraction2 = 0
endif

! Given theta = mcell0(icirc,ilong)%theta, how to determine isector = 1,..,NR?
! The sector boundaries are sector(isector,2)
theta = mcell0(icirc,ilong)%theta
yfraction = real(ilong)/nlong
do isector = 1,NR
	lo = sector(isector,1)
	hi = sector(isector,2)
	if (lo > hi) then
		if (theta < 180) then
			lo = lo - 360
		else
			hi = hi + 360
		endif
	endif
	if (lo <= theta .and. theta <= hi) then
		theta0 = (lo + hi)/2
		exit
	endif
enddo
do isection = 1,NL
	lo = section(isection,1)
	hi = section(isection,2)
	if (lo <= yfraction .and. yfraction <= hi) exit
enddo

! Need to interpolate in space as well
!!!!!!!!!!!!!!!!!!!!!!!!  THIS IS ALL SCREWED UP !!!!!!!!!!!!!!!!!!!!
!theta0 = (sector(isector,1) + sector(isector,2))/2
if (theta > theta0) then
	isector1 = isector
	theta1 = theta0
	isector2 = isector1 + 1
	if (isector2 > NR) isector2 = 1
	lo = sector(isector2,1)
	hi = sector(isector2,2)
	if (lo > hi) then
		if (theta < 180) then
			lo = lo - 360
		else
			hi = hi + 360
		endif
	endif
	theta2 = (lo + hi)/2
else
	isector2 = isector
	theta2 = theta0
	isector1 = isector - 1
	if (isector1 < 1) isector1 = NR
	lo = sector(isector1,1)
	hi = sector(isector1,2)
	if (lo > hi) then
		if (theta < 180) then
			lo = lo - 360
		else
			hi = hi + 360
		endif
	endif	
	theta1 = (lo + hi)/2
endif
alfa = (theta - theta1)/(theta2 - theta1)

!if (ilong == nlong/4) then
!!	write(*,'(3i6,f8.4)') istep,icirc,isector,rate
!	write(*,'(4i4,4f8.1,f8.4)') icirc,isector,isector1,isector2,theta0,theta1,theta,theta2,alfa
!	if (icirc == ncirc) stop
!endif    

yfraction0 = (section(isection,1) + section(isection,2))/2
if (yfraction > yfraction0) then
	isection1 = isection
	yfraction1 = yfraction0
	isection2 = min(isection1 + 1, NL)
	yfraction2 = (section(isection2,1) + section(isection2,2))/2
else
	isection2 = isection
	yfraction2 = yfraction0
	isection1 = max(isection2-1,1)
	yfraction1 = (section(isection1,1) + section(isection1,2))/2
endif
if (isection1 /= isection2) then
	beta = (yfraction - yfraction1)/(yfraction2 - yfraction1)
else
	beta = 1
endif

rate = tfraction1*((1-alfa)*(1-beta)*growth_rate(isector1,isection1,itime1) &
                 + (1-alfa)*beta*growth_rate(isector1,isection2,itime1) &
                 + alfa*(1-beta)*growth_rate(isector2,isection1,itime1) &
                 + alfa*beta*growth_rate(isector2,isection2,itime1)) &
     + tfraction2*((1-alfa)*(1-beta)*growth_rate(isector1,isection1,itime2) &
                 + (1-alfa)*beta*growth_rate(isector1,isection2,itime2) &
                 + alfa*(1-beta)*growth_rate(isector2,isection1,itime2) &
                 + alfa*beta*growth_rate(isector2,isection2,itime2))
rate = rate_factor*rate	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,'(a,2e12.3)') 'rates: ',growth_rate(isector,isection,itime1),growth_rate(isector,isection,itime2)
!write(*,'(5i4,f8.4)') istep,ilong,icirc,isector,isection,rate
if (abs(rate) > 0.1) then
	write(*,'(a,i3,f8.1)') 'istep,tnow: ',istep, tnow
	write(*,'(a,2i4,2f8.4)') 'itime1, itime2: ',itime1,itime2,tfraction1,tfraction2
	write(*,'(a,5f8.1)') 'growth_time: ',growth_time(:)
	write(*,'(a,f8.4,i3)') 'theta,isector: ',theta,isector
	write(*,'(a,f8.4,i3)') 'yfraction,isection: ',yfraction,isection
	write(*,'(a,2i3,2e12.3)') 'rates: ',isector,isection,growth_rate(isector,isection,itime1),growth_rate(isector,isection,itime2)
	write(*,'(3f8.4)') growth_rate(:,:,itime2)
	stop
endif
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
real(REAL_KIND) function get_growthrate1(ilong,icirc) result(rate)
integer :: ilong, icirc
integer :: isection	! 1 = bottom half, 2 = top half
real(REAL_KIND) :: x, fx, fy, ymid

ymid = (1 + Nlong)/2.
if (ilong <= ymid) then
	isection = 1
	fy = (ilong-1)/(ymid-1)
else
	isection = 2
	fy = (ilong-ymid)/(Nlong - ymid)
endif
x = mcell0(icirc,ilong)%centre(1)
fx = (x + Rinitial)/(2*Rinitial)
rate = (1-fx)*(1-fy)*growth_rate(2,isection,1) + fx*(1-fy)*growth_rate(1,isection,1)  &
    + (1-fx)*fy*growth_rate(2,isection+1,1) + fx*fy*growth_rate(1,isection+1,1)
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_scene(Nhex, hex_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: Nhex
type(hexahedron) :: hex_list(*)
integer :: icirc, ilong, k
!type(mcell_type), pointer :: p
!real(REAL_KIND) :: a, b, vx, vy, vz, s, c, s2
!integer :: kcell, j, ix, iy

Nhex = 0
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		Nhex = Nhex + 1
		do k = 1,8
			hex_list(Nhex)%vertex(k)%x(:) = mcell(icirc,ilong)%vert(k,:)
		enddo
		hex_list(Nhex)%centre = mcell(icirc,ilong)%centre
		hex_list(Nhex)%width = mcell(icirc,ilong)%width
		hex_list(Nhex)%vx = mcell(icirc,ilong)%vx
		hex_list(Nhex)%vy = mcell(icirc,ilong)%vy
		hex_list(Nhex)%vz = mcell(icirc,ilong)%vz
	enddo
enddo
!nEC_list = ncells
!kcell = 0
!do ix = 1,Ncirc
!do iy = 1,Nlong
!    p => mcell(ix,iy)
!    kcell = kcell + 1
!    j = (kcell-1)*12
!    EC_list(j+1:j+3) = p%centre
!    EC_list(j+4:j+12) = 0
!enddo
!enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_fieldinfo(nxx, axis, fraction, ns, nc, cused, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fieldinfo
use, intrinsic :: iso_c_binding
integer(c_int) :: nxx, axis, ns, nc, cused(*), res
real(c_double) :: fraction

end subroutine

!--------------------------------------------------------------------------------
! Need to transmit medium concentration data.  This could be a separate subroutine.
!--------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, nfdata, nc, fdata, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
real(c_double) :: fraction
integer(c_int) :: axis, nc, nfdata, res
type(FIELD_DATA) :: fdata(*)

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: ires, it, nt, k, kd, icirc, ilong
real(REAL_KIND) :: dt, obj, width(3)
!integer :: kcell
real(REAL_KIND), pointer :: pv(:)
logical :: ok

res = 0
istep = istep + 1
if (mod(istep,10) == 0) then
	write(logmsg,*) 'simulate_step: ',istep
	call logger(logmsg)
endif
ok = .true.

call updateWidths
if (Ndim == 2) then
	pv => v2D
else
	pv => v3D
endif
k = 0
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		do kd = 1,Ndim
			k = k+1
			pv(k) = mcell(icirc,ilong)%centre(kd)		! initial guess
		enddo
!		mcell(icirc,ilong)%width(1) = (1 + delta)*mcell(icirc,ilong)%width(1)	! target widths
!		mcell(icirc,ilong)%width(2) = (1 + delta)*mcell(icirc,ilong)%width(2)
!		mcell(icirc,ilong)%width(3) = (1 + delta)*mcell(icirc,ilong)%width(3)
!		write(*,'(4i4,2f8.4)') i,j,k-1,k,v2D(k-1),v2D(k)	!,mcell(i,j)%centre(1:2)
	enddo
enddo
!vbase = (1 + delta)*vbase	! grow the base only if width(1) and width(2) are growing

!write(*,*) 'v3D'
!write(*,'(3f8.3)') v3D

!write(nflog,*)
!write(nflog,*) 'Centres'
!call newtonsolve(pv,ires)
call nitsolve(pv,obj)
k = 0
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		do kd = 1,Ndim
			k = k+1
!			width(1) = pv(k) - mcell(i,j)%centre(kd)
			mcell(icirc,ilong)%centre(kd) = pv(k)
		enddo
!		k = k+1
!		width(2) = p(k) - mcell(i,j)%centre(2)
!		mcell(i,j)%centre(2) = p(k)
!		write(nflog,'(2i4,6f8.4)') ilong,icirc,mcell(icirc,ilong)%centre(1:Ndim)	!,width(1:Ndim)
	enddo
enddo

!call logger('makeVertices')
!write(nflog,*) 'makeVertices'
!write(nflog,*) 'Vertices'
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		call makeVertices(icirc,ilong)
!		do k = 1,8
!			write(nflog,'(3i4,3f8.4)') ilong,icirc,k,mcell(icirc,ilong)%vert(k,:)
!		enddo
	enddo
enddo

if (ok) then
    res = 0
else
    res = 1
endif
!call logger('end simulate_step')
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NstepsGUI, NcircGUI, NlongGUI, deltat) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NstepsGUI, NcircGUI, NlongGUI
real(c_double) :: deltat

NstepsGUI = nsteps
NcircGUI = Ncirc
NlongGUI = Nlong
deltat = DELTA_T

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen, res
character*(128) :: infile, outfile
integer :: i
logical :: ok, isopen

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

awp_0%is_open = .false.
awp_1%is_open = .false.
par_zig_init = .false.
logfile = 'mcell.log'
inquire(unit=nflog,OPENED=isopen)
if (.not.isopen) then
    open(nflog,file=logfile,status='replace')
endif
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
	write(nflog,*) 'call connecter'
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call Setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	res = 0
else
	call logger('=== Setup failed ===')
	res = 1
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
!call CheckAngles
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
else
	call logger('  === Execution failed ===')
endif
!write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine


end module
