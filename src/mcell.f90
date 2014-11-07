! ABM for mouse heart-tube cells, off-lattice 
module mcell_mod
use global
use geometry
!use newton

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
real(REAL_KIND) :: hours
integer :: nb0, nt_anim
character*(256) :: cellmlfile

open(nfin,file=inputfile)
!read(nfin,*) nb0
!read(nfin,*) igrow
!read(nfin,*) t_div_median
!read(nfin,*) t_div_shape
read(nfin,*) Ncirc
read(nfin,*) Nlong
read(nfin,*) hours
read(nfin,*) DELTA_T
read(nfin,*) dx_init
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
read(nfin,*) cellmlfile
close(nfin)

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
if (allocated(Jac)) then
	deallocate(Jac)
endif
allocate(mcell(Ncirc,Nlong))
allocate(mcell0(Ncirc,Nlong))
allocate(vbase(Ncirc,3))
!allocate(v1(Ncirc,3))
if (Ndim == 2) then
	allocate(v2D(2*Ncirc*Nlong))
	allocate(Jac(2*Ncirc*Nlong,2*Ncirc*Nlong))
else
	allocate(v3D(3*Ncirc*Nlong))
	allocate(Jac(3*Ncirc*Nlong,3*Ncirc*Nlong))
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
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! xcur is the array containing the current x value, fcur 
! is f(xcur) on output, rpar and ipar are, respectively, real 
! and integer parameter/work arrays for use by the subroutine,
! and itrmf is an integer termination flag.  The meaning of
! itrmf is as follows:
!   0 => normal termination; desired function value calculated.
!   1 => failure to produce f(xcur).
!-----------------------------------------------------------------------------------------
subroutine fnit(n, xcur, fcur, rpar, ipar, itrmf)
integer :: n, ipar(*), itrmf
real(REAL_KIND) :: xcur(*), fcur(*), rpar(*)
integer :: icirc, ilong, nd, k, kd, nvars, i
real(REAL_KIND) :: F(3), fmax(3)

nd = Ndim
nvars = nd*Ncirc*Nlong
if (n /= nvars) then
	write(*,*) 'ERROR: fnit: n != nvars: ',n,nvars
	stop
endif
nd = Ndim
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		call getForce(nd,xcur,icirc,ilong,F)
		do kd = 1,nd
			k = (ilong-1)*nd*Ncirc + (icirc-1)*nd + kd 
			fcur(k) = F(kd)
		enddo
	enddo
enddo
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
! Dummy
!-----------------------------------------------------------------------------------------
subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
integer :: n, ijob, ipar(*), itrmjv
real(REAL_KIND) :: xcur(*), fcur(*), v, z, rpar(*)

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
!write(*,*) 'v'
!write(*,'(10f7.4)') v(1:n)
!write(*,*) 'f'
!write(*,'(10f7.4)') rwork(1:n)
!deallocate(rwork)
!write(nflog,*) 'nitsolved: iterm: ',iterm
end subroutine

!--------------------------------------------------------------------------------
! Programmed test case for cell shape change.
! Case 1:
! width(2) at x = Rinitial  -> 1.5*width0
!          at x = -Rinitial -> 0.5*width0
! as istep -> nstep_var
!--------------------------------------------------------------------------------
subroutine updateWidths
integer :: ilong, icirc, ivar
real(REAL_KIND) :: x, y, z, ymid, tfactor, xfactor, yfactorx, yfactorz, zfactor, deltay, deltaz, cyz
real(REAL_KIND) :: deltay_max = 0.6
real(REAL_KIND) :: deltaz_max = 0.35
real(REAL_KIND) :: t_start = 60		! time until bending start (min)
real(REAL_KIND) :: t_var = 20*60	! duration of bending (min)
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
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
integer :: i, res
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
	stop
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
