module cellml
use global
use, intrinsic :: iso_c_binding
use csim_abm_mp

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CellML_load_file(file_array, buflen, nvar, ncomp) BIND(C,NAME='CellML_load_file')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_load_file
character(c_char),dimension(*) :: file_array
integer(c_int),value :: buflen
integer(c_int) :: nvar, ncomp
integer :: k

cellmlfile = ' '
do k = 1,buflen
	cellmlfile(k:k) = file_array(k)
enddo
write(*,*) cellmlfile
call loadCellML
nvar = nvariables
ncomp = ncomponents
write(nflog,*) 'Loaded CellML model successfully'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine loadCellML
character(64) :: ID, comp_name
logical :: hit
integer :: nlen, k, i, clen, ksim

write(nflog,'(a)') 'loadCellML: cellmlfile:'
write(nflog,'(a)') cellmlfile
call freeCellmlSimulators
ksim = 0
if (allocated(IDlist)) deallocate(IDlist)
if (allocated(csim_state)) deallocate(csim_state)
if (allocated(csim_savestate)) deallocate(csim_savestate)
if (allocated(component)) deallocate(component)
!call setupCellml(cellmlfile)
call setupCellmlSimulator(ksim,cellmlfile)
call getNvariables(ksim,nvariables)
write(nflog,*) 'nvariables: ',nvariables
allocate(IDlist(0:nvariables-1))
allocate(csim_state(0:nvariables-1))
allocate(csim_savestate(0:nvariables-1))
allocate(component(0:nvariables-1))

ncomponents = 0
do k = 0,nvariables-1
	call getID(ksim,k,ID,clen)
!	IDlist(k)%name = ID(1:clen)
	i = index(ID,'.')
	comp_name = ID(1:i-1)
	IDlist(k)%comp_name = comp_name
	IDlist(k)%var_name = ID(i+1:clen)
	if (ncomponents == 0) then
		component(ncomponents) = comp_name
		IDlist(k)%comp_index = ncomponents
		ncomponents = ncomponents + 1
	else
		hit = .false.
		do i = 0,ncomponents-1
			if (comp_name == component(i)) then
				hit = .true.
				IDlist(k)%comp_index = i
				exit
			endif
		enddo
		if (.not.hit) then
			component(ncomponents) = comp_name
			IDlist(k)%comp_index = ncomponents
			ncomponents = ncomponents + 1
		endif
	endif		
enddo
call getState(ksim,csim_state)
write(nflog,*) 'Components: #: ',ncomponents
write(nflog,'(a)') component(0:ncomponents-1)
write(nflog,*) 'Variables: '
do k = 0,nvariables-1
	write(nflog,'(i4,2x,2a30,i3,f10.6)') k,IDlist(k)%comp_name(1:30),IDlist(k)%var_name(1:30),IDlist(k)%comp_index,csim_state(k)
enddo
end subroutine

!--------------------------------------------------------------------------------
! Get component name and number of variables
!--------------------------------------------------------------------------------
subroutine CellML_get_component_info(icomp,component_name, nvars) BIND(C,name='CellML_get_component_info')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_component_info
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp
integer(c_int) :: nvars
character(c_char),dimension(1) :: component_name
integer :: nlen, i

nlen = index(component(icomp),' ')
do i = 1,nlen-1
	component_name(i) = component(icomp)(i:i)
enddo
component_name(nlen) = char(0)
nvars = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		nvars = nvars + 1
	endif
enddo
write(nflog,*) icomp,nvars,component(icomp)(1:nlen-1)
end subroutine

!--------------------------------------------------------------------------------
! Get variable name
!--------------------------------------------------------------------------------
subroutine CellML_get_variable_name(icomp,ivar,variable_name) BIND(C,name='CellML_get_variable_name')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_variable_name
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
character(c_char),dimension(1) :: variable_name
character(64) :: varname
integer :: nv, i, k, nlen

nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			nlen = index(IDlist(i)%var_name,' ')
			do k = 1,nlen-1
				variable_name(k) = IDlist(i)%var_name(k:k)
			enddo
			variable_name(nlen) = char(0)
			write(nflog,*) icomp,ivar,IDlist(i)%var_name(1:nlen-1)
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Get variable value
!--------------------------------------------------------------------------------
subroutine CellML_get_variable_value(icomp,ivar,variable_value) BIND(C,name='CellML_get_variable_value')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_variable_value
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
real(c_double) :: variable_value
integer :: nv, i

write(nflog,*) 'CellML_get_variable_value: ',icomp,ivar
call getState(0,csim_state)
!write(nflog,'(10f7.3)') csim_state
nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			variable_value = csim_state(i)
			write(nflog,*) icomp,ivar,variable_value
!			if (icomp == 2 .and. ivar == 5) stop
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Set variable value
!--------------------------------------------------------------------------------
subroutine CellML_set_variable_value(icomp,ivar,variable_value) BIND(C,name='CellML_set_variable_value')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_set_variable_value
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
real(c_double),value :: variable_value
integer :: nv, i

nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			call setStateValue(0,i,variable_value)
			write(nflog,*) icomp,ivar,variable_value
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine


end module