!-------------------------------------------------------------------------------------
! To determine corners of a cell that is in a sheet of such cells making up the wall of the mouse heart tube
! Cells exist (in general) left and right and above and below the focus cell (i,j)
!
!                           U o (i,j+1)
!                             |
!                             |
!                             |
!               (i-1,j)       |         (i+1,j)
!                 o-----------O-----------o
!                 L           |           R
!                             |
!                             |
!                             |
!                           D o (i,j-1) 
!
! We can define six planes:
!    Normal to the line OL, intersecting it at the mid-point
!    Normal to the line OR, intersecting it at the mid-point
!    Normal to the line OD, intersecting it at the mid-point
!    Normal to the line OU, intersecting it at the mid-point
! The last two planes are a bit more complicated to define.
! First we estimate the plane P0 that is the best fit to all five points: O, L, R, D, U
! The last two planes are parallel to this plane and displaced from O by +- w(z)/2
!
! Plane normal to a line N and passing through a point (xm,ym,zm):
! Any point (x,y,z) in the plane satisfies: V (x-xm,y-ym,z-zm) is normal to N: V.N = 0
! If the line N is OL, i.e. (xL-xO,yL-yO,zL-zO) then we get:
!     (x-xm)(xL-xO) + (y-ym)(yL-yO) + (z-zm)(zL-zO) = 0
! or  (x-xm)v_N(1) + (y-ym)v_N(2) + (z-zm)v_N(3)
! where
!     xm = (xL+xO)/2
!     ym = (yL+yO)/2
!     zm = (zL+zO)/2
!
!     x(xL-xO) + y(yL-yO) + z(zL-zO) = (xL+xO)(xL-xO)/2 + (yL+yO)(yL-yO)/2 + (zL+zO)(zL-zO)/2
!                                    = (xL^2 + yL^2 + zL^2)/2 - (xO^2 + yO^2 + zO^2)/2
! or  x.v_N(1) + y.v_N(2) + z.v_N(3) = xm.v_N(1) + ym.v_N(2) + zm.v_N(3)
! => a(k).x + b(k).y + c(k).z = d(k)
!    a(k) = v_N(1)
!    b(k) = v_N(2)
!    c(k) = v(N(3)
!    d(k) = xm.v_N(1) + ym.v_N(2) + zm.v_N(3)
! k=1 line OL
! k=2 line OR
! k=3 line OD
! k=4 line OU
!
! The general plane through O is:
!    P: a(x-xO) + b(y-yO) + c(z-zO) = 0 => ax + by + cz + d = 0 where d = -(a.x0 + b.y0 + c.z0)
! We want to choose a,b,c such that Q = dist(LP)^2 + dist(RP)^2 + dist(DP)^2 + dist(UP)^2
! is a minimum.
! The distance of p1 = (x1,y1,z1) from the plane: ax + by + cz + d = 0 is:
! |a.x1 + b.y1 + c.z1 + d|/sqrt(a^2 + b^2 + c^2)
! Therefore dist^2 = (a.x1 + b.y1 + c.z1 + d)^2/(a^2 + b^2 + c^2)
! and Q = { (a.xL + b.yL + c.zL + d)^2 + (a.xR + b.yR + c.zR + d)^2
!         + (a.xD + b.yD + c.zD + d)^2 + (a.xU + b.yU + c.zU + d)^2 }/((a^2 + b^2 + c^2)
! 
! Best fit plane through N points:
! math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
! Subtract out the centroid, form a 3xN matrix X out of the resulting coordinates and calculate its 
! singular value decomposition.  
! en.wikipedia.org/wiki/Singular_value_decomposition#Calculating_the_SVD
! The normal vector of the best-fitting plane is the left singular vector corresponding to the least singular value.


! Simple method: Just use three points.
!   either (i,j+1), (i-1,j-1) and (i+1,j-1)
!   or     (i,j-1), (i-1,j+1) and (i+1,j+1)
! or better, use mid-points p(1), p(2), p(3):
!   p(1) = (O+U)/2
!   p(2) = (O+L+D+LD)/4
!   p(3) = (O+R+D+RD)/4
! where LD = (i-1,j-1) and RD = (i+1,j-1)
!
! Vector normal to the plane: N = (p(2)-p(1))x(p(3)-p(1))
! This vector points outwards.  The two planes for i=5,6 are normal to n and pass through the two
! points offset from O along N:
! k=5 (xm,ym,zm) = O + (w/2)N/|N| is outside
! k=6 (xm,ym,zm) = O - (w/2)N/|N| is inside
! where w = the cell 'z' width.
! Then with v_N = N the plane is generated as for k=1,2,3,4
!
! (Need to decide how to treat cells at the edges: j=1 and j=Nlong
!-------------------------------------------------------------------------------------

module geometry
use global

implicit none

contains

!-------------------------------------------------------------------------------------
! Make adjacent cell vertices coincident.  
!                   inext    icirc
! Bottom vertices:  2   6    1   5   1
! Top vertices:     3   7    4   8   Nlong
! Middle vertices:  3   7    4   8   ilong
!                   2   6    1   5   ilong + 1
!                  out  in  out  in
!-------------------------------------------------------------------------------------
subroutine smoother
integer :: icirc, ilong, inext, k
real(REAL_KIND) :: v(4,3), vave(3), d(5), vside(3)
real(REAL_KIND) :: a, b, c, s

do ilong = 0,Nlong
	do icirc = 1,Ncirc
		inext = icirc - 1
		if (inext == 0) inext = Ncirc
		if (ilong == 0) then			! bottom vertices, 2 cells
			do k = 1,5,4
				v(1,:) = mcell(icirc,1)%vert(k,:)
				v(2,:) = mcell(inext,1)%vert(k+1,:)
				vave = (v(1,:) + v(2,:))/2
				mcell(icirc,1)%vert(k,:) = vave
				mcell(inext,1)%vert(k+1,:) = vave
			enddo
		elseif (ilong == Nlong) then	! top vertices, 2 cells
			do k = 4,8,4
				v(1,:) = mcell(icirc,Nlong)%vert(k,:)
				v(2,:) = mcell(inext,Nlong)%vert(k-1,:)
				vave = (v(1,:) + v(2,:))/2
				mcell(icirc,Nlong)%vert(k,:) = vave
				mcell(inext,Nlong)%vert(k-1,:) = vave
			enddo		
		else							! middle vertices, 4 cells
			do k = 4,8,4
				v(1,:) = mcell(icirc,ilong)%vert(k,:)
				v(2,:) = mcell(inext,ilong)%vert(k-1,:)
				v(3,:) = mcell(icirc,ilong+1)%vert(k-3,:)
				v(4,:) = mcell(inext,ilong+1)%vert(k-2,:)
				vave = (v(1,:) + v(2,:) + v(3,:) + v(4,:))/4
				mcell(icirc,ilong)%vert(k,:) = vave
				mcell(inext,ilong)%vert(k-1,:) = vave
				mcell(icirc,ilong+1)%vert(k-3,:) = vave
				mcell(inext,ilong+1)%vert(k-2,:) = vave
			enddo
		endif
	enddo	
enddo
! Compute area of internal faces (5,6,7,8)
! Triangle area A = sqrt(s(s-a)(s-b)(s-c)) where s = (a+b+c)/2
! Two triangles: (5,6,7) and (5,7,8)
do ilong = 1,Nlong
	do icirc = 1,Ncirc
		do k = 1,4
			if (k < 4) then
				vside = mcell(icirc,ilong)%vert(4+k+1,:) - mcell(icirc,ilong)%vert(4+k,:)
			else
				vside = mcell(icirc,ilong)%vert(4+1,:) - mcell(icirc,ilong)%vert(4+k,:)
			endif
			d(k) = sqrt(dot_product(vside,vside))
		enddo
		vside = mcell(icirc,ilong)%vert(7,:) - mcell(icirc,ilong)%vert(5,:)
		d(5) = sqrt(dot_product(vside,vside))
		! triangle (5,6,7)
		a = d(1)
		b = d(2)
		c = d(5)
		s = (a+b+c)/2
		mcell(icirc,ilong)%area = sqrt(s*(s-a)*(s-b)*(s-c))
		! triangle (5,7,8)
		a = d(5)
		b = d(3)
		c = d(4)
		s = (a+b+c)/2
		mcell(icirc,ilong)%area = mcell(icirc,ilong)%area + sqrt(s*(s-a)*(s-b)*(s-c))
	enddo
enddo
end subroutine

!-------------------------------------------------------------------------------------
! i = icirc, j = ilong
!-------------------------------------------------------------------------------------
subroutine makeVertices(i,j)
integer :: i,j
real(REAL_KIND) :: p_O(3), p_L(3), p_R(3), p_D(3), p_U(3), p_LD(3), p_RD(3)
real(REAL_KIND) :: p_M_L(3), p_M_R(3), p_M_D(3), p_M_U(3), p_M_O(3), p_M_I(3)
real(REAL_KIND) :: v_N(3), a(6), b(6), c(6), d(6), p_T(3,3), vx(3), vy(3), va, w
real(REAL_KIND) :: AV(3,3), rhs(3), sol(3)
integer :: k, row, kv(8,3)
logical :: dbug = .false.

p_O = mcell(i,j)%centre
if (i > 1) then
	p_L = mcell(i-1,j)%centre
else
	p_L = mcell(Ncirc,j)%centre
endif
if (i < Ncirc) then
	p_R = mcell(i+1,j)%centre
else
	p_R = mcell(1,j)%centre
endif
if (j > 1) then
	p_D = mcell(i,j-1)%centre	
endif
if (j < Nlong) then
	p_U = mcell(i,j+1)%centre
endif

k = 1	! line OL
p_M_L = (p_O + p_L)/2
v_N = p_L - p_O
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_L: ',k,p_M_L
d(k) = dot_product(p_M_L,v_N)
k = 2	! line OR
p_M_R = (p_O + p_R)/2
v_N = p_R - p_O
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_R: ',k,p_M_R
d(k) = dot_product(p_M_R,v_N)
k = 3	! line OD
if (j > 1) then
	p_M_D = (p_O + p_D)/2
	v_N = p_D - p_O
else
	p_M_D = 1.5*P_O  - 0.5*P_U
	v_N = P_O - p_U
endif
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_D: ',k,p_M_D
d(k) = dot_product(p_M_D,v_N)
k = 4	! line OU
if (j < Nlong) then
	p_M_U = (p_O + p_U)/2
	v_N = p_U - p_O
else
	p_M_U = 1.5*P_O  - 0.5*P_D
	v_N = P_O - p_D
endif
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_U: ',k,p_M_U
d(k) = dot_product(p_M_U,v_N)

vx = p_M_R - p_M_L
va = sqrt(dot_product(vx,vx))
vx = vx/va
vy = p_M_U - p_M_D
va = sqrt(dot_product(vy,vy))
vy = vy/va
if (dbug) write(*,'(a,3f8.3)') 'vx: ',vx
if (dbug) write(*,'(a,3f8.3)') 'vy: ',vy
if (dbug) write(*,'(a,f8.4)') 'vx.vy: ',dot_product(vx,vy)
call cross_product(vx,vy,v_N)
if (dbug) write(*,'(a,3f8.3)') 'v_N: ',v_N
va = sqrt(dot_product(v_N,v_N))
v_N = v_N/va
mcell(i,j)%vx = vx
mcell(i,j)%vy = vy
mcell(i,j)%vz = v_N
if (dbug) write(*,'(a,3f8.3)') 'v_N: ',v_N
w = mcell(i,j)%width(3)
if (dbug) write(*,'(a,f8.3)') 'w: ',w
k=5		! (xm,ym,zm) = O + (w/2)N/|N| is outside
p_M_O = p_O + (w/2)*v_N
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_O: ',k,p_M_O
d(k) = dot_product(p_M_O,v_N)
k=6		! (xm,ym,zm) = O - (w/2)N/|N| is inside
p_M_I = p_O - (w/2)*v_N
a(k) = v_N(1)
b(k) = v_N(2)
c(k) = v_N(3)
if (dbug) write(*,'(a,i4,3f8.4)') 'p_M_I: ',k,p_M_I
d(k) = dot_product(p_M_I,v_N)

if (dbug) then
	write(*,*) 'i,j: ',i,j
	write(*,'(a,3f8.4)') 'p_O: ',p_O
	write(*,'(a,3f8.4)') 'p_L: ',p_L
	write(*,'(a,3f8.4)') 'p_R: ',p_R
	write(*,'(a,3f8.4)') 'p_D: ',p_D
	write(*,'(a,3f8.4)') 'p_U: ',p_U
	write(*,'(a,3f8.4)') 'p_LD: ',p_LD
	write(*,'(a,3f8.4)') 'p_RD: ',p_RD
	do k = 1,6
		write(*,'(a,i2,4f8.4)') 'k,a,b,c,d: ', k,a(k),b(k),c(k),d(k)
	enddo
endif

! The next step is to deduce the 8 intersection points of the 6 planes: L, R, U, D, O, I
! Number the vertices anticlockwise starting from outside.
! The outside 'O' vertices are: OLD, ORD, ORU, OLU
! The inside 'I' vertices are: ILD, IRD, IRU, ILU
! The cases are:
!   vertex  planes 
!      1    5 1 3 OLD
!      2    5 2 3 ORD
!      3    5 2 4 ORU
!      4    5 1 4 OLU
!      5    6 1 3 ILD
!      6    6 2 3 IRD
!      7    6 2 4 IRU
!      8    6 1 4 ILU
kv(1,:) = [5, 1, 3]
kv(2,:) = [5, 2, 3]
kv(3,:) = [5, 2, 4]
kv(4,:) = [5, 1, 4]
kv(5,:) = [6, 1, 3]
kv(6,:) = [6, 2, 3]
kv(7,:) = [6, 2, 4]
kv(8,:) = [6, 1, 4]
! Now for each vertex need to solve 3 eqtns in 3 unknowns
! In matrix form:
!  | a1 b1 c1 | | x | = | d1 |
!  | a2 b2 c2 | | y | = | d2 |
!  | a3 b3 c3 | | z | = | d3 |
do k = 1,8
!	write(*,*) 'AV: ',k
	do row = 1,3
		AV(row,1) = a(kv(k,row))
		AV(row,2) = b(kv(k,row))
		AV(row,3) = c(kv(k,row))
		rhs(row) = d(kv(k,row))
!		write(*,'(2i4,4f8.3)') row,kv(k,row),AV(row,:),rhs(row)
	enddo
	call solveNxN(3,AV,rhs,sol)
	mcell(i,j)%vert(k,:) = sol
	if (dbug) write(*,'(a,i2,3f8.3)') 'vertex: ',k,sol
enddo
call smoother
end subroutine

!-------------------------------------------------------------------------------------
! Solve 3x3 linear system of eqtns: A.sol = rhs
! Cramer's rule:
! D = det(A), Dx = det(Ax), Dy = det(Ay), Dz = det(Az)
! where Ax has column 1 replaced by rhs, Ay has column 2 replaced by rhs, Az has column 3 replaced by rhs
! Then solution is: x = Dx/D, y = Dy/D, z = Dz/D
!-------------------------------------------------------------------------------------
subroutine solveNxN(n,A,rhs,sol)
integer :: n
real(REAL_KIND) :: A(n,n), rhs(n), sol(n)
real(REAL_KIND) :: AA(n,n), detA, detAA, eps
integer :: i

eps = 1.0e-10
if (n == 3) then
	detA = det3(A)
else
	detA = dmgt(eps,n,A)
endif
!write(*,*) 'detA: ',detA
do i = 1,n
	AA = A
	AA(:,i) = rhs	! Ax, Ay, Az
	if (n == 3) then
		detAA = det3(AA)
	else
		detAA = dmgt(eps,n,AA)
	endif
!	write(*,*) 'detAA: ',i,detAA
	sol(i) = detAA/detA
enddo
end subroutine

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
real(REAL_KIND) function det3(A)
real(REAL_KIND) :: A(3,3)

det3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
      -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
      +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function


end module
