! Interface to min_dist() in ellipsoid_dist.dll

interface

	subroutine min_dist(a1,b1,centre1,orient1,a2,b2,centre2,orient2,tol,s1,s2,r1,r2,d,res) BIND(C,NAME='min_dist')
	use :: iso_c_binding
	real(c_double),value :: a1, b1, a2, b2, tol
	real(c_double) :: centre1(*), orient1(*), centre2(*), orient2(*)
	real(c_double) :: s1, s2, r1, r2, d
	integer(c_int) :: res
	end subroutine

end interface
