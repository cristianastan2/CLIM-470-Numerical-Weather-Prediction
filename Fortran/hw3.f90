program main
	implicit none
	real, parameter :: omega = 7.292e-5, pi = 3.14159, Rd = 287.0
	real :: lat, f, v, T, p1, p2, dx, rho, dpdx, dudx, u, dt, time
	integer :: n, nsteps
	
	lat = 35.0*pi/180.0
	f = 2*omega*sin(lat)
	v = 5.0
	T = 285.0
	p1 =1010.0e2    !in pascals, west
	p2 = 1006.0e2   !in pascals, east
	dx = 5000.0     !meters
	rho = 0.5*(p1+p2)/(Rd*T)
	dpdx = (p2-p1)/dx
	dudx = (4.0-3.0)/dx
	u = 0.5*(3.0+4.0)
	dt = 1.0
	nsteps = 3600

	open(10, file='u.dat')
	write(10,*) 0.0, u
	do n = 1, nsteps
		u = u + dt*(f*v - dpdx/rho - u*dudx)
		time = n*dt
		write(10,*) time, u
	end do
	close(10)
	write(*,*) 'u after 1 hour = ', u
end program main
