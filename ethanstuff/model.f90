program lowresmodel
implicit none

!variable declaration
real, parameter :: Lx = 6000, Ly = 2000                 !domain in km
real, parameter :: Dx = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)        !create grid points for array
integer :: i, j, t, latitude, degrees, nstep            !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, f, hsurf, lh, lu, lv, lq   !array creation with internal variables for each point

!showcase model info
print *, "Grid Type: C"
print *, "Domain (km):", Lx, Ly
print *, "Grid Resolution (km):", Dx
print *, "X-Direction Grid Point #:", Nx
print *, "Y-Direction Grid Point #:", Ny

!create Adams-Bashforth time differencing loop
time = 0
nstep = 100

do t = 1, nstep
	if (t == 1) then 
		call firststep()
	else if (t == 2) then
		call secondstep()
	else 
		call thirdstep()
	end if
end do

call initgrid(h, u, v, q, hsurf)
open(unit=10, file='topo.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*Nx*Ny)  !scary thing
write(10, rec=1) hsurf


end program lowresmodel




subroutine initgrid(h, u, v, q, hsurf)	!creates grid, ridge, and topography.
implicit none

real, parameter :: Lx = 6000, Ly = 2000
real, parameter :: Dx = 500
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)

real, dimension(Nx,Ny), intent(out) :: h, u, v, q, hsurf
integer :: i, j

hsurf(:,:)=0

!grid resolution thing (topography)
if (Dx==500) then 
	hsurf(Nx/2,:)=2000

elseif (Dx==250) then
	hsurf(Nx/2-1,:)=1000
	hsurf(Nx/2,:)=2000
	hsurf(Nx/2+1,:)=1000

elseif (Dx==125) then
	hsurf(Nx/2-2,:)=500
	hsurf(Nx/2-1,:)=1000
	hsurf(Nx/2,:)=2000
	hsurf(Nx/2+1,:)=1000
	hsurf(Nx/2+2,:)=500
endif
	
!initial conditions set to 0
h = 0
u = 20
v = 0 ! needs to be changed
q = 0

end subroutine initgrid

subroutine firststep()		!(ethan) - My thoughts with the steps as subroutines is that we can use the same forumlas each time around
end subroutine firststep	! 	   It seems as though realisitically the only thing changing is the variables, and the steps just need
				!	   to calculate initial than the "next step" and then subtract them, so if we leave the subroutines down there and
				!          focus on these ones being the "complex" stuff with just calling the math we should be fine.
subroutine secondstep()
end subroutine secondstep

subroutine thirdstep()
end subroutine thirdstep

subroutine math() !big subroutine to call all three main math nodes at the same time.
end subroutine math

subroutine momentum()	!math for velocity (WIP)
end subroutine momentum

subroutine height()	!math for determining scale height
endsubroutine height




subroutine calculate_potential_vorticity(u, v, h, f, q, nx, ny, dx)	!math for vorticity
implicit none
real, intent(in) :: u(nx,ny), v(nx,ny), h(nx,ny), f(nx,ny), dx
integer, intent(inout) :: nx, ny
real, intent(out) :: q(nx,ny)
real :: zeta(nx,ny) !relative vorticity

!initialize to zero
q = 0
zeta = 0

!relative vorticity calculation
do ny=2, 1
	do nx=2, 1
		zeta(nx,ny) = (v(nx+1,ny)-v(nx-1,ny))/(dx) - (u(nx,ny+1)-u(nx,ny-1))/(dx)
	end do
end do

!potential vorticity calculations
do ny=2, 1
	do nx=2, 1
		if (h(nx,ny)>0.0) then !the >0 avoids divison by zero (just in case)
			q(nx,ny) = (zeta(nx,ny) + f(nx,ny)) / h(nx,ny)
		else
			q(nx,ny) = 0.0
		end if
	end do
end do
end subroutine calculate_potential_vorticity
