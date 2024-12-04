program lowresmodel
implicit none

!variable declaration
real, parameter :: Lx = 6000, Ly = 2000                 !domain in km
real, parameter :: Dx = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)        !create grid points for array
integer :: i, j, t, latitude, degrees, nstep            !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, f, hsurf, lq   !array creation with internal variables for each point
real, dimension(Nx, Ny) :: hprev1, hprev2 ! Previous h
real, dimension(Nx, Ny) :: uprev1, uprev2 ! Previous u
real, dimension(Nx, Ny) :: vprev1, vprev2 ! Previous v

!Initial conditions, t = 0
 
!showcase model info
print *, "Grid Type: C"
print *, "Domain (km):", Lx, Ly
print *, "Grid Resolution (km):", Dx
print *, "X-Direction Grid Point #:", Nx
print *, "Y-Direction Grid Point #:", Ny

!calculate f
degrees = 45.0
f = 10e-04             	!f is going to be constant f = 10e-04

call initgrid(h, u, v, q, hsurf)


hprev1 = h
hprev2 = h
uprev1 = u
uprev2 = u
vprev1 = v
vprev2 = v

!create Adams-Bashforth time differencing loop
time = 0
nstep = 100

do t = 1, nstep
	if (t == 1) then 
		call firststep()
        if (t == 2) then
		call secondstep()
	else 
		call thirdstep()
	end if

!Update conditions
time = time + dT
h_prev2 = h_prev1
h_prev1 = h
u_prev2 = u_prev1
u_prev1 = u
v_prev2 = v_prev1
v_prev1 = v

end do 

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

subroutine firststep()
implicit none 
!

		
end subroutine firststep	
	                
			
subroutine secondstep()
end subroutine secondstep

subroutine thirdstep()
end subroutine thirdstep

subroutine math() !big subroutine to call all three main math nodes at the same time.
end subroutine math()

subroutine momentum()	!math for velocity (WIP)
end subroutine momentum

subroutine height()	!math for determining scale height
endsubroutine height




subroutine calculate_potential_vorticity(u, v, h, f, q, nx, ny)	!math for vorticity
implicit none
real, intent(in) :: u(nx,ny), v(nx,ny), h(nx,ny), f(nx,ny)
integer, intent(in) :: nx, ny
real, intent(out) :: q(nx,ny)
real :: zeta(nx,ny) !relative vorticity

!initialize to zero
q = 0
zeta = 0

!relative vorticity calculation
do j=2, ny-1
	do i=2, nx-1
		zeta(i,j) = (v(i+1,j)-v(i-1,j))/(dx) - (u(i,j+1)-u(i,j-1))/(dx)
	end do
end do

!potential vorticity calculations
do j=2, ny-1
	do i=2, nx-1
		if (h(i,j)>0.0) then !the >0 avoids divison by zero (just in case)
			q(i,j) = (zeta(i,j) + f(i,j)) / h(i,j)
		else
			q(i,j) = 0.0
		end if
	end do
end do
end subroutine calculate_potential_vorticity
