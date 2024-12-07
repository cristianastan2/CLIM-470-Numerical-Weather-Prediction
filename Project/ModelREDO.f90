program lowresmodel
implicit none

!variable declaration
real, parameter :: Lx = 6000, Ly = 2000                 !domain in km
real, parameter :: Dx = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)        !create grid points for array
integer :: i, j, t, n                !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, hsurf, lq, zeta, ght, ken   !array creation with internal variables for each point
real, dimension(Nx, Ny, 3) :: hu0, hv0, hq0, us0, vs0, h0 ! Variables for time differencing
real, dimension(Nx, Ny) :: hu1, hu2, hu3, hv1, hv2, hv3, us1, us2, us3, vs1, vs2, vs3
real, dimension(Nx, Ny, 3) :: alp0, bet0, gam0, del0, eps0, ken0, ght0, q0, z0, phi0
real, parameter :: g = 9.81, f = 1.0e-4 !Constants



!showcase model info
print *, "Grid Type: C"
print *, "Domain (km):", Lx, Ly
print *, "Grid Resolution (km):", Dx
print *, "X-Direction Grid Point #:", Nx
print *, "Y-Direction Grid Point #:", Ny


! Set initial conditions
u(1:Nx,1:Ny) = 20 ! uniform zonal wind
v(1:Nx,1) = 0.0
v(1:Nx, Ny) = 0 ! rigid computational boundary
v(1:Nx,2:Ny-1) = 0.1

do i = 1,Nx
	h(i,:) = 5000-hsurf(i, :) ! horizontal free surface
end do

!Initialize hu0 and hv0 arrays
hu0(1,:,1) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,1) = (h(Nx-1,:) + h(1,:))/2.0

do i = 2, Nx-1
	hu0(i,:,1) = (h(i-1,:) + h(i+1,:))/2.0
end do

hv0(:,1,1) = (h(:,Ny) + h(:,2))/2.0
hv0(:,Ny,1) = (h(:,Ny-1) + h(:,1))/2.0

do j = 2, Ny-1
	hv0(:,j,1) = (h(:,j-1) + h(:,j+1))/2.0
end do

! Time differencing for u and v velocities
us0(:,:,1)= hu0(:,:,1)*u(:,:)
vs0(:,:,1) = hv0(:,:,1)*v(:,:)

!Initial conditions, t = 0
hq0(1,1,1)=(h(1,1) + h(Nx,1))/2.0
hq0(1,Ny,1)=(h(1,Ny-1) + h(Nx,Ny-2))/2.0

do i = 2,Nx
	hq0(i,1,1)=(h(i,1) + h(i-1,1))/2.0
	hq0(i,Ny,1)=(h(i,Ny-1) + h(i-1,Ny-1))/2.0
end do

do j = 2,Ny-1
	hq0(1,j,1)=(h(1,j) + h(i-1,j) + h(Nx,j) + h(Nx-1,j-1))/4.0
end do

do j = 2,Ny-1
	do i = 2,Nx
		hq0(i,j,1)=(h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
	end do
end do

do j = 1,Ny
	do i = 1,Nx
		q0(i,j,1)=(z0(i,j,1) + f)/hq0(i,j,1)
	end do
end do

do i = 1,Nx
	do j = 1, Ny
		ght0(i,j,1) = g*(hsurf(i,j) + h(i,j))
	end do
end do

do i = 1, Nx
	ken0(i,Ny,1)=0.
end do

do j = 1,Ny-1
	ken0(Nx,j,1)=(u(Nx-1,j)**2 + u(1,j)**2 + v(Nx-1,j)**2 + v(1,Ny)**2)/4.
end do

do i = 1, Nx-1
	do j = 1, Ny-1
		ken0(i,j,1)=(u(i,j)**2 + u(i+1,j)**2 + v(i,j)**2 + v(1,j+1)**2)/4.
	end do
end do

!!! We can print out the initial conditions


! Time-stepping loop
do n = 2, 3
	do i = 2, Nx-1
		do j = 2, Ny-1
		
			h(i,j) = h(i,j) - dT * (us0(i+1,j+1,n-1) - us0(i,j+1,n-1) + vs0(i+1,j+1,n-1) - vs0(i+1,j,n-1))/Dx
			
			u(i,j) = u(i,j) + dT * (alp0(i,j+1,n-1) * vs0(i,j,n-1) + bet0(i,j+1,n-1) * vs0(i-1,j+1,n-1) + &
                        gam0(i,j+1,n-1) * vs0(i-1,j,n-1) + del0(i,j+1,n-1) * vs0(i+1,j,n-1) - eps0(i+1,j+1,n-1) * &
			us0(i+1,j+1,n-1) + eps0(i-1,j+1,n-1) * us0(i-1,j+1,n-1) - (ken0(i+1,j+1,n-1) + &
			ght0(i+1,j+1,n-1) - ken0(i-1,j+1,n-1) - ght0(i-1,j+1,n-1)) / Dx)
			
                        v(i,j) = v(i,j) - dT * (gam0(i+1,j+1,n-1) * us0(i+1,j+1,n-1) + del0(i,j+1,n-1) * us0(i,j+1,n-1) + &
                        alp0(i,j-1, n-1) * us0(i,j-1,n-1) + bet0(i+1,j-1,n-1) * us0(i+1,j-1,n-1) + phi0(i+1,j+1,n-1) * &
			vs0(i+1,j+1,n-1) - phi0(i+1,j-1,n-1) * vs0(i+1,j-1,n-1) - (ken0(i+1,j+1,n-1) + &
			ght0(i+1,j+1,n-1) - ken0(i+1,j-1,n-1) - phi0(i+1,j-1,n-1)) / Dx)
			
			! Update zeta (z)
			zeta(i,j) = (u(i,j-1) - u(i,j+1) + v(i+1,j) - v(i-1,j)) / Dx
			
			! Update Potential Vorticity
			q(i,j) = (f + zeta(i,j)) / h(i,j)
			
			! Compute Greek letters (alp0, bet0, gam0, del0, eps0)
			alp0(i, j, n) = (1.0 / 24.0) * (2.0 * q(i+1, j+1) + q(i, j+1) + 2.0 * q(i, j) + q(i+1, j))
			bet0(i, j, n) = (1.0 / 24.0) * (q(i, j+1) + 2.0 * q(i-1, j+1) + q(i-1, j) + 2.0 * q(i, j))
			gam0(i, j, n) = (1.0 / 24.0) * (2.0 * q(i, j+1) + q(i-1, j+1) + 2.0 * q(i-1, j) + q(i, j))
			del0(i, j, n) = (1.0 / 24.0) * (q(i+1, j+1) + 2.0 * q(i, j+1) + q(i, j) + 2.0 * q(i+1, j))
			eps0(i+1, j+1, n) = (1.0 / 24.0) * (q(i+1, j+1) + q(i, j+1) - q(i, j+1) - q(i+1, j))
	  
		end do
	end do


! Update hu0 and hv0 for next time step 
hu0(1,:,n) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,n) = (h(Nx-1,:) + h(1,:))/2.0

	do i = 2, Nx-1
	
	hu0(i,:,n) = (h(i-1,:) + h(i+1,:))/2.0
	
	end do
	
hv0(:,1,n) = (h(:,Ny) + h(:,2))/2.0
hv0(:,Ny,n) = (h(:,Ny-1) + h(:,1))/2.0

	do j = 2, Ny-1
	hv0(:,j,n) = ( h(:,j-1) + h(:,j+1) )/2.0
	end do
	
us0(:,:,n)= hu0(:,:,n)*u(:,:)
vs0(:,:,n) = hv0(:,:,n)*v(:,:)

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
