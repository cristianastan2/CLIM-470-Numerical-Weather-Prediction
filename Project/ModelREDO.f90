program lowresmodel
implicit none

!variable declaration
integer, parameter :: Lx = 6000, Ly = 2000                 !domain in km
integer, parameter :: Dx = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)        !create grid points for array
integer :: i, j, t, n, irec                !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, hsurf, lq, zeta, ght, ken   !array creation with internal variables for each point
real, dimension(Nx, Ny, 3) :: hu0, hv0, hq0, us0, vs0, h0 ! Variables for time differencing
real, dimension(Nx, Ny) :: hu1, hu2, hu3, hv1, hv2, hv3, us1, us2, us3, vs1, vs2, vs3
real, dimension(Nx, Ny, 3) :: alp0, bet0, gam0, del0, eps0, ken0, ght0, q0, z0, phi0
real, dimension(Nx, Ny, 3) :: u0, v0

!Variables for further time-stepping
real, dimension(Nx, Ny) :: alp1, alp2, alp3, bet1, bet2, bet3
real, dimension(Nx, Ny) :: gam1, gam2, gam3, del1, del2, del3
real, dimension(Nx, Ny) :: eps1, eps2, eps3, ken1, ken2, ken3
real, dimension(Nx, Ny) :: ght1, ght2, ght3
real, parameter :: g = 9.81, f = 1.0e-4 !Constants

!showcase model info
print *, "Grid Type: C"
print *, "Domain (km):", Lx, Ly
print *, "Grid Resolution (km):", Dx
print *, "X-Direction Grid Point #:", Nx
print *, "Y-Direction Grid Point #:", Ny


! Set initial conditions
u(1:Nx,1:Ny) = 20 ! uniform zonal wind
u(1:Nx,1) = 0.0
v(1:Nx, Ny) = 0.0 ! rigid computational boundary
v(1:Nx,2:Ny-1) = 0.1

!Just for writing out
u0(:,:,1) = u
v0(:,:,1) = v

call initgrid(hsurf,Dx,Nx,Ny)

print*, "Calculated hsurf"

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

do i = 1, Nx
	z0(i,1,1)=0.
	z0(i,Ny,1)=0.
end do

do j=2,Ny-1
	z0(1,j,1)=(u(1,j-1)-u(1,j+1)+v(2,j)-v(Nx,j))/Dx
	z0(Nx,j,1)=(u(Nx,j-1)-u(1,j+1)+v(1,j)-v(Nx-1,j))/Dx
end do

do i=2,Nx-1
	do j=2,Ny-1
		z0(i,j,1)=(u(i,j-1)-u(i,j+1)+v(i+1,j)-v(i-1,j))/Dx
	end do
end do

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

! We can print out the initial conditions


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

u0(:,:,n) = u
v0(:,:,n) = v

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

!Update us0 and vs0 for the next time step
us0(:,:,n)= hu0(:,:,n)*u(:,:)
vs0(:,:,n) = hv0(:,:,n)*v(:,:)

! Update z0 for the next time step
	do j=2,Ny-1

	z0(1,j,n)=(u(1,j-1)-u(1,j+1)+v(2,j)-v(Nx,j))/Dx
	z0(Nx,j,n)=(u(Nx,j-1)-u(1,j+1)+v(1,j)-v(Nx-1,j))/Dx

	end do
	
! Update h10 for the next time step
do j = 2, Ny-1
   do i = 2, Nx
      hq0(i, j, n) = (h(i, j) + h(i-1, j) + h(i-1, j-1) + h(i, j-1)) / 4.0
   end do
end do

! Update q0 for the next time step 

q0(:,:,n) = (f+z0(:,:,n))/hq0(:,:,n)

! update ght0 for the next time step
	do i = 1, Nx
	ght0(i,:,n) = g*(h(i,:)+hsurf(i,:))
	end do
	
! Update ken0 for the next time step
ken0(:,:,n)=(u(:,:)*u(:,:) + v(:,:)*v(:,:))/2


end do ! Time loop 

!Update previous values prior to next time step
us1(:,:) = us0(:,:,1)
us2(:,:) = us0(:,:,2)
us3(:,:) = us0(:,:,3)

hu1(:,:) = h0(:,:,1)
hu2(:,:) = h0(:,:,2)
hu3(:,:) = h0(:,:,3)

hv1(:,:) = hv0(:,:,1)
hv2(:,:)= hv0(:,:,2)
hv3(:,:) = hv0(:,:,3)

alp1(:,:) = alp0(:,:,1)
alp2(:,:) = alp0(:,:,2)
alp3(:,:)= alp0(:,:,3)

bet1(:,:) = bet0(:,:,1)
bet2(:,:) = bet0(:,:,2)
bet3(:,:)= bet0(:,:,3)

gam1(:,:) = gam0(:,:,1)
gam2(:,:) = gam0(:,:,2)
gam3(:,:)= gam0(:,:,3)

del1(:,:) = del0(:,:,1)
del2(:,:) = del0(:,:,2)
del3(:,:)= del0(:,:,3)

eps1(:,:) = eps0(:,:,1)
eps2(:,:) = eps0(:,:,2)
eps3(:,:)= eps0(:,:,3)

ken1(:,:) = ken0(:,:,1)
ken2(:,:) = ken0(:,:,2)
ken3(:,:)= ken0(:,:,3)

ght1(:,:) = ght0(:,:,1)
ght2(:,:) = ght0(:,:,2)
ght3(:,:)= ght0(:,:,3)



open(unit=10, file='topo.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*(Nx-2)*(Ny-2))  !scary thing
write(10, rec=1) hsurf(2:Nx-1,2:Ny-1)

open(unit=15, file='geopot.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*(Nx-2)*(Ny-2))
irec=0
do n=1,3
	irec=irec+1
	write(15, rec=irec) ght0(2:Nx-1,2:Ny-1,n)
end do

open(unit=20, file='uwind.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*(Nx-2)*(Ny-2))
irec=0 
do n=1,3
        irec=irec+1
        write(20, rec=irec) u0(2:Nx-1,2:Ny-1,n)
end do

open(unit=25, file='vwind.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*(Nx-2)*(Ny-2))
irec=0 
do n=1,3
        irec=irec+1
        write(25, rec=irec) v0(2:Nx-1,2:Ny-1,n)
end do

end program lowresmodel



subroutine initgrid(hsurf, Dx, Nx, Ny)	!creates grid, ridge, and topography.
    implicit none
  
    integer, parameter :: Lx = 6000, Ly = 2000
    integer, intent(in) :: Dx, Nx, Ny 
   

    real, dimension(Nx,Ny), intent(out) :: hsurf
    integer :: i, j
  
    print*, "hsurf =", hsurf
    hsurf(:,:)=0.0
    print*, "Dx =", Dx 

    !grid resolution thing (topography)

    hsurf((Nx)/2 + 1,:) = 2000
    if (Dx==250) then
        hsurf((Nx)/2,:) = 1000
        hsurf((Nx)/2 + 2,:) = 1000
    endif
    if (Dx==125) then
        hsurf((Nx)/2 - 1,:) = 500
        hsurf((Nx)/2 + 3,:) = 500
    endif
	

    end subroutine initgrid
