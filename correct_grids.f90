program lowresmodel
implicit none

!variable declaration
integer, parameter :: Lx = 6000, Ly = 2000                 !domain in km
integer, parameter :: Dx = 500, Dxm = 500000                   !grid resolution 
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx) + 3, Ny = (Ly/Dx) + 3        !create grid points for array
integer :: i, j, t, n                !define other potentially important variables
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

call initgrid(hsurf, Dx, Nx, Ny)

! Set initial conditions 
u(:,:) = 0.0 
u(1:Nx,2:Ny-1) = 20 ! uniform zonal wind
v(:,:) = 0.0 
v(1:Nx,2:Ny-1) = 0.1

u0(:,:,1) = u(:,:)
v0(:,:,1) = v(:,:)

h(:,:) = 0.0 
do i = 1,Nx
	h(i,2:Ny-1) = 5000-hsurf(i, 2:Ny-1) ! horizontal free surface
end do

h0(:,:,1) = h(:,:)

hu0(:,:,1) = 0.0
hv0(:,:,1) = 0.0
us0(:,:,1) = 0.0 
vs0(:,:,1) = 0.0 
z0(:,:,1) = 0.0
hq0(:,:,1) = 0.0 
ght0(:,:,1) = 0.0
ken0(:,:,1) = 0.0
q0(:,:,1) = 0.0 

do i = 2, Nx-1
    do j = 2, Ny-1
        hu0(i,j,1) = (h(i-1,j) + h(i+1,j))/2.0
        hv0(i,j,1) = (h(i,j-1) + h(i,j+1))/2.0
        z0(i,j,1)=(u(i,j-1)-u(i,j+1)+v(i+1,j)-v(i-1,j))/Dxm
        hq0(i,j,1)=(h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
        ght0(i,j,1) = g*(hsurf(i,j) + h(i,j))
        ken0(i,j,1)=(u(i,j)**2 + u(i+1,j)**2 + v(i,j)**2 + v(1,j+1)**2)/4.
    end do
end do 

us0(2:Nx-1,2:Ny-1,1)= hu0(2:Nx-1,2:Ny-1,1)*u(2:Nx-1,2:Ny-1)
vs0(2:Nx-1,2:Ny-1,1)= hv0(2:Nx-1,2:Ny-1,1)*v(2:Nx-1,2:Ny-1)
      

u(1,:) = u(Nx-1,:)
u(Nx,:) = u(2,:)

v(1,:) = v(Nx-1,:)
v(Nx,:) = v(2,:)

hu0(1,:,1) = hu0(Nx-1,:,1)
hu0(Nx,:,1) = hu0(2,:,1)

hv0(1,:,1) = hv0(Nx-1,:,1)
hv0(Nx,:,1) = hv0(2,:,1)

z0(1,:,1) = z0(Nx-1,:,1) 
z0(Nx,:,1) = z0(2,:,1) 

hq0(1,:,1) = hq0(Nx-1,:,1) 
hq0(Nx,:,1) = hq0(2,:,1) 

ght0(1,:,1) = ght0(Nx-1,:,1) 
ght0(Nx,:,1) = ght0(2,:,1) 

ght0(1,:,1) = ght0(Nx-1,:,1) 
ght0(Nx,:,1) = ght0(2,:,1) 

! We can print out the initial conditions


! Time-stepping loop
do n = 2, 3
	do i = 2, Nx-1
		do j = 2, Ny-1
		
			h0(i,j,n) = h0(i,j,n-1) - dT * (us0(i+1,j,n-1) - us0(i,j,n-1) + vs0(i,j+1,n-1) - vs0(i,j,n-1))/Dxm
			
			! Compute Greek letters
                        alp0(i, j, n) = ((1.0/24.0)*(2*q0(i+1,j+1,n-1) + q0(i,j+1,n-1) + 2*q0(i,j,n-1) + q0(i+1,j,n-1)))
        		bet0(i, j, n) = ((1.0/24.0)*(q0(i,j+1,n-1) + 2*q0(i-1,j+1,n-1) + q0(i-1,j,n-1) + 2*q0(i,j,n-1)))
        		gam0(i, j, n) = ((1.0/24.0)*(2*q0(i,j+1,n-1) + q0(i-1,j+1,n-1) + 2*q0(i-1,j,n-1) + q0(i,j,n-1)))
        		del0(i, j, n) = ((1.0/24.0)*(q0(i+1,j+1,n-1) + 2*q0(i,j+1,n-1) + q0(i,j,n-1) + 2*q0(i+1,j,n-1)))
        		eps0(i, j, n) = ((1.0/24.0)*(q0(i+1,j+1,n-1) + q0(i,j+1,n-1) - q0(i,j,n-1) - q0(i+1,j,n-1)))
        		phi0(i, j, n) = ((1.0/24.0)*(-q0(i+1,j+1,n-1) + q0(i,j+1,n-1) + q0(i,j,n-1) - q0(i+1,j,n-1))) 
			
			! Calculate u
   			u0(i,j,n) = u0(i,j,n-1) + dT * (alp0(i,j,n-1) * vs0(i,j+1,n-1) + bet0(i,j,n-1) * vs0(i-1,j+1,n-1) + &
                        gam0(i,j,n-1) * vs0(i-1,j,n-1) + del0(i,j,n-1) * vs0(i,j,n-1) - eps0(i,j,n-1) * &
			us0(i+1,j,n-1) + eps0(i-1,j,n-1) * us0(i-1,j,n-1) - (ken0(i,j,n-1) - &
			ght0(i,j,n-1) + ken0(i,j,n-1) + ght0(i,j,n-1)) / Dxm)
			
			! Calculate v 
                        v0(i,j,n) = v0(i,j,n-1) - dT * (gam0(i+1,j,n-1) * us0(i+1,j,n-1) + del0(i,j,n-1) * us0(i,j,n-1) + &
                        alp0(i,j-1, n-1) * us0(i,j-1,n-1) + bet0(i+1,j-1,n-1) * us0(i+1,j-1,n-1) + phi0(i,j,n-1) * &
			vs0(i,j+1,n-1) - phi0(i,j-1,n-1) * vs0(i,j-1,n-1) - (ken0(i,j,n-1) + &
			ght0(i,j,n-1) - ken0(i,j-1,n-1) - ght0(i,j-1,n-1)) / Dxm)
			
			! Calculate zeta (z0)
			z0(i,j,n) = (u0(i,j-1,n-1) - u0(i,j,n-1) + v0(i,j,n-1) - v0(i-1,j,n-1)) / Dxm
			
			! Calculate Potential Vorticity
			q0(i,j,n) = (f + z0(i,j,n-1)) / hq0(i,j,n-1)

			! Update everything else
			hu0(i,j,n) = (h0(i-1,j,n-1) + h0(i+1,j,n-1))/2.0
        		hv0(i,j,n) = (h0(i,j-1,n-1) + h0(i,j+1,n-1))/2.0
        		z0(i,j,n)=(u0(i,j-1,n-1)-u0(i,j+1,n-1)+v0(i+1,j,n-1)-v0(i-1,j,n-1))/Dxm
        		hq0(i,j,n)=(h0(i,j,n-1) + h0(i-1,j,n-1) + h0(i-1,j-1,n-1) + h0(i,j-1,n-1))/4.0
        		ght0(i,j,n) = g*(hsurf(i,j) + h0(i,j,n-1))
        		ken0(i,j,n)=(u0(i,j,n-1)**2 + u0(i+1,j,n-1)**2 + v0(i,j,n-1)**2 + v0(1,j+1,n-1)**2)/4.
	  
		end do
	end do

	us0(2:Nx-1,2:Ny-1,n)= hu0(2:Nx-1,2:Ny-1,n)*u(2:Nx-1,2:Ny-1)
	vs0(2:Nx-1,2:Ny-1,n)= hv0(2:Nx-1,2:Ny-1,n)*v(2:Nx-1,2:Ny-1)


	alp0(1,:,1) = alp0(Nx-1,:,1)
	alp0(Nx,:,1) = alp0(2,:,1)
	bet0(1,:,1) = bet0(Nx-1,:,1)
	bet0(Nx,:,1) = bet0(2,:,1)
	gam0(1,:,1) = gam0(Nx-1,:,1)
	gam0(Nx,:,1) = gam0(2,:,1)
	del0(1,:,1) = del0(Nx-1,:,1)
	del0(Nx,:,1) = del0(2,:,1)
	eps0(1,:,1) = eps0(Nx-1,:,1)
	eps0(Nx,:,1) = eps0(2,:,1)
	phi0(1,:,1) = phi0(Nx-1,:,1)
	phi0(Nx,:,1) = phi0(2,:,1)

	hu0(1,:,1) = hu0(Nx-1,:,1)
	hu0(Nx,:,1) = hu0(2,:,1)

	hv0(1,:,1) = hv0(Nx-1,:,1)
	hv0(Nx,:,1) = hv0(2,:,1)

	z0(1,:,1) = z0(Nx-1,:,1) 
	z0(Nx,:,1) = z0(2,:,1) 

	hq0(1,:,1) = hq0(Nx-1,:,1) 
	hq0(Nx,:,1) = hq0(2,:,1) 

	ght0(1,:,1) = ght0(Nx-1,:,1) 
	ght0(Nx,:,1) = ght0(2,:,1) 

	ght0(1,:,1) = ght0(Nx-1,:,1) 
	ght0(Nx,:,1) = ght0(2,:,1) 

	q0(1,:,1) = q0(Nx-1,:,1) 
	q0(Nx,:,1) = q0(2,:,1) 

	h0(1,:,1) = h0(Nx-1,:,1) 
	h0(Nx,:,1) = h0(2,:,1) 

end do ! time loop

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


open(unit=10, file='topo.dat', access='direct', form='unformatted', status='unknown', action='write', recl=4*Nx*Ny)  !scary thing
write(10, rec=1) hsurf


end program lowresmodel


subroutine initgrid(hsurf, Dx, Nx, Ny)	!creates grid, ridge, and topography.
    implicit none
  
    integer, parameter :: Lx = 6000, Ly = 2000
    integer, intent(in) :: Dx, Nx, Ny

    real, dimension(Nx,Ny), intent(out) :: hsurf
    integer :: i, j 

    hsurf(:,:)=0.0

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


