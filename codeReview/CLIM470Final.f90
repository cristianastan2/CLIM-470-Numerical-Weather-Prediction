
PROGRAM finalproject

!Make actual variable and parameters
! need to define Nx,Ny for each vector
! define ntime, how many dT 
! compile code, look for syntax errors, run code (as of 11/21)
! one subroutine involved in this project 
! use print statements to help with debugging the code 

IMPLICIT NONE 

INTEGER :: d,Lx,Ly,Nx,Ny,h0,nstep,i,j,n,ntime,ierror        
REAL :: delt, f_1, f_2, f_3, hs_t
REAL, DIMENSION(:,:,:), ALLOCATABLE :: alp0,bet0,gam0,del0,eps0,ken0,phi0,alp1,bet1,gam1,del1,eps1,&
    ken1,phi1,alp2,bet2,gam2,del2,eps2,ken2,phi2,alp3,bet3,gam3,del3,eps3,ken3,phi3,q0,z0,hu0,hu1,&
    hu2,hu3,hv0,hv1,hv2,hv3,hq0,us0,vs0,phi,ght0,vs1, vs2, vs3,us1, us2, us3, q, ght1, ght2, ght3
REAL, DIMENSION(:,:), ALLOCATABLE :: u, v, h, z, ken, ght
REAL, Dimension(:), ALLOCATABLE :: hs
REAL, parameter:: f_cor=10e-04, g=9.8



open (unit = 10, file = '470dat1.txt', status = 'new', action = 'readwrite', iostat = ierror)
!!!!! ierror is not defined
open(unit = 99, file = 'info.dat', status = 'new', action = 'readwrite', iostat = ierror)
!!!!!! same comment as above
!!!!!! after adding ierror to the list of variables both files have been created
!!!!!! I wrote the initial conditions for h,u and v in info.dat. See the write statement below. 

!Initialize domain and resolution
    Lx = 6e+06                                  !domain size in x direction, real numbers
    Ly = 2e+06                                  !domain size in y direction, real numbers

    !d = 5e+05                                   !Resolution (changes for each run)
    !d = 2.5e+05
    d = 1.25e+05

    Nx = Lx/d + 1                               !number of grid points in the x direction (13, 25, 49), integers
    Ny = Ly/d + 1                               !number of grid points in the y direction (5, 9, 17), integers
    
    
    hs_t = 2e+03                                !height of the topography
    ! need to initialize h0, cannot equal/exceed top of the model or 5000! 
   

    !delt = 5*60                                !time step (10 min) in seconds (changes for each run)
    !delt = 5*60
    delt = 2.5*60

    f_1= (23*delt)/12                           !alphas for 3rd order adams-bashforth schemes, double check these 
    f_2 = (4*delt)/3
    f_3 = (5*delt)/12	
    

allocate(hs(Nx))
hs(1:Nx)=0.

!Establish topography

SELECT CASE (d)
	CASE (500000)                       !topography variable for first resolution
        	hs(Nx/2) = hs_t
	CASE (250000)                    !topography variables for second resolution
        	hs(Nx/2-1) = 1e+03
        	hs(Nx/2) = hs_t
        	hs(Nx/2+1) = 1e+03

	CASE (125000)                    !topography variable for third resolution
        	hs(Nx/2-2) = 0.5e+03
        	hs(Nx/2-1) = 1.5e+03
        	hs(Nx/2) = hs_t
        	hs(Nx/2+1) = 0.5e+03
        	hs(Nx/2+2) = 0.5e+03
 END SELECT

                                                !Time discretization using 3rd order adams-bashforth scheme 

do i = 1,Nx
 print*, 'hs(',i,') = ',hs(i)
end do

!!!!!!!!!! This printout shows that hs is not properly initialized. See the output file
!!!!!!!!! after allocate (hs(Nx)) the following statement must be added: hs(1:Nx)=0.

!initialize variables 



allocate(hu0(Nx,Ny,3)) 
allocate(hu1(Nx,Ny,3)) 
allocate(hu2(Nx,Ny,3)) 
allocate(hu3(Nx,Ny,3))
allocate(hq0(Nx,Ny,3))

allocate(hv0(Nx,Ny,3))
allocate(hv1(Nx,Ny,3))
allocate(hv2(Nx,Ny,3))
allocate(hv3(Nx,Ny,3))

allocate(us0(Nx,Ny,3))
allocate(us1(Nx,Ny,3))
allocate(us2(Nx,Ny,3))
allocate(us3(Nx,Ny,3))

allocate(vs0(Nx,Ny,3)) 
allocate(vs1(Nx,Ny,3)) 
allocate(vs2(Nx,Ny,3)) 
allocate(vs3(Nx,Ny,3)) 

allocate(alp0(Nx,Ny,3)) 
allocate(alp1(Nx,Ny,3)) 
allocate(alp2(Nx,Ny,3)) 
allocate(alp3(Nx,Ny,3)) 

allocate(bet0(Nx,Ny,3))
allocate(bet1(Nx,Ny,3))
allocate(bet2(Nx,Ny,3))
allocate(bet3(Nx,Ny,3))

allocate(gam0(Nx,Ny,3))
allocate(gam1(Nx,Ny,3))
allocate(gam2(Nx,Ny,3))
allocate(gam3(Nx,Ny,3))

allocate(del0(Nx,Ny,3))
allocate(del1(Nx,Ny,3))
allocate(del2(Nx,Ny,3))
allocate(del3(Nx,Ny,3))

allocate(eps0(Nx,Ny,3))
allocate(eps1(Nx,Ny,3))
allocate(eps2(Nx,Ny,3))
allocate(eps3(Nx,Ny,3))

allocate(ken0(Nx,Ny,3))

allocate(phi(Nx,Ny,3))
allocate(phi0(Nx,Ny,3))
allocate(phi1(Nx,Ny,3))
allocate(phi2(Nx,Ny,3))
allocate(phi3(Nx,Ny,3))

allocate(ght0(Nx,Ny,3))
allocate(ght1(Nx,Ny,3))
allocate(ght2(Nx,Ny,3))
allocate(ght3(Nx,Ny,3))
 
allocate(q0(Nx,Ny,3))
allocate(z0(Nx,Ny,3))
!allocate(z(Nx,Ny))
allocate(q(Nx,Ny,3))                                       ! absolute potential vorticity                                       
allocate(ght(Nx,Ny))                                     ! geopotential 
allocate(ken(Nx,Ny)) 					 !kinetic energy 
allocate(ken1(Nx,Ny,3)) 
allocate(ken2(Nx,Ny,3)) 
allocate(ken3(Nx,Ny,3)) 

allocate(h(Nx,Ny)) !Error: Shape specification for allocatable scalar at (1)
allocate(u(Nx,Ny))
allocate(v(Nx,Ny))                                        ! vertical velocity  !!!!!!! The shallow water model is a 2D model. It does not have vertical velocity; only zonal (u) and meridional (v) velocity!!!!!!!!
allocate(z(Nx,Ny))
h0 = 4999		                                     !top of fluid
!add remaining new variables and initial conditions for q etc (Done?)

!initial conditions, t=0 or n=1

u(1:Nx,1:Ny) = 0.5                          !value of u
v(1:Nx,1) = 0                               !value of v
v(1:Nx,Ny) = 0                              !value of v
V(1:Nx,2:Ny-1) = 0.1                        !horizontal velocity 


do i = 1, Nx
    h(i,:) = h0 - hs(i)   !we only have the height of the topography so we account for fluid above
end do 

write(99,*)h
write(99,*)u
write(99,*)v
stop

hu0(1,:,1) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,1) = (h(Nx-1,:) + h(1,:))/2.0

do i = 2, Nx-1
    hu0(i,:,1) = (h(i-1,:) + h(i+1,:))/2.0
end do

hv0(:,1,1) = (h(:,Ny) +h(:,2))/2.0
hv0(:,Ny,1) = (h(:,Ny-1) + h(:,1))/2.0

do j = 2, Ny-1
    hv0(:,j,1) = (h(:,j-1) + h(:,j+1))/2.0
end do 

us0(:,:,1) = hu0(:,:,1)*u(:,:)
vs0(:,:,1) = hv0(:,:,1)*v(:,:) 

do n = 2,3

do i = 2, Nx-1
    do j = 2, Ny-1 
        h(i,j) = h(i,j)-delt*(us0(i+1,j+1,n-1)-us0(i,j+1,n-1)+vs0(i+1,j+1,n-1)-vs0(i+1,j,n-1))/d
    end do
end do

hu0(1,:,n) = (h(Nx,:) + h(2,:))/2.0
hu0(Nx,:,n) = (h(Nx-1,:) + h(1,:))/2.0

do i = 2, Nx-1
    hu0(i,:,n) = (h(i-1,:) + h(i+1,:))/2.0
end do 

hv0(:,1,n) = (h(:,Ny) +h(:,2))/2.0
hv0(:,Ny,n) = (h(:,Ny-1) + h(:,1))/2.0

do j = 2, Ny-1
    hv0(:,j,n) = (h(:,j-1) + h(:,j+1))/2.0
end do 

us0(:,:,n) = hu0(:,:,n)*u(:,:)          ! I included these two again, will double check for errors 
vs0(:,:,n) = hv0(:,:,n)*v(:,:)

end do
!begin momentum :mike:

do i = 1, Nx
    z0(i,1,:)=0.
    z0(i,Ny,:)=0.
end do

do j=2,Ny-1
    z0(1,j,1)=(u(1,j-1)-u(1,j+1)+v(2,j)-v(Nx,j))/d
    z0(Nx,j,1)=(u(Nx,j-1)-u(1,j+1)+v(1,j)-v(Nx-1,j))/d
end do

do i=2,Nx-1
    z0(i,1,1)=(-u(i,2)+v(i+1,1)-v(i-1,1))/d
    z0(i,Ny,1)=(u(i,Ny-1)+v(i+1,Ny)-v(i-1,Ny))/d
end do

do i=2,Nx-1
	do j=2,Ny-1
    		z0(i,j,1)=(u(i,j-1)-u(i,j+1)+v(i+1,j)-v(i-1,j))/d
	end do
end do

hq0(1,:,1)=(h(1,1) + h(Nx,1))/2.0                !What are these two lines? Are these in the correct spot? Yes, they're part of the initial conditions of the model  
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
    		hq0(i,j,0)=(h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
	end do
end do

do j = 1,Ny
	do i = 1,Nx
    		q0(i,j,1)=(z0(i,j,0) + f_cor)/hq0(i,j,0)
	end do
end do

do i = 1,Nx
	do j = 1, Ny
    		ght0(i,j,1) = g*(hs(i) + h(i,j))
	end do
end do

do i = 1, Nx
    ken0(i,Ny,1)=0.
end do

do j = 1, Ny-1
    ken0(Nx,j,1)=(u(Nx-1,j)**2 + u(1,j)**2 + v(Nx-1,j)**2 + v(1,Ny+1)**2)/4.0    
end do

do i = 1, Nx-1
	do j = 1, Ny-1
   		 ken0(i,j,1) =(u(i,j)**2 + u(i+1,j)**2 + v(i,j)**2 + v(1,j+1)**2)/4.0     
	end do
end do

!**if this works, do we need the three different **# vars for hu and hv? can we just assign them to the first three spaces of the array?
                        !as so:
                        !hu(:,:,1)= hu0 (:,:,1)
                        !hu(:,:,2)= hu0 (:,:,2)
                        !hu(:,:,3)= hu0 (:,:,3)
                        !hv(:,:,1) = hv0(:,:,1) 
                        !hv(:,:,2) = hv0(:,:,2)
                        !hv(:,:,3) = hv0(:,:,3) !**it cuts down on vars, look through the rest of the code to see if it could be implemented
			
! define momentum coefficients from eq. 3.13 as completely specified in 3.34 (***do we need to do this still?)
eps0(:,:,n) = (q(i+1,j+1,n)+q(i,j+1,n)-q(i,j,n)-q(i+1,j,n))/24.0

phi0(:,:,n) = (-q(i+1,j+1,n)+q(i,j+1,n)+q(i,j,n)-q(i+1,j,n))/24.0

alp0(:,:,n) = (2*q(i+1,j+1,n)+q(i,j+1,n)+2*q(i,j,n)+q(i+1,j,n))/24.0

bet0(:,:,n) = (q(i,j+1,n)+2*q(i-1,j+1,n)+q(i-1,j,n)+2*q(i,j,n))/24.0

gam0(:,:,n) = (2*q(i,j+1,n)+q(i-1,j+1,n)+2*q(i-1,j,n)+q(i,j,n))/24.0

del0(:,:,n) = (q(i+1,j+1,n)+2*q(i,j+1,n)+q(i,j,n)+2*q(i+1,j,n))/24.0 

do n = 2,3
    do i = 2, Nx-1 
        do j = 2, Ny-1
            h(i,j) = h(i,j)-delt*(us0(i+1,j+1,n-1)) - (us0(i,j+1,n-1)) + (vs0(i+1,j+1,n-1)) - ((vs0(i+1,j,n-1))/d)
	    
            u(i,j) = u(i,j)+delt*(alp0(i,j+1,n-1)*vs0(i,j,n-1)+bet0(i,j+1,n-1)*vs0(i-1,j+1,n-1)+ gam0(i,j+1,n-1)*vs0(i-1,j,n-1) &
	    +del0(i,j+1,n-1)*vs0(i+1,j,n-1)-eps0(i+1,j+1,n-1)*us0(i+1,j+1,n-1)+eps0(i-1,j+1,n-1)*us0(i-1,j+1,n-1) &
            -(ken0(i+1,j+1,n-1)+ght0(i+1,j+1,n-1)-ken0(i-1,j+1,n-1)-ght0(i-1,j+1,n-1))/d) 

            v(i,j) = (v(i,j))-(delt*(gam0(i+1,j+1,n-1)))*(us0(i+1,j+1,n-1))+(del0(i,j+1,n-1))*(us0(i,j+1,n-1)) &
             +(alp0(i,j-1,n-1))*(us0(i,j-1,n-1))+(bet0(i+1,j-1,n-1))*(us0(i+1,j-1,n-1))+ (phi0(i+1,j+1,n-1))*(vs0(i+1,j+1,n-1)) &
	         -(phi0(i+1,j-1,n-1))*(vs0(i+1,j-1,n-1))- (ken0(i+1,j+1,n-1))+(ght0(i+1,j+1,n-1)) &
	         -(ken0(i+1,j-1,n-1))-((phi(i+1,j-1,n-1))/d)
  
            z0(1,j,n)=(u(1,j-1)-u(1,j+1)+v(2,j)-v(Nx,j))/d 
            z0(Nx,j,n)=(u(Nx,j-1)-u(1,j+1)+v(1,j)-v(Nx-1,j))/d 
        end do 
    end do
end do

! what is hq0??????? add in parentheses
!hq0(:,:,n) = (h(i+1/2,j+1/2,n)+h(i-1/2,j+1/2,n)+h(i-1/2,j-1/2,n)+h(i+1/2,j-1/2,n))/4.0            !Added eqn 3.15 from Arakawa for the def of hq0 for a square grid 
q0(:,:,n) = (f_cor+z0(:,:,n))/hq0(:,:,n) 

do i = 1,Nx
   ght0(i,:,n) = g*(h(i,:)+hs(i))
end do

ken0(:,:,n)=(u(:,:)*u(:,:) + v(:,:)*v(:,:))/2.0

us1(:,:,1) = us0(:,:,1)  !Might need to move these to line 340 and change back to original 2D arrays, will explain 
us2(:,:,2) = us0(:,:,2)
us3(:,:,3) = us0(:,:,3)
vs1(:,:,1) = vs0(:,:,1)
vs2(:,:,2) = vs0(:,:,2)
vs3(:,:,3) = vs0(:,:,3)
hu1(:,:,1) = hu0(:,:,1) !***Tried to assign 2 dimensional array to an allocated 3 dimensional array :: hotfix; give 2d arrays (hu#,hv#) a static third dim.
hu2(:,:,2) = hu0(:,:,2)
hu3(:,:,3) = hu0(:,:,3)
hv1(:,:,1) = hv0(:,:,1) 
hv2(:,:,2) = hv0(:,:,2)
hv3(:,:,3) = hv0(:,:,3)
alp1(:,:,1) = alp0(:,:,1) 
alp2(:,:,2) = alp0(:,:,2)
alp3(:,:,3) = alp0(:,:,3)
bet1(:,:,1) = bet0(:,:,1)
bet2(:,:,2) = bet0(:,:,2)
bet3(:,:,3) = bet0(:,:,3)
gam1(:,:,1) = gam0(:,:,1)
gam2(:,:,2) = gam0(:,:,2)
gam3(:,:,3) = gam0(:,:,3)
del1(:,:,1) = del0(:,:,1)
del2(:,:,2) = del0(:,:,2)
del3(:,:,3) = del0(:,:,3)
eps1(:,:,1) = eps0(:,:,1)
eps2(:,:,2) = eps0(:,:,2)
eps3(:,:,3) = eps0(:,:,3)
ken1(:,:,1) = ken0(:,:,1)
ken2(:,:,2) = ken0(:,:,2)
ken3(:,:,3) = ken0(:,:,3)
ght1(:,:,1) = ght0(:,:,1)
ght2(:,:,2) = ght0(:,:,2)
ght3(:,:,3) = ght0(:,:,3)       !Added these back, will need them for the momentum equations next 

! added in u(i,j) and v(i,j) for momentum
nstep = 4

do n = 4, 1440
    nstep = nstep + 1
    do i = 2, Nx-1
        do j = 2, Ny-1
            h(i,j) = h(i,j)-f_1*(us1(i+1,j+1,1)-us1(i,j+1,1)+vs1(i+1,j+1,1)-vs1(i+1,j,1)) &
                     +f_2*(us2(i+1,j+1,2)-us2(i,j+1,2)+vs2(i+1,j+1,2)-vs2(i+1,j,2)) &
                     -f_3*(us3(i+1,j+1,3)-us3(i,j+1,3)+vs3(i+1,j+1,3)-vs3(i+1,j,3))
            
            u(i,j) = u(i,j) + f_1*(alp1(i,j+1,1)*vs1(i,j,1)+bet1(i,j+1,1)*vs1(i-1,j+1,1) &
	    	+gam1(i,j+1,1)*vs1(i-1,j,1)+del1(i,j+1,1)*vs1(i+1,j,1) &
		-eps1(i+1,j+1,1)*us1(i+1,j+1,1)+eps1(i-1,j+1,1)*us1(i-1,j+1,1) &
		-(ken1(i+1,j+1,1)+ght1(i+1,j+1,1)-ken1(i-1,j+1,1)-ght1(i-1,j+1,1))/d) &
		-f_2*(alp2(i,j+1,2)*vs2(i,j,2)+bet2(i,j+1,2)*vs2(i-1,j+1,2)&
	    	+gam2(i,j+1,2)*vs2(i-1,j,2)+del2(i,j+1,2)*vs2(i+1,j,2)&
		-eps2(i+1,j+1,2)*us2(i+1,j+1,2)+eps2(i-1,j+1,2)*us2(i-1,j+1,2)&
		-(ken2(i+1,j+1,2)+ght2(i+1,j+1,2)-ken2(i-1,j+1,2)-ght2(i-1,j+1,2))/d)&      !***I think we can optomize this code: have all of the **# coefficients 
		+f_3*(alp3(i,j+1,3)*vs3(i,j,3)+bet3(i,j+1,3)*vs3(i-1,j+1,3)&                !be 3 dimension arrays, and just call (i*,j*,1), (i*,j*,2),(i*,j*,3) of each array. (can we also add the 0?) 
	    	+gam3(i,j+1,3)*vs3(i-1,j,3)+del3(i,j+1,3)*vs3(i+1,j,3)&                 !*Also, we may need to parenthesize the arrays in order for the compiler to run without error (see 266-277)
		-eps3(i+1,j+1,3)*us3(i+1,j+1,3)+eps3(i-1,j+1,3)*us3(i-1,j+1,3)&
		-(ken3(i+1,j+1,3)+ght3(i+1,j+1,3)-ken3(i-1,j+1,3)-ght3(i-1,j+1,3))/d)
            
	    v(i,j) = v(i,j) - f_1*(gam1(i+1,j+1,1)*us1(i+1,j+1,1)+del1(i,j+1,1)*us1(i,j+1,1)&
	    		+alp1(i,j-1,1)*us1(i,j-1,1)+bet1(i+1,j-1,1)*us1(i+1,j-1,1)&
			+phi1(i+1,j+1,1)*vs1(i+1,j+1,1)-phi1(i+1,j-1,1)*vs1(i+1,j-1,1)&
			-(ken1(i+1,j+1,1)+ght1(i+1,j+1,1)-ken1(i+1,j-1,1)-ght1(i+1,j-1,1))/d)&
			+f_2*(gam2(i+1,j+1,2)*us2(i+1,j+1,2)+del2(i,j+1,2)*us2(i,j+1,2)&
	    		+alp2(i,j-1,2)*us2(i,j-1,2)+bet2(i+1,j-1,2)*us2(i+1,j-1,2)&
			+phi2(i+1,j+1,2)*vs2(i+1,j+1,2)-phi2(i+1,j-1,2)*vs2(i+1,j-1,2)&
			-(ken2(i+1,j+1,2)+ght2(i+1,j+1,2)-ken2(i+1,j-1,2)-ght2(i+1,j-1,2))/d)&
			-f_3*(gam3(i+1,j+1,3)*us3(i+1,j+1,3)+del3(i,j+1,3)*us3(i,j+1,3)&
	    		+alp3(i,j-1,3)*us3(i,j-1,3)+bet3(i+1,j-1,3)*us3(i+1,j-1,3)&
			+phi3(i+1,j+1,3)*vs3(i+1,j+1,3)-phi3(i+1,j-1,3)*vs3(i+1,j-1,3)&
			-(ken3(i+1,j+1,3)+ght3(i+1,j+1,3)-ken3(i+1,j-1,3)-ght3(i+1,j-1,3))/d)
        end do
    end do
    
    us0(:,:,1) = hu0(:,:,1)*u(:,:)
    vs0(:,:,1) = hv0(:,:,1)*v(:,:)

    us1 = us2
    us2 = us3
    us3 = us0
    vs1 = vs2
    vs2 = vs3
    vs3 = vs0
    hu1 = hu2
    hu2 = hu3
    hu3 = hu0
    hv1 = hv2 
    hv2 = hv3
    hv3 = hv0 
    
    if (nstep == 1440) then
    	!write(10,*) "h",h
	!write(10,*) "u",u
	!write(10,*) "v",v
    	print*, "h",h
	print*, "u",u
	print*, "v",v
        !write(10,*) h, u, V
        nstep = 0
    end if 
end do

do i = 2, Nx-1
        hu0(i,:,1) = (h(i-1,:) + h(i+1,:))/2.0 !***need specified third dim in hu0 to assign the two dimensions in the eq to
end do

do j = 2, Ny-1
        hv0(:,j,1) = (h(:,j-1) + h(:,j+1))/2.0 !***need specified third dim in hv0 to assign the two dimensions in the eq to (same for the last few errors)
end do
    
us0(:,:,n) = hu0(:,:,n)*u(:,:)          ! I included these two again, will double check for errors 
vs0(:,:,n) = hv0(:,:,n)*v(:,:)

!if (((h(ken3+((1/2)*(g)*(h))+((g)*(hs)))-(h(ken0+((1/2)*(g)*(h))+((g)*(hs))))/(delt)) = 0) then 
!write(*,*) 'The scheme conserves the total energy.'
!else
!write(*,*) 'The scheme does not conserve the total energy.'


!do i=1,1440
    !if (nstep == 1440) then
    	!print*, "h",h
	!print*, "u",u
	!print*, "V",V
        !write(10,*) h, u, V
        !nstep = 0
   ! end if 
!end do !time loop 










end program finalproject
