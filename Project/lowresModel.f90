program lowresmodel
implicit none

!variable declaration
real, parameter :: Lx = 6000, Ly = 2000                 !domain in km
real, parameter :: Dx = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dx)        !create grid points for array
integer :: i, j, t, nstep                               !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, hsurf, lh, lu, lv, lq   !array creation with internal variables for each point


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



subroutine firststep()
end subroutine firststep

subroutine secondstep()
end subroutine secondstep

subroutine thirdstep()
end subroutine thirdstep

subroutine momentum()	!math for velocity (WIP)
end subroutine momentum

subroutine height()	!math for determining scale height
endsubroutine height

subroutine vorticity()	!math for vorticity
end subroutine vorticity
