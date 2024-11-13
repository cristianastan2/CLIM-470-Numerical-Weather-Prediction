program lowresmodel

implicit none

!variable declaration
real, parameter :: Lx = 6000, Ly = 2000                 !domain in km
real, parameter :: Dx = 500, Dy = 500                   !grid resolution in km
real, parameter :: dT = 100                             !time step, 100 seconds
integer, parameter :: Nx = (Lx/Dx), Ny = (Ly/Dy)        !create grid points for array
integer :: i, j, t, nstep                               !define other potentially important variables
real :: time                                            !TIME

real, dimension(Nx, Ny) :: h, u, v, q, hsurf, lh, lu, lv, lq   !array creation with internal variables for each point


!showcase model info
print *, "Grid Type: C"
print *, "Domain (km):", Lx, Ly
print *, "Grid Resolution (km):", Dx, Dy
print *, "X-Direction Grid Point #:", Nx
print *, "Y-Direction Grid Point #:", Ny

!initialize main grid and topopgraphy
call initgrid()

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

end program lowresmodel

subroutine initgrid()	!creates grid, ridge, and topography.

real, dimension(:,:), intent(out) :: h, u, v, q, hsurf
integer :: i, j

!initial conditions set to 0
h = 0
u = 0
v = 0
q = 0

!create first perturbations
h(Nx/2, Ny/2) = 2		!this creates a small bump in scale height towards the middle of the model

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
