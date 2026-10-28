!meridional wind (v) for shallow water model

program meridional_wind_v
implicit none

!variables and integers

integer :: i, j
integer :: nx, ny

real, parameter :: Lx = 6000.0 !x length (km)
real, parameter :: Ly = 2000.0 !y length (km)
real, parameter :: h_max = 2.0 !height of ridge (km)
real, parameter :: x_center = 3000.0 !ridge peak in x value (km)
real, parameter :: g = 9.8 !gravity, m/s^2

real, dimension (:,:), allocatable :: v
real, dimension (:,:), allocatable :: v_new
real, dimension (:,:), allocatable :: dvdt
real, dimension (:,:), allocatable :: h !for topography height (km)
!begin blueprint

!determine grid points (nx,ny)

!allocate a 2-d array calculating v at each grid point

!loops to progress through grid in x & y direction

!initial conditions

!map out topography

!identify topography effects on meridional wind & calculate

!update & store values as model progresses over time 

!deallocate arrays

end program meridional_wind_v
