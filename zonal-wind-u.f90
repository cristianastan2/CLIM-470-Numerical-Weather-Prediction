
! shallow water model, zonal (u), from Arakawa and Lamb (1981)

program zonal_wind_u

implicit none

!declare variables

integer :: i, j
integer :: nx, ny

real, parameter :: u_initial = 20.0
real, dimension(:,:), allocatable :: u
real, dimension(:,:), allocatable :: u_new
rean, dimension(:,:), allocatable :: dudt


!concept blueprint

!define number of grid points in x and y (nx, ny above)

!allocate a two-dimensional array to store zonal at each point

!set initial to give 20 m/s in paper across entire domain

!use loops to move through each grid in x and y dir

!assign initial to appropriate grid point

!identify topographical interactions

!determind how topo affects zonal field

!calculate zonal wind using shallow water eqns

!update zonal values as model progesses through time

!store updated field for later calcs and plots

!deallocate arrays when done

end program zonal_wind_u
