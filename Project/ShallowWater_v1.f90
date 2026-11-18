!
! this program creates the grid for the bottom topography for
! CLIM470 Project.  The Plus version adds more model features,
! including initial values and 2 timesteps for u, v and h.
! Note that the ridge has not been included
!
! the overall grid is Lx = 6000km, and Ly = 2000km
! grid resolutions (for different simulation runs) are 500, 250 and 125km
!
! in the center is a ridge.  It starts at x = 2500km,
! has a 2km peak at x = 3000km and the other edge is at 3500km.

! the software will create a grid with each cell populated by Hs value
! number of grid cells created is a function of resolution
!
program ShallowWater_v1

    ! integer, parameter  :: nlon=144, nlat=72, ntime=20
    integer, parameter :: i_Lx = 6000, i_Ly = 2000   ! grid size in km
    ! integer, parameter :: i_resolution = 125        ! grid resolution
            ! alternatives are 500, 250, 125
    integer :: i_nlon, i_nlat       ! (i_Lx/i_resolution)+1
    ! integer, parameter :: i_nlat = (i_Ly/i_resolution)+1
    integer :: allocation_status
    integer, parameter :: i_mtnstart = 2500, i_mtnend = 3500 ! ridge parms

    integer, parameter :: i_numvars = 1    ! number variables in output
            ! will only do surface height Hs initially
    real :: r_Hs = 0
    integer :: i_i, i_outer      ! for do loop

    integer, dimension(3) :: a_resolution, a_timestep_sec
    integer :: i_recl, i_recl_steps      ! record length for write
    ! idata is an array for the bottom surface topography
    ! n_zero, n_one, n_two, and n_steps are arrays for timesteps
    real*4, dimension(:, :, :), allocatable :: idata    ! bottom surface topo
    real*4, dimension(:, :, :), allocatable :: n_zero, n_one, n_two, n_steps
    real :: r_uwind, r_vwind, r_h, r_f      ! variables used by model
    real :: r_utend, r_vtend, r_htend    ! diagnostic variables for time steps
    integer :: i_deltat, i_stepnum          ! values for time steps

! mel - variable for u_tend
real :: dx, dy
real, parameter :: grav = 9.81

    integer :: col_u = 1        ! u value (fortran is 1-based for array columns)
    integer :: col_v = 2        ! v value
    integer :: col_h = 3        ! height of water column above surface
    integer :: col_utend = 4   ! tendency for this u
    integer :: col_vtend = 5    ! tendency for v
    integer :: col_htend = 6   ! perhaps useful for geopotential height
                ! note that not all 6 columns need to be written
                ! columns to be written are selected during the write process
    logical :: b_runsteps = .false.
    character(len=50) :: outfile    ! used for file name in open statements
    ! *****  end of variable definitions

    ! ******  variable initialization follows
    a_resolution(1:3) =  (/ 500, 250, 125 /)
    a_timestep_sec(1:3) = (/ 10*60 , 5*60, 25*60/)  ! from Arakawa

    !****
    ! this outerloop wraps around the inner steps.  it uses a do loop
    ! to determine which resolutions to process
    ! PART 1 produces the initial topo, creating
    ! the idata array for later use in the timesteps.  Idata is also
    ! written out to a file for debugging and plotting
    ! PART 2 builds arrays for timesteps 1, 2 3
    !********
    outerdo: do i_outer = 1, 1 ! outer loop... one loop per resolution
        i_nlon = (i_Lx/a_resolution(i_outer))+1
        i_nlat = (i_Ly/a_resolution(i_outer))+1
        ! Note that arrays are for cell edges, not centers, hence +1
        allocate (idata(i_nlon, i_nlat, 1), stat = allocation_status)
        allocate (n_zero(i_nlon, i_nlat, 6))   ! 6 values associated with each lon, lat
        allocate (n_one(i_nlon, i_nlat, 6))
        allocate (n_two(i_nlon, i_nlat, 6))
        allocate (n_steps(i_nlon, i_nlat, 6))    ! for step 3 and beyond
        print*, " status of idata allocation: ", allocation_status
        i_recl = i_nlon*i_nlat*4
        i_deltat = a_timestep_sec(i_outer)

! mel grid spacing for u_tend
dx = real(a_resolution(i_outer))*1000.0
dy = dx

        ! *** now set up for writing to the output
        ! outfile="./Hs" + a_resolution(i_outer) + ".dat"
        ! it's BEST to follow this convention so the plotting program can find
        ! your file....
        write(outfile, '(A,I0,A)') './Hs', a_resolution(i_outer), '.dat'
        open(10,file=outfile,form='unformatted',action='write', &
                access='direct',recl=i_recl)
        print*, "outfile for idata setup: ", outfile, i_recl
        ! recl is record length.  intel fortran has recl unit = 4-byte words
        ! for gfortran it is bytes....
        ! so i_nlon*inlat is right for intel.  i_nlon*i_nlat*4 for gfortran

        ! *************  end of initialization
        !
        ! **** PART 1 - CALCULATE TOPO  ***
        !
        !    now populate the array with problems values, and write out data file
        dotopo: do i_i = 1, i_nlon
            print*, "processing for resolution: ", a_resolution(i_outer)
            print*, "currently i_i: ", i_i, " and i_nlon", i_nlon
            r_Hs = 0.0
            i_temp = (i_i - 1) * a_resolution(i_outer)
            print*,"i_temp after init for this do (i_i): ", i_i, i_temp
            if ((i_temp <= i_mtnstart).or.(i_temp >= i_mtnend)) then
                    idata(i_i, :, 1) = r_Hs
            else
                CALL calc_Hs(i_temp, a_resolution(i_outer), i_mtnstart, i_mtnend, r_Hs)
                print*, "longitude and Mountain height: ", i_temp, r_Hs
                idata(i_i, :, 1) = r_Hs
            end if
            print*, "i_temp (0-6000), r_Hs: ", i_temp, r_Hs
        end do dotopo

        n = 1
        ! the topo data is one real value per grid cell.  All are written
        ! at once.  The record lenght is set in the OPEN statement above (i_recl)
        ! each value is a single precision real number (4 bytes).
        write(10,rec=n)idata(:,:,1)
        print*, "output file ", outfile, "written"
        close(10)
        ! deallocate (idata)   ! need to keep idata for remaining calcs

        ! ********* PART2
        ! calculate initial values and create arrays for timesteps 0, 1 and 2
        ! timestep 0 array is populated with initial values
        ! timestep 1 and 2 arrays are calculate with Euler function
        !
        ! now setup initial values, from arakawa paper
! initialize model arrays
n_zero = 0.0
n_one = 0.0
n_two = 0.0
        r_uwind = 20.0
        r_vwind = 0.0
        r_h = 4000.0            ! top surface of the fluid

        r_f = 0.0000937        ! this Coriolis f is for lat 40N.  later we will add variable lat
        ! now calculate tendency values for first time-step
        i_stepnum = 0

        ! THE REMAINDER OF THE CODE FOR THIS STEP IS TO BE PROVIDED BY TEAM MEMBERS

! zonal (u-direction), by Mel (however, may be extra since Nathan initialized)
! maybe delete first n_zero and then keep verify init zonal
! Using Arakawa and Lamb, init zonal uniform at 20 m/s
! Store initial u values at timestep 0

n_zero(:,:,col_u) = r_uwind

! verify init zonal

print*, "Initial zonal wind (m/s): ", n_zero(1,1,col_u)

! end zonal wind init mel

! start lorenzo
! meridional wind (v) initialize
! v = 0 initially for horizontal free surface (arakawa)
! timestep 0 for v

n_zero(:,:,col_v) = r_vwind

!verify init meridional

print*, "Initial Meridional wind (m/s):", n_zero(1,1,col_v)

!end meridional wind init lorenzo

! mel calculate initial u_tend
! du/dt = q*vstar - d(K + Phi)/dx

call calculate_u_tendency( &
n_zero(:,:,col_u), &
n_zero(:,:,col_v), &
n_zero(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_zero(:,:,col_utend))

! lorenzo calculate initial v_tend
!dv/dt = -q*ustar - d(k + Phi)/dy

call calculate_v_tendency( &
n_zero(:,:,col_u), &
n_zero(:,:,col_v), &
n_zero(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_zero(:,:,col_vtend))

! process timestep 1

i_stepnum = i_stepnum +1

! mel timestep 1 zonal wind
! fwd Euler: u(1) = u(0) + dt * du/dt(0)

n_one(:,:,col_u) = n_zero(:,:,col_u) + &
real(i_deltat)*n_zero(:,:,col_utend)

! lorenzo timestep 1 meridional wind
! fwd Euler: v(1) = v(0) + dt * dv/dt(0)

n_one(:,:,col_v) = n_zero(:,:,col_v) + & 
real(i_deltat)*n_zero(:,:,col_vtend)

! calculate u_tend at timestep 1

call calculate_u_tendency( &
n_one(:,:,col_u), &
n_one(:,:,col_v), &
n_one(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_one(:,:,col_utend))

call calculate_v_tendency( &
n_one(:,:,col_u), &
n_one(:,:,col_v), &
n_one(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_one(:,:,col_vtend))

! end mel timestep 1
! end lorenzo timestep 1

print*, "timestep for resolution:", a_resolution(i_outer), "=", i_deltat
print*, "output file", outfile, "NOT written for timestep 1"
print*, "array:", n_zero(1,1,:)

        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat
        print*, "timestep 0 u, v, h: TBD"

        ! This code is here as an example of how the data files for the steps
        ! can (should?) be written.  b_runsteps is false, so the code is skipped
        if (b_runsteps) then
            write(outfile, '(A,I0,A)') './timesteps', a_resolution(i_outer), '.dat'
            open(12,file=outfile,form='unformatted',action='write', &
                    access='direct',recl=i_recl_steps)
            print*, "outfile for timesteps setup: ", outfile, i_recl_steps
            n = 0       ! each write will add 1

            ! **** write out arrays
            n = n+1
            write(12,rec=n)n_zero(:,:,col_u)
            n = n + 1
            write(12,rec=n)n_zero(:,:,col_v)
            n = n + 1
            write(12,rec=n)n_zero(:,:,col_h)
            print*, "output file ", outfile, "written for timestep 0"
        end if
        ! ****************************************
        ! DEBUGGING stuff - remove when working
        ! print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think
        !print*, "exiting"
        !close(12)
        !exit outerdo

        ! now process for timestep 1

        i_stepnum = i_stepnum +1
        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat
        print*, "output file ", outfile, "NOT written for timestep 1"
        print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

        ! now process for timestep 2
        i_stepnum = i_stepnum + 1
        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat

        ! TEAM MEMBERS INSERT CODE HERE

! mel zonal timestep 2
! Fwd Euler 2: u(2) = u(1) + dt * du/dt(1)

n_two(:,:,col_u) = n_one(:,:,col_u) + &
real(i_deltat)*n_one(:,:,col_utend)

! lorenzo meridional timestep 2
! fwd Euler 2: v(2) = v(1) + dt * dv/dt(1)

n_two(:,:,col_v) = n_one(:,:,col_v) + &
real(i_deltat)*n_one(:,:,col_vtend)

! calculate u_tend at timestamp 2

call calculate_u_tendency( &
n_two(:,:,col_u), &
n_two(:,:,col_v), &
n_two(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_two(:,:,col_utend))

! calculate v_tend at timestep 2

call calculate_v_tendency( &
n_two(:,:,col_u), &
n_two(:,:,col_v), &
n_two(:,:,col_h), &
idata(:,:,1), &
r_f, grav, dx, dy, &
n_two(:,:,col_vtend))


! end mel zonal wind timestep 2
! end lorenzo meridional wind timestep 2

        print*, "output file ", outfile, "NOT written for timestep 2"
        print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

        ! notes on multiple variables in one array:
        ! Fortran stores arrays column-major, so in n_zero(i_nlon, i_nlat, 6)
        ! the first index (lon) varies fastest in the file, then lat.
        ! the varialbles from one column are printed in one record
        ! Written with rec=n, one record holds all 3 (or more) layers back to back:
        ! the full lon x lat grid for variable 1, then variable 2, and so on.

        ! now have three timesteps written out
        ! close(12)
        ! leave file open .. will continue writing as we go through timesteps

        ! ******** PART 3 ********
        ! this part will calculate values for timesteps 3 through n
        !
        ! At the endcopy array data backward to support Adams

        do_all_steps: do
            ! now process for timestep 3 and until end
            i_stepnum = i_stepnum + 1
            print*, " timestep", i_stepnum, " for resolution: ", a_resolution(i_outer), " = ", i_deltat

            ! TEAM MEMBERS INSERT CODE HERE

            print*, "output file ", outfile, "NOT written for timestep ", i_stepnum
            print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

            if (i_stepnum > 20) exit do_all_steps

            ! move arrays back in time
        end do do_all_steps


        ! ******* WRAP UP
        close(12)
        deallocate(n_zero, n_one, n_two, n_steps, idata)
        print*, " **** finished outer loop for iouter: ", i_outer
    end do outerdo         ! outer do loop

end program ShallowWater_v1


SUBROUTINE calc_Hs(long, n, i_mtnstart, i_mtnend, r_Hs)
    implicit none
    integer, intent(in) :: long
    integer, intent(in) :: n, i_mtnstart, i_mtnend
    real, intent(out)  :: r_Hs
    real :: r_deltalong
    integer :: r_slopewidth
    integer :: i_ridge_long = 3000    ! longitude for ridge
    integer :: i_temp

    r_slopewidth = real(((i_mtnend - i_mtnstart)/2))
    if (long <= i_ridge_long) then
        i_temp = long-i_mtnstart
        print*, "long-i_mtnstart: ", i_temp
        !i_temp = i_temp / i_resolution
        !print*, "long - i_mtnstart) / i_resolution: ", i_temp
        r_Hs = ((real(long) - real(i_mtnstart)) / r_slopewidth) * 2000
        print*, "r_Hs: ", r_Hs
    else
        r_Hs = ((real(i_mtnend) - real(long)) / r_slopewidth) * 2000
    end if
END SUBROUTINE calc_Hs

! mel - provisional zonal tendency
! du/dt = q*vstat - d(K + Phi)/dx

subroutine calculate_u_tendency(u, v, h, hs, f, g, dx, dy, utend)
implicit none
real, intent(in) :: u(:,:), v(:,:), h(:,:), hs(:,:)
real, intent(in) :: f, g, dx, dy
real, intent(out) :: utend(:,:)

integer :: i, j, ip, im
integer :: nx, ny
real :: zeta, q, vstar
real :: k_left, k_right
real :: phi_left, phi_right
real :: dvdx, dudy
real :: grad_energy

nx = size(u,1)
ny = size(u,2)

utend = 0.0

! last x-column dupicates first periodic point

do j = 1, ny
do i = 1, nx-1

ip = i+1
im = i-1

if(i == 1) im = nx-1

! relative vort
dvdx = (v(ip,j)-v(im,j))/(2.0*dx)

if (j == 1) then
dudy = (u(i,j+1)-u(i,j))/dy
else if (j == ny) then
dudy = (u(i,j)-u(i,j-1))/dy
else
dudy = (u(i,j+1)-u(i,j-1))/(2.0*dy)
end if

zeta = dvdx-dudy

! potential vort and mass flux
q = (zeta+f)/h(i,j)
vstar = h(i,j)*v(i,j)

! kinetic energy
k_left = 0.5*(u(im,j)**2+v(im,j)**2)
k_right = 0.5*(u(ip,j)**2+v(ip,j)**2)

! geopotential
phi_left = g*(h(im,j)+hs(im,j))
phi_right = g*(h(ip,j)+hs(ip,j))

! zonal gradient of K + Phi
grad_energy = (k_right+phi_right &
-k_left-phi_left)/(2.0*dx)

! zonal wind tendency
utend(i,j) = q*vstar-grad_energy

end do
end do

! periodic x boundary
utend(nx,:) = utend(1,:)

end subroutine calculate_u_tendency

!lorenzo

subroutine calculate_v_tendency (u,v,h,hs,f,g,dx,dy,vtend)
implicit none
real, intent(in) :: u(:,:), v(:,:), h(:,:), hs(:,:)
real, intent(in) :: f, g, dx, dy
real, intent(out) :: vtend(:,:)

integer :: i, j, ip, im, jp, jm
integer :: nx, ny
real :: zeta, q, ustar
real :: k_down, k_up
real :: phi_down, phi_up
real :: dvdx, dudy
real :: grad_energy

nx = size(u,1)
ny = size(u,2)

vtend = 0.0

!north and south boundary j = 1 and j = ny, start j = 2
!last x point wraps around to starting point

do j = 2, ny-1
do i = 1, nx-1

ip = i + 1
im = i - 1

if(i==1) im = nx - 1

jp = j + 1
jm = j - 1

!relative vorticity
dvdx = (v(ip,j)-v(im,j))/(2.0*dx)
dudy = (u(i,jp) -u(i,jm))/(2.0*dy)

zeta = dvdx-dudy

!potential vorticity and mass flux
q = (zeta + f)/h(i,j)
ustar = h(i,j)*u(i,j)

!kinetic energy
k_down = 0.5*(u(i,jm)**2+v(i,jm)**2)
k_up = 0.5*(u(i,jp)**2+v(i,jp)**2)

!geopotential
phi_down = g*(h(i,jm)+hs(i,jm))
phi_up = g*(h(i,jp)+hs(i,jp))

!meridional gradient of K + phi
grad_energy = (k_up + phi_up - k_down - phi_down)/(2.0*dy)

!meridional wind tendency
vtend(i,j) = -q*ustar-grad_energy

end do
end do

!periodic x boundary
vtend(nx,:) = vtend(1,:)

end subroutine calculate_v_tendency
