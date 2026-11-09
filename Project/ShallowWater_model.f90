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
!
program ShallowWater_model

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
    ! real, dimension(i_Lx / i_resolution) :: xpos
    ! real, dimension(i_Ly / i_resolution) :: ypos
    ! real*4, dimension(i_nlon,i_nlat,i_numvars), allocatable :: idata
    ! n_zero, n_one, n_two are arrays for timesteps
    real*4, dimension(:, :, :), allocatable :: idata
    real*4, dimension(:, :, :), allocatable :: n_zero, n_one, n_two, n_steps
    real :: r_uwind, r_vwind, r_h, r_f      ! variables used by model
    real :: r_utend, r_vtend, r_htend    ! diagnostic variables for time steps
    integer :: i_deltat, i_stepnum          ! values for time steps

    integer :: col_u = 1        ! fortran is 0-based for array columns
    integer :: col_v = 2        ! v value
    integer :: col_h = 3        ! height of water column above surface
    integer :: col_utend = 4   ! tendency for this u
    integer :: col_vtend = 5
    integer :: col_htend = 6   !

    character(len=50) :: outfile
    ! *****  end of variable definitions

    ! ******  variable initialization follows
    a_resolution(1:3) = (/ 500, 250, 125 /)
    a_timestep_sec(1:3) = (/ 10*60 , 5*60, 25*60/)

    !****
    ! this outerloop wraps around the inner steps.  it uses a do loop
    ! to determine which resolutions to process
    ! the individual parts process (1) the initial topo, creating
    ! the idata array for later use in the timesteps.  Idata is also
    ! written out to a file for debugging and plotting
    !********
    outerdo: do i_outer = 1, 3  ! outer loop... one loop per resolution
        i_nlon = (i_Lx/a_resolution(i_outer))+1
        i_nlat = (i_Ly/a_resolution(i_outer))+1
        ! Note that arrays are for cell edges, not centers, hence +1
        allocate (idata(i_nlon, i_nlat, 1), stat = allocation_status)
        allocate (n_zero(i_nlon, i_nlat, 6))
        allocate (n_one(i_nlon, i_nlat, 6))
        allocate (n_two(i_nlon, i_nlat, 6))
        allocate (n_steps(i_nlon, i_nlat, 6))    ! for step 3 and beyond
        print*, " status of idata allocation: ", allocation_status
        i_recl = i_nlon*i_nlat*4
        i_deltat = a_timestep_sec(i_outer)

        ! *** now set up for writing to the output
        ! outfile="./Hs" + a_resolution(i_outer) + ".dat"
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
        !    now populate the array with problems values, and writ
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
        r_uwind = 20.0
        r_vwind = 0.0
        r_h = 4000.0            ! top surface of the fluid

        ! now calculate initial values. and add to n_zero array
        !do i_j = 1, i_nlat          ! process across each lattitude
        !    print*, "processing by lat for resolution: ", a_resolution(i_outer)
        !    print*, " and nlat: ", i_nlat, " currently i_j: ", i_j

        !    i_temp = (i_i - 1) * a_resolution(i_outer)
        !    n_zero(:, i_j, col_u) = r_uwind       ! set to initial constant from Arakawa
        !   n_zero(:, i_j, col_v) = r_vwind       ! u is specified as zonal wind, so v starts as 0
        !    n_zero(:, i_j, col_h) = r_h -  idata(1, i_j, 1)    ! h varies by long across the grid
        !end do

        r_f = 0.0000937        ! this Coriolis f is for lat 40N.  later we will add variable lat
        ! now calculate tendency values for first time-step
        i_stepnum = 0
        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat
        r_vtend = 0.0
        r_vwind = 0.0
        r_utend = 0.0
        r_uwind = 20.0
        r_htend = 0.0      ! h doesn't change over time yet, but h_s changes over lon
        r_h = r_h - r_htend
        print*, "timestep 0 tendencies u, v, h: ", r_utend, r_vtend, r_htend
        print*, "timestep 0 u, v, h: ", r_uwind, r_vwind, r_h

        dostep0: do i_i = 1, i_nlon          ! process across each lattitude
            ! print*, " processing timestep 1 at lat: ", i_j
            ! print*, "u, v, h values: ", r_uwind, r_vwind, r_h
            n_zero(i_i, :, col_u) = r_uwind
            n_zero(i_i, :, col_v) = r_vwind
            n_zero(i_i, :, col_h) = r_h -  idata(i_i, 1, 1)    ! h varies by long across the grid
            print*, "DEBUG in GRID: r_h and r_htend and i_i and idata val: ", &
                r_h, r_htend, i_i, idata(i_i, 1, 1)
            ! n_zero(:, i_j, col_h) = r_h
            ! and for diagnostic purposes
            n_zero(i_i, :, col_utend) = a_resolution(i_outer)     ! r_utend
            n_zero(i_i, :, col_vtend) = i_stepnum        ! r_vtend
            n_zero(i_i, :, col_htend) = r_htend
        end do dostep0    ! finished timestep 0 - initial values
        ! ********
        ! do initial set up for writing
        i_recl_steps = i_recl       ! timesteps have more data than initial,
                        ! but each record has the same amount, one variable
                        ! for all grid locations
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
        n = n + 1
        write(12,rec=n)n_zero(:,:,col_utend)
        n = n + 1
        write(12,rec=n)n_zero(:,:,col_vtend)
        print*, "output file ", outfile, "written for timestep 0"
        ! ****************************************
        ! DEBUGGING stuff - remove when working
        print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think
        !print*, "exiting"
        !close(12)
        !exit outerdo

        ! now process for timestep 1
        i_stepnum = i_stepnum +1
        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat
        r_utend = 0.0
        r_uwind = r_uwind - r_utend
        r_vtend = - (r_uwind * r_f * i_deltat)    ! tendency is neg or positive
        r_vwind = r_vwind + (r_vtend)
        r_htend = 0.0
        r_h = r_h - r_htend    ! h doesn't change over time yet, but h_s changes over lon
        print*, "timestep 1 tendencies u, v, h: ", r_utend, r_vtend, r_htend
        print*, "timestep 1 u, v, h: ", r_uwind, r_vwind, r_h

        dostep1: do i_i = 1, i_nlon          ! process across each lattitude
            !print*, " processing timestep 1 at lat: ", i_j
            !print*, "u, v, h values: ", r_uwind, r_vwind, r_h
            n_one(i_i, :, col_u) = r_uwind
            n_one(i_i, :, col_v) = r_vwind
            n_one(i_i, :, col_h) = r_h -  idata(i_i, 1, 1)
            ! and for diagnostic purposes
            n_one(i_i, :, col_utend) = a_resolution(i_outer)        !r_utend
            n_one(i_i, :, col_vtend) = i_stepnum            ! r_vtend
            n_one(i_i, :, col_htend) = r_htend
        end do dostep1! finished timestep 1
        ! **** write out arrays
        n = n+1
        write(12,rec=n)n_one(:,:,col_u)
        n = n + 1
        write(12,rec=n)n_one(:,:,col_v)
        n = n + 1
        write(12,rec=n)n_one(:,:,col_h)
        n = n + 1
        write(12,rec=n)n_one(:,:,col_utend)
        n = n + 1
        write(12,rec=n)n_one(:,:,col_vtend)
        print*, "output file ", outfile, "written for timestep 1"
        print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

        ! now process for timestep 2
        i_stepnum = i_stepnum + 1
        print*, " timestep for resolution: ", a_resolution(i_outer), " = ", i_deltat
        r_utend = (r_vtend)      ! use v tendency from previous step
        r_uwind = r_uwind + r_utend
        r_vtend = - (r_uwind * r_f * i_deltat)
        r_vwind = r_vwind + r_vtend
        r_htend = 0.0
        r_h = r_h - r_htend    ! h doesn't change over time yet, but h_s changes over lon
        print*, "timestep 2 tendencies u, v, h: ", r_utend, r_vtend, r_htend
        print*, "timestep 2 u, v, h: ", r_uwind, r_vwind, r_h

        dostep2: do i_i = 1, i_nlon          ! process across each lattitude
            !print*, " processing timestep 2 at lat: ", i_j
            !print*, "u, v, h values: ", r_uwind, r_vwind, r_h
            n_two(i_i, :, col_u) = r_uwind
            n_two(i_i, :, col_v) = r_vwind
            n_two(i_i, :, col_h) = r_h -  idata(i_i, 1, 1)
            ! and for diagnostic purposes
            n_two(i_i, :, col_utend) = i_stepnum            !r_utend
            n_two(i_i, :, col_vtend) = a_resolution(i_outer)        !r_vtend
            n_two(i_i, :, col_htend) = r_htend
        end do dostep2
        ! **** write out arrays
        n = n+1
        write(12,rec=n)n_two(:,:,col_u)
        n = n + 1
        write(12,rec=n)n_two(:,:,col_v)
        n = n + 1
        write(12,rec=n)n_two(:,:,col_h)
        n = n + 1
        write(12,rec=n)n_two(:,:,col_utend)
        n = n + 1
        write(12,rec=n)n_two(:,:,col_vtend)
        print*, "output file ", outfile, "written for timestep 2"
        print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

        ! notes on multiple variables in one array:
        ! Fortran stores arrays column-major, so in n_zero(i_nlon, i_nlat, 6)
        ! the first index (lon) varies fastest in the file, then lat, then the variable.
        ! Written with rec=n, one record holds all 6 layers back to back:
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
            r_utend = (r_vtend)      ! use v tendency from previous step
            r_uwind = r_uwind + r_utend
            r_vtend = - (r_uwind * r_f * i_deltat)
            r_vwind = r_vwind + r_vtend
            r_htend = 0.0
            r_h = r_h - r_htend    ! h doesn't change over time yet, but h_s changes over lon
            print*, "timestep", i_stepnum, " tendencies u, v, h: ", r_utend, r_vtend, r_htend
            print*, "timestep", i_stepnum, " u, v, h: ", r_uwind, r_vwind, r_h

            dostepn: do i_i = 1, i_nlon          ! process across each lattitude
                !print*, " processing timestep 2 at lat: ", i_j
                !print*, "u, v, h values: ", r_uwind, r_vwind, r_h
                n_steps(i_i, :, col_u) = r_uwind
                n_steps(i_i, :, col_v) = r_vwind
                n_steps(i_i, :, col_h) = r_h -  idata(i_i, 1, 1)
                ! and for diagnostic purposes
                n_steps(i_i, :, col_utend) = i_stepnum
                n_steps(i_i, :, col_vtend) = a_resolution(i_outer)
                n_steps(i_i, :, col_htend) = r_htend
            end do dostepn

            ! *****  write out data for timestep n
            n = n+1
            write(12,rec=n)n_steps(:,:,col_u)
            n = n + 1
            write(12,rec=n)n_steps(:,:,col_v)
            n = n + 1
            write(12,rec=n)n_steps(:,:,col_h)
            n = n + 1
            write(12,rec=n)n_steps(:,:,col_utend)
            n = n + 1
            write(12,rec=n)n_steps(:,:,col_vtend)
            print*, "output file ", outfile, "written for timestep ", i_stepnum
            print*, "array: ", n_zero(1,1, :)   ! print first row, all columns, I think

            if (i_stepnum > 20) exit do_all_steps

            ! move arrays back in time
            n_zero = n_one
            n_one = n_two
            n_two = n_steps

        end do do_all_steps


        ! ******* WRAP UP
        close(12)
        deallocate(n_zero, n_one, n_two, n_steps, idata)
        print*, " **** finished outer loop for iouter: ", i_outer
    end do outerdo         ! outer do loop

end program ShallowWater_model


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
