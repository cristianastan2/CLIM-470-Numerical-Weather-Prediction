!
! this program creates the grid for the bottom topography for
! CLIM470 Project.
!
! the overall grid is Lx = 6000km, and Ly = 2000km
! grid resolutions (for different simulation runs) are 500, 250 and 125km
!
! in the center is a ridge.  It starts at x = 2500km,
! has a 2km peak at x = 3000km and the other edge is at 3500km.

! the software will create a grid with each cell populated by Hs value
!
program ProjectTopo

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

    integer, dimension(3) :: a_resolution
    integer :: i_recl       ! record length for write
    ! real, dimension(i_Lx / i_resolution) :: xpos
    ! real, dimension(i_Ly / i_resolution) :: ypos
    ! real*4, dimension(i_nlon,i_nlat,i_numvars), allocatable :: idata
    real*4, dimension(:, :, :), allocatable :: idata

    character(len=50) :: outfile

    ! ******  variable initialization follows
    a_resolution(1:3) = (/ 125, 250, 500 /)

    do i_outer = 1, 3  ! outer loop... one loop per resolution
        i_nlon = (i_Lx/a_resolution(i_outer))+1
        i_nlat = (i_Ly/a_resolution(i_outer))+1
        allocate (idata(i_nlon, i_nlat, i_numvars), stat = allocation_status)
        print*, " status of idata allocation: ", allocation_status
        i_recl = i_nlon*i_nlat*4

        ! *** now set up for writing to the output
        ! outfile="./Hs" + a_resolution(i_outer) + ".dat"
        write(outfile, '(A,I0,A)') './Hs', a_resolution(i_outer), '.dat'
        open(10,file=outfile,form='unformatted',action='write', &
                access='direct',recl=i_recl)
        ! recl is record length.  intel fortran has recl unit = 4-byte words
        ! for gfortran it is bytes....
        ! so i_nlon*inlat is right for intel.  i_nlon*i_nlat*4 for gfortran

        ! *************  end of initialization
        !    now populate the array with problems values, and writ
        do i_i = 1, i_nlon
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
        end do

        n = 1
        write(10,rec=n)idata(:,:,1)
        print*, "output file ", outfile, "written"
        close(10)
        deallocate (idata)
    end do          ! outer do loop

end program ProjectTopo


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
