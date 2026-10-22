
program hw4p5
    implicit none

    integer :: i, nsteps
    real(8) :: dx, dt, v, lat, temp
    real(8) :: p1, p2, pavg, rho
    real(8) :: u1, u2, u, dudx, dpdx
    real(8) :: omega, rd, pi, f
    real(8) :: dudt, time

    ! values given
    dx = 5000.0d0
    dt = 1.0d0
    nsteps = 3600
    v = 5.0d0
    lat = 35.0d0
    temp = 285.0d0

    p1 = 101000.0d0
    p2 = 100600.0d0
    u1 = 3.0d0
    u2 = 4.0d0

    ! constants
    pi = acos(-1.0d0)
    omega = 7.2921d-5
    rd = 287.0d0

    ! coriolis param
    f = 2.0d0 * omega * sin(lat * pi / 180.0d0)

    ! density air
    pavg = (p1 + p2) / 2.0d0
    rho = pavg / (rd * temp)

    ! pressure and wind gradients
    dpdx = (p2 - p1) / dx
    dudx = (u2 - u1) / dx

    ! initial wind at center of cell
    u = (u1 + u2) / 2.0d0

    ! create output file
    open(unit=10, file='zonal_wind.dat', status='replace')
    write(10,*) 0.0d0, u

    ! wind speed every second for one hour
    do i = 1, nsteps
        dudt = -u * dudx + f * v - dpdx / rho

        u = u + dt * dudt
        time = i * dt

        write(10,*) time, u
    end do

    close(10)

    print *, 'Simulation finished!'
    print *, 'Final zonal wind:', u, 'm/s'

end program hw4p5
