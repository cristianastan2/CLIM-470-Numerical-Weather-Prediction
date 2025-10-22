program zonal_momentum
    
    real, parameter :: PI = 3.141592653589793
    real, parameter :: OMEGA = 7.2921e-5    ! Earth's rotation rate
    real, parameter :: R = 287.05           ! Gas constant for dry air
    real, parameter :: dx = 5000.0          ! Grid cell size (5 km)
    real, parameter :: dt = 1.0             ! time step of 1 second
    real, parameter :: total_time = 3600.0  ! 1 hour in seconds

    real, parameter :: phi = 35.0 * PI/180.0
    real, parameter :: T = 285.0             ! Temperature in kelvin
    real, parameter :: v = 5.0
    real, parameter :: pa1 = 101000.0        ! West boundary atmospheric pressure in pascals
    real, parameter :: pa2 = 100600.0        ! East boundary atmospheric pressure in pascals
    real, parameter :: u1 = 3.0              ! West u-wind
    real, parameter :: u2 = 4.0              ! East u-wind
    
    real :: f, rho_a, dudx, dpdx, u_center
    real :: time, local_accel, total_deriv
    integer :: i
    real, allocatable :: u_history(:), time_history(:)
    
    f = 2.0 * OMEGA * sin(phi)              ! Coriolis parameter
    rho_a = (pa1 + pa2)/2.0 / (R * T)       ! Air density using average pressure
    dudx = (u2 - u1) / dx                  
    dpdx = (pa2 - pa1) / dx                 
    
    u_center = (u1 + u2) / 2.0              ! u_center is the average between east and west wind pressure
    
    ! Allocate arrays for storing results
    allocate(u_history(0:int(total_time)))
    allocate(time_history(0:int(total_time)))
    
    ! initial conditions
    u_history(0) = u_center
    time_history(0) = 0.0
    
    ! Time loop
    do i = 1, int(total_time)
        time = i
        
        ! du/dt = -u*du/dx - v*du/dy - w*du/dz + f*v - (1/?)*dp/dx + Fr
        ! du/dy, dy/dz and Fr are all 0, so left them out of the equation

        local_accel = f*v - (1.0/rho_a)*dpdx - u_center*dudx
        
        ! forward Euler method
        u_center = u_center + local_accel * dt
        
        u_history(i) = u_center
        time_history(i) = time
    end do
    
    ! Write results to data file
    open(unit=10, file='u_change.dat', status='replace')
    write(10, *) '# Time(s) Time(min) u(m/s)'
    do i = 0, int(total_time)
        write(10, '(3F12.4)') time_history(i), time_history(i)/60.0, u_history(i)
    end do
    close(10)
    
    deallocate(u_history, time_history)
    
end program zonal_momentum
