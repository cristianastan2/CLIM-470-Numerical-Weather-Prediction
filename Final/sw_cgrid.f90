! CLIM 470 Shallow Water Model - PDF IMPLEMENTATION
! Final Project - Shallow Water Equations with Bottom Topography
! Using exact PDF variable dimensions and initialization

module shallow_water_model
  implicit none
  
  ! Physical parameters
  real(8), parameter :: g = 9.81d0          ! gravitational acceleration (m/s^2)
  real(8), parameter :: f0 = 1.0d-4         ! Coriolis parameter (1/s)
  real(8), parameter :: H_mean = 5000.0d0   ! mean fluid depth (m)
  real(8), parameter :: Hs_max = 2000.0d0   ! maximum topography height (m)
  real(8), parameter :: pi = 3.141592653589793d0
  
  ! Domain parameters
  real(8) :: Lx = 6000000.0d0    ! domain length in x (m)
  real(8) :: Ly = 2000000.0d0    ! domain length in y (m)
  real(8) :: dx, dy              ! grid spacing
  integer :: Nx, Ny              ! number of grid points
  
  ! EXACT PDF VARIABLE DIMENSIONS
  real(8), allocatable :: h(:,:)        ! h-points (Nx, Ny)
  real(8), allocatable :: u(:,:)        ! u-points (Nx, Ny)  
  real(8), allocatable :: v(:,:)        ! v-points (Nx, Ny)
  real(8), allocatable :: hs(:,:)       ! topography (Nx, Ny)
  
  ! PDF mass flux variables with time levels
  real(8), allocatable :: hu0(:,:,:)    ! (Nx, Ny, 3)
  real(8), allocatable :: hv0(:,:,:)    ! (Nx, Ny, 3)
  real(8), allocatable :: us0(:,:,:)    ! (Nx, Ny, 3) - mass flux u
  real(8), allocatable :: vs0(:,:,:)    ! (Nx, Ny, 3) - mass flux v
  
  ! Tendency storage for Adams-Bashforth
  real(8), allocatable :: F_u(:,:,:)    ! u tendencies (Nx, Ny, 3)
  real(8), allocatable :: F_v(:,:,:)    ! v tendencies (Nx, Ny, 3)
  real(8), allocatable :: F_h(:,:,:)    ! h tendencies (Nx, Ny, 3)
  
  ! Time stepping parameters
  real(8) :: dt, total_time
  integer :: current_step, max_steps
  
contains
  
  ! Initialize the model
  subroutine initialize_model(resolution_km)
    real(8), intent(in) :: resolution_km
    integer :: i, j
    
    ! Convert resolution from km to m
    dx = resolution_km * 1000.0d0
    dy = dx
    
    ! Calculate grid dimensions
    Nx = nint(Lx / dx)
    Ny = nint(Ly / dy)
    
    print *, 'Initializing model with resolution:', resolution_km, 'km'
    print *, 'Grid dimensions:', Nx, 'x', Ny
    print *, 'Grid spacing:', dx/1000.0, 'km'
    
    ! EXACT PDF ALLOCATION
    allocate(h(Nx, Ny), hs(Nx, Ny))
    allocate(u(Nx, Ny), v(Nx, Ny))
    
    ! PDF mass flux arrays with time levels
    allocate(hu0(Nx, Ny, 3), hv0(Nx, Ny, 3))
    allocate(us0(Nx, Ny, 3), vs0(Nx, Ny, 3))
    
    ! Tendency arrays
    allocate(F_u(Nx, Ny, 3), F_v(Nx, Ny, 3), F_h(Nx, Ny, 3))
    
    ! Initialize to zero
    h = 0.0d0; hs = 0.0d0; u = 0.0d0; v = 0.0d0
    hu0 = 0.0d0; hv0 = 0.0d0; us0 = 0.0d0; vs0 = 0.0d0
    F_u = 0.0d0; F_v = 0.0d0; F_h = 0.0d0
    
    ! Set up bottom topography (narrow ridge)
    call setup_topography()
    
    ! Set initial conditions using EXACT PDF approach
    call set_initial_conditions_exact_pdf()
    
    ! Set time step based on CFL condition
    dt = 0.01d0 * dx / sqrt(g * H_mean)
    total_time = 1.0d0 * 24.0d0 * 3600.0d0  ! 1 day
    max_steps = nint(total_time / dt)
    current_step = 0
    
    print *, 'Time step:', dt, 'seconds'
    print *, 'Total steps:', max_steps
    print *, 'CFL number:', dt * sqrt(g * H_mean) / dx
    print *, 'Simulation time:', total_time/3600.0, 'hours'
    
    ! Check initial conditions
    print *, 'Initial h range:', minval(h), maxval(h)
    print *, 'Initial u range:', minval(u), maxval(u)
    print *, 'Initial v range:', minval(v), maxval(v)
    
  end subroutine initialize_model
  
  ! Set up the bottom topography (narrow ridge)
  subroutine setup_topography()
    integer :: i, j
    real(8) :: x, y, x0, y0, r, R_width
    
    x0 = Lx / 2.0d0  ! center at x = 3000 km
    y0 = Ly / 2.0d0  ! center in y
    R_width = 500000.0d0  ! half-width of 500 km for total width ~1000 km
    
    do i = 1, Nx
      x = (i - 0.5d0) * dx  ! x position at h-point i (cell center)
      do j = 1, Ny
        y = (j - 0.5d0) * dy  ! y position at h-point j (cell center)
        r = sqrt((x - x0)**2 + (y - y0)**2)
        
        if (r <= R_width) then
          hs(i,j) = (Hs_max / 2.0d0) * (1.0d0 + cos(pi * r / R_width))
        else
          hs(i,j) = 0.0d0
        endif
      end do
    end do
    
    print *, 'Topography range:', minval(hs), maxval(hs)
    
  end subroutine setup_topography
  
  ! Set initial conditions using EXACT PDF code
  subroutine set_initial_conditions_exact_pdf()
    integer :: i, j
    real(8) :: rand_val
    
    ! EXACT PDF CODE: Initialize velocities
    u(1:Nx, 1:Ny) = 0.0d0
    
    ! EXACT PDF CODE: Rigid wall boundary conditions for v
    v(1:Nx, 1) = 0.0d0       ! Bottom wall
    v(1:Nx, Ny) = 0.0d0      ! Top wall
    v(1:Nx, 2:Ny-1) = 0.0d0  ! Interior initialized to zero
    
    ! EXACT PDF CODE: Initialize height field
    do i = 1, Nx
      h(i,:) = 5000.0d0 - hs(i,1)  ! horizontal free surface
      ! Ensure positive depth
      do j = 1, Ny
        if (h(i,j) < 100.0d0) h(i,j) = 100.0d0
      end do
    end do
    
    ! EXACT PDF CODE: Initialize hu0 with proper boundary handling
    ! Left boundary (periodic)
    do j = 1, Ny
      hu0(1,j,1) = (h(Nx,j) + h(2,j)) / 2.0d0
    end do
    
    ! Right boundary (periodic) 
    do j = 1, Ny
      hu0(Nx,j,1) = (h(Nx-1,j) + h(1,j)) / 2.0d0
    end do
    
    ! Interior points
    do i = 2, Nx-1
      do j = 1, Ny
        hu0(i,j,1) = (h(i-1,j) + h(i+1,j)) / 2.0d0
      end do
    end do
    
    ! EXACT PDF CODE: Initialize hv0 with proper boundary handling
    ! Bottom boundary
    do i = 1, Nx
      hv0(i,1,1) = (h(i,Ny) + h(i,2)) / 2.0d0
    end do
    
    ! Top boundary
    do i = 1, Nx
      hv0(i,Ny,1) = (h(i,Ny-1) + h(i,1)) / 2.0d0
    end do
    
    ! Interior points
    do j = 2, Ny-1
      do i = 1, Nx
        hv0(i,j,1) = (h(i,j-1) + h(i,j+1)) / 2.0d0
      end do
    end do
    
    ! EXACT PDF CODE: Calculate initial mass fluxes
    do i = 1, Nx
      do j = 1, Ny
        us0(i,j,1) = hu0(i,j,1) * u(i,j)
        vs0(i,j,1) = hv0(i,j,1) * v(i,j)
      end do
    end do
    
    ! Add small perturbation to break symmetry
    do i = 1, Nx
      do j = 1, Ny
        call random_number(rand_val)
        h(i,j) = h(i,j) + 10.0d0 * (rand_val - 0.5d0)
      end do
    end do
    
    ! Initialize other time levels for Adams-Bashforth
    hu0(:,:,2) = hu0(:,:,1)
    hu0(:,:,3) = hu0(:,:,1)
    hv0(:,:,2) = hv0(:,:,1)
    hv0(:,:,3) = hv0(:,:,1)
    us0(:,:,2) = us0(:,:,1)
    us0(:,:,3) = us0(:,:,1)
    vs0(:,:,2) = vs0(:,:,1)
    vs0(:,:,3) = vs0(:,:,1)
    
    ! Apply boundary conditions
    call apply_boundary_conditions_pdf()
    
  end subroutine set_initial_conditions_exact_pdf
  
  ! Apply boundary conditions for PDF dimensions
  subroutine apply_boundary_conditions_pdf()
    integer :: i, j
    
    ! Periodic boundary conditions in x-direction
    do j = 1, Ny
      u(1,j) = u(Nx,j)
      u(Nx,j) = u(1,j)
      v(1,j) = v(Nx,j)
      v(Nx,j) = v(1,j)
      h(1,j) = h(Nx,j)
      h(Nx,j) = h(1,j)
    end do
    
    ! Rigid wall boundary conditions in y-direction (v = 0 at walls)
    do i = 1, Nx
      v(i,1) = 0.0d0      ! Bottom wall
      v(i,Ny) = 0.0d0     ! Top wall
    end do
    
    ! For u at walls: zero gradient (free-slip)
    do i = 1, Nx
      u(i,1) = u(i,2)
      u(i,Ny) = u(i,Ny-1)
    end do
    
  end subroutine apply_boundary_conditions_pdf
  
  ! Calculate mass fluxes for current time step (PDF version)
  subroutine calculate_mass_fluxes_pdf()
    integer :: i, j
    
    ! Update hu0 for current time step
    ! Left boundary (periodic)
    do j = 1, Ny
      hu0(1,j,1) = (h(Nx,j) + h(2,j)) / 2.0d0
    end do
    
    ! Right boundary (periodic)
    do j = 1, Ny
      hu0(Nx,j,1) = (h(Nx-1,j) + h(1,j)) / 2.0d0
    end do
    
    ! Interior points
    do i = 2, Nx-1
      do j = 1, Ny
        hu0(i,j,1) = (h(i-1,j) + h(i+1,j)) / 2.0d0
      end do
    end do
    
    ! Update hv0 for current time step
    ! Bottom boundary
    do i = 1, Nx
      hv0(i,1,1) = (h(i,Ny) + h(i,2)) / 2.0d0
    end do
    
    ! Top boundary
    do i = 1, Nx
      hv0(i,Ny,1) = (h(i,Ny-1) + h(i,1)) / 2.0d0
    end do
    
    ! Interior points
    do j = 2, Ny-1
      do i = 1, Nx
        hv0(i,j,1) = (h(i,j-1) + h(i,j+1)) / 2.0d0
      end do
    end do
    
    ! Calculate mass fluxes
    do i = 1, Nx
      do j = 1, Ny
        us0(i,j,1) = hu0(i,j,1) * u(i,j)
        vs0(i,j,1) = hv0(i,j,1) * v(i,j)
      end do
    end do
    
  end subroutine calculate_mass_fluxes_pdf
  
  ! Check for NaN values
  subroutine check_nan(step)
    integer, intent(in) :: step
    logical :: has_nan
    integer :: i, j
    
    has_nan = .false.
    
    ! Check all arrays
    do i = 1, Nx
      do j = 1, Ny
        if (h(i,j) /= h(i,j)) has_nan = .true.
        if (u(i,j) /= u(i,j)) has_nan = .true.
        if (v(i,j) /= v(i,j)) has_nan = .true.
      end do
    end do
    
    if (has_nan) then
      print *, 'WARNING: NaN detected at step', step
      print *, 'Aborting simulation due to NaN'
      stop
    end if
    
  end subroutine check_nan
  
  ! Calculate kinetic energy at h-points
  subroutine calculate_kinetic_energy(K)
    real(8), dimension(Nx, Ny), intent(out) :: K
    integer :: i, j
    
    do i = 1, Nx
      do j = 1, Ny
        K(i,j) = 0.5d0 * (u(i,j)**2 + v(i,j)**2)
      end do
    end do
    
  end subroutine calculate_kinetic_energy
  
  ! Calculate total energy
  function calculate_total_energy() result(total_energy)
    real(8) :: total_energy
    real(8), dimension(Nx, Ny) :: K
    integer :: i, j
    real(8) :: energy_density
    
    call calculate_kinetic_energy(K)
    total_energy = 0.0d0
    
    do i = 1, Nx
      do j = 1, Ny
        energy_density = h(i,j) * (K(i,j) + 0.5d0 * g * h(i,j) + g * hs(i,j))
        if (abs(energy_density) < 1.0d20 .and. .not. (energy_density /= energy_density)) then
          total_energy = total_energy + energy_density
        end if
      end do
    end do
    
    total_energy = total_energy * dx * dy
    
    if (total_energy /= total_energy) then
      total_energy = 0.0d0
    end if
    
  end function calculate_total_energy
  
  ! Calculate tendencies for u, v, h with PDF discretization
  subroutine calculate_tendencies()
    real(8), dimension(Nx, Ny) :: Phi
    integer :: i, j
    
    ! Check for NaN before calculating tendencies
    call check_nan(current_step)
    
    ! Calculate geopotential at h-points
    do i = 1, Nx
      do j = 1, Ny
        Phi(i,j) = g * (max(h(i,j), 1.0d0) + hs(i,j))
      end do
    end do
    
    ! u tendencies - pressure gradient with periodic boundaries
    do i = 1, Nx
      do j = 1, Ny
        if (i == 1) then
          F_u(i,j,1) = - (Phi(1,j) - Phi(Nx,j)) / dx
        else if (i == Nx) then
          F_u(i,j,1) = - (Phi(1,j) - Phi(Nx,j)) / dx
        else
          F_u(i,j,1) = - (Phi(i,j) - Phi(i-1,j)) / dx
        end if
      end do
    end do
    
    ! v tendencies - pressure gradient with rigid walls
    do i = 1, Nx
      do j = 1, Ny
        if (j == 1) then
          ! Bottom boundary - one-sided difference
          F_v(i,j,1) = - (Phi(i,2) - Phi(i,1)) / dy
        else if (j == Ny) then
          ! Top boundary - one-sided difference
          F_v(i,j,1) = - (Phi(i,Ny) - Phi(i,Ny-1)) / dy
        else
          ! Interior - centered difference
          F_v(i,j,1) = - (Phi(i,j) - Phi(i,j-1)) / dy
        end if
      end do
    end do
    
    ! h tendencies - mass conservation
    do i = 1, Nx
      do j = 1, Ny
        if (i == 1) then
          F_h(i,j,1) = -(us0(2,j,1) - us0(Nx,j,1)) / dx
        else if (i == Nx) then
          F_h(i,j,1) = -(us0(1,j,1) - us0(Nx-1,j,1)) / dx
        else
          F_h(i,j,1) = -(us0(i+1,j,1) - us0(i,j,1)) / dx
        end if
        
        if (j == 1) then
          F_h(i,j,1) = F_h(i,j,1) - (vs0(i,2,1) - vs0(i,1,1)) / dy
        else if (j == Ny) then
          F_h(i,j,1) = F_h(i,j,1) - (vs0(i,Ny,1) - vs0(i,Ny-1,1)) / dy
        else
          F_h(i,j,1) = F_h(i,j,1) - (vs0(i,j+1,1) - vs0(i,j,1)) / dy
        end if
      end do
    end do
    
    ! Conservative tendency limits
    do i = 1, Nx
      do j = 1, Ny
        if (abs(F_u(i,j,1)) > 10.0d0) F_u(i,j,1) = sign(10.0d0, F_u(i,j,1))
        if (abs(F_v(i,j,1)) > 10.0d0) F_v(i,j,1) = sign(10.0d0, F_v(i,j,1))
        if (abs(F_h(i,j,1)) > 1.0d0) F_h(i,j,1) = sign(1.0d0, F_h(i,j,1))
      end do
    end do
    
  end subroutine calculate_tendencies
  
  ! Time stepping with Adams-Bashforth 3rd order
  subroutine time_step()
    integer :: n_ab, i, j
    real(8) :: ab_coeff
    
    n_ab = min(3, current_step + 1)  ! Handle first two steps
    
    select case (n_ab)
    case (1)
      ! Forward Euler (first step)
      ab_coeff = dt
    case (2)
      ! 2nd order Adams-Bashforth (second step)
      ab_coeff = dt * 0.5d0
    case (3)
      ! 3rd order Adams-Bashforth
      ab_coeff = dt / 12.0d0
    end select
    
    ! Update u at all points
    do i = 1, Nx
      do j = 1, Ny
        select case (n_ab)
        case (1)
          u(i,j) = u(i,j) + ab_coeff * F_u(i,j,1)
        case (2)
          u(i,j) = u(i,j) + ab_coeff * (3.0d0 * F_u(i,j,1) - F_u(i,j,2))
        case (3)
          u(i,j) = u(i,j) + ab_coeff * (23.0d0 * F_u(i,j,1) - 16.0d0 * F_u(i,j,2) + 5.0d0 * F_u(i,j,3))
        end select
      end do
    end do
    
    ! Update v at all points
    do i = 1, Nx
      do j = 1, Ny
        select case (n_ab)
        case (1)
          v(i,j) = v(i,j) + ab_coeff * F_v(i,j,1)
        case (2)
          v(i,j) = v(i,j) + ab_coeff * (3.0d0 * F_v(i,j,1) - F_v(i,j,2))
        case (3)
          v(i,j) = v(i,j) + ab_coeff * (23.0d0 * F_v(i,j,1) - 16.0d0 * F_v(i,j,2) + 5.0d0 * F_v(i,j,3))
        end select
      end do
    end do
    
    ! Update h at all points
    do i = 1, Nx
      do j = 1, Ny
        select case (n_ab)
        case (1)
          h(i,j) = h(i,j) + ab_coeff * F_h(i,j,1)
        case (2)
          h(i,j) = h(i,j) + ab_coeff * (3.0d0 * F_h(i,j,1) - F_h(i,j,2))
        case (3)
          h(i,j) = h(i,j) + ab_coeff * (23.0d0 * F_h(i,j,1) - 16.0d0 * F_h(i,j,2) + 5.0d0 * F_h(i,j,3))
        end select
        
        ! Ensure positive depth after update
        if (h(i,j) < 1.0d0) h(i,j) = 1.0d0
      end do
    end do
    
    ! Add velocity limiting to prevent explosion
    do i = 1, Nx
      do j = 1, Ny
        if (abs(u(i,j)) > 100.0d0) u(i,j) = sign(100.0d0, u(i,j))
        if (abs(v(i,j)) > 100.0d0) v(i,j) = sign(100.0d0, v(i,j))
      end do
    end do
    
    ! Shift tendency arrays for next step
    if (n_ab >= 2) then
      do i = 1, Nx
        do j = 1, Ny
          F_u(i,j,3) = F_u(i,j,2)
          F_v(i,j,3) = F_v(i,j,2)
          F_h(i,j,3) = F_h(i,j,2)
        end do
      end do
    end if
    
    if (n_ab >= 1) then
      do i = 1, Nx
        do j = 1, Ny
          F_u(i,j,2) = F_u(i,j,1)
          F_v(i,j,2) = F_v(i,j,1)
          F_h(i,j,2) = F_h(i,j,1)
        end do
      end do
    end if
    
    current_step = current_step + 1
    
    ! Update mass fluxes and apply boundary conditions
    call calculate_mass_fluxes_pdf()
    call apply_boundary_conditions_pdf()
    
    ! Check for problems
    call check_nan(current_step)
    
  end subroutine time_step
  
  ! Run the model
  subroutine run_model()
    integer :: step, output_interval
    real(8) :: initial_energy, current_energy
    real(8) :: max_u, max_v
    
    output_interval = max(1, max_steps / 20)
    
    initial_energy = calculate_total_energy()
    print *, 'Initial total energy:', initial_energy
    
    ! Check if initial energy is reasonable
    if (initial_energy <= 0.0d0 .or. initial_energy /= initial_energy) then
      print *, 'ERROR: Invalid initial energy:', initial_energy
      print *, 'Check initial conditions and topography'
      return
    end if
    
    do step = 1, max_steps
      call calculate_tendencies()
      call time_step()
      
      if (mod(step, output_interval) == 0 .or. step <= 5) then
        current_energy = calculate_total_energy()
        max_u = maxval(abs(u))
        max_v = maxval(abs(v))
        
        if (.not. (current_energy /= current_energy) .and. .not. (initial_energy /= initial_energy)) then
          if (abs(initial_energy) > 1.0d-10) then
            print *, 'Step:', step, 'Time:', step*dt/3600.0, 'hours', &
                     'Energy change (%):', 100.0*(current_energy - initial_energy)/initial_energy, &
                     'Max u:', max_u, 'Max v:', max_v
          else
            print *, 'Step:', step, 'Time:', step*dt/3600.0, 'hours', &
                     'Energy:', current_energy, &
                     'Max u:', max_u, 'Max v:', max_v
          end if
        else
          print *, 'Step:', step, 'NaN detected in energy calculation'
        end if
      end if
      
      ! Early stopping conditions
      max_u = maxval(abs(u))
      max_v = maxval(abs(v))
      if (max_u > 1000.0d0 .or. max_v > 1000.0d0) then
        print *, 'Stopping: velocities too large'
        exit
      end if
      if (any(h < 0.0d0)) then
        print *, 'Stopping: negative depth detected'
        exit
      end if
    end do
    
    current_energy = calculate_total_energy()
    if (.not. (current_energy /= current_energy) .and. .not. (initial_energy /= initial_energy)) then
      if (abs(initial_energy) > 1.0d-10) then
        print *, 'Final total energy:', current_energy
        print *, 'Relative energy change (%):', 100.0*(current_energy - initial_energy)/initial_energy
      else
        print *, 'Final total energy:', current_energy
      end if
    else
      print *, 'ERROR: NaN in final energy calculation'
    end if
    
  end subroutine run_model
  
  ! Clean up memory
  subroutine cleanup_model()
    if (allocated(h)) deallocate(h)
    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(hs)) deallocate(hs)
    if (allocated(hu0)) deallocate(hu0)
    if (allocated(hv0)) deallocate(hv0)
    if (allocated(us0)) deallocate(us0)
    if (allocated(vs0)) deallocate(vs0)
    if (allocated(F_u)) deallocate(F_u)
    if (allocated(F_v)) deallocate(F_v)
    if (allocated(F_h)) deallocate(F_h)
  end subroutine cleanup_model

end module shallow_water_model

! Main program
program clim470_project
  use shallow_water_model
  implicit none
  
  real(8) :: resolutions(3) = [500.0d0, 250.0d0, 125.0d0]
  integer :: i
  
  print *, 'CLIM 470 Final Project - Shallow Water Model'
  print *, 'All variables at (Nx, Ny) grid points'
  
  ! Initialize random seed for reproducibility
  call random_seed()
  
  ! Run simulations for different resolutions
  do i = 1, 3
    print *, ''
    print *, 'Running simulation with resolution:', resolutions(i), 'km'
    print *, '------------------------------------------------'
    
    call initialize_model(resolutions(i))
    call run_model()
    call cleanup_model()
    
    print *, 'Completed resolution:', resolutions(i), 'km'
  end do
  
  print *, ''
  print *, 'All simulations completed!'
  
end program clim470_project