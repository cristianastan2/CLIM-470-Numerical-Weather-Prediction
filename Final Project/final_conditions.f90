program shallow_water_model
      implicit none
      
      !parameters! 
      integer, parameter::Lx = 6e+06 !domain size in x-direction (m)
      integer, parameter::Ly = 2e+06 !domain size in y-direction (m)
      real, parameter::hs_top = 2e+03 !height of the mountain (m)
      real, parameter::g = 9.8 !the acceleration of gravity (m/s^2)
      real, parameter::f = 1e-04 !the Coriolis parameter (s^-1)
      integer, parameter::dt = 1 !time step (s)
      integer, parameter::ntime = 1440 !(s)
      integer::nstep
      real, parameter::f1 = 23.0/12.0
      real, parameter::f2 = -4.0/3.0
      real, parameter::f3 = 5.0/12.0
      
      !resolution!
      real:: d !model resolution, ie delta_x, delta_y

      !variables!
      integer:: Nx !number of grid points in x-direction (13, 25, 49)
      integer:: Ny !number of grid points in y-direction (5, 9, 17)
      integer:: i, j, ii, jj, ierr, n !working variables
      real, allocatable::vor(:,:) !vorticity
      real, allocatable::q(:,:) !potential vorticity
      real, allocatable::pe(:,:) !potential energy
      real, allocatable::ke(:,:) !kinetic energy
      real, allocatable::hs(:,:), u(:,:), v(:,:), h(:,:)

      !store the results from forward scheme(Euler Scheme) for each time steps(n=2, n=3)
      real, allocatable::hu0(:,:,:), hv0(:,:,:), us0(:,:,:), vs0(:,:,:), hq0(:,:,:)
      real, allocatable::vor0(:,:,:), q0(:,:,:), pe0(:,:,:), ke0(:,:,:)
      real, allocatable::alp0(:,:,:), bet0(:,:,:), gam0(:,:,:), del0(:,:,:), eps0(:,:,:), pi0(:,:,:) !constant

      !set up for progressing Adams-Bashforth third order scheme!
      real, allocatable::us1(:, :), us2(:, :), us3(:, :), us(:, :)
      real, allocatable::vs1(:, :), vs2(:, :), vs3(:, :), vs(:, :)
      real, allocatable::hu(:, :)
      real, allocatable::hv(:, :)
      real, allocatable::alp1(:, :), alp2(:, :), alp3(:, :), alp(:, :)
      real, allocatable::bet1(:, :), bet2(:, :), bet3(:, :), bet(:, :)
      real, allocatable::gam1(:, :), gam2(:, :), gam3(:, :), gam(:, :)
      real, allocatable::del1(:, :), del2(:, :), del3(:, :), del(:, :)
      real, allocatable::eps1(:, :), eps2(:, :), eps3(:, :), eps(:, :)
      real, allocatable::pi1(:, :), pi2(:, :), pi3(:, :), pi(:, :)
      real, allocatable::pe1(:, :), pe2(:, :), pe3(:, :)
      real, allocatable::ke1(:, :), ke2(:, :), ke3(:, :)
      real, allocatable::hq(:, :)

      !resolution!
      !d = 5e+05
      !d = 2.5e+05
      d = 1.25e+05

      Nx = Lx/d + 1
      Ny = Ly/d + 1

      !allocate hs variable!
      allocate(hs(Nx, Ny))
      hs(1:Nx, 1:Ny) = 0.0
      hs((Nx+1)/2, 1:Ny) = hs_top

      if (d==2.5e+05) then
              hs((Nx+1)/2-1, 1:Ny) = 1e+03
              hs((Nx+1)/2+1, 1:Ny) = 1e+03
      else if (d==1.25e+05) then
              hs((Nx+1)/2-3, 1:Ny) = 0.5e+03
              hs((Nx+1)/2-2, 1:Ny) = 1e+03
              hs((Nx+1)/2-1, 1:Ny) = 1.5e+03
              hs((Nx+1)/2+1, 1:Ny) = 1.5e+03
              hs((Nx+1)/2+2, 1:Ny) = 1e+03
              hs((Nx+1)/2+3, 1:Ny) = 0.5e+03
      endif

      !!!============initial conditions============!!!
      !allocate u variable!
      allocate(u(Nx,Ny))
      u(1:Nx, 1:Ny) = 20.0 !m/s, initial condition on reference paper

      !allocate v variable!
      allocate(v(Nx, Ny))
      v(1:Nx, 1) = 0.0 !rigid computational boundary
      v(1:Nx, Ny) = 0.0 !rigid computational boundary
      v(1:Nx, 2:Ny-1) = 10.0

      !allocate h variable!
      allocate(h(Nx,Ny))
      do i = 1,Nx
      h(i,:) = 5e+03 - hs(i, :)!in m, initial height "hzero" defined in Arakawa and Lamb 1981, vertical extent fluid column above the bottom surface
      end do
      
      !allocate hu0 vairable!
      allocate(hu0(Nx,Ny,3))
      hu0(1,:,1) = (h(Nx,:) + h(2,:))/2.0 !arithmetic average described in paper for h^u
      hu0(Nx,:,1) = (h(Nx-1,:) + h(1,:))/2.0
      do i = 2, Nx-1
      hu0(i,:,1) = (h(i-1,:) + h(i+1,:))/2.0
      end do
      
      !allocate hv0 variable!
      allocate(hv0(Nx,Ny,3))
      hv0(:,1,1) = (h(:,Ny) + h(:,2))/2.0 !arithmetic average described in paper for h^v
      hv0(:,Ny,1) = (h(:,Ny-1) + h(:,1))/2.0
      do j = 2, Ny-1
      hv0(:,j,1) = (h(:,j-1) + h(:,j+1))/2.0
      end do
      
      !allocate us0 variable!
      allocate(us0(Nx, Ny, 3))
      us0(1:Nx, 1:Ny, 1) = hu0(1:Nx, 1:Ny, 1)*u(1:Nx, 1:Ny)

      !allocate vs0 variable!
      allocate(vs0(Nx, Ny, 3))
      vs0(1:Nx, 1:Ny, 1) = hv0(1:Nx, 1:Ny, 1)*v(1:Nx, 1:Ny)

      !create the initial data (u, v, h) files for each resolutions!
      open(11, file='u_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(11, rec=1)u
      close(11)
      
      open(12, file='v_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(12, rec=1)v
      close(12)

      open(13, file='h_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(13, rec=1)h
      close(13)

      !allocate vor0 variable!
      allocate(vor0(Nx,Ny,3))
      do i = 1, Nx
      vor0(i,1,1) = 0.0
      vor0(i,Ny,1) = 0.0
      end do

      do j = 2, Ny-1
      vor0(1,j,1) = (u(1,j-1) - u(1,j+1) + v(2,j) - v(Nx, j))/d
      vor0(Nx,j,1) = (u(Nx,j-1) - u(1,j+1) + v(1,j) - v(Nx-1,j))/d
      end do

      do i = 2, Nx-1
       do j = 2, Ny-1
        vor0(i,j,1) = (u(i,j-1) - u(i,j+1) + v(i+1,j) - v(i-1,j))/d
       end do
      end do

      !allocate hq0 variable!
      allocate(hq0(Nx,Ny,3))
      hq0(1,1,1) = (h(1,1) + h(Nx,1))/2.0
      hq0(1,Ny,1) = (h(1,Ny-1) + h(Nx,Ny-1))/2.0

      do i = 2, Nx
       hq0(i,1,1) = (h(i,1) + h(i-1, 1))/2.0
       hq0(i,Ny,1) = (h(i,Ny-1) + h(i-1,Ny-1))/2.0
      end do

      do j = 2, Ny-1
       hq0(1,j,1) = (h(1,j) + h(Nx,j) + h(Nx, j-1) + h(1,j-1))/4.0
      end do

      do j = 2, Ny-1
       do i = 2, Nx
        hq0(i,j,1) = (h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
       end do
      end do

      !allocate q0 variable!
      allocate(q0(Nx,Ny,3))
      do j = 1, Ny
       do i = 1, Nx
        q0(i,j,1) = (vor0(i,j,1) + f)/hq0(i,j,1)
       end do
      end do

      !allocate pe0 variable!
      allocate(pe0(Nx,Ny,3))
      do i = 1, Nx
       do j = 1, Ny
        pe0(i,j,1) = g * (hs(i,j) + h(i,j))
       end do
      end do

      !allocate ke0 variable!
      allocate(ke0(Nx,Ny,3))
      do i = 1, Nx
       ke0(i,Ny,1) = 0.0
      end do

      do j = 1, Ny-1
       ke0(Nx,j,1) = (u(Nx-1,j)**2 + u(1,j)**2 + v(Nx-1,j)**2 + v(1,j+1)**2)/4.0
      end do

      do i = 1, Nx-1
       do j = 1, Ny-1
        ke0(i, j, 1) = (u(i, j)**2 + u(i+1, j)**2 + v(i, j)**2 + v(i, j+1)**2)/4.0
       end do
      end do

      !allocate alp0, bet0, gam0, del0, eps0, pi0 forcing coefficients!
      allocate(alp0(Nx, Ny, 3), bet0(Nx, Ny, 3), gam0(Nx, Ny, 3), del0(Nx, Ny, 3), eps0(Nx, Ny, 3), pi0(Nx, Ny, 3))
      alp0(:, :, 1) = 0.0
      bet0(:, :, 1) = 0.0
      gam0(:, :, 1) = 0.0
      del0(:, :, 1) = 0.0
      eps0(:, :, 1) = 0.0
      pi0(:, :, 1) = 0.0
      
      do i = 1, Nx-1
       do j = 1, Ny-1
        alp0(i, j, 1) = (2.0*q0(i+1, j+1, 1) + q0(i, j+1, 1) + 2.0*q0(i, j, 1) + q0(i+1, j, 1))/24.0
        del0(i, j, 1) = (q0(i+1, j+1, 1) + 2.0*q0(i, j+1, 1) + q0(i, j, 1) + 2.0*q0(i+1, j, 1))/24.0
        eps0(i, j, 1) = (q0(i+1, j+1, 1) + q0(i, j+1, 1) - q0(i, j, 1) - q0(i+1, j, 1))/24.0
        pi0(i, j, 1) = (-q0(i+1, j+1, 1) + q0(i, j+1, 1) + q0(i, j, 1) - q0(i+1, j, 1))/24.0
       end do
      end do

      do j = 1, Ny-1
       alp0(Nx, j, 1) = (2.0*q0(1, j+1, 1) + q0(Nx, j+1, 1) + 2.0*q0(Nx, j, 1) + q0(1, j, 1))/24.0
       del0(Nx, j, 1) = (q0(1, j+1, 1) + 2.0*q0(Nx, j+1, 1) + q0(Nx, j, 1) + 2.0*q0(1, j, 1))/24.0
       bet0(1, j, 1) = (q0(1, j+1, 1) + 2.0*q0(Nx, j+1, 1) + q0(Nx, j, 1) + 2.0*q0(1, j, 1))/24.0
       gam0(1, j, 1) = (2.0*q0(1, j+1, 1) + q0(Nx, j+1, 1) + 2.0*q0(Nx, j, 1) + q0(1, j, 1))/24.0
      end do

      do i = 2, Nx
       do j = 1, Ny-1
        bet0(i, j, 1) = (q0(i, j+1, 1) + 2.0*q0(i-1, j+1, 1) + q0(i-1, j, 1) + 2.0*q0(i, j, 1))/24.0
        gam0(i, j, 1) = (2.0*q0(i, j+1, 1) + q0(i-1, j+1, 1) + 2.0*q0(i-1, j, 1) + q0(i, j, 1))/24.0
       end do
      end do

      !create the initial data (vor, pe0, ke0) files for each resolutions!
      open(14, file='vor_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(14, rec=1)vor0
      close(14)
      
      open(15, file='pe0_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(15, rec=1)pe0
      close(15)

      open(16, file='ke0_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(16, rec=1)ke0
      close(16)
      
      !!!============Forward (Euler) scheme============!!!
      do n = 2, 3

      !h update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        h(i, j) = h(i, j) - (dt*(us0(i+1, j, n-1) - us0(i, j, n-1) + vs0(i, j+1, n-1) - vs0(i, j, n-1)))/d   
       end do
      end do
      
      !hu0, hv0, us0, vs0 update!
      hu0(1, :, n) = (h(Nx, :) + h(2, :))/2.0
      hu0(Nx, :, n) = (h(Nx-1, :) + h(1, :))/2.0

      do i = 2, Nx-1
       hu0(i, :, n) = (h(i-1, :) + h(i+1, :))/2.0
      end do

      hv0(:, 1, n) = (h(:, Ny) + h(:, 2))/2.0
      hv0(:, Ny, n) = (h(:, Ny-1) + h(:, 1))/2.0

      do j = 2, Ny-1
       hv0(:, j, n) = (h(:, j-1) + h(:, j+1))/2.0
      end do

      us0(:, :, n) = hu0(:, :, n) * u(:, :)
      vs0(:, :, n) = hv0(:, :, n) * v(:, :)

      !u update!
      do i = 2, Nx-1
       do j = 1, Ny-1
        u(i, j) = u(i, j) + dt*(alp0(i, j, n-1)*vs0(i, j+1, n-1) + bet0(i, j, n-1)*vs0(i-1, j+1, n-1) + gam0(i, j, n-1)*vs0(i-1, j, n-1) & 
                            + del0(i, j, n-1)*vs0(i, j, n-1) - eps0(i, j, n-1)*us0(i+1, j, n-1) + eps0(i-1, j, n-1)*us0(i-1, j, n-1)) &
                            - (dt*(ke0(i, j, n-1) + pe0(i, j, n-1) - ke0(i-1, j, n-1) - pe0(i-1, j, n-1)))/d
       end do
      end do

      do j = 1, Ny-1
       u(1, j) = u(1, j) + dt*(alp0(1, j, n-1)*vs0(1, j+1, n-1) + bet0(1, j, n-1)*vs0(Nx-1, j+1, n-1) + gam0(1, j, n-1)*vs0(Nx-1, j, n-1) & 
                            + del0(1, j, n-1)*vs0(1, j, n-1) - eps0(1, j, n-1)*us0(2, j, n-1) + eps0(Nx-1, j, n-1)*us0(Nx, j, n-1)) &
                            - (dt*(ke0(1, j, n-1) + pe0(1, j, n-1) - ke0(Nx-1, j, n-1) - pe0(Nx-1, j, n-1)))/d
       u(Nx, j) = u(Nx, j) + dt*(alp0(Nx, j, n-1)*vs0(1, j+1, n-1) + bet0(Nx, j, n-1)*vs0(Nx-1, j+1, n-1) + gam0(Nx, j, n-1)*vs0(Nx-1, j, n-1) & 
                            + del0(Nx, j, n-1)*vs0(1, j, n-1) - eps0(1, j, n-1)*us0(1, j, n-1) + eps0(Nx-1, j, n-1)*us0(Nx-1, j, n-1)) &
                            - (dt*(ke0(1, j, n-1) + pe0(1, j, n-1) - ke0(Nx-1, j, n-1) - pe0(Nx-1, j, n-1)))/d
      end do

      !v update!
      do i = 1, Nx-1
       do j = 1, Ny
        v(i, j) = v(i, j) - dt*(gam0(i+1, j, n-1)*us0(i+1, j, n-1) + del0(i, j, n-1)*us0(i, j, n-1) + alp0(i, j-1, n-1)*us0(i, j-1, n-1) &
                              + bet0(i+1, j-1, n-1)*us0(i+1, j-1, n-1) + pi0(i, j, n-1)*vs0(i, j+1, n-1) - pi0(i, j-1, n-1)*vs0(i, j-1, n-1)) &
                              -(dt*(ke0(i, j, n-1) + pe0(i, j, n-1) - ke0(i, j-1, n-1) - pe0(i, j-1, n-1)))/d
       end do
      end do
      
      !vorticity update!
      vor0(1, 1, n) = (-u(1, 1) + v(1, 1) - v(Nx-1, 1))/d
      vor0(1, Ny, n) = (u(1, Ny-1) + v(1, Ny) - v(Nx-1, Ny))/d
      vor0(Nx, 1, n) = (-u(Nx, 1) + v(1, 1) - v(Nx-1, 1))/d
      vor0(Nx, Ny, n) = (u(Nx, Ny-1) + v(1, Ny) - v(Nx-1, Ny))/d

      do j = 2, Ny-1
       vor0(1, j, n) = (u(1, j-1) - u(1, j) + v(1, j) - v(Nx-1, j))/d
       vor0(Nx, j, n) = (u(Nx, j-1) - u(Nx, j) + v(1, j) - v(Nx-1, j))/d
      end do
      
      do i = 2, Nx-1
       vor0(i, 1, n) = (-u(i, 1) + v(i, 1) - v(i-1, 1))/d
       vor0(i, Ny, n) = (u(i, Ny-1) + v(i, Ny) - v(i-1, Ny))/d
      end do
      
      do i = 2, Nx-1
       do j = 2, Ny-1
        vor0(i, j, n) = (u(i, j-1) - u(i, j) + v(i, j) - v(i-1, j))/d
       end do
      end do

      !hq0 update!
      hq0(1,1,n) = (h(1,1) + h(Nx-1,1))/2.0
      hq0(1,Ny,n) = (h(1,Ny-1) + h(Nx-1,Ny-1))/2.0
      hq0(Nx,1,n) = (h(Nx-1,1) + h(1,1))/2.0
      hq0(Nx,Ny,n) = (h(1,Ny-1) + h(Nx-1,Ny-1))/2.0

      do i = 2, Nx-1
       hq0(i,1,n) = (h(i,1) + h(i-1, 1))/2.0
       hq0(i,Ny,n) = (h(i,Ny-1) + h(i-1,Ny-1))/2.0
      end do

      do j = 2, Ny-1
       hq0(1,j,n) = (h(1,j) + h(Nx-1,j) + h(Nx-1, j-1) + h(1,j-1))/4.0
       hq0(Nx,j,n) = (h(1,j) + h(Nx-1,j) + h(Nx-1, j-1) + h(1,j-1))/4.0
      end do

      do j = 2, Ny-1
       do i = 2, Nx-1
        hq0(i,j,n) = (h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
       end do
      end do

      !q0 update!
      q0(:, :, n) = (vor0(:, :, n) + f)/hq0(:, :, n)

      !pe0 update!
      do i = 1, Nx
       pe0(i, :, n) = g*(h(i, :) + hs(i, :))
      end do

      !ke0 update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        ke0(i, j, n) = (u(i, j)**2 + u(i+1, j)**2 + v(i, j)**2 + v(i, j+1)**2)/4.0
       end do
      end do

      !alp0, bet0, del0, gam0, eps0, pi0 update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        alp0(i, j, n) = (2.0*q0(i+1, j+1, n) + q0(i, j+1, n) + 2.0*q0(i, j, n) + q0(i+1, j, n))/24.0
        del0(i, j, n) = (q0(i+1, j+1, n) + 2.0*q0(i, j+1, n) + q0(i, j, n) + 2.0*q0(i+1, j, n))/24.0
        eps0(i, j, n) = (q0(i+1, j+1, n) + q0(i, j+1, n) - q0(i, j, n) - q0(i+1, j, n))/24.0
        pi0(i, j, n) = (-q0(i+1, j+1, n) + q0(i, j+1, n) + q0(i, j, n) - q0(i+1, j, n))/24.0
       end do
      end do

      do j = 1, Ny-1
       bet0(1, j, n) = (q0(1, j+1, n) + 2.0*q0(Nx, j+1, n) + q0(Nx, j, n) + 2.0*q0(1, j, n))/24.0
       gam0(1, j, n) = (2.0*q0(1, j+1, n) + q0(Nx, j+1, n) + 2.0*q0(Nx, j, n) + q0(1, j, n))/24.0
      end do

      do i = 2, Nx
       do j = 1, Ny-1
        bet0(i, j, n) = (q0(i, j+1, n) + 2.0*q0(i-1, j+1, n) + q0(i-1, j, n) + 2.0*q0(i, j, n))/24.0
        gam0(i, j, n) = (2.0*q0(i, j+1, n) + q0(i-1, j+1, n) + 2.0*q0(i-1, j, n) + q0(i, j, n))/24.0
       end do
      end do
      end do
      
      !!!============Adams-bashforth third order scheme============!!!
      !allocate variables for Adams-bashforth scheme!
      allocate(us1(Nx, Ny), us2(Nx, Ny), us3(Nx, Ny), vs1(Nx, Ny), vs2(Nx, Ny), vs3(Nx, Ny), hu(Nx, Ny), hv(Nx, Ny))
      allocate(alp1(Nx, Ny), alp2(Nx, Ny), alp3(Nx, Ny), bet1(Nx, Ny), bet2(Nx, Ny), bet3(Nx, Ny), gam1(Nx, Ny), gam2(Nx, Ny), gam3(Nx, Ny))
      allocate(del1(Nx, Ny), del2(Nx, Ny), del3(Nx, Ny), eps1(Nx, Ny), eps2(Nx, Ny), eps3(Nx, Ny), pi1(Nx, Ny), pi2(Nx, Ny), pi3(Nx, Ny), pe1(Nx, Ny), pe2(Nx, Ny), pe3(Nx, Ny), ke1(Nx, Ny), ke2(Nx, Ny), ke3(Nx, Ny))
      allocate(hq(Nx, Ny), vor(Nx, Ny), us(Nx, Ny), vs(Nx, Ny), q(Nx, Ny), pe(Nx, Ny), ke(Nx, Ny))
      allocate(alp(Nx,Ny), bet(Nx,Ny), del(Nx,Ny), gam(Nx,Ny), eps(Nx,Ny), pi(Nx,Ny))
      
      us1(:,:) = us0(:,:,1)
      us2(:,:) = us0(:,:,2)
      us3(:,:) = us0(:,:,3)
      vs1(:,:) = vs0(:,:,1)
      vs2(:,:) = vs0(:,:,2)
      vs3(:,:) = vs0(:,:,3)
      alp1(:,:) = alp0(:,:,1)
      alp2(:,:) = alp0(:,:,2)
      alp3(:,:) = alp0(:,:,3)
      bet1(:,:) = bet0(:,:,1)
      bet2(:,:) = bet0(:,:,2)
      bet3(:,:) = bet0(:,:,3)
      gam1(:,:) = gam0(:,:,1)
      gam2(:,:) = gam0(:,:,2)
      gam3(:,:) = gam0(:,:,3)
      del1(:,:) = del0(:,:,1)
      del2(:,:) = del0(:,:,2)
      del3(:,:) = del0(:,:,3)
      eps1(:,:) = eps0(:,:,1)
      eps2(:,:) = eps0(:,:,2)
      eps3(:,:) = eps0(:,:,3)
      pi1(:,:) = pi0(:,:,1)
      pi2(:,:) = pi0(:,:,2)
      pi3(:,:) = pi0(:,:,3)
      ke1(:,:) = ke0(:,:,1)
      ke2(:,:) = ke0(:,:,2)
      ke3(:,:) = ke0(:,:,3)
      pe1(:,:) = pe0(:,:,1)
      pe2(:,:) = pe0(:,:,2)
      pe3(:,:) = pe0(:,:,3)

      nstep = 4
      do n = 4, ntime
       nstep = nstep + 1
      !h update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        h(i, j) = h(i, j) - (f1*dt*(us1(i+1, j) - us1(i, j) + vs1(i, j+1) - vs1(i, j)))/d - &
                           (f2*dt*(us2(i+1, j) - us2(i, j) + vs2(i, j+1) - vs2(i, j)))/d - &
                           (f3*dt*(us3(i+1, j) - us3(i, j) + vs3(i, j+1) - vs3(i, j)))/d
       end do
      end do

      !hu, hv, us, vs update!
      do i = 2, Nx-1
       hu(i, :) = (h(i-1, :) + h(i, :))/2.0
      end do

      do j = 2, Ny-1
       hv(:, j) = (h(:, j-1) + h(:, j))/2.0
      end do

      us(:, :) = hu(:, :) * u(:, :)
      vs(:, :) = hv(:, :) * v(:, :)

      !u update!
      do i = 2, Nx-1
       do j = 1, Ny-1
        u(i, j) = u(i, j) + f1*(dt*(alp1(i, j)*vs1(i, j+1) + bet1(i, j)*vs1(i-1, j+1) + gam1(i, j)*vs1(i-1, j) & 
                            + del1(i, j)*vs1(i, j) - eps1(i, j)*us1(i+1, j) + eps1(i-1, j)*us1(i-1, j)) &
                            - (dt*(ke1(i, j) + pe1(i, j) - ke1(i-1, j) - pe1(i-1, j)))/d) &
                          + f2*(dt*(alp2(i, j)*vs2(i, j+1) + bet2(i, j)*vs2(i-1, j+1) + gam2(i, j)*vs2(i-1, j) & 
                            + del2(i, j)*vs2(i, j) - eps2(i, j)*us2(i+1, j) + eps2(i-1, j)*us2(i-1, j)) &
                            - (dt*(ke2(i, j) + pe2(i, j) - ke2(i-1, j) - pe2(i-1, j)))/d) &
                          + f3*(dt*(alp3(i, j)*vs3(i, j+1) + bet3(i, j)*vs3(i-1, j+1) + gam3(i, j)*vs3(i-1, j) & 
                            + del3(i, j)*vs3(i, j) - eps3(i, j)*us3(i+1, j) + eps3(i-1, j)*us3(i-1, j)) &
                            - (dt*(ke3(i, j) + pe3(i, j) - ke3(i-1, j) - pe3(i-1, j)))/d)
       end do
      end do

      do j = 1, Ny-1
       u(1, j) = u(1, j) + f1*(dt*(alp1(1, j)*vs1(1, j+1) + bet1(1, j)*vs1(Nx-1, j+1) + gam1(1, j)*vs1(Nx-1, j) & 
                            + del1(1, j)*vs1(1, j) - eps1(1, j)*us1(2, j) + eps1(Nx-1, j)*us1(Nx, j)) &
                            - (dt*(ke1(1, j) + pe1(1, j) - ke1(Nx-1, j) - pe1(Nx-1, j)))/d) &
                         + f2*(dt*(alp2(1, j)*vs2(1, j+1) + bet2(1, j)*vs2(Nx-1, j+1) + gam2(1, j)*vs2(Nx-1, j) & 
                            + del2(1, j)*vs2(1, j) - eps2(1, j)*us2(2, j) + eps2(Nx-1, j)*us2(Nx, j)) &
                            - (dt*(ke2(1, j) + pe2(1, j) - ke2(Nx-1, j) - pe2(Nx-1, j)))/d) &
                         + f3*(dt*(alp3(1, j)*vs3(1, j+1) + bet3(1, j)*vs3(Nx-1, j+1) + gam3(1, j)*vs3(Nx-1, j) & 
                            + del3(1, j)*vs3(1, j) - eps3(1, j)*us3(2, j) + eps3(Nx-1, j)*us3(Nx, j)) &
                            - (dt*(ke3(1, j) + pe3(1, j) - ke3(Nx-1, j) - pe3(Nx-1, j)))/d)
                            
       u(Nx, j) = u(Nx, j) + f1*(dt*(alp1(Nx, j)*vs1(1, j+1) + bet1(Nx, j)*vs1(Nx-1, j+1) + gam1(Nx, j)*vs1(Nx-1, j) & 
                             + del1(Nx, j)*vs1(1, j) - eps1(1, j)*us1(1, j) + eps1(Nx-1, j)*us1(Nx-1, j)) &
                             - (dt*(ke1(1, j) + pe1(1, j) - ke1(Nx-1, j) - pe1(Nx-1, j)))/d) &
                           + f2*(dt*(alp2(Nx, j)*vs2(1, j+1) + bet2(Nx, j)*vs2(Nx-1, j+1) + gam2(Nx, j)*vs2(Nx-1, j) & 
                             + del2(Nx, j)*vs2(1, j) - eps2(1, j)*us2(1, j) + eps2(Nx-1, j)*us2(Nx-1, j)) &
                             - (dt*(ke2(1, j) + pe2(1, j) - ke2(Nx-1, j) - pe2(Nx-1, j)))/d) &
                           + f3*(dt*(alp3(Nx, j)*vs3(1, j+1) + bet3(Nx, j)*vs3(Nx-1, j+1) + gam3(Nx, j)*vs3(Nx-1, j) & 
                             + del3(Nx, j)*vs3(1, j) - eps3(1, j)*us3(1, j) + eps3(Nx-1, j)*us3(Nx-1, j)) &
                             - (dt*(ke3(1, j) + pe3(1, j) - ke3(Nx-1, j) - pe3(Nx-1, j)))/d)  
      end do

      !v update!
      do i = 1, Nx-1
       do j = 2, Ny-1
        v(i, j) = v(i, j) - f1*(dt*(gam1(i+1, j)*us1(i+1, j) + del1(i, j)*us1(i, j) + alp1(i, j-1)*us1(i, j-1) &
                              + bet1(i+1, j-1)*us1(i+1, j-1) + pi1(i, j)*vs1(i, j+1) - pi1(i, j-1)*vs1(i, j-1)) &
                              -(dt*(ke1(i, j) + pe1(i, j) - ke1(i, j-1) - pe1(i, j-1)))/d) &
                          - f2*(dt*(gam2(i+1, j)*us2(i+1, j) + del2(i, j)*us2(i, j) + alp2(i, j-1)*us2(i, j-1) &
                              + bet2(i+1, j-1)*us2(i+1, j-1) + pi2(i, j)*vs2(i, j+1) - pi2(i, j-1)*vs2(i, j-1)) &
                              -(dt*(ke2(i, j) + pe2(i, j) - ke2(i, j-1) - pe2(i, j-1)))/d) &
                          - f3*(dt*(gam3(i+1, j)*us3(i+1, j) + del3(i, j)*us3(i, j) + alp3(i, j-1)*us3(i, j-1) &
                              + bet3(i+1, j-1)*us3(i+1, j-1) + pi3(i, j)*vs3(i, j+1) - pi3(i, j-1)*vs3(i, j-1)) &
                              -(dt*(ke3(i, j) + pe3(i, j) - ke3(i, j-1) - pe3(i, j-1)))/d)
       end do
      end do

      !vorticity update!
      vor(1, 1) = (-u(1, 1) + v(1, 1) - v(Nx-1, 1))/d
      vor(1, Ny) = (u(1, Ny-1) + v(1, Ny) - v(Nx-1, Ny))/d
      vor(Nx, 1) = (-u(Nx, 1) + v(1, 1) - v(Nx-1, 1))/d
      vor(Nx, Ny) = (u(Nx, Ny-1) + v(1, Ny) - v(Nx-1, Ny))/d

      do j = 2, Ny-1
       vor(1, j) = (u(1, j-1) - u(1, j) + v(1, j) - v(Nx-1, j))/d
       vor(Nx, j) = (u(Nx, j-1) - u(Nx, j) + v(1, j) - v(Nx-1, j))/d
      end do
      
      do i = 2, Nx-1
       vor(i, 1) = (-u(i, 1) + v(i, 1) - v(i-1, 1))/d
       vor(i, Ny) = (u(i, Ny-1) + v(i, Ny) - v(i-1, Ny))/d
      end do
      
      do i = 2, Nx-1
       do j = 2, Ny-1
        vor(i, j) = (u(i, j-1) - u(i, j) + v(i, j) - v(i-1, j))/d
       end do
      end do
      
      !hq update!
      hq(1,1) = (h(1,1) + h(Nx-1,1))/2.0
      hq(1,Ny) = (h(1,Ny-1) + h(Nx-1,Ny-1))/2.0
      hq(Nx,1) = (h(Nx-1,1) + h(1,1))/2.0
      hq(Nx,Ny) = (h(1,Ny-1) + h(Nx-1,Ny-1))/2.0

      do i = 2, Nx-1
       hq(i,1) = (h(i,1) + h(i-1, 1))/2.0
       hq(i,Ny) = (h(i,Ny-1) + h(i-1,Ny-1))/2.0
      end do

      do j = 2, Ny-1
       hq(1,j) = (h(1,j) + h(Nx-1,j) + h(Nx-1, j-1) + h(1,j-1))/4.0
       hq(Nx,j) = (h(1,j) + h(Nx-1,j) + h(Nx-1, j-1) + h(1,j-1))/4.0
      end do

      do j = 2, Ny-1
       do i = 2, Nx-1
        hq(i,j) = (h(i,j) + h(i-1,j) + h(i-1,j-1) + h(i,j-1))/4.0
       end do
      end do

      !q update!
      q(:, :) = (vor(:, :) + f)/hq(:, :)

      !pe update!
      do i = 1, Nx
       pe(i, :) = g*(h(i, :) + hs(i, :))
      end do

      !ke update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        ke(i, j) = (u(i, j)**2 + u(i+1, j)**2 + v(i, j)**2 + v(i, j+1)**2)/4.0
       end do
      end do

      !update alp, bet, del, gam, eps, pi forcing coefficients!
      do i = 1, Nx-1
       do j = 1, Ny-1
        alp(i, j) = (2.0*q(i+1, j+1) + q(i, j+1) + 2.0*q(i, j) + q(i+1, j))/24.0
        del(i, j) = (q(i+1, j+1) + 2.0*q(i, j+1) + q(i, j) + 2.0*q(i+1, j))/24.0
        eps(i, j) = (q(i+1, j+1) + q(i, j+1) - q(i, j) - q(i+1, j))/24.0
        pi(i, j) = (-q(i+1, j+1) + q(i, j+1) + q(i, j) - q(i+1, j))/24.0
       end do
      end do

      do j = 1, Ny-1
       bet(1, j) = (q(1, j+1) + 2.0*q(Nx, j+1) + q(Nx, j) + 2.0*q(1, j))/24.0
       gam(1, j) = (2.0*q(1, j+1) + q(Nx, j+1) + 2.0*q(Nx, j) + q(1, j))/24.0
      end do

      do i = 2, Nx
       do j = 1, Ny-1
        bet(i, j) = (q(i, j+1) + 2.0*q(i-1, j+1) + q(i-1, j) + 2.0*q(i, j))/24.0
        gam(i, j) = (2.0*q(i, j+1) + q(i-1, j+1) + 2.0*q(i-1, j) + q(i, j))/24.0
       end do
      end do

      us1(:, :) = us2(:, :)
      us2(:, :) = us3(:, :)
      us3(:, :) = us(:, :)
      vs1(:, :) = vs2(:, :)
      vs2(:, :) = vs3(:, :)
      vs3(:, :) = vs(:, :)
      alp1(:, :) = alp2(:, :)
      alp2(:, :) = alp3(:, :)
      alp3(:, :) = alp(:, :)
      bet1(:, :) = bet2(:, :)
      bet2(:, :) = bet3(:, :)
      bet3(:, :) = bet(:, :)
      gam1(:, :) = gam2(:, :)
      gam2(:, :) = gam3(:, :)
      gam3(:, :) = gam(:, :)
      del1(:, :) = del2(:, :)
      del2(:, :) = del3(:, :)
      del3(:, :) = del(:, :)
      eps1(:, :) = eps2(:, :)
      eps2(:, :) = eps3(:, :)
      eps3(:, :) = eps(:, :)
      pi1(:, :) = pi2(:, :)
      pi2(:, :) = pi3(:, :)
      pi3(:, :) = pi(:, :)
      pe1(:, :) = pe2(:, :)
      pe2(:, :) = pe3(:, :)
      pe3(:, :) = pe(:, :)
      ke1(:, :) = ke2(:, :)
      ke2(:, :) = ke3(:, :)
      ke3(:, :) = ke(:, :)
      
      !create the final data (u, v, h) files for each resolutions!
      if (nstep == 1440) then
       open(20, file='u_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(20, rec=1)u
       close(20)

       open(21, file='v_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(21, rec=1)v
       close(21)

       open(22, file='h_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(22, rec=1)h
       close(22)

       open(23, file='pe_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(23, rec=1)pe
       close(23)

       open(24, file='ke_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(24, rec=1)ke
       close(24)

       open(25, file='vor_final_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
       write(25, rec=1)vor
       close(25)
      endif
      end do

end program shallow_water_model
