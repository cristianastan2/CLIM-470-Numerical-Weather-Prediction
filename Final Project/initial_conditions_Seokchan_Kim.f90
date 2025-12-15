program shallow_water_model
      implicit none
      
      !parameters! 
      integer, parameter::Lx = 6e+06 !domain size in x-direction (m)
      integer, parameter::Ly = 2e+06 !domain size in y-direction (m)
      real, parameter::hs_top = 2e+03 !height of the mountain (m)
      real, parameter::g = 9.8 !the acceleration of gravity (m/s^2)
      real, parameter::f = 1e-04 !the Coriolis parameter (s^-1)
      integer, parameter::t = 1 !time step (s)
      
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

      !allocate alp0, bet0, gam0, del0, eps0, pi0 variables!
      allocate(alp0(Nx, Ny, 3), bet0(Nx, Ny, 3), gam0(Nx, Ny, 3), del0(Nx, Ny, 3), eps0(Nx, Ny, 3), pi0(Nx, Ny, 3))
      alp0(:, :, 1) = 0.0
      bet0(:, :, 1) = 0.0
      gam0(:, :, 1) = 0.0
      del0(:, :, 1) = 0.0
      eps0(:, :, 1) = 0.0
      pi0(:, :, 1) = 0.0
      
      do i = 2, Nx-1
       do j = 2, Ny-1
        alp0(i, j, 1) = (2.0*q0(i+1, j+1, 1) + q0(i, j+1, 1) + 2.0*q0(i, j) + q0(i+1, j, 1))/24.0
        del0(i, j, 1) = (q0(i+1, j+1, 1) + 2.0*q0(i, j+1, 1) + q0(i, j) + 2.0*q0(i+1, j, 1))/24.0
        eps0(i, j, 1) = (q0(i+1, j+1, 1) + q0(i, j+1, 1) - q0(i, j) - q0(i+1, j, 1))/24.0
        pi0(i, j, 1) = (-q0(i+1, j+1, 1) + q0(i, j+1, 1) + q0(i, j) - q0(i+1, j, 1))/24.0
       !(should be modified)
       end do
      end do

      do j = 1, Ny-1
       alp0(Nx, j, 1) = (1/24) * (2*q0(2, j+1, 1) + q0(1, j+1, 1) + 2*q0(1, j) + q0(2, j, 1))
       del0(Nx, j, 1) = (1/24) * (q0(2, j+1, 1) + 2*q0(1, j+1, 1) + q0(1, j) + 2*q0(2, j, 1))
       eps0(Nx, j, 1) = 0
       pi0(Nx, j, 1) = 0
       bet0(1, j, 1) = (1/24) * (q0(1, j+1) + 2*q0(Nx, j+1) + q0(Nx, j) + 2*q0(1, j))
       gam0(1, j, 1) = (1/24) * (2*q0(1, j+1) + q0(Nx, j+1) + 2*q0(Nx, j) + q0(1, j))
      end do

      do i = 2, Nx
       do j = 1, Ny-1
        bet0(i, j, 1) = (1/24) * (q0(i, j+1) + 2*q0(i-1, j+1) + q0(i-1, j) + 2*q0(i, j))
        gam0(i, j, 1) = (1/24) * (2*q0(i, j+1) + q0(i-1, j+1) + 2*q0(i-1, j) + q0(i, j))
       end do
      end do

      !create the initial data (vor, pe0, ke0) files for each resolutions!
      open(14, file='vor_initial_high_res.dat', status='unknown',form='unformatted', action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(14, rec=1)vor
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

      !h, u, v update!
      do i = 1, Nx-1
       do j = 1, Ny-1
        do ii = 2, Nx-1
         do jj = 2, Ny-1
          h(ii, jj) = h(ii, jj) - (t*(us0(i+1, jj, n-1) - us0(i, jj, n-1) + vs0(ii, j+1, n-1) - vs0(ii, j, n-1)))/d
          u(i, jj) = u(i, jj) + t*(alp0(i, jj+1, n-1)*vs0(ii, j, n-1) + bet0(i, jj+1, n-1)*vs0(ii-1, j+1, n-1) + gam0(i, jj+1, n-1)*vs0(ii-1, j, n-1) + del0(i, jj+1, n-1)*vs0(ii+1, j, n-1) - eps0(ii+1, jj+1, n-1)*us0(i+1, jj+1, n-1) + eps0(ii-1, jj+1, n-1)*us0(i-1, jj+1, n-1)) - (t*(ke0(ii+1, jj+1, n-1) + pe0(ii+1, jj+1, n-1) - ke0(ii-1, jj+1, n-1) - pe0(ii-1, jj+1, n-1)))/d
          v(ii, j) = v(ii, j) - t*(gam0(i+1, jj+1, n-1)*us0(i+1, jj+1, n-1) + del0(i, jj+1, n-1)*us0(i, jj+1, n-1) + alp0(i, jj-1, n-1)*us0(i, jj-1, n-1) + bet0(i+1, jj-1, n-1)*us0(i+1, jj-1, n-1) + pi0(ii+1, jj+1, n-1)*vs0(ii+1, j+1, n-1) - pi0(ii+1, jj-1, n-1)*vs0(ii+1, jj-1, n-1) -(t*(ke0(ii+1, jj+1, n-1) + pe0(ii+1, jj+1, n-1) - ke0(ii+1, jj-1, n-1) - pe0(ii+1, jj-1, n-1)))/d
         end do
        end do
       end do 
      end do

      !hu0, hv0, us0, vs0 update!
      hu0(1, :, n) = (h(Nx-1, :) + h(1, :))/2.0

      do ii = 2, Nx-1
       hu0(ii, :, n) = (h(ii-1, :) + h(ii, :))/2.0
      end do

      hv0(:, 1, n) = (h(:, Ny-1) + h(:, 1))/2.0

      do jj = 2, Ny-1
       hv0(:, jj, n) = (h(:, jj-1) + h(:, jj))/2.0
      end do

      us0(:, :, n) = hu0(:, :, n) * u(:, :)
      vs0(:, :, n) = hv0(:, :, n) * v(:, :)

      !vorticity update!
      do i = 1, Nx
       do j = 1, Ny
        do ii = 2, Nx-1
         do jj = 2, Ny-1
          vor0(i, j, n) = (u(i, jj-1, n) - u(i, jj+1, n) + v(ii+1, j, n) - v(ii-1, j, n))/d
         end do
        end do
       end do
      end do

      !hq0 update!
      do i = 1, Nx
       do j = 1, Ny
        do ii = 2, Nx-1
         do jj = 2, Ny-1
          hq0(i, j, n) = (h(ii, jj) + h(ii-1, jj) + h(ii-1, jj-1) + h(ii, jj-1))/d
         end do
        end do
       end do
      end do

      !q0 update!
      q0(:, :, n) = (vor0(:, :, n) + f)/hq0(:, :, n)

      !pe0 update!
      do i = 1, Nx
       pe0(i, :, n) = g*(h(i, :) + hs(i, :))
      end do

      !ke0 update!
      ke0(:, :, n) = (u(:, :)*u(:, :) + v(:, :)*v(:, :))/2

      !alp0, bet0, del0, gam0, eps0, pi0 update!
      
         
      

      end do
      
      !set up for progressing Adams-Bashforth third order scheme!
      integer::ntime = 1440
      integer::nstep = 4
      
      real, allocatable::us1(:, :), us2(:, :), us3(:, :)
      real, allocatable::vs1(:, :), vs2(:, :), vs3(:, :)
      real, allocatable::hu1(:, :), hu2(:, :), hu3(:, :)
      real, allocatable::hv1(:, :), hv2(:, :), hv3(:, :)
      real, allocatable::alp1(:, :), alp2(:, :), alp3(:, :)
      real, allocatable::bet1(:, :), bet2(:, :), bet3(:, :)
      real, allocatable::gam1(:, :), gam2(:, :), gam3(:, :)
      real, allocatable::del1(:, :), del2(:, :), del3(:, :)
      real, allocatable::eps1(:, :), eps2(:, :), eps3(:, :)
      real, allocatable::pi1(:, :), pi2(:, :), pi3(:, :)
      real, allocatable::pe1(:, :), pe2(:, :), pe3(:, :)
      real, allocatable::ke1(:, :), ke2(:, :), ke3(:, :)

      allocate(us1(Nx, Ny), us2(Nx, Ny), us3(Nx, Ny), vs1(Nx, Ny), vs2(Nx, Ny), vs3(Nx, Ny), hu1(Nx, Ny), hu2(Nx, Ny), hu3(Nx, Ny), hv1(Nx, Ny), hv2(Nx, Ny), hv3(Nx, Ny))
      allocate(alp1(Nx, Ny), alp2(Nx, Ny), alp3(Nx, Ny), bet1(Nx, Ny), bet2(Nx, Ny), bet3(Nx, Ny), gam1(Nx, Ny), gam2(Nx, Ny), gam3(Nx, Ny))
      allocate(del1(Nx, Ny), del2(Nx, Ny), del3(Nx, Ny), eps1(Nx, Ny), eps2(Nx, Ny), eps3(Nx, Ny), pi1(Nx, Ny), pi2(Nx, Ny), pi3(Nx, Ny), pe1(Nx, Ny), pe2(Nx, Ny), pe3(Nx, Ny), ke1(Nx, Ny), ke2(Nx, Ny), ke3(Nx, Ny))

      us1(:,:) = us0(:,:,1)
      us2(:,:) = us0(:,:,2)
      us3(:,:) = us0(:,:,3)
      vs1(:,:) = vs0(:,:,1)
      vs2(:,:) = vs0(:,:,2)
      vs3(:,:) = vs0(:,:,3)
      hu1(:,:) = hu0(:,:,1)
      hu2(:,:) = hu0(:,:,2)
      hu3(:,:) = hu0(:,:,3)
      hv1(:,:) = hv0(:,:,1)
      hv2(:,:) = hv0(:,:,2)
      hv3(:,:) = hv0(:,:,3)
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

      !!!============Adams-bashforth third order scheme============!!!
      do n = 4, ntime
       nstep = nstep + 1
       do 

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
      endif
      
      
      end do

end program shallow_water_model

