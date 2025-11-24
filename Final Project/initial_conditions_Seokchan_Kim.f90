program initial_conditions
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
      integer:: i, j, ii, jj, ierr !working variables
      real, allocatable::vor(:,:) !vorticity
      real, allocatable::q(:,:) !potential vorticity
      real, allocatable::pe(:,:) !potential energy
      real, allocatable::ke(:,:) !kinetic energy
      real, allocatable::hs(:,:), u(:,:), v(:,:), h(:,:)

      !store the results from forward scheme(Euler Scheme) for each time steps(n=2, n=3)
      real, allocatable::hu0(:,:,:), hv0(:,:,:), us0(:,:,:), vs0(:,:,:), hq0(:,:,:)
      real, allocatable::vor0(:,:,:), q0(:,:,:), pe(:,:,:), ke(:,:,:)
      real, allocatable::alp0(:,:,:), bet0(:,:,:), gam0(:,:,:), del0(:,:,:), eps(:,:,:) !constant

      !resolution!
      !d = 5e+05
      !d = 2.5e+05
      d = 1.25e+05

      Nx = Lx/d + 1
      Ny = Ly/d + 1

      !allocate hs variable!
      allocate(hs(Nx, Ny))
      hs(1:Nx, 1:Ny) = 0
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

      !allocate u variable!
      allocate(u(Nx,Ny))
      u(1:Nx, 1:Ny) = 20.0 !m/s, initial condition on reference paper

      !allocate v variable!
      allocate(v(Nx, Ny))
      v(1:Nx, 1) = 0 !rigid computational boundary
      v(1:Nx, Ny) = 0 !rigid computational boundary
      v(1:Nx, 2:Ny-1) = 10

      !allocate h variable!
      allocate(h(Nx,Ny))
      do i = 1,Nx
      h(i,:) = 5e+03 - hs(i, :)!in m, initial height "hzero" defined in Arakawa and Lamb 1981
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
      vs0(1:Nx, 1:Ny, 1) = vs0(1:Nx, 1:Ny, 1)*v(1:Nx, 1:Ny)

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
      vor0(i,1,1) = 0
      vor0(i,Ny,,1) = 0
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
       hq0(1,j,1) = (h(1,j) + h(i-1,j) + h(Nx, j) + h(Nx-1,j-1))/4.0
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
        q0(i,j,1) = (z0(i,j,1) + f)/hq0(i,j,1)
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
       ke0(i,Ny,1) = 0
      end do

      do j = 1, Ny-1
       ke0(Nx,j,1) = (u(Nx-1,j)**2 + u(1,j)**2 + v(Nx-1,j)**2 + v(1,Ny+1)**2)/4
      

end program initial_conditions

