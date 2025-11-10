program initial_conditions
      implicit none
      
      !parameters! 
      integer, parameter::Lx = 6e+06 !domain size in x-direction
      integer, parameter::Ly = 2e+06 !domain size in y-direction
      real, parameter::hs_top = 2e+03 !height of the mountain
      integer, parameter::t = 1 !time step (s)
      
      !resolution!
      real:: d !model resolution, ie delta_x, delta_y

      !variables!
      integer:: Nx !number of grid points in x-direction (13, 25, 49)
      integer:: Ny !number of grid points in y-direction (5, 9, 17)
      integer:: i, j, ii, jj, ierr !working variables
      real, allocatable::hs(:,:), u(:,:), v(:,:), h(:,:)
      real, allocatable::hu0(:,:,3), hv0(:,:,3), us0(:,:,3),vs0(:,:,3) !store the results from forward scheme for each time steps

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
      u(1:Nx, 1:Ny) = 20 !m/s, initial condition on reference paper

      !allocate v variable!

      !allocate h variable!
      allocate(h(Nx,Ny))
      do i = 1,Nx
      h(i:) = 5e+03 !in m, initial height "hzero" defined in Arakawa and Lamb 1981
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

end program initial_conditions

