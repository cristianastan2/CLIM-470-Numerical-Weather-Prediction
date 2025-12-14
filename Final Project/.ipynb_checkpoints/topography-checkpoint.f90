program topography
      implicit none

      integer, parameter::Lx = 6e+06 !domain size in x-direction
      integer, parameter::Ly = 2e+06 !domain size in y-direction
      real, parameter::hs_top = 2e+03 !height of the mountain

      real:: d !model resolution, ie delta_x, delta_y

      integer:: Nx !number of grid points in x-direction (13,25,49)
      integer:: Ny !number of grid points in y-direction (5,9,17)
      integer:: i, j ,ierr !working variables
      real, allocatable::hs(:,:)

      !d = 5e+05
      !d = 2.5e+05
      d = 1.25e+05

      Nx = Lx/d + 1
      Ny = Ly/d + 1

      write(*,100)'The number of points in x-direction for d=',d/1000,'km is',Nx
      write(*,100)'The number of points in y-direction for d=',d/1000,'km is',Ny

      100 format(A,f6.2,A,I3)

      allocate(hs(Nx,Ny))
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


      do i = 1, Nx
       do j = 1, 1
       print*,'hs(',i,',',j,')=',hs(i,j)
       enddo
      enddo

      open(10,file='topography_high_res.dat',status='unknown', form='unformatted',action='write',&
              access='direct',recl=4*Nx*Ny,iostat=ierr)
      write(10, rec=1)hs

      close(10)
      end program topography

      
