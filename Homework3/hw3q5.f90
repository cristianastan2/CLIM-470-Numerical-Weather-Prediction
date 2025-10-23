program hw3q5

      !just leaving this here
      !Du/Dt = du/dt + u(du/dx) + v(du/dy) + w(du/dz)
      != fv - (1/rhoa)(dpa/dx) + Fr
      !Note Fr=0, du/dy=0, du/dz=0
      !so really Du/Dt = du/dt + u(du/dx) = fv - (1/rhoa)(dpa/dx)
      !and fv-(1/rhoa)(dpa/dx) is constant
      !and du/dx is constant
      !so du/dt = 0.065 -u/5000 

      implicit none
      real::u,dudt
      integer::i
      u=3.5 !initial zonal wind

      open(unit=40,file='zonalwind.dat',status='old',action='write')
      write(40,'(2F9.4)') u

      do i = 1,3600,1
      dudt = 0.065 -(u/5000)
      u = u + dudt
      
      write(40,'(2F9.4)') u

      end do

      close(40)

      end program hw3q5
