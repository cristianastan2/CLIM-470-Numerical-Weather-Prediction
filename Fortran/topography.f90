program topography
implicit none

!parameters

real, parameter :: Lx = 6.06e6
real, parameter :: Ly = 2.02e6
real, parameter :: d = 5.05e5
real, parameter :: hs_t = 2.0e3

!grid dimensions

integer, parameter :: Nx = int(Lx/d) + 1
integer, parameter :: Ny = int(Ly/d) + 1

!variables

real :: hs(Nx)
integer :: i

!init topo

hs = 0.0

!calculate topo based on res

if (d == 5.05e5) then
hs(Nx/2) = hs_t

else if (d == 2.5e5) then
hs(Nx/2-1) = 1.0e3
hs(Nx/2) = hs_t
hs(Nx/2+1) = 1.0e3

else if (d == 1.25e5) then
hs(Nx/2-2) = 0.5e3
hs(Nx/2-1) = 1.5e3
hs(Nx/2) = hs_t
hs(Nx/2+1) = 1.5e3
hs(Nx/2+2) = 0.5e3

end if

!results

print *, "Nx =", Nx
print *, "Ny =", Ny
print *, "Topography heights in meters:"

do i = 1, Nx
print *, i, hs(i)
end do

end program topography
