program meridional

! initialize all variables
implicit none
integer, parameter :: Nx = 13 + 2
integer, parameter :: Ny = 5 + 2
integer :: i,j
real, dimension(Nx, Ny) :: q
real, dimension(Nx, Ny) :: alpha
real, dimension(Nx, Ny) :: beta
real, dimension(Nx, Ny) :: gamma
real, dimension(Nx, Ny) :: delta
real, dimension(Nx, Ny) :: epsilon
real, dimension(Nx, Ny) :: phi
real, dimension(Nx, Ny) :: u
real, dimension(Nx, Ny) :: v
real, dimension(Nx, Ny) :: h
real, dimension(Nx, Ny) :: hsurf

! make all arrays all zeros
q = 0
alpha = 0
beta = 0
gamma = 0
delta = 0
epsilon = 0
phi = 0

! make q grid all 2s just as an example
do i=2, Nx - 1
    do j=2, Ny - 1
        q(i,j) = 2
    end do
end do

! set first column to second to last column and last column to second column (for q) 
q(1,:) = q(Nx - 1,:)
q(Nx,:) = q(2,:)

! calculate grid values for greek letters 
do i=2, Nx - 1
    do j=2, Ny - 1
        alpha(i,j) = ((1/24)*(2*q(i+1,j+1) + q(i,j+1) + 2*q(i,j) + q(i+1,j)))
        beta(i,j) = ((1/24)*(q(i,j+1) + 2*q(i-1,j+1) + q(i-1,j) + 2*q(i,j)))
        gamma(i,j) = ((1/24)*(2*q(i,j+1) + q(i-1,j+1) + 2*q(i-1,j) + q(i,j)))
        delta(i,j) = ((1/24)*(q(i+1,j+1) + 2*q(i,j+1) + q(i,j) + 2*q(i+1,j)))
        epsilon(i,j) = ((1/24)*(q(i+1,j+1) + q(i,j+1) - q(i,j) - q(i+1,j)))
        phi(i,j) = ((1/24)*(-q(i+1,j+1) + q(i,j+1) + q(i,j) - q(i+1,j))) 
    end do
end do 

! set first columns to second to last columns and last columns to second columns
alpha(1,:) = alpha(Nx - 1,:)
alpha(Nx,:) = alpha(2,:)
beta(1,:) = beta(Nx - 1,:)
beta(Nx,:) = beta(2,:)
gamma(1,:) = gamma(Nx - 1,:)
gamma(Nx,:) = gamma(2,:)
delta(1,:) = delta(Nx - 1,:)
delta(Nx,:) = delta(2,:)
epsilon(1,:) = epsilon(Nx - 1,:)
epsilon(Nx,:) = epsilon(2,:)
phi(1,:) = phi(Nx - 1,:)
phi(Nx,:) = phi(2,:)


call initgrid(h, u, v, q, hsurf)

end program meridional

subroutine initgrid(h, u, v, q, hsurf)	!creates grid, ridge, and topography.
    implicit none
  
    real, parameter :: Lx = 6000, Ly = 2000
    real, parameter :: Dx = 500
    integer, parameter :: Nx = (Lx/Dx) + 3, Ny = (Ly/Dx) + 3

    real, dimension(Nx,Ny), intent(out) :: h, u, v, q, hsurf
    integer :: i, j

    hsurf(:,:)=0

    !grid resolution thing (topography)

    hsurf((Nx)/2 + 1,:) = 2000
    if (Dx==250) then
        hsurf((Nx)/2,:) = 1000
        hsurf((Nx)/2 + 2,:) = 1000
    endif
    if (Dx==125) then
        hsurf((Nx)/2 - 1,:) = 500
        hsurf((Nx)/2 + 3,:) = 500
    endif
	
    !initial conditions set to 0
    h = 0
    u = 20 
    v = 0.5 
    q = 0

    u(:,Ny) = 0
    u(:,1) = 0
    v(:,Ny) = 0
    v(:,1) = 0

    end subroutine initgrid
