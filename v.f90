program v

! initialize all variables
implicit none
integer, parameter :: Nx = 12
integer, parameter :: Ny = 4
integer :: i,j
real, dimension(Nx + 2, Ny + 2) :: q
real, dimension(Nx + 2, Ny + 2) :: alpha
real, dimension(Nx + 2, Ny + 2) :: beta
real, dimension(Nx + 2, Ny + 2) :: gamma
real, dimension(Nx + 2, Ny + 2) :: delta
real, dimension(Nx + 2, Ny + 2) :: epsilon
real, dimension(Nx + 2, Ny + 2) :: phi

! make all arrays all zeros
q = 0
alpha = 0
beta = 0
gamma = 0
delta = 0
epsilon = 0
phi = 0

! make q grid all 2s just as an example
do i=2, Nx + 1
    do j=2, Ny + 1
        q(i,j) = 2
    end do
end do

! set first column to second to last column and last column to second column (for q) 
q(1,:) = q(Nx + 1,:)
q(Nx + 2,:) = q(2,:)

! calculate grid values for greek letters 
do i=2, Nx + 1
    do j=2, Ny + 1
        alpha(i,j) = ((1.0/24.0)*(2*q(i+1,j+1) + q(i,j+1) + 2*q(i,j) + q(i+1,j)))
        beta(i,j) = ((1.0/24.0)*(q(i,j+1) + 2*q(i-1,j+1) + q(i-1,j) + 2*q(i,j)))
        gamma(i,j) = ((1.0/24.0)*(2*q(i,j+1) + q(i-1,j+1) + 2*q(i-1,j) + q(i,j)))
        delta(i,j) = ((1.0/24.0)*(q(i+1,j+1) + 2*q(i,j+1) + q(i,j) + 2*q(i+1,j)))
        epsilon(i,j) = ((1.0/24.0)*(q(i+1,j+1) + q(i,j+1) - q(i,j) - q(i+1,j)))
        phi(i,j) = ((1.0/24.0)*(-q(i+1,j+1) + q(i,j+1) + q(i,j) - q(i+1,j))) 
    end do
end do 

! set first columns to second to last columns and last columns to second columns
alpha(1,:) = alpha(Nx + 1,:)
alpha(Nx + 2,:) = alpha(2,:)
beta(1,:) = beta(Nx + 1,:)
beta(Nx + 2,:) = beta(2,:)
gamma(1,:) = gamma(Nx + 1,:)
gamma(Nx + 2,:) = gamma(2,:)
delta(1,:) = delta(Nx + 1,:)
delta(Nx + 2,:) = delta(2,:)
epsilon(1,:) = epsilon(Nx + 1,:)
epsilon(Nx + 2,:) = epsilon(2,:)
phi(1,:) = phi(Nx + 1,:)
phi(Nx + 2,:) = phi(2,:)

end program v
