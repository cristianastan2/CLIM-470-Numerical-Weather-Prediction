!   Purpose:
!   To illustrate some of the basic features of a Fortran program
!

 program example2 
!   Declare the variables used in this program
INTEGER:: i, j, k		! All variables are integers

!   Get two values to store variables i and j
write(*,*) 'Enter numbers to multiply:'
read(*,*) i, j

! Multiply the numbers together
k = i * j

stop

! Write out the result
write(*,*) 'Result=  ', k

!stop
end program example2


