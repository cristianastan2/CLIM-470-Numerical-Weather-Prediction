!Rachel's Code
program triangle
implicit none

real::a,b,c,theta,thetar
real:: pi = 3.14

write(*,*)'Enter the length of the hypotenuse C:'
read(*,*)c
write(*,*)'Enter the angle theta in degrees'
read(*,*)theta

thetar = theta * 180 / pi

a=c*cos(thetar)
b=c*sin(thetar)

write(*,*)'The length of adjacent side is', a
write(*,*)'The length of opposite side is', b

! Alex's fortran code
program triangle
real:: a,b,c,theta

write(*,*) 'Enter the length of hypotenuse C:'
read(*,*) c
write(*,*)'Enter the angle THETA  in degrees'
read(*,*) theta

a=c*cos(theta)
b=c*sin(theta)

write(*,*) 'The length of adjacent side is ', a
write(*,*) 'The length of opposite side is ', b

! Matthew's code
real :: a,b,c,theta
real, parameter :: pi = 3.1415926536

write(*,*) "Enter the length of hypotenuse C:"
read(*,*) c
write(*,*) "Enter the angle theta in degrees:"
read(*,*) theta

theta = theta*180/pi

a = c*cos(theta)
b = c*sin(theta)

write(*,*) "The length of the adjacent side is",a
write(*,*) "The length of the opposite side is",b


! Ethan's code
real :: a, b, c, theta, degree
real, parameter :: pi = 3.14159265358979323846
write(*,*) 'enter the length of hypotenuse'
read(*,*) c

write(*,*) 'enter the angle in degrees'
read(*,*) theta

degree = theta * 180/pi

a=c * cos(degree)
b=c * sin(degree)

write(*,*) 'adjacent side is', a
write(*,*) 'opposite side is', b

stop


end program triangle
