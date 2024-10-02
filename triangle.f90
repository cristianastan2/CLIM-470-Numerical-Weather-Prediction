program triangle

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

end program triangle
