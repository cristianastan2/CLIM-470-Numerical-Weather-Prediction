program triangle

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
