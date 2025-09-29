program triangle

real :: a, b, c, theta

write(*,*) 'Enter the length of hypoteneuse C:'
read(*,*) c
write(*,*) 'Enter the angle THETA in degrees:'
read(*,*) theta

a=c*cos(theta)
b=c*sin(theta)

write(*,*) 'The length of adjacent side is ', a
write(*,*) 'The length of opposite side is ', b

stop
end program triangle
