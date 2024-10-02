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

end program triangle
