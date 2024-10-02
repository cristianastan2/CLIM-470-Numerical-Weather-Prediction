program triangle

real:: a,b,c,theta


! Alex's fortran code
write(*,*) 'Enter the length of hypotenuse C:'
read(*,*) c
write(*,*)'Enter the angle THETA  in degrees'
read(*,*) theta

a=c*cos(theta)
b=c*sin(theta)

write(*,*) 'The length of adjacent side is ', a
write(*,*) 'The length of opposite side is ', b

end program triangle
