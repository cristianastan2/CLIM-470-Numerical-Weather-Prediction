program var_example

implicit none

<<<<<<< HEAD
integer:: i
real:: j
integer, parameter:: k=4

!i = k**2
i = 25
j = float(i)/k
k=i+j
print*, 'k=', k
print*, 'i=', i
print*, 'j=', j
=======
integer:: i, j
integer, parameter:: k=4

i = k**2

print*,'i=', i
>>>>>>> bc83b8543cbf6b0b5c1c6d534f354215aa6e8338

stop

end program var_example
<<<<<<< HEAD
=======

>>>>>>> bc83b8543cbf6b0b5c1c6d534f354215aa6e8338
