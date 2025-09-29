program var_example

implicit none

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

stop

end program var_example
