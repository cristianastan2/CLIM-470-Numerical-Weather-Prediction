program precision_demo
!including both parts of this slide in one code here

implicit none
integer, parameter :: a = 3**3**2
integer, parameter :: b = (3**3)**2
integer, parameter :: c = 3**(3**2)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: long_int=selected_int_kind(15)
real :: single_val
double precision :: double_val
!integer :: big_number cannot get this line to work
single_val = 1.0_sp/3.0_sp
double_val = 1.0_dp/3.0_dp
!big_number = 123456789012345_long_int cannot get this line to work

print*, 'a=', a
print*, 'b=', b
print*, 'c=', c
print*, '(A,F10.6)', 'Single Precision:', single_val
print*, '(A,F20.15)', 'Double Precision:', double_val
!print '(A,10)', 'Long Integer:', big_number cannot get this line to work

stop
end program precision_demo
