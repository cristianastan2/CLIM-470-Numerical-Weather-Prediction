program constants_demo

implicit none
integer, parameter :: MAX_SIZE = 1000
real, parameter :: GRAVITY = 9.81
character(12), parameter :: SCHOOL_NAME = 'George Mason'
logical, parameter :: DEBUG_MODE = .false.
! MAX_SIZE = 2000
print*, 'School:', SCHOOL_NAME
print*, 'Debug Mode:', DEBUG_MODE

stop
end program constants_demo
