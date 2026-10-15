program number_checker
      implicit none

      real:: number
      real:: square, sqroot    ! for gary's code	

! Prompt user for input
      print*, "Enter a number:"
      read*, number

! --- Part a: positive, negative or zero ---
!lorenzos part
if (number >0.0) then
	print*, "positive"
else if (number <0.0) then
	print*, "negative"
else 
	print*, "zero"
end if

! --- Part b: odd or even ---
    if (mod(number, 2.0) == 0.0) then
        print*, number, "is even."
    else
        print*, number, "is odd."
    end if


! --- Part c: divisible by 3 and/or 5 ---
if (modulo(number, 3.0)== 0) then
	print*,number, " is divisible by 3."
else if (modulo(number, 5.0) == 0) then
	print*,number, " is divisible by 5."
else
	print*,number, " is not divisible by 3 or 5."
endif


! --- Part d: is the square root a positive integer? ---
! i decided to analyze both the square of the entered number and it
! square root
	square = number*number
	sqroot = sqrt(number)
	print*, "number squared = ", square
	print*, "square root = ", sqroot
	if (int(square) == square) then
		print*, "number squared is integer ", number, square
	else
		print*, "number squared is NOT integer ", number, square
	end if

	if (int(sqroot) == sqroot) then
		print*, "square root of number is integer ", sqroot, number
	else
		print*, "square root of number is NOT integer", sqroot, number
	end if


end program number_checker

