program number_checker
      implicit none

      real:: number

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




! --- Part c: divisible by 3 and/or 5 ---
if (modulo(number, 3.0)== 0) then
	print*,number, " is divisible by 3."
else if (modulo(number, 5.0) == 0) then
	print*,number, " is divisible by 5."
else
	print*,number, " is not divisible by 3 or 5."
endif


! --- Part d: is the square root a positive integer? ---


end program number_checker

