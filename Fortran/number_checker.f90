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




! --- Part d: is the square root a positive integer? ---


end program number_checker

