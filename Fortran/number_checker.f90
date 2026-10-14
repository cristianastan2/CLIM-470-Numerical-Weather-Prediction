program number_checker
      implicit none

      real:: number

! Prompt user for input
      print*, "Enter a number:"
      read*, number

! --- Part a: positive, negative or zero ---




! --- Part b: odd or even ---
    if (mod(number, 2.0) == 0.0) then
        print*, number, "is even."
    else
        print*, number, "is odd."
    end if


! --- Part c: divisible by 3 and/or 5 ---




! --- Part d: is the square root a positive integer? ---


end program number_checker

