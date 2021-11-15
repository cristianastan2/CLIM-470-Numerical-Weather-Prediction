program format_example

implicit none
integer :: index = -12, junk = 4, number = -12345

       write (*,200) index, index+12, junk, number
       write(*,210) index, index+12, junk, number
       write(*,220) index, index+12, junk, number

  200  format ('', 2I5, I6, I10)
  210  format ('', 2I5.0, I6, I10.8)
  220  format ('', 2I5.3, I6, I5)

end program format_example
