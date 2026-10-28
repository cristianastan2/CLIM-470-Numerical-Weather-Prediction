program main
implicit none
integer :: ierr
real:: a, b, c
open(20,file ='output.dat',status='new',  action='write', iostat=ierr)
if ( ierr==0 ) then
write (*,*) 'File successfully open. Error code =',ierr
read(10,*) a, b, c
write(20,100)a, b, c
100 format(1X, 3(F5.3,'-----'))
else
write(*,*) 'File not open. Error code=',ierr
end if
end program main

