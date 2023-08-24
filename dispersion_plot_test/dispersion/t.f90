program test


implicit none

integer n, ifreq
character*2 ifnum2

n=1
ifreq=2

!open(10,file='/home/william019/dispersion/2Dstable/test.dat') 
!close(10)

write(ifnum2, '(i2.2)') ifreq
open(11,file='test'//ifnum2//'.dat')
close(11)


end program
