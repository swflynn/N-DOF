module tester_mod
    implicit none
    contains

    subroutine lam(a,b,c)
    implicit none
    integer :: a,b,c
    c = a+b
    write(*,*) 'result of lam ', c
    end subroutine lam
   
end module

program tester
use tester_mod
implicit none
integer :: i,j,k,l

i =1
j = 2
write(*,*) 'i', i
write(*,*) 'j', j
call lam(i,j,l)
k=1+ l
write(*,*) k
end program tester
