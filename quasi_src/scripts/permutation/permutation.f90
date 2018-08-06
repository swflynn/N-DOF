!=============================================================================80
!                               Permutations 
!=============================================================================80
!           Discussion:
!   Fortran 90 implementation for computing total number of permutations 
!   associated with a d-dimensional wavefunction, and a total maximum 
!   excitation allowed in the system starting from the ground state.
!==============================================================================!
!           Note:
!   The code assumes a total dimensionality of less than 10. 
!==============================================================================!
!           Modified:
!   18 January 2018
!           Author:
!   Shane Flynn
!==============================================================================!
program perms
implicit none
!==============================================================================!
!           Discussion: 
!   dimen  ==> The spatial dimension of the wavefunction 
!   Vmax_T ==> The largest total excitation allowed in the system.
!   Jmax   ==> The total number of permutations given dimen, Vmax_T.
!==============================================================================!
integer:: dimen,Vmax,Vmax_T
integer :: i,j,k,v1,v2,v3,v4,v5,v6,v7,v8,v9,Jmax
integer, allocatable :: Vm(:),vv(:)
!==============================================================================!
!                             Read Input File
!==============================================================================!
read(*,*) dimen 
read(*,*) Vmax_T
write(*,*) 'Spatial Dimension ==> ', dimen
write(*,*) 'Maximum Total Excitation ==> ', Vmax_T
if(dimen>9) stop 'Spatial Dimension > 9 See Source Code'
allocate(Vm(dimen), vv(dimen))
write(*,*) 'Input Confirmed, Computing Permutations'
!==============================================================================!
!                       Loop Over Total Excitations
!==============================================================================!
do i=0,Vmax_T
    Vmax = i
    Jmax = 0
    j=0
    Vm(1)=Vmax
    do v1=0,Vm(1)
        vv(1)=v1
        if(dimen>= 2) then
            Vm(2)=Vm(1)-vv(1)
            do v2=0,Vm(2)
                vv(2)=v2
                if(dimen>= 3) then
                    Vm(3)=Vm(2)-vv(2)
                    do v3=0,Vm(3)
                        vv(3)=v3
                        if(dimen>= 4) then
                            Vm(4)=Vm(3)-vv(3)
                            do v4=0,Vm(4)
                                vv(4)=v4
                                if(dimen>= 5) then
                                    Vm(5)=Vm(4)-vv(4)
                                    do v5=0,Vm(5)
                                        vv(5)=v5
                                        if(dimen>= 6) then
                                            Vm(6)=Vm(5)-vv(5)
                                            do v6=0,Vm(6)
                                                vv(6)=v6
                                                if(dimen>= 7) then
                                                    Vm(7)=Vm(6)-vv(6)
                                                    do v7=0,Vm(7)
                                                        vv(7)=v7
                                                        if(dimen>= 8) then
                                                            Vm(8)=Vm(7)-vv(7)
                                                            do v8=0,Vm(8)
                                                                vv(8)=v8
                                                                if(dimen== 9) then
                                                                    Vm(9)=Vm(8)-vv(8)
                                                                    do v9=0,Vm(9)
                                                                        vv(9)=v9
                                                                        j=j+1
                                                                    enddo
                                                                else
                                                                    j=j+1
                                                                endif
                                                            enddo
                                                        else
                                                            j=j+1
                                                        endif
                                                    enddo
                                                else
                                                    j=j+1
                                                endif
                                            enddo
                                        else
                                            j=j+1
                                        endif
                                    enddo
                                else 
                                    j=j+1 
                                endif
                            enddo
                        else
                            j=j+1
                        endif
                    enddo
                else
                    j=j+1
                endif
            enddo
        else
            j=j+1
        endif
    enddo 
    Jmax = j
!    write(*,*) 'Allowed Total Excitation ==> ', Vmax
!    write(*,*) 'Permutations for Given Excitation ==> ', Jmax
end do
write(*,*) 'Jmax ==>', Jmax

End Program perms
