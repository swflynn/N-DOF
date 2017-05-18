! 5-17-17 Module For N-Dimensional Scrambled Sobol Calculation
! Subroutine to determine Jmax, v(d,Jmax) works

MODULE NDOF_s_module
IMPLICIT NONE
INTEGER, PARAMETER :: d=3           ! Spatial Dimension
INTEGER, PARAMETER :: Vmax=9        ! Maximum Excitation
INTEGER :: Jmax
INTEGER, ALLOCATABLE :: v(:,:)

CONTAINS


!=================================Permutation Subroutine=================================!
!========================================================================================!
! Subroutine Uses the spatial dimesnion d (global parameter) and the maximum excitation
! Vmax (global Parameter) to first determine Jmax, the number of possible combinations. 
! Next the code writes out an array v(d,Jmax) for our matrix calculation. 
!========================================================================================!
!========================================================================================!

SUBROUTINE permutation(d,Vmax)

Implicit none
INTEGER :: d, Vmax, Jmax
INTEGER :: j,Vm(d),vv(d),k,v1,v2,v3,v4,v5,v6,v7,v8,v9

Jmax = 0
j=0
       Vm(1)=Vmax
       if(d>9) stop 'Spatial_dim>9'
        do v1=0,Vm(1)
           vv(1)=v1
           if(d>= 2) then
              Vm(2)=Vm(1)-vv(1)
              do v2=0,Vm(2)
                 vv(2)=v2
                 if(d>= 3) then
                    Vm(3)=Vm(2)-vv(2)
                    do v3=0,Vm(3)
                       vv(3)=v3
                       if(d>= 4) then
                          Vm(4)=Vm(3)-vv(3)
                          do v4=0,Vm(4)
                             vv(4)=v4
                             if(d>= 5) then
                                Vm(5)=Vm(4)-vv(4)
                                do v5=0,Vm(5)
                                   vv(5)=v5
                                   if(d>= 6) then
                                      Vm(6)=Vm(5)-vv(5)
                                      do v6=0,Vm(6)
                                         vv(6)=v6
                                         if(d>= 7) then
                                            Vm(7)=Vm(6)-vv(6)
                                            do v7=0,Vm(7)
                                               vv(7)=v7
                                               if(d>= 8) then
                                                  Vm(8)=Vm(7)-vv(7)
                                                  do v8=0,Vm(8)
                                                     vv(8)=v8
                                                     if(d== 9) then
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
      WRITE(*,*) 'Spatial Dimension = ', d 
      WRITE(*,*) 'max excitation = ', Vmax
      WRITE(*,*) 'Jmax = ', Jmax
!==================With Jmax, Run again to determine v(d,Jmax)===========================!
ALLOCATE(v(d,Jmax))

j=0
 Vm(1)=Vmax
 if(d>9) stop 'Spatial_dim>9'
  do v1=0,Vm(1)
     vv(1)=v1
     if(d>= 2) then
        Vm(2)=Vm(1)-vv(1)
        do v2=0,Vm(2)
           vv(2)=v2
           if(d>= 3) then
              Vm(3)=Vm(2)-vv(2)
              do v3=0,Vm(3)
                 vv(3)=v3
                 if(d>= 4) then
                    Vm(4)=Vm(3)-vv(3)
                    do v4=0,Vm(4)
                       vv(4)=v4
                       if(d>= 5) then
                          Vm(5)=Vm(4)-vv(4)
                          do v5=0,Vm(5)
                             vv(5)=v5
                             if(d>= 6) then
                                Vm(6)=Vm(5)-vv(5)
                                do v6=0,Vm(6)
                                   vv(6)=v6
                                   if(d>= 7) then
                                      Vm(7)=Vm(6)-vv(6)
                                      do v7=0,Vm(7)
                                         vv(7)=v7
                                         if(d>= 8) then
                                            Vm(8)=Vm(7)-vv(7)
                                            do v8=0,Vm(8)
                                               vv(8)=v8
                                               if(d== 9) then
                                                  Vm(9)=Vm(8)-vv(8)
                                                  do v9=0,Vm(9)
                                                     vv(9)=v9
                                                     j=j+1
                                                     do k=1,d
                                                        v(k,j)=vv(k)       
                                                     enddo
                                                  enddo
                                               else
                                                  j=j+1
                                                  do k=1,d
                                                    v(k,j)=vv(k)           
                                                  enddo
                                               endif
                                            enddo
                                         else
                                            j=j+1
                                            do k=1,d
                                               v(k,j)=vv(k)                
                                            enddo
                                         endif
                                      enddo
                                   else
                                      j=j+1
                                      do k=1,d
                                         v(k,j)=vv(k)                      
                                      enddo
                                   endif
                                enddo
                             else
                                j=j+1
                                do k=1,d
                                   v(k,j)=vv(k)                            
                                enddo
                             endif
                          enddo
                       else 
                          j=j+1 
                          do k=1,d
                             v(k,j)=vv(k)                                  
                          enddo
                       endif
                    enddo
                 else
                    j=j+1
                    do k=1,d
                       v(k,j)=vv(k)                                       
                    enddo
                 endif
              enddo
           else
              j=j+1
              do k=1,d
                 v(k,j)=vv(k)                                              
              enddo
           endif
        enddo
     else
        j=j+1
        do k=1,d
           v(k,j)=vv(k)                                                    
        enddo
     endif
  enddo 

!write(*,*) v  

DEALLOCATE(v)

END SUBROUTINE permutation







END MODULE NDOF_s_module


PROGRAM NDOF_s

USE NDOF_s_module

implicit none

call permutation(d,Vmax)

END PROGRAM NDOF_s
