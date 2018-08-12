MODULE perm2_mod
IMPLICIT NONE

INTEGER :: Dim, Vtot, Vmax(9), Jmax
INTEGER, ALLOCATABLE :: v(:,:)

CONTAINS


SUBROUTINE permutation(d,Vtot,Vmax)

    Implicit none
    INTEGER :: d, Vtot,Vmax(d)
    INTEGER :: j,vv(d),k,v1,v2,v3,v4,v5,v6,v7,v8,v9

    j=0
    if(d>9) stop 'Spatial_dim>9'
    do v1=0,Vmax(1)
       vv(1)=v1
       if(d>= 2) then
          do v2=0,min(Vtot-vv(1),Vmax(2))
             vv(2)=v2
             if(d>= 3) then
                do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
                   vv(3)=v3
                   if(d>= 4) then
                      do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                         vv(4)=v4
                         if(d>= 5) then
                            do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                               vv(5)=v5
                               if(d>= 6) then
                                  do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                                     vv(6)=v6
                                     if(d>= 7) then
                                        do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                           vv(7)=v7
                                           if(d>= 8) then
                                              do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                                 vv(8)=v8
                                                 if(d== 9) then
                                                    do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
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
    !================With Jmax, Run again to determine v(d,Jmax)===========================!
    !WRITE(*,*) 'Jmax = ', Jmax
    ALLOCATE(v(d,Jmax))
    
    j=0
    if(d>9) stop 'Spatial_dim>9'
    do v1=0,Vmax(1)
       vv(1)=v1
       if(d>= 2) then
          do v2=0,min(Vtot-vv(1),Vmax(2))
             vv(2)=v2
             if(d>= 3) then
                do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
                   vv(3)=v3
                   if(d>= 4) then
                      do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                         vv(4)=v4
                         if(d>= 5) then
                            do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                               vv(5)=v5
                               if(d>= 6) then
                                  do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                                     vv(6)=v6
                                     if(d>= 7) then
                                        do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                           vv(7)=v7
                                           if(d>= 8) then
                                              do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                                 vv(8)=v8
                                                 if(d== 9) then
                                                    do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
                                                       vv(9)=v9
                                                       j=j+1
                                                       v(:,j)=vv 
                                                    enddo
                                                 else
                                                    j=j+1
                                                    v(:,j)=vv           
                                                 endif
                                              enddo
                                           else
                                              j=j+1
                                              v(:,j)=vv           
                                           endif
                                        enddo
                                     else
                                        j=j+1
                                        v(:,j)=vv           
                                     endif
                                  enddo
                               else
                                  j=j+1
                                  v(:,j)=vv                            
                               endif
                            enddo
                         else 
                            j=j+1 
                            v(:,j)=vv
                         endif
                      enddo
                   else
                      j=j+1
                      v(:,j)=vv           
                   endif
                enddo
             else
                j=j+1
                v(:,j)=vv           
             endif
          enddo
       else
          j=j+1
          v(:,j)=vv           
       endif
    enddo
    !WRITE(*,*) v! (This contains all of the permutations given Vtot, dim). 
  END SUBROUTINE permutation

END MODULE perm2_mod

PROGRAM main
USE perm2_mod
IMPLICIT NONE

DIM = 9
Vmax = [0,0,0,0,0,0,3,3,3]
Vtot = 4

CALL permutation(dim,Vtot,Vmax)
write(*,*) v

END PROGRAM main
