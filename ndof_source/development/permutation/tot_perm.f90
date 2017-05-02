Program perms

IMPLICIT NONE
INTEGER, PARAMETER :: spatial_dim = 9           !======spatial dimension !!!!code only works spatial_dim = 9======!
INTEGER, PARAMETER :: vmax = 3                 !Highest Excitation to consider
integer :: j,Vm(Spatial_dim),vv(Spatial_dim),k,v1,v2,v3,v4,v5,v6,v7,v8,v9
 
j=0
 Vm(1)=Vmax
 if(Spatial_dim>9) stop 'Spatial_dim>9'
  do v1=0,Vm(1)
     vv(1)=v1
     if(Spatial_dim >= 2) then
        Vm(2)=Vm(1)-vv(1)
        do v2=0,Vm(2)
           vv(2)=v2
           if(Spatial_dim >= 3) then
              Vm(3)=Vm(2)-vv(2)
              do v3=0,Vm(3)
                 vv(3)=v3
                 if(Spatial_dim >= 4) then
                    Vm(4)=Vm(3)-vv(3)
                    do v4=0,Vm(4)
                       vv(4)=v4
                       if(Spatial_dim >= 5) then
                          Vm(5)=Vm(4)-vv(4)
                          do v5=0,Vm(5)
                             vv(5)=v5
                             if(Spatial_dim >= 6) then
                                Vm(6)=Vm(5)-vv(5)
                                do v6=0,Vm(6)
                                   vv(6)=v6
                                   if(Spatial_dim >= 7) then
                                      Vm(7)=Vm(6)-vv(6)
                                      do v7=0,Vm(7)
                                         vv(7)=v7
                                         if(Spatial_dim >= 8) then
                                            Vm(8)=Vm(7)-vv(7)
                                            do v8=0,Vm(8)
                                               vv(8)=v8
                                               if(Spatial_dim == 9) then
                                                  Vm(9)=Vm(8)-vv(8)
                                                  do v9=0,Vm(9)
                                                     vv(9)=v9
                                                     j=j+1
                                                  write(*,*) 'j=',j,vv !here
                                                     do k=1,Spatial_dim
!                                                        v(k,j)=vv(k) 
                                                     enddo
                                                  enddo
                                               else
                                                  j=j+1
                                                  do k=1,Spatial_dim
!                                                    v(k,j)=vv(k) 
                                                  enddo
                                               endif
                                            enddo
                                         else
                                            j=j+1
                                            do k=1,Spatial_dim
!                                               v(k,j)=vv(k) 
                                            enddo
                                         endif
                                      enddo
                                   else
                                      j=j+1
                                      do k=1,Spatial_dim
!                                         v(k,j)=vv(k) 
                                      enddo
                                   endif
                                enddo
                             else
                                j=j+1
                                do k=1,Spatial_dim
!                                   v(k,j)=vv(k)  
                                enddo
                             endif
                          enddo
                       else 
                          j=j+1 
                          do k=1,Spatial_dim
!                             v(k,j)=vv(k)  
                          enddo
                       endif
                    enddo
                 else
                    j=j+1
                    do k=1,Spatial_dim
!                       v(k,j)=vv(k)   
                    enddo
                 endif
              enddo
           else
              j=j+1
              do k=1,Spatial_dim
!                 v(k,j)=vv(k)   
              enddo
           endif
        enddo
     else
        j=j+1
        do k=1,Spatial_dim
!           v(k,j)=vv(k)   
        enddo
     endif
  enddo 

 write(*,*) 'Jmax=', j
   
end Program perms
