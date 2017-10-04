PROGRAM perm
IMPLICIT NONE

INTEGER, PARAMETER :: spatial_dim=3
INTEGER, PARAMETER :: Vmax=2
INTEGER :: j,Vm(Spatial_dim),vv(Spatial_dim),k,v1,v2,v3,v4,v5,v6,v7,v8,v9, Jmax
INTEGER, ALLOCATABLE :: v(:,:)
!====================variables for knonicker delta matrix====================================!
double precision, allocatable :: kd_mat(:,:), omega(:)
double precision :: enum
integer :: i, i1, qnum

write(*,*) spatial_dim
write(*,*) Vmax

allocate(omega(spatial_dim))
do i=0,spatial_dim
  omega(i) = i
end do

Jmax = 0
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
      WRITE(*,*) 'Spatial Dimension = ', spatial_dim
      WRITE(*,*) 'max excitation = ', Vmax
      WRITE(*,*) 'Jmax = ', Jmax

!==================With Jmax, Run again to determine v(d,Jmax)===========================!
ALLOCATE(v(spatial_dim,Jmax))

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
                                                     do k=1,Spatial_dim
                                                        v(k,j)=vv(k)       
                                                     enddo
                                                  enddo
                                               else
                                                  j=j+1
                                                  do k=1,Spatial_dim
                                                    v(k,j)=vv(k)           
                                                  enddo
                                               endif
                                            enddo
                                         else
                                            j=j+1
                                            do k=1,Spatial_dim
                                               v(k,j)=vv(k)                
                                            enddo
                                         endif
                                      enddo
                                   else
                                      j=j+1
                                      do k=1,Spatial_dim
                                         v(k,j)=vv(k)                      
                                      enddo
                                   endif
                                enddo
                             else
                                j=j+1
                                do k=1,Spatial_dim
                                   v(k,j)=vv(k)                            
                                enddo
                             endif
                          enddo
                       else 
                          j=j+1 
                          do k=1,Spatial_dim
                             v(k,j)=vv(k)                                  
                          enddo
                       endif
                    enddo
                 else
                    j=j+1
                    do k=1,Spatial_dim
                       v(k,j)=vv(k)                                       
                    enddo
                 endif
              enddo
           else
              j=j+1
              do k=1,Spatial_dim
                 v(k,j)=vv(k)                                              
              enddo
           endif
        enddo
     else
        j=j+1
        do k=1,Spatial_dim
           v(k,j)=vv(k)                                                    
        enddo
     endif
  enddo 
write(*,*) 'v'
write(*,*) v  
write(*,*)'Jmax'
write(*,*) Jmax
Allocate(kd_mat(Jmax,Jmax))
kd_mat = 0d0
qnum = 0
enum = 0d0
do i=1,Jmax
  do i1=1,Jmax
    if(i.ne.i1) THEN
      kd_mat(i,i1)=0d0
    else
      qnum = 0d0
      do j=1,spatial_dim
        write(*,*) 'here is the element of v'
       write(*,*) v(j,i)
       qnum = qnum + v(j,i)
      end do
      write(*,*) 'here is qnum'
      write(*,*) qnum
      enum = 0d0
      do k=0,qnum       ! hermite index form 0, so are omega
          enum = enum + omega(k)*(k+.5)
          write(*,*) 'here is enum'
          write(*,*) enum
      end do
      write(*,*) 'final enum'
      write(*,*) enum
      kd_mat(i,i1)= enum
    end if
  end do
end do
write(*,*) 'kd_mat'
write(*,*) kd_mat

END Program perm

! Code is working as of 5-17-17 
