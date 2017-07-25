MODULE NDOF_module
IMPLICIT NONE
INTEGER, PARAMETER :: d = 2           ! Spatial Dimension
INTEGER, PARAMETER :: Vmax = 2        ! Maximum Excitation
INTEGER :: Jmax
INTEGER, ALLOCATABLE :: v(:,:)

CONTAINS

SUBROUTINE permutation(d,Vmax)

Implicit none
INTEGER :: d, Vmax
INTEGER :: j,Vm(d),vv(d),k,v1,v2,v3,v4,v5,v6,v7,v8,v9

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
!      WRITE(*,*) 'Spatial Dimension = ', d 
!      WRITE(*,*) 'max excitation = ', Vmax
!      WRITE(*,*) 'Jmax = ', Jmax
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

END SUBROUTINE permutation

END MODULE NDOF_module

PROGRAM NDOF
USE NDOF_module

IMPLICIT NONE
REAL :: initial_time, final_time
DOUBLE PRECISION :: B
DOUBLE PRECISION, ALLOCATABLE :: scrambled_z(:), herm(:,:), A(:,:,:), U(:,:), U1(:,:)
DOUBLE PRECISION, ALLOCATABLE :: C(:,:), FV1(:), FV2(:), eigenvalues(:)
INTEGER :: j, j1, k, m, n, IERR, Nsobol
INTEGER*8 :: i, skip

CALL CPU_TIME(initial_time)

Nsobol = 1000000
skip = 1000000 
ALLOCATE(scrambled_z(d), herm(0:Vmax,d), A(0:Vmax,0:Vmax,d))

!========================Jmax and v(d,Jmax)=============================================!
CALL permutation(d,Vmax)
! Make sure we are calculating all the permutation v(d,Jmax) 
!WRITE(*,*) 'spacial dimensions = ', d
!WRITE(*,*) 'Maximum Excitation = ', Vmax
!WRITE(*,*) 'Jmax = ', Jmax                      ! remove once done testing
!WRITE(*,*) 'Test v'                         ! remove these once done testing, v is large
!WRITE(*,*) v                                ! remove once done testing
ALLOCATE(U(Jmax,Jmax), U1(Jmax,Jmax), C(Jmax,Jmax), FV1(Jmax), FV2(Jmax))
ALLOCATE(eigenvalues(Jmax))                      !Allocate our matrix elements
U=0d0
U1=0d0
!========================Jmax and v(d,Jmax)=============================================!

!=================================Sobol Points=================================!
!========================================================================================!
! Make herm for each point up to deg we have redefined the Hermite Polynomial
! Therefore we no longer have a coeficient (our definition or,alized them)
!========================================================================================!
OPEN(UNIT=80, FILE='matrix.dat')
OPEN(UNIT=81, FILE='eigenvalues.dat')
DO i = 1, Nsobol
  CALL sobol_stdnormal(d,skip,scrambled_z)
!  WRITE(*,*) 'Here is the transformed point'                            ! Just testing remove this after
!  WRITE(*,*) scrambled_z                                    ! Just testing remove this after
  scrambled_z = scrambled_z/SQRT(2.)    ! This factor is from definition of normal distribution
  herm(0,:) = 1.0                       ! Re-normalized hermite polynomial now
  herm(1,:) = SQRT(2.)*scrambled_z(:)       
  DO j = 2,Vmax      
    herm(j,:) = (SQRT(2./j)*scrambled_z(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
  END DO
!  WRITE(*,*) 'Herm test'                           ! Make sure we have proper hermite evaluation
!  WRITE(*,*) herm                                           ! remove once done testing
!=================================Evaluate Herm * Herm =================================!
! Make a Matrix A that evaluates herm(deg)*herm(deg) 
! A is a colllection of 1D calculations from before it contains all of the polynomial products
! the only difference is that we repeat this same calculation for each spatial dimension
!========================================================================================!
  DO m=0,Vmax
    DO n=0,Vmax
      A(m,n,:) = herm(m,:)*herm(n,:)
    END DO
  END DO
!  WRITE(*,*) 'Test A'                                     ! check HO * HO evaluations for each dimension
!  WRITE(*,*) A                                            ! remove this line after we get it working

!==============================Evaluate Matrix Elements =================================!
! Our matrix U will contains the matrix elements U1 for partial average
!========================================================================================!
  DO j=1,Jmax         ! Loop over jmax
    DO j1=j,Jmax 
      B=A(v(1,j),v(1,j1),1)
      DO k=2,d
        B=B*A(v(k,j),v(k,j1),k)
      END DO
      U(j1,j) = U(j1,j) + B
    END DO
  END DO
! Write out partial average and flush
  if(mod(i,10000)==0) then
    U = U + U1
    U1 = 0d0
    C=U/i
    WRITE(80,*) i, C !Writes entire matrix out
    call RS(Jmax,Jmax,C,eigenvalues,0,B,FV1,FV2,IERR) 
    write(81,*) i, ABS(1-eigenvalues(:))
    flush(80)
    flush(81)
  endif

! This end do closes off the loop over all your sobol points    
END DO
CLOSE(UNIT=80)
CLOSE(UNIT=81)

!DEALLOCATE(scrambled_z, herm, A, v, U, U1, C)
CALL CPU_TIME(final_time)
WRITE(*,*) 'Total Time:', final_time - initial_time

OPEN(UNIT=83, FILE='output.dat')
WRITE(83,*) 'Sobol Numers = ', Nsobol
WRITE(83,*) 'spacial dimensions = ', d
WRITE(83,*) 'Maximum Excitation = ', Vmax
WRITE(83,*) 'Jmax = ', Jmax
WRITE(83,*) 'This calculation ran for (s): ', final_time - initial_time
CLOSE(UNIT=83)

END PROGRAM NDOF
