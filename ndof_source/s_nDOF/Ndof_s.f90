MODULE NDOF_s_module
IMPLICIT NONE
INTEGER, PARAMETER :: d = 2           ! Spatial Dimension
INTEGER, PARAMETER :: Vmax = 2        ! Maximum Excitation
INTEGER :: Jmax
INTEGER, ALLOCATABLE :: v(:,:)

CONTAINS

!=================================Transformation Function=================================!
!========================================================================================!
!> Computes the inverse cumulative density function (CDF), i.e., the quantile,
! of the standard normal distribution given u uniform on the unit hypercube.
! This is used to transform the uniformly distributed points to normally distributed
!========================================================================================!
!========================================================================================!

FUNCTION beasley_springer_moro(u)
    IMPLICIT NONE

    INTEGER :: i, j
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u

    DOUBLE PRECISION :: r
    DOUBLE PRECISION, DIMENSION(SIZE(u)) :: beasley_springer_moro, y


    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: a = (/ &
            2.50662823884, &
            -18.61500062529, &
            41.39119773534, &
            -25.44106049637 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: b = (/ &
            -8.47351093090, &
            23.08336743743, &
            -21.06224101826, &
            3.13082909833 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:8) :: c = (/ &
            0.3374754822726147, &
            0.9761690190917186, &
            0.1607979714918209, &
            0.0276438810333863, &
            0.0038405729373609, &
            0.0003951896511919, &
            0.0000321767881768, &
            0.0000002888167364, &
            0.0000003960315187 /)


    y = u - 0.5D0

    DO j = 1, SIZE(u)
        IF (ABS(y(j)) < 0.42) THEN
            r = y(j)*y(j)
            beasley_springer_moro(j) = y(j)*(((a(3)*r + a(2))*r + a(1))*r + a(0))/((((b(3)*r + b(2))*r + b(1))*r + b(0))*r + 1)
        ELSE
            IF (y(j) > 0) THEN
                r = LOG(-LOG(1-u(j)))
            ELSE IF (y(j) < 0) THEN
                r = LOG(-LOG(u(j)))
            END IF
            beasley_springer_moro(j) = c(0) + r*(c(1) + r*(c(2) + r*(c(3) + r*(c(4) + r*(c(5) + r*(c(6) + r*(c(7) + r*c(8))))))))
            IF (y(j) < 0) THEN
                beasley_springer_moro(j) = -beasley_springer_moro(j)
            END IF
        END IF
    END DO

END FUNCTION beasley_springer_moro

!=================================Permutation Subroutine=================================!
!========================================================================================!
! Subroutine Uses the spatial dimesnion d and the maximum excitation
! Vmax to first determine Jmax, the number of possible combinations. 
! Next the code writes out an array v(d,Jmax) for our matrix calculation, which contains 
! a list of polynomial degrees to satisfy permutations. 
!========================================================================================!
!========================================================================================!

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

END MODULE NDOF_s_module

PROGRAM NDOF_s
USE NDOF_s_module

IMPLICIT NONE
REAL :: initial_time, final_time
DOUBLE PRECISION :: B
DOUBLE PRECISION, ALLOCATABLE :: scrambled_u(:), scrambled_z(:), herm(:,:), A(:,:,:), U(:,:), C(:,:),FV1(:),FV2(:),eigenvalues(:)
INTEGER :: i, j, j1, k, m, n, IERR
INTEGER Nsobol

CALL CPU_TIME(initial_time)

!=================================Read input file========================================!
OPEN(60,FILE='input.dat')
READ(60,*) Nsobol
CLOSE(60)
!=================================Read input file========================================!

ALLOCATE(scrambled_u(d), scrambled_z(d), herm(0:Vmax,d), A(0:Vmax,0:Vmax,d))

!========================Jmax and v(d,Jmax)=============================================!
CALL permutation(d,Vmax)
! Make sure we are calculating all the permutation v(d,Jmax) 
WRITE(*,*) 'spacial dimensions = ', d
WRITE(*,*) 'Maximum Excitation = ', Vmax
WRITE(*,*) 'Jmax = ', Jmax                      ! remove once done testing
WRITE(*,*) 'Test v'                         ! remove these once done testing, v is large
WRITE(*,*) v                                ! remove once done testing
ALLOCATE(U(Jmax,Jmax),C(Jmax,Jmax),FV1(Jmax),FV2(Jmax),eigenvalues(Jmax))                      !Allocate our matrix elements
U=0d0
!========================Jmax and v(d,Jmax)=============================================!

!=================================Sobol Points=================================!
!========================================================================================!
! The Scrambled Sobol Points are generated seperatly through a matlab code.
! We need to open the file containing our points: s_sobol_unif.dat
! For each point use beasley_springer_moro function to convert to a normal distribution
! Make herm for each point up to deg we have redefined the Hermite Polynomial
! Therefore we no longer have a coeficient (our definition or,alized them)
!========================================================================================!
OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
OPEN(UNIT=80, FILE='matrix.dat')
OPEN(UNIT=81, FILE='eigenvalues.dat')
DO i = 1, Nsobol              
  READ(70,*) scrambled_u
  !WRITE(*,*) 'Here is the point'                            ! Just testing remove this after
  !WRITE(*,*) scrambled_u                                    ! Just testing remove this after
  scrambled_z(:)=beasley_springer_moro(scrambled_u)
  scrambled_z = scrambled_z/SQRT(2.)    ! This factor is from definition of normal distribution
!  WRITE(*,*) 'Here is the transformed point'                ! Just testing remove this after
!  WRITE(*,*) scrambled_z                                    ! Just testing remove this after
  herm(0,:) = 1.0                       ! Re-normalized hermite polynomial now
  herm(1,:) = SQRT(2.)*scrambled_z(:)       
  DO j = 2,Vmax      
    herm(j,:) = (SQRT(2./j)*scrambled_z(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
  END DO
!  WRITE(*,*) 'Herm test'                                    ! Make sure we have proper hermite evaluation
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
  !WRITE(*,*) 'Test A'                     ! check HO * HO evaluations for each dimension
  !WRITE(*,*) A                       ! remove this line after we get it working

!=================================Evaluate Matrix Elements =================================!
! Our matrix U will contains the matrix elements 
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

!Evaluation of convergence
  if(mod(i,1000000)==0) then
     C=U/i
!   write(80,*) '******* Nsobol=',i,'*****' not convenient for analysis
!   DO j =1,Jmax          Not sure why we have a do loop here - Shane 7-1-17
     WRITE(80,*) i, C !Writes entire matrix out
      !  WRITE(80,*) U(:,i)
!   END DO
     write(80,*)
     flush(80)
     call RS(Jmax,Jmax,C,eigenvalues,0,B,FV1,FV2,IERR) 
     write(81,*) i, ABS(1-eigenvalues(:))       ! should be 1 get 1-eigenvalue to plot 
     flush(81)
  endif

! This end do closes off the loop over all your sobol points    
END DO
CLOSE(70)
CLOSE(UNIT=80)
CLOSE(UNIT=81)

!DEALLOCATE(scrambled_u, scrambled_z, herm, A, v, U)
CALL CPU_TIME(final_time)
WRITE(*,*) 'Total Time:', final_time - initial_time

OPEN(UNIT=83, FILE='output.dat')
WRITE(83,*) 'Sobol Numers = ', Nsobol
WRITE(83,*) 'spacial dimensions = ', d
WRITE(83,*) 'Maximum Excitation = ', Vmax
WRITE(83,*) 'Jmax = ', Jmax                    
WRITE(83,*) 'This calculation ran for (s): ', final_time - initial_time
CLOSE(UNIT=83)

END PROGRAM NDOF_s
