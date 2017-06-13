MODULE NDOF_s_module
IMPLICIT NONE
INTEGER, PARAMETER :: d = 2           ! Spatial Dimension
INTEGER, PARAMETER :: Vmax = 2        ! Maximum Excitation
INTEGER :: Jmax
INTEGER, ALLOCATABLE :: v(:,:)
INTEGER, PARAMETER :: deg = 2         ! Highest Degree Polynomial

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

implicit none
REAL :: initial_time, final_time
DOUBLE PRECISION :: B
DOUBLE PRECISION, ALLOCATABLE :: coef(:), scrambled_u(:), scrambled_z(:), herm(:,:), A(:,:,:), U(:,:) 
INTEGER :: i, j, j1, k, m
INTEGER Nsobol

CALL CPU_TIME(initial_time)

!=================================Read input file========================================!
OPEN(60,FILE='input.dat')
READ(60,*) Nsobol
CLOSE(60)
!=================================Read input file========================================!

ALLOCATE(coef(deg), scrambled_u(d), scrambled_z(d), herm(deg,d), A(deg,deg,d))

!=================================Coeficient Generation=================================!
! Start By generating the Coeficients for the Hermite Polynomials 
! for however many spatial dimensions you will be calcualtion. 
! P170 McQuarrie the coef = 1 / (2^v * v!)^.5
! Coef * Herm = HO Wavefunction
!========================================================================================!
coef(1) = 1.0
coef(2) = 1.0 / SQRT(2.0)
DO i = 3, deg
  coef(i) = coef(i-1) * (1 / SQRT(2.*(i-1)))
END DO

! Make Sure the coeficients are being generated, this can be removed
WRITE(*,*) 'Coef Test'
WRITE(*,*) coef
! Make Sure the coeficients are being generated, this can be removed
!=================================Coeficient Generation=================================!

!========================Jmax and v(d,Jmax)=============================================!
CALL permutation(d,Vmax)
! Make sure we are calculation all the permutation v(d,Jmax)
WRITE(*,*) 'Test v'
WRITE(*,*) v
WRITE(*,*) 'Test Jmax'
WRITE(*,*) Jmax
ALLOCATE(U(Jmax,Jmax))
U=0d0
!========================Jmax and v(d,Jmax)=============================================!

!=================================Sobol Points=================================!
!========================================================================================!
! The Scrambled Sobol Points are generated seperatly through a matlab code.
! We need to open the file containing our points: s_sobol_unif.dat
! For each point use beasley_springer_moro function to convert to a normal distribution
! Make herm for each point up to deg and multiply by coef to get the polynomial for each point
!========================================================================================!
OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
DO i = 1, Nsobol              
  READ(70,*) scrambled_u
  WRITE(*,*) 'Here is the point'                            ! Just testing remove this after
  WRITE(*,*) scrambled_u                                    ! Jsut testing remove this after
  scrambled_z(:)=beasley_springer_moro(scrambled_u)
  WRITE(*,*) 'Here is the transformed point'                ! Just testing remove this after
  WRITE(*,*) scrambled_z                                    ! Just testing remove this after
  scrambled_z = scrambled_z/SQRT(2.)    ! This factor is from definition of normal distribution 
  herm(1,:) = 1.0             
  herm(2,:) = 2.0*scrambled_z(:)       
  DO j = 3,deg      
    herm(j,:) =(2.*scrambled_z(:)*herm(j-1,:)) - (2.*(j-2)*herm(j-2,:))
  END DO
  WRITE(*,*) 'Herm test'                                    ! Make sure we have proper hermite evaluation
  WRITE(*,*) herm                                           ! remove once done testing
  do j=1,deg
     herm(j,:)=herm(j,:)*coef(j)
  enddo
! Make Sure the Polynomial works for each dimension/point can be deleted
  WRITE(*,*) 'Harmonic oscilator test'
  WRITE(*,*) herm
! Make Sure the Polynomial works for each dimension/point can be deleted

!=================================Evaluate Herm * Herm =================================!
! Make a Matrix A that evaluates herm(deg)*herm(deg) 
! A is a colllection of 1D calculations from before it contains all of the polynomial products
! the only difference is that we repeat this same calculation for each spatial dimension
!========================================================================================!
  DO j=1,deg
    DO k=1,deg
      A(j,k,:) = herm(j,:)*herm(k,:)
    END DO
  END DO
  WRITE(*,*) 'Test A'                                     ! check HO * HO evaluations for each dimension
  WRITE(*,*) A                                            ! remove this line after we get it working

!=================================Evaluate Matrix Elements =================================!
! Our matrix U will contains the matrix elements 
!========================================================================================!

  DO j=1,Jmax
    DO j1=j,Jmax
       B = A(v(1,j),v(1,j1),1)
      DO m=2,d
        B = B * A(v(m,j),v(m,j1),m)
      END DO
      U(j,j1) = U(j,j1) + B
    END DO
  END DO


! This end do closes off the loop over all your sobol points    
END DO
CLOSE(70)
U = U/Nsobol

OPEN(UNIT=80, FILE='matrix.dat')
DO i =1,Jmax
  WRITE(80,*) U(1:Jmax,i)
END DO
CLOSE(UNIT=80)

DEALLOCATE(coef, scrambled_u, scrambled_z, herm, A, v, U)
CALL CPU_TIME(final_time)
WRITE(*,*) 'Total Time:', final_time - initial_time

END PROGRAM NDOF_s
