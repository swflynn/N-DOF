! Fortran-90 program to determine the first 10 Hermite Polynomials Coefficients

!The code asks what degree polynomial you are interested in, and provides the coeficients for each polynomial up to the specified degree. 

PROGRAM Hermite
  IMPLICIT NONE

  REAL*8 :: A(0:10)      !Final coeficient for each degree
  REAL*8 :: B(0:10,0:10) 
  INTEGER :: degree,k    
  CHARACTER :: x

  PRINT *, "What degree polynomial would you like to consider" 
  READ (*,*) degree           !which Hermite polynomial you want
     WRITE(*,25) degree !write what degree polynomial you are calculating
     CALL Hermite_Coeff(degree,A,B) 
     DO k = 0, degree
       WRITE(*,50) k, A(k)  ! write out all of our polynomial coefficients
     END DO 

!Formatting for writing coefficients to screen 
25 format('Coefficients for Hermite Polynomial',i2/)
50 format('A(',i2,') = ',f10.0)

END PROGRAM HERMITE

SUBROUTINE Hermite_Coeff(degree,A,B)   !Calculate our coefficients
  IMPLICIT NONE
  INTEGER i,j,degree
  REAL*8  A(0:10), B(0:10,0:10)  
  B(0,0)=1.D0 ; B(1,0)=0.D0 ; B(1,1)=2.D0 ! Recursion requires first terms

 IF (degree == 0) THEN  ! Zeroth polynomial H_0 = 1
   A(0) = 1
   RETURN 

 ELSE IF (degree == 1) THEN ! First polynomial H_1 = 2x + 0
   A(0) = 0
   A(1) = 2
   RETURN 

  ELSE IF (degree>1) THEN  ! Calculate others recursively  starting at H_2
    DO i = 2, degree       ! loop from 2-degree
      B(i,0)=-2.D0*(i-1)*B(i-2,0) !Determines x**0 term (constant term)
      DO j = 1, i    !Determines the value for each power of x
        B(i,j)=2.D0*B(i-1,j-1)-2.D0*(i-1)*B(i-2,j)
      END dO 
    END DO 

    DO i = 0, degree ! Set the coefficients for each degree
      A(i)=B(degree,i)
    END DO
  END IF
  RETURN 

END SUBROUTINE Hermite_Coeff

! Code is working as of 2-17-17. -Shane
