! Intend for this to be the main program for running my ndof module. 
! will put all functions/subroutines within the module, then run this. 


PROGRAM main
  USE ndof_module
  IMPLICIT NONE

  REAL*8 :: A(0:10)
  REAL*8 :: B(0:10,0:10)
  INTEGER :: degree,k

  WRITE(*,25) degree

  CALL Hermite_Coeff(5,A,B)

  DO k = 0, degree
      WRITE(*,50) k, A(k)
  END DO 

  stop

	25 format('Coefficients for Hermite Polynomial',i2/)
	50 format('A(',i2,') = ',f10.0)



END PROGRAM main
