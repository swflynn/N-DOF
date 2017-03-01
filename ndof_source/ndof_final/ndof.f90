! Intend for this to be the main program for running my ndof module. 
! will put all functions/subroutines within the module, then run this. 


PROGRAM Main
  USE ndof_mod
  IMPLICIT NONE

  REAL*8 :: A(0:10)
  REAL*8 :: B(0:10,0:10)
  INTEGER :: degree,k, herm_1

  herm_1 = 0
  print *, herm_1

  PRINT *, "What degree polynomial would you like to consider?"
  READ (*,*) degree

  CALL Hermite_Coeff(degree,A,B)

  DO k = 0, degree
      WRITE(*,50) k, A(k)
  END DO 

  stop

	50 format('A(',i2,') = ',f10.0)


END PROGRAM Main
