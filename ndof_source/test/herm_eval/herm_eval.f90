PROGRAM herm_poly_eval
      IMPLICIT NONE

      REAL, DIMENSION(0:1) :: chi 
      INTEGER :: i, j
      DOUBLE PRECISION, DIMENSION(0:9, 0:9) :: herm

      ! make an array with values (this would be our sobol numbers)
      chi(0) = 3.0
      chi(1) = 5.0

      ! initialize the first 2 polynomials for recursion 
      !I should go back and do this for a general case
      !For now the specific definitions work
      herm(0,0) = 1.0
      herm(1,0) = 1.0

      herm(0,1) = 2.0*chi(0)
      herm(1,1) = 2.0*chi(1)
      
     
      DO i = 0, 1
      PRINT *, herm(i,0)
      PRINT *, herm(i,1)

         DO j = 2, 9

             herm(i,j) = 2.0*chi(i)*herm(i, j-1) - 2.0*(j-1)*herm(i, j-2)
             PRINT *, herm(i,j)

         END DO

      END DO
    

END PROGRAM herm_poly_eval
