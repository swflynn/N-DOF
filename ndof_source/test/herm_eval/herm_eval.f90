PROGRAM herm_poly_eval
      IMPLICIT NONE

      REAL, DIMENSION(0:2) :: chi 
      INTEGER :: i, j, size_chi, k
      DOUBLE PRECISION, DIMENSION(0:9, 0:9) :: herm

      ! make an array with values (this would be our sobol numbers)
      chi(0) = 0.5
      chi(1) = 1.0
      chi(2) = 5.0


      size_chi = size(chi)
      size_chi = size_chi - 1         !subtract 1 because index is 0 for polynomials

      DO i = 0, size_chi
          herm(i,0) = 1.0             !initialize first 2 polynomials for recursion
          herm(i,1) = 2.0*chi(i)

          PRINT *, 'Evaluating the polynomials at: ', chi(i)

          PRINT *, herm(i,0)
          PRINT *, herm(i,1)

         DO j = 2, 9

             herm(i,j) = 2.0*chi(i)*herm(i, j-1) - 2.0*(j-1)*herm(i, j-2)
             PRINT *, herm(i,j)

         END DO

      END DO

END PROGRAM herm_poly_eval
