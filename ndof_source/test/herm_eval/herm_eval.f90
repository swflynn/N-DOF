PROGRAM herm_poly_eval
      IMPLICIT NONE

      REAL, PARAMETER :: chi  = 5.0
      INTEGER :: i
      DOUBLE PRECISION, DIMENSION(0:9) :: herm

      herm(0) = 1
      herm(1) = 2.0*chi
      
       DO i = 2, 9

        herm(i) = 2.0*chi*herm(i-1) - 2.0*(i-1)*herm(i-2)

      END DO

      PRINT *, herm

END PROGRAM herm_poly_eval
