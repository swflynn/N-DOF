! Quick code to call sobol.f90 and generate sobol points
! Once we generate the points we integrate using MC

PROGRAM unif_sobol_points
      USE  sobol
      IMPLICIT NONE

      INTEGER (kind=8) :: m, n, skip
      Double Precision, allocatable :: r(:,:), q(:,:)
      REAL :: end_val, start_val, int_range, integral
      INTEGER :: i
      REAL*8 :: f

      m = 1                   !spatial dimension
      n = 10000               !number of points to generate
      skip = 0                !starting sobol point
      end_val = 10.0          !Integral bounds
      start_val = 0.0
      int_range = end_val - start_val
      ALLOCATE (r(m,n))       !allocate space for arrays
      ALLOCATE (q(m,n))

      CALL i8_sobol_generate(m, n, skip, r)
      q = r*int_range    ! we now have random points sampling whole integral bound

     
      f = 0.0D0
      integral = 0.0
      DO i=1,n
          f = f+integrand(q(m,i))
      END DO
      f = f/n

      integral = (end_val-start_val)*f

      PRINT *, integral
        


CONTAINS

      FUNCTION integrand(r) RESULT (values)
        IMPLICIT NONE
        DOUBLE PRECISION :: r
        REAL :: values
        values = (r**2)*EXP(-r)
      END FUNCTION integrand

END PROGRAM unif_sobol_points
