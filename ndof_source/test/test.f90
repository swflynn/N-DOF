! attempting to calculate an integral using sobol numbers

! first we need to generate the sobol sequence

PROGRAM main
      USE  sobol

      IMPLICIT NONE
      INTEGER (kind=8) :: m
      INTEGER (kind=8) :: n
      INTEGER (kind=8) :: skip
      REAL (kind=8), DIMENSION(:) :: r

      m = 1
      n = 10
      skip = 0

      CALL i8_sobol_generate(m, n, skip, r)

END PROGRAM main
