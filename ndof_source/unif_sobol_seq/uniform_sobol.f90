! Quick code to call sobol.f90 and generate sobol points
!This makes a uniform dist. We want normal so use the scp_lmon module, not just sobol

PROGRAM unif_sobol_points
      USE  sobol
      IMPLICIT NONE

      INTEGER (kind=8) :: m, n, skip
      Double Precision, allocatable :: r(:,:)

      m = 1                   !spatial dimension
      n = 10                  !number of points to generate
      skip = 0                !starting sobol point
      ALLOCATE (r(m,n))

      CALL i8_sobol_generate(m, n, skip, r)
      PRINT *, r              ! print out point to see if it worked

END PROGRAM unif_sobol_points
