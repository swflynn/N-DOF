! This is the program main for ndof. 
!Please see the test directory for documentation. 

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER (kind=8) :: m, n, skip, i, j
  Double Precision, allocatable :: r(:,:)
  REAL :: end_val, start_val, int_range, integral
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: herm
! Hard code number of sobol points to evaluate here. 
! Could generalize to n points, however not good for memory, as n gets massive.

!===================taken from unif_sobol===================================!

   m = 1                   !spatial dimension only works for 1 currently
   n = 2               !number of points to generate
   skip = 2                !starting sobol point

   ALLOCATE (r(m,n))       !allocate space for sobol points array

   CALL i8_sobol_generate(m, n, skip, r) !populates r(m,n) 

   PRINT *, 'here are your sobol points'
   PRINT *, r
   Print *, 'Next are the hermite evaluations'
!===================taken from unif_sobol===================================!

!===================taken from herm_eval===================================!
  DO i = 1, n              
    herm(i,1) = 1.0             !initialize first 2 polynomials for recursion
    herm(i,2) = 2.0*r(m,i)        !This is correct

    PRINT *, 'next point'

    PRINT *, herm(i,1)
    PRINT *, herm(i,2)

    DO j = 3,10       

      herm(i,j) = (2.0*r(m,i)*herm(i, j-1)) - (2.0*(j-2)*herm(i, j-2))
      PRINT *, herm(i,j)

    END DO

  END DO
!===================taken from herm_eval===================================!

END PROGRAM main
