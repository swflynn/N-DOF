! Generate normal sobol sequence and evaluate hermite polynomials at each point
! Code is written for 1 dimension only, would need to modify for more. 
! This is a modified version of unif_sob_herm

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER (kind=8) :: m, Nsobol, skip, i, j
  Double Precision, allocatable :: r(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

   m = 1                   !spatial dimension only works for 1 currently
   Nsobol = 100               !number of points to generate
   skip = 2                !Seed starting point for sobol sequence

   OPEN(unit=10, file='data.dat')


   ALLOCATE (r(m,Nsobol))       !allocate space for sobol points array

   CALL i8_sobol_generate(m, Nsobol, skip, r) !populates r(m,n) 

   Write(10,*), 'here are your sobol points'
   WRITE(10,*), r
   WRITE(10,*), 'Next are the hermite evaluations'

  DO i = 1, Nsobol              
    herm(1) = 1.0             !initialize first 2 polynomials for recursion
    herm(2) = 2.0*r(m,i)       

    WRITE(10,*), 'next point', r(m,i)

    WRITE(10,*), herm(1)
    WRITE(10,*), herm(2)

    DO j = 3,10       

      herm(j) = (2.0*r(m,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
      WRITE(10,*), herm(j)

    END DO

  END DO
CLOSE(10)
END PROGRAM main
