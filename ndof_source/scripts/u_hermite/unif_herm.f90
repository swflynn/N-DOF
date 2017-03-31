! Generate normal sobol sequence and evaluate hermite polynomials H-0-H-9 at each point in the sequence (1D only). 
! Generate sequence of uniformly distributed sobol points 
! Evaluates each  polynomial for a given point iteratively
! To use:
! Reguires sobol.f90 module, set Nsobol (m=1 only) and skip. 
! Output is data.dat, has each sobol point, and H_0-H_9 evaluated at that point recursively
! Code is working fine 3/11/17 -Shane

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER (kind=8) :: m, Nsobol, skip, i, j
  Double Precision, allocatable :: r(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

   m = 1                   !spatial dimension only works for 1 currently
   Nsobol = 100               !number of points to generate
   skip = 100                !Seed starting point for sobol sequence

   OPEN(unit=10, file='data.dat')


   ALLOCATE (r(m,Nsobol))       !allocate space for sobol points array

   CALL i8_sobol_generate(m, Nsobol, skip, r) !populates r(m,n) 

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
