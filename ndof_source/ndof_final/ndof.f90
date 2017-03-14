! Please see the test directory for documentation. 
! Generate normal sobol sequence and evaluate hermite polynomials at each point
! Code is written for 1 dimension only, would need to modify for more. 

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER (kind=8) :: m, Nsobol, skip, i, j
  Double Precision, allocatable :: r(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

   m = 1                   
   Nsobol = 100           
   skip = 2              

   OPEN(unit=10, file='data.dat')

   ALLOCATE (r(m,Nsobol)) 

   CALL i8_sobol_generate(m, Nsobol, skip, r) 

   Write(10,*), 'here are your sobol points'
   WRITE(10,*), r
   WRITE(10,*), 'Next are the hermite evaluations'

   DO i = 1, Nsobol              
    herm(1) = 1.0             
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
