! Code to calculate integral of a function using uniform distribution of sobol points

! Call the sobol.f90 module to generate a sequence of sobol points
!Then do MC integration numerically over a function defined in the program. 

!======================================TO USE================================================! 
!Specify spatial dimension: m, Number of Sobol Points to generate: n, Skip when generating: seed (set the same as n)
! Output prints integral of interest numerical solution
!======================================TO USE================================================! 

PROGRAM unif_sobol_points
      USE  sobol
      IMPLICIT NONE

      INTEGER (kind=8) :: m, n, skip
      Double Precision, allocatable :: r(:,:), q(:,:)
      REAL :: end_val, start_val, int_range, integral
      INTEGER :: i
      REAL*8 :: f


      m = 1                                         !spatial dimension
      n = 100                                       !number of sobol points generated
      skip = 100                                    !starting sobol point
      end_val = 10.0                                !Integral bounds
      start_val = 0.0
      int_range = end_val - start_val               !range of integral to normalize MC
      f = 0.0D0                                     ! initialize values for integral
      integral = 0.0

      ALLOCATE (r(m,n))       
      ALLOCATE (q(m,n))       

!Call the module to generate a uniform sobol sequence
      CALL i8_sobol_generate(m, n, skip, r) 

! sobol numbers are [0,1] we need to sample the whole function range
      q = r*int_range    
    
! for each value within q(m,n) call function and add up values
      DO i=1,n
          f = f+integrand(q(m,i))
      END DO

! normalize your MC summation
      f = f/n
      integral = (end_val-start_val)*f

! Here is the numerical approximation of the integral!
      WRITE(*,*) 'The numerical Integral is'
      WRITE(*,*) integral
        
! Function that we are trying to integrate
CONTAINS

      FUNCTION integrand(r) RESULT (values)
        IMPLICIT NONE
        DOUBLE PRECISION :: r
        REAL :: values
        values = (r**2)*EXP(-r)
      END FUNCTION integrand

END PROGRAM unif_sobol_points

!This code is working properly 3/8/17 -Shane
