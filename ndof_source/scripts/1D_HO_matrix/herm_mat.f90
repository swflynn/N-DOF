!Code to calculate 1D-HO matrix elements using normal sobol sequence 
!Matrix Element convergence analysis 

!=========================================To Use=============================================!
!Set what degree (deg) polynomial they would like to calculate up to
!===============  NOTE: Indexing starts at 1 therefore deg = 10 := H_9 ======================!
!Set the Nsobol number of points to use for the integration (in general set skip = Nsobol )
!===============  NOTE: Matrix Elements for Spatial Dimension (d) = 1 only===================!
! User should also set an appropriate number of convergence evaluations
!Code generates a sequence of NORMALLY distributed sobol points 'norm' using sobol.f90 and sobol_stdnormal.f90
! Recursion is used to calculate the HO Wavefunction Coefficients
! Recursion is used to calculate the Hermite Poynomials 
! The HO Wavefucntion is contructed as the product of the Coeficients and Polynomials
!=========================================To Use=============================================!

PROGRAM mat_eval
  
  USE sobol
  IMPLICIT NONE
  INTEGER :: Nsobol                           ! Number of Sobol Points
  INTEGER :: n, i, j, k, m, o, p         
  INTEGER, PARAMETER :: d = 1                 ! Code for 1D case Only
  INTEGER:: deg = 10                          !polynomial we want to calculate up to
  INTEGER*8 :: skip                           !seed, set = Nsobol
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)   !vector of normal sobol points
  DOUBLE PRECISION, DIMENSION(1:10) :: herm, coef    !hermite polynomial and coef vectors
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A     ! matrix elements
	REAL :: initial_time, final_time
!=============================Variables to be set by User==================================!
  deg = 10                                    ! Indexing starts at 1 not 0!
  Nsobol = 10000000
  skip = 10000000
!=============================Variables to be set by User==================================!

CALL CPU_TIME(initial_time)

  ALLOCATE (norm(d,Nsobol))
  A=0d0
!==============================Get each sobol pointpoint===================================!
  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
  norm=norm/SQRT(2.)      
!=============================Get each sobol pointpoint====================================!

!=============================Get each wavefn coef.========================================!
  coef(1) = 1
  coef(2) = 1.0 / (SQRT(2.0))
  DO p = 3,deg
    coef(p) = coef(p-1)*(1 /SQRT(2.0*REAL(p-1)))
  END DO
  !=========================Get each wavefn coef.=========================!
  OPEN(UNIT=9, FILE='converge.dat')
!=========================evaluate each sobol point=========================!
DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(1,i)       
!=============evaluate each herm for a single sobol point===================!
        DO j = 3,deg      
             herm(j) =(2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO
!=========evaluate each herm for a single sobol point======================!

!==============evaluate each matrix element for a single point===================!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + coef(k)*herm(k)*coef(m)*herm(m)
             END DO
        END DO
 !=============evaluate each matrix element for a single point======================!

! add in function to print results as a function of N
  IF (mod(i,1000000)==0) THEN
      WRITE(9,*) (A) / REAL(i)
  END IF

  END DO

  CLOSE(9)
!=========================evaluate each sobol point=========================!
  A = A / Nsobol
!=========================write out Final Matrix =========================!
  OPEN(UNIT=10, FILE='final_matrix.dat')
  do o=1,deg
     Write(10,*) A(1:deg,o)
  enddo
  CLOSE(10)
!=========================write out matrix elements=========================!

!============================Output File====================================!
OPEN(UNIT=83, FILE='output.dat')
WRITE(83,*) 'Sobol Numers = ', Nsobol
WRITE(83,*) 'Polynomial Degree= ', deg
WRITE(83,*) 'This calculation ran for (s): ', final_time - initial_time
CLOSE(UNIT=83)

CALL CPU_TIME(final_time)
WRITE(*,*) 'Total Time:', final_time - initial_time

END PROGRAM mat_eval

!Everything is working as of 3-20-17   -Shane 
