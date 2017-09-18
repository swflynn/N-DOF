PROGRAM s_mat_eval
 IMPLICIT NONE
 
 INTEGER ::  Nsobol
 INTEGER :: deg                               
 INTEGER, PARAMETER :: d = 1                      ! Code for 1D case only
 INTEGER :: i, j, k, m
 REAL :: initial_time, final_time
 DOUBLE PRECISION, ALLOCATABLE :: scrambled_u(:,:), scrambled_z(:), herm(:), coef(:), A(:,:)

!==========================================================================================!
!=================================Variables to set=========================================!
!======================deg: highest polynomial, Nsobol: # iterations=======================!
!==========================================================================================!
 CALL CPU_TIME(initial_time)
 Nsobol=100
 deg = 10
 ALLOCATE(scrambled_u(d, Nsobol), scrambled_z(d), herm(deg), coef(deg), A(deg,deg))
 A = 0d0
!==========================================================================================!
!==============================Read in Scrambled Sequence==================================!
!==========================================================================================!
 OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
 READ(70,*) scrambled_u
 CLOSE(UNIT=70)
!==========================================================================================!
!=============================Coeficients each Polynomial==================================!
!==========================================================================================!
  coef(1) = 1.0
  coef(2) = 1.0 / SQRT(2.)
  DO i = 3,deg
    coef(i) = coef(i-1)*(1 / SQRT(2.*(i-1)))
  END DO
  OPEN(UNIT=75, FILE='converge.dat') 
!==========================================================================================!
!=================================Evaluate Polynomials=====================================!
!==========================================================================================!
  DO i = 1, Nsobol              
    CALL scrambled_sobol_stdnormal(d, scrambled_u(:,i), scrambled_z(d))
        scrambled_z = scrambled_z/SQRT(2.)
        herm(1) = 1.0             
        herm(2) = 2.0*scrambled_z(d)       
        DO j = 3,deg      
             herm(j) =(2.*scrambled_z(d)*herm(j-1)) - (2.*(j-2)*herm(j-2))
        END DO
        herm(:)=herm(:)*coef(:)
!==========================================================================================!
!===============================Potential Energy Matrix====================================!
!==========================================================================================!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + herm(k)*herm(m)
             END DO
        END DO
!==========================================================================================!
!=================================Convergence Analysis=====================================!
!==========================================================================================!
        IF (mod(i,100)==0) THEN
            WRITE(75,*) (A) / i
        END IF
  END DO ! sobol point loop
  CLOSE(UNIT=75)
  A = A / Nsobol 
!==========================================================================================!
!======================================Final Matrix========================================!
!==========================================================================================! 
  OPEN(UNIT=80, FILE='final_matrix.dat')
  DO i=1,deg
     Write(80,*) A(1:deg,i)
  END DO
  CLOSE(UNIT=80)
  DEALLOCATE(scrambled_u, scrambled_z, herm, coef, A)
  CALL CPU_TIME(final_time)
!==========================================================================================!
!=======================================Output File========================================!
!==========================================================================================!
OPEN(UNIT=83, FILE='output.dat')
WRITE(83,*) 'Sobol Numers = ', Nsobol
WRITE(83,*) 'Polynomial Degree= ', deg
WRITE(83,*) 'This calculation ran for (s): ', final_time - initial_time
CLOSE(UNIT=83)
WRITE(*,*) 'TOTAL TIME (s): ', final_time - initial_time
 
END PROGRAM s_mat_eval
