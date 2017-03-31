PROGRAM s_mat_eval

 IMPLICIT NONE
 INTEGER :: d, Nsobol
 INTEGER, PARAMETER :: deg = 10
 INTEGER :: i, j, k, m
 REAL :: initial_time, final_time
 DOUBLE PRECISION, ALLOCATABLE :: scrambled_u(:,:), scrambled_z(:), herm(:), coef(:), A(:,:)

  CALL CPU_TIME(initial_time)

 d=1
 Nsobol=100000
 ALLOCATE(scrambled_u(d, Nsobol), scrambled_z(d), herm(deg), coef(deg), A(deg,deg))

 A = 0d0
 
 !=========================Read in Scrambled Sequence=========================!
 OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
 READ(70,*) scrambled_u
 CLOSE(UNIT=70)

 !=========================Convert to normal dist=========================!

!=========================Get each wavefn coef.=========================!
  coef(1) = 1.0
  coef(2) = 1.0 / SQRT(2.)
  DO i = 3,deg
    coef(i) = coef(i-1)*(1 / SQRT(2.*(i-1)))
  END DO

  OPEN(UNIT=75, FILE='converge.dat') 
!=========================evaluate each sobol point=========================!
  DO i = 1, Nsobol              
    CALL scrambled_sobol_stdnormal(d, scrambled_u(:,i), scrambled_z(d))
        scrambled_z = scrambled_z/SQRT(2.)
        herm(1) = 1.0             
        herm(2) = 2.0*scrambled_z(d)       
!=============evaluate each herm for a single sobol point===================!
        DO j = 3,deg      
             herm(j) =(2.*scrambled_z(d)*herm(j-1)) - (2.*(j-2)*herm(j-2))
        END DO
        herm(:)=herm(:)*coef(:)
!==============evaluate each matrix element for a single point===================!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + herm(k)*herm(m)
             END DO
        END DO

!==============Matrix Convergence as a function of N===================!
  IF (mod(i,100)==0) THEN
      WRITE(75,*) (A) / i
  END IF

  END DO

  CLOSE(UNIT=75)

  A = A / Nsobol 
 
!=========================write out final matrix elements=========================!
  OPEN(UNIT=80, FILE='final_matrix.dat')
  DO i=1,deg
     Write(80,*) A(1:deg,i)
  END DO

  CLOSE(UNIT=80)
  DEALLOCATE(scrambled_u, scrambled_z, herm, coef, A)

   CALL CPU_TIME(final_time)
   WRITE(*,*) 'TOTAL TIME: ', final_time - initial_time
 
END PROGRAM s_mat_eval

