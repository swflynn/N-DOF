PROGRAM s_mat_eval

 IMPLICIT NONE
 INTEGER :: d, Nsobol
 INTEGER, PARAMETER :: deg = 3
 INTEGER :: i, j, k, m
 DOUBLE PRECISION, ALLOCATABLE :: scrambled_u(:,:), scrambled_z(:,:), herm(:), coef(:), A(:,:)

 d=1
 Nsobol=100000
 ALLOCATE(scrambled_u(d, Nsobol), scrambled_z(d, Nsobol), herm(deg), coef(deg), A(deg,deg))

 herm = 0d0
 coef = 0d0
 A = 0d0
 scrambled_u = 0d0
 scrambled_z = 0d0

 OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
 READ(70,*) scrambled_u
 CLOSE(UNIT=70)

 !=========================Get each sobol point normal_dist=========================!
 DO i = 1, Nsobol
    CALL scrambled_sobol_stdnormal(d, scrambled_u(:,i), scrambled_z(:,i))
 END DO
 scrambled_z = scrambled_z/SQRT(2.)
 !=========================Get each sobol point normal_dist=========================!

!=========================Get each wavefn coef.=========================!
  coef(1) = 1.0
  coef(2) = 1.0 / SQRT(2.)
  DO i = 3,deg
    coef(i) = coef(i-1)*(1 / SQRT(2.*(i-1)))
  END DO
!=========================Get each wavefn coef.=========================!
  OPEN(UNIT=75, FILE='converge.dat') 
!=========================evaluate each sobol point=========================!
  DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*scrambled_z(1,i)       
!=============evaluate each herm for a single sobol point===================!
        DO j = 3,deg      
             herm(j) =(2.*scrambled_z(d,i)*herm(j-1)) - (2.*(j-2)*herm(j-2))
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
  IF (mod(i,1000)==0) THEN
      WRITE(75,*) (A) / i
  END IF

  END DO

!=========================evaluate each sobol point=========================!
  CLOSE(UNIT=75)

  A = A / Nsobol 
 
!=========================write out final matrix elements=========================!
  OPEN(UNIT=80, FILE='final_matrix.dat')
  DO i=1,deg
     Write(80,*) A(1:deg,i)
  END DO
!=========================write out matrix elements=========================! 
  CLOSE(UNIT=80)
  DEALLOCATE(scrambled_u, scrambled_z, herm, coef, A)
 
END PROGRAM s_mat_eval
