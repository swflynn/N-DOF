PROGRAM s_mat_eval
 IMPLICIT NONE
  
 INTEGER :: d, Nsobol
 INTEGER, PARAMETER :: deg = 10
 INTEGER :: n, p, i, j, k, m, q
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: scrambled_u
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: scrambled_z
 DOUBLE PRECISION, DIMENSION(1:10) :: herm, coef 
 DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A

 d=1
 Nsobol=10000
 
 ALLOCATE(scrambled_u(d, Nsobol))
 ALLOCATE(scrambled_z(d, Nsobol))
 
 OPEN(UNIT=6, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
 READ(6,*) scrambled_u
 CLOSE(UNIT=6)

 !=========================Get each sobol point normal_dist=========================!
 DO n = 1, Nsobol
    CALL scrambled_sobol_stdnormal(d, scrambled_u(:,n), scrambled_z(:,n))
 END DO
 scrambled_z = scrambled_z/SQRT(2.)
 !=========================Get each sobol point normal_dist=========================!

!=========================Get each wavefn coef.=========================!
  coef(1) = 1.0
  coef(2) = 1.0 / (SQRT(2.0))
  DO p = 3,deg
    coef(p) = coef(p-1)*(1 /SQRT(2.0*REAL(p-1)))
  END DO
!=========================Get each wavefn coef.=========================!
  OPEN(UNIT=9, FILE='converge.dat') 
!=========================evaluate each sobol point=========================!
DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*scrambled_z(1,i)       
!=============evaluate each herm for a single sobol point===================!
        DO j = 3,deg      
             herm(j) =(2.0*scrambled_z(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
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
  IF (mod(i,100)==0) THEN
      WRITE(9,*) (A) / REAL(i)
  END IF

  END DO

  CLOSE(UNIT=9)
!=========================evaluate each sobol point=========================!
  A = A / Nsobol 
 
!=========================write out matrix elements=========================!
  OPEN(UNIT=10, FILE='final_matrix.dat')
  do q=1,deg
     Write(10,*) A(1:deg,q)
  enddo
  CLOSE(UNIT=10)
!=========================write out matrix elements=========================! 
 
 
 DEALLOCATE(scrambled_u)
 DEALLOCATE(scrambled_z)

 
END PROGRAM s_mat_eval
