! develop code to evaluate any hermite polynomial for a value of x. 

PROGRAM herm_fn_1
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 2/28/17

  DO i = 0, n
    u = 3.0*i/n
    WRITE (*,*) u, herm_fn_1_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_1_eval(x) RESULT (value_herm_fn_1)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_1

    value_herm_fn_1 = 2.0*x

   END FUNCTION herm_fn_1_eval

END PROGRAM herm_fn_1
