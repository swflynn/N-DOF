! develop code to evaluate any hermite polynomial for a value of x. 

PROGRAM herm_fn_5
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 3/2/17

  DO i = 0, n
    u = 5.0*i/n
    WRITE (*,*) u, herm_fn_5_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_5_eval(x) RESULT (value_herm_fn_5)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_5

    value_herm_fn_5 = 32.0*(x**5) - 160.0*(x**3) + 120.0*(x)

   END FUNCTION herm_fn_5_eval

END PROGRAM herm_fn_5
