! develop code to evaluate any hermite polynomial for a value of x. 

PROGRAM herm_fn_6
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 3/2/17

  DO i = 0, n
    u = 6.0*i/n
    WRITE (*,*) u, herm_fn_6_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_6_eval(x) RESULT (value_herm_fn_6)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_6

    value_herm_fn_6 = 64.0*(x**6) - 480.0*(x**4) + 720.0*(x**2) - 120.0

   END FUNCTION herm_fn_6_eval

END PROGRAM herm_fn_6
