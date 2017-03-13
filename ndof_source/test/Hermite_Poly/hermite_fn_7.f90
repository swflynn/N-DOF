! Evaluate H_7 for a given x

PROGRAM herm_fn_7
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 3/2/17

  DO i = 0, n
    u = 2.0*i/n
    WRITE (*,*) u, herm_fn_7_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_7_eval(x) RESULT (value_herm_fn_7)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_7

    value_herm_fn_7 = 128.0*(x**7) - 13440.0*(x**5) + 3360.0*(x**3) - 1680.0*(x)

   END FUNCTION herm_fn_7_eval

END PROGRAM herm_fn_7
