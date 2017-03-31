! Evaluate H_3 for a given x

PROGRAM herm_fn_3
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 3/3/17

  DO i = 0, n
    u = 3.0*i/n
    WRITE (*,*) u, herm_fn_3_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_3_eval(x) RESULT (value_herm_fn_3)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_3

    value_herm_fn_3 = 8.0*(x**3) - 12.0*x

   END FUNCTION herm_fn_3_eval

END PROGRAM herm_fn_3
