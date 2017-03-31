! Evaluate H_4 for a given x

PROGRAM herm_fn_4
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 4/4/17

  DO i = 0, n
    u = 4.0*i/n
    WRITE (*,*) u, herm_fn_4_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_4_eval(x) RESULT (value_herm_fn_4)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_4

    value_herm_fn_4 = 16.0*(x**4) - 48.0*(x**2) + 12.0

   END FUNCTION herm_fn_4_eval

END PROGRAM herm_fn_4
