! Evaluate H_2 for a given x.

PROGRAM herm_fn_2
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 3/2/17

  DO i = 0, n
    u = 3.0*i/n
    WRITE (*,*) u, herm_fn_2_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_2_eval(x) RESULT (value_herm_fn_2)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_2

    value_herm_fn_2 = 4.0*(x**2) - 2.0

   END FUNCTION herm_fn_2_eval

END PROGRAM herm_fn_2
