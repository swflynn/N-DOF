!==============================================================================!
!                               DSYEV
! DSYEV computes all eigenvalues and, optionally, eigenvectors of a
! real symmetric matrix A.
!==============================================================================!
program mats
implicit none
integer :: dimen
double precision, allocatable :: A(:,:), eigenvalues(:)
!==============================================================================!
!                         Variables to run dsyev                               !
!==============================================================================!
integer :: info, lwork 
double precision, allocatable :: work(:)
!==============================================================================!
!==============================================================================!

dimen = 3
allocate(A(dimen,dimen), eigenvalues(dimen))
!==============================================================================!
! Allocations are suggested by llapack developers and should not be changed
lwork = max(1,3*dimen-1)
allocate(work(max(1,lwork)))
work=0d0
info = 0
!==============================================================================!

A(1,1) = 1
A(1,2) = 2
A(1,3) = 3
A(2,1) = 2
A(2,2) = 5
A(2,3) = 3
A(3,1) = 3
A(3,2) = 3
A(3,3) = 4

write(*,*) A

!==============================================================================!
!                                 Dsyev
!           Eigenvalues and Eigenvectors of Real Symmetic matrix A
!==============================================================================!
! arguments: DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
! 1(Jobz): (input) N/V eigenvalues,eigenvetors+eigenvalues
! 2(Uplo): (Input) U/L upper or lower triangle of matrix A is stored
! 3(N):    (Input)The order of the Matrix
! 4(A): (Input/Output) array A(LDA,N) A should be symmetric
! 5(LDA): (Input) Leading Dimension of the array 
! 6(W): (Output) Array(N) eigenvalues in ascending order
! 7(WORK): (Output) Returns optimal Lwork
! 8(LWORK): (Input) Efficiency Array
! 9(INFO): (Output) Success Flag
!==============================================================================!
CALL dsyev('n', 'u', dimen, A, dimen, eigenvalues, work, Lwork, info)
write(*,*) 'eigenvalues'
write(*,*) eigenvalues

end program mats
