!==============================================================================!
!                               Dsygv
! DSYGV computes all the eigenvalues, and optionally, the eigenvectors
! of a real generalized symmetric-definite eigenproblem, of the form
! A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
! Here A and B are assumed to be symmetric and B is also
! positive definite. 
!==============================================================================!
! A would be the hamiltonian, and B the overlap matrix for QM problems
!==============================================================================!
program mats2
implicit none
integer :: dimen, i
double precision, allocatable :: A(:,:), B(:,:), eigenvalues(:)
!==============================================================================!
!                         Variables to run dsygv                               !
!==============================================================================!
integer :: itype, info, lwork
double precision, allocatable :: work(:)
!==============================================================================!
itype = 1
!==============================================================================!
dimen = 3
allocate(A(dimen,dimen), B(dimen,dimen), eigenvalues(dimen))
!==============================================================================!
! Allocations are suggested by llapack developers and should not be changed
lwork = max(1,3*dimen-1)
allocate(work(max(1,lwork)))
work=0d0
!==============================================================================!
B=0d0
!do i=1,dimen
!B(i,i) = 1
!end do

B(1,1) = 1
B(1,2) = .5
B(1,3) = 0
B(2,1) = .5
B(2,2) = 1
B(2,3) = .5
B(3,1) = 0
B(3,2) = .5
B(3,3) = 1

A(1,1) = -1
A(1,2) = 1
A(1,3) = 1
A(2,1) = 1
A(2,2) = 7
A(2,3) = 1
A(3,1) = 1
A(3,2) = 1
A(3,3) = .5
!==============================================================================!
!                                 Dsygv
!==============================================================================!
! arguments: DSYGV(Itype, JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
! 0(Itype): (input) Specified generalized eigenvalue problem to solve
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

CALL dsygv(itype, 'n', 'u', dimen, A, dimen, B, dimen, eigenvalues, work, Lwork, info)
! if info neq there was an error, go check, if larger than 1 you are not positive definite
! You will get garbage if B is not POSITIVE DEFINITE
write(*,*) 'info ==>', info
write(*,*) 'eigenvalues'
write(*,*) eigenvalues

end program mats2
