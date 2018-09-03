!=============================================================================80
!                       1D Harmonic Oscillator Scrambled
!=============================================================================80
!		Discussion:
! Quasi Monte Carlo Integration, computing the Potential Energy Matrix for a 
! 1-Dimensional Harmonic Oscillator (Hermite Polynomials)
! Scrambled Sobol Sequence used for QMC.
! Basis functions are orthonormal, therefore off-diagonal=0,diagonal=1.
!==============================================================================!
! Code requires a scrambled sobol sequences, see matlab script (s_sobol.m)
!==============================================================================!
!		Modified:
! 20 March 2017
!		Author:
! Shane Flynn
!==============================================================================!
Program HO_1D_s
!==============================================================================!
!                                 Variables
!==============================================================================!
!               Integer:
! dimen         ==> Spatial Dimension (parameter=1 for this code)
! exc           ==> maximum excitation (herminte polynomial)
! Nsobol        ==> Number of Sobol Points for MC Integration
! data_freq     ==> Interval for Convergence Analysis
!               Double Precision:
! scrambled_u   ==> (dimen,Nsobol)  Scrambled Sobol Points (unif dist)
! scrambled_z   ==> (dimen)         Scrambled Sobol Points (norm dist)
! herm          ==> (exc)           Hermite polynomials
! coef          ==> (exc)           Hermite Polynomial Coefficients
! Vmat          ==> (exc,exc)       Potential Matrix 
!               Character:
! seq_in        ==>                 Scrambled Sobol Sequence File Name
!==============================================================================!
!==============================================================================!
implicit none
integer,parameter:: dimen=1
integer :: exc,Nsobol,data_freq
integer :: i,j,k
double precision,allocatable :: scrambled_u(:,:),scrambled_z(:),herm(:),coef(:)
double precision,allocatable :: Vmat(:,:)
character(len=50) :: seq_in
real :: initial_time,final_time
!==============================================================================!
!                               Read Input File
!==============================================================================!
call cpu_time(initial_time)
read(*,*) exc
read(*,*) Nsobol
read(*,*) data_freq
read(*,*) seq_in
allocate(scrambled_u(dimen,Nsobol),scrambled_z(dimen),herm(exc),coef(exc))
allocate(Vmat(exc,exc))
Vmat=0d0
write(*,*) 'Test 1; Successfully Read Input File!'
!==============================================================================!
!                           Read Scrambled Sequence
!==============================================================================!
open(unit=70,file=seq_in)
    read(70,*) scrambled_u
close(70)
!==============================================================================!
!                     Generate Wavefunction Coefficients
!==============================================================================!
coef(1)=1.
coef(2)=1./sqrt(2.)
do i=3,exc
    coef(i)=coef(i-1)*(1./sqrt(2.*(i-1)))
enddo
open(unit=75,file='converge.dat') 
!==============================================================================!
!                           Evaluate Sobol Points 
!==============================================================================!
do i=1,Nsobol              
    call scrambled_sobol_stdnormal(dimen,scrambled_u(:,i),scrambled_z(dimen))
    scrambled_z=scrambled_z/sqrt(2.)      ! factor from normal distribution
!==============================================================================!
!                       Generate Hermite Polynomials
!==============================================================================!
    herm(1)=1.             
    herm(2)=2.*scrambled_z(dimen)       
    do j=3,exc
        herm(j)=(2.*scrambled_z(dimen)*herm(j-1))-(2.*(j-2)*herm(j-2))
    enddo
    herm(:)=herm(:)*coef(:)
!==============================================================================!
!                         Evaluate Matrix Elements
!==============================================================================!
    do j=1,exc 
        do k=1,exc 
            Vmat(j,k)=Vmat(j,k)+herm(j)*herm(k)
        enddo
    enddo
!==============================================================================!
!                           Convergence Analysis
!==============================================================================!
    if (mod(i,data_freq)==0) then
      write(75,*) i,Vmat/i
    endif
enddo
close(75)
Vmat=Vmat/Nsobol 
write(*,*) 'Test 2; Successfully Computed Potential Matrix!'
!==============================================================================!
!                             Final Matrix
!==============================================================================! 
open(unit=80,file='final_matrix.dat')
do i=1,exc
    write(80,*) Vmat(1:exc,i)
enddo
close(80)
call cpu_time(final_time)
write(*,*) 'Total Time:', final_time-initial_time
write(*,*) 'Final Test, Hello Universe!'
!==============================================================================!
!                              output.dat
!==============================================================================!
open(unit=83,file='output.dat')
write(83,*) 'Sobol Points ==> ', Nsobol
write(83,*) 'Maximum Excitation ==> ',exc 
write(83,*) 'Calculation Time (s) ==> ',final_time-initial_time
close(83)
End Program HO_1D_s
