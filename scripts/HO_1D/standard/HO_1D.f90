!=============================================================================80
!                           1D Harmonic Oscillator
!=============================================================================80
!		Discussion:
! Quasi Monte Carlo Integration, computing the Potential Energy Matrix 
! for a 1-Dimensional Harmonic Oscillator (Hermite Polynomials)
! The basis functions are orthonormal, therefore off-diagonal=0,diagonal=1.
!==============================================================================!
! Code requires the sobol module to generate the sobol sequence (sobol.f90)
!==============================================================================!
!		Modified:
! 20 March 2017
!		Author:
! Shane Flynn
!==============================================================================!
Program HO_1D
use sobol
!==============================================================================!
!                                 Variables
!==============================================================================!
!               Integer:
! d         ==> Spatial Dimension (parameter=1 for this code)
! exc       ==> maximum excitation (herminte polynomial)
! Nsobol    ==> Number of Sobol Points for MC Integration
! data_freq ==> Interval for Convergence Analysis
! skip      ==> Initialize Sobol Sequence (suggested=Nsobol)
!               Double Precision:
! norm      ==> (d,Nsobol)  Sobol Sequence
! herm      ==> (exc)       Hermite polynomials
! coef      ==> (exc)       Hermite Polynomial Coefficients
! Vmat      ==> (exc,exc)   Potential Matrix 
!==============================================================================!
implicit none
integer,parameter:: d=1
integer :: exc,Nsobol,data_freq
integer :: i,j,k
integer*8 :: skip                       ! Must be *8 for sobol.f90 module 
double precision,allocatable :: norm(:,:),herm(:),coef(:),Vmat(:,:)
real :: initial_time,final_time
!==============================================================================!
!                               Read Input File
!==============================================================================!
call cpu_time(initial_time)
read(*,*) exc
read(*,*) Nsobol
read(*,*) data_freq
skip=Nsobol
allocate(norm(d,Nsobol),herm(exc),coef(exc),Vmat(exc,exc))
Vmat=0d0
write(*,*) 'Test 1; Successfully Read Input File!'
!==============================================================================!
!                           Generate Sobol Sequence 
!==============================================================================!
do i=1,Nsobol                
    call sobol_stdnormal(d,skip,norm(:,i))
enddo
norm=norm/sqrt(2.)      ! factor from the normal distribution
!==============================================================================!
!                     Generate Wavefunction Coefficients 
!==============================================================================!
coef(1)=1
coef(2)=1./(sqrt(2.))
do i=3,exc
    coef(i)=coef(i-1)*(1 /sqrt(2.*real(i-1)))
enddo
open(unit=9,file='converge.dat')
!==============================================================================!
!                           Evaluate Sobol Points 
!==============================================================================!
do i=1,Nsobol              
!==============================================================================!
!                       Generate Hermite Polynomials
!==============================================================================!
    herm(1)=1.             
    herm(2)=2.*norm(1,i)       
    do j=3,exc
         herm(j)=(2.*norm(d,i)*herm(j-1))-(2.*(j-2)*herm(j-2))
    enddo
!==============================================================================!
!                         Evaluate Matrix Elements
!==============================================================================!
    do j=1,exc     
        do k=1,exc 
            Vmat(j,k)=Vmat(j,k)+coef(j)*herm(j)*coef(k)*herm(k)
        enddo
    enddo
!==============================================================================!
!                           Convergence Analysis
!==============================================================================!
    if(mod(i,data_freq)==0) then
        write(9,*) i, Vmat/real(i)
    endif
enddo
close(9)
Vmat=Vmat/Nsobol
write(*,*) 'Test 2; Successfully Computed Potential Matrix!'
!==============================================================================!
!                             Final Matrix
!==============================================================================!
open(unit=10,file='final_matrix.dat')
do i=1,exc
    write(10,*) Vmat(1:exc,i)
enddo
close(10)
call cpu_time(final_time)
WRITE(*,*) 'Total Time:', final_time-initial_time
write(*,*) 'Final Test, Hello Universe!'
!==============================================================================!
!                              output.dat
!==============================================================================!
open(unit=83,file='output.dat')
write(83,*) 'Sobol Points = ', Nsobol
write(83,*) 'Maximum Excitation= ',exc 
write(83,*) 'Calculation Time (s): ', final_time-initial_time
close(unit=83)
End Program HO_1D
