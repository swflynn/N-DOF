!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!       02
!==============================================================================!
!    Discussion:
!Distributed Gaussian Basis Set (DGB) 
!Calculation for a 2D system only, single atom x,y coordinates, mass=1
!==============================================================================!
!    Modified:
!       28 August 2018
!    Author:
!       Shane Flynn & Vladimir Mandelshtam
!==============================================================================!
!   Things To Do:
! compute potential matrix
! fix prefactor for gaussian centers
! Compute Kinetic on the fly 
! Compute Lambda on the Fly
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                           Global Paramaters
!==============================================================================!
!       Discussion:
!Atomic Properties
!alpha0                 ==> Scaling Parameter for Gaussians
!c1,c2          constants for 2D potetnial, c1 > c2 to have positive Hessian
!==============================================================================!
double precision, parameter :: Hmass=1d0
double precision, parameter :: alpha0=1d0
double precision, parameter :: pi=4.0*atan(1d0)
double precision, parameter :: c1=2d0
double precision, parameter :: c2=1d0
!==============================================================================!
!                               Global Variables 
!==============================================================================!
!       Discussion:
!Natoms             ==> Number of atoms 
!Dimen              ==> Dimensionality of Input System 
!atom_type          ==> Element Abbreviation
!mass(Natoms)       ==> Atom Masses
!sqrt_mass(Dimen)   ==> Square Root Mass
!==============================================================================!
integer :: Natoms, Dimen
character(len=2), allocatable :: atom_type(:)
double precision, allocatable :: mass(:),sqrt_mass(:)
!==============================================================================!
!                           Begin Module 
!==============================================================================!
contains
!==============================================================================!
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!       Discussion:
!Computes the mass for each atom 
!Currently defined for water only
!==============================================================================!
implicit none
double precision :: Atom_Mass
character(len=2) :: atom
if (atom=='H' .or. atom=='h') then 
    Atom_mass=Hmass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    write(*,*) 'Check Atom_Mass Function in the dgb module'
    STOP 
end if
end function Atom_Mass
!==============================================================================!
subroutine Toy_Potential(x,energies)
!==============================================================================!
!       Discussion:
! V := c1(x+y)^2 + c2(x-y)^2
! 2D potential hard-coded
!==============================================================================!
implicit none
double precision :: x(Dimen),energies
energies=0d0
energies = c1*(x(1)+x(2))**2 + c2*(x(1)-x(2))**2
!write(*,*) 'Toy Energy from potential subroutine'
!write(*,*) energies
end subroutine Toy_Potential
!==============================================================================!
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
! Toy Forces, hardcoded for 2D potential
!==============================================================================!
implicit none
integer :: i
double precision :: x(Dimen),forces(Dimen)
forces(1) = -(2.*x(1)*(c1+c2) + 2.*x(2)*(c1-c2))
forces(2) = -(2.*x(1)*(c1-c2) + 2.*x(2)*(c1+c2))
write(*,*) 'Toy Forces from force subroutine'
write(*,*) forces
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
! Computes Hessian From Forces
!==============================================================================!
implicit none 
integer :: i,j
double precision :: Hess_Mat(Dimen,Dimen),x(Dimen),r(Dimen),force(Dimen)
double precision :: force0(Dimen)
double precision, parameter :: s=1d-6
r=x
call Toy_Force(r, force0)
write(*,*) 'Force0 Hessian ==>'
write(*,*) force0
do i=1,Dimen
    r(i)=x(i)+s
    call Toy_Force(r, force)
    r(i)=x(i)
    do j=1,Dimen
        Hess_Mat(i,j)=(force0(j)-force(j))/s
    end do
end do
write(*,*) 'Hessian matrix subroutine ==>'
write(*,*) Hess_Mat
! symmetrize and mass scale the Hessian
do i=1,Dimen
    do j=1,i
        if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
        Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    end do
end do
!write(*,*) 'symmetrized-mass-scaled Hessian from subroutine'
!write(*,*) Hess_Mat
end subroutine Toy_Hessian
!==============================================================================!
subroutine Frequencies_From_Hess(Dimen,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Computes the Eigenvalues and Eigenvectors of the Hessian
!dsyev subroutine (real symmetrix eigen-solver) from LAPACK is used.
!       Arguments:
!Dimen                  ==> Dimensionality of System
!Hess(Dimen,Dimen)      ==> Hessian 
!omega(Dimen)           ==> Eigenvalues of the Hessian
!U(Dimen,Dimen)         ==> Eigenvectors of the Hessian
!       Output:
!omega                  ==> Frequencies (Eigenvalues)
!U                      ==> Normal Modes (Eigenvectors)
!       Other:
!dsyev                  ==> v: eigenvalues and eigenvectors
!                           u: upper triangle of matrix
!==============================================================================!
implicit none
integer :: i, info, lwork, Dimen
double precision :: Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen)
double precision, allocatable :: work(:) 
lwork = max(1,3*Dimen-1)        !suggested by LAPACK Developers
allocate(work(max(1,lwork)))    !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',Dimen,U,Dimen,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=Dimen,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?',sum(U(:,i)**2)
end do
!write(*,*) 'Frequencies from hessian subroutine'
!write(*,*) omega
end subroutine frequencies_from_Hess
!==============================================================================!
!==============================================================================!
!==============================================================================!
end module dgb_groundstate
!==============================================================================!
!==============================================================================!
program DGB_ground_state
use dgb_groundstate
!==============================================================================!
!==============================================================================!
!               Discussion:
! coor_in           ==> Input Water Geometry 
! NG                ==> Number of Gaussians 
! Nsobol            ==> Number of Sobol Points for numerical integration
! ii, skip          ==> Integer; (kind=8) necessary for using Sobol Module
! q0(Dimen)         ==> Input Water Geometry x,y,z coordinates 
!		        Assumes Input Coordinates in Angstroms
! q(dimen, NG)      ==> Gaussian Centers (distributed according to Ground State)
! q2(dimen, Nsobol) ==> Sequence for potential evaluation (distrbuted to ground state)
! q3(dimen, Nsobol) ==> Sequence for potential sacaled by alpha factors for Vij
! force(Dimen)      ==> Forces associated with atoms
! omega(Dimen)      ==> Eigenvalues of the Hessian (frequencies)
! U(Dimen,Dimen)    ==> Normal Modes from Hessian
! y(Dimen,NG)       ==> Sobol points for Gaussian Centers
! y2(Dimen,integ)   ==> Sobol Sequence for computing Matrices
! S(NG,NG)          ==> Overlap matrix for Gaussians
! S1(NG,NG)         ==> Overlap matrix for eigenvalues (destroyed each iteration)
! Vmat(NG,NG)       ==> Potential Energy Matrix
! V1mat(NG,NG)      ==> Potential Energy Matrix Partial Average
! Hmat(NG,NG)       ==> Hamiltonian Matrix (V+T) for eigenvalue problem
! eigenvalues(NG)   ==> Eigenvalues of the Hamiltonian Matrix
!==============================================================================!
implicit none
character(len=50) :: coord_in
integer :: NG,Nsobol,data_freq
integer :: i,j,k,n,counter
integer*8 :: ii,skip,skip2
double precision :: E0,detO,lsum,prod_omega,Tsum,pot_ene,alpha_par
!integer :: data_freq
!==============================================================================!
double precision, allocatable :: q0(:),force(:),q(:,:),q2(:,:)
double precision, allocatable :: Hess(:,:),omega(:),U(:,:),z(:,:),z2(:,:),alpha(:)
double precision, allocatable :: lmat(:,:),Smat(:,:),S1mat(:,:),eigenvalues(:)
double precision, allocatable :: Tmat(:,:),Vmat(:,:),Hmat(:,:)
!double precision, allocatable :: q3(:,:)
!double precision, allocatable :: V1mat(:,:),ltest(:)
!double precision, allocatable :: val
!==============================================================================!
!                           dsygv variables                                    !
integer :: itype, info, lwork
double precision, allocatable :: work(:)
!==============================================================================!
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) coord_in
read(*,*) NG
read(*,*) Nsobol
read(*,*) alpha_par
read(*,*) data_freq
!read(*,*) data_freq
write(*,*) 'Number of Gaussians         ==> ', NG
write(*,*) 'Number of Sobol points      ==> ', Nsobol
!write(*,*) 'Partial Average Frequency   ==> ', data_freq
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File!'
!==============================================================================!
!                         Set Input Water Geometry 
! set dim=2*Natoms for x-y
!==============================================================================!
open(66,File=coord_in)
read(66,*) Natoms
read(66,*) 
Dimen=2*Natoms 
write(*,*) 'Dimensionality ==> ', Dimen
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen))
allocate(q(Dimen,NG),z(Dimen,NG),z2(Dimen,Nsobol),alpha(NG),lmat(NG,NG),Smat(NG,NG),S1mat(NG,NG))
allocate(eigenvalues(NG),Tmat(NG,NG),Vmat(NG,NG),Hmat(NG,NG),q2(Dimen,Nsobol))
!allocate(q3(Dimen,Nsobol, V1mat(NG,NG))
do i=1,Natoms
    read(66,*) atom_type(i), q0(2*i-1:2*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(2*i-1:2*i)=SQRT(mass(i))
end do
close(66)
write(*,*) 'q0 ==> '
write(*,*)  q0
write(*,*) 'Test 2; Successfully Read Input Geometry File!'
!==============================================================================!
!                   Input Configuration Energy
!==============================================================================!
call toy_potential(q0,E0)
write(*,*) 'E0 ==> ', E0
write(*,*) 'Test 3; Successfully Evaluate Input Geometry Energy!'
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0, Hess)
write(*,*) 'Mass-Scaled Hessian ==>'
write(*,*)  Hess
call Frequencies_From_Hess(Dimen,Hess,omega,U)
write(*,*) 'Hessian Eigenvalues (omega)==> '
write(*,*) omega
write(*,*) 'Hessian Eigenvectors ==> '
write(*,*) U
!==============================================================================!
! Scale Eigenvectors for std normal distribution
!==============================================================================!
do i=1,Dimen
    U(:,i) = U(:,i) / sqrt_mass(:) / sqrt(omega(i))
end do
write(*,*) 'Scaled Eigenvectors U ==> '
write(*,*)  U
write(*,*) 'Test 3; Successfully Determined Normal Modes'
!==============================================================================!
!                       Generate Gaussian Centers (q^i)
!==============================================================================!
z=0d0
q=0d0
do ii=1,NG
    call sobol_stdnormal(Dimen,skip,z(:,ii))
    q(:,ii)=q0(:)
    do i=1,Dimen
        q(:,ii)=q(:,ii)+z(i,ii)*U(:,i)
    end do
end do
!==============================================================================!
!                   Compute Ground State Pre-factor
! Not sure if this is correct, need to check this
!==============================================================================! 
detO=1d0
do i=1,Dimen
    detO=detO*(Hess(i,i))
enddo
q = q*sqrt(detO/(2*pi)**dimen)
write(*,*) 'Test 4; Successfully Generated Gaussian Centers'
!==============================================================================!
!                       Generate Alpha Scaling 
! A single paramater for each gaussian (the same across all dimensions)
! set to a value for now, no transformations, can generalize later, fo rnow just 
! use more gaussians
!==============================================================================!
alpha=alpha_par
write(*,*) 'alpha scaling'
write(*,*) alpha
!==============================================================================!
!                       Lambda matrix lmat(NG,NG)
! using for testing
! once done re-write this to compute on the fly
! Only scale by lambda at the end when you are passing in values to Hamiltonian
!==============================================================================!
do i=1,NG
    do j=i,NG
        lsum=0d0
        do k=1,dimen
            lsum=lsum+omega(k)*(q(k,i)-q(k,j))**2
        end do
        lmat(i,j) = alpha(i)*alpha(j)/(2*(alpha(i)+alpha(j)))*lsum
        lmat(j,i) = lmat(i,j)
    end do
end do
lmat=exp(-lmat)
write(*,*) 'Here is lambda'
write(*,*) lmat
!==============================================================================!
!                           Overlap Matrix (S)
! scale by lambda when solving generalized eigenvalue problem
!==============================================================================!
Smat=1d0
prod_omega=1d0
do k=1,Dimen
    prod_omega=prod_omega*omega(k)
end do
prod_omega=prod_omega**(-.5)
do i=1,NG
    do j=i,NG
        Smat(i,j) = prod_omega*sqrt(2.0*pi/(alpha(i)+alpha(j)))
        Smat(j,i) = Smat(i,j)
    end do
end do
!==============================================================================!
!!!!!!!!!!!!!!!Check to see if S is positive definite!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing purposes only dont need this here
!==============================================================================!
S1mat=Smat
S1mat=S1mat*Lmat
eigenvalues=0d0
lwork = max(1,3*NG-1)
allocate(work(max(1,lwork)))
write(*,*) 'info ==>', info
CALL dsyev('n', 'u', NG, S1mat, NG, eigenvalues, work, Lwork, info)
write(*,*) 'info ==>', info
write(*,*) 'eigenvalues for Overlap matrix'
write(*,*) eigenvalues
!==============================================================================!
!                           Kinetic Matrix (T)
!==============================================================================!
Tmat = 0d0
do i=1,NG
    do j=i,NG
        Tsum=0d0
        do k=1,Dimen
            Tsum=Tsum+(omega(k)*(1.-(2.*alpha(i)*alpha(j)/(alpha(i)+alpha(j)))*omega(k)*(q(k,i)-q(k,j))**2))
        end do
        Tmat(i,j) = -prod_omega*Tsum*(alpha(i)*(alpha(j)*(2*pi)**(dimen/2.))/(4.*(alpha(i)+alpha(j))**(1.+dimen/2.)))
        Tmat(j,i) = Tmat(i,j)
    end do
end do
write(*,*) 'Kinetic Matrix'
write(*,*) Tmat
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
z2=0d0
q2=0d0
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
    do i=1,Dimen
        q2(:,ii)=q2(:,ii)+z2(i,ii)*U(:,i)
    end do
end do
!==============================================================================!
!                              Evaluate Potential 
!==============================================================================!
Vmat=0d0
do i=1,NG
    do j=i,NG
        pot_ene=0d0
        do n=1,Nsobol
            call Toy_Potential(q2(:,n)+((alpha(i)*q(:,i)+alpha(j)*q(:,j))/(alpha(i)+alpha(j))),pot_ene)
            Vmat(i,j)=Vmat(i,j)*(1/sqrt((alpha(i)+alpha(j))))+pot_ene
            Vmat(j,i)=Vmat(i,j)
        end do
    end do
end do
Vmat=Vmat/(Nsobol*prod_omega)
!==============================================================================!
Hmat=Vmat
Hmat=Hmat+Tmat
Hmat=Hmat*Lmat
S1mat=Smat
S1mat=S1mat*Lmat
write(*,*) 'Overlap matrix'
write(*,*) S1mat

itype=1
eigenvalues=0d0
lwork = max(1,3*NG-1)
!allocate(work(max(1,lwork)))
write(*,*) 'info ==>', info
CALL dsygv(itype, 'n', 'u', NG, Hmat, NG, S1mat, NG, eigenvalues, work, Lwork, info)
write(*,*) 'info after ==>', info
write(*,*) 'eigenvalues'
write(*,*) eigenvalues

write(*,*) '2D E0'
write(*,*) .5*(omega(1)+omega(2))
!!==============================================================================!
!!                           Partial Average and Energy
!! Check 'INFO' if you have an error, info>1 then overlap matrix isn't positive definite
!! scale matricies by lambda at the end
!!==============================================================================!
!do i=1,NG
    !do  j=1,NG
        !Vmat(j,i) = Vmat(i,j)
    !end do
!end do
!do i=1,NG
!    do j=i,NG
!    Hmat(i,j)=Vmat(i,j)/Nsobol
!    Hmat(j,i) = Hmat(i,j)
!    end do
!end do
!Hmat=Hmat+Tmat
!S1=S
!do i=1,NG
!    do j=i,NG
!        Hmat(i,j)=Hmat(i,j)*exp(-lmat(i,j))
!        Hmat(j,i) = Hmat(i,j)
!        S1(i,j)=S1(i,j)*exp(-lmat(i,j))
!        S1(j,i) = S1(i,j)
!    end do 
!end do
!write(*,*) 'Hamiltonian ===> '
!write(*,*) Hmat
!write(*,*) 'Scaled S1 ===> '
!write(*,*) S1
!
!!val = S1(1,1)
!!write(*,*) 'testing S remove when done'
!!do i=1,NG
!!do j=1,NG
!!S1(i,j) = S1(i,j) / val
!!end do 
!!end do
!!write(*,*) 'test S1'
!!write(*,*) S1
!
!! DSYGV allocations
!itype=1
!eigenvalues=0d0
!lwork = max(1,3*NG-1)
!allocate(work(max(1,lwork)))
!write(*,*) 'info ==>', info
!CALL dsygv(itype, 'n', 'u', NG, Hmat, NG, S1, NG, eigenvalues, work, Lwork, info)
!write(*,*) 'info after ==>', info
!write(*,*) 'eigenvalues'
!write(*,*) eigenvalues
!!==============================================================================!
!!                               output file                                    !
!!==============================================================================!
!open(90,file='simulation.dat')
!write(90,*) 'Natoms ==>', Natoms
!write(90,*) 'Dimensionality Input System ==>', Dimen !write(90,*) 'N_gauss ==>', NG
!close(90)
end program DGB_ground_state
