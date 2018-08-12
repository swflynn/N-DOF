!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!==============================================================================!
!    Discussion:
!Distributed Gaussian Basis Set (DGB) 
!Computing ground state energy/properties of water
!==============================================================================!
!    Modified:
!       25 July 2018
!    Author:
!       Shane Flynn & Vladimir Mandelshtam
!==============================================================================!
! compute Potential Matrix
!                           Things to Do
! compute matricies on the fly (lambda,kinetic)
! parallel implementation
!==============================================================================!
!==============================================================================!
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                           Global Paramaters
!==============================================================================!
!       Discussion:
!Atomic Properties
!alpha0                 ==> Scaling Parameter for Gaussians
!==============================================================================!
double precision, parameter :: melectron=1822.88839
double precision, parameter :: Hmass=1.00782503223*melectron
double precision, parameter :: Omass=15.99491461957*melectron
double precision, parameter :: bohr=0.52917721092
double precision, parameter :: autocm=2.194746313D5
double precision, parameter :: alpha0=1d0
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
if (atom=='O'.or. atom=='o') then  
    Atom_mass=Omass
else if (atom=='H' .or. atom=='h') then 
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
!Toy Potential Energy for Configuration
!Assume Harmonic Potential: V =  0.5*k*(isum[x(:))^2]
!Assume force constants are all the same (1d0). 
!Assumes you pass in normal mode coordinates, so all dimensions are indepedent
!returns a single value for the potential energy
!==============================================================================!
implicit none
integer :: i
double precision :: x(Dimen),force_con(Dimen),energies
force_con = 1d0
energies=0d0
do i=1,Dimen
    energies=energies+(x(i)**2)*force_con(i)
end do 
energies=energies*.5
!write(*,*) 'Toy Energy from potential'
!write(*,*) energies
end subroutine Toy_Potential
!==============================================================================!
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
!Toy Forces for atoms
!Assume Harmonic Foces: F(:) = - k*x(:)
!Assume force constants are all the same (1d0). 
!==============================================================================!
implicit none
integer :: i
double precision :: x(Dimen),forces(Dimen),force_con(Dimen)
force_con = 1d0
do i=1,Dimen
    forces(i) = -x(i)*force_con(i)
end do 
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
!Assume Harmonic Foces: F(:) = - k*x(:)
!==============================================================================!
implicit none 
integer :: i,j
double precision :: Hess_Mat(Dimen,Dimen),x(Dimen),r(Dimen),force(Dimen)
double precision :: force0(Dimen)
double precision, parameter :: s=1d-6
r=x
call Toy_Force(r, force0)
write(*,*) 'Force0 ==>'
write(*,*) force0
do i=1,Dimen
    r(i)=x(i)+s
    call Toy_Force(r, force)
    r(i)=x(i)
    do j=1,Dimen
        Hess_Mat(i,j)=(force0(j)-force(j))/s
    end do
end do
! symmetrize and mass scale the Hessian
do i=1,Dimen
    do j=1,i
        if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
        Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    end do
end do
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
    write(*,*) omega(i)*autocm, 'normalized = 1?',sum(U(:,i)**2)
end do
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
integer :: i,j,k,n,NG,Nsobol,data_freq
double precision, parameter :: pi=4.0*atan(1d0)
double precision :: lsum,prod_omega,Tsum,pot_ene
integer*8 :: ii,skip,skip2
!==============================================================================!
double precision, allocatable :: q0(:),force(:),force_con(:),q(:,:),q2(:,:),q3(:,:)
double precision, allocatable :: Hess(:,:),omega(:),U(:,:),alpha(:),y(:,:)
double precision, allocatable :: y2(:,:),S(:,:),S1(:,:),T(:,:),Vmat(:,:)
double precision, allocatable :: V1mat(:,:),Hmat(:,:),eigenvalues(:),lmat(:,:),ltest(:)
double precision, allocatable :: Tmat(:,:), val
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
read(*,*) data_freq
write(*,*) 'Number of Gaussians         ==> ', NG
write(*,*) 'Number of Sobol points      ==> ', Nsobol
write(*,*) 'Partial Average Frequency   ==> ', data_freq
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File!'
!==============================================================================!
!                         Set Input Water Geometry 
!==============================================================================!
open(66,File=coord_in)
read(66,*) Natoms
read(66,*) 
Dimen=3*Natoms 
write(*,*) 'Dimensionality ==> ', Dimen
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(force_con(Dimen), Hess(Dimen,Dimen), omega(Dimen),U(Dimen,Dimen))
allocate(alpha(NG),q(Dimen,NG),q2(Dimen,Nsobol),q3(Dimen,Nsobol),y(Dimen,NG),y2(Dimen,Nsobol),S(NG,NG))
allocate(S1(NG,NG),Vmat(NG,NG),V1mat(NG,NG),Hmat(NG,NG),eigenvalues(NG),lmat(NG,NG))
allocate(Tmat(NG,NG))
do i=1,Natoms
    read(66,*) atom_type(i), q0(3*i-2:3*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(3*i-2:3*i)=SQRT(mass(i))
end do
close(66)
q0=q0/bohr                            !Convert xyz Coordinates to Atomic Units
write(*,*) 'q0 ==> '
write(*,*)  q0
write(*,*) 'Test 2; Successfully Read Input Geometry File!'
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0, Hess)
write(*,*) 'Mass-Scaled Hessian ==>'
write(*,*)  Hess
call Frequencies_From_Hess(Dimen,Hess,omega,U)
write(*,*) 'Hessian Eigenvalues (omega)==> '
write(*,*) omega
write(*,*) 'Hessian Eigenvectors (U)==> '
write(*,*) U
!==============================================================================!
do i=1,Dimen
    U(:,i) = U(:,i) / sqrt_mass(:) / sqrt(omega(i))
end do
write(*,*) 'Hessian Normal Modes (Scaled U)> '
write(*,*)  U
write(*,*) 'Test 3; Successfully Determined Normal Modes'
!==============================================================================!
!                       Generate Gaussian Centers (q)
!==============================================================================!
q=0d0
do ii=1,NG
    call sobol_stdnormal(Dimen,skip,y(:,ii))
    q(:,ii)=q0(:)
    do k=1,Dimen
        q(:,ii)=q(:,ii)+y(k,ii)*U(:,k)
    end do
end do
write(*,*) 'Gaussian Centers ==>'
do i=1,NG
write(*,*) q(:,i)
end do
write(*,*) 'Test 4; Successfully Generated Gaussian Centers'
!==============================================================================!
!                       Generate Alpha Scaling 
! A single paramater for each gaussian (the same across all dimensions)
! alpha := omega ===> sqrt{Hess}
!==============================================================================!
alpha=0d0
do i=1,NG
write(*,*) 'omega', omega(:)
write(*,*) 'vec'
write(*,*) q(:,i)
write(*,*) 'vec*omega'
write(*,*)  omega(:)*q(:,i)
write(*,*) 'dot product', dot_product(q(:,i),(omega(:)*q(:,i)))
    alpha(i)=alpha0*NG**(2./dimen)*exp((-1./dimen)*dot_product(q(:,i),(omega(:)*q(:,i))))
write(*,*) 'alpha i', alpha(i)
end do
write(*,*) 'Alpha Scaling ==> '
write(*,*)  alpha
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
write(*,*) 'lmat ===> '
write(*,*) lmat
!==============================================================================!
!                           Overlap Matrix (S)
! scale by lambda when solving generalized eigenvalue problem
!==============================================================================!
S=1d0
prod_omega=1d0
do k=1,Dimen
!    write(*,*) 'omega ', omega(k)
    prod_omega=prod_omega*omega(k)
end do
prod_omega=prod_omega**(-.5)
write(*,*) 'prod omega', prod_omega
do i=1,NG
    do j=i,NG
        S(i,j) = prod_omega*sqrt(2.0*pi/(alpha(i)+alpha(j)))
        S(j,i) = S(i,j)
    end do
end do
write(*,*) 'S ===> '
write(*,*) S

!!!!!!!!!!!!!!!Check to see if S is positive definite!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eigenvalues=0d0
!lwork = max(1,3*NG-1)
!allocate(work(max(1,lwork)))
!write(*,*) 'info ==>', info
!CALL dsyev('n', 'u', NG, S, NG, eigenvalues, work, Lwork, info)
!write(*,*) 'info ==>', info
!write(*,*) 'eigenvalues'
!write(*,*) eigenvalues
!!!!!!!!!!!!!!!Check to see if S is positive definite!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
write(*,*) 'T ===> '
write(*,*) Tmat
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
q2=0d0
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,y2(:,ii))
    q2(:,ii) = q0(:)
    do k=1,Dimen
        q2(:,ii)=q2(:,ii)+y2(k,ii)*U(:,k)
    end do
end do
!write(*,*) 'q2'
!write(*,*) q2
!==============================================================================!
!                              Evaluate Potential 
! loop over each matrix element and run 1 ---> Nsobol
! q0 is over U/omega, need q2 over omega*alpha scaling
! q2 is sequence of len. Nsobol, and distributed to q0
! q3 takes q2 and scales by alpha
!==============================================================================!
Vmat=0d0
!V1mat=0d0
do i=1,NG
    do j=i,NG
        pot_ene=0d0
        do n=1,Nsobol
            q3 = q
!            write(*,*) 'q3 test'
!            write(*,*) q3
            q3(:,n) = q3(:,n) * (alpha(i) + alpha(j))**(-.5)* ((alpha(i)*q(:,i) + alpha(j)*q(:,j)) / (alpha(i) + alpha(j)))
!            write(*,*) 'q3 test again'
!            write(*,*) q3
            call Toy_Potential(q3(:,n),pot_ene)
            Vmat(i,j)=Vmat(i,j)+pot_ene
!            Vmat(j,i)=Vmat(i,j)
        end do
    end do
end do
!==============================================================================!
!                           Partial Average and Energy
! Check 'INFO' if you have an error, info>1 then overlap matrix isn't positive definite
! scale matricies by lambda at the end
!==============================================================================!
do i=1,NG
do  j=1,NG
Vmat(j,i) = Vmat(i,j)
end do
end do
write(*,*) 'Vmat'
write(*,*) Vmat
do i=1,NG
    do j=i,NG
    Hmat(i,j)=Vmat(i,j)/Nsobol
    Hmat(j,i) = Hmat(i,j)
    end do
end do
Hmat=Hmat+Tmat
S1=S
do i=1,NG
    do j=i,NG
        Hmat(i,j)=Hmat(i,j)*exp(-lmat(i,j))
        Hmat(j,i) = Hmat(i,j)
        S1(i,j)=S1(i,j)*exp(-lmat(i,j))
        S1(j,i) = S1(i,j)
    end do 
end do
write(*,*) 'Hamiltonian ===> '
write(*,*) Hmat
write(*,*) 'Scaled S1 ===> '
write(*,*) S1

!val = S1(1,1)
!write(*,*) 'testing S remove when done'
!do i=1,NG
!do j=1,NG
!S1(i,j) = S1(i,j) / val
!end do 
!end do
!write(*,*) 'test S1'
!write(*,*) S1

! DSYGV allocations
itype=1
eigenvalues=0d0
lwork = max(1,3*NG-1)
allocate(work(max(1,lwork)))
write(*,*) 'info ==>', info
CALL dsygv(itype, 'n', 'u', NG, Hmat, NG, S1, NG, eigenvalues, work, Lwork, info)
write(*,*) 'info after ==>', info
write(*,*) 'eigenvalues'
write(*,*) eigenvalues
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Natoms ==>', Natoms
write(90,*) 'Dimensionality Input System ==>', Dimen
write(90,*) 'N_gauss ==>', NG
close(90)
end program DGB_ground_state
