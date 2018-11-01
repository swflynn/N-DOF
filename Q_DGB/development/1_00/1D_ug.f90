!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!==============================================================================!
!    Discussion:
!Calculation for a 1D system only, single atom x coordinates, mass=1
!uniform grid, new derivation using normal mode coordinates
!==============================================================================!
!    Modified:
!       25 October 2018
!    Author:
!       Shane Flynn 
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
double precision, parameter :: Hmass=1d0
double precision, parameter :: alpha0=1d0
double precision, parameter :: pi=4.*atan(1d0)
double precision, parameter :: c1=1d0
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
! V := 1/2 c1(x)^2 
! 1D potential hard-coded
!==============================================================================!
implicit none
double precision :: x(Dimen),energies
energies=0d0
energies = c1*.5*(x(1))**2 
!write(*,*) 'Toy Energy from potential subroutine'
!write(*,*) energies
end subroutine Toy_Potential
!==============================================================================!
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
! Toy Forces, hardcoded for 1D potential
!==============================================================================!
implicit none
integer :: i
double precision :: x(Dimen),forces(Dimen)
forces(1) = -(x(1)*c1) 
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
integer :: NG,Nsobol
integer :: i,j,k,l,n,counter,iter
integer*8 :: ii,skip,skip2
double precision :: E0,detO,Asum,prod_omega,Tsum,pot_ene,alpha_par,Tpre
double precision :: low_bound, renorm
!==============================================================================!
double precision, allocatable :: q0(:),force(:),r(:,:),r2(:,:)
double precision, allocatable :: Hess(:,:),omega(:),U(:,:),z(:,:),z2(:,:),alpha(:)
double precision, allocatable :: Smat(:,:),S1mat(:,:),eigenvalues(:)
double precision, allocatable :: Tmat(:,:),Vmat(:,:),Hmat(:,:)
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
read(*,*) low_bound
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
Dimen=1*Natoms 
write(*,*) 'Dimensionality ==> ', Dimen
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen),z2(Dimen,Nsobol))
allocate(r(Dimen,NG),z(Dimen,Nsobol),alpha(NG),Smat(NG,NG),S1mat(NG,NG))
allocate(eigenvalues(NG),Tmat(NG,NG),Vmat(NG,NG),Hmat(NG,NG),r2(Dimen,Nsobol))
do i=1,Natoms
    read(66,*) atom_type(i), q0(1)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(1)=sqrt(mass(i))
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
! Scale Eigenvectors by mass
!==============================================================================!
do i=1,Dimen
    U(:,i)=U(:,i)*sqrt_mass(:) 
end do
!==============================================================================!
!                       Generate Gaussian Centers (q^i) use uniform grid
! generate coordinates in coordinate space q, then scale to normal modes r  
!==============================================================================!
r=0d0
open(unit=17,file='centers.dat')
r(:,1) = low_bound
write(17,*) r(:,1), 1
do ii=2,NG
    r(:,ii) = r(:,ii-1) + (-low_bound - low_bound)/NG
    write(17,*) r(:,ii), 1
end do
!==============================================================================!
! replace q with r, go to normal mode coordinates for analysis
! r = U^Tsqrt(M)(q-q0)
!==============================================================================!
do ii=1,NG
    do k=1,Dimen
        r(:,ii) = U(:,k)*(r(:,ii)-q0(:))
    enddo
enddo
close(17)
write(*,*) 'r='
write(*,*) r
!==============================================================================!
!                       Generate Alpha Scaling 
! A single paramater for each gaussian (the same across all dimensions)
! set to a value for now, no transformations, can generalize later
!==============================================================================!
alpha=alpha_par
write(*,*) 'alpha scaling'
write(*,*) alpha
!==============================================================================!
!                           Overlap Matrix (S)
! Checked numerically for 1D case, this is working
!==============================================================================!
Smat=1d0
do k=1,dimen
    do i=1,NG
        do j=i,NG
            Smat(i,j) = smat(i,j) * sqrt(2.*pi/((alpha(i)+alpha(j))*omega(k)))*&
                exp(-alpha(i)*alpha(j)*omega(k)*(r(k,i)-r(k,j))**2./(2.*(alpha(i)+alpha(j))))
            Smat(j,i)=Smat(i,j)
        enddo
    enddo
enddo
write(*,*) 'smat'
write(*,*) smat
!==============================================================================!
!                   Check to see if S is positive definite
! If this is removed, you need to allocate llapack arrays at hamiltonian matrix
!==============================================================================!
S1mat=Smat
eigenvalues=0d0
lwork = max(1,3*NG-1)
allocate(work(max(1,lwork)))
write(*,*) 'info ==>', info
CALL dsyev('n', 'u', NG, S1mat, NG, eigenvalues, work, Lwork, info)
write(*,*) 'info ==>', info
write(*,*) 'eigenvalues for Overlap matrix'
do i=1,NG
    write(*,*) eigenvalues(i)
end do
write(*,*) 'working test'
!==============================================================================!
!                           Kinetic Matrix (T)
! Checked numerically for 1D case, this is working
!==============================================================================!
Tmat=0d0
Tpre=0d0
Tsum=0d0
do i=1,NG
    do j=i,NG
        Tpre = alpha(i)*alpha(j)/(2.*(alpha(i)+alpha(j)))
        Tsum=0d0
        do k=1,Dimen
            Tsum=Tsum + (omega(k) - (alpha(i)*alpha(j)*omega(k)**2.*(r(k,i)-r(k,j))**2./(alpha(i)+alpha(j))))
        end do
        Tmat(i,j) = Smat(i,j)*Tpre*Tsum
        Tmat(j,i) = Tmat(i,j)
    end do
end do
write(*,*) 'Kinetic Matrix ==>'
write(*,*) Tmat
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
! Want to generate a single sequence and then scale for each Gaussian
!==============================================================================!
z=0d0
r2=0d0
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
end do
!==============================================================================!
!                              Evaluate Potential 
! I am here
! add indicies, check by hand with 1 d case and 2 gaussians
!==============================================================================!
Vmat=0d0
do i=1,NG
    do j=i,NG
        do l=1,Nsobol
            r2(:,l) = z2(:,l) * (omega(:)*(alpha(i)+alpha(j)))**(-0.5) + &
                (alpha(i)*r(:,i)+alpha(j)*r(:,j)/(alpha(i)+alpha(j)))
            call Toy_Potential((q0+(matmul(U,r2(:,l))/sqrt_mass(:))),pot_ene)
            Vmat(i,j)=Vmat(i,j)+pot_ene
            Vmat(j,i)=Vmat(i,j)
        enddo
    enddo
enddo
Vmat=Vmat*Smat/Nsobol
write(*,*) 'vmat'
write(*,*) Vmat
!==============================================================================!
!                     Solve Generalized Eigenvalue Problem
!==============================================================================!
Hmat=Vmat
Hmat=Hmat+Tmat
S1mat=Smat
itype=1
eigenvalues=0d0
!==============================================================================!
! Allocations needed if overlap is not diagonalized above
!==============================================================================!
!lwork = max(1,3*NG-1)
!!allocate(work(max(1,lwork)))
!==============================================================================!
write(*,*) 'info ==>', info
CALL dsygv(itype, 'n', 'u', NG, Hmat, NG, S1mat, NG, eigenvalues, work, Lwork, info)
write(*,*) 'info after ==>', info
write(*,*) 'eigenvalues'
write(*,*) eigenvalues
write(*,*) 'Computed,                   True Energy,                    %Error'
do i=1,NG
    write(*,*) eigenvalues(i), ((i-1)+.5)*omega(1), ((eigenvalues(i) - ((i-1)+.5)*omega(1))/(((i-1)+.5)*omega(1)))*100
end do 
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Natoms ==>', Natoms
write(90,*) 'Dimensionality Input System ==>', Dimen 
write(90,*) 'N_gauss ==>', NG
write(*,*) 'Alpha Parameter ==>', alpha_par
close(90)
end program DGB_ground_state
