!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!==============================================================================!
!       Discussion:
!DGB analysis for TIP4P single water atom
!==============================================================================!
!       Current Developments:
! implement subspace for calculations  (dimen)
! fix/understnad alpha(r) implementation (currently getting -overlap matrix
!==============================================================================!
!       Modified:
!   27 November 2018
!       Author:
!   Shane Flynn 
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                            Global Paramaters
!==============================================================================!
double precision, parameter:: bohr=0.52917721092
double precision, parameter:: autocm=2.194746313D5
double precision, parameter:: autokcalmol=627.5096 
double precision, parameter:: melectron=1822.88839
double precision, parameter:: Hmass=1.00782503223*melectron
double precision, parameter:: Omass=15.99491461957*melectron
double precision, parameter:: alpha0=1d0
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!       Discussion:
!Natoms             ==> Number of Atoms 
!Dimen              ==> System Dimensionality 
!dimen1             ==> Dimensionality of subspace
!==============================================================================!
integer,parameter::dimen=3
integer :: Natoms, Dimen
character(len=2),allocatable::atom_type(:)
double precision,allocatable::mass(:),sqrt_mass(:)
character(len=5)potential
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
if(atom=='H'.or.atom=='h')then 
    Atom_mass=Hmass
elseif(atom=='O'.or.atom=='o')then 
    Atom_mass=Omass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    stop 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine Hessian(q,H)
!==============================================================================!
!       Discussion:
!Numerically computed Hessian 
!       Variables:
!s          ==> Perturbation parameter 
!H	    ==> (Dimen,Dimen); Symmetrized Mass-Scaled Hessian
!q          ==> (Dimen); XYZ Configuration 
!force	    ==> (Dimen); Forces from Water Potential
!==============================================================================!
implicit none
integer::i,j
double precision::H(dimen,dimen),q(Dimen),r(Dimen),E,force(Dimen),force0(Dimen)
double precision,parameter::s=1d-6
r=q
call water_potential(Natoms/3,r,E,force0)
do i=1,dimen
   r(i)=q(i)+s
   call water_potential(Natoms/3, r, E, force)
   r(i)=q(i)
   do j=1,Dimen
      H(i,j)=(force0(j)-force(j))/s
   enddo
enddo
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
do i=1,Dimen
   do j=1,i
      if(i.ne.j) H(i,j)=(H(i,j)+H(j,i))/2.
      H(i,j)=H(i,j)/(sqrt_mass(i)*sqrt_mass(j)) 
      if(i.ne.j)  H(j,i)=H(i,j)
   enddo
enddo
! write(*,*) 'Mass-Scaled Hessian', H
end subroutine Hessian
!==============================================================================!
subroutine Freq_Hess(Dimen,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!       Variables:
!Hess   ==> (Dimen,Dimen); Hessian Matrix
!omega  ==> (Dimen); Eigenvalues of the Hessian
!U      ==> (Dimen,Dimen); Hessian Eigenvectors
!       LLAPACK (dsyev):
!v      ==> Compute both Eigenvalues and Eigenvectors
!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer :: i,info,lwork,Dimen
double precision :: Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen)
double precision, allocatable :: work(:) 
lwork=max(1,3*Dimen-1)                           !suggested by LAPACK Developers
allocate(work(max(1,lwork)))                     !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',Dimen,U,Dimen,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=Dimen,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
end do
!write(*,*) 'Frequencies From Hess ==> ', omega
end subroutine Freq_Hess
!==============================================================================!
subroutine water_potential(NO,q,energy,force)
use iso_c_binding
use TIP4P_module
!==============================================================================!
!       Discussion:
!External Potential Energy module for water
!       Variables:
!NO         ==> number of water molecules
!q          ==> (9*NO); coordinates
!force      ==> (9*NO) Forces computed in external water potential
!energy     ==> Potential Energy
!==============================================================================!
implicit none
integer,intent(in) :: NO                              
double precision,dimension(9*NO),intent(in)::q   
double precision,dimension(9*NO),intent(inout)::force
double precision,intent(inout)::energy
if(potential=='tip4p') then
   call TIP4P(NO,q,energy,force)
else
   stop 'Check water_potential subroutine'
endif
end subroutine water_potential
!==============================================================================!
end module dgb_groundstate
!==============================================================================!
!==============================================================================!
program main
use dgb_groundstate
!==============================================================================!
!==============================================================================!
!               Discussion:
!coord_in       ==> Input Geometry File-Name (xyz file)
!NG             ==> Number of Gaussians 
!Nsobol         ==> Number of Sobol Points for numerical integration
!ii,skip        ==> Integer(kind=8): necessary for using Sobol Module
!alpha          ==> (NG): Width of Gaussian
!q0             ==> (Dimen): Input Geometry coordinates 
!r              ==> (Dimen,NG): Gaussian Centers 
!r2             ==> (Dimen): ith Gaussian coordinate for integration
!force          ==> (Dimen): Forces. See "Toy_Force" 
!Hess           ==> (Dimen,Dimen): Mass-Scaled Hessian. See "Toy_Hessian"
!omega          ==> (Dimen): Frequencies (Eigenvalues). See "Freq_Hess"
!U              ==> (Dimen,Dimen): Hess Eigenvectors. See "Freq_Hess"
!z              ==> (Dimen,NG): Quasi-Random Sequence for Gaussian Centers
!z2             ==> (Dimen,Nsobol): Quasi-Random Sequence for Integration 
!Smat           ==> (NG,NG): Overlap Matrix
!Tmat           ==> (NG,NG): Kinetic Energy Matrix
!Vmat           ==> (NG,NG): Potential Energy Matrix
!Hmat           ==> (NG,NG): Hamiltonian Matrix (V+T) for eigenvalue problem
!eigenvalues    ==> (NG): Eigenvalues of the Hamiltonian Matrix
!lambda         ==> (NG): Eigenvalues of the Overlap Matrix
!r_ij           ==> center of the i,j matrix element 
!E0             ==> Energy Evaluation at the minimum configuration
!pot_ene        ==> Potential Energy evaluated in q space
!==============================================================================!
implicit none
character(len=50) :: coord_in
integer::NG,Nsobol
integer::i,j,k,l,n,counter
integer*8::ii,skip,skip2
double precision::E0,Tsum,pot_ene,s_sum
!==============================================================================!
double precision, allocatable :: q0(:),force(:),r(:,:),r2(:),eigenvalues(:)
double precision, allocatable :: Hess(:,:),omega(:),U(:,:),z(:,:),z2(:,:)
double precision, allocatable :: Smat(:,:),alpha(:),Tmat(:,:),Vmat(:,:)
double precision, allocatable :: Hmat(:,:),S1mat(:,:),lambda(:),r_ij(:)
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable::work(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) potential
read(*,*) coord_in
read(*,*) NG
read(*,*) Nsobol
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File!'
!==============================================================================!
!                         Set Input Water Geometry 
!Set dimen=2Natoms (2D Potential Hard-Coded) 
!==============================================================================!
open(16,File=coord_in)
read(16,*) Natoms
read(16,*) 
Dimen=3*Natoms 
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen),z(Dimen,NG),z2(Dimen,Nsobol))
allocate(r(Dimen,NG),alpha(NG),Smat(NG,NG),eigenvalues(NG))
allocate(Tmat(NG,NG),Vmat(NG,NG),Hmat(NG,NG),r2(Dimen))
allocate(S1mat(NG,NG),lambda(NG),r_ij(Dimen))
!==============================================================================!
!                         Input Configuration Energy
!Assumes Input Configuration is in Angstroms
!==============================================================================!
do i=1,Natoms
    read(16,*) atom_type(i),q0(3*i-2:3*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(16)
q0=q0/bohr
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Hessian(q0,Hess)
write(*,*) 'Mass-Scaled Hessian ==>'
write(*,*)  Hess
call Freq_Hess(Dimen,Hess,omega,U)
write(*,*) 'Hessian Eigenvalues (omega)==> '
write(*,*) omega
write(*,*) 'Hessian Eigenvectors ==> '
write(*,*) U
write(*,*) 'Test 2; Successfully computed Hessian'
!==============================================================================!
!                 Generate Gaussian Centers with a quasi grid
! Assume centers are in normal-mode space r (not coordinate space q)
!==============================================================================!
do ii=1,NG
    call sobol_stdnormal(Dimen,skip,z(:,ii))
    r(:,ii)=q0(:)
    do i=1,Dimen
        r(:,ii)=r(:,ii)+z(i,ii)*U(:,i)
    enddo
enddo
!==============================================================================!
!                       Generate Alpha Scaling 
!==============================================================================!
do i=1,NG
    alpha(i) = alpha0*exp(-1./dimen*sum(r(:,i)*omega(:)*r(:,i)))
enddo
!open(unit=17,file='centers.dat')
!do i=1,NG
!    write(17,*) r(:,i)
!enddo
!close(17)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
        s_sum=sum(omega(:)*(r(:,i)-r(:,j))**2)
        Smat(i,j)=sqrt(alpha(i)+alpha(j))**(-dimen)&
                 *exp(-0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*s_sum)
        Smat(j,i)=Smat(i,j)
    enddo
enddo
write(*,*) 'Test 4; Successfully computed Overlap Matrix'
!==============================================================================!
!                   Check to see if S is positive definite
! If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
S1mat=Smat
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v', 'u', NG, S1mat, NG, lambda, work, Lwork, info)
write(*,*) 'Info (Overlap Matrix) ===>', info
write(*,*) 'Overlap Matrix Eigenvalues'
open(unit=17,file='overlap_eigenvalues.dat')
do i=1,NG
    write(17,*) lambda(i)
enddo
close(17)
!==============================================================================!
!                           Kinetic Matrix (T)
!==============================================================================!
do i=1,NG
    do j=i,NG
        Tmat(i,j)=Smat(i,j)*0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*&
        sum(omega(:)-(alpha(i)*alpha(j)*(omega(:)**2*(r(:,i)-r(:,j))**2)/(alpha(i)+alpha(j))))
        Tmat(j,i)=Tmat(i,j)
    end do
end do
write(*,*) 'Test 6; Successfully computed Kinetic Matrix'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
! Want to generate a single sequence and then scale for each Gaussian
!==============================================================================!
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
end do
write(*,*) 'Test 7; Successfully generated integration sequence'
!==============================================================================!
!                              Evaluate Potential 
! will need to be replaced with sliding scale
!==============================================================================!
Vmat=0d0
do i=1,NG
   do j=i,NG
      r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
        do l=1,Nsobol
            r2(:)=r_ij(:)+z2(:,l)/(sign(sqrt(abs(omega(:))),omega(:))*(alpha(i)+alpha(j)))
            call water_potential(Natoms/3,(q0+matmul(U,r2)/sqrt_mass(:)),pot_ene,force)
            Vmat(i,j)=Vmat(i,j)+pot_ene
        enddo
            Vmat(j,i)=Vmat(i,j)
    enddo
enddo
Vmat=Vmat*Smat/Nsobol
write(*,*) 'Test 8; Successfully computed Potential Matrix'
!==============================================================================!
!                     Solve Generalized Eigenvalue Problem
!==============================================================================!
Hmat=Vmat
Hmat=Hmat+Tmat
itype=1
eigenvalues=0d0
!==============================================================================!
! Allocations needed if overlap is not diagonalized above
!==============================================================================!
!lwork = max(1,3*NG-1)                                                         !
!!allocate(work(max(1,lwork)))                                                 !
!==============================================================================!
write(*,*) 'info ==>', info
CALL dsygv(itype, 'n', 'u', NG, Hmat, NG, Smat, NG, eigenvalues, work, Lwork, info) 
open(unit=19,file='eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Natoms ==>', Natoms
write(90,*) 'Dimensionality Input System ==>', Dimen 
write(90,*) 'N_gauss ==>', NG
write(90,*) 'N_Sobol ==>', Nsobol
write(90,*) 'omega ==>', omega(1), omega(2)
close(90)
end program main
