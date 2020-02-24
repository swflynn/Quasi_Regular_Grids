!=============================================================================80
!                       2D Morse with Quasi-Monte Carlo
!==============================================================================!
!       Discussion:
!Quasi-Regular Gaussian Basis
!using optimized grid solve generalized eigenvalue problem
!this code computes the potential using quasi Monte Carlo (not GHQ)
!note, this has not been tested for 3D morse, may not have correct dimen-scaling
!==============================================================================!
!       Modified:
!   1 May 2019
!       Author:
!   Shane Flynn
!==============================================================================!
module QRGB_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==>i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!integral_P     ==>Normalization factor for P(x) (from Grid Generation)
!E_cut          ==>Energy Cutoff Contour    (be consistent with Grid Generation)
!c_LJ           ==>Lennard-Jones parameter  (be consistent with Grid Generation)
!==============================================================================!
integer::d,NG
double precision::E_cut,c_LJ,integral_P
!==============================================================================!
contains
!==============================================================================!
function P(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Delta          ==>Delta parameter for P(x) distribution, :=10% of Ecut
!==============================================================================!
implicit none
double precision::Delta,x(d),P
Delta=0.01*E_cut
if(Potential(x)<E_cut) then
   P=(E_cut+Delta-Potential(x))**(d/2.)/integral_P
else                                              !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
!==============================================================================!
function Potential(x)
!==============================================================================!
!Hard-coded Morse Potential Energy
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Potential      ==>evaluate Potential(x)
!D_morse        ==>Parameter for Morse Potential
!c_x/c_y        ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),Potential
double precision,parameter::D_morse=12.
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
Potential=D_morse*((exp(-c_x*x(1))-1)**2+(exp(-c_y*x(2))-1)**2)
end function Potential
!==============================================================================!
end module QRGB_mod
!==============================================================================!
!==============================================================================!
program main
use QRGB_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!grid_in        ==>Filename Containing Gridpoints, see qlj_morse2d_grid.f90
!theory_in      ==>Filename Containing Analytic Eigenvalues, see theory.f90
!Nsobol         ==>Number of sobol points for qMC integration
!skip           ==>Seed for sobol generator
!NG             ==>Number of Gaussian Basis Functions (gridpoints)
!alpha0         ==>Flat Scaling Parameter for Gaussian Widths
!alpha          ==>(d) Gaussian Widths
!RCN            ==>Recriprical Convergence Number, stability of Overlap Matrix
!x              ==>(d) ith atoms coordinates
!x_ij           ==>i-jth gaussian center (product of Gaussians is a Gaussian)
!Smat           ==>(NG,NG) Overlap Matrix
!Vmat           ==>(NG,NG) Potential Matrix
!Hmat           ==>(NG,NG) Hamiltonian Matrix
!z              ==>sequence for matrix elements integration
!data_freq      ==>Interval to evaluate qMC integration
!==============================================================================!
implicit none
character(len=50)::grid_in,theory_in
integer::Nsobol,data_freq,counter,i,j,l,n
integer*8::skip
double precision::aij,r2,alpha0
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,theory
double precision,allocatable,dimension(:,:)::x,Smat,z,Vmat,Hmat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) d
read(*,*) NG
read(*,*) Nsobol
read(*,*) alpha0
read(*,*) integral_P
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) grid_in
read(*,*) theory_in
read(*,*) data_freq
Nsobol=Nsobol/data_freq
Nsobol=Nsobol*data_freq
write(*,*) 'data_freq=',data_freq, '  Nsobol=',Nsobol
skip=Nsobol
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),alpha(NG),eigenvalues(NG),Smat(NG,NG),Vmat(NG,NG))
allocate(z(d,Nsobol),Hmat(NG,NG),theory(NG))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Input GridPoints x(d,NG)
!==============================================================================!
open(17,File=grid_in)
do n=1,NG
    read(17,*) x(:,n)
enddo
close(17)
!==============================================================================!
!                          Generate Alpha Scaling
!==============================================================================!
open(unit=18,file='alphas.dat')
do i=1,NG
    alpha(i)=alpha0/(c_LJ*(P(x(:,i))*NG)**(-1./d))**2
    write(18,*) alpha(i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Generated Gaussian Widths'
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(alpha(i)*alpha(j))**(d/4.)*(alpha(i)+alpha(j))**(-0.5*d)&
             *exp(-0.5*aij*r2)
       Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
open(unit=19,file='overlap_eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'RCN = ', eigenvalues(1)/eigenvalues(NG)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
do l=1,Nsobol
    call sobol_stdnormal(d,skip,z(:,l))
enddo
write(*,*) 'Test 5; Successfully Generated Integration Sequence'
!==============================================================================!
!                              Theory Eigenvalues Value
!==============================================================================!
open(21,File=theory_in)
do i=1,NG
    read(21,*) theory(i)
enddo
close(21)
!==============================================================================!
!                             Evaluate Potential
!==============================================================================!
Vmat=0d0
open(unit=22,file='eigenval_conv.dat')
open(unit=23,file='abs_error_conv.dat')
do counter=1,Nsobol/data_freq
   do i=1,NG
      do j=i,NG
         x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
         do l=(counter-1)*data_freq+1,counter*data_freq
            Vmat(i,j)=Vmat(i,j)+Potential(x_ij(:)+z(:,l)/sqrt(alpha(i)+alpha(j)))
         enddo
      enddo
   enddo
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
   write(*,*) 'Generalized Eigenvalue Iteration', counter*data_freq
   do i=1,NG
      do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(alpha(i)*alpha(j))**(d/4.)*(alpha(i)+alpha(j))**(-0.5*d)&
             *exp(-0.5*aij*r2)
         Smat(j,i)=Smat(i,j)
!==============================================================================!
!                          Kinetic Energy Matrix
!==============================================================================!
         Hmat(i,j)=0.5*aij*(d-aij*r2)
!==============================================================================!
!                      Hamiltonian = Kinetic + Potential
!==============================================================================!
         Hmat(i,j)=(Hmat(i,j)+Vmat(i,j)/(counter*data_freq))*Smat(i,j)
         Hmat(j,i)=Hmat(i,j)
      enddo
   enddo
   itype=1
   eigenvalues=0d0
   call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
   write(22,*) alpha0, eigenvalues(:)
   write(23,*) alpha0, abs(theory(:)-eigenvalues(:))
   write(*,*) 'info ==> ', info
enddo
close(22)
close(23)
open(unit=24,file='rel_error.dat')
open(unit=25,file='alpha_rel.dat')
do i=1,NG
    write(24,*) i, (eigenvalues(i)-theory(i))/theory(i)
    write(25,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(23)
close(25)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'particle dsionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'Nsobol==>', Nsobol
write(90,*) 'Data Frequency==>', Data_Freq
write(90,*) 'alpha 0==>', alpha0
write(90,*) 'integral_P==>', integral_P
write(90,*) 'Ecut==>', E_cut
write(90,*) 'c_LJ ==>', c_LJ
close(90)
write(*,*) 'Hello Universe!'
end program main
