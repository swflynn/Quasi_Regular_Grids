!=================20==================40==================60==================80
!                         nD Vibrational EigenSpectra
!==============================================================================!
!Vibrational EigenSpectra calculations for a given potential and grid.
!All equations written for general nD case
!Basis Functions:
!             phi(r):=(2*alpha_i/pi)^(d/4)*exp[-alpha_i(r-r_i)^2]
!Generates alpha (inverse gaussian width i.e. small alpha=broad gaussian)
!constant (uniform grids), based on the distribtion (pseudo/quasi random)
!or using nearest neighbor (quasi-regular grids)
!Needs gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture) for potential eval.
!Needs LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!We have also implemented 2D basis functions proposed by Garaschuk and light
!S. Garashchuk and J. C. Light, “Quasirandom distributed Gaussian bases for
!bound problems,” J. Chem. Phys. 114, 3929–3939 (2001).
!These basis functions do not scale with d, however, given a few orders of
!magnitude higher accuracy for the 2D morse case.
!==============================================================================!
!       Modified:
!   2 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module spectra_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine read_grid(d,NG,x)
!==============================================================================!
!Reads in grid points, each line contains coordinates for a d-dimesional point
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!==============================================================================!
implicit none
integer::d,NG,i
double precision::x(d,NG)
open(17,File='grid.dat')
do i=1,NG
  read(17,*) x(:,i)
enddo
close(17)
end subroutine read_grid
!==============================================================================!
subroutine generate_alphas(alpha_method,potential,d,NG,alpha0,alpha,x,V_i,E_cut&
  ,c_LJ,integral_P)
!==============================================================================!
!Determines which alpha generation method to call
!Symmetric gaussians (alpha_i is the same for each dimension of the gaussian)
!QRG==>nearest neighbor, uniform grid==>constant, random grid==>distribution
!==============================================================================!
!alpha_method       ==>Alpha generation method (constant, nearest, distribution)
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::d,NG
double precision::alpha0,alpha(NG),x(d,NG),V_i,E_cut,c_LJ,integral_P
character(len=20)::alpha_method,potential
!==============================================================================!
if(alpha_method=='constant') then
  call const_alpha(NG,alpha0,alpha)
elseif(alpha_method=='nearest') then
  call near_alpha(d,NG,alpha0,alpha,x)
elseif(alpha_method=='distribution') then
  call dist_alpha(potential,d,NG,alpha0,alpha,x,V_i,E_cut,c_LJ,integral_P)
else
  stop 'Cannot Identify Alpha Generation Method, Check "alpha_type" Subroutine'
endif
end subroutine generate_alphas
!==============================================================================!
subroutine const_alpha(NG,alpha0,alpha)
!==============================================================================!
!Sets the alpha parameter to a constant
!This should be used for uniform (square) grids
!==============================================================================!
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!==============================================================================!
implicit none
integer::NG
double precision::alpha0,alpha(NG)
!==============================================================================!
alpha=alpha0
call write_alphas(NG,alpha)
end subroutine const_alpha
!==============================================================================!
subroutine near_alpha(d,NG,alpha0,alpha,x)
!==============================================================================!
!Uses the nearest neighbor to sets the alpha parameter
!This should be used for Quasi-Regular Grids
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!==============================================================================!
implicit none
integer::d,NG,i,j
double precision::alpha0,alpha(NG),x(d,NG),r2
!==============================================================================!
do i=1,NG
  alpha(i)=1d20                          !large distance for initial placeholder
  do j=1,NG
    if(j.ne.i) then
      r2=sum((x(:,i)-x(:,j))**2)                    !distance between gridpoints
      if(r2<alpha(i)) alpha(i)=r2
    endif
  enddo
  alpha(i)=alpha0/alpha(i)
enddo
call write_alphas(NG,alpha)
end subroutine near_alpha
!==============================================================================!
subroutine dist_alpha(potential,d,NG,alpha0,alpha,x,V_i,E_cut,c_LJ,integral_P)
!==============================================================================!
!Uses the Target Distribution Function to generate the alpha parameter
!This should be used for pseudo-random and quasi-random grids
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,NG,i
double precision::x(d,NG),E_cut,c_LJ,alpha0,alpha(NG),integral_P,V_i
!==============================================================================!
if(d==2)then
  do i=1,NG
    alpha(i)=alpha0/(c_LJ*(P_i(potential,d,x(:,i),V_i,E_cut,integral_P)*NG)&
    **(-1./d))**2
  enddo
  call write_alphas(NG,alpha)
else
  do i=1,NG
    alpha(i)=.5*(alpha0/(c_LJ*(P_i(potential,d,x(:,i),V_i,E_cut,integral_P)*NG)&
    **(-1./d))**2)
  enddo
  call write_alphas(NG,alpha)
endif
end subroutine dist_alpha
!==============================================================================!
subroutine write_alphas(NG,alpha)
!==============================================================================!
!Write alphas to file
!==============================================================================!
!NG                 ==>Number of grid points
!alpha(NG)          ==>Inverse Gaussian Widths
!==============================================================================!
integer::NG,i
double precision::alpha(NG)
open(unit=18,file='alphas.dat')
do i=1,NG
  write(18,*) alpha(i)
enddo
close(18)
end subroutine write_alphas
!==============================================================================!
subroutine potentials(potential,d,x_i,V_i)
!==============================================================================!
!Determines which Potential Energy Function to call
!See the potentials_mod.f90 for available potentials, and update this
!subroutine accordingly.
!All equations derived for a d-dimensional system, any potential is valid.
!==============================================================================!
!potential        ==>potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i              ==>Potential Energy evaluation V(x_i)
!==============================================================================!
use potential_mod
implicit none
integer::d
double precision::x_i(d),V_i
character(len=10)::potential
!==============================================================================!
if(potential=='morse_pot') then
  call morse(d,x_i,v_i)
elseif(potential=='henon_pot') then
  call henon(d,x_i,v_i)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine potentials
!==============================================================================!
function P_i(potential,d,x_i,V_i,E_cut,integral_P)
!==============================================================================!
!Target Distribution Function, !defined according to a semi-classical argument:
!B. Poirier, “Algebraically self-consistent quasiclassical approximation on
!phase space,” Found. Phys. 30, 1191–1226 (2000).
!This subroutine computes P(x_i)
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)             ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!integral_P         ==>Normalization constant for the distribtion P(x)
!P_i                ==>evaluate P(x_i)
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),E_cut,P_i,integral_P,V_i
character(len=10)::potential
!==============================================================================!
call potentials(potential,d,x_i,V_i)
if(V_i<E_cut) P_i=(E_cut-V_i)**(d/2.)/integral_P
if(V_i>=E_cut) P_i=1d-20                      !Define distribution=0 beyond Ecut
end function P_i
!==============================================================================!
subroutine overlap_elements(d,x_i,x_j,alpha_i,alpha_j,S_ij)
!==============================================================================!
!Compute i-jth element of the overlap matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameters
!Sij              ==>i-j element of the overlap matrix
!==============================================================================!
integer::d
double precision::alpha_i,alpha_j,S_ij,x_i(d),x_j(d),aij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
if(d==2)then
  S_ij=(sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d)*exp(-0.5*aij*r2)
else
  S_ij=(2*sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d)*exp(-aij*r2)
endif
end subroutine overlap_elements
!==============================================================================!
subroutine overlap_matrix(d,NG,x,alpha,Smat)
!==============================================================================!
!Compute the overlap matrix (symmetric, positive-definite)
!Used to check the stability of the basis; see overlap_eigenvalues subroutine
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!==============================================================================!
integer::d,NG,i,j
double precision::alpha(NG),Smat(NG,NG),x(d,NG)
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    Smat(j,i)=Smat(i,j)
  enddo
enddo
end subroutine overlap_matrix
!==============================================================================!
subroutine overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
!==============================================================================!
!Compute the Eigenvalues for the overlap matrix to check numerical stability
!Must be positive definite; check the Recriprocal Condition Number
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!eigenvalues(NG)  ==>Eigenvalues of the matrix
!==============================================================================!
implicit none
integer::i,NG,info,Lwork
double precision::eigenvalues(NG),Smat(NG,NG),work(max(1,Lwork))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)              !LLAPACK
write(*,*) 'Info (Overlap Matrix) ==> ', info
write(*,*) 'RCN ==> ', eigenvalues(1)/eigenvalues(NG)
open(unit=19,file='overlap.dat')
do i=1,NG
  write(19,*) eigenvalues(i)
enddo
close(19)
end subroutine overlap_eigenvalues
!==============================================================================!
subroutine kinetic_elements(d,x_i,x_j,alpha_i,alpha_j,T_ij)
!==============================================================================!
!Compute i-jth element of the kinetic matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameter
!T_ij             ==>i-j element of the kinetic matrix
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),x_j(d),alpha_i,alpha_j,aij,r2,T_ij
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
if(d==2)then
  T_ij=0.5*aij*(d-aij*r2)
else
  T_ij=aij*(d-2.*aij*r2)
endif
end subroutine kinetic_elements
!==============================================================================!
subroutine get_hamiltonian(potential,d,NG,x,alpha,Smat,GH_order,Hmat)
!==============================================================================!
!Compute the Hamiltonian Matrix
!Use Gauss-Hermite quadrature for the potential
!==============================================================================!
!potential        ==>Potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!GH_order         ==>Gauss-Hermite quadrature order
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!z(GH_order)      ==>Quadrature points
!w(GH_order)      ==>Quadrature points weights
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,NG,GH_order,i,j,k,ll
double precision,parameter::pi=4.*atan(1d0)
double precision::x(d,NG),alpha(NG),Smat(NG,NG),Hmat(NG,NG),z(GH_order)
double precision::w(GH_order),x_ij(d),V_ij,l(d),rr(d),r2
if(d==2)then
  call cgqf(GH_order,6,0d0,0d0,0d0,0.5d0,z,w)                           !LLAPACK
  w=w/sqrt(2.*pi)
else
  call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                             !LLAPACK
  w=w/sqrt(pi)
endif
!==============================================================================!
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    call kinetic_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Hmat(i,j))
!==============================================================================!
!                    Compute Potential Energy with Quadriture
!==============================================================================!
    x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
    V_ij=0d0
    l(:)=1
    do ll=1,GH_order**d
      do k=1,d
        rr(k)=z(l(k))
      enddo
      rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
      call potentials(potential,d,rr,r2)
      do k=1,d
        r2=r2*w(l(k))
      enddo
      V_ij=V_ij+r2
      do k=1,d
        l(k)=mod(l(k),float(GH_order))+1
        if(l(k).ne.1) exit
      enddo
    enddo
    Hmat(i,j)=(Hmat(i,j)+V_ij)*Smat(i,j)
    Hmat(j,i)=Hmat(i,j)
  enddo
enddo
end subroutine get_hamiltonian
!==============================================================================!
subroutine hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
!==============================================================================!
!Compute the Eigenvalues for the Hamiltonian matrix
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!eigenvalues(NG)  ==>Matrix eigenvalues
!==============================================================================!
implicit none
integer::NG,info,Lwork,itype,i
double precision::Smat(NG,NG),Hmat(NG,NG),eigenvalues(NG),work(max(1,lwork))
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
open(unit=20,file='eigenvalues.dat')
do i=1,NG
  write(20,*) eigenvalues(i)
enddo
close(20)
end subroutine hamiltonian_eigenvalues
!==============================================================================!
subroutine write_out(alpha_method,potential,d,NG,alpha0,GH_order)
!==============================================================================!
!Simulation details
!==============================================================================!
!alpha_method       ==>Alpha generation method (constant, nearest, distribution)
!potential          ==>Potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!alpha0           ==>Flat scaling parameter for optimizing widths
!GH_order         ==>Gauss-Hermite quadrature order
!==============================================================================!
implicit none
integer::d,NG,GH_order
double precision::alpha0
character(len=20)::alpha_method,potential
open(unit=99,file='out')
write(99,*) 'alpha method ==> ', alpha_method
write(99,*) 'potential ==> ', potential
write(99,*) 'dimensionality ==> ', d
write(99,*) 'Number of  Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
close(99)
end subroutine write_out
!==============================================================================!
end module spectra_mod
