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
!       Modified:
!   2 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use spectra_mod
!==============================================================================!
!alpha_method       ==>Alpha generation method (constant, nearest, distribution)
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!GH_order           ==>Order for the Gauss-Hermite Quadriture
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!eigenvalues(NG)    ==>Eigenvalues of the matrix
!Smat(NG,NG)        ==>Overlap Matrix
!Hmat(NG,NG)        ==>Hamiltonian Matrix
!==============================================================================!
implicit none
character(len=20)::alpha_method,potential
integer::d,NG,GH_order
double precision::alpha0,V_i,E_cut,c_LJ,integral_P
double precision,allocatable,dimension(:)::alpha,eigenvalues
double precision,allocatable,dimension(:,:)::x,Smat,Hmat
integer::Lwork                                                    !LLAPACK dsygv
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
!information on the potential is only needed for alpha generation using the
!distribtion itself: (potential,Ecut,c_LJ,integral_P)
!==============================================================================!
read(*,*) d
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
read(*,*) alpha_method
read(*,*) potential
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) integral_P
!==============================================================================!
allocate(x(d,NG),alpha(NG))
allocate(eigenvalues(NG),Smat(NG,NG),Hmat(NG,NG))
Lwork=max(1,3*NG-1)                                          !LLAPACK Suggestion
!==============================================================================!
call read_grid(d,NG,x)
call generate_alphas(alpha_method,potential,d,NG,alpha0,alpha,x,V_i,E_cut,&
  c_LJ,integral_P)
call overlap_matrix(d,NG,x,alpha,Smat)
call overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
call get_hamiltonian(potential,d,NG,x,alpha,Smat,GH_order,Hmat)
call hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
call write_out(alpha_method,potential,d,NG,alpha0,GH_order)
!==============================================================================!
end program main
