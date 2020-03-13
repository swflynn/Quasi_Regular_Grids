!=================20==================40==================60==================80
!                           nD-Morse EigenSpectra
!==============================================================================!
!Compute the EigenSpectra for the nD Morse Potential
!generate alpha (Gaussian Widths) based on nearest neghbor
!This code needs the gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture)
!This code uses LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   13 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use spectra_mod
!==============================================================================!
implicit none
integer::d,NG,GH_order
double precision::alpha0
double precision,allocatable,dimension(:)::alpha,eigenvalues
double precision,allocatable,dimension(:,:)::x,Smat,Hmat
!==============================================================================!
!                            LLAPACK dsygv variables
!==============================================================================!
integer::itype,info,Lwork
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) d
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
!==============================================================================!
!                                 Allocations
!==============================================================================!
allocate(x(d,NG),alpha(NG),eigenvalues(NG),Smat(NG,NG),Hmat(NG,NG))
Lwork=max(1,3*NG-1)
!==============================================================================!
call read_grid(d,NG,x)
call gaussian_widths(d,NG,x,alpha0,alpha)
call overlap_matrix(d,NG,x,alpha,Smat)
call overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
call get_hamiltonian(d,NG,x,alpha,Smat,GH_order,Hmat)
call hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
call write_out(d,NG,alpha0,GH_order)
!==============================================================================!
end program main
