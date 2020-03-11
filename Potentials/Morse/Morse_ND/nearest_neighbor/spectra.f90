!=================20==================40==================60==================80
!                           nD-Morse EigenSpectra
!==============================================================================!
!Compute the EigenSpectra for the nD Morse Potential
!generate alpha (Gaussian Widths) based on nearest neghbor
!This code needs the gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture)
!This code uses LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   8 March 2020
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
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
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
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
!==============================================================================!
call read_grid(d,NG,x)
call gaussian_widths(d,NG,alpha0,alpha,x)
call overlap_matrix(d,NG,alpha,Smat,x)
call diagonalize_overlap(NG,Smat,eigenvalues,lwork)
call get_hamiltonian(d,NG,GH_order,x,alpha,Smat,Hmat)
call hamiltonian_eigenvalues(NG,Hmat,Smat,eigenvalues,work,Lwork)
call write_out(d,NG,alpha0,GH_order)
!==============================================================================!
end program main
