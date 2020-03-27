!=================20==================40==================60==================80
!                     Quasi-Regular Grid Generation Main                       !
!==============================================================================!
!Generate a Quasi-Regular Grid for any d-dimensional system of interest.
!Requires the potentials_mod.f90 (see Main)
!==============================================================================!
!       Modified:
!   26 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use qrg_grid_mod
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,N_MMC_grid,MMC_freq
double precision::E_cut,c_LJ,integral_P,V_i
double precision,allocatable,dimension(:)::xmin,xmax
double precision,allocatable,dimension(:,:)::x,U
character(len=20)::potential
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
read(*,*) Npoints
read(*,*) d
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) c_LJ
read(*,*) N_MMC_grid
read(*,*) MMC_freq
read(*,*) potential
!==============================================================================!
allocate(x(d,Npoints),U(Npoints,Npoints),xmin(d),xmax(d))
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
call initial_distribution(potential,Npoints,d,x,V_i,E_cut,xmin,xmax)
call initial_pair_energy(potential,Npoints,d,x,U,V_i,E_cut,c_LJ,integral_P)
call Q_Reg(potential,Npoints,d,x,U,V_i,E_cut,c_LJ,N_MMC_grid,MMC_freq,&
  integral_P)
call write_out(potential,Npoints,d,E_cut,xmin,xmax,N_1D,N_MMC_box,c_LJ,&
  N_MMC_grid,MMC_freq,integral_P)
end program main
