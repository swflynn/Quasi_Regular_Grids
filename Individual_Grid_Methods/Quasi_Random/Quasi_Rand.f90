!=================20==================40==================60==================80
!                     Quasi-Random Grid Generation Main                       !
!==============================================================================!
!Quasi-Random Grid for a d-dimensional system using a cutoff or rejection.
!Requires the potentials_mod.f90 and sobol.f90 (see Github Main)
!==============================================================================!
!       Modified:
!   26 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program qrand_grid
use qrand_grid_mod
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D
double precision::E_cut,integral_P,V_i
double precision,allocatable,dimension(:)::xmin,xmax
double precision,allocatable,dimension(:,:)::x
character(len=20)::potential
logical::reject
!==============================================================================!
!potential          ==>Potential name
!reject             ==>Rejection logical (pseudo-random/quasi-random grids only)
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
read(*,*) Npoints
read(*,*) d
read(*,*) reject
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) potential
!==============================================================================!
allocate(x(d,Npoints),xmin(d),xmax(d))
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
call Q_Rand(potential,reject,Npoints,d,x,V_i,E_cut,xmin,xmax,&
  integral_P)
call write_out(potential,reject,Npoints,d,E_cut,xmin,xmax,N_1D,N_MMC_box,&
  integral_P)
end program qrand_grid
