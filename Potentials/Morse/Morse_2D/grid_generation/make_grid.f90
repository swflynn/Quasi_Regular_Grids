!=================20==================40==================60==================80
!                   2D Morse Quasi-Reglar Grid Generation
!==============================================================================!
!Generate gridpoints distributed via the 2D Morse Oscillator
!This Quasi-Regular gird is optimized using a quasi-Lennard Jones Potential
!Force minimization, accept any trial moves that reduces the system's energy
!See MMC grid for metropolis optimization (acceptance criteria)
!==============================================================================!
!Start developing 2D module want to make a module for each code then remove the
!various seperate directories
!the 2D case can use more efficient basis functions, these are not generalizable
!==============================================================================!
!       Modified:
!   23 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use make_grid_mod
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,N_MMC_grid,MMC_freq
double precision::E_cut,c_LJ,D_morse,integral_P
double precision,allocatable,dimension(:)::omega
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
read(*,*) Npoints
read(*,*) d
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) Npoints
read(*,*) c_LJ
allocate(omega(d))
read(*,*) omega
read(*,*) D_morse
read(*,*) N_MMC_grid
read(*,*) MMC_freq
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(xmin(d),xmax(d),x(d,Npoints),U(Npoints,Npoints))
!==============================================================================!
call box_size_P(d,N_MMC_box,E_cut,omega,D_morse,integral_P)
call compute_integral_P(d,N_1D,E_cut,omega,D_morse,integral_P)
call initial_distribution(Npoints,d,x,E_cut,omega,D_morse)
call initial_pairwise_energy(Npoints,d,x,U,E_cut,c_LJ,omega,D_morse,integral_P)
call generate_grid(Npoints,d,x,U,E_cut,c_LJ,omega,D_morse,N_MMC_Grid,MMC_freq,integral_P)
call write_out(Npoints,d,E_cut,N_1D,c_LJ,omega,D_morse,N_MMC_grid,MMC_freq,integral_P)
end program main
