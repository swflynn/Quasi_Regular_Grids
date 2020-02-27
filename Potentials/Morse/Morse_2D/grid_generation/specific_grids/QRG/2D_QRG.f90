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
use Morse2D_QRG_mod
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D,N_MMC_grid,MMC_freq
double precision::E_cut,c_LJ
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) Npoints
read(*,*) c_LJ
read(*,*) N_MMC_grid
read(*,*) MMC_freq
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(xmin(d),xmax(d),x(d,Npoints),U(Npoints,Npoints))
!==============================================================================!
call box_size_P(N_MMC_box,E_cut)
call compute_integral_P(N_1D,E_cut)
call initial_distribution(x,E_cut)
call initial_pairwise_energy(x,U,E_cut,c_LJ)
call generate_grid(x,U,E_cut,c_LJ,N_MMC_Grid,MMC_freq)
call write_out(E_cut,N_1D,c_LJ,N_MMC_grid,MMC_freq)
end program main
