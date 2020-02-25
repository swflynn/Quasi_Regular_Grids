!=================20==================40==================60==================80
!                   2D Morse Quasi-Reglar Grid Generation
!==============================================================================!
!This grid is generated using a sobol sequence and the rejection method or Ecut
!==============================================================================!
!       Modified:
!   23 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use Morse2D_sobol_mod
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D
double precision::E_cut
double precision,allocatable,dimension(:,:)::x
logical::rej
!==============================================================================!
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) Npoints
read(*,*) rej
!==============================================================================!
!                                 Allocations
!==============================================================================!
allocate(xmin(d),xmax(d),x(d,Npoints))
!==============================================================================!
call box_size_P(N_MMC_box,E_cut)
call compute_integral_P(N_1D,E_cut)
call generate_grid(rej,E_cut,x)
call write_out(E_cut,N_1D,rej,N_MMC_box)
end program main
