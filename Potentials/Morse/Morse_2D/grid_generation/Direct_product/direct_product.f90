!=================20==================40==================60==================80
!                   2D Morse Direct Product Grid Generation
!==============================================================================!
!       Modified:
!   25 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use Morse2D_direct_product_mod
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D
double precision::E_cut
double precision,allocatable,dimension(:,:)::x
!==============================================================================!
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
!==============================================================================!
!                                 Allocations
!==============================================================================!
allocate(xmin(d),xmax(d),x(d,Npoints))
!==============================================================================!
call box_size_P(N_MMC_box,E_cut)
call generate_grid(N_1D,E_cut)
call write_out(E_cut,N_1D,N_MMC_box)
end program main
