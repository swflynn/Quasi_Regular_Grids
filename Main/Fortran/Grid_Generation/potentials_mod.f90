!=================20==================40==================60==================80
!                              Potentials Module
!==============================================================================!
!Any user defined potential is valid assuming a single point x_i = (x_1,..x_d)
!results in a single energy evaluation [V_i(x_i) = #]
!After defining a potential add the appropriate call statement to the
!"potentials" subroutine in grids_mod.f90
!All equations in grid_generation and vibrational_spectra are derived 
!and implemented as d-dimensional.
!==============================================================================!
!       Modified:
!   24 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module potential_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine morse(d,x_i,V_i)
!==============================================================================!
!Potential Energy: Sum of 1D Morse potentials
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!D_morse          ==>Parameter for Morse Potential
!omega(d)         ==>Parameter for Morse Potential
!V_i              ==>Potential Energy evaluation V(x_i)
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),V_i
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
!double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
V_i=D_morse*sum((exp(-omega(:)*x_i(:))-1.)**2 )
end subroutine morse
!==============================================================================!
subroutine henon(d,x_i,V_i)
!==============================================================================!
!Potential Energy: 2d Henon-Heiles
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i              ==>Potential Energy evaluation V(x_i)
!lambda           ==>Parameter for Henon-Heiles Potential
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),V_i
double precision,parameter::lambda=sqrt(0.0125)
V_i=0.5*(x_i(1)**2+x_i(2)**2)+lambda*(x_i(1)**2*x_i(2)-x_i(2)**3/3.)
end subroutine henon
!==============================================================================!
end module potential_mod
