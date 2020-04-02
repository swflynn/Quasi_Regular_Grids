!=============================================================================80
!                           Gaussian Sobol Grid
!==============================================================================!
!       Discussion:
!Generate a d-dimensional Multivariate Gaussian Distribution of points using 
!the a sobol sequence and the Beasley-Springer-Moro transformation. 
!Requires the sobol.f90 fortran module
!==============================================================================!
!       Modified:                                                              !
!   16 April 2019                                                              !
!       Author:                                                                !
!   Shane Flynn                                                                !
!==============================================================================!
program main
!==============================================================================!
!               Discussion:
!==============================================================================!
!d              ==> Dimensionality of the ith point (x^i=x^i_1,x^i_2,..,x^i_d)
!Npoints        ==> Total Number of Points
!x              ==>(d,Npoints) points coordinates
!==============================================================================!
implicit none
integer::d,Npoints
integer*8::ii,skip                           !=need *8 for the sobol generator=!
double precision,allocatable,dimension(:,:)::x
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) Npoints
read(*,*) d
skip=Npoints
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints))
!==============================================================================!
!                              Generate Grid 
!==============================================================================!
open(unit=17,file='grid.dat')
do ii=1,Npoints
    call sobol_stdnormal(d,skip,x(:,ii))
    write(17,*) x(:,ii)
enddo
close(17)
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
open(99,file='simulation.dat')
write(99,*) 'particle dimensionality ==> ', d
write(99,*) 'Npoints ==> ', Npoints
close(99)
write(*,*) 'Hello Universe!'
end program main
