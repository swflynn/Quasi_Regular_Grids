!=============================================================================80
!                      2D Gaussian Quasi Regular Grid
!!!!!Need to finish this!!!!
!==============================================================================!
!       Discussion:
!Generate a Quasi-Regular (QR) 2D Standard Normal Grid using a pair-wise 
!interaction between grid points (modified Lennard-Jones Potential).
!MMC is used to update the grid points.
!This code is for demonstration purposes, start with a 2D-Standard-Normal 
!distribution generated with sobol points for convenience (requires sobol.f90)
!==============================================================================!
!       Modified:
!   16 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module QR_std_norm
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==>Dimensionality of the grid points
!Npoints        ==>Number of Grid-Points
!c_LJ           ==>Lennard-Jones Parameter
!==============================================================================!
integer::d,Npoints
double precision::c_LJ
double precision,parameter::pi=acos(-1d0)
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the Nmin to Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!random_integer ==>integer returned
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function P_x(x)
!==============================================================================!
!       Discussion:
!Target distribution function
!==============================================================================!
!P_x            ==>Evaluate P_x
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none 
double precision::x(d),P_x
P_x=(2.*pi)**(-d/2.)*exp(-0.5*sum(x(:)**2))
end function P_x
!==============================================================================!
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!       Discussion:
!Modifies form of the LJ-12-6 Pair-wise interaction between particles
!Generates a global distribution according to P_x, while maintaining local 
!uniformity.
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!Pair_LJ_NRG    ==>Energy due to ith,jth particle interaction
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_x(x1)*Npoints)**(-1./d)    
sigma2=c_LJ*(P_x(x2)*Npoints)**(-1./d)    
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_LJ_NRG=(a**2-a+b**2-b)
end function Pair_LJ_NRG
!==============================================================================!
end module QR_std_norm
!==============================================================================!
!==============================================================================!
!==============================================================================!
program main
use QR_std_norm
!==============================================================================!
!       Discussion:
!==============================================================================!
!x              ==>(d,Npoins) Coordinates of the Grid Points
!N_iter         ==>Number of MMC Iterations for Optimizing the Grid 
!mv_cutoff      ==>Cutoff Distance for MMC Trial Move
!freq           ==>Frequency at which the acceptance ratio is adjusted
!x0             ==>(d) store initial grid points coordinate before trial move
!U              ==>(Npoints,Npoints) Store all i,j pair-wise energies
!U_move         ==>(Npoints) All effected pair-wise energies due to trial move 


!acc_coef       ==> scaling coefficient for making acceptance rate ~50%
!Delta_E        ==> change in energy due to potential
!t1             ==> random number for accepting higher energy movement
!accept         ==> number of accepted trial moves, for acceptance~50%  
!counter        ==> total number of moves, for acceptance~50%
!moment_i/f     ==> the second moment of the initial and final distributions
!==============================================================================!
implicit none
integer::i,j,k,n,N_iter,accept,counter,freq
integer*8::ii,skip                          !need *8 for sobol generator!
double precision::t1,Delta_E,mv_cutoff,t_i,t_f,moment_i,moment_f,error
double precision,allocatable,dimension(:)::x0,s,U_move
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(t_i)
read(*,*) Npoints
read(*,*) d
read(*,*) N_iter
read(*,*) c_LJ
read(*,*) mv_cutoff
read(*,*) freq
skip=Npoints                            !skip is required for the sobol genertor
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints),x0(d),s(d),U(Npoints,Npoints),U_move(Npoints))
!==============================================================================!
!                           Generate Initial Guess
!Any initial guess is sufficient, for convenience use a 2D sobol standard normal
!This requires sobol.f90 code, and sobol_stdnormal.f90
!==============================================================================!
open(unit=16,file='coor_ini.dat')
do ii=1,Npoints
    call sobol_stdnormal(d,skip,x(:,ii))
    write(16,*) x(:,ii)
enddo
close(16)
!==============================================================================!
!               Compute second moment of sobol distribution   <|x|2> 
!moments are a simple metric to assess if the distribution is correct
!==============================================================================!
moment_i=0d0
do i=1,Npoints
    moment_i=moment_i+sum(x(:,i)**2)
enddo
moment_i=moment_i/Npoints
write(*,*) 'moment ', moment_i
!==============================================================================!
!                               Compute U[x_ij] 
!compute pairwise energies for the entire set of initial particle 
!==============================================================================!
do i=2,Npoints
    do j=1,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
enddo
!==============================================================================!
!                               Begin MMC 
!==============================================================================!
accept=0
counter=0
do n=1,N_iter
!==============================================================================!
!                       Randomly Select Atom to Move
!==============================================================================!
    k=random_integer(1,Npoints)
!==============================================================================!
!                   Generate coordinates for Random move trial
!random numbers generated [0,1], make it [-1,1] ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    counter=counter+1
!==============================================================================!
!                           make trial move
!==============================================================================!
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Compute Energy Change due to Trial Move
!first compute Delta V = Vnew[x_i] - V[x_i]
!Then compute delta U and add together for total change in energy
!==============================================================================!
    Delta_E=0d0
    do j=1,Npoints
        if(j.ne.k) then
           U_move(j)=Pair_LJ_NRG(x(:,j),x0)
           Delta_E=Delta_E+U(j,k)-U_move(j)
        endif
    enddo
!==============================================================================!
!               Test to see if you should accept higher energy move
!==============================================================================!
    call random_number(t1) 
    if(exp(Delta_E).ge.t1)then
       U(:,k)=U_move(:)
       U(k,:)=U_move(:)
       accept=accept+1
       x(:,k)=x0(:)
    else
!==============================================================================!
!                   Otherwise Reject Configuration
!==============================================================================!
    endif
!==============================================================================!
!                           Update Cutoff Paramater
!acceptance rate ~50%, adjust random movement displacement length accordingly
!==============================================================================!
    if(mod(n,freq)==0)then
        write(*,*) 'iteration', n
        if(dble(accept)/counter<0.5)then 
            mv_cutoff=mv_cutoff*0.9
        else
            mv_cutoff=mv_cutoff*1.1
        endif
        accept=0
        counter=0
    endif
enddo
!==============================================================================!
!                   Write Final Configuration to file 
!may want to make this as a function of MMC Iteration in the future
!==============================================================================!
open(unit=17,file='qlj_grid.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(t_f)
open(90,file='simulation.dat')
write(90,*) 'Npoints ==> ', Npoints
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Number of MMC Iterations ==> ', N_iter
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'Final Move Cutoff ==> ', mv_cutoff
write(90,*) 'Moment Sobol Dist.', moment_i
write(90,*) 'Moment final Dist.', moment_f
write(90,*) 'Moment Error', error
write(90,*) 'Total Time ==> ', t_f-t_i
close(90)
write(*,*) 'Hello Universe!'
end program main
