!=============================================================================80
!                 2D Std. Norm. Quais-Lennard-Jones Grid
!==============================================================================!
!       Discussion:
!Generate gridpoints (Particles) distributed via the 2D Std. Norm. Dist.
!This grid is 'optimized' using a quasi-Lennard Jones Potential
!Gridpoints have a quasi-regular distribution locally
!Minimization implemented analagous to traditional MMC (beta,epsilon). 
!==============================================================================!
!       Modified:
!   16 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module qlj_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==>Dimensionality of the ith point (x^i=x^i_1,x^i_2,..,x^i_d)
!epsilon_lj     ==>epsilon parameter for LJ 
!c_LJ           ==>parameter for LJ
!==============================================================================!
integer::d,Npoints
double precision::epsilon_LJ,c_LJ
double precision,parameter::pi=acos(-1d0)
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the range 1-Npoints for MMC Algorithm
!==============================================================================!
!Nmin           ==>minimum index value (1)
!Nmax           ==>maximum index value (Npoints)
!random_integer ==>integer returned
!a              ==>Fortran intrinsic random number [0,1]
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
!Target Distribution Function (Hard-Coded)
!==============================================================================!
!P_x            ==>evaluates P(x)
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
!quasi-LJ energy between 2 grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate LJ
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_x(x1)*Npoints)**(-1./d)    
sigma2=c_LJ*(P_x(x2)*Npoints)**(-1./d)    
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_LJ_NRG=4.*epsilon_LJ*(a**2-a+b**2-b)
end function Pair_LJ_NRG
!==============================================================================!
end module qlj_mod
!==============================================================================!
!==============================================================================!
program main
use qlj_mod
!==============================================================================!
!               Discussion:
!==============================================================================!
!Npoints     ==> Number of Particles
!N_iter         ==>Number of MMC Iterations
!beta           ==>Inverse Temperature
!x              ==>(dimen) all atom coordinates
!dimen          ==> total system dimensionality (configuration space)
!d              ==> i-th particle dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
!V              ==>(Npoints) All energies due to distribution V(x)=-ln[P(x)]
!U              ==>(Npoints,Npoints) All i,j pair-wise energies (LJ-12-6)
!U_move         ==>(Npoints)  LJ-Energy associated with trial movement
!V_move         ==> V(x) associated with trial move V(x) = -ln[P(x)]
!mv_cutoff      ==> maximum displacement parameter for trial displacement
!x0             ==>(d) store previous coordinate before trial move
!acc_coef       ==> scaling coefficient for making acceptance rate ~50%
!max_move       ==> cutoff for local movement
!Delta_E        ==> change in energy due to potential
!t1             ==> random number for accepting higher energy movement
!s              ==>(d) random number to move coordinate by
!freq           ==> Interval to update mv_cutoff size
!==============================================================================!
implicit none
integer::N_iter,accept,counter,freq,i,j,k,n
integer*8::ii,skip                          !need *8 for sobol generator!
double precision::t1,Delta_E,mv_cutoff,beta
double precision,allocatable,dimension(:)::x0,s,U_move
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) Npoints
read(*,*) d
read(*,*) N_iter
read(*,*) beta
read(*,*) epsilon_LJ
read(*,*) c_LJ
read(*,*) mv_cutoff
read(*,*) freq
skip=Npoints
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints),x0(d),s(d),U(Npoints,Npoints),U_move(Npoints))
!==============================================================================!
!                       Generate Initial Distribution
!               Initally accept any point where Potential<Ecut
!==============================================================================!
i=1
do while(i.le.Npoints)   
    call random_number(s)
    x(:,i)=s(:)
    i=i+1            
enddo
!==============================================================================!
!                          Write Initial Coordinates
!==============================================================================!
open(unit=17,file='coor_ini.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
close(17)
write(*,*) 'Test 2; Successfully generated initial grid' 
!==============================================================================!
!                             Compute U[x_ij] 
!compute pairwise energies for the entire set of initial particle 
!==============================================================================!
do i=2,Npoints
    do j=1,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
enddo
!==============================================================================!
!                               Compute V[x_i] 
!compute pairwise energies for the entire set of initial particle 
!==============================================================================!
do i=1,Npoints
    U(i,i)=-log(P_x(x(:,i)))
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
    U_move(k)=-log(P_x(x0))
    Delta_E=U(K,k)-U_move(k) 
    do j=1,Npoints
        if(j.ne.k) then
           U_move(j)=Pair_LJ_NRG(x(:,j),x0)
           Delta_E=Delta_E+U(j,k)-U_move(j)
        endif
    enddo
!==============================================================================!
!               Test to see if you should accept higher energy move
!criteria depends on how you define inequality
!if beta > 0 always, if delta e > 0 than always larger than 1, always accept
!==============================================================================!
    call random_number(t1) 
    if(exp(beta*Delta_E).ge.t1)then
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
open(90,file='simulation.dat')
write(90,*) 'Npoints ==> ', Npoints
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Number of MMC Iterations ==> ', N_iter
write(90,*) 'epsilon LJ ==> ', epsilon_LJ 
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'Final Move Cutoff ==> ', mv_cutoff
close(90)
write(*,*) 'Hello Universe!'
end program main
