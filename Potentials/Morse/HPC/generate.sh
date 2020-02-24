#=============================================================================80
#                              Morse HPC Scripts
#=============================================================================80
#Bash script for generating calculations as a function of alpha0
#Makes a new directory for each calcualtion
#==============================================================================#
#       Modified:
#   19 May 2019
#       Author:
#   Shane Flynn
#==============================================================================#
#grid_in        ==>Filename Containing Gridpoints
#theory_in      ==>Filename Containing Analytic Eigenvalues
#d              ==>dimensionality
#NG             ==>Number of Gaussian Basis Functions (gridpoints)
#GH_order       ==>Number of Points for evaluating the potential (Gauss-Hermite)
#alpha0         ==>scaling parameters, vary with start/inc counters
#n              ==>place each calculation in a new directory (counter)
#max            ==>Number of calculations to generate
#start          ==>starting value for alpha0 parameter
#inc            ==>increment to increase alpha0 parameter
#               Directory for Calculations as a function of alpha0
#==============================================================================#
mkdir alpha
cd alpha
cp ../grid.dat .
cp ../all_theory.dat .
cp ../submit.sh .
cp ../runner.sh .
#==============================================================================#
#                         Set Calculation Parameters
#==============================================================================#
d=2;
NG=300;
GH_order=10;
grid_in=grid.dat;
theory_in=all_theory.dat;
n=1;
max=100;
start=0.0;
inc=0.1;
#==============================================================================#
#                       New directory n for each calculation
#==============================================================================#
while [ "$n" -le "$max" ]; do
  mkdir "$n"
  cp ../submit.sh "$n"
  cp ../grid.dat "$n"
  cp ../all_theory.dat "$n"
  cd "$n"
#single arrow makes a file, double arrow appends to an existing file
  echo "$d" > input
  echo "$NG" >> input
  echo "$GH_order" >> input
  echo "$grid_in" >> input
  echo "$theory_in" >> input
  echo "$start + $n*$inc" | bc  >> input
  cd ..
  n=`expr "$n" + 1`;
done
