# Quasi-Regular Grids; Applications to the Vibrational Spectra Calculations
Fortran and Python implementations for generating a Quasi-Regular Grid (QRG)
using model Potential Energy Functions.
See our [manuscript](https://doi.org/10.1063/1.5134677) (JCP Communication) 
for details.

Given a general distribution of interest the Quasi-Regular algorithm generates a
set of points (the grid) that optimally samples the distribution.
We demonstrate the utility of these grids by constructing a Distributed Gaussian
Basis (DGB) and computing the associated Vibrational Eigenspectra.

## Main
Fortran and Python implementations for generating a grid and evaluating the
associated Vibrational Eigenspectra.
All equations are d-dimensional, therefore any potential of interest can be
defined by the user.

### Grid_Generation
Consistent with the manuscript the following grid generation methods are
available:
* Quasi-Regular
* Pseudo-Random with a cutoff
* Pseudo-Random with Rejection
* Quasi-Random (Sobol Sequence) with a cutoff
* Quasi-Random (Sobol Sequence) with Rejection
* Direct-Product

As discussed in the manuscript the Distribution Function (P) we use is defined
as:
P(x):= c^{-1}(E_{cut}-V(x))^{d/2}

## Analysis
Various scripts I have generated for studying/analyzing the accuracy of these
methods.

## miscellaneous
Various scripts I have used for understanding the project.

## External Codes
External Fortran modules were used during this project. 
A detailed list is provided below along with links to the appropriate source code.

#### Quasi-Random Sequences:
The [Sobol Sequence Generator](https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html)
was used for the Quasi-Random Calculations (sobol.f90).
This code is under a [GNU LGPL License](https://www.gnu.org/licenses/lgpl-3.0.en.html)
and has been included here for convenience.

Special thanks to John Burkardt for taking the time to discuss various
quasi-random number generators and their implementation.

#### Evaluation of the Potential Energy Matrix Elements:
We used Gauss-Hermite Quadrature to evaluate the Potential Energy Matrix
Elements.
Again we have used code made available by
[John Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/gen_hermite_rule/gen_hermite_rule.html)
(gen_hermite_rule.f90) which is also under a [GNU LGPL License](https://www.gnu.org/licenses/lgpl-3.0.en.html)
and has been included here for convenience.

The Monte Carlo and Quasi-Monte Carlo Wiki page
([MCQMC](http://roth.cs.kuleuven.be/wiki/Main_Page)) has been a useful resource
for generating and applying various Monte Carlo schemes.

## Scientific Context
Accurate RoVibrational Energy calculations have a long history in computational
chemistry/physics.
The following references could be helpful for understanding the context for 
this project.

#### Distributed Gaussian Basis Sets
The project was motivated in part by the work of
[Garashchuk and Light](https://aip.scitation.org/doi/abs/10.1063/1.1348022).
Special thanks to Sophya Garashcuk for discussing her previous results with us.

The work by [Poirier and Light](https://aip.scitation.org/doi/abs/10.1063/1.481787)
is a good reference for using Distributed Gaussian Basis sets in the context of
RoVibrational spectroscopy.

Interested readers should also consider modern work in the field.
[Bill Poirier](https://aip.scitation.org/doi/full/10.1063/1.4769402) and
[Tucker Carrington's](https://aip.scitation.org/doi/full/10.1063/1.3246593) work
is a good place to start.

Special thanks to both Bill Poirier and Tucker Carrington Jr. for their insight
and suggestions during the early phases of this project.

## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2019. UCI
