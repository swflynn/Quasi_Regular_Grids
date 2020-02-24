# Program Main:
3D Morse Oscillator Implementation

A Quasi-Regular grid is generated according to our target distribution P(x).

## Overview:
A quick summary of the various pieces of code available. 

#### Makefile
A simple Makefile for compiling the program.

#### sobol_stdnormal.f90:
This fortran program takes in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform it to a Gaussian Distribution. 
