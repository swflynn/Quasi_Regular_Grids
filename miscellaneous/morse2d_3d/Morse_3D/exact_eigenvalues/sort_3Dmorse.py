#=============================================================================80
#                   Sort 3D Morse Analytic Eigenvalues
#==============================================================================#
#Script for sorting (by magnitude) eigenvalues from 3Dmorse_theory.f90
#==============================================================================#
#       Modified:
#   15 May 2019
#       Author:
#   Shane Flynn 
#==============================================================================#
import pandas as pd
import numpy as np
from sys import argv
#==============================================================================#
#               Discussion:
#script             ==> Name of THIS Python script 
#theory_values      ==> datafile containing 3DMorse eigenvalue combinatorics
#df_theory          ==> dataframe for theory_values
#df_sort            ==> dataframe with theory_values sorted by magnitude
#==============================================================================!
#                           Read In Theory Values
#==============================================================================!
script,theory_values=argv
df_theory=pd.read_csv(theory_values,delim_whitespace=True,header=None,dtype=np.float64)
print('theory input')
print(df_theory.head())
#==============================================================================!
#                           Sort Theory Eigenvalues
#==============================================================================!
df_sort=df_theory.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
print('Sorted Theory')
print(df_sort.head())
#==============================================================================!
#                           Write Data to CSV File
#==============================================================================!
df_sort.to_csv('3Dmorse_sorted.dat',sep='\t',index=False,header=None)
print('Hello Universe')
