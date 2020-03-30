#=============================================================================80
#               Distributed Gaussian Basis (Eigen Error Analysis)
#==============================================================================#
#       Discussion:
#Python script for computing relative error for the final iteration of eigenvalues
#'Theory' values are passed in (computed elsewhere, already sorted)
#Error is taken as (calc-theory)/theory
#first column is assumed to be x axis (alpha0 or some other value, not an 
#eigenvalue
#also assume only theory resutls is a column of data
#==============================================================================#
#       Modified:
#   8 May 2019
#       Author:
#   Shane Flynn 
#==============================================================================#
import pandas as pd
import numpy as np
from sys import argv
#==============================================================================#
#               Discussion:
#script             ==> Name of this script (error.py)
#eigenvalues        ==> datafile containing computed eigenvalues
#theory_values      ==> datafile containing theoretical eigenvalues
#df_eigen           ==> python dataframe for eigenvalues
#df_theory          ==> python dataframe for theory_values
#df_sort            ==> python dataframe with sorted theory_values
#errors             ==> difference between computed/theory eigenvalues
#index              ==> eigenvalue index (for plotting convenience)
#df_analysis        ==> store all data to write to file (for convenience)
#==============================================================================!
#                    Read In Eigenvalues and Theory Values
#first value is alpha 0 not an eigenvalue want to keep this
#==============================================================================!
script,eigenvalues,theory_values=argv
df_eigen=pd.read_csv(eigenvalues,delim_whitespace=True,header=None,dtype=np.float64)
df_theory=pd.read_csv(theory_values, delim_whitespace=True,header=None,dtype=np.float64)
print 'Input Eigenvalues'
print df_eigen.head()
print 'theory input'
print df_theory.head()
#==============================================================================!
#           Determine the Number of Iterations and Eigenvalues Computed 
#first column is alpha0 not na eigenvalue
#==============================================================================!
num_iter = int(df_eigen.shape[0])  #number of rows
num_eig_val = int(df_eigen.shape[1])
print'Number of Iterations => ', num_iter
print'Number of Eigenvalues => ', num_eig_val
size_theory = int(df_theory.shape[0])
#==============================================================================!
#                   Remove extra sorted theory values
#Fortran code computes all combinatorics for i,j index
#Only compute i eigenvalues (true=i^2, but eigenvalues=i)
#So remove tail (theory - number of eigenvalues computed)
#==============================================================================!
df_theory.drop(df_theory.tail((size_theory-num_eig_val)).index,inplace=True)
print'sorted and truncated data frame'
print df_theory.head()
#==============================================================================!
#                    Compute Error ===> abs(computed-theory)
#want to compute error for each eigenvalue and each iteration
#==============================================================================!
errors = []
errors.append(df_eigen[0][0])
for i in range(0,num_iter):
#    for j in range(0,num_eig_val):
#start on second column since first marks alpha0
    for j in range(1,num_eig_val):
        errors.append((df_eigen[j][i] - df_theory[0][j-1])/df_theory[0][j-1])
#        errors.append(abs(df_eigen[j][i] - df_sort[0][j]))
#==============================================================================!
#                       Write Data to CSV File
#==============================================================================!
df_error=pd.DataFrame(np.array(errors).reshape(num_iter,num_eig_val))
print df_error.head()
df_error.to_csv('rel_error.dat',sep='\t',index=False)
#df_analysis=pd.concat([df_eigen, df_sort, df_error],axis=1)
#df_analysis=pd.concat([df_index, df_eigen, df_sort, df_error],axis=1)
#df_analysis.to_csv('error.dat',sep='\t',index=False)
