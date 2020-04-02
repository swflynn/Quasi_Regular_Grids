# 2D and 3D Morse Hard-coded Potential Indexing
Before generalizing the code we hard-coded the 2D and 3D morse cases for testing. 
The results are usually consistent with the generalized code for the first 9 
digits. 
The difference come in computing the potential energy matrix in how we do our 
indexing for rr and the subsequent Vij value.

The hard-coded cases do seem to be slighly better upon comparison, however, the 
error is so small in every test case that I am going to proceed with the 
general code. 

These legacy codes are here for reference, the paper actually used these codes for 
the calculations. 
We generalized to begin working on the QRG-Lmon work.
