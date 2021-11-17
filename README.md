There are 3 files:

1) Gauss.f90 : here a random matrix is defined inside the program, the program applies gauss reduction algorithm to obtain a up-triangular matrix.

2) GaussSub.f90 : is the same algorithm but is expected to be used as a subroutine to be called in the main program (called prova_gauss, prova=try in Ita) 
There you can just cancel "call rand matrix" and put your own LinSis to be Gaussed

3) GJSub2.f90 : is the Gauss-Jordan elimination on a rectangular matrix
(usually n x n+1) where the last column is the known values column)
The output will be the identity matrix merged with the vector solution on the last column. 
Gauss-sub is called as subroutine, it's analogue to the gauss_sub.f90 file. 
"Jordaning" process is very non optimized, don't use it for long computations 
(e.g. 5000x5000)+
I guess it is useful for beginners to understand how it works. 
Exist better optimized versions of this code.
