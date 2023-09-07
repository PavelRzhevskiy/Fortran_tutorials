### Eigens 

- The program reads the matrix from "input.dat" file, solves eigenvalue problem for the matrix and prints eigenvalues in the form (real, imag)
- sample_output.dat contains the program output for some matrix

### Eigens_Makefile

- The same program with subroutines reading matrix from file, finding eigenvalues, and printing eigenvalues in a seperate module
- The Makefile contains two builds with release and debug modes

### MPI_Integral

- The program calculates a particulat 3D integral in rectangular area with openMPI 
- Domain is divided between processes in chunks parallel to x-planes

### MPI_Integral_V2

- Another implementation of the program with the most uniform distribution of domain among processes
- Each point of the domain is referenced through the one index only. But i,j,k are recalculated on each iteration