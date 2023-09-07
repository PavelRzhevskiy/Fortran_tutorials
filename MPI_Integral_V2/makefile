FC=mpif90
FFLAGS=-Wall -Wextra -Wimplicit-interface -fmax-errors=1
RELFLAGS=-O3 -funroll-loops -ffast-math -march=native 
SRC=IntegralModule.f90 IntegralUniform.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) $(RELFLAGS) -o $@ -c $<

release: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm *.o *.mod 


