FC=mpif90
FFLAGS=-Wall -Wextra -Wimplicit-interface -fmax-errors=1
SRC=Integral.f90

#debug settings
DBGDIR=debug
DBGTRG=$(DBGDIR)/integral_debug
DBGOBJ=${SRC:%.f90=$(DBGDIR)/%.o}
DBGFLAGS=-Og -g -fbacktrace -fcheck=all -ffpe-trap=zero,overflow -finit-real=NaN

#release settings
RELDIR=release
RELTRG=$(RELDIR)/integral_release
RELOBJ=${SRC:%.f90=$(RELDIR)/%.o}
RELFLAGS=-O3 -funroll-loops -ffast-math -march=native 

#debig build
debug: dirs $(DBGTRG)

$(DBGTRG): $(DBGOBJ)
	$(FC) $(FFLAGS) $(DBGFLAGS) $^ -o $@ $(FLIBS) 

$(DBGOBJ): $(DBGDIR)/%.o: %.f90 
	$(FC) $(FFLAGS) $(DBGFLAGS) -c $< -o $@ 

#release build
release: dirs $(RELTRG)

$(RELTRG): $(RELOBJ)
	$(FC) $(FFLAGS) $(RELFLAGS) $^ -o $@ $(FLIBS) 

$(RELOBJ): $(RELDIR)/%.o: %.f90 
	$(FC) $(FFLAGS) $(RELFLAGS) -c $< -o $@ 




.PHONY: dirs clean
dirs:
	@mkdir -p $(DBGDIR) $(RELDIR)

clean:
	@rm -f $(DBGOBJ) $(DBGTRG) $(RELOBJ) $(RELTRG) 
