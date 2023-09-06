FC=gfortran
FFLAGS=-Wall -Wextra -Wimplicit-interface -fmax-errors=1
FLIBS=-llapack

SRC=workmodule.f90 main.f90

#debug settings
DBGDIR=debug
DBGTRG=$(DBGDIR)/eigen_dbg
DBGOBJ=${SRC:%.f90=debug/%.o}
DBGFLAGS=-Og -g -fbacktrace -fcheck=all -ffpe-trap=zero,overflow -finit-real=NaN

#release settings
RELDIR=release
RELTRG=$(RELDIR)/eigen_rel
RELOBJ=${SRC:%.f90=release/%.o}
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
	@rm -f $(DBGOBJ) $(DBGTRG) $(RELOBJ) $(RELTRG) ./**/*.mod
