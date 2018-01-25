FF     = gfortran
IDIR   = ./src
ODIR   = ./src
_FOBJS = pop_synth.o gaia_errors.o tgas_errors.o
FOBJS  = $(patsubst %,$(ODIR)/%,$(_FOBJS))
FFILES = src/pop_synth.F src/gaia_errors.F src/tgas_errors.F
FFLAGS = -O3 -I$(IDIR)
_DEPS  = const_math.h const_ast.h
DEPS   = $(patsubst %,$(IDIR)/%,$(_DEPS))

$(ODIR)/%.o: $(ODIR)/%.F $(DEPS)
	$(FF) -c -o $@ $< $(FFLAGS)

snapdragons: $(FOBJS)
	gfortran -o $@ $^

clean:
	rm -rf *.o rm *~ src/*.o src/*~ snapdragons
