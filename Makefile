this= .
PROG = snapdragons
O = ./src

SRCS =  $(O)/PopSynth.F $(O)/Coord_trans.F $(O)/Gaia-errors.F

OBJS = $(O)/PopSynth.o $(O)/Coord_trans.o $(O)/Gaia-errors.o

LIBS =

FC = gfortran

#FFLAGS  = -O3 -i-dynamic -convert big_endian -mcmodel=large -auto-scalar
#LDFLAGS = -O3 -i-dynamic -convert big_endian -mcmodel=large -auto-scalar
FFLAGS = -O3
LDFLAGS = -O3

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(this)/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

LIB =

$(O)/Main.o: PopSynth.F parameters.h Mockinioptions.h $(depall)
	$(FC) $(FFLAGS) -c PopSynth.F $(LIB)

$(O)/READ_HART_cosmo.o: Gaia-errors.F const_math.h const_astro.h $(depall)
	$(FC) $(FFLAGS) -c Gaia-errors.F $(LIB)

$(O)/trans3.o: Coord_trans.F parameters.h const_math.h const_astro.h $(depall)
	$(FC) $(FFLAGS) -c Coord_trans.F $(LIB)
