#$Id: makefile,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $

include ../makefile.inc

.SUFFIXES: .o .f90 .C .for

.f90.o:
	$(f90) $(f90appflags) -c $< -o $@

.C.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@

.for.o:
	$(f77) -c $< -o $@

# Input files

TwoPunctureFILES = TwoPunctureABE.o TwoPunctures.o write_Pablo_file.o \
		   ABp_to_xyz.o 

# file dependences

$(TwoPunctureFILES): TwoPunctures.h comm.h

# projects

All: TwoPunctureABE PN_orbit

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

PN_orbit: PN_orbit.o misc.o lagpolint.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ \
	PN_orbit.o misc.o lagpolint.o $(LDLIBS)

clean:
	rm -f TwoPunctureABE $(TwoPunctureFILES)
