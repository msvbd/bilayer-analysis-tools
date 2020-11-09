# Fortran compiler
#FC = ifort
FC = gfortran
# Fortran compiler flags
FF = -O3 -std=f2008
FFDBG = -std=f2008 -g -fcheck=all

# C compiler
CC = gcc
# C compiler flags
CF = 

.PHONY: clean all

all: clean \
	 config_mod.o \
	 lmp_traject_reader_mod.o \
	 gnuplot_mod.o \
	 single_patt_mod.o \
	 mol_pattern_mod.o

%.o: %.f90
	$(FC) $(FFDBG) -c -o $@ $<
	
clean:
	rm -f *.mod *.o
