 # path to library
OBJ := $(shell find "$$(cd ..; pwd)" -type d -name "fortran_mods")

LIBS =  $(OBJ)/config_mod.o \
		$(OBJ)/lmp_traject_reader_mod.o \
		$(OBJ)/gnuplot_mod.o \
		$(OBJ)/mol_pattern_mod.o \
		$(OBJ)/single_patt_mod.o

# Fortran compiler
#FC = ifort
FC = gfortran
# Fortran compiler flags
FF = -O3 -std=f2008 -J$(OBJ)
FFDBG = -std=f2008 -g -fcheck=all -J$(OBJ)

# C compiler
CC = gcc
# C compiler flags
CF = 

.PHONY: clean all

all: clean \
	 bilayer_Dc_Lz_Doverlap \
	 bilayer_Dol_Dc_Lz_N_width \
	 bilayer_realSpace_rigidity \
	 bilayer_realSpace_rigidity_doubleTail

#bilayer_realSpace_rigidity_doubleTail: bilayer_realSpace_rigidity_doubleTail.f90
#	$(FC) $(FF) -o $@ $< $(LIBS)
#	./$@ asdf lam_vt_diCTAC.cs

%: %.f90
	$(FC) $(FF) -o $@ $< $(LIBS)
#	./$@ asdf lam_vt_CTAC.cs

clean:
	rm -f *.mod *.o \
	bilayer_Dol_Dc_Lz_N_width \
	bilayer_Dc_Lz_Doverlap \
	bilayer_realSpace_rigidity \
	bilayer_realSpace_rigidity_doubleTail
