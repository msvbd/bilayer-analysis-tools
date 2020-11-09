.PHONY: clean all bilayerTools fortranMods

all: clean \
	 fortranMods \
	 bilayerTools \
	 
fortranMods: fortranMods
	cd fortran_mods && make all

bilayerTools:
	cd bilayer_tools &&	make all
	
clean:
	cd bilayer_tools &&	make clean
	cd fortran_mods  && make clean
