# C compiler:
CC = gcc
#CC = cc
#CC = icc


C_FOLDER = dream/C_functions/



.DEFAULT_GOAL := _dream_


_dream_:
	$(CC) -fPIC -shared -o $(C_FOLDER)numerical/num_c.so $(C_FOLDER)numerical/num_c.c
	$(CC) -fPIC -shared -o $(C_FOLDER)dark_matter/mergers.so $(C_FOLDER)dark_matter/mergers.c
	$(CC) -fPIC -shared -o $(C_FOLDER)dark_matter/centrals.so $(C_FOLDER)dark_matter/centrals.c
	$(CC) -fPIC -shared -o $(C_FOLDER)analyze.so $(C_FOLDER)analyze.c
	$(CC) -fPIC -shared -o $(C_FOLDER)dark_matter/Mhalo_to_Mstar.so $(C_FOLDER)dark_matter/Mhalo_to_Mstar.c
	$(CC) -fPIC -shared -o $(C_FOLDER)baryonic_matter/central_galaxies.so $(C_FOLDER)baryonic_matter/central_galaxies.c
	$(CC) -fPIC -shared -o $(C_FOLDER)baryonic_matter/satellite_galaxies.so $(C_FOLDER)baryonic_matter/satellite_galaxies.c
	$(CC) -fPIC -shared -o $(C_FOLDER)star_formation/star_formation.so $(C_FOLDER)star_formation/star_formation.c
	$(CC) -fPIC -shared -o $(C_FOLDER)star_formation/star_formation_satellites.so $(C_FOLDER)star_formation/star_formation_satellites.c


clean: 	
	rm $(C_FOLDER)*.so
	rm $(C_FOLDER)numerical/*.so
	rm $(C_FOLDER)dark_matter/*.so
