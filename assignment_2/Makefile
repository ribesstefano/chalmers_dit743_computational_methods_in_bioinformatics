# Author: Stefano Ribes
CXX = gcc
CFLAGS = -g -std=c99 -O3 -fopenmp # -Wall 
LDFLAGS = -lm
SRC = atom.c residue.c pdb_handler.c segment.c

all: pdb_io.exe atom_array.exe residue_array.exe make_distance_map.exe domak_partition.exe multi_domak_partition.exe

pdb_io.exe: pdb_io.c
	$(CXX) $(CFLAGS) -o pdb_io.exe pdb_io.c $(SRC) $(LDFLAGS)

atom_array.exe: atom_array.c
	$(CXX) $(CFLAGS) -o atom_array.exe atom_array.c $(SRC) $(LDFLAGS)

residue_array.exe: residue_array.c
	$(CXX) $(CFLAGS) -o residue_array.exe residue_array.c $(SRC) $(LDFLAGS)

make_distance_map.exe: make_distance_map.c
	$(CXX) $(CFLAGS) -o make_distance_map.exe make_distance_map.c $(SRC) $(LDFLAGS)

domak_partition.exe: domak_partition.c
	$(CXX) $(CFLAGS) -o domak_partition.exe domak_partition.c $(SRC) $(LDFLAGS)

multi_domak_partition.exe: multi_domak_partition.c
	$(CXX) $(CFLAGS) -DNO_LOOKUP_TABLE=1 -o multi_domak_partition.exe multi_domak_partition.c $(SRC) $(LDFLAGS)

clean:
	rm -f *.exe *.o *.ps