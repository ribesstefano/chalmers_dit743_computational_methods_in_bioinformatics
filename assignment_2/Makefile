#
# Instructions
# ------------
#
# To compile the programs 'pdb_io', 'atom_array' and 'residue_array', type:
#
#	make
# or
#	make all
# or
#	make pdb_io.exe #	make atom_array.exe #	make residue_array.exe #
# To remove those files that can be recompiled from the source code, type:
#
#	make clean
#
#
#
# Description
# -----------
#
# This 'Makefile' contains instructions for creating the executable programs
# 'pdb_io', 'atom_array' and 'residue_array' from the C source file 'pdb_io.c',
# 'atom_array.c' and 'residue_array.c' using the GNU C compiler 'gcc'.
#
# When the user types 'make' at the UNIX prompt, the 'make' program looks
# for the first target in the 'Makefile', and builds that target.
# Targets are given at the start of a line, and are followed by a colon (':')
# then a tab character(s), then a list of files on which the target depends
# (the "dependencies").  The target will be built from the dependencies using
# the command(s) on the following line(s).
# If the target already exists, and is newer than all of the files on which
# it depends, then the target will not be rebuilt.
# For more information, see the UNIX manual page for 'make' (type 'man make').
#
# In this case, target 'all' depends on the executable files 'pdb_io',
# 'atom_array' and 'residue_array'. These in turn depend on the C source files.
#
# Compiler options used:
#
#	-Wall	report possible programming errors in the C source
#	-g	produce debugging information
#	-o	put the output (the executable file) in named file
#

CXX = gcc
CFLAGS = -Wall -g -std=c99
LDFLAGS = -lm
SRC = atom.c residue.c pdb_handler.c

all: pdb_io.exe atom_array.exe residue_array.exe make_distance_map.exe domak_partition.exe

pdb_io.exe:	pdb_io.c
	$(CXX) $(CFLAGS) -o pdb_io.exe pdb_io.c $(SRC) $(LDFLAGS)

atom_array.exe:	atom_array.c
	$(CXX) $(CFLAGS) -o atom_array.exe atom_array.c $(SRC) $(LDFLAGS)

residue_array.exe:	residue_array.c
	$(CXX) $(CFLAGS) -o residue_array.exe residue_array.c $(SRC) $(LDFLAGS)

make_distance_map.exe:	make_distance_map.c
	$(CXX) $(CFLAGS) -o make_distance_map.exe make_distance_map.c $(SRC) $(LDFLAGS)

domak_partition.exe:	domak_partition.c
	$(CXX) $(CFLAGS) -o domak_partition.exe domak_partition.c $(SRC) $(LDFLAGS)

clean:
	rm -f *.exe *.o *.ps