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
#	make pdb_io
#	make atom_array
#	make residue_array
#
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
CFLAGS= -Wall -g -std=c99

all:	pdb_io atom_array residue_array

pdb_io:	pdb_io.c
	$(CXX) $(CFLAGS) -o pdb_io pdb_io.c

atom_array:	atom_array.c
	$(CXX) $(CFLAGS) -o atom_array atom_array.c

residue_array:	residue_array.c
	$(CXX) $(CFLAGS) -o residue_array residue_array.c

clean:
	rm -f pdb_io
	rm -f atom_array
	rm -f residue_array
