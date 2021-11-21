# Assignment 2: Domain assignment

## Purpose

The C programs `pdb_io.c`, `atom_array.c` and `residue_array.c`  demonstrate file handling and input/output.  


## Description

The program `pdb_io` takes the name of a Protein Data Bank (PDB) file as a
command line argument, reads each line of that file, splits up each ATOM line
into individual fields (placing their values in separate variables), then
writes out the contents of each ATOM line in three ways:

	1. (most of) the individual fields, one field per line;
	2. the single string containing the whole ATOM line (as read);
	3. (most of) the individual fields, in PDB format.

The program `atom_array` first reads PDB atom records into an array of
"atom" structures, then writes out the content of this array in PDB format.

The program `residue_array` first reads PDB atom records into an array of
"residue" structures, then writes out the content of this array in PDB format.

Later in this course you will write programs that have to read/write
Protein Data Bank files, and you might find it useful to use these
programs as a starting point.


These programs can be compiled using the Makefile (see below).
You don't *need* to understand (or even look at) the contents of the
Makefile, or how the `make` utility works.  But, if you are interested,
there are some explanatory comments in the Makefile.


## Instructions

Copy the files from this directory into a directory of your own.
I recommend that you use a new directory for each practical.

To compile the programs, type:

```bash
make
```

To run the programs with data file '1CRN.pdb', type:

```bash
./pdb_io 1CRN.pdb
./atom_array 1CRN.pdb
./residue_array 1CRN.pdb
```

To remove those files that can be recompiled from the source code, type:

```bash
make clean
```

## Improvements

List of modifications:

* Cleaned up formatting and style
* Created a separate header file with all data structure definitions
* Removed all global variables
* Array indexing now starts from zero in `atom_array.c`
* Written generic `read_data` function: it now accepts a function pointer to a callback function. Such callback expects the PDB line index, the read PDB entry and any user data. The prototype of the callback function looks something like the following:

```c
typedef void (*callback_ptr)(const PdbEntry*, int*, void* user_data); // Type definition

void atom_callback(const PdbEntry* entry, int* line_idx, void* user_data); // Example of definition

int main(int argc, char **argv) {
  int numAtoms, i;
  if (argc < 2) {
    (void) fprintf(stderr, "usage: atom_array file.pdb\n");
    exit(0);
  }
  numAtoms = read_data(argv[1], &atom_callback, (void*)atom);
  for (i = 0; i < numAtoms; ++i) {
    write_pdb_atom(
      atom[i].serial,
      atom[i].atomName,
      atom[i].altLoc,
      atom[i].resName,
      atom[i].chainID,
      atom[i].resSeq,
      atom[i].iCode,
      atom[i].centre);
  }
  return 0;
}

```
