Author: Stefano Ribes ([ribes@chalmers.se](ribes@chalmers.se))

# Assignment 5: Detecting Steric Overlap

The assignment has been coded in C++, with some libraries from Assignment 2 written in C99.

## Instructions

Download the two PDB files into this directory. The `main` can be found in `detect_steric_overlap.cc`. To compile the program, type:

```bash
make
```

To run the program, type:
```bash
./detect_steric_overlap.exe 1cdh.pdb 2csn.pdb
```

## Changelog

* The function `read_data` in `pdb_handler.c` has been modified to include _HETATM_ entries in the PDB file.

## Implementation Details

## Outputs
