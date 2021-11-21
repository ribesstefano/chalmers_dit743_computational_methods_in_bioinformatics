#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "atom.h"
#include <stdbool.h>

#define MAX_ATOMS_PER_RESIDUE 14

typedef struct {
  int numAtoms;
  char resName[4];
  char chainID[2];
  int resSeq;
  char iCode[2];
  Atom atom[MAX_ATOMS_PER_RESIDUE+1];
} Residue;

bool get_atom_from_residue(const Residue residue, const char* atom_name,
  Atom* atom);

#endif // end RESIDUE_H_