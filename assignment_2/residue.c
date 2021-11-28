/*
 * File:  residue.c
 * Author: Stefano Ribes
 */
#include "residue.h"
#include <string.h>

bool get_atom_from_residue(const Residue residue, const char* atom_name,
    Atom* atom) {
  bool atom_found = false;
  for (int i = 1; i <= residue.numAtoms; ++i) {
    if (strcmp(residue.atom[i].atomName, atom_name) == 0) {
      *atom = residue.atom[i];
      atom_found = true;
      break;
    }
  }
  return atom_found;
}