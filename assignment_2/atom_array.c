/*
 * File:  atom_array.c
 * Purpose:  Read PDB atom records into an array of "atom" structures.
 * Author: Stefano Ribes
 */
#include "pdb_handler.h"
#include "atom.h"
#include "residue.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void atom_callback(const PdbEntry* entry, int* line_idx, void* user_data) {
  // printf("I'm atom_callback!\n");
  Atom* atoms = (Atom*)user_data;
  int serial = atoi(entry->s_serial);
  int resSeq = atoi(entry->s_resSeq);
  int i = *(line_idx);
  double x = atof(entry->s_x);
  double y = atof(entry->s_y);
  double z = atof(entry->s_z);
  /*
   * Copy values to the next element in the atom array.
   */
  if (i > MAX_ATOMS - 1) {
    (void) fprintf(stderr, "Too many elements read\n");
    exit(0);
  }
  /*
   * Convert the numeric fields to integers or doubles.
   * The library functions atoi() and atof() are
   * described in the UNIX manual pages ('man atoi' and
   * 'man atof').
   */
  atoms[i].serial = serial;
  strcpy(atoms[i].atomName, entry->s_name);
  strcpy(atoms[i].altLoc, entry->s_altLoc);
  strcpy(atoms[i].resName, entry->s_resName);
  strcpy(atoms[i].chainID, entry->s_chainID);
  atoms[i].resSeq = resSeq;
  strcpy(atoms[i].iCode, entry->s_iCode);
  atoms[i].centre.x = x;
  atoms[i].centre.y = y;
  atoms[i].centre.z = z;
  (*(line_idx))++;
}

int main(int argc, char **argv) {
  int  numAtoms;
  int  i;

  if (argc < 2) {
    (void) fprintf(stderr, "usage: atom_array file.pdb\n");
    exit(0);
  }
  /*
   * Declare an array to hold data read from the ATOM records of a PDB file.
   */
  Atom atoms[MAX_ATOMS + 1];
  numAtoms = read_data(argv[1], &atom_callback, (void*)atoms);
  for (i = 0; i < numAtoms; ++i) {
    print_pdb_atom (
      atoms[i].serial,
      atoms[i].atomName,
      atoms[i].altLoc,
      atoms[i].resName,
      atoms[i].chainID,
      atoms[i].resSeq,
      atoms[i].iCode,
      atoms[i].centre);
  }
  return 0;
}
