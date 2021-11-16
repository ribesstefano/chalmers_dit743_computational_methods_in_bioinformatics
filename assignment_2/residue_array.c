/*
 * File:  residue_array.c
 * Purpose:  Read PDB atom records into an array of "residue" structures.
 */
#include "data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  Residue* residues;
  int j;
  int previousSeq;
  char previousID[2];
  char previousICode[2];
} user_data_t;

void residue_callback(const PdbEntry* entry, int* line_idx, void* data) {
  // printf("I'm residue_callback!\n");
  int i = *line_idx;
  user_data_t* d = (user_data_t*)data;
  Residue* residues = d->residues;
  char* previousID = d->previousID;
  char* previousICode = d->previousICode;

  int serial = atoi(entry->s_serial);
  int resSeq = atoi(entry->s_resSeq);
  double x = atof(entry->s_x);
  double y = atof(entry->s_y);
  double z = atof(entry->s_z);
  /*
   * Copy values to the next element in the array.
   */
  if (resSeq != d->previousSeq || strcmp(entry->s_chainID, previousID) != 0 ||
       strcmp(entry->s_iCode, previousICode) != 0) {
    if (i > MAX_RESIDUES) {
      (void) fprintf(stderr, "Too many residues\n");
      exit(0);
    }
    (*(line_idx))++;
    ++i;
    d->previousSeq = resSeq;
    strcpy(previousID, entry->s_chainID);
    strcpy(previousICode, entry->s_iCode);
    residues[i].numAtoms = 0;
    strcpy(residues[i].resName, entry->s_resName);
    strcpy(residues[i].chainID, entry->s_chainID);
    residues[i].resSeq = resSeq;
    strcpy(residues[i].iCode, entry->s_iCode);
  }
  if (++(residues[i].numAtoms) > MAX_ATOMS_PER_RESIDUE) {
    (void) fprintf(stderr, "Too many atoms in residues %d\n", i);
    exit(0);
  }
  int j = residues[i].numAtoms;
  residues[i].atom[j].serial = serial;
  strcpy(residues[i].atom[j].atomName, entry->s_name);
  strcpy(residues[i].atom[j].altLoc, entry->s_altLoc);
  residues[i].atom[j].centre.x = x;
  residues[i].atom[j].centre.y = y;
  residues[i].atom[j].centre.z = z;
  d->j = j;
}

int main(int argc, char **argv) {
  int numResidues;
  int i, j;
  if (argc < 2) {
    (void) fprintf(stderr, "usage: residue_array file.pdb\n");
    exit(0);
  }
  Residue residues[MAX_RESIDUES + 1];
  user_data_t data;
  data.residues = residues;
  data.j = 0;
  data.previousSeq = 0;
  strcpy(data.previousID, "");
  strcpy(data.previousICode, "");
  numResidues = read_data(argv[1], &residue_callback, (void*)&data);
  for (i = 1; i <= numResidues; ++i) {
    for (j = 1; j <= residues[i].numAtoms; ++j) {
      print_pdb_atom (
          residues[i].atom[j].serial,
          residues[i].atom[j].atomName,
          residues[i].atom[j].altLoc,
          residues[i].resName,
          residues[i].chainID,
          residues[i].resSeq,
          residues[i].iCode,
          residues[i].atom[j].centre);
    }
  }

  return 0;
}
