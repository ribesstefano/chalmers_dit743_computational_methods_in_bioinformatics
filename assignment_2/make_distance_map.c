/*
 * File:  residue_array.c
 * Purpose:  Read PDB atom records into an array of "residue" structures.
 */
#include "pdb_handler.h"
#include "atom.h"
#include "residue.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
  Residue* residues;
  int previousSeq;
  char previousID[2];
  char previousICode[2];
} user_data_t;

Atom get_ca_from_residue(const Residue residue);
void residue_callback(const PdbEntry* entry, int* line_idx, void* data);

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: residue_array file.pdb [threshold=7]\n");
    exit(1);
  }
  double dist_threshold = 7.0;
  if (argc == 3) {
    dist_threshold = atof(argv[2]);
  }
  int num_residues;
  Residue residues[MAX_RESIDUES + 1];
  user_data_t data;
  data.residues = residues;
  data.previousSeq = 0;
  strcpy(data.previousID, "");
  strcpy(data.previousICode, "");
  num_residues = read_data(argv[1], &residue_callback, (void*)&data);
  for (int i = 1; i <= num_residues; ++i) {
    const Atom a = get_ca_from_residue(residues[i]);
    for (int j = 1; j <= num_residues; ++j) {
      const Atom b = get_ca_from_residue(residues[j]);
      const double distance = get_atoms_distance(a, b);
      if (distance < dist_threshold) {
        printf("%d %d\n", i, j);
      }
    }
  }
  return 0;
}

Atom get_ca_from_residue(const Residue residue) {
  const char ca_str[5] = " CA ";
  for (int i = 1; i <= residue.numAtoms; ++i) {
    if (strcmp(residue.atom[i].atomName, ca_str) == 0) {
      return residue.atom[i];
    }
  }
  // Should never get here.
  fprintf(stderr, "ERROR. Unable to find CA atom in given residue.\n");
  fprintf(stderr, "ERROR. List of atoms in the residue follows:\n");
  for (int i = 1; i <= residue.numAtoms; ++i) {
    fprintf(stderr, "ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n",
      residue.atom[i].serial, residue.atom[i].atomName, residue.atom[i].altLoc,
      residue.resName, residue.chainID, residue.resSeq, residue.iCode,
      residue.atom[i].centre.x, residue.atom[i].centre.y,
      residue.atom[i].centre.z);
  }
  fprintf(stderr, "ERROR. Exiting.\n");
  exit(2);
}

void residue_callback(const PdbEntry* entry, int* line_idx, void* data) {
  // printf("I'm residue_callback!\n");
  int i = *line_idx;
  user_data_t* d = (user_data_t*)data;
  Residue* residues = d->residues;
  char* previousID = d->previousID;
  char* previousICode = d->previousICode;
  /*
   * Copy values to the next element in the array.
   */
  int serial = atoi(entry->s_serial);
  int resSeq = atoi(entry->s_resSeq);
  double x = atof(entry->s_x);
  double y = atof(entry->s_y);
  double z = atof(entry->s_z);
  const bool kNewSequence = resSeq != d->previousSeq;
  const bool kNewChanID = strcmp(entry->s_chainID, previousID) != 0;
  const bool kNewICode = strcmp(entry->s_iCode, previousICode) != 0;
  if (kNewSequence || kNewChanID || kNewICode) {
    if (i > MAX_RESIDUES) {
      fprintf(stderr, "Too many residues\n");
      exit(0);
    }
    ++(*(line_idx));
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
  /*
   * Add C_alpha atoms only.
   */
  if (strcmp(entry->s_name, " CA ") == 0) {
    ++(residues[i].numAtoms);
    if (residues[i].numAtoms > MAX_ATOMS_PER_RESIDUE) {
      fprintf(stderr, "ERROR. Too many atoms in residues %d. Exiting.\n", i);
      exit(0);
    }
    int j = residues[i].numAtoms;
    residues[i].atom[j].serial = serial;
    strcpy(residues[i].atom[j].atomName, entry->s_name);
    strcpy(residues[i].atom[j].altLoc, entry->s_altLoc);
    residues[i].atom[j].centre.x = x;
    residues[i].atom[j].centre.y = y;
    residues[i].atom[j].centre.z = z;
  }
}