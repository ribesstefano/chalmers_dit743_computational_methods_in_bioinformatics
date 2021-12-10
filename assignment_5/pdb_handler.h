/*
 * File:  pdb_handler.h
 * Author: Stefano Ribes
 */
#ifndef PDB_HANDLER_H_
#define PDB_HANDLER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "atom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ATOMS 10000
#define MAX_RESIDUES 1000
#define MAX_ATOMS_PER_RESIDUE 14
#define LINE_LENGTH 81

typedef struct {
  char s_serial[6];
  char s_name[5]; // Alpha-carbon is " CA "; calcium is "CA  " (different space)
  char s_altLoc[2]; // Usually " "
  char s_resName[4];
  char s_chainID[2];
  char s_resSeq[5];
  char s_iCode[2]; // Usually " "
  char s_x[9];
  char s_y[9];
  char s_z[9];
} PdbEntry;

typedef void (*callback_ptr)(const PdbEntry*, int*, void* user_data);

int read_data(const char *filename, const callback_ptr callback, void* user_data);

#ifdef __cplusplus
}
#endif
#endif // end PDB_HANDLER_H_