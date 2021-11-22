/*
 * File:  residue_array.c
 * Purpose:  Read PDB atom records into an array of "residue" structures.
 */
#include "pdb_handler.h"
#include "atom.h"
#include "residue.h"
#include "segment.h"

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

void residue_callback(const PdbEntry* entry, int* line_idx, void* data);

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "ERROR. Usage: residue_array file.pdb [threshold=7]\n");
    exit(1);
  }
  double dist_threshold = 5.0;
  if (argc == 3) {
    dist_threshold = atof(argv[2]);
  }
  const int kMaxNumSegments = 2;
  int num_residues;
  Residue residues[MAX_RESIDUES + 1];
  user_data_t data;
  data.residues = residues;
  data.previousSeq = 0;
  strcpy(data.previousID, "");
  strcpy(data.previousICode, "");
  num_residues = read_data(argv[1], &residue_callback, (void*)&data);
  for (int i = 1; i <= num_residues; ++i) {
    if (residues[i].numAtoms == 0) {
      fprintf(stderr, "ERROR. Residue n.%d doesn't contain any heavy atom (CA). Exiting.\n", i);
      exit(3);   
    }
  }
  double* split_values = malloc((num_residues+1) * sizeof(double));
  int num_exterior_contacts[kMaxNumSegments];
  // NOTE: The lookup table is kept symmetric to exploit data locality as much
  // as possible.
  const int kLookupSize = (num_residues+1) * (num_residues+1);
  double* dist_lookup = malloc(kLookupSize * sizeof(double));
  memset(dist_lookup, -1, kLookupSize * sizeof(double));
  Segment segments[kMaxNumSegments];
  for (int i = 2; i <= num_residues - 1; ++i) {
    segments[0].start = 1;
    segments[0].end = i;
    segments[1].start = i + 1;
    segments[1].end = num_residues;
    segments[0].num_internal_contacts = 0;
    segments[1].num_internal_contacts = 0;
    set_int_cnt_for_atom(dist_threshold, num_residues, residues, " CA ",
      &segments[0], dist_lookup);
    set_int_cnt_for_atom(dist_threshold, num_residues, residues, " CA ",
      &segments[1], dist_lookup);
    num_exterior_contacts[0] = get_ext_cnt_for_atom(dist_threshold,
      num_residues, residues, " CA ", segments[0], segments[1], dist_lookup);
    split_values[i] = (double)segments[0].num_internal_contacts * (double)segments[1].num_internal_contacts / ((double)num_exterior_contacts[0] * (double)num_exterior_contacts[0]);
  }
  int split_idx = 0;
  double max_split = -1;
  for (int i = 2; i <= num_residues - 1; ++i) {
    if (split_values[i] > max_split) {
      max_split = split_values[i];
      split_idx = i;
    }
  }
  printf("[INFO] Maximum split value: %.3f, corresponding to index %d.\n", max_split, split_idx);
  printf("[INFO] Bar plot with normalized values:\n");
  for (int i = 2; i <= num_residues - 1; ++i) {
    double norm = split_values[i] / max_split;
    printf("%5d ", i);
    for (int j = 0; j < (int)(norm * 100.0); ++j) {
      // printf("%c", (char)254u);
      printf("*");
    }
    printf(" %.2f", norm * 100);
    printf("\n");
  }
  free(dist_lookup);
  free(split_values);
  return 0;
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
   * Add "heavy" atoms only.
   */
  if (is_heavy_atom(entry->s_name)) {
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