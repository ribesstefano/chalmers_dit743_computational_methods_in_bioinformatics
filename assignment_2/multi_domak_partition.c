/*
 * File:  residue_array.c
 * Purpose:  Read PDB atom records into an array of "residue" structures.
 * Author: Stefano Ribes
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
#include <omp.h>

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
  const int kMaxNumDomains = 40;
  int num_residues;
  Residue residues[MAX_RESIDUES + 1];
  user_data_t data;
  data.residues = residues;
  data.previousSeq = 0;
  strcpy(data.previousID, "");
  strcpy(data.previousICode, "");
  num_residues = read_data(argv[1], &residue_callback, (void*)&data);
  for (int i = 1; i <= num_residues; ++i) {
    if (!residues[i].numAtoms) {
      fprintf(stderr, "ERROR. Residue n.%d doesn't contain any heavy atom (CA). Exiting.\n", i);
      exit(3);   
    }
  }
  // NOTE: The lookup table is kept symmetric to exploit data locality as much
  // as possible.
  const int kLookupSize = (num_residues+1) * (num_residues+1);
  double* dist_lookup = malloc(kLookupSize * sizeof(double));
  memset(dist_lookup, -1, kLookupSize * sizeof(double));
  // Init domains with the whole list of residues.
  Domain domains[kMaxNumDomains];
  int num_domains = 0; // starting indexing from zero.
  domains[0].segments[0].start = 1;
  domains[0].segments[0].end = num_residues;
  domains[0].segments[1] = null_segment;
  domains[0].num_segments = 1;
  int num_splitted_domains = 1;
  int iter_cnt = 0;
  while (num_splitted_domains > 0) {
    printf("[DEBUG] Iteration n.%d.\n", ++iter_cnt);
    num_splitted_domains = 0;
    const int curr_num_domains = num_domains+1;
    for (int i = 0; i < curr_num_domains; ++i) {
      const int total_len = len(domains[i].segments[0]) + len(domains[i].segments[1]);
      if (total_len >= DOMAK_MDSP) {
        if (domains[i].num_segments == 1) {
          printf("[DEBUG] Calling single_segment_scan.\n");
          single_segment_scan(dist_threshold, num_residues,
            residues, i, domains, &num_domains, dist_lookup);
        } else {
          printf("[DEBUG] Calling two_segment_scan_of_two_segment_domain.\n");
          two_segment_scan_of_two_segment_domain(dist_threshold, num_residues,
            residues, i, domains, &num_domains, dist_lookup);
        }
        if (num_domains+1 > curr_num_domains) {
          // If the above methods have extended the domains size, then repeat.
          ++num_splitted_domains;
        }
      }
    }
  }
  printf("[INFO] Found %d domain. Each segment is defined as: (residue start index, residue end index)\n", num_domains+1);
  for (int i = 0; i < num_domains+1; ++i) {
    printf("Domain n.%d\n", i+1);
    for (int j = 0; j < domains[i].num_segments; ++j) {
      printf("\tsegment n.%d: (%d, %d)\n", j+1, domains[i].segments[j].start,
        domains[i].segments[j].end);
    }
  }
  free(dist_lookup);
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