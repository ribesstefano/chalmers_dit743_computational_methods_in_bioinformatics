/*
 * File:  residue_array.c
 * Purpose:  Read PDB atom records into an array of "residue" structures.
 */
#include "data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
  Residue* residues;
  int j;
  int previousSeq;
  char previousID[2];
  char previousICode[2];
} user_data_t;

double get_distance(const Point a, const Point b);
double get_atoms_distance(const Atom a, const Atom b);
bool get_atom_from_residue(const Residue residue, const char* atom_name, Atom* atom);
void residue_callback(const PdbEntry* entry, int* line_idx, void* data);

bool is_heavy(Atom a) {
  if (strcmp(a.atomName, " CA ") == 0) {
    return true;
  }
  return false;
}

typedef struct {
  int start;
  int end;
  int num_internal_contacts;
} segment_t;

void set_interior_count(const double dist_threshold, const int num_residues,
    const Residue* residues, segment_t* segment, double* dist_lookup) {
  const char* atom_name = " CA ";
  double dist;
  Atom a, b;
  for (int i = segment->start; i <= segment->end; ++i) {
    if (get_atom_from_residue(residues[i], atom_name, &a)) {
      for (int j = segment->start; j <= segment->end; ++j) {
        if (get_atom_from_residue(residues[j], atom_name, &b)) {
          // Check if the distance is already in the lookup table (distances
          // cannot be negative). Compute it otherwise (the lookup matrix is
          // symmetric).
          if (dist_lookup[i * num_residues + j] > 0) {
            dist = dist_lookup[i * num_residues + j];
          } else {
            dist = get_atoms_distance(a, b);
            dist_lookup[i * num_residues + j] = dist;
            dist_lookup[j * num_residues + i] = dist;
          }
          if (dist < dist_threshold) {
            (segment->num_internal_contacts)++;
          }
        }
      }
    }
  }
}

int get_num_exterior_count(const double dist_threshold, const int num_residues,
    const Residue* residues, const segment_t a, const segment_t b,
    double* dist_lookup) {
  const char* atom_name = " CA ";
  double dist;
  Atom x, y;
  int num_exterior_contacts = 0;
  for (int i = a.start; i <= a.end; ++i) {
    if (get_atom_from_residue(residues[i], atom_name, &x)) {
      for (int j = b.start; j <= b.end; ++j) {
        if (get_atom_from_residue(residues[j], atom_name, &y)) {
          // Check if the distance is already in the lookup table (distances
          // cannot be negative). Compute it otherwise (the lookup matrix is
          // symmetric).
          if (dist_lookup[i * num_residues + j] > 0) {
            dist = dist_lookup[i * num_residues + j];
          } else {
            dist = get_atoms_distance(x, y);
            dist_lookup[i * num_residues + j] = dist;
            dist_lookup[j * num_residues + i] = dist;
          }
          if (dist < dist_threshold) {
            ++num_exterior_contacts;
          }
        }
      }
    }
  }
  return num_exterior_contacts;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    (void) fprintf(stderr, "ERROR. Usage: residue_array file.pdb [threshold=7]\n");
    exit(1);
  }
  double dist_threshold = 5.0;
  if (argc == 3) {
    dist_threshold = atof(argv[2]);
  }
  int num_residues;
  Residue residues[MAX_RESIDUES + 1];
  user_data_t data;
  data.residues = residues;
  data.j = 0;
  data.previousSeq = 0;
  strcpy(data.previousID, "");
  strcpy(data.previousICode, "");
  // TODO: Modify the callback above to only store the "heavy atoms" of the
  // residues.
  num_residues = read_data(argv[1], &residue_callback, (void*)&data);
  double* split_values = malloc((num_residues+1) * sizeof(double));
#define MAX_NUM_SEGMENTS 2
  int num_exterior_contacts[MAX_NUM_SEGMENTS];
  // NOTE: The lookup table is kept symmetric to exploit data locality as much
  // as possible.
  const int kLookupSize = (num_residues+1) * (num_residues+1);
  double* dist_lookup = malloc(kLookupSize * sizeof(double));
  memset(dist_lookup, -1, kLookupSize * sizeof(double));
  segment_t segments[MAX_NUM_SEGMENTS];
  for (int i = 2; i <= num_residues - 1; ++i) {
    segments[0].start = 1;
    segments[0].end = i;
    segments[1].start = i + 1;
    segments[1].end = num_residues;
    segments[0].num_internal_contacts = 0;
    segments[1].num_internal_contacts = 0;
    set_interior_count(dist_threshold, num_residues, residues, &segments[0], dist_lookup);
    set_interior_count(dist_threshold, num_residues, residues, &segments[1], dist_lookup);
    num_exterior_contacts[0] = get_num_exterior_count(dist_threshold, num_residues, residues,
      segments[0], segments[1], dist_lookup);
    split_values[i] = (double)segments[0].num_internal_contacts * (double)segments[1].num_internal_contacts / (double)num_exterior_contacts[0];
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

double get_distance(const Point a, const Point b) {
  const Point diff = {a.x - b.x, a.y - b.y, a.z - b.z};
  return sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
}

double get_atoms_distance(const Atom a, const Atom b) {
  return get_distance(a.centre, b.centre);
}

bool get_atom_from_residue(const Residue residue, const char* atom_name, Atom* atom) {
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