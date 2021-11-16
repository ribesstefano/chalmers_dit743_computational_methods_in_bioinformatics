#ifndef DATA_H_
#define DATA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ATOMS       10000
#define LINE_LENGTH     81

#define MAX_RESIDUES    1000
#define MAX_ATOMS_PER_RESIDUE  14

typedef struct {
  double  x, y, z;
} Point;

typedef struct {
  int serial;
  char atomName[5];
  char altLoc[2];
  char resName[4];
  char chainID[2];
  int resSeq;
  char iCode[2];
  Point centre;
} Atom;

typedef struct {
  int numAtoms;
  char resName[4];
  char chainID[2];
  int resSeq;
  char iCode[2];
  Atom atom[MAX_ATOMS_PER_RESIDUE+1];
} Residue;

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

void print_pdb_atom (
    const int serial,
    const char* s_name,
    const char* s_altLoc,
    const char* s_resName,
    const char* s_chainID,
    const int resSeq,
    const char* s_iCode,
    const Point centre) {
  printf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n", serial, s_name,
    s_altLoc, s_resName, s_chainID, resSeq, s_iCode, centre.x, centre.y,
    centre.z);
}

typedef void (*callback_ptr)(const PdbEntry*, int*, void* user_data);

int read_data(const char *filename, const callback_ptr callback, void* user_data) {
  FILE  *stream;
  char line[LINE_LENGTH];
  PdbEntry entry;
  int i = 0;
  if ((stream = fopen(filename, "r")) == NULL) {
    (void) fprintf(stderr, "Unable to open %s\n", filename);
    exit(0);
  }
  while (fgets(line, LINE_LENGTH, stream)) {
    if (strncmp(line, "ATOM  ", 6) == 0) {
      /*
       * Split the line into its constituent fields.
       * We are only interested in columns 1-54.
       */
      strncpy(entry.s_serial,  &line[6],  5); entry.s_serial[5] = '\0';
      strncpy(entry.s_name,    &line[12], 4); entry.s_name[4] = '\0';
      strncpy(entry.s_altLoc,  &line[16], 1); entry.s_altLoc[1] = '\0';
      strncpy(entry.s_resName, &line[17], 3); entry.s_resName[3] = '\0';
      strncpy(entry.s_chainID, &line[21], 1); entry.s_chainID[1] = '\0';
      strncpy(entry.s_resSeq,  &line[22], 4); entry.s_resSeq[4] = '\0';
      strncpy(entry.s_iCode,   &line[26], 1); entry.s_iCode[1] = '\0';
      strncpy(entry.s_x,       &line[30], 8); entry.s_x[8] = '\0';
      strncpy(entry.s_y,       &line[38], 8); entry.s_y[8] = '\0';
      strncpy(entry.s_z,       &line[46], 8); entry.s_z[8] = '\0';
      /*
       * Call given callback function on each entry: each program will have its
       * own definition.
       */
      callback(&entry, &i, user_data);
    }
  }
  return i;
}




#endif // end DATA_H_