#include "pdb_handler.h"

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