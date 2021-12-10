#ifndef ATOM_H_
#define ATOM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

typedef struct {
  double x, y, z;
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

double get_distance(const Point a, const Point b);
double get_atoms_distance(const Atom a, const Atom b);
bool is_heavy_atom(const char* atom_name);
void print_pdb_atom (const int serial, const char* s_name, const char* s_altLoc,
  const char* s_resName, const char* s_chainID, const int resSeq,
  const char* s_iCode, const Point centre);

#ifdef __cplusplus
}
#endif
#endif // end ATOM_H_