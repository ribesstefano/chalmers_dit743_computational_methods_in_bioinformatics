#ifndef ATOM_H_
#define ATOM_H_

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

#endif // end ATOM_H_