/*
 * File:  atom.c
 * Author: Stefano Ribes
 */
#include "atom.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

double get_distance(const Point a, const Point b) {
  const Point diff = {a.x - b.x, a.y - b.y, a.z - b.z};
  return sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
}

double get_atoms_distance(const Atom a, const Atom b) {
  return get_distance(a.centre, b.centre);
}

bool is_heavy_atom(const char* atom_name) {
  if (strcmp(atom_name, " CA ") == 0) {
    return true;
  }
  return false;
}

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