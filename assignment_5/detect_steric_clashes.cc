/*
 * File:  detec_steric_overlap.cc
 * Author: Stefano Ribes
 */
extern "C" {
#include "pdb_handler.h"
#include "atom.h"
}

#include <iostream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <map>
#include <stdio.h>
#include <cmath>
#include <tuple>
#include <utility>
#include <sstream>

/**
 * @brief      This struct wraps an hash function for Point coordinates class.
 */
struct PointHashFunc {
  size_t operator() (const Point &k) const {
    size_t h1 = std::hash<double>()(k.x);
    size_t h2 = std::hash<double>()(k.y);
    size_t h3 = std::hash<double>()(k.z);
    return (h1 ^ (h2 << 1)) ^ h3;
  }
};

/**
 * @brief      This struct wraps a function for determining if two Point objects
 *             are equal. Used in hash table for having a Point class as key
 *             type.
 */
struct PointEqualsFunc {
  bool operator() (const Point& lhs, const Point& rhs) const {
    return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
  }
};

/**
 * @brief      This class describes a protein.
 */
class Protein {
public:
  typedef std::unordered_map<Point, std::vector<Atom>, PointHashFunc,
    PointEqualsFunc> HashTableType;
  static const int kAtomRadius = 2;

  /**
   * @brief      Callback function to pass to read_data(), which is used for
   *             reading PDB files. Populate a vector of atoms.
   *
   * @param[in]  entry      The read PDB entry
   * @param      line_idx   The line index of the PDB file
   * @param      user_data  The user data
   */
  static void atom_callback(const PdbEntry* entry, int* line_idx,
      void* user_data) {
    std::vector<Atom>* atoms = static_cast<std::vector<Atom>*>(user_data);
    // Convert the numeric fields to integers or doubles.
    Atom a;
    a.serial = atoi(entry->s_serial);
    a.resSeq = atoi(entry->s_resSeq);
    strcpy(a.atomName, entry->s_name);
    strcpy(a.altLoc, entry->s_altLoc);
    strcpy(a.resName, entry->s_resName);
    strcpy(a.chainID, entry->s_chainID);
    strcpy(a.iCode, entry->s_iCode);
    a.centre.x = double(atof(entry->s_x));
    a.centre.y = double(atof(entry->s_y));
    a.centre.z = double(atof(entry->s_z));
    atoms->push_back(a);
    (*line_idx)++; // Advance to the next line in the PDB file.
  }

  /**
   * @brief      Gets the cube coordinates of an Atom, scaled by the minimum
   *             dimensions of a given Protein.
   *
   * @param[in]  a     The input atom
   * @param[in]  min   The minimum coordinates to scale to
   *
   * @return     The cube coordinate as a Point class.
   */
  static Point GetCubeCoord(const Atom& a, const Point& min) {
    // Use the minimum coordinates of the given protein as a sort of origin to
    // scale the atom original coordinates. Then divide by the cube size (i.e.
    // the sphere diameter). Finally round up if positive coordinates, round
    // down otherwise.
    const double kCubeSize = double(Protein::kAtomRadius * 2);
    Point cube;
    // NOTE: Add one to avoid mapping wrong atoms to the same cube index, i.e.
    // avoind rounding errors.
    cube.x = (a.centre.x - min.x) / kCubeSize + 1;
    cube.y = (a.centre.y - min.y) / kCubeSize + 1;
    cube.z = (a.centre.z - min.z) / kCubeSize + 1;
    cube.x = (cube.x > 0) ? std::ceil(cube.x) : std::floor(cube.x);
    cube.y = (cube.y > 0) ? std::ceil(cube.y) : std::floor(cube.y);
    cube.z = (cube.z > 0) ? std::ceil(cube.z) : std::floor(cube.z);
    return cube;
  }

  /**
   * @brief      Detect and print out number of steric overlaps and number of
   *             atom-atom comparisons.
   *
   * @param      a_protein  The A protein
   * @param      b_protein  The B protein
   * @param[in]  use_hash   Whether to use an hash table. Brute force otherwise
   */
  static void DetectStericOverlaps(Protein& a_protein, Protein& b_protein,
      bool use_hash) {
    const double kCubeSize = Protein::kAtomRadius * 2;
    std::map<int, Atom> clashes;
    int num_clashes = 0;
    int num_comparisons = 0;
    if (use_hash) {
      const auto a_atoms = a_protein.get_atoms();
      Point min_coords = b_protein.get_min_coods();
      // Use protein B min dimensions to store the cubes in the hash tables.
      b_protein.MapAtomsToCubes(min_coords);
      auto hashmap = b_protein.get_hash_table();
      // Loop over all atoms in Protein A.
      for (auto a = a_atoms.begin(); a != a_atoms.end(); ++a) {
        // Get the coordinates of the cube containing atom A.
        auto atom_coords = Protein::GetCubeCoord(*a, min_coords);
        // Check all cubes around the cube containing atom A and the cube
        // itself. Check in the range -1 to 1 for each of the thee dimensions.
        for (int i = -1; i < 2; ++i) {
          for (int j = -1; j < 2; ++j) {
            for (int k = -1; k < 2; ++k) {
              Point coords = {
                atom_coords.x + i,
                atom_coords.y + j,
                atom_coords.z + k
              };
              if (hashmap.find(coords) == hashmap.end()) {
                // Cube coordinates of neighbouring atom of Protein A are not in
                // the hashmap of Protein B. Continue.
                continue;
              } else {
                const auto b_atoms = hashmap[coords];
                for (auto b = b_atoms.begin(); b != b_atoms.end(); ++b) {
                  ++num_comparisons;
                  if (get_atoms_distance(*a, *b) < kCubeSize) {
                    ++num_clashes;
                    clashes[b->serial] = *b;
                  }
                }
              }
            }
          }
        }
      }
    } else {
      const auto a_atoms = a_protein.get_atoms();
      const auto b_atoms = b_protein.get_atoms();
      for (auto a = a_atoms.begin(); a != a_atoms.end(); ++a) {
        for (auto b = b_atoms.begin(); b != b_atoms.end(); ++b) {
          ++num_comparisons;
          if (get_atoms_distance(*a, *b) < kCubeSize) {
            ++num_clashes;
            clashes[b->serial] = *b;
          }
        }
      }
    }
    for (auto i = clashes.begin(); i != clashes.end(); ++i) {
      const auto atom = i->second;
      std::cout << atom.serial << " " << atom.resName << " "
                << atom.resSeq << " " << atom.atomName
                << std::endl;
    }
    std::cout << "[INFO] Number of clashing atoms:   " << clashes.size()
              << std::endl;
    std::cout << "[INFO] Number of comparisons made: " << num_comparisons
              << std::endl;
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  pdb_file  The pdb file
   */
  Protein(const char* pdb_file) {
    this->max_coords_ = {-1e308, -1e308, -1e308};
    this->min_coords_ = {1e308, 1e308, 1e308};
    const int num_atoms = read_data(pdb_file, &this->atom_callback,
      static_cast<void*>(&this->atoms_));
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      this->UpdateMaxCoordinates(a->centre);
      this->UpdateMinCoordinates(a->centre);
    }
  }

  ~Protein() {}

  /**
   * @brief      For each atom in the protein, get the coordinates of its
   *             surrounding cube. Then store all atoms into the hash table
   *             using the coordinates as keys.
   *
   * @param[in]  min_coords  The minimum coordinates of a given Protein
   */
  void MapAtomsToCubes(const Point& min_coords) {
    // For each atom in the protein, get the coordinates of its surrounding
    // cube. Then store all atoms into the hash table using the coordinates as
    // keys.
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      const Point cube_coord = Protein::GetCubeCoord(*a, min_coords);
      this->hash_table_[cube_coord].push_back(*a);
    }
  }

  /**
   * @brief      Given a Point, update the Point member keeping track of the
   *             maximum coordinates of the Protein.
   *
   * @param[in]  p     The centre coordinates of an Atom
   */
  void UpdateMaxCoordinates(const Point p) {
    if (p.x > this->max_coords_.x) {
      this->max_coords_.x = p.x;
    }
    if (p.y > this->max_coords_.y) {
      this->max_coords_.y = p.y;
    }
    if (p.z > this->max_coords_.z) {
      this->max_coords_.z = p.z;
    }
  }

  /**
   * @brief      Given a Point, update the Point member keeping track of the
   *             minimum coordinates of the Protein.
   *
   * @param[in]  p     The centre coordinates of an Atom
   */
  void UpdateMinCoordinates(const Point p) {
    if (p.x < this->min_coords_.x) {
      this->min_coords_.x = p.x;
    }
    if (p.y < this->min_coords_.y) {
      this->min_coords_.y = p.y;
    }
    if (p.z < this->min_coords_.z) {
      this->min_coords_.z = p.z;
    }
  }

  /**
   * @brief      Prints all atoms in the protein.
   */
  void PrintAtoms() {
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      print_pdb_atom(a->serial, a->atomName, a->altLoc, a->resName, a->chainID,
        a->resSeq, a->iCode, a->centre);
    }
  }

  HashTableType get_hash_table() {
    return this->hash_table_;
  }

  Point get_min_coods() {
    return this->min_coords_;
  }

  Point get_max_coods() {
    return this->max_coords_;
  }

  std::vector<Atom> get_atoms() {
    return this->atoms_;
  }

private:
  std::vector<Atom> atoms_;
  Point max_coords_;
  Point min_coords_;
  HashTableType hash_table_;
};


int main(int argc, char const *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "ERROR. Usage: detect_steric_clashes.exe file1.pdb file2.pdb [use_bruteforce]\n");
    exit(1);
  }
  bool use_hash = true;
  if (argc >= 4) {
    use_hash = false;
  }
  Protein a = Protein(argv[1]);
  Protein b = Protein(argv[2]);
  Protein::DetectStericOverlaps(a, b, use_hash);
  return 0;
}