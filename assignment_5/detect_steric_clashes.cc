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
#include <stdio.h>
#include <cmath>
#include <tuple>
#include <utility>
#include <sstream>

std::string PrintPoint(const Point a) {
  std::stringstream str;
  str << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return str.str();
}

struct PointHashFunc {
    size_t operator() (const Point &k) const {
    size_t h1 = std::hash<double>()(k.x);
    size_t h2 = std::hash<double>()(k.y);
    size_t h3 = std::hash<double>()(k.z);
    return (h1 ^ (h2 << 1)) ^ h3;
    }
};

struct PointEqualsFunc {
  bool operator() (const Point& lhs, const Point& rhs) const {
    return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
  }
};

class Protein {
public:
  static const int kAtomRadius = 2;

  static void atom_callback(const PdbEntry* entry, int* line_idx, void* user_data) {
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

  static Point GetCubeCoord(const Atom& a, const Point& min) {
    const double kCubeSize = Protein::kAtomRadius * 2;
    Point cube;
    cube.x = (a.centre.x - min.x) / kCubeSize;
    cube.y = (a.centre.y - min.y) / kCubeSize;
    cube.z = (a.centre.z - min.z) / kCubeSize;
    cube.x = (cube.x > 0) ? std::ceil(cube.x) : std::floor(cube.x);
    cube.y = (cube.y > 0) ? std::ceil(cube.y) : std::floor(cube.y);
    cube.z = (cube.z > 0) ? std::ceil(cube.z) : std::floor(cube.z);
    return cube;
  }

  static void DetectStericOverlaps(Protein& a_protein, Protein& b_protein, bool use_hash) {
    const int kNumDims = 3;
    const double kCubeSize = Protein::kAtomRadius * 2;
    // a.MapAtomsToCubes(b);
    // a.ResetComparisonCnt();
    // b.ResetComparisonCnt();
    // auto a_hashtable = a.get_hash_table();
    Point min_coords = b_protein.get_min_coods();
    // Use protein B dimensions to label the cubes in the hash tables.
    b_protein.MapAtomsToCubes(min_coords);
    auto hashmap = b_protein.get_hash_table();
    std::vector<std::pair<Atom, Atom> > clashing_atoms;
    std::unordered_map<int, Atom> clashes;
    int num_clashes = 0;
    int num_comparisons = 0;
    if (use_hash) {
      // Loop over all atoms in Protein A.
      const auto a_atoms = a_protein.get_atoms();
      for (auto a = a_atoms.begin(); a != a_atoms.end(); ++a) {
        // Get the coordinates of the cube containing atom A.
        auto atom_coords = Protein::GetCubeCoord(*a, min_coords);
        // Check all cubes around the cube containing atom A.
        // std::cout << PrintPoint(atom_coords) << ":" << std::endl;
        for (int i = -1; i < 2; ++i) {
          for (int j = -1; j < 2; ++j) {
            for (int k = -1; k < 2; ++k) {
              Point coords = {
                atom_coords.x + i,
                atom_coords.y + j,
                atom_coords.z + k
              };
              // std::cout << "\t" << PrintPoint(coords) << std::endl;
              if (hashmap.find(coords) == hashmap.end()) {
                continue;
              } else {
                const auto b_atoms = hashmap[coords];
                for (auto b = b_atoms.begin(); b != b_atoms.end(); ++b) {
                  ++num_comparisons;
                  if (get_atoms_distance(*a, *b) < kCubeSize) {
                    ++num_clashes;
                    clashing_atoms.push_back(std::make_pair(*a, *b));
                    clashes[a->serial] = *a;
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
            clashing_atoms.push_back(std::make_pair(*a, *b));
            clashes[a->serial] = *a;
          }
        }
      }
    }
    for (auto i = clashing_atoms.begin(); i != clashing_atoms.end(); ++i) {
      const auto atom_a = i->first;
      const auto atom_b = i->second;
      // std::cout << atom_a.serial << " " << atom_a.resName << "\t"
      //           << atom_b.serial << " " << atom_b.resName
      //           << std::endl;
      // std::cout << atom_a.serial << " " << atom_a.resName << " "
      //           << atom_a.resSeq << " " << atom_a.atomName
      //           << std::endl;
    }
    // std::cout << "[INFO] Number of clashing atoms:   " << num_clashes << std::endl;
    std::cout << "[INFO] Number of clashing atoms:   " << clashes.size() << std::endl;
    std::cout << "[INFO] Number of comparisons made: " << num_comparisons
              << std::endl;
  }

  Protein(const char* pdb_file) {
    this->max_coords_ = {-1e308, -1e308, -1e308};
    this->min_coords_ = {1e308, 1e308, 1e308};
    const int num_atoms = read_data(pdb_file, &this->atom_callback,
      static_cast<void*>(&this->atoms_));
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      this->UpdateMaxCoordinates(a->centre);
      this->UpdateMinCoordinates(a->centre);
    }
    this->comparisons_cnt_ = 0;
  }

  ~Protein() {}

  void MapAtomsToCubes(const Point& min_coords) {
    // // TODO: Don't I actually need the min of the maxima and the max of the minima??
    // // Afterall, I can skip the parts where the two proteins do not overlap at
    // // all...
    // const double max_x = other.max_coords_.x; // std::min(this->max_coords_.x, other.max_coords_.x);
    // const double max_y = other.max_coords_.y; // std::min(this->max_coords_.y, other.max_coords_.y);
    // const double max_z = other.max_coords_.z; // std::min(this->max_coords_.z, other.max_coords_.z);
    // const double min_x = other.min_coords_.x; // std::max(this->min_coords_.x, other.min_coords_.x);
    // const double min_y = other.min_coords_.y; // std::max(this->min_coords_.y, other.min_coords_.y);
    // const double min_z = other.min_coords_.z; // std::max(this->min_coords_.z, other.min_coords_.z);
    // const double kCubeSize = Protein::kAtomRadius * 2;
    // this->num_cubes_.x = std::ceil((max_x - min_x) / kCubeSize);
    // this->num_cubes_.y = std::ceil((max_y - min_y) / kCubeSize);
    // this->num_cubes_.z = std::ceil((max_z - min_z) / kCubeSize);
    // // NOTE: Having multiple atoms within a cube with a protein is fine (the
    // // clashes arise when we compare ANOTHER protein).
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      Point cube_coord = Protein::GetCubeCoord(*a, min_coords);

      // Point cube_coord;
      // cube_coord.x = (a->centre.x - min_x) / kCubeSize;
      // cube_coord.y = (a->centre.y - min_y) / kCubeSize;
      // cube_coord.z = (a->centre.z - min_z) / kCubeSize;
      // cube_coord.x = (cube_coord.x > 0) ? std::ceil(cube_coord.x) : std::floor(cube_coord.x);
      // cube_coord.y = (cube_coord.y > 0) ? std::ceil(cube_coord.y) : std::floor(cube_coord.y);
      // cube_coord.z = (cube_coord.z > 0) ? std::ceil(cube_coord.z) : std::floor(cube_coord.z);
      // std::cout << "(" << a->centre.x << ", " << a->centre.y << ", " << a->centre.z << ") -> ("
      //           << cube_coord.x << ", " << cube_coord.y << ", " << cube_coord.z << ")"
      //           << std::endl;
      // int id = cube_coord.z * num_cubes_y * num_cubes_x + cube_coord.y * num_cubes_x + cube_coord.x;
      // if (this->ht_.find(id) != this->ht_.end()) {
      //   std::cout << "ERROR. Atom already in HT!" << std::endl;
      // }
      // this->ht_[id] = *a;

      // if (this->hash_table_.find(cube_coord) != this->hash_table_.end()) {
      //   std::cout << "ERROR. Atom already in hash table!" << std::endl;
      //   std::cout << "(" << cube_coord.x << ", " << cube_coord.y << ", " << cube_coord.z << ") "
      //             << "(" << a->centre.x << ", " << a->centre.y << ", " << a->centre.z << ")"
      //             << "(" << this->hash_table_[cube_coord].centre.x << ", " << this->hash_table_[cube_coord].centre.y << ", " << this->hash_table_[cube_coord].centre.z << ")"
      //             << std::endl;
      // }
      this->hash_table_[cube_coord] = *a;
      this->ht_[cube_coord].push_back(*a);
    }
  }

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

  bool IsAtomClashing(const Point& atom_centre, const Point& min, const Point& max) {
    const double kCubeSize = Protein::kAtomRadius * 2;
    this->comparisons_cnt_++;
    // // Check if atom from protein A is within min-max of protein B.
    // Point tmp = {
    //   (1 + max.x - min.x) * kCubeSize,
    //   (1 + max.y - min.y) * kCubeSize,
    //   (1 + max.z - min.z) * kCubeSize
    // };
    // if (atom_centre.x < min.x || atom_centre.y < min.y || atom_centre.z < min.z) {
    //   return false;
    // }
    // if (atom_centre.x > tmp.x || atom_centre.y > tmp.y || atom_centre.z > tmp.z) {
    //   return false;
    // }
    if (this->hash_table_.find(atom_centre) != this->hash_table_.end()) {
      return true;
    } else {
      return false;
    }
  }

  Atom GetAtomFromHTable(const Point& coords) {
    return this->hash_table_[coords];
  }

  void PrintAtoms() {
    for (auto a = this->atoms_.begin(); a != this->atoms_.end(); ++a) {
      print_pdb_atom(a->serial, a->atomName, a->altLoc, a->resName, a->chainID,
        a->resSeq, a->iCode, a->centre);
    }
  }

  void ResetComparisonCnt() {
    this->comparisons_cnt_ = 0;
  }

  // std::unordered_map<Point, Atom, PointHashFunc, PointEqualsFunc> get_hash_table() {
  //   return this->hash_table_;
  // }
  std::unordered_map<Point, std::vector<Atom>, PointHashFunc, PointEqualsFunc> get_hash_table() {
    return this->ht_;
  }

  int get_comparisons_cnt() {
    return this->comparisons_cnt_;
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
  std::unordered_map<Point, Atom, PointHashFunc, PointEqualsFunc> hash_table_;
  std::unordered_map<Point, std::vector<Atom>, PointHashFunc, PointEqualsFunc> ht_;
  Point num_cubes_;
  int comparisons_cnt_;
};


int main(int argc, char const *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "ERROR. Usage: detect_steric_clashes.exe file1.pdb file2.pdb\n");
    exit(1);
  }
  Protein a = Protein(argv[1]);
  Protein b = Protein(argv[2]);
  Protein::DetectStericOverlaps(a, b, true);
  Protein::DetectStericOverlaps(b, a, true);
  Protein::DetectStericOverlaps(a, b, false);
  Protein::DetectStericOverlaps(b, a, false);

  // a.MapAtomsToCubes(b);
  // b.MapAtomsToCubes(a);
  // a.PrintAtoms();
  // b.PrintAtoms();
  return 0;
}