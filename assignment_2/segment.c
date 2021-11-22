#include "segment.h"

int len(const Segment s) {
  return s.end - s.start;
}

void set_int_cnt_for_atom(const double dist_threshold, const int num_residues,
    const Residue* residues, const char* atom_name, Segment* segment,
    double* dist_lookup) {
  double dist;
  Atom a, b;
  for (int i = segment->start; i <= segment->end; ++i) {
    if (get_atom_from_residue(residues[i], atom_name, &a)) {
      for (int j = segment->start; j <= segment->end; ++j) {
        if (get_atom_from_residue(residues[j], atom_name, &b)) {
          // Check if the distance is already in the lookup table (distances
          // cannot be negative). Otherwise, compute it (the lookup matrix is
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

/**
 * @brief      Sets the interior count. Only considering the first atom in the
 *             residue, i.e. only the heavy atom is stored in the residue.
 *
 * @param[in]  dist_threshold  The distance threshold
 * @param[in]  num_residues    The number of residues
 * @param[in]  residues        The residues
 * @param      segment         The segment to update
 * @param      dist_lookup     The distance lookup table
 */
void set_int_cnt(const double dist_threshold, const int num_residues,
    const Residue* residues, Segment* segment, double* dist_lookup) {
  double dist;
  Atom a, b;
  for (int i = segment->start; i <= segment->end; ++i) {
    a = residues[i].atom[0];
    for (int j = segment->start; j <= segment->end; ++j) {
      b = residues[j].atom[0];
      // Check if the distance is already in the lookup table (distances
      // cannot be negative). Otherwise, compute it (the lookup matrix is
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

int get_ext_cnt_for_atom(const double dist_threshold, const int num_residues,
    const Residue* residues, const char* atom_name, const Segment a,
    const Segment b, double* dist_lookup) {
  double dist;
  Atom x, y;
  int num_exterior_contacts = 0;
  for (int i = a.start; i <= a.end; ++i) {
    if (get_atom_from_residue(residues[i], atom_name, &x)) {
      for (int j = b.start; j <= b.end; ++j) {
        if (get_atom_from_residue(residues[j], atom_name, &y)) {
          // Check if the distance is already in the lookup table (distances
          // cannot be negative). Otherwise, compute it (the lookup matrix is
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

int get_ext_cnt(const double dist_threshold, const int num_residues,
    const Residue* residues, const Segment a, const Segment b,
    double* dist_lookup) {
  double dist;
  Atom x, y;
  int num_exterior_contacts = 0;
  for (int i = a.start; i <= a.end; ++i) {
    x = residues[i].atom[0];
    for (int j = b.start; j <= b.end; ++j) {
      y = residues[j].atom[0];
      // Check if the distance is already in the lookup table (distances
      // cannot be negative). Otherwise, compute it (the lookup matrix is
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
  return num_exterior_contacts;
}

double get_split_val(const double dist_threshold, const int num_residues,
    const Residue* residues, Segment* a, Segment* b, double* dist_lookup) {
  if (len(*a) <= DOMAK_MDS || len(*b) <= DOMAK_MDS) {
    return 0;
  }
  set_int_cnt(dist_threshold, num_residues, residues, a, dist_lookup);
  set_int_cnt(dist_threshold, num_residues, residues, b, dist_lookup);
  const double int_a = (double)a->num_internal_contacts;
  const double int_b = (double)b->num_internal_contacts;
  const double ext_ab = (double)get_ext_cnt(dist_threshold, num_residues,
    residues, *a, *b, dist_lookup);
  return (int_a / ext_ab) * (int_b / ext_ab);
}

void single_segment_scan(const double dist_threshold, const int num_residues,
    const Residue* residues, const int curr_domain_idx, Domain* domains,
    int* num_domains, double* dist_lookup) {
  Segment a_max, a, b1, b2;
  Segment b = domains[curr_domain_idx].segments[0];
  double max_split_val = 0;
  // TODO: Change for indeces to account for the DOMAK_MDS constraint.
  for (int i = b.start; i <= b.end; ++i) {
    for (int j = i; j <= b.end; ++j) {
      a.start = i;
      a.end = j;
      b1.start = b.start;
      b1.end = i;
      b2.start = j;
      b2.end = b.end;
      const double split_b1 = get_split_val(dist_threshold, num_residues,
        residues, &a, &b1, dist_lookup);
      const double split_b2 = get_split_val(dist_threshold, num_residues,
        residues, &a, &b2, dist_lookup);
      if (split_b1 > max_split_val) {
        max_split_val = split_b1;
        a_max = a;
      }
      if (split_b2 > max_split_val) {
        max_split_val = split_b2;
        a_max = a;
      }
    }
  }
  if (max_split_val > DOMAK_MSV) {
    // A_max and B_max not correlated: extract A as a new domain and update the
    // current domain (eventually with B1 and B2).
    ++(*num_domains);
    domains[*num_domains].segments[0] = a_max;
    domains[*num_domains].segments[1] = null_segment;
    domains[*num_domains].num_segments = 1;
    if (a_max.end == b.end) {
      domains[curr_domain_idx].segments[0].end = a_max.start;
      domains[curr_domain_idx].segments[1] = null_segment;
      domains[curr_domain_idx].num_segments = 1;
    } else if (a_max.start == b.start) {
      domains[curr_domain_idx].segments[0].start = a_max.start;
      domains[curr_domain_idx].segments[1] = null_segment;
      domains[curr_domain_idx].num_segments = 1;
    } else {
      b1.start = b.start;
      b1.end = a_max.start;
      b2.start = a_max.end;
      b2.end = b.end;
      const double split_b = get_split_val(dist_threshold, num_residues,
        residues, &b1, &b2, dist_lookup);
      if (max_split_val > DOMAK_MSV) {
        // B1 and B2 not correlated: generate another domain
        domains[curr_domain_idx].segments[0] = b1;
        domains[curr_domain_idx].segments[1] = null_segment;
        domains[curr_domain_idx].num_segments = 1;
        ++(*num_domains);
        domains[*num_domains].segments[0] = b2;
        domains[*num_domains].segments[1] = null_segment;
        domains[*num_domains].num_segments = 1;
      } else {
        domains[curr_domain_idx].segments[0] = b1;
        domains[curr_domain_idx].segments[1] = b2;
        domains[curr_domain_idx].num_segments = 2;
      }
    }
  }
}

void two_segment_scan_of_two_segment_domain(const double dist_threshold,
    const int num_residues, const Residue* residues, const int curr_domain_idx,
    Domain* domains, int* num_domains, double* dist_lookup) {
  double max_split_a = 0;
  double max_split_b = 0;
  double max_split_val = 0;
  Segment a1, a2, b1, b2, a1_max, a2_max, b1_max, b2_max;
  Segment seg_lo = domains[*num_domains].segments[0];
  Segment seg_hi = domains[*num_domains].segments[1];
  a1.end = seg_lo.end;
  a2.start = seg_hi.start;
  b1.start = seg_lo.start;
  b2.end = seg_hi.end;
  for (int i = seg_lo.start; i <= seg_lo.end; ++i) {
    for (int j = seg_hi.start; j <= seg_hi.end; ++j) {
      a1.start = i;
      a2.end = j;
      b1.end = i;
      b2.start = j;
      const double split_a1b1 = get_split_val(dist_threshold, num_residues,
        residues, &a1, &b1, dist_lookup);
      const double split_a1b2 = get_split_val(dist_threshold, num_residues,
        residues, &a1, &b2, dist_lookup);
      const double split_a2b1 = get_split_val(dist_threshold, num_residues,
        residues, &a2, &b1, dist_lookup);
      const double split_a2b2 = get_split_val(dist_threshold, num_residues,
        residues, &a2, &b2, dist_lookup);
      if (split_a1b1 > max_split_val) {
        max_split_val = split_a1b1;
      } else if (split_a1b2 > max_split_val) {
        max_split_val = split_a1b2;
      } else if (split_a2b1 > max_split_val) {
        max_split_val = split_a2b1;
      } else if (split_a2b2 > max_split_val) {
        max_split_val = split_a2b2;
      }
      if (split_a1b1 > max_split_val || split_a1b2 > max_split_val ||
          split_a2b1 > max_split_val || split_a2b2 > max_split_val) {
        a1_max = a1;
        a2_max = a2;
        b1_max = b1;
        b2_max = b2;
      }
    }
  }
  if (max_split_val > DOMAK_MSV) {
    // A_max and B_max not correlated: extract A as a new two-segment domain and
    // update the current domain (eventually with B1 and B2).
    ++(*num_domains);
    domains[*num_domains].segments[0] = a1_max;
    domains[*num_domains].segments[1] = a2_max;
    domains[*num_domains].num_segments = 2;
    if (a1_max.start == seg_lo.start) {
      domains[curr_domain_idx].segments[0] = b2_max;
      domains[curr_domain_idx].segments[1] = null_segment;
      domains[curr_domain_idx].num_segments = 1;
    } else if (a2_max.end == seg_hi.end) {
      domains[curr_domain_idx].segments[0] = b1_max;
      domains[curr_domain_idx].segments[1] = null_segment;
      domains[curr_domain_idx].num_segments = 1;
    } else {
      const double split_b = get_split_val(dist_threshold, num_residues,
        residues, &b1_max, &b2_max, dist_lookup);
      if (max_split_val > DOMAK_MSV) {
        // B1 and B2 not correlated: generate another domain
        domains[curr_domain_idx].segments[0] = b1_max;
        domains[curr_domain_idx].segments[1] = null_segment;
        domains[curr_domain_idx].num_segments = 1;
        ++(*num_domains);
        domains[*num_domains].segments[0] = b2_max;
        domains[*num_domains].segments[1] = null_segment;
        domains[*num_domains].num_segments = 1;
      } else {
        domains[curr_domain_idx].segments[0] = b1_max;
        domains[curr_domain_idx].segments[1] = b2_max;
        domains[curr_domain_idx].num_segments = 2;
      }
    }
  }
}