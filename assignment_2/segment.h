/*
 * File:  segment.h
 * Author: Stefano Ribes
 */
#ifndef SEGMENT_H_
#define SEGMENT_H_

#include "atom.h"
#include "residue.h"

#define DOMAK_MAX_NUM_DOMAINS 40
#define DOMAK_MAX_SEGMENTS_PER_DOMAIN 2
#define DOMAK_MDS 40 // Minimum domain size
#define DOMAK_MNCCm 30 // Minimum no contact cut-off middle of chain
#define DOMAK_MNCCe 10 // Minimum no contact cut-off end of chain
#define DOMAK_MSSm 25 // Minimum segment size middle of chain
#define DOMAK_MSSe 5 // Minimum segment size end of chain
#define DOMAK_sso 0.57 // If percentage of secondary structure is greater than
                     // this only use secondary structure contacts
#define DOMAK_MSV 9.5 // Minimum split value
#define DOMAK_MSVsso 17.05 // Minimum split value using only secondary structure
                        // contacts
#define DOMAK_MSVcs 60.0 // Minimum split value for chopped segments
#define DOMAK_MDSP 120 // Minimum size of segment for a double split
#define DOMAK_ID 250 // Increment divider

typedef struct {
  int start;
  int end;
  int num_internal_contacts;
  double max_split_val;
} Segment;

static const Segment null_segment = {0, 0, 0};

typedef struct {
  Segment segments[DOMAK_MAX_SEGMENTS_PER_DOMAIN];
  int num_segments;
} Domain;

int len(const Segment s);

void set_int_cnt_for_atom(const double dist_threshold, const int num_residues,
    const Residue* residues, const char* atom_name, Segment* segment,
    double* dist_lookup);

void set_int_cnt(const double dist_threshold, const int num_residues,
    const Residue* residues, Segment* segment, double* dist_lookup);

int get_ext_cnt_for_atom(const double dist_threshold, const int num_residues,
    const Residue* residues, const char* atom_name, const Segment a,
    const Segment b, double* dist_lookup);

int get_ext_cnt(const double dist_threshold, const int num_residues,
    const Residue* residues, const Segment a, const Segment b,
    double* dist_lookup);

void single_segment_scan(const double dist_threshold, const int num_residues,
    const Residue* residues, const int curr_domain_idx, Domain* domains,
    int* num_domains, double* dist_lookup);

void two_segment_scan_of_two_segment_domain(const double dist_threshold,
    const int num_residues, const Residue* residues, const int curr_domain_idx,
    Domain* domains, int* num_domains, double* dist_lookup);

#endif // end SEGMENT_H_