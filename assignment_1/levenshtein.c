/*
 * @author     Stefano Ribes
 *
 * @brief      Global alignment with levenshtein function.
 *
 * @details    To compile this C program, type:
 *
 *             gcc -O3 -std=c99 levenshtein.c -o levenshtein.exe
 *
 *             To run the program, type:
 *
 *             ./levenshtein.exe
 */
#include <stdio.h>

#define MAX_LENGTH 100

#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_PENALTY 2

#define STOP 0
#define UP 1
#define LEFT 2
#define DIAG 3

/**
 * @brief      Calculate Levenshtein distance between two sequences.
 *
 * @param[in]  length_a  The length of sequence a
 * @param[in]  length_b  The length of sequence b
 * @param[in]  start_a   The start index of a
 * @param[in]  start_b   The start index of b
 * @param      a         The a sequence
 * @param      b         The b sequence
 *
 * @return     The Levenshtein distance
 */
int levenshtein(const int length_a, const int length_b, const int start_a,
    const int start_b, char* const a, char* const b) {
  // NOTE: Increasing the start index effectively translates to taking the tail
  // of the string.
  if (start_b == length_b) {
    return length_a - start_a;
  }
  if (start_a == length_a) {
    return length_b - start_b;
  }
  if (a[start_a] == b[start_b]) {
    return levenshtein(length_a, length_b, start_a + 1, start_b + 1, a, b);
  }
  int min = levenshtein(length_a, length_b, start_a + 1, start_b, a, b);
  int tmp = levenshtein(length_a, length_b, start_a, start_b + 1, a, b);
  min = (tmp < min) ? tmp : min;
  tmp = levenshtein(length_a, length_b, start_a + 1, start_b + 1, a, b);
  min = (tmp < min) ? tmp : min;
  return 1 + min;
}

int main() {
  int i, j;
  int m, n;
  int alignment_length, score, tmp;
  char X[MAX_LENGTH+1] = "ATCGAT";
  char Y[MAX_LENGTH+1] = "ATACGT";

  int F[MAX_LENGTH+1][MAX_LENGTH+1]; // score matrix
  int trace[MAX_LENGTH+1][MAX_LENGTH+1]; // trace matrix
  char alignX[MAX_LENGTH*2]; // aligned X sequence
  char alignY[MAX_LENGTH*2]; // aligned Y sequence
  /*
   * Find lengths of (null-terminated) strings X and Y
   */
  m = 0;
  n = 0;
  while (X[m] != 0) {
    m++;
  }
  while (Y[n] != 0) {
    n++;
  }
  /*
   * Initialise matrices
   */
  F[0][0] = 0;
  trace[0][0] = STOP;
  for (int i = 1; i <= m; ++i) {
    F[i][0] = F[i-1][0] - GAP_PENALTY;
    trace[i][0] = STOP;
  }
  for (int j = 1; j <= n; j++) {
    F[0][j] = F[0][j-1] - GAP_PENALTY;
    trace[0][j] = STOP;
  }
  /*
   * Fill matrices
   */
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; j++) {
      if (X[i-1] == Y[j-1]) {
        score = F[i-1][j-1] + MATCH_SCORE;
      } else {
        score = F[i-1][j-1] + MISMATCH_SCORE;
      }
      trace[i][j] = DIAG;
      tmp = F[i-1][j] - GAP_PENALTY;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = UP;
      }
      tmp = F[i][j-1] - GAP_PENALTY;
      if(tmp > score) {
        score = tmp;
        trace[i][j] = LEFT;
      }
      F[i][j] = score;
     }
  }
  /*
   * Print score matrix
   */
#define PRINT_MATRIX(m, n, x, y, f) printf("      "); \
  for (int j = 0; j < n; ++j) { \
    printf("%5c", y[j]); \
  } \
  printf("\n"); \
  for (int i = 0; i <= m ; ++i) { \
    if (i == 0) { \
      printf(" "); \
    } else { \
      printf("%c", x[i-1]); \
    } \
    for (int j = 0; j <= n; j++) { \
      printf("%5d", f[i][j]); \
    } \
    printf("\n"); \
  } \
  printf("\n");

  printf("[INFO] Score matrix:\n");
  PRINT_MATRIX(m, n, X, Y, F);
  printf("[INFO] Trace matrix:\n");
  PRINT_MATRIX(m, n, X, Y, trace);
  /*
   * Trace back from the lower-right corner of the matrix
   */
  i = m;
  j = n;
  alignment_length = 0;
  while (trace[i][j] != STOP) {
    switch (trace[i][j]) {
      case DIAG:
        alignX[alignment_length] = X[i-1];
        alignY[alignment_length] = Y[j-1];
        --i;
        --j;
        ++alignment_length;
        break;
      case LEFT:
        alignX[alignment_length] = '-';
        alignY[alignment_length] = Y[j-1];
        --j;
        ++alignment_length;
        break;
      case UP:
        alignX[alignment_length] = X[i-1];
        alignY[alignment_length] = '-';
        --i;
        ++alignment_length;
    }
  }
  /*
   * Unaligned beginning
   */
  while (i > 0) {
    alignX[alignment_length] = X[i-1];
    alignY[alignment_length] = '-';
    --i;
    ++alignment_length;
  }
  while (j > 0) {
    alignX[alignment_length] = '-';
    alignY[alignment_length] = Y[j-1];
    --j;
    ++alignment_length;
  }
  /*
   * Print alignment
   */
  for (i = alignment_length - 1; i >= 0; --i) {
    printf("%c", alignX[i]);
  }
  printf("\n");
  int match_cnt = 0;
  for (i = alignment_length - 1; i >= 0; --i) {
    if (alignX[i] == alignY[i]) {
      printf("|");
      ++match_cnt;
    } else {
      printf(" ");
    }
  }
  printf("\n");
  for (i = alignment_length - 1; i >= 0; --i) {
    printf("%c", alignY[i]);
  }
  printf("\n");
  const float perc_identity = (float)match_cnt / (float)alignment_length * 100.;
  const int lev_dist = levenshtein(alignment_length, alignment_length, 0, 0,
                                   alignX, alignY);
  printf("[INFO] Percent identity: %.2f%\n", perc_identity);
  printf("[INFO] Levenshtein dist: %d\n", lev_dist);
  return 0;
}
