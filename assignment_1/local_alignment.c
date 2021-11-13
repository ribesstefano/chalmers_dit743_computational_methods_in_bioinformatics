/*
 * @author     Stefano Ribes
 *
 * @brief      Local alignment.
 *
 * @details    To compile this C program, type:
 *
 *             gcc -O3 -std=c99 local_alignment.c -o local_alignment.exe
 *
 *             To run the program, type:
 *
 *             ./local_alignment.exe
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

int main(int argc, char** argv) {
  int i, j;
  int m, n;
  int alignment_length, score;
  char X[MAX_LENGTH + 1] = "PAWHEAE";
  char Y[MAX_LENGTH + 1] = "HDAGAWGHEQ";

  int F[MAX_LENGTH + 1][MAX_LENGTH + 1]; // score matrix
  int trace[MAX_LENGTH + 1][MAX_LENGTH + 1]; // trace matrix
  char alignX[MAX_LENGTH * 2]; // aligned X sequence
  char alignY[MAX_LENGTH * 2]; // aligned Y sequence
  /*
   * Find lengths of (null-terminated) strings X and Y
   */
  m = 0;
  n = 0;
  while (X[m] != 0) {
    ++m;
  }
  while (Y[n] != 0) {
    ++n;
  }
  /*
   * Initialise matrices
   */
  F[0][0] = 0;
  trace[0][0] = STOP;
  for (int i = 1; i <= m; ++i) {
    F[i][0] = F[i-1][0];
    trace[i][0] = STOP;
  }
  for (int j = 1; j <= n; j++) {
    F[0][j] = F[0][j-1];
    trace[0][j] = STOP;
  }

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

  /*
   * Fill matrices
   */
  int max_i = 0;
  int max_j = 0;
  int max_score = -((1 << 30) - 1); // Fairly small number
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; j++) {
      if (X[i-1] == Y[j-1]) {
        score = F[i-1][j-1] + MATCH_SCORE;
      } else {
        score = F[i-1][j-1] + MISMATCH_SCORE;
      }
      trace[i][j] = DIAG;
      int tmp = F[i-1][j] - GAP_PENALTY;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = UP;
      }
      tmp = F[i][j-1] - GAP_PENALTY;
      if(tmp > score) {
        score = tmp;
        trace[i][j] = LEFT;
      }
      score = score > 0 ? score : 0;
      F[i][j] = score;
      // Keep track of highest score.
      if (score > max_score) {
        max_score = score;
        max_i = i;
        max_j = j;
      }
      // Stop trace if score is equal to zero.
      if (score == 0) {
        trace[i][j] = STOP;
      }
    }
  }
  /*
   * Print score matrix
   */
  printf("[INFO] Best score %d at position (%d, %d)\n", max_score, max_i, max_j);
  printf("[INFO] Score matrix:\n");
  PRINT_MATRIX(m, n, X, Y, F);
  printf("[INFO] Trace matrix:\n");
  PRINT_MATRIX(m, n, X, Y, trace);
  /*
   * Trace back from the maximum score coordinates of the matrix
   */
  i = max_i;
  j = max_j;
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
  printf("\n\n");
  printf("[INFO] Percent identity: %.2f%\n", (float)match_cnt / (float)alignment_length * 100.0);
  printf("[INFO] Hamming distance: %d\n", alignment_length - match_cnt);
  return 0;
}
