/*
 * @author     Stefano Ribes
 *
 * @brief      Global alignment.
 *
 * @details    To compile this C program, type:
 *
 *             gcc -O3 -std=c99 global_alignment.c -o global_alignment.exe
 *
 *             To run the program, type:
 *
 *             ./global_alignment.exe
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define MAX_LENGTH 100
#define MAX_PATHS 100

#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_PENALTY 2

#define STOP 0
#define UP 1
#define LEFT 2
#define DIAG 3

/*
 * @brief      Data structure to store the up, left and diagonal scores of a
 *             cell.
 */
typedef struct {
  int up;
  int diag;
  int left;
} score_t;

/*
 * @brief      Index data structure to store cells coordinates.
 */
typedef struct {
  int x;
  int y;
} idx_t;

/*
 * @brief      Data structure holding information for performing backtracking
 *             search of the optimal paths.
 */
typedef struct {
  char* X;
  char* Y;
  int m;
  int n;
  idx_t dst;
  idx_t* path;
  bool* visited_cells;
  int path_index;
  int num_optimal_paths;
} search_t;

/**
 * @brief      Prints an alignment sequence.
 *
 * @param[in]  alignment_length  The alignment length
 * @param[in]  alignX            The aligned string of the x sequence
 * @param[in]  alignY            The aligned string of the y sequence
 */
void print_alignment(const int alignment_length,
    const char alignX[MAX_LENGTH * 2], const char alignY[MAX_LENGTH * 2]);

/**
 * @brief      Prints all optimal paths given two sequences.
 *
 * @param[in]  X       The X input sequence
 * @param[in]  Y       The Y input sequence
 * @param[in]  m       The length of the X sequence
 * @param[in]  n       The length of the Y sequence
 * @param[in]  trace   The trace matrix
 * @param[in]  scores  The scores matrix
 * @param[in]  src     The starting cell from which starting the alignment
 * @param[in]  dst     The destination cell to end the alignment
 */
void print_all_paths(
    const char* X, //[MAX_LENGTH + 1],
    const char* Y, //[MAX_LENGTH + 1],
    const int m,
    const int n,
    const int trace[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const score_t scores[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const idx_t src,
    const idx_t dst);

/**
 * @brief      Prints all paths utility (recursive function).
 *
 * @param[in]  trace      The trace
 * @param[in]  scores     The scores
 * @param[in]  curr_cell  The current cell to visit
 * @param      search     The search data structure
 */
void print_all_paths_util(
  const int trace[MAX_LENGTH + 1][MAX_LENGTH + 1],
  const score_t scores[MAX_LENGTH + 1][MAX_LENGTH + 1],
  idx_t curr_cell,
  search_t* search);

int main(int argc, char** argv) {
  int i, j;
  int m, n;
  int alignment_length, score, tmp;
  char X[MAX_LENGTH + 1] = "ATCGAT"; // "ATTA";
  char Y[MAX_LENGTH + 1] = "ATACGT"; // "ATTTTA";

  if (argc == 3) {
    strcpy(X, argv[1]);
    strcpy(Y, argv[2]);
  }

  int F[MAX_LENGTH + 1][MAX_LENGTH + 1]; // Score matrix
  int trace[MAX_LENGTH + 1][MAX_LENGTH + 1]; // Trace matrix
  char alignX[MAX_LENGTH * 2]; // Aligned X sequence
  char alignY[MAX_LENGTH * 2]; // Aligned Y sequence

  int max_score = -((1 << 30) - 1); // Fairly small number
  const score_t no_match = {0};
  score_t scores[MAX_LENGTH + 1][MAX_LENGTH + 1] = {no_match}; // Scores matrix
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
      scores[i][j].diag = score;
      scores[i][j].up = tmp;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = UP;
      }
      tmp = F[i][j-1] - GAP_PENALTY;
      scores[i][j].left = tmp;
      if(tmp > score) {
        score = tmp;
        trace[i][j] = LEFT;
      }
      F[i][j] = score;
      max_score = score > max_score ? score : max_score;
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
        break;
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
  print_alignment(alignment_length, alignX, alignY);
  idx_t start, end;
  start.x = m;
  start.y = n;
  end.x = 0;
  end.y = 0;
  print_all_paths(X, Y, m, n, trace, scores, start, end);
  return 0;
}

void print_alignment(const int alignment_length,
    const char alignX[MAX_LENGTH * 2],
    const char alignY[MAX_LENGTH * 2]) {
  int match_cnt = 0;
  printf("* Alignment Sequence:\n");
  for (int i = alignment_length - 1; i >= 0; --i) {
    printf("%c", alignX[i]);
  }
  printf("\n");
  for (int i = alignment_length - 1; i >= 0; --i) {
    if (alignX[i] == alignY[i]) {
      printf("|");
      ++match_cnt;
    } else {
      printf(" ");
    }
  }
  printf("\n");
  for (int i = alignment_length - 1; i >= 0; --i) {
    printf("%c", alignY[i]);
  }
  printf("\n");
  // NOTE: The percentage identity is calculated based on the total alignment
  // length, i.e. including all the indel.
  const float perc_identity = (float)match_cnt / (float)alignment_length * 100.;
  printf("[INFO] Percent identity: %.2f%\n", perc_identity);
  // NOTE: The Hamming distance only works for sequences of same length. We
  // calculate the Hamming distance by getting the difference between the
  // aligned sequences and the number of matches.
  printf("[INFO] Hamming distance: %d\n", alignment_length - match_cnt);
}

void print_all_paths(
    const char X[MAX_LENGTH + 1],
    const char Y[MAX_LENGTH + 1],
    const int m,
    const int n,
    const int trace[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const score_t scores[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const idx_t src,
    const idx_t dst) {
  // Setup search data structure
  const int kNumCells = (m + 1) * (n + 1);
  search_t* search = malloc(sizeof(search_t));
  search->visited_cells = malloc(sizeof(bool) * kNumCells);
  search->path = malloc(sizeof(idx_t) * kNumCells);
  search->path_index = 0;
  search->num_optimal_paths = 0;
  search->dst = dst;
  search->X = malloc(sizeof(m+1));
  search->Y = malloc(sizeof(n+1));
  search->m = m;
  search->n = n;
  memcpy(search->X, X, m+1);
  memcpy(search->Y, Y, n+1);
  for (int i = 0; i < kNumCells; ++i) {
    search->visited_cells[i] = false;
  }
  printf("[INFO] All optimal paths:\n");
  print_all_paths_util(trace, scores, src, search);
  printf("[INFO] Number of optimal paths found: %d\n", search->num_optimal_paths);
  free(search->X);
  free(search->Y);
  free(search->visited_cells);
  free(search->path);
  free(search);
}

int get_cell_id(const idx_t p, const int size) {
  return p.x * size + p.y;
}

void print_all_paths_util(
    const int trace[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const score_t scores[MAX_LENGTH + 1][MAX_LENGTH + 1],
    const idx_t c,
    search_t* search) {
  // If cell with negative coordinates, end recursion.
  if (c.x < 0 || c.y < 0) {
    return;
  }
  // Update visited cells and append current cell c to path.
  search->visited_cells[get_cell_id(c, search->n)] = true;
  search->path[search->path_index] = c;
  search->path_index++;
  const idx_t* path = search->path;
  const bool* visited_cells = search->visited_cells;
  // If current cell is same as destination or STOP reached, then print path
  if ((c.x == search->dst.x && c.y == search->dst.y) || trace[c.x][c.y] == STOP) {
    const char* X = search->X;
    const char* Y = search->Y;
    char alignX[MAX_LENGTH * 2];
    char alignY[MAX_LENGTH * 2];
    int alignment_length = 0;
    int x_idx = path[0].x;
    int y_idx = path[0].y;
    // In order to trace back the path, compare the COORDINATES of the cells in
    // the path.
    for (int i = 0; i < search->path_index - 1; ++i) {
      if (path[i].x-1 == path[i+1].x && path[i].y-1 == path[i+1].y) {
        // diag
        alignX[alignment_length] = X[path[i].x-1];
        alignY[alignment_length] = Y[path[i].y-1];
        --x_idx;
        --y_idx;
      } else if (path[i].x-1 == path[i+1].x && path[i].y == path[i+1].y) {
        // up
        alignX[alignment_length] = X[path[i].x-1];
        alignY[alignment_length] = '-';
        --x_idx;
      } else if (path[i].x == path[i+1].x && path[i].y-1 == path[i+1].y) {
        // left
        alignX[alignment_length] = '-';
        alignY[alignment_length] = Y[path[i].y-1];
        --y_idx;
      }
      ++alignment_length;
    }
    // Unaligned beginning
    while (x_idx > 0) {
      alignX[alignment_length] = X[x_idx-1];
      alignY[alignment_length] = '-';
      --x_idx;
      ++alignment_length;
    }
    while (y_idx > 0) {
      alignX[alignment_length] = '-';
      alignY[alignment_length] = Y[y_idx-1];
      --y_idx;
      ++alignment_length;
    }
    search->num_optimal_paths++;
    printf("* Path n.%d: List of coordinates: ", search->num_optimal_paths);
    for (int i = 0; i < search->path_index; ++i) {
      printf("(%d, %d) ", path[i].x, path[i].y);
    }
    printf("\n");
    print_alignment(alignment_length, alignX, alignY);
    printf("\n");
  } else {
    // If current cell is not the final destination, call
    // recursion/backtracking/visit on cells with same score directions and
    // which haven't been visited yet.
    idx_t up, diag, left;
    up.x = c.x - 1;
    up.y = c.y;
    diag.x = c.x - 1;
    diag.y = c.y - 1;
    left.x = c.x;
    left.y = c.y - 1;
    switch (trace[c.x][c.y]) {
      case UP:
        if (!visited_cells[get_cell_id(up, search->n)]) {
          print_all_paths_util(trace, scores, up, search); // Visit Up cell.
        }
        if (scores[c.x][c.y].up == scores[c.x][c.y].diag &&
            !visited_cells[get_cell_id(diag, search->n)]) {
          print_all_paths_util(trace, scores, diag, search); // Visit Diag cell.
        }
        if (scores[c.x][c.y].up == scores[c.x][c.y].left &&
            !visited_cells[get_cell_id(left, search->n)]) {
          print_all_paths_util(trace, scores, left, search); // Visit Left cell.
        }
        break;
      case LEFT:
        if (!visited_cells[get_cell_id(left, search->n)]) {
          print_all_paths_util(trace, scores, left, search); // Visit Left cell.
        }
        if (scores[c.x][c.y].left == scores[c.x][c.y].diag &&
            !visited_cells[get_cell_id(diag, search->n)]) {
          print_all_paths_util(trace, scores, diag, search); // Visit Diag cell.
        }
        if (scores[c.x][c.y].left == scores[c.x][c.y].up &&
            !visited_cells[get_cell_id(up, search->n)]) {
          print_all_paths_util(trace, scores, up, search); // Visit Up cell.
        }
        break;
      case DIAG:
        if (!visited_cells[get_cell_id(diag, search->n)]) {
          print_all_paths_util(trace, scores, diag, search); // Visit Diag cell.
        }
        if (scores[c.x][c.y].diag == scores[c.x][c.y].up &&
            !visited_cells[get_cell_id(up, search->n)]) {
          print_all_paths_util(trace, scores, up, search); // Visit Up cell.
        }
        if (scores[c.x][c.y].diag == scores[c.x][c.y].left &&
            !visited_cells[get_cell_id(left, search->n)]) {
          print_all_paths_util(trace, scores, left, search); // Visit Left cell.
        }
        break;
    }
  }
  // Remove current cell from path and set it as not visited
  search->path_index--;
  search->visited_cells[get_cell_id(c, search->n)] = false;
}