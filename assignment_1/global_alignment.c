/*
 * To compile this C program, placing the executable file in 'global', type:
 *
 *      gcc -o global global_alignment.c
 *
 * To run the program, type:
 *
 *      ./global
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_LENGTH 100
#define MAX_PATHS 100

#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_PENALTY 2

#define STOP 0
#define UP 1
#define LEFT 2
#define DIAG 3

typedef struct {
  int up;
  int diag;
  int left;
} match_t;

typedef struct {
  int x;
  int y;
} point_t;

int get_cell_id(const point_t p, const int size) {
  return p.x * size + p.y;
}

int get_max_direction(match_t match) {
  int max = UP;
  if (match.left > max) {
    max = LEFT;
  }
  if (match.diag > max) {
    max = DIAG;
  }
  return max;
}

void print_alignment(const int alignment_length,
    const char alignX[MAX_LENGTH*2],
    const char alignY[MAX_LENGTH*2]) {
  int match_cnt = 0;
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
  printf("[INFO] Percent identity: %.2f%\n", (float)match_cnt / (float)alignment_length * 100.0);
  printf("[INFO] Hamming distance: %d\n", alignment_length - match_cnt); // Only for sequences of same length 
}

void print_all_paths_util(
  const char X[MAX_LENGTH+1],
  const char Y[MAX_LENGTH+1],
  const int n, const int trace[MAX_LENGTH+1][MAX_LENGTH+1], const match_t scores[MAX_LENGTH+1][MAX_LENGTH+1], point_t u, point_t d, bool* visited, point_t* path, int* path_index, int* num_optimal_paths);

void print_all_paths(
    const char X[MAX_LENGTH+1],
    const char Y[MAX_LENGTH+1],
    const int m, const int n, const int trace[MAX_LENGTH+1][MAX_LENGTH+1], match_t scores[MAX_LENGTH+1][MAX_LENGTH+1], const point_t s, const point_t d) {
  const int kNumCells = (m + 1) * (n + 1);
  bool* visited = malloc(sizeof(bool) * kNumCells);
  point_t* path = malloc(sizeof(point_t) * kNumCells);
  int path_index = 0;
  int num_optimal_paths = 0;
  for (int i = 0; i < kNumCells; ++i) {
    visited[i] = false;
  }
  printf("[INFO] All optimal paths:\n");
  print_all_paths_util(X, Y, n, trace, scores, s, d, visited, path, &path_index, &num_optimal_paths);
  printf("[INFO] Number of optimal paths found: %d\n", num_optimal_paths);
  free(visited);
  free(path);
}

void print_all_paths_util(
  const char X[MAX_LENGTH+1],
  const char Y[MAX_LENGTH+1],
  const int n,
    const int trace[MAX_LENGTH+1][MAX_LENGTH+1],
    const match_t scores[MAX_LENGTH+1][MAX_LENGTH+1],
    const point_t u, const point_t d,
    bool* visited,
    point_t* path,
    int* path_index,
    int* num_optimal_paths) {
  if (u.x < 0 || u.y < 0) {
    return;
  }
  // Mark the current node and store it in path[]
  visited[get_cell_id(u, n)] = true;
  path[*path_index] = u;
  (*path_index)++;

  // If current cell is same as destination or STOP reached, then print path
  if ((u.x == d.x && u.y == d.y) || trace[u.x][u.y] == STOP) {
    char alignX[MAX_LENGTH*2]; // aligned X sequence
    char alignY[MAX_LENGTH*2]; // aligned Y sequence
    int alignment_length = 0;
    int x_idx = path[0].x;
    int y_idx = path[0].y;

    // TODO: The printing is NOT correct yet.
    for (int i = 0; i < (*path_index) - 1; ++i) {
      if (path[i].x-1 == path[i+1].x && path[i].y-1 == path[i+1].y) {
        // diag
        alignX[alignment_length] = X[path[i].x-1];
        alignY[alignment_length] = Y[path[i].y-1];
        --x_idx;
        --y_idx;
      } else if (path[i].x-1 == path[i+1].x && path[i].y == path[i+1].y) {
        // left
        alignX[alignment_length] = '-';
        alignY[alignment_length] = Y[path[i].y-1];
        --y_idx;
      } else if (path[i].x == path[i+1].x && path[i].y-1 == path[i+1].y) {
        // up
        alignX[alignment_length] = X[path[i].x-1];
        alignY[alignment_length] = '-';
        --x_idx;
      } 
      ++alignment_length;
    }
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

    (*num_optimal_paths)++;
    printf("[path n.%d]: ", *num_optimal_paths);
    for (int i = 0; i < *path_index; ++i) {
      // printf("(%s: %d, %d) ", 
      //   (trace[path[i].x][path[i].y] == DIAG) ? "DIAG" : ((trace[path[i].x][path[i].y] == LEFT) ? "LEFT" : "UP"),
      //   path[i].x, path[i].y);
      printf("(%d, %d) ", path[i].x, path[i].y);
    }
    printf("\n");
    print_alignment(alignment_length, alignX, alignY);
    printf("\n");
  } else {
    // If current cell is not the destination, call recursion on nodes with
    // same score directions and which haven't been visited yet.
    point_t up, diag, left;
    up.x = u.x - 1;
    up.y = u.y;
    diag.x = u.x - 1;
    diag.y = u.y - 1;
    left.x = u.x;
    left.y = u.y - 1;
    switch (trace[u.x][u.y]) {
      case UP:
        if (!visited[get_cell_id(up, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, up, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].up == scores[u.x][u.y].diag && !visited[get_cell_id(diag, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, diag, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].up == scores[u.x][u.y].left && !visited[get_cell_id(left, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, left, d, visited, path, path_index, num_optimal_paths);
        }
        break;
      case LEFT:
        if (!visited[get_cell_id(left, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, left, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].left == scores[u.x][u.y].diag && !visited[get_cell_id(diag, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, diag, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].left == scores[u.x][u.y].up && !visited[get_cell_id(up, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, up, d, visited, path, path_index, num_optimal_paths);
        }
        break;
      case DIAG:
        if (!visited[get_cell_id(diag, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, diag, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].diag == scores[u.x][u.y].up && !visited[get_cell_id(up, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, up, d, visited, path, path_index, num_optimal_paths);
        }
        if (scores[u.x][u.y].diag == scores[u.x][u.y].left && !visited[get_cell_id(left, n)]) {
          print_all_paths_util(X, Y, n, trace, scores, left, d, visited, path, path_index, num_optimal_paths);
        }
        break;
    }
  }
  // Remove current cell from path and set it as not visited
  (*path_index)--;
  visited[get_cell_id(u, n)] = false;
}

int main() {
  int i, j;
  int m, n;
  int alignment_length, score, tmp, up_score, left_score, diag_score;
  char X[MAX_LENGTH+1] = "ATTA"; // "ATCGAT";
  char Y[MAX_LENGTH+1] = "ATTTTA"; // "ATACGT";

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
  int max_score = -((1 << 30) - 1); // Fairly small number

  const match_t no_match = {0, 0, 0};
  match_t matches[MAX_LENGTH+1][MAX_LENGTH+1] = {no_match}; // branch matrix

  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; j++) {
      if (X[i-1] == Y[j-1]) {
        score = F[i-1][j-1] + MATCH_SCORE;
      } else {
        score = F[i-1][j-1] + MISMATCH_SCORE;
      }
      trace[i][j] = DIAG;
      tmp = F[i-1][j] - GAP_PENALTY;
      matches[i][j].diag = score;
      matches[i][j].up = tmp;
      if (tmp > score) {
        score = tmp;
        trace[i][j] = UP;
      }
      tmp = F[i][j-1] - GAP_PENALTY;
      matches[i][j].left = tmp;
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

  point_t start, end;
  start.x = m;
  start.y = n;
  end.x = 0;
  end.y = 0;
  print_all_paths(X, Y, m, n, trace, matches, start, end);
  return 0;
}
