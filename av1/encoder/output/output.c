#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#include "av1/encoder/output/output.h"

#define BUF_SIZE 500
static char PATH[BUF_SIZE] = "";
static FILE *writer = NULL;

int get_iter() {
  const char *s = getenv("ITER_NUM");
  if (s == NULL || strlen(s) == 0) {
    return -1;
  }
  return atoi(s);
}

#define CHECK_OK(s)                                                \
  do {                                                             \
    srecordio_status stsjse = (s);                                 \
    if (stsjse != SRECORDIO_OK) {                                  \
      printf(                                                      \
          "Check failed in function %s, line number %d, file %s\n" \
          "Error: %d\n",                                           \
          __func__, __LINE__, __FILE__, stsjse);                   \
      exit(2);                                                     \
    }                                                              \
  } while (0);

static const char EXT[] = "srio";

static int find_last(const char *s, char c) {
  size_t len = strlen(s);
  for (int i = len - 1; i >= 0; --i) {
    if (s[i] == c) {
      return i;
    }
  }
  return -1;
}

void srecordio_path(const char *path) {
  // Find the last slash and use anything after it.
  int last_slash = find_last(path, '/');
  int last_period = find_last(path, '.');
  if (last_period == -1 || last_period > BUF_SIZE - (int)strlen(EXT) - 2 ||
      last_period < last_slash) {
    printf("Error: unable to find extension\n");
    exit(23);
  }
  strncpy(PATH, path, last_period + 1);
  strncpy(PATH + last_period + 1, EXT, strlen(EXT) + 1);
  printf("Path is now: %s\n", PATH);
}

void srecordio_output(const Elliottk *elliottk) {
  if (strlen(PATH) == 0) {
    printf("PATH NOT SET, ERROR\n");
    exit(3);
  }
  if (writer == NULL) {
    writer = fopen(PATH, "wb");
    int r = fprintf(writer,
                    "bsize,source_variance,mi_row,mi_col,rd_cost,exit_point,"
                    "best_mode_index\n");
    CHECK(r > 0);
  }
  int r = fprintf(writer, "%d,%d,%d,%d,%d,%d,%d\n", elliottk->bsize,
                  elliottk->source_variance, elliottk->mi_row, elliottk->mi_col,
                  elliottk->rd_cost, elliottk->exit_point,
                  elliottk->best_mode_index);
  CHECK(r > 0);
}

char *get_basename() {
  size_t l = strlen(PATH);
  if (l == 0) {
    printf("PATH NOT SET, ERROR\n");
    exit(3);
  }
  int slash_loc;
  for (slash_loc = l - 1; slash_loc >= 0; --slash_loc) {
    if (PATH[slash_loc] == '/') {
      break;
    }
  }
  return PATH + slash_loc + 1;
}

void srecordio_output_close() {
  if (writer != NULL) {
    fclose(writer);
    writer = NULL;
  }
}
