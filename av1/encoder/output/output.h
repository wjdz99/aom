#ifndef AOM_AV1_ENCODER_OUTPUT_H_
#define AOM_AV1_ENCODER_OUTPUT_H_

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CHECK(b)                                       \
  do {                                                 \
    if (!(b)) {                                        \
      printf("ERROR at %s, %d\n", __FILE__, __LINE__); \
      exit(1);                                         \
    }                                                  \
  } while (0);

typedef struct Elliottk {
  int bsize;
  int source_variance;
  int mi_row;
  int mi_col;
  int64_t rd_cost;
  int exit_point;
  int best_mode_index;
} Elliottk;

void srecordio_output(const Elliottk *elliottk);
void srecordio_path(const char *path);
char *get_path();
void srecordio_output_close();
int get_iter();

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_OUTPUT_H_
