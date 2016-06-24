#include <assert.h>
#include <stdio.h>
#include "aom_util/debug_util.h"
#if CONFIG_BITSTREAM_DEBUG
int result_queue[1000000];
int prob_queue[1000000];
int queue_max_size = 1000000;
int queue_r = 0;
int queue_w = 0;
int queue_prev_w = -1;
int pop_cnt = 0;
int push_cnt = 0;

void bitstream_queue_record() {
  queue_prev_w = queue_w;
  push_cnt = 0;
}

void bitstream_queue_reset() {
  queue_w = queue_prev_w;
  push_cnt = 0;
}

void bitstream_queue_pop(int* result, int* prob) {
  if (queue_w == queue_r) {
    printf("buffer underflow queue_w %d queue_r %d\n", queue_w, queue_r);
    assert(0);
  }
  *result = result_queue[queue_r];
  *prob = prob_queue[queue_r];
  queue_r = (queue_r + 1) % queue_max_size;
}

void bitstream_queue_push(int result, int prob) {
  result_queue[queue_w] = result;
  prob_queue[queue_w] = prob;
  queue_w = (queue_w + 1) % queue_max_size;
  if (queue_w == queue_r) {
    printf("buffer overflow queue_w %d queue_r %d\n", queue_w, queue_r);
    assert(0);
  }
}
#endif
