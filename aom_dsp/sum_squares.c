/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>

#include "config/aom_dsp_rtcd.h"

uint64_t aom_sum_squares_2d_i16_c(const int16_t *src, int src_stride, int width,
                                  int height, double *log_var) {
  int r, c;
  uint64_t ss = 0;
  int64_t sum = 0;
  double temp_var;

  for (r = 0; r < height; r++) {
    for (c = 0; c < width; c++) {
      const int16_t v = src[c];
      ss += v * v;
      sum += v;
    }
    src += src_stride;
  }
  if (log_var != NULL) {
    temp_var = (double)(ss - ((double)(sum * sum) / (double)(width * height)));
    *log_var = (((double)256 * temp_var) / (double)(width * height));
  }
  return ss;
}

uint64_t aom_sum_squares_i16_c(const int16_t *src, uint32_t n) {
  uint64_t ss = 0;
  do {
    const int16_t v = *src++;
    ss += v * v;
  } while (--n);

  return ss;
}
