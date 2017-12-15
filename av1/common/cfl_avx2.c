/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <immintrin.h>

#include "./av1_rtcd.h"

#include "av1/common/cfl.h"

void av1_subtract_average_avx2(int16_t *pred_buf_q3, int width, int height,
                               int16_t avg_q3) {
  const __m256i avg = _mm256_set1_epi16(avg_q3);
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i += 16) {
      // The CfL prediction buffer is 32x32, and CfL is capped at 32x32.
      // Writing outside of the boundary of a smaller block is harmless.
      // We don't need to manage borders
      __m256i vals = _mm256_loadu_si256((__m256i *)(pred_buf_q3 + i));
      _mm256_storeu_si256((__m256i *)(pred_buf_q3 + i),
                          _mm256_sub_epi16(vals, avg));
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}
