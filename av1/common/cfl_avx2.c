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

void av1_cfl_subtract_avx2(int16_t *pred_buf_q3, int width, int height,
                           int16_t avg_q3) {
  const __m256i avg = _mm256_set1_epi16(avg_q3);
  // The CfL prediction buffer is 32x32, and CfL is capped at 32x32.
  // Writing outside of the boundary of a smaller block is harmless.
  // We don't need to manage borders
  //
  // If the width of the block is less than or equal to 16, the whole row will
  // fit in 1 register. If the width of the block is 32, the row will fit
  // exactly in two registers.
  //
  // If the whole row fits in 1 register, we jump to the next row (stride = 32),
  // else we jump to the end of the register (stride == 16).
  const int stride = CFL_BUF_LINE - ((width > 16) << 4);
  // If the whole row does not fit in 1 register, we will have twice as many
  // iterations to perform.
  const int16_t *end = pred_buf_q3 + (height << (width > 16)) * stride;
  do {
    __m256i vals0 = _mm256_loadu_si256((__m256i *)pred_buf_q3);
    _mm256_storeu_si256((__m256i *)pred_buf_q3, _mm256_sub_epi16(vals0, avg));
  } while ((pred_buf_q3 += stride) < end);
}
