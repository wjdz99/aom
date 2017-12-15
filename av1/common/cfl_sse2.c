
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

#include <emmintrin.h>

#include "./av1_rtcd.h"

#include "av1/common/cfl.h"

void av1_subtract_average_sse2(int16_t *pred_buf_q3, int width, int height,
                               int16_t avg_q3) {
  const __m128i avg = _mm_set1_epi16(avg_q3);
  __m128i *pred_buf;
  for (int j = 0; j < height; j++) {
    pred_buf = (__m128i *)pred_buf_q3;
    for (int i = 0; i < width; i += 8) {
      __m128i vals = _mm_load_si128(pred_buf);
      _mm_store_si128(pred_buf, _mm_sub_epi16(vals, avg));
      pred_buf++;
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}
