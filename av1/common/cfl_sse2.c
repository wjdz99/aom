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
  // The CfL prediction buffer is 32x32, and CfL is capped at 32x32.
  // Writing outside of the boundary of a smaller block is harmless.
  // We don't need to manage borders
  //
  // If the width of the block is less than or equal to 8, the whole row will
  // fit in 1 register. If the width of the block is greater than 8, the row
  // will fit in four registers.
  //
  // If the whole row fits in 1 register, we jump to the next row (stride = 32),
  // else we jump to the end of the register (stride == 8).
  const int stride = CFL_BUF_LINE >> (2 * (width > 8));
  // If the whole row does not fit in 1 register, we will have four times as
  // many iterations to perform.
  const int16_t *end = pred_buf_q3 + (height << (2 * (width > 8))) * stride;
  do {
    __m128i vals0 = _mm_load_si128((__m128i *)pred_buf_q3);
    __m128i vals1 = _mm_load_si128((__m128i *)(pred_buf_q3 + stride));
    __m128i vals2 = _mm_load_si128((__m128i *)(pred_buf_q3 + stride * 2));
    __m128i vals3 = _mm_load_si128((__m128i *)(pred_buf_q3 + stride * 3));

    _mm_store_si128((__m128i *)pred_buf_q3, _mm_sub_epi16(vals0, avg));
    _mm_store_si128((__m128i *)(pred_buf_q3 + stride),
                    _mm_sub_epi16(vals1, avg));
    _mm_store_si128((__m128i *)(pred_buf_q3 + stride * 2),
                    _mm_sub_epi16(vals2, avg));
    _mm_store_si128((__m128i *)(pred_buf_q3 + stride * 3),
                    _mm_sub_epi16(vals3, avg));

    pred_buf_q3 += stride << 2;
  } while (pred_buf_q3 < end);
}
