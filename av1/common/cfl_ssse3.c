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

#include <tmmintrin.h>

#include "./av1_rtcd.h"

#include "av1/common/cfl.h"

void av1_cfl_build_prediction_lbd_ssse3(const int16_t *pred_buf_q3,
                                        uint8_t *dst, int dst_stride, int width,
                                        int height, int alpha_q3) {
  const __m128i zeros = _mm_setzero_si128();
  const __m128i alpha_q12 = _mm_set1_epi16(abs(alpha_q3) * (1 << 9));
  const __m128i alpha_sign = alpha_q3 < 0 ? _mm_set1_epi16(-1) : zeros;
  const __m128i dc_packed = _mm_loadu_si128((__m128i *)(dst));
  const __m128i dc_q0 = _mm_unpacklo_epi8(dc_packed, zeros);

  uint8_t *row_end = dst + height * dst_stride;
  do {
    for (int m = 0; m < width; m += 8) {
      __m128i ac_q3 = _mm_loadu_si128((__m128i *)(pred_buf_q3 + m));
      __m128i ac_sign = _mm_srai_epi16(ac_q3, 15);
      ac_q3 = _mm_xor_si128(ac_q3, ac_sign);
      ac_q3 = _mm_sub_epi16(ac_q3, ac_sign);
      ac_sign = _mm_xor_si128(ac_sign, alpha_sign);
      __m128i scaled_luma_q0 = _mm_mulhrs_epi16(ac_q3, alpha_q12);
      scaled_luma_q0 = _mm_xor_si128(scaled_luma_q0, ac_sign);
      scaled_luma_q0 = _mm_sub_epi16(scaled_luma_q0, ac_sign);
      __m128i tmp = _mm_add_epi16(scaled_luma_q0, dc_q0);
      __m128i res = _mm_packus_epi16(tmp, tmp);
      _mm_storel_epi64((__m128i *)(dst + m), res);
    }
    dst += dst_stride;
    pred_buf_q3 += CFL_BUF_LINE;
  } while (dst < row_end);
}
