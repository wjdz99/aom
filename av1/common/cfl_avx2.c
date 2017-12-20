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

static void cfl_luma_subsampling_420_lbd_avx2(const uint8_t *input,
                                              int input_stride,
                                              int16_t *output_q3, int width,
                                              int height) {
  const __m256i twos = _mm256_set1_epi8(2);
  const int16_t *row_end = output_q3 + height * CFL_BUF_LINE;
  const int luma_stride = input_stride << 1;
  (void)width;
  do {
    __m256i top = _mm256_loadu_si256((__m256i *)(input));
    __m256i bot = _mm256_loadu_si256((__m256i *)(input + input_stride));

    top = _mm256_maddubs_epi16(top, twos);
    bot = _mm256_maddubs_epi16(bot, twos);

    _mm256_storeu_si256((__m256i *)output_q3, _mm256_add_epi16(top, bot));

    input += luma_stride;
    output_q3 += CFL_BUF_LINE;
  } while (output_q3 < row_end);
}

cfl_subsample_lbd_fn get_subsample_lbd_fn_avx2(int sub_x, int sub_y) {
  static const cfl_subsample_lbd_fn subsample_lbd[2][2] = {
    //  (sub_y == 0, sub_x == 0)       (sub_y == 0, sub_x == 1)
    //  (sub_y == 1, sub_x == 0)       (sub_y == 1, sub_x == 1)
    { cfl_luma_subsampling_444_lbd, cfl_luma_subsampling_422_lbd },
    { cfl_luma_subsampling_440_lbd, cfl_luma_subsampling_420_lbd_avx2 },
  };
  return subsample_lbd[sub_y & 1][sub_x & 1];
}

static INLINE __m256i get_scaled_luma_q0_avx2(const __m256i *input,
                                              __m256i alpha_q12,
                                              __m256i alpha_sign) {
  __m256i ac_q3 = _mm256_loadu_si256(input);
  __m256i ac_sign = _mm256_srai_epi16(ac_q3, 15);
  ac_q3 = _mm256_xor_si256(ac_q3, ac_sign);
  ac_q3 = _mm256_sub_epi16(ac_q3, ac_sign);
  ac_sign = _mm256_xor_si256(ac_sign, alpha_sign);
  __m256i scaled_luma_q0 = _mm256_mulhrs_epi16(ac_q3, alpha_q12);
  scaled_luma_q0 = _mm256_xor_si256(scaled_luma_q0, ac_sign);
  scaled_luma_q0 = _mm256_sub_epi16(scaled_luma_q0, ac_sign);
  return scaled_luma_q0;
}

void av1_cfl_build_prediction_lbd_avx2(const int16_t *pred_buf_q3, uint8_t *dst,
                                       int dst_stride, int width, int height,
                                       int alpha_q3) {
  const __m256i zeros = _mm256_setzero_si256();
  const __m256i alpha_q12 = _mm256_set1_epi16(abs(alpha_q3) * (1 << 9));
  const __m256i alpha_sign = alpha_q3 < 0 ? _mm256_set1_epi16(-1) : zeros;
  const __m256i dc_packed = _mm256_loadu_si256((__m256i *)dst);
  const __m256i dc_q0 = _mm256_unpacklo_epi8(dc_packed, zeros);

  const int16_t *row_end = pred_buf_q3 + height * CFL_BUF_LINE;
  (void)width;
  do {
    __m256i scaled_luma_q0 =
        get_scaled_luma_q0_avx2((__m256i *)pred_buf_q3, alpha_q12, alpha_sign);
    __m256i tmp0 = _mm256_add_epi16(scaled_luma_q0, dc_q0);
    if (width >= 16)
      scaled_luma_q0 = get_scaled_luma_q0_avx2((__m256i *)(pred_buf_q3 + 8),
                                               alpha_q12, alpha_sign);
    __m256i tmp1 = _mm256_add_epi16(scaled_luma_q0, dc_q0);
    __m256i res = _mm256_packus_epi16(tmp0, tmp1);
    _mm256_storeu_si256((__m256i *)dst, res);
    if (width == 32) {
      scaled_luma_q0 = get_scaled_luma_q0_avx2((__m256i *)(pred_buf_q3 + 16),
                                               alpha_q12, alpha_sign);
      tmp0 = _mm256_add_epi16(scaled_luma_q0, dc_q0);
      scaled_luma_q0 = get_scaled_luma_q0_avx2((__m256i *)(pred_buf_q3 + 24),
                                               alpha_q12, alpha_sign);
      tmp1 = _mm256_add_epi16(scaled_luma_q0, dc_q0);
      res = _mm256_packus_epi16(tmp0, tmp1);
      _mm256_storeu_si256((__m256i *)(dst + 16), res);
    }
    dst += dst_stride;
    pred_buf_q3 += CFL_BUF_LINE;
  } while (pred_buf_q3 < row_end);
}
