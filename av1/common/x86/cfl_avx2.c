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

#include "av1/common/x86/cfl_simd.h"

/**
 * Adds 4 pixels (in a 2x2 grid) and multiplies them by 2. Resulting in a more
 * precise version of a box filter 4:2:0 pixel subsampling in Q3.
 *
 * The CfL prediction buffer is always of size CFL_BUF_SQUARE. However, the
 * active area is specified using width and height.
 *
 * Note: We don't need to worry about going over the active area, as long as we
 * stay inside the CfL prediction buffer.
 *
 * Note: For 4:2:0 luma subsampling, the width will never be greater than 16.
 */
static void cfl_luma_subsampling_420_lbd_avx2(const uint8_t *input,
                                              int input_stride,
                                              int16_t *pred_buf_q3, int width,
                                              int height) {
  (void)width;                               // Forever 32
  const __m256i twos = _mm256_set1_epi8(2);  // Thirty two twos
  const int luma_stride = input_stride << 1;
  __m256i *row = (__m256i *)pred_buf_q3;
  const __m256i *row_end = row + (height >> 1) * CFL_BUF_LINE_I256;
  do {
    __m256i top = _mm256_loadu_si256((__m256i *)input);
    __m256i bot = _mm256_loadu_si256((__m256i *)(input + input_stride));

    __m256i top_16x16 = _mm256_maddubs_epi16(top, twos);
    __m256i bot_16x16 = _mm256_maddubs_epi16(bot, twos);
    __m256i sum_16x16 = _mm256_add_epi16(top_16x16, bot_16x16);

    _mm256_storeu_si256(row, sum_16x16);

    input += luma_stride;
  } while ((row += CFL_BUF_LINE_I256) < row_end);
}

CFL_SUBSAMPLE(avx2, 420, lbd, 32, 32)
CFL_SUBSAMPLE(avx2, 420, lbd, 32, 16)
CFL_SUBSAMPLE(avx2, 420, lbd, 32, 8)

cfl_subsample_lbd_fn cfl_get_luma_subsampling_420_lbd_avx2(TX_SIZE tx_size) {
  static const cfl_subsample_lbd_fn subfn_420[TX_SIZES_ALL] = {
    subsample_lbd_420_4x4_ssse3,   /* 4x4 */
    subsample_lbd_420_8x8_ssse3,   /* 8x8 */
    subsample_lbd_420_16x16_ssse3, /* 16x16 */
    subsample_lbd_420_32x32_avx2,  /* 32x32 */
    cfl_subsample_lbd_null,        /* 64x64 (invalid CFL size) */
    subsample_lbd_420_4x8_ssse3,   /* 4x8 */
    subsample_lbd_420_8x4_ssse3,   /* 8x4 */
    subsample_lbd_420_8x16_ssse3,  /* 8x16 */
    subsample_lbd_420_16x8_ssse3,  /* 16x8 */
    subsample_lbd_420_16x32_ssse3, /* 16x32 */
    subsample_lbd_420_32x16_avx2,  /* 32x16 */
    cfl_subsample_lbd_null,        /* 32x64 (invalid CFL size) */
    cfl_subsample_lbd_null,        /* 64x32 (invalid CFL size) */
    subsample_lbd_420_4x16_ssse3,  /* 4x16  */
    subsample_lbd_420_16x4_ssse3,  /* 16x4  */
    subsample_lbd_420_8x32_ssse3,  /* 8x32  */
    subsample_lbd_420_32x8_avx2,   /* 32x8  */
    cfl_subsample_lbd_null,        /* 16x64 (invalid CFL size) */
    cfl_subsample_lbd_null,        /* 64x16 (invalid CFL size) */
  };
  return subfn_420[tx_size];
}

static INLINE __m256i predict_unclipped(const __m256i *input, __m256i alpha_q12,
                                        __m256i alpha_sign, __m256i dc_q0) {
  __m256i ac_q3 = _mm256_loadu_si256(input);
  __m256i ac_sign = _mm256_sign_epi16(alpha_sign, ac_q3);
  __m256i scaled_luma_q0 =
      _mm256_mulhrs_epi16(_mm256_abs_epi16(ac_q3), alpha_q12);
  scaled_luma_q0 = _mm256_sign_epi16(scaled_luma_q0, ac_sign);
  return _mm256_add_epi16(scaled_luma_q0, dc_q0);
}

static INLINE void cfl_predict_lbd_avx2(const int16_t *pred_buf_q3,
                                        uint8_t *dst, int dst_stride,
                                        int alpha_q3, int width, int height) {
  (void)width;
  const __m256i alpha_sign = _mm256_set1_epi16(alpha_q3);
  const __m256i alpha_q12 = _mm256_slli_epi16(_mm256_abs_epi16(alpha_sign), 9);
  const __m256i dc_q0 = _mm256_set1_epi16(*dst);
  __m256i *row = (__m256i *)pred_buf_q3;
  const __m256i *row_end = row + height * CFL_BUF_LINE_I256;

  do {
    __m256i res = predict_unclipped(row, alpha_q12, alpha_sign, dc_q0);
    __m256i next = predict_unclipped(row + 1, alpha_q12, alpha_sign, dc_q0);
    res = _mm256_packus_epi16(res, next);
    res = _mm256_permute4x64_epi64(res, _MM_SHUFFLE(3, 1, 2, 0));
    _mm256_storeu_si256((__m256i *)dst, res);
    dst += dst_stride;
  } while ((row += CFL_BUF_LINE_I256) < row_end);
}

CFL_PREDICT_X(avx2, 32, 8, lbd);
CFL_PREDICT_X(avx2, 32, 16, lbd);
CFL_PREDICT_X(avx2, 32, 32, lbd);

cfl_predict_lbd_fn get_predict_lbd_fn_avx2(TX_SIZE tx_size) {
  static const cfl_predict_lbd_fn pred[TX_SIZES_ALL] = {
    predict_lbd_4x4_ssse3,   /* 4x4 */
    predict_lbd_8x8_ssse3,   /* 8x8 */
    predict_lbd_16x16_ssse3, /* 16x16 */
    predict_lbd_32x32_avx2,  /* 32x32 */
    cfl_predict_lbd_null,    /* 64x64 (invalid CFL size) */
    predict_lbd_4x8_ssse3,   /* 4x8 */
    predict_lbd_8x4_ssse3,   /* 8x4 */
    predict_lbd_8x16_ssse3,  /* 8x16 */
    predict_lbd_16x8_ssse3,  /* 16x8 */
    predict_lbd_16x32_ssse3, /* 16x32 */
    predict_lbd_32x16_avx2,  /* 32x16 */
    cfl_predict_lbd_null,    /* 32x64 (invalid CFL size) */
    cfl_predict_lbd_null,    /* 64x32 (invalid CFL size) */
    predict_lbd_4x16_ssse3,  /* 4x16  */
    predict_lbd_16x4_ssse3,  /* 16x4  */
    predict_lbd_8x32_ssse3,  /* 8x32  */
    predict_lbd_32x8_avx2,   /* 32x8  */
    cfl_predict_lbd_null,    /* 16x64 (invalid CFL size) */
    cfl_predict_lbd_null,    /* 64x16 (invalid CFL size) */
  };
  /* Modulo TX_SIZES_ALL to ensure that an attacker won't be able to
              */ /* index the function pointer array out of bounds. */
  return pred[tx_size % TX_SIZES_ALL];
}

static __m256i highbd_max_epi16(int bd) {
  const __m256i neg_one = _mm256_set1_epi16(-1);
  // (1 << bd) - 1 => -(-1 << bd) -1 => -1 - (-1 << bd) => -1 ^ (-1 << bd)
  return _mm256_xor_si256(_mm256_slli_epi16(neg_one, bd), neg_one);
}

static __m256i highbd_clamp_epi16(__m256i u, __m256i zero, __m256i max) {
  return _mm256_max_epi16(_mm256_min_epi16(u, max), zero);
}

static INLINE void cfl_predict_hbd_avx2(const int16_t *pred_buf_q3,
                                        uint16_t *dst, int dst_stride,
                                        int alpha_q3, int bd, int width,
                                        int height) {
  // Use SSSE3 version for smaller widths
  assert(width == 16 || width == 32);
  const __m256i alpha_sign = _mm256_set1_epi16(alpha_q3);
  const __m256i alpha_q12 = _mm256_slli_epi16(_mm256_abs_epi16(alpha_sign), 9);
  const __m256i dc_q0 = _mm256_loadu_si256((__m256i *)dst);
  const __m256i max = highbd_max_epi16(bd);

  __m256i *row = (__m256i *)pred_buf_q3;
  const __m256i *row_end = row + height * CFL_BUF_LINE_I256;
  do {
    __m256i res = predict_unclipped(row, alpha_q12, alpha_sign, dc_q0);
    _mm256_storeu_si256((__m256i *)dst,
                        highbd_clamp_epi16(res, _mm256_setzero_si256(), max));
    if (width == 32) {
      res = predict_unclipped(row + 1, alpha_q12, alpha_sign, dc_q0);
      _mm256_storeu_si256((__m256i *)(dst + 16),
                          highbd_clamp_epi16(res, _mm256_setzero_si256(), max));
    }
    dst += dst_stride;
  } while ((row += CFL_BUF_LINE_I256) < row_end);
}

CFL_PREDICT_X(avx2, 16, 4, hbd)
CFL_PREDICT_X(avx2, 16, 8, hbd)
CFL_PREDICT_X(avx2, 16, 16, hbd)
CFL_PREDICT_X(avx2, 16, 32, hbd)
CFL_PREDICT_X(avx2, 32, 8, hbd)
CFL_PREDICT_X(avx2, 32, 16, hbd)
CFL_PREDICT_X(avx2, 32, 32, hbd)

cfl_predict_hbd_fn get_predict_hbd_fn_avx2(TX_SIZE tx_size) {
  static const cfl_predict_hbd_fn pred[TX_SIZES_ALL] = {
    predict_hbd_4x4_ssse3,  /* 4x4 */
    predict_hbd_8x8_ssse3,  /* 8x8 */
    predict_hbd_16x16_avx2, /* 16x16 */
    predict_hbd_32x32_avx2, /* 32x32 */
    cfl_predict_hbd_null,   /* 64x64 (invalid CFL size) */
    predict_hbd_4x8_ssse3,  /* 4x8 */
    predict_hbd_8x4_ssse3,  /* 8x4 */
    predict_hbd_8x16_ssse3, /* 8x16 */
    predict_hbd_16x8_avx2,  /* 16x8 */
    predict_hbd_16x32_avx2, /* 16x32 */
    predict_hbd_32x16_avx2, /* 32x16 */
    cfl_predict_hbd_null,   /* 32x64 (invalid CFL size) */
    cfl_predict_hbd_null,   /* 64x32 (invalid CFL size) */
    predict_hbd_4x16_ssse3, /* 4x16  */
    predict_hbd_16x4_avx2,  /* 16x4  */
    predict_hbd_8x32_ssse3, /* 8x32  */
    predict_hbd_32x8_avx2,  /* 32x8  */
    cfl_predict_hbd_null,   /* 16x64 (invalid CFL size) */
    cfl_predict_hbd_null,   /* 64x16 (invalid CFL size) */
  };
  /* Modulo TX_SIZES_ALL to ensure that an attacker won't be able to
    */ /* index the function pointer array out of bounds. */
  return pred[tx_size % TX_SIZES_ALL];
}

// Returns a vector where all the (32-bits) elements are the sum of all the
// lanes in a.
static INLINE __m256i fill_sum_epi32(__m256i a) {
  // Given that a == [A, B, C, D, E, F, G, H]
  a = _mm256_hadd_epi32(a, a);
  // Given that A' == A + B, C' == C + D, E' == E + F, G' == G + H
  // a == [A', C', A', C', E', G', E', G']
  a = _mm256_permute4x64_epi64(a, _MM_SHUFFLE(3, 1, 2, 0));
  // a == [A', C', E', G', A', C', E', G']
  a = _mm256_hadd_epi32(a, a);
  // Given that A'' == A' + C' and E'' == E' + G'
  // a == [A'', E'', A'', E'', A'', E'', A'', E'']
  return _mm256_hadd_epi32(a, a);
  // Given that A''' == A'' + E''
  // a == [A''', A''', A''', A''', A''', A''', A''', A''']
}

static INLINE __m256i _mm256_addl_epi16(__m256i a) {
  return _mm256_add_epi32(_mm256_unpacklo_epi16(a, _mm256_setzero_si256()),
                          _mm256_unpackhi_epi16(a, _mm256_setzero_si256()));
}

static INLINE void subtract_average_avx2(int16_t *pred_buf, int width,
                                         int height, int round_offset,
                                         int num_pel_log2) {
  // Use SSE2 version for smaller widths
  assert(width == 16 || width == 32);
  __m256i *row = (__m256i *)pred_buf;
  const __m256i *const end = row + height * CFL_BUF_LINE_I256;
  // To maximize usage of the AVX2 registers, we sum two rows per loop
  // iteration
  const int step = 2 * CFL_BUF_LINE_I256;
  __m256i sum = _mm256_setzero_si256();

  // For width 32, we use a second sum accumulator to reduce accumulator
  // dependencies in the loop.
  __m256i sum2;
  if (width == 32) sum2 = _mm256_setzero_si256();
  do {
    // Add top row to the bottom row
    __m256i l0 = _mm256_add_epi16(_mm256_loadu_si256(row),
                                  _mm256_loadu_si256(row + CFL_BUF_LINE_I256));
    sum = _mm256_add_epi32(sum, _mm256_addl_epi16(l0));
    if (width == 32) { /* Don't worry, this if it gets optimized out. */
      // Add the second part of the top row to the second part of the bottom row
      __m256i l1 =
          _mm256_add_epi16(_mm256_loadu_si256(row + 1),
                           _mm256_loadu_si256(row + 1 + CFL_BUF_LINE_I256));
      // Store the sum of the second part in the same accumulator as the first
      // part
      sum2 = _mm256_add_epi32(sum2, _mm256_addl_epi16(l1));
    }
  } while ((row += step) < end);
  // Combine both sum accumulator
  if (width == 32) sum = _mm256_add_epi32(sum, sum2);

  // The sum accumulator now contains the 8 lanes
  __m256i fill = fill_sum_epi32(sum);

  __m256i avg_epi16 = _mm256_srli_epi32(
      _mm256_add_epi32(fill, _mm256_set1_epi32(round_offset)), num_pel_log2);
  avg_epi16 = _mm256_packs_epi32(avg_epi16, avg_epi16);

  // Store and subtract loop
  row = (__m256i *)pred_buf;
  do {
    _mm256_storeu_si256(row,
                        _mm256_sub_epi16(_mm256_loadu_si256(row), avg_epi16));
    if (width == 32) {
      _mm256_storeu_si256(
          row + 1, _mm256_sub_epi16(_mm256_loadu_si256(row + 1), avg_epi16));
    }
  } while ((row += CFL_BUF_LINE_I256) < end);
}

// Declare wrappers for AVX2 sizes
CFL_SUB_AVG_X(avx2, 16, 4, 32, 6)
CFL_SUB_AVG_X(avx2, 16, 8, 64, 7)
CFL_SUB_AVG_X(avx2, 16, 16, 128, 8)
CFL_SUB_AVG_X(avx2, 16, 32, 256, 9)
CFL_SUB_AVG_X(avx2, 32, 8, 128, 8)
CFL_SUB_AVG_X(avx2, 32, 16, 256, 9)
CFL_SUB_AVG_X(avx2, 32, 32, 512, 10)

// Based on the observation that for small blocks AVX2 does not outperform
// SSE2, we call the SSE2 code for block widths 4 and 8.
cfl_subtract_average_fn get_subtract_average_fn_avx2(TX_SIZE tx_size) {
  static const cfl_subtract_average_fn sub_avg[TX_SIZES_ALL] = {
    subtract_average_4x4_sse2,   /* 4x4 */
    subtract_average_8x8_sse2,   /* 8x8 */
    subtract_average_16x16_avx2, /* 16x16 */
    subtract_average_32x32_avx2, /* 32x32 */
    cfl_subtract_average_null,   /* 64x64 (invalid CFL size) */
    subtract_average_4x8_sse2,   /* 4x8 */
    subtract_average_8x4_sse2,   /* 8x4 */
    subtract_average_8x16_sse2,  /* 8x16 */
    subtract_average_16x8_avx2,  /* 16x8 */
    subtract_average_16x32_avx2, /* 16x32 */
    subtract_average_32x16_avx2, /* 32x16 */
    cfl_subtract_average_null,   /* 32x64 (invalid CFL size) */
    cfl_subtract_average_null,   /* 64x32 (invalid CFL size) */
    subtract_average_4x16_sse2,  /* 4x16 */
    subtract_average_16x4_avx2,  /* 16x4 */
    subtract_average_8x32_sse2,  /* 8x32 */
    subtract_average_32x8_avx2,  /* 32x8 */
    cfl_subtract_average_null,   /* 16x64 (invalid CFL size) */
    cfl_subtract_average_null,   /* 64x16 (invalid CFL size) */
  };
  // Modulo TX_SIZES_ALL to ensure that an attacker won't be able to
  // index the function pointer array out of bounds.
  return sub_avg[tx_size % TX_SIZES_ALL];
}
