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
static void cfl_luma_subsampling_420_lbd_ssse3(const uint8_t *input,
                                               int input_stride,
                                               int16_t *pred_buf_q3, int width,
                                               int height) {
  const __m128i twos = _mm_set1_epi8(2);  // Sixteen twos

  // Sixteen int8 values fit in one __m128i register. If this is enough to do
  // the entire row, the next value is two rows down, otherwise we move to the
  // next sixteen values.
  const int next = (width == 16) ? 16 : input_stride << 1;

  // Values in the prediction buffer are subsampled, so we only need to move
  // down one row or forward by eight values.
  const int next_chroma = (width == 16) ? 8 : CFL_BUF_LINE;

  // When the width is less than 16, we double the stride, because we process
  // four lines by iteration (instead of two).
  const int luma_stride = input_stride << (1 + (width < 16));
  const int chroma_stride = CFL_BUF_LINE << (width < 16);

  const int16_t *end = pred_buf_q3 + height * CFL_BUF_LINE;
  do {
    // Load 16 values for the top and bottom rows.
    // t_0, t_1, ... t_15
    __m128i top = _mm_loadu_si128((__m128i *)(input));
    // b_0, b_1, ... b_15
    __m128i bot = _mm_loadu_si128((__m128i *)(input + input_stride));

    // Load either the next line or the next 16 values
    __m128i next_top = _mm_loadu_si128((__m128i *)(input + next));
    __m128i next_bot =
        _mm_loadu_si128((__m128i *)(input + next + input_stride));

    // Horizontal add of the 16 values into 8 values that are multiplied by 2
    // (t_0 + t_1) * 2, (t_2 + t_3) * 2, ... (t_14 + t_15) *2
    top = _mm_maddubs_epi16(top, twos);
    next_top = _mm_maddubs_epi16(next_top, twos);
    // (b_0 + b_1) * 2, (b_2 + b_3) * 2, ... (b_14 + b_15) *2
    bot = _mm_maddubs_epi16(bot, twos);
    next_bot = _mm_maddubs_epi16(next_bot, twos);

    // Add the 8 values in top with the 8 values in bottom
    _mm_storeu_si128((__m128i *)pred_buf_q3, _mm_add_epi16(top, bot));
    _mm_storeu_si128((__m128i *)(pred_buf_q3 + next_chroma),
                     _mm_add_epi16(next_top, next_bot));

    input += luma_stride;
    pred_buf_q3 += chroma_stride;
  } while (pred_buf_q3 < end);
}

cfl_subsample_lbd_fn get_subsample_lbd_fn_ssse3(int sub_x, int sub_y) {
  static const cfl_subsample_lbd_fn subsample_lbd[2][2] = {
    //  (sub_y == 0, sub_x == 0)       (sub_y == 0, sub_x == 1)
    //  (sub_y == 1, sub_x == 0)       (sub_y == 1, sub_x == 1)
    { cfl_luma_subsampling_444_lbd, cfl_luma_subsampling_422_lbd },
    { cfl_luma_subsampling_440_lbd, cfl_luma_subsampling_420_lbd_ssse3 },
  };
  // AND sub_x and sub_y with 1 to ensures that an attacker won't be able to
  // index the function pointer array out of bounds.
  return subsample_lbd[sub_y & 1][sub_x & 1];
}

static INLINE __m128i predict_unclipped(const __m128i *input, __m128i alpha_q12,
                                        __m128i alpha_sign, __m128i dc_q0) {
  __m128i ac_q3 = _mm_loadu_si128(input);
  __m128i ac_sign = _mm_sign_epi16(alpha_sign, ac_q3);
  __m128i scaled_luma_q0 = _mm_mulhrs_epi16(_mm_abs_epi16(ac_q3), alpha_q12);
  scaled_luma_q0 = _mm_sign_epi16(scaled_luma_q0, ac_sign);
  return _mm_add_epi16(scaled_luma_q0, dc_q0);
}

static INLINE void cfl_predict_lbd_x(const int16_t *pred_buf_q3, uint8_t *dst,
                                     int dst_stride, TX_SIZE tx_size,
                                     int alpha_q3, int width) {
  uint8_t *row_end = dst + tx_size_high[tx_size] * dst_stride;
  const __m128i alpha_sign = _mm_set1_epi16(alpha_q3);
  const __m128i alpha_q12 = _mm_slli_epi16(_mm_abs_epi16(alpha_sign), 9);
  const __m128i dc_q0 = _mm_set1_epi16(*dst);
  do {
    __m128i res = predict_unclipped((__m128i *)(pred_buf_q3), alpha_q12,
                                    alpha_sign, dc_q0);
    if (width < 16) {
      res = _mm_packus_epi16(res, res);
      if (width == 4)
        *(uint32_t *)dst = _mm_cvtsi128_si32(res);
      else
        _mm_storel_epi64((__m128i *)dst, res);
    } else {
      __m128i next = predict_unclipped((__m128i *)(pred_buf_q3 + 8), alpha_q12,
                                       alpha_sign, dc_q0);
      res = _mm_packus_epi16(res, next);
      _mm_storeu_si128((__m128i *)dst, res);
      if (width == 32) {
        res = predict_unclipped((__m128i *)(pred_buf_q3 + 16), alpha_q12,
                                alpha_sign, dc_q0);
        next = predict_unclipped((__m128i *)(pred_buf_q3 + 24), alpha_q12,
                                 alpha_sign, dc_q0);
        res = _mm_packus_epi16(res, next);
        _mm_storeu_si128((__m128i *)(dst + 16), res);
      }
    }
    dst += dst_stride;
    pred_buf_q3 += CFL_BUF_LINE;
  } while (dst < row_end);
}

static INLINE __m128i highbd_max_epi16(int bd) {
  const __m128i neg_one = _mm_set1_epi16(-1);
  // (1 << bd) - 1 => -(-1 << bd) -1 => -1 - (-1 << bd) => -1 ^ (-1 << bd)
  return _mm_xor_si128(_mm_slli_epi16(neg_one, bd), neg_one);
}

static INLINE __m128i highbd_clamp_epi16(__m128i u, __m128i zero, __m128i max) {
  return _mm_max_epi16(_mm_min_epi16(u, max), zero);
}

static INLINE void cfl_predict_hbd(__m128i *dst, __m128i *src,
                                   __m128i alpha_q12, __m128i alpha_sign,
                                   __m128i dc_q0, __m128i zero, __m128i max) {
  __m128i res = predict_unclipped(src, alpha_q12, alpha_sign, dc_q0);
  _mm_storeu_si128(dst, highbd_clamp_epi16(res, zero, max));
}

static INLINE void cfl_predict_hbd_x(const int16_t *pred_buf_q3, uint16_t *dst,
                                     int dst_stride, TX_SIZE tx_size,
                                     int alpha_q3, int bd, int width) {
  uint16_t *row_end = dst + tx_size_high[tx_size] * dst_stride;
  const __m128i alpha_sign = _mm_set1_epi16(alpha_q3);
  const __m128i alpha_q12 = _mm_slli_epi16(_mm_abs_epi16(alpha_sign), 9);
  const __m128i dc_q0 = width == 4 ? _mm_loadl_epi64((__m128i *)dst)
                                   : _mm_load_si128((__m128i *)dst);
  const __m128i max = highbd_max_epi16(bd);
  const __m128i zero = _mm_setzero_si128();
  do {
    if (width == 4) {
      __m128i res = predict_unclipped((__m128i *)(pred_buf_q3), alpha_q12,
                                      alpha_sign, dc_q0);
      _mm_storel_epi64((__m128i *)dst, highbd_clamp_epi16(res, zero, max));
    } else {
      cfl_predict_hbd((__m128i *)dst, (__m128i *)pred_buf_q3, alpha_q12,
                      alpha_sign, dc_q0, zero, max);
    }
    if (width >= 16)
      cfl_predict_hbd((__m128i *)(dst + 8), (__m128i *)(pred_buf_q3 + 8),
                      alpha_q12, alpha_sign, dc_q0, zero, max);
    if (width == 32) {
      cfl_predict_hbd((__m128i *)(dst + 16), (__m128i *)(pred_buf_q3 + 16),
                      alpha_q12, alpha_sign, dc_q0, zero, max);
      cfl_predict_hbd((__m128i *)(dst + 24), (__m128i *)(pred_buf_q3 + 24),
                      alpha_q12, alpha_sign, dc_q0, zero, max);
    }
    dst += dst_stride;
    pred_buf_q3 += CFL_BUF_LINE;
  } while (dst < row_end);
}

#define CFL_PREDICT_LBD_X(width)                                               \
  static void cfl_predict_lbd_##width(const int16_t *pred_buf_q3,              \
                                      uint8_t *dst, int dst_stride,            \
                                      TX_SIZE tx_size, int alpha_q3) {         \
    cfl_predict_lbd_x(pred_buf_q3, dst, dst_stride, tx_size, alpha_q3, width); \
  }

CFL_PREDICT_LBD_X(4)
CFL_PREDICT_LBD_X(8)
CFL_PREDICT_LBD_X(16)
CFL_PREDICT_LBD_X(32)

#define CFL_PREDICT_HBD_X(width)                                               \
  static void cfl_predict_hbd_##width(const int16_t *pred_buf_q3,              \
                                      uint16_t *dst, int dst_stride,           \
                                      TX_SIZE tx_size, int alpha_q3, int bd) { \
    cfl_predict_hbd_x(pred_buf_q3, dst, dst_stride, tx_size, alpha_q3, bd,     \
                      width);                                                  \
  }

CFL_PREDICT_HBD_X(4)
CFL_PREDICT_HBD_X(8)
CFL_PREDICT_HBD_X(16)
CFL_PREDICT_HBD_X(32)

cfl_predict_lbd_fn get_predict_lbd_fn_ssse3(TX_SIZE tx_size) {
  static const cfl_predict_lbd_fn predict_lbd[4] = {
    cfl_predict_lbd_4, cfl_predict_lbd_8, cfl_predict_lbd_16, cfl_predict_lbd_32
  };
  return predict_lbd[(tx_size_wide_log2[tx_size] - tx_size_wide_log2[0]) & 3];
}

cfl_predict_hbd_fn get_predict_hbd_fn_ssse3(TX_SIZE tx_size) {
  static const cfl_predict_hbd_fn predict_hbd[4] = {
    cfl_predict_hbd_4, cfl_predict_hbd_8, cfl_predict_hbd_16, cfl_predict_hbd_32
  };
  return predict_hbd[(tx_size_wide_log2[tx_size] - tx_size_wide_log2[0]) & 3];
}

static INLINE int sum_epi32(__m128i l0) {
  l0 = _mm_add_epi32(l0, _mm_shuffle_epi32(l0, _MM_SHUFFLE(1, 0, 3, 2)));
  return _mm_cvtsi128_si32(
      _mm_add_epi32(l0, _mm_shuffle_epi32(l0, _MM_SHUFFLE(2, 3, 0, 1))));
}

int sum_block_4x4_ssse3(int16_t *pred_buf) {
  const __m128i zeros = _mm_setzero_si128();
  const __m128i *pred_buf_i128 = (__m128i *)pred_buf;
  const int stride_i128 = CFL_BUF_LINE >> 3;

  // At this stage in CfL, the maximum value in the CfL prediction buffer can
  // only reach 15 bit, so it is safe to do one addition inside 16 bit.
  // r0 = r0_0 + r1_0, r2_0 + r3_0, ..., r0_3 + r1_3, r2_3 + r3_3
  __m128i l0 = _mm_add_epi16(_mm_loadl_epi64(pred_buf_i128),
                             _mm_loadl_epi64(pred_buf_i128 + stride_i128));
  __m128i l1 = _mm_add_epi16(_mm_loadl_epi64(pred_buf_i128 + 2 * stride_i128),
                             _mm_loadl_epi64(pred_buf_i128 + 3 * stride_i128));

  // For the other additions, we need to convert to 32 bits.
  // To do so, we add the low part with the high part.
  // r0 = r0_0 + r0_4, r0_1 + r0_5, r0_2 + r0_6, r0_3 + r0_7
  l0 = _mm_add_epi32(_mm_unpacklo_epi16(l0, zeros),
                     _mm_unpacklo_epi16(l1, zeros));

  return sum_epi32(l0);
}

static __m128i load_epi16_add_epi32(const __m128i *pred_buf, int inc_si128) {
  const __m128i zeros = _mm_setzero_si128();
  const __m128i r0 = _mm_add_epi16(_mm_loadu_si128(pred_buf),
                                   _mm_loadu_si128(pred_buf + inc_si128));
  return _mm_add_epi32(_mm_unpacklo_epi16(r0, zeros),
                       _mm_unpackhi_epi16(r0, zeros));
}

// TODO(ltrudeau) const pred_buf
int sum_block_8x8_ssse3(int16_t *pred_buf) {
  const int stride_si128 = CFL_BUF_LINE >> 3;
  __m128i *pred_buf_si128 = (__m128i *)pred_buf;

  // (Row0 + Row1) + (Row2 + Row3)
  __m128i l0 = _mm_add_epi32(
      load_epi16_add_epi32(pred_buf_si128, stride_si128),
      load_epi16_add_epi32(pred_buf_si128 + 2 * stride_si128, stride_si128));

  // (Row4 + Row5) + (Row6 + Row7)
  __m128i l1 = _mm_add_epi32(
      load_epi16_add_epi32(pred_buf_si128 + 4 * stride_si128, stride_si128),
      load_epi16_add_epi32(pred_buf_si128 + 6 * stride_si128, stride_si128));

  // ((Row0 + Row1) + (Row2 + Row3)) + ((Row4 + Row5) + (Row6 + Row7))
  return sum_epi32(_mm_add_epi32(l0, l1));
}

int sum_block_16x16_ssse3(int16_t *pred_buf) {
  const __m128i *pred_buf_i128 = (__m128i *)pred_buf;
  const int stride_i128 = CFL_BUF_LINE >> 3;

  // At this stage in CfL, the maximum value in the CfL prediction buffer
  // can only reach 15 bit, so it is safe to do one addition inside 16 bit.
  // r0 = r0_0 + r1_0, r0_1 + r1_1, ..., r0_7 + r1_7
  __m128i l0 = load_epi16_add_epi32(pred_buf_i128, 1);
  __m128i l1 = load_epi16_add_epi32(pred_buf_i128 + stride_i128, 1);
  l0 = _mm_add_epi32(l0, l1);  // Row0 + Row1

  l1 = load_epi16_add_epi32(pred_buf_i128 + 2 * stride_i128, 1);
  __m128i l2 = load_epi16_add_epi32(pred_buf_i128 + 3 * stride_i128, 1);
  l1 = _mm_add_epi32(l1, l2);  // Row2 + Row3
  l0 = _mm_add_epi32(l0, l1);  // (Row0 + Row1) + (Row2 + Row3)

  l1 = load_epi16_add_epi32(pred_buf_i128 + 4 * stride_i128, 1);
  l2 = load_epi16_add_epi32(pred_buf_i128 + 5 * stride_i128, 1);
  l1 = _mm_add_epi32(l1, l2);  // Row4 + Row5

  l2 = load_epi16_add_epi32(pred_buf_i128 + 6 * stride_i128, 1);
  __m128i l3 = load_epi16_add_epi32(pred_buf_i128 + 7 * stride_i128, 1);
  l2 = _mm_add_epi32(l2, l3);  // Row6 + Row7
  l1 = _mm_add_epi32(l1, l2);  // (Row4 + Row5) + (Row6 + Row7)
  // ((Row0 + Row1) + (Row2 + Row3)) + ((Row4 + Row5) + (Row6 + Row7))
  l0 = _mm_add_epi32(l0, l1);

  l1 = load_epi16_add_epi32(pred_buf_i128 + 8 * stride_i128, 1);
  l2 = load_epi16_add_epi32(pred_buf_i128 + 9 * stride_i128, 1);
  l1 = _mm_add_epi32(l1, l2);  // Row8 + Row9

  l2 = load_epi16_add_epi32(pred_buf_i128 + 10 * stride_i128, 1);
  l3 = load_epi16_add_epi32(pred_buf_i128 + 11 * stride_i128, 1);
  l2 = _mm_add_epi32(l2, l3);  // Row10 + Row11
  l1 = _mm_add_epi32(l1, l2);  // (Row8 + Row9) + (Row10 + Row11)

  l2 = load_epi16_add_epi32(pred_buf_i128 + 12 * stride_i128, 1);
  l3 = load_epi16_add_epi32(pred_buf_i128 + 13 * stride_i128, 1);
  l2 = _mm_add_epi32(l2, l3);  // Row12 + Row13

  l3 = load_epi16_add_epi32(pred_buf_i128 + 14 * stride_i128, 1);
  __m128i l4 = load_epi16_add_epi32(pred_buf_i128 + 15 * stride_i128, 1);
  l3 = _mm_add_epi32(l3, l4);  // Row14 + Row15
  l2 = _mm_add_epi32(l2, l3);  // (Row12 + Row13) + (Row14 + Row15)
  // ((Row8 + Row9) + (Row10 + Row11)) + ((Row12 + Row13) + (Row14 + Row15))
  l1 = _mm_add_epi32(l1, l2);

  l0 = _mm_add_epi32(l0, l1);

  // r0 = r0_0 + r0_1, ...
  return sum_epi32(l0);
}

cfl_sum_block_fn get_sum_block_fn_ssse3(TX_SIZE tx_size) {
  static const cfl_sum_block_fn sum_block[TX_SIZES_ALL] = {
    sum_block_4x4_ssse3,    // 4x4
    sum_block_8x8_ssse3,    // 8x8
    sum_block_16x16_ssse3,  // 16x16
    cfl_sum_block_null,     // 32x32
#if CONFIG_TX64X64
    cfl_sum_block_null,  // 64x64 (invalid CFL size)
#endif                   // CONFIG_TX64X64
    cfl_sum_block_null,  // 4x8
    cfl_sum_block_null,  // 8x4
    cfl_sum_block_null,  // 8x16
    cfl_sum_block_null,  // 16x8
    cfl_sum_block_null,  // 16x32
    cfl_sum_block_null,  // 32x16
#if CONFIG_TX64X64
    cfl_sum_block_null,  // 32x64 (invalid CFL size)
    cfl_sum_block_null,  // 64x32 (invalid CFL size)
#endif                   // CONFIG_TX64X64
    cfl_sum_block_null,  // 4x16 (invalid CFL size)
    cfl_sum_block_null,  // 16x4 (invalid CFL size)
    cfl_sum_block_null,  // 8x32 (invalid CFL size)
    cfl_sum_block_null,  // 32x8 (invalid CFL size)
#if CONFIG_TX64X64
    cfl_sum_block_null,  // 16x64 (invalid CFL size)
    cfl_sum_block_null,  // 64x16 (invalid CFL size)
#endif                   // CONFIG_TX64X64
  };
  // Modulo TX_SIZES_ALL to ensure that an attacker won't be able to
  // index the function pointer array out of bounds.
  return sum_block[tx_size % TX_SIZES_ALL];
}
