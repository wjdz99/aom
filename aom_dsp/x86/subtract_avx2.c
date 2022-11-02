/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <immintrin.h>

#include "config/aom_dsp_rtcd.h"

static INLINE void subtract_wd32_2rows_avx2(
    int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  const __m256i set_one_minusone = _mm256_set1_epi16((short)0xff01);

  const __m256i s_r0 = _mm256_lddqu_si256((__m256i *)(src_ptr));
  const __m256i p_r0 = _mm256_lddqu_si256((__m256i *)(pred_ptr));
  const __m256i s_r1 = _mm256_lddqu_si256((__m256i *)(src_ptr + src_stride));
  const __m256i p_r1 = _mm256_lddqu_si256((__m256i *)(pred_ptr + pred_stride));

  const __m256i s_r0_temp = _mm256_permute4x64_epi64(s_r0, 0xd8);
  const __m256i p_r0_temp = _mm256_permute4x64_epi64(p_r0, 0xd8);
  const __m256i s_r1_temp = _mm256_permute4x64_epi64(s_r1, 0xd8);
  const __m256i p_r1_temp = _mm256_permute4x64_epi64(p_r1, 0xd8);

  const __m256i r0_low_256 = _mm256_unpacklo_epi8(s_r0_temp, p_r0_temp);
  const __m256i r0_high_256 = _mm256_unpackhi_epi8(s_r0_temp, p_r0_temp);
  const __m256i r1_low_256 = _mm256_unpacklo_epi8(s_r1_temp, p_r1_temp);
  const __m256i r1_high_256 = _mm256_unpackhi_epi8(s_r1_temp, p_r1_temp);

  const __m256i d_r0_00 = _mm256_maddubs_epi16(r0_low_256, set_one_minusone);
  const __m256i d_r0_01 = _mm256_maddubs_epi16(r0_high_256, set_one_minusone);
  const __m256i d_r1_00 = _mm256_maddubs_epi16(r1_low_256, set_one_minusone);
  const __m256i d_r1_01 = _mm256_maddubs_epi16(r1_high_256, set_one_minusone);

  _mm256_store_si256((__m256i *)(diff_ptr), d_r0_00);
  _mm256_store_si256((__m256i *)(diff_ptr + 16), d_r0_01);
  _mm256_store_si256((__m256i *)(diff_ptr + diff_stride), d_r1_00);
  _mm256_store_si256((__m256i *)(diff_ptr + diff_stride + 16), d_r1_01);
}

static INLINE void subtract_block_16xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  for (int32_t j = 0; j < rows; j += 4) {
    const __m128i s_r0 = _mm_lddqu_si128((__m128i *)(src_ptr));
    const __m128i p_r0 = _mm_lddqu_si128((__m128i *)(pred_ptr));
    const __m128i s_r1 = _mm_lddqu_si128((__m128i *)(src_ptr += src_stride));
    const __m128i p_r1 = _mm_lddqu_si128((__m128i *)(pred_ptr += pred_stride));
    const __m128i s_r2 = _mm_lddqu_si128((__m128i *)(src_ptr += src_stride));
    const __m128i p_r2 = _mm_lddqu_si128((__m128i *)(pred_ptr += pred_stride));
    const __m128i s_r3 = _mm_lddqu_si128((__m128i *)(src_ptr += src_stride));
    const __m128i p_r3 = _mm_lddqu_si128((__m128i *)(pred_ptr += pred_stride));

    const __m256i s_0 = _mm256_cvtepu8_epi16(s_r0);
    const __m256i p_0 = _mm256_cvtepu8_epi16(p_r0);
    const __m256i s_1 = _mm256_cvtepu8_epi16(s_r1);
    const __m256i p_1 = _mm256_cvtepu8_epi16(p_r1);
    const __m256i s_2 = _mm256_cvtepu8_epi16(s_r2);
    const __m256i p_2 = _mm256_cvtepu8_epi16(p_r2);
    const __m256i s_3 = _mm256_cvtepu8_epi16(s_r3);
    const __m256i p_3 = _mm256_cvtepu8_epi16(p_r3);

    const __m256i d_0 = _mm256_sub_epi16(s_0, p_0);
    const __m256i d_1 = _mm256_sub_epi16(s_1, p_1);
    const __m256i d_2 = _mm256_sub_epi16(s_2, p_2);
    const __m256i d_3 = _mm256_sub_epi16(s_3, p_3);

    _mm256_store_si256((__m256i *)(diff_ptr), d_0);
    _mm256_store_si256((__m256i *)(diff_ptr += diff_stride), d_1);
    _mm256_store_si256((__m256i *)(diff_ptr += diff_stride), d_2);
    _mm256_store_si256((__m256i *)(diff_ptr += diff_stride), d_3);

    src_ptr += src_stride;
    pred_ptr += pred_stride;
    diff_ptr += diff_stride;
  }
}

static INLINE void subtract_block_32xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  for (int32_t j = 0; j < rows; j += 2) {
    subtract_wd32_2rows_avx2(diff_ptr, diff_stride, src_ptr, src_stride,
                             pred_ptr, pred_stride);
    src_ptr += 2 * src_stride;
    pred_ptr += 2 * pred_stride;
    diff_ptr += 2 * diff_stride;
  }
}

static INLINE void subtract_block_64xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  for (int32_t j = 0; j < rows; j += 2) {
    subtract_wd32_2rows_avx2(diff_ptr, diff_stride, src_ptr, src_stride,
                             pred_ptr, pred_stride);
    subtract_wd32_2rows_avx2(diff_ptr + 32, diff_stride, src_ptr + 32,
                             src_stride, pred_ptr + 32, pred_stride);
    src_ptr += 2 * src_stride;
    pred_ptr += 2 * pred_stride;
    diff_ptr += 2 * diff_stride;
  }
}

static INLINE void subtract_block_128xn_avx2(
    int rows, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr,
    ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride) {
  for (int32_t j = 0; j < rows; j += 2) {
    subtract_wd32_2rows_avx2(diff_ptr, diff_stride, src_ptr, src_stride,
                             pred_ptr, pred_stride);
    subtract_wd32_2rows_avx2(diff_ptr + 32, diff_stride, src_ptr + 32,
                             src_stride, pred_ptr + 32, pred_stride);
    subtract_wd32_2rows_avx2(diff_ptr + 64, diff_stride, src_ptr + 64,
                             src_stride, pred_ptr + 64, pred_stride);
    subtract_wd32_2rows_avx2(diff_ptr + 96, diff_stride, src_ptr + 96,
                             src_stride, pred_ptr + 96, pred_stride);

    src_ptr += 2 * src_stride;
    pred_ptr += 2 * pred_stride;
    diff_ptr += 2 * diff_stride;
  }
}

void aom_subtract_block_avx2(int rows, int cols, int16_t *diff_ptr,
                             ptrdiff_t diff_stride, const uint8_t *src_ptr,
                             ptrdiff_t src_stride, const uint8_t *pred_ptr,
                             ptrdiff_t pred_stride) {
  switch (cols) {
    case 16:
      subtract_block_16xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    case 32:
      subtract_block_32xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    case 64:
      subtract_block_64xn_avx2(rows, diff_ptr, diff_stride, src_ptr, src_stride,
                               pred_ptr, pred_stride);
      break;
    case 128:
      subtract_block_128xn_avx2(rows, diff_ptr, diff_stride, src_ptr,
                                src_stride, pred_ptr, pred_stride);
      break;
    default:
      aom_subtract_block_sse2(rows, cols, diff_ptr, diff_stride, src_ptr,
                              src_stride, pred_ptr, pred_stride);
      break;
  }
}
