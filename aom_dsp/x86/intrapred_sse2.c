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

#include "./aom_dsp_rtcd.h"

// -----------------------------------------------------------------------------
// DC_PRED

void aom_dc_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  __m128i sum_left = _mm_load_si128((const __m128i *)left);
  __m128i sum_above = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();

  sum_left = _mm_sad_epu8(sum_left, zero);
  sum_above = _mm_unpacklo_epi8(sum_above, zero);
  sum_above = _mm_sad_epu8(sum_above, zero);
  sum_above = _mm_add_epi16(sum_above, sum_left);
  uint32_t sum = _mm_cvtsi128_si32(sum_above);
  sum += 6;
  sum /= 12;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);
  const uint32_t pred = _mm_cvtsi128_si32(row);

  int i;
  for (i = 0; i < 2; ++i) {
    *(uint32_t *)dst = pred;
    dst += stride;
    *(uint32_t *)dst = pred;
    dst += stride;
    *(uint32_t *)dst = pred;
    dst += stride;
    *(uint32_t *)dst = pred;
    dst += stride;
  }
}

void aom_dc_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t stride,
                               const uint8_t *above, const uint8_t *left) {
  __m128i sum_left = _mm_load_si128((const __m128i *)left);
  __m128i sum_above = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();

  sum_above = _mm_sad_epu8(sum_above, zero);
  sum_left = _mm_unpacklo_epi8(sum_left, zero);
  sum_left = _mm_sad_epu8(sum_left, zero);
  sum_above = _mm_add_epi16(sum_above, sum_left);
  uint32_t sum = _mm_cvtsi128_si32(sum_above);
  sum += 6;
  sum /= 12;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);
  _mm_storel_epi64((__m128i *)dst, row);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row);
  dst += stride;
  _mm_storel_epi64((__m128i *)dst, row);
}

void aom_dc_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  __m128i sum_left = _mm_load_si128((const __m128i *)left);
  __m128i sum_above = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();

  sum_above = _mm_sad_epu8(sum_above, zero);
  sum_left = _mm_sad_epu8(sum_left, zero);
  const __m128i sum_left_high = _mm_unpackhi_epi64(sum_left, sum_left);
  sum_above = _mm_add_epi16(sum_above, sum_left);
  sum_above = _mm_add_epi16(sum_above, sum_left_high);

  uint32_t sum = _mm_cvtsi128_si32(sum_above);
  sum += 12;
  sum /= 24;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);

  int i;
  for (i = 0; i < 4; ++i) {
    _mm_storel_epi64((__m128i *)dst, row);
    dst += stride;
    _mm_storel_epi64((__m128i *)dst, row);
    dst += stride;
    _mm_storel_epi64((__m128i *)dst, row);
    dst += stride;
    _mm_storel_epi64((__m128i *)dst, row);
    dst += stride;
  }
}

void aom_dc_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  __m128i sum_left = _mm_load_si128((const __m128i *)left);
  __m128i sum_above = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();

  sum_above = _mm_sad_epu8(sum_above, zero);
  sum_left = _mm_sad_epu8(sum_left, zero);
  const __m128i sum_above_high = _mm_unpackhi_epi64(sum_above, sum_above);
  sum_above = _mm_add_epi16(sum_above, sum_left);
  sum_above = _mm_add_epi16(sum_above, sum_above_high);

  uint32_t sum = _mm_cvtsi128_si32(sum_above);
  sum += 12;
  sum /= 24;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);

  int i;
  for (i = 0; i < 2; ++i) {
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
  }
}

void aom_dc_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  __m128i sum_left0 = _mm_load_si128((const __m128i *)left);
  __m128i sum_left1 = _mm_load_si128((const __m128i *)(left + 16));
  __m128i sum_above = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();

  sum_above = _mm_sad_epu8(sum_above, zero);
  sum_left0 = _mm_sad_epu8(sum_left0, zero);
  sum_left1 = _mm_sad_epu8(sum_left1, zero);
  sum_left0 = _mm_add_epi64(sum_left0, sum_left1);
  sum_above = _mm_add_epi64(sum_left0, sum_above);

  __m128i sum_high = _mm_unpackhi_epi64(sum_above, sum_above);
  sum_above = _mm_add_epi32(sum_above, sum_high);

  uint32_t sum = _mm_cvtsi128_si32(sum_above);
  sum += 24;
  sum /= 48;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);

  int i;
  for (i = 0; i < 8; ++i) {
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
  }
}

void aom_dc_predictor_32x16_sse2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  __m128i sum_above0 = _mm_load_si128((const __m128i *)above);
  __m128i sum_above1 = _mm_load_si128((const __m128i *)(above + 16));
  __m128i sum_left = _mm_load_si128((const __m128i *)left);
  const __m128i zero = _mm_setzero_si128();

  sum_left = _mm_sad_epu8(sum_left, zero);
  sum_above0 = _mm_sad_epu8(sum_above0, zero);
  sum_above1 = _mm_sad_epu8(sum_above1, zero);

  sum_above0 = _mm_add_epi64(sum_above0, sum_above1);
  sum_above0 = _mm_add_epi64(sum_left, sum_above0);

  __m128i sum_high = _mm_unpackhi_epi64(sum_above0, sum_above0);
  sum_above0 = _mm_add_epi32(sum_above0, sum_high);

  uint32_t sum = _mm_cvtsi128_si32(sum_above0);
  sum += 24;
  sum /= 48;
  const __m128i row = _mm_set1_epi8((uint8_t)sum);

  int i;
  for (i = 0; i < 4; ++i) {
    _mm_store_si128((__m128i *)dst, row);
    _mm_store_si128((__m128i *)(dst + 16), row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    _mm_store_si128((__m128i *)(dst + 16), row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    _mm_store_si128((__m128i *)(dst + 16), row);
    dst += stride;
    _mm_store_si128((__m128i *)dst, row);
    _mm_store_si128((__m128i *)(dst + 16), row);
    dst += stride;
  }
}
