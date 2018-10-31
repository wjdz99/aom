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

#include "config/aom_dsp_rtcd.h"
#include "aom_dsp/x86/lpf_common_sse2.h"

static INLINE __m256i dc_sum_64(const uint8_t *ref) {
  const __m256i x0 = _mm256_loadu_si256((const __m256i *)ref);
  const __m256i x1 = _mm256_loadu_si256((const __m256i *)(ref + 32));
  const __m256i zero = _mm256_setzero_si256();
  __m256i y0 = _mm256_sad_epu8(x0, zero);
  __m256i y1 = _mm256_sad_epu8(x1, zero);
  y0 = _mm256_add_epi64(y0, y1);
  __m256i u0 = _mm256_permute2x128_si256(y0, y0, 1);
  y0 = _mm256_add_epi64(u0, y0);
  u0 = _mm256_unpackhi_epi64(y0, y0);
  return _mm256_add_epi16(y0, u0);
}

static INLINE __m256i dc_sum_32(const uint8_t *ref) {
  const __m256i x = _mm256_loadu_si256((const __m256i *)ref);
  const __m256i zero = _mm256_setzero_si256();
  __m256i y = _mm256_sad_epu8(x, zero);
  __m256i u = _mm256_permute2x128_si256(y, y, 1);
  y = _mm256_add_epi64(u, y);
  u = _mm256_unpackhi_epi64(y, y);
  return _mm256_add_epi16(y, u);
}

static INLINE void row_store_32xh(const __m256i *r, int height, uint8_t *dst,
                                  ptrdiff_t stride) {
  for (int i = 0; i < height; ++i) {
    _mm256_storeu_si256((__m256i *)dst, *r);
    dst += stride;
  }
}

static INLINE void row_store_32x2xh(const __m256i *r0, const __m256i *r1,
                                    int height, uint8_t *dst,
                                    ptrdiff_t stride) {
  for (int i = 0; i < height; ++i) {
    _mm256_storeu_si256((__m256i *)dst, *r0);
    _mm256_storeu_si256((__m256i *)(dst + 32), *r1);
    dst += stride;
  }
}

static INLINE void row_store_64xh(const __m256i *r, int height, uint8_t *dst,
                                  ptrdiff_t stride) {
  for (int i = 0; i < height; ++i) {
    _mm256_storeu_si256((__m256i *)dst, *r);
    _mm256_storeu_si256((__m256i *)(dst + 32), *r);
    dst += stride;
  }
}

static INLINE void highbd_transpose16x4_8x8_sse2(__m128i *x, __m128i *d) {
  __m128i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;

  r0 = _mm_unpacklo_epi16(x[0], x[1]);
  r1 = _mm_unpacklo_epi16(x[2], x[3]);
  r2 = _mm_unpacklo_epi16(x[4], x[5]);
  r3 = _mm_unpacklo_epi16(x[6], x[7]);

  r4 = _mm_unpacklo_epi16(x[8], x[9]);
  r5 = _mm_unpacklo_epi16(x[10], x[11]);
  r6 = _mm_unpacklo_epi16(x[12], x[13]);
  r7 = _mm_unpacklo_epi16(x[14], x[15]);

  r8 = _mm_unpacklo_epi32(r0, r1);
  r9 = _mm_unpackhi_epi32(r0, r1);
  r10 = _mm_unpacklo_epi32(r2, r3);
  r11 = _mm_unpackhi_epi32(r2, r3);

  r12 = _mm_unpacklo_epi32(r4, r5);
  r13 = _mm_unpackhi_epi32(r4, r5);
  r14 = _mm_unpacklo_epi32(r6, r7);
  r15 = _mm_unpackhi_epi32(r6, r7);

  r0 = _mm_unpacklo_epi64(r8, r9);
  r1 = _mm_unpackhi_epi64(r8, r9);
  r2 = _mm_unpacklo_epi64(r10, r11);
  r3 = _mm_unpackhi_epi64(r10, r11);

  r4 = _mm_unpacklo_epi64(r12, r13);
  r5 = _mm_unpackhi_epi64(r12, r13);
  r6 = _mm_unpacklo_epi64(r14, r15);
  r7 = _mm_unpackhi_epi64(r14, r15);

  d[0] = _mm_unpacklo_epi64(r0, r2);
  d[1] = _mm_unpacklo_epi64(r4, r6);
  d[2] = _mm_unpacklo_epi64(r1, r3);
  d[3] = _mm_unpacklo_epi64(r5, r7);

  d[4] = _mm_unpackhi_epi64(r0, r2);
  d[5] = _mm_unpackhi_epi64(r4, r6);
  d[6] = _mm_unpackhi_epi64(r1, r3);
  d[7] = _mm_unpackhi_epi64(r5, r7);
}

static INLINE void highbd_transpose4x16_avx2(__m256i *x, __m256i *d) {
  __m256i w0, w1, w2, w3, ww0, ww1;

  w0 = _mm256_unpacklo_epi16(x[0], x[1]);  // 00 10 01 11 02 12 03 13
  w1 = _mm256_unpacklo_epi16(x[2], x[3]);  // 20 30 21 31 22 32 23 33
  w2 = _mm256_unpackhi_epi16(x[0], x[1]);  // 40 50 41 51 42 52 43 53
  w3 = _mm256_unpackhi_epi16(x[2], x[3]);  // 60 70 61 71 62 72 63 73

  ww0 = _mm256_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
  ww1 = _mm256_unpacklo_epi32(w2, w3);  // 40 50 60 70 41 51 61 71

  d[0] = _mm256_unpacklo_epi64(ww0, ww1);  // 00 10 20 30 40 50 60 70
  d[1] = _mm256_unpackhi_epi64(ww0, ww1);  // 01 11 21 31 41 51 61 71

  ww0 = _mm256_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
  ww1 = _mm256_unpackhi_epi32(w2, w3);  // 42 52 62 72 43 53 63 73

  d[2] = _mm256_unpacklo_epi64(ww0, ww1);  // 02 12 22 32 42 52 62 72
  d[3] = _mm256_unpackhi_epi64(ww0, ww1);  // 03 13 23 33 43 53 63 73
}

static INLINE void highbd_transpose8x16_16x8_avx2(__m256i *x, __m256i *d) {
  __m256i w0, w1, w2, w3, ww0, ww1;

  w0 = _mm256_unpacklo_epi16(x[0], x[1]);  // 00 10 01 11 02 12 03 13
  w1 = _mm256_unpacklo_epi16(x[2], x[3]);  // 20 30 21 31 22 32 23 33
  w2 = _mm256_unpacklo_epi16(x[4], x[5]);  // 40 50 41 51 42 52 43 53
  w3 = _mm256_unpacklo_epi16(x[6], x[7]);  // 60 70 61 71 62 72 63 73

  ww0 = _mm256_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
  ww1 = _mm256_unpacklo_epi32(w2, w3);  // 40 50 60 70 41 51 61 71

  d[0] = _mm256_unpacklo_epi64(ww0, ww1);  // 00 10 20 30 40 50 60 70
  d[1] = _mm256_unpackhi_epi64(ww0, ww1);  // 01 11 21 31 41 51 61 71

  ww0 = _mm256_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
  ww1 = _mm256_unpackhi_epi32(w2, w3);  // 42 52 62 72 43 53 63 73

  d[2] = _mm256_unpacklo_epi64(ww0, ww1);  // 02 12 22 32 42 52 62 72
  d[3] = _mm256_unpackhi_epi64(ww0, ww1);  // 03 13 23 33 43 53 63 73

  w0 = _mm256_unpackhi_epi16(x[0], x[1]);  // 04 14 05 15 06 16 07 17
  w1 = _mm256_unpackhi_epi16(x[2], x[3]);  // 24 34 25 35 26 36 27 37
  w2 = _mm256_unpackhi_epi16(x[4], x[5]);  // 44 54 45 55 46 56 47 57
  w3 = _mm256_unpackhi_epi16(x[6], x[7]);  // 64 74 65 75 66 76 67 77

  ww0 = _mm256_unpacklo_epi32(w0, w1);  // 04 14 24 34 05 15 25 35
  ww1 = _mm256_unpacklo_epi32(w2, w3);  // 44 54 64 74 45 55 65 75

  d[4] = _mm256_unpacklo_epi64(ww0, ww1);  // 04 14 24 34 44 54 64 74
  d[5] = _mm256_unpackhi_epi64(ww0, ww1);  // 05 15 25 35 45 55 65 75

  ww0 = _mm256_unpackhi_epi32(w0, w1);  // 06 16 26 36 07 17 27 37
  ww1 = _mm256_unpackhi_epi32(w2, w3);  // 46 56 66 76 47 57 67 77

  d[6] = _mm256_unpacklo_epi64(ww0, ww1);  // 06 16 26 36 46 56 66 76
  d[7] = _mm256_unpackhi_epi64(ww0, ww1);  // 07 17 27 37 47 57 67 77
}

static INLINE void highbd_transpose16x16_avx2(__m256i *x, __m256i *d) {
  __m256i w0, w1, w2, w3, ww0, ww1;
  __m256i dd[16];
  w0 = _mm256_unpacklo_epi16(x[0], x[1]);
  w1 = _mm256_unpacklo_epi16(x[2], x[3]);
  w2 = _mm256_unpacklo_epi16(x[4], x[5]);
  w3 = _mm256_unpacklo_epi16(x[6], x[7]);

  ww0 = _mm256_unpacklo_epi32(w0, w1);  //
  ww1 = _mm256_unpacklo_epi32(w2, w3);  //

  dd[0] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[1] = _mm256_unpackhi_epi64(ww0, ww1);

  ww0 = _mm256_unpackhi_epi32(w0, w1);  //
  ww1 = _mm256_unpackhi_epi32(w2, w3);  //

  dd[2] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[3] = _mm256_unpackhi_epi64(ww0, ww1);

  w0 = _mm256_unpackhi_epi16(x[0], x[1]);
  w1 = _mm256_unpackhi_epi16(x[2], x[3]);
  w2 = _mm256_unpackhi_epi16(x[4], x[5]);
  w3 = _mm256_unpackhi_epi16(x[6], x[7]);

  ww0 = _mm256_unpacklo_epi32(w0, w1);  //
  ww1 = _mm256_unpacklo_epi32(w2, w3);  //

  dd[4] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[5] = _mm256_unpackhi_epi64(ww0, ww1);

  ww0 = _mm256_unpackhi_epi32(w0, w1);  //
  ww1 = _mm256_unpackhi_epi32(w2, w3);  //

  dd[6] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[7] = _mm256_unpackhi_epi64(ww0, ww1);

  w0 = _mm256_unpacklo_epi16(x[8], x[9]);
  w1 = _mm256_unpacklo_epi16(x[10], x[11]);
  w2 = _mm256_unpacklo_epi16(x[12], x[13]);
  w3 = _mm256_unpacklo_epi16(x[14], x[15]);

  ww0 = _mm256_unpacklo_epi32(w0, w1);
  ww1 = _mm256_unpacklo_epi32(w2, w3);

  dd[8] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[9] = _mm256_unpackhi_epi64(ww0, ww1);

  ww0 = _mm256_unpackhi_epi32(w0, w1);
  ww1 = _mm256_unpackhi_epi32(w2, w3);

  dd[10] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[11] = _mm256_unpackhi_epi64(ww0, ww1);

  w0 = _mm256_unpackhi_epi16(x[8], x[9]);
  w1 = _mm256_unpackhi_epi16(x[10], x[11]);
  w2 = _mm256_unpackhi_epi16(x[12], x[13]);
  w3 = _mm256_unpackhi_epi16(x[14], x[15]);

  ww0 = _mm256_unpacklo_epi32(w0, w1);
  ww1 = _mm256_unpacklo_epi32(w2, w3);

  dd[12] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[13] = _mm256_unpackhi_epi64(ww0, ww1);

  ww0 = _mm256_unpackhi_epi32(w0, w1);
  ww1 = _mm256_unpackhi_epi32(w2, w3);

  dd[14] = _mm256_unpacklo_epi64(ww0, ww1);
  dd[15] = _mm256_unpackhi_epi64(ww0, ww1);

  for (int i = 0; i < 8; i++) {
    d[i] = _mm256_insertf128_si256(dd[i], _mm256_castsi256_si128(dd[i + 8]), 1);
    d[i + 8] = _mm256_insertf128_si256(dd[i + 8],
                                       _mm256_extracti128_si256(dd[i], 1), 0);
  }
}

void aom_dc_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m256i sum_above = dc_sum_32(above);
  __m256i sum_left = dc_sum_32(left);
  sum_left = _mm256_add_epi16(sum_left, sum_above);
  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum_left = _mm256_add_epi16(sum_left, thirtytwo);
  sum_left = _mm256_srai_epi16(sum_left, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum_left, zero);
  row_store_32xh(&row, 32, dst, stride);
}

void aom_dc_top_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_32(above);
  (void)left;

  const __m256i sixteen = _mm256_set1_epi16(16);
  sum = _mm256_add_epi16(sum, sixteen);
  sum = _mm256_srai_epi16(sum, 5);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_32xh(&row, 32, dst, stride);
}

void aom_dc_left_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m256i sum = dc_sum_32(left);
  (void)above;

  const __m256i sixteen = _mm256_set1_epi16(16);
  sum = _mm256_add_epi16(sum, sixteen);
  sum = _mm256_srai_epi16(sum, 5);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_32xh(&row, 32, dst, stride);
}

void aom_dc_128_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_32xh(&row, 32, dst, stride);
}

void aom_v_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row = _mm256_loadu_si256((const __m256i *)above);
  (void)left;
  row_store_32xh(&row, 32, dst, stride);
}

// There are 32 rows togeter. This function does line:
// 0,1,2,3, and 16,17,18,19. The next call would do
// 4,5,6,7, and 20,21,22,23. So 4 times of calling
// would finish 32 rows.
static INLINE void h_predictor_32x8line(const __m256i *row, uint8_t *dst,
                                        ptrdiff_t stride) {
  __m256i t[4];
  __m256i m = _mm256_setzero_si256();
  const __m256i inc = _mm256_set1_epi8(4);
  int i;

  for (i = 0; i < 4; i++) {
    t[i] = _mm256_shuffle_epi8(*row, m);
    __m256i r0 = _mm256_permute2x128_si256(t[i], t[i], 0);
    __m256i r1 = _mm256_permute2x128_si256(t[i], t[i], 0x11);
    _mm256_storeu_si256((__m256i *)dst, r0);
    _mm256_storeu_si256((__m256i *)(dst + (stride << 4)), r1);
    dst += stride;
    m = _mm256_add_epi8(m, inc);
  }
}

void aom_h_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  (void)above;
  const __m256i left_col = _mm256_loadu_si256((__m256i const *)left);

  __m256i u = _mm256_unpacklo_epi8(left_col, left_col);

  __m256i v = _mm256_unpacklo_epi8(u, u);
  h_predictor_32x8line(&v, dst, stride);
  dst += stride << 2;

  v = _mm256_unpackhi_epi8(u, u);
  h_predictor_32x8line(&v, dst, stride);
  dst += stride << 2;

  u = _mm256_unpackhi_epi8(left_col, left_col);

  v = _mm256_unpacklo_epi8(u, u);
  h_predictor_32x8line(&v, dst, stride);
  dst += stride << 2;

  v = _mm256_unpackhi_epi8(u, u);
  h_predictor_32x8line(&v, dst, stride);
}

// -----------------------------------------------------------------------------
// Rectangle

// TODO(luoyi) The following two functions are shared with intrapred_sse2.c.
// Use a header file, intrapred_common_x86.h
static INLINE __m128i dc_sum_16_sse2(const uint8_t *ref) {
  __m128i x = _mm_load_si128((__m128i const *)ref);
  const __m128i zero = _mm_setzero_si128();
  x = _mm_sad_epu8(x, zero);
  const __m128i high = _mm_unpackhi_epi64(x, x);
  return _mm_add_epi16(x, high);
}

static INLINE __m128i dc_sum_32_sse2(const uint8_t *ref) {
  __m128i x0 = _mm_load_si128((__m128i const *)ref);
  __m128i x1 = _mm_load_si128((__m128i const *)(ref + 16));
  const __m128i zero = _mm_setzero_si128();
  x0 = _mm_sad_epu8(x0, zero);
  x1 = _mm_sad_epu8(x1, zero);
  x0 = _mm_add_epi16(x0, x1);
  const __m128i high = _mm_unpackhi_epi64(x0, x0);
  return _mm_add_epi16(x0, high);
}

void aom_dc_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m128i top_sum = dc_sum_32_sse2(above);
  __m128i left_sum = dc_sum_16_sse2(left);
  left_sum = _mm_add_epi16(top_sum, left_sum);
  uint16_t sum = _mm_cvtsi128_si32(left_sum);
  sum += 24;
  sum /= 48;
  const __m256i row = _mm256_set1_epi8((uint8_t)sum);
  row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m256i sum_above = dc_sum_32(above);
  __m256i sum_left = dc_sum_64(left);
  sum_left = _mm256_add_epi16(sum_left, sum_above);
  uint16_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
  sum += 48;
  sum /= 96;
  const __m256i row = _mm256_set1_epi8((uint8_t)sum);
  row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m256i sum_above = dc_sum_64(above);
  __m256i sum_left = dc_sum_64(left);
  sum_left = _mm256_add_epi16(sum_left, sum_above);
  uint16_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
  sum += 64;
  sum /= 128;
  const __m256i row = _mm256_set1_epi8((uint8_t)sum);
  row_store_64xh(&row, 64, dst, stride);
}

void aom_dc_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m256i sum_above = dc_sum_64(above);
  __m256i sum_left = dc_sum_32(left);
  sum_left = _mm256_add_epi16(sum_left, sum_above);
  uint16_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
  sum += 48;
  sum /= 96;
  const __m256i row = _mm256_set1_epi8((uint8_t)sum);
  row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                 const uint8_t *above, const uint8_t *left) {
  const __m256i sum_above = dc_sum_64(above);
  __m256i sum_left = _mm256_castsi128_si256(dc_sum_16_sse2(left));
  sum_left = _mm256_add_epi16(sum_left, sum_above);
  uint16_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
  sum += 40;
  sum /= 80;
  const __m256i row = _mm256_set1_epi8((uint8_t)sum);
  row_store_64xh(&row, 16, dst, stride);
}

void aom_dc_top_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_32(above);
  (void)left;

  const __m256i sixteen = _mm256_set1_epi16(16);
  sum = _mm256_add_epi16(sum, sixteen);
  sum = _mm256_srai_epi16(sum, 5);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_top_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_32(above);
  (void)left;

  const __m256i sixteen = _mm256_set1_epi16(16);
  sum = _mm256_add_epi16(sum, sixteen);
  sum = _mm256_srai_epi16(sum, 5);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_top_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_64(above);
  (void)left;

  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum = _mm256_add_epi16(sum, thirtytwo);
  sum = _mm256_srai_epi16(sum, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_64xh(&row, 64, dst, stride);
}

void aom_dc_top_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_64(above);
  (void)left;

  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum = _mm256_add_epi16(sum, thirtytwo);
  sum = _mm256_srai_epi16(sum, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_top_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  __m256i sum = dc_sum_64(above);
  (void)left;

  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum = _mm256_add_epi16(sum, thirtytwo);
  sum = _mm256_srai_epi16(sum, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_64xh(&row, 16, dst, stride);
}

void aom_dc_left_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m128i sum = dc_sum_16_sse2(left);
  (void)above;

  const __m128i eight = _mm_set1_epi16(8);
  sum = _mm_add_epi16(sum, eight);
  sum = _mm_srai_epi16(sum, 4);
  const __m128i zero = _mm_setzero_si128();
  const __m128i r = _mm_shuffle_epi8(sum, zero);
  const __m256i row = _mm256_inserti128_si256(_mm256_castsi128_si256(r), r, 1);
  row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_left_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m256i sum = dc_sum_64(left);
  (void)above;

  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum = _mm256_add_epi16(sum, thirtytwo);
  sum = _mm256_srai_epi16(sum, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_left_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m256i sum = dc_sum_64(left);
  (void)above;

  const __m256i thirtytwo = _mm256_set1_epi16(32);
  sum = _mm256_add_epi16(sum, thirtytwo);
  sum = _mm256_srai_epi16(sum, 6);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_64xh(&row, 64, dst, stride);
}

void aom_dc_left_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m256i sum = dc_sum_32(left);
  (void)above;

  const __m256i sixteen = _mm256_set1_epi16(16);
  sum = _mm256_add_epi16(sum, sixteen);
  sum = _mm256_srai_epi16(sum, 5);
  const __m256i zero = _mm256_setzero_si256();
  __m256i row = _mm256_shuffle_epi8(sum, zero);
  row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_left_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                      const uint8_t *above,
                                      const uint8_t *left) {
  __m128i sum = dc_sum_16_sse2(left);
  (void)above;

  const __m128i eight = _mm_set1_epi16(8);
  sum = _mm_add_epi16(sum, eight);
  sum = _mm_srai_epi16(sum, 4);
  const __m128i zero = _mm_setzero_si128();
  const __m128i r = _mm_shuffle_epi8(sum, zero);
  const __m256i row = _mm256_inserti128_si256(_mm256_castsi128_si256(r), r, 1);
  row_store_64xh(&row, 16, dst, stride);
}

void aom_dc_128_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_128_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_128_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_64xh(&row, 64, dst, stride);
}

void aom_dc_128_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_128_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  (void)above;
  (void)left;
  const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
  row_store_64xh(&row, 16, dst, stride);
}

void aom_v_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row = _mm256_loadu_si256((const __m256i *)above);
  (void)left;
  row_store_32xh(&row, 16, dst, stride);
}

void aom_v_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row = _mm256_loadu_si256((const __m256i *)above);
  (void)left;
  row_store_32xh(&row, 64, dst, stride);
}

void aom_v_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
  const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
  (void)left;
  row_store_32x2xh(&row0, &row1, 64, dst, stride);
}

void aom_v_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
  const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
  (void)left;
  row_store_32x2xh(&row0, &row1, 32, dst, stride);
}

void aom_v_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                const uint8_t *above, const uint8_t *left) {
  const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
  const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
  (void)left;
  row_store_32x2xh(&row0, &row1, 16, dst, stride);
}

// -----------------------------------------------------------------------------
// PAETH_PRED

// Return 16 16-bit pixels in one row (__m256i)
static INLINE __m256i paeth_pred(const __m256i *left, const __m256i *top,
                                 const __m256i *topleft) {
  const __m256i base =
      _mm256_sub_epi16(_mm256_add_epi16(*top, *left), *topleft);

  __m256i pl = _mm256_abs_epi16(_mm256_sub_epi16(base, *left));
  __m256i pt = _mm256_abs_epi16(_mm256_sub_epi16(base, *top));
  __m256i ptl = _mm256_abs_epi16(_mm256_sub_epi16(base, *topleft));

  __m256i mask1 = _mm256_cmpgt_epi16(pl, pt);
  mask1 = _mm256_or_si256(mask1, _mm256_cmpgt_epi16(pl, ptl));
  __m256i mask2 = _mm256_cmpgt_epi16(pt, ptl);

  pl = _mm256_andnot_si256(mask1, *left);

  ptl = _mm256_and_si256(mask2, *topleft);
  pt = _mm256_andnot_si256(mask2, *top);
  pt = _mm256_or_si256(pt, ptl);
  pt = _mm256_and_si256(mask1, pt);

  return _mm256_or_si256(pt, pl);
}

// Return 16 8-bit pixels in one row (__m128i)
static INLINE __m128i paeth_16x1_pred(const __m256i *left, const __m256i *top,
                                      const __m256i *topleft) {
  const __m256i p0 = paeth_pred(left, top, topleft);
  const __m256i p1 = _mm256_permute4x64_epi64(p0, 0xe);
  const __m256i p = _mm256_packus_epi16(p0, p1);
  return _mm256_castsi256_si128(p);
}

static INLINE __m256i get_top_vector(const uint8_t *above) {
  const __m128i x = _mm_load_si128((const __m128i *)above);
  const __m128i zero = _mm_setzero_si128();
  const __m128i t0 = _mm_unpacklo_epi8(x, zero);
  const __m128i t1 = _mm_unpackhi_epi8(x, zero);
  return _mm256_inserti128_si256(_mm256_castsi128_si256(t0), t1, 1);
}

void aom_paeth_predictor_16x8_avx2(uint8_t *dst, ptrdiff_t stride,
                                   const uint8_t *above, const uint8_t *left) {
  __m128i x = _mm_loadl_epi64((const __m128i *)left);
  const __m256i l = _mm256_inserti128_si256(_mm256_castsi128_si256(x), x, 1);
  const __m256i tl16 = _mm256_set1_epi16((uint16_t)above[-1]);
  __m256i rep = _mm256_set1_epi16(0x8000);
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i top = get_top_vector(above);

  int i;
  for (i = 0; i < 8; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);
    const __m128i row = paeth_16x1_pred(&l16, &top, &tl16);

    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

static INLINE __m256i get_left_vector(const uint8_t *left) {
  const __m128i x = _mm_load_si128((const __m128i *)left);
  return _mm256_inserti128_si256(_mm256_castsi128_si256(x), x, 1);
}

void aom_paeth_predictor_16x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i l = get_left_vector(left);
  const __m256i tl16 = _mm256_set1_epi16((uint16_t)above[-1]);
  __m256i rep = _mm256_set1_epi16(0x8000);
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i top = get_top_vector(above);

  int i;
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);
    const __m128i row = paeth_16x1_pred(&l16, &top, &tl16);

    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

void aom_paeth_predictor_16x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  __m256i l = get_left_vector(left);
  const __m256i tl16 = _mm256_set1_epi16((uint16_t)above[-1]);
  __m256i rep = _mm256_set1_epi16(0x8000);
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i top = get_top_vector(above);

  int i;
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);
    const __m128i row = paeth_16x1_pred(&l16, &top, &tl16);

    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }

  l = get_left_vector(left + 16);
  rep = _mm256_set1_epi16(0x8000);
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);
    const __m128i row = paeth_16x1_pred(&l16, &top, &tl16);

    _mm_store_si128((__m128i *)dst, row);
    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

void aom_paeth_predictor_16x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i tl16 = _mm256_set1_epi16((uint16_t)above[-1]);
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i top = get_top_vector(above);

  for (int j = 0; j < 4; ++j) {
    const __m256i l = get_left_vector(left + j * 16);
    __m256i rep = _mm256_set1_epi16(0x8000);
    for (int i = 0; i < 16; ++i) {
      const __m256i l16 = _mm256_shuffle_epi8(l, rep);
      const __m128i row = paeth_16x1_pred(&l16, &top, &tl16);

      _mm_store_si128((__m128i *)dst, row);
      dst += stride;
      rep = _mm256_add_epi16(rep, one);
    }
  }
}

// Return 32 8-bit pixels in one row (__m256i)
static INLINE __m256i paeth_32x1_pred(const __m256i *left, const __m256i *top0,
                                      const __m256i *top1,
                                      const __m256i *topleft) {
  __m256i p0 = paeth_pred(left, top0, topleft);
  __m256i p1 = _mm256_permute4x64_epi64(p0, 0xe);
  const __m256i x0 = _mm256_packus_epi16(p0, p1);

  p0 = paeth_pred(left, top1, topleft);
  p1 = _mm256_permute4x64_epi64(p0, 0xe);
  const __m256i x1 = _mm256_packus_epi16(p0, p1);

  return _mm256_permute2x128_si256(x0, x1, 0x20);
}

void aom_paeth_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i l = get_left_vector(left);
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  __m256i rep = _mm256_set1_epi16(0x8000);
  const __m256i one = _mm256_set1_epi16(1);

  int i;
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);

    const __m256i r = paeth_32x1_pred(&l16, &t0, &t1, &tl);

    _mm256_storeu_si256((__m256i *)dst, r);

    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

void aom_paeth_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  __m256i l = get_left_vector(left);
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  __m256i rep = _mm256_set1_epi16(0x8000);
  const __m256i one = _mm256_set1_epi16(1);

  int i;
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);

    const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
    const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);

    _mm_store_si128((__m128i *)dst, r0);
    _mm_store_si128((__m128i *)(dst + 16), r1);

    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }

  l = get_left_vector(left + 16);
  rep = _mm256_set1_epi16(0x8000);
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);

    const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
    const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);

    _mm_store_si128((__m128i *)dst, r0);
    _mm_store_si128((__m128i *)(dst + 16), r1);

    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

void aom_paeth_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  const __m256i one = _mm256_set1_epi16(1);

  int i, j;
  for (j = 0; j < 4; ++j) {
    const __m256i l = get_left_vector(left + j * 16);
    __m256i rep = _mm256_set1_epi16(0x8000);
    for (i = 0; i < 16; ++i) {
      const __m256i l16 = _mm256_shuffle_epi8(l, rep);

      const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
      const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);

      _mm_store_si128((__m128i *)dst, r0);
      _mm_store_si128((__m128i *)(dst + 16), r1);

      dst += stride;
      rep = _mm256_add_epi16(rep, one);
    }
  }
}

void aom_paeth_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i t2 = get_top_vector(above + 32);
  const __m256i t3 = get_top_vector(above + 48);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  const __m256i one = _mm256_set1_epi16(1);

  int i, j;
  for (j = 0; j < 2; ++j) {
    const __m256i l = get_left_vector(left + j * 16);
    __m256i rep = _mm256_set1_epi16(0x8000);
    for (i = 0; i < 16; ++i) {
      const __m256i l16 = _mm256_shuffle_epi8(l, rep);

      const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
      const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);
      const __m128i r2 = paeth_16x1_pred(&l16, &t2, &tl);
      const __m128i r3 = paeth_16x1_pred(&l16, &t3, &tl);

      _mm_store_si128((__m128i *)dst, r0);
      _mm_store_si128((__m128i *)(dst + 16), r1);
      _mm_store_si128((__m128i *)(dst + 32), r2);
      _mm_store_si128((__m128i *)(dst + 48), r3);

      dst += stride;
      rep = _mm256_add_epi16(rep, one);
    }
  }
}

void aom_paeth_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i t2 = get_top_vector(above + 32);
  const __m256i t3 = get_top_vector(above + 48);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  const __m256i one = _mm256_set1_epi16(1);

  int i, j;
  for (j = 0; j < 4; ++j) {
    const __m256i l = get_left_vector(left + j * 16);
    __m256i rep = _mm256_set1_epi16(0x8000);
    for (i = 0; i < 16; ++i) {
      const __m256i l16 = _mm256_shuffle_epi8(l, rep);

      const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
      const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);
      const __m128i r2 = paeth_16x1_pred(&l16, &t2, &tl);
      const __m128i r3 = paeth_16x1_pred(&l16, &t3, &tl);

      _mm_store_si128((__m128i *)dst, r0);
      _mm_store_si128((__m128i *)(dst + 16), r1);
      _mm_store_si128((__m128i *)(dst + 32), r2);
      _mm_store_si128((__m128i *)(dst + 48), r3);

      dst += stride;
      rep = _mm256_add_epi16(rep, one);
    }
  }
}

void aom_paeth_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
                                    const uint8_t *above, const uint8_t *left) {
  const __m256i t0 = get_top_vector(above);
  const __m256i t1 = get_top_vector(above + 16);
  const __m256i t2 = get_top_vector(above + 32);
  const __m256i t3 = get_top_vector(above + 48);
  const __m256i tl = _mm256_set1_epi16((uint16_t)above[-1]);
  const __m256i one = _mm256_set1_epi16(1);

  int i;
  const __m256i l = get_left_vector(left);
  __m256i rep = _mm256_set1_epi16(0x8000);
  for (i = 0; i < 16; ++i) {
    const __m256i l16 = _mm256_shuffle_epi8(l, rep);

    const __m128i r0 = paeth_16x1_pred(&l16, &t0, &tl);
    const __m128i r1 = paeth_16x1_pred(&l16, &t1, &tl);
    const __m128i r2 = paeth_16x1_pred(&l16, &t2, &tl);
    const __m128i r3 = paeth_16x1_pred(&l16, &t3, &tl);

    _mm_store_si128((__m128i *)dst, r0);
    _mm_store_si128((__m128i *)(dst + 16), r1);
    _mm_store_si128((__m128i *)(dst + 32), r2);
    _mm_store_si128((__m128i *)(dst + 48), r3);

    dst += stride;
    rep = _mm256_add_epi16(rep, one);
  }
}

#define PERM4x64(c0, c1, c2, c3) c0 + (c1 << 2) + (c2 << 4) + (c3 << 6)
#define PERM2x128(c0, c1) c0 + (c1 << 4)

static AOM_FORCE_INLINE void av1_highbd_dr_prediction_z1_4xN_internal_avx2(
    int N, __m128i *dst, const uint16_t *above, int upsample_above, int dx) {
  const int frac_bits = 6 - upsample_above;
  const int max_base_x = ((N + 4) - 1) << upsample_above;
  int x;
  // a assert(dx > 0);
  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be caluculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m256i a0, a1, a32, a16;
  __m256i diff;
  __m128i a_mbase_x, max_base_x128, base_inc128, mask128;

  a16 = _mm256_set1_epi32(16);
  a_mbase_x = _mm_set1_epi16(above[max_base_x]);
  max_base_x128 = _mm_set1_epi32(max_base_x);

  x = dx;
  for (int r = 0; r < N; r++) {
    __m256i b, res, shift;
    __m128i res1;

    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        dst[i] = a_mbase_x;  // save 4 values
      }
      return;
    }

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

    if (upsample_above) {
      a0 = _mm256_permutevar8x32_epi32(
          a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
      a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));
      base_inc128 = _mm_setr_epi32(base, base + 2, base + 4, base + 6);
      shift = _mm256_srli_epi32(
          _mm256_and_si256(
              _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
              _mm256_set1_epi32(0x3f)),
          1);
    } else {
      base_inc128 = _mm_setr_epi32(base, base + 1, base + 2, base + 3);
      shift = _mm256_srli_epi32(
          _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
    }

    diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

    b = _mm256_mullo_epi32(diff, shift);
    res = _mm256_add_epi32(a32, b);
    res = _mm256_srli_epi32(res, 5);

    res1 = _mm256_castsi256_si128(res);
    res1 = _mm_packus_epi32(res1, res1);

    mask128 = _mm_cmpgt_epi32(max_base_x128, base_inc128);
    mask128 = _mm_packs_epi32(mask128, mask128);  // goto 16 bit
    dst[r] = _mm_blendv_epi8(a_mbase_x, res1, mask128);
    x += dx;
  }
}

void av1_highbd_dr_prediction_z1_4xN_avx2(int N, uint16_t *dst,
                                          ptrdiff_t stride,
                                          const uint16_t *above,
                                          int upsample_above, int dx, int bd) {
  __m128i dstvec[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_4xN_internal_avx2(N, dstvec, above,
                                                upsample_above, dx);
  for (int i = 0; i < N; i++) {
    _mm_storel_epi64((__m128i *)(dst + stride * i), dstvec[i]);
  }
}

static AOM_FORCE_INLINE void av1_highbd_dr_prediction_z1_8xN_internal_avx2(
    int N, __m128i *dst, const uint16_t *above, int upsample_above, int dx) {
  const int frac_bits = 6 - upsample_above;
  const int max_base_x = ((8 + N) - 1) << upsample_above;

  int x;
  // a assert(dx > 0);
  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be caluculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m256i a0, a1, a0_1, a1_1, a32, a16;
  __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

  a16 = _mm256_set1_epi32(16);
  a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
  max_base_x256 = _mm256_set1_epi32(max_base_x);

  x = dx;
  for (int r = 0; r < N; r++) {
    __m256i b, res, res1, shift;

    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        dst[i] = _mm256_castsi256_si128(a_mbase_x);  // save 8 values
      }
      return;
    }

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

    if (upsample_above) {
      a0 = _mm256_permutevar8x32_epi32(
          a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
      a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));

      a0_1 =
          _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
      a0_1 = _mm256_permutevar8x32_epi32(
          a0_1, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
      a1_1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0_1, 1));

      a0 = _mm256_inserti128_si256(a0, _mm256_castsi256_si128(a0_1), 1);
      a1 = _mm256_inserti128_si256(a1, _mm256_castsi256_si128(a1_1), 1);
      base_inc256 =
          _mm256_setr_epi32(base, base + 2, base + 4, base + 6, base + 8,
                            base + 10, base + 12, base + 14);
      shift = _mm256_srli_epi32(
          _mm256_and_si256(
              _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
              _mm256_set1_epi32(0x3f)),
          1);
    } else {
      base_inc256 = _mm256_setr_epi32(base, base + 1, base + 2, base + 3,
                                      base + 4, base + 5, base + 6, base + 7);
      shift = _mm256_srli_epi32(
          _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
    }

    diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

    b = _mm256_mullo_epi32(diff, shift);
    res = _mm256_add_epi32(a32, b);
    res = _mm256_srli_epi32(res, 5);

    res1 = _mm256_packus_epi32(
        res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));

    mask256 = _mm256_cmpgt_epi32(max_base_x256, base_inc256);
    mask256 = _mm256_packs_epi32(
        mask256, _mm256_castsi128_si256(
                     _mm256_extracti128_si256(mask256, 1)));  // goto 16 bit
    res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
    dst[r] = _mm256_castsi256_si128(res1);
    x += dx;
  }
}

void av1_highbd_dr_prediction_z1_8xN_avx2(int N, uint16_t *dst,
                                          ptrdiff_t stride,
                                          const uint16_t *above,
                                          int upsample_above, int dx, int bd) {
  __m128i dstvec[32];
  (void)bd;
  av1_highbd_dr_prediction_z1_8xN_internal_avx2(N, dstvec, above,
                                                upsample_above, dx);
  for (int i = 0; i < N; i++) {
    _mm_storeu_si128((__m128i *)(dst + stride * i), dstvec[i]);
  }
}

static AOM_FORCE_INLINE void av1_highbd_dr_prediction_z1_16xN_internal_avx2(
    int N, __m256i *dstvec, const uint16_t *above, int upsample_above, int dx) {
  int x;
  // here upsample_above is 0 by design of av1_use_intra_edge_upsample
  (void)upsample_above;
  const int frac_bits = 6;
  const int max_base_x = ((16 + N) - 1);

  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be caluculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m256i a0, a0_1, a1, a1_1, a32, a16;
  __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

  a16 = _mm256_set1_epi32(16);
  a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
  max_base_x256 = _mm256_set1_epi16(max_base_x);

  x = dx;
  for (int r = 0; r < N; r++) {
    __m256i b, res[2], res1;

    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        dstvec[i] = a_mbase_x;  // save 16 values
      }
      return;
    }

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
    a0_1 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));
    a1_1 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 9)));
    base_inc256 = _mm256_setr_epi16(base, base + 1, base + 2, base + 3,
                                    base + 4, base + 5, base + 6, base + 7,
                                    base + 8, base + 9, base + 10, base + 11,
                                    base + 12, base + 13, base + 14, base + 15);

    __m256i shift = _mm256_srli_epi32(
        _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

    diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
    b = _mm256_mullo_epi32(diff, shift);

    res[0] = _mm256_add_epi32(a32, b);
    res[0] = _mm256_srli_epi32(res[0], 5);

    diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
    b = _mm256_mullo_epi32(diff, shift);

    res[1] = _mm256_add_epi32(a32, b);
    res[1] = _mm256_srli_epi32(res[1], 5);

    res[0] = _mm256_packus_epi32(
        res[0], _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
    res[1] = _mm256_packus_epi32(
        res[1], _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));

    res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
                                   1);  // 16 16bit values

    mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);  
    dstvec[r] = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
    x += dx;
  }
}

void av1_highbd_dr_prediction_z1_16xN_avx2(int N, uint16_t *dst,
                                           ptrdiff_t stride,
                                           const uint16_t *above,
                                           int upsample_above, int dx, int bd) {
  __m256i dstvec[64];
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(N, dstvec, above,
                                                 upsample_above, dx);
  for (int i = 0; i < N; i++) {
    _mm256_storeu_si256((__m256i *)(dst + stride * i), dstvec[i]);
  }
}

static AOM_FORCE_INLINE void av1_highbd_dr_prediction_z1_32xN_internal_avx2(
    int N, __m256i *dstvec, const uint16_t *above, int upsample_above, int dx) {
  int x;
  // here upsample_above is 0 by design of av1_use_intra_edge_upsample
  (void)upsample_above;
  const int frac_bits = 6;
  const int max_base_x = ((32 + N) - 1);

  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be caluculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m256i a0, a0_1, a1, a1_1, a32, a16;
  __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

  a16 = _mm256_set1_epi32(16);
  a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
  max_base_x256 = _mm256_set1_epi16(max_base_x);

  x = dx;
  for (int r = 0; r < N; r++) {
    __m256i b, res[2], res1;

    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        dstvec[i] = a_mbase_x;  // save 32 values
        dstvec[i + N] = a_mbase_x;
      }
      return;
    }

    __m256i shift = _mm256_srli_epi32(
        _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

    base_inc256 = _mm256_setr_epi16(base, base + 1, base + 2, base + 3,
                                    base + 4, base + 5, base + 6, base + 7,
                                    base + 8, base + 9, base + 10, base + 11,
                                    base + 12, base + 13, base + 14, base + 15);
    for (int j = 0; j < 32; j += 16) {
      a0 =
          _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + j)));
      a0_1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 8 + j)));
      a1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 1 + j)));
      a1_1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 9 + j)));

      diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
      a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
      a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
      b = _mm256_mullo_epi32(diff, shift);

      res[0] = _mm256_add_epi32(a32, b);
      res[0] = _mm256_srli_epi32(res[0], 5);

      diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
      a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
      a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
      b = _mm256_mullo_epi32(diff, shift);

      res[1] = _mm256_add_epi32(a32, b);
      res[1] = _mm256_srli_epi32(res[1], 5);

      res[0] = _mm256_packus_epi32(
          res[0], _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
      res[1] = _mm256_packus_epi32(
          res[1], _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));

      res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
                                     1);  // 16 16bit values
      base_inc256 = _mm256_setr_epi16(
          base + j, base + j + 1, base + j + 2, base + j + 3, base + j + 4,
          base + j + 5, base + j + 6, base + j + 7, base + j + 8, base + j + 9,
          base + j + 10, base + j + 11, base + j + 12, base + j + 13,
          base + j + 14, base + j + 15);

      mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);
      res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
      if (!j)
        dstvec[r] = res1;
      else
        dstvec[r + N] = res1;
    }
    x += dx;
  }
}

void av1_highbd_dr_prediction_z1_32xN_avx2(int N, uint16_t *dst,
                                           ptrdiff_t stride,
                                           const uint16_t *above,
                                           int upsample_above, int dx, int bd) {
  __m256i dstvec[128];
  (void)bd;
  av1_highbd_dr_prediction_z1_32xN_internal_avx2(N, dstvec, above,
                                                 upsample_above, dx);
  for (int i = 0; i < N; i++) {
    _mm256_storeu_si256((__m256i *)(dst + stride * i), dstvec[i]);
    _mm256_storeu_si256((__m256i *)(dst + stride * i + 16), dstvec[i + N]);
  }
}

void av1_highbd_dr_prediction_z1_64xN_avx2(int N, uint16_t *dst,
                                           ptrdiff_t stride,
                                           const uint16_t *above,
                                           int upsample_above, int dx, int bd) {
  int x;
  (void)bd;
  // here upsample_above is 0 by design of av1_use_intra_edge_upsample
  (void)upsample_above;
  const int frac_bits = 6;
  const int max_base_x = ((64 + N) - 1);

  // pre-filter above pixels
  // store in temp buffers:
  //   above[x] * 32 + 16
  //   above[x+1] - above[x]
  // final pixels will be caluculated as:
  //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
  __m256i a0, a0_1, a1, a1_1, a32, a16;
  __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

  a16 = _mm256_set1_epi32(16);
  a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
  max_base_x256 = _mm256_set1_epi16(max_base_x);

  x = dx;
  for (int r = 0; r < N; r++, dst += stride) {
    __m256i b, res[2], res1;

    int base = x >> frac_bits;
    if (base >= max_base_x) {
      for (int i = r; i < N; ++i) {
        _mm256_storeu_si256((__m256i *)dst, a_mbase_x);  // save 32 values
        _mm256_storeu_si256((__m256i *)(dst + 16), a_mbase_x);
        _mm256_storeu_si256((__m256i *)(dst + 32), a_mbase_x);
        _mm256_storeu_si256((__m256i *)(dst + 48), a_mbase_x);
        dst += stride;
      }
      return;
    }

    __m256i shift = _mm256_srli_epi32(
        _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

    base_inc256 = _mm256_setr_epi16(base, base + 1, base + 2, base + 3,
                                    base + 4, base + 5, base + 6, base + 7,
                                    base + 8, base + 9, base + 10, base + 11,
                                    base + 12, base + 13, base + 14, base + 15);
    for (int j = 0; j < 64; j += 16) {
      a0 =
          _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + j)));
      a0_1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 8 + j)));
      a1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 1 + j)));
      a1_1 = _mm256_cvtepu16_epi32(
          _mm_loadu_si128((__m128i *)(above + base + 9 + j)));

      diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
      a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
      a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
      b = _mm256_mullo_epi32(diff, shift);

      res[0] = _mm256_add_epi32(a32, b);
      res[0] = _mm256_srli_epi32(res[0], 5);

      diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
      a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
      a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
      b = _mm256_mullo_epi32(diff, shift);

      res[1] = _mm256_add_epi32(a32, b);
      res[1] = _mm256_srli_epi32(res[1], 5);

      res[0] = _mm256_packus_epi32(
          res[0], _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
      res[1] = _mm256_packus_epi32(
          res[1], _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));

      res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
                                     1);  // 16 16bit values
      base_inc256 = _mm256_setr_epi16(
          base + j, base + j + 1, base + j + 2, base + j + 3, base + j + 4,
          base + j + 5, base + j + 6, base + j + 7, base + j + 8, base + j + 9,
          base + j + 10, base + j + 11, base + j + 12, base + j + 13,
          base + j + 14, base + j + 15);

      mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);
      res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
      _mm256_storeu_si256((__m256i *)(dst + j), res1);
    }
    x += dx;
  }
}

// Directional prediction, zone 1: 0 < angle < 90
void av1_highbd_dr_prediction_z1_avx2(uint16_t *dst, ptrdiff_t stride, int bw,
                                      int bh, const uint16_t *above,
                                      const uint16_t *left, int upsample_above,
                                      int dx, int dy, int bd) {
  (void)left;
  (void)dy;
  (void)bd;

  switch (bw) {
    case 4:
      av1_highbd_dr_prediction_z1_4xN_avx2(bh, dst, stride, above,
                                           upsample_above, dx, bd);
      break;
    case 8:
      av1_highbd_dr_prediction_z1_8xN_avx2(bh, dst, stride, above,
                                           upsample_above, dx, bd);
      break;
    case 16:
      av1_highbd_dr_prediction_z1_16xN_avx2(bh, dst, stride, above,
                                            upsample_above, dx, bd);
      break;
    case 32:
      av1_highbd_dr_prediction_z1_32xN_avx2(bh, dst, stride, above,
                                            upsample_above, dx, bd);
      break;
    case 64:
      av1_highbd_dr_prediction_z1_64xN_avx2(bh, dst, stride, above,
                                            upsample_above, dx, bd);
      break;
    default: break;
  }
  return;
}


uint16_t z2BlendMaskabove[] = {
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF
};
uint16_t z2BlendMaskleft[] = {
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
  0xFFFF, 0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0,      0,      0,      0,      0,      0,      0,      0,
  0,      0
};

void transpose_TX_8X8(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                      ptrdiff_t pitchDst) {
  __m128i r0, r1, r2, r3, r4, r5, r6, r7, r0_Lo, r1_Lo, r2_Lo, r3_Lo, r4_Lo,
      r5_Lo, r6_Lo;
  r0 = _mm_load_si128(
      (__m128i *)(src + 0 * pitchSrc));  // 07,06,05,04,03,02,01,00
  r1 = _mm_load_si128(
      (__m128i *)(src + 1 * pitchSrc));  // 17,16,15,14,13,12,11,10
  r2 = _mm_load_si128(
      (__m128i *)(src + 2 * pitchSrc));  // 27,26,25,24,23,22,21,20
  r3 = _mm_load_si128(
      (__m128i *)(src + 3 * pitchSrc));  // 37,36,35,34,33,32,31,30
  r4 = _mm_load_si128(
      (__m128i *)(src + 4 * pitchSrc));  // 47,46,45,44,43,42,41,40
  r5 = _mm_load_si128(
      (__m128i *)(src + 5 * pitchSrc));  // 57,56,55,54,53,52,51,50
  r6 = _mm_load_si128(
      (__m128i *)(src + 6 * pitchSrc));  // 67,66,65,64,63,62,61,60
  r7 = _mm_load_si128(
      (__m128i *)(src + 7 * pitchSrc));  // 77,76,75,74,73,72,71,70

  r0_Lo = _mm_unpacklo_epi16(r0, r1);
  r2_Lo = _mm_unpacklo_epi16(r2, r3);
  r4_Lo = _mm_unpacklo_epi16(r4, r5);
  r6_Lo = _mm_unpacklo_epi16(r6, r7);

  r1_Lo = r0_Lo;
  r0_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
  r1_Lo = _mm_unpackhi_epi32(r1_Lo, r2_Lo);
  r5_Lo = r4_Lo;
  r4_Lo = _mm_unpacklo_epi32(r4_Lo, r6_Lo);
  r5_Lo = _mm_unpackhi_epi32(r5_Lo, r6_Lo);
  r2_Lo = r0_Lo;
  r0_Lo = _mm_unpacklo_epi64(r0_Lo, r4_Lo);  // 64
  r2_Lo = _mm_unpackhi_epi64(r2_Lo, r4_Lo);
  r3_Lo = r1_Lo;
  r1_Lo = _mm_unpacklo_epi64(r1_Lo, r5_Lo);
  r3_Lo = _mm_unpackhi_epi64(r3_Lo, r5_Lo);

  _mm_storeu_si128((__m128i *)(dst + 0 * pitchDst), r0_Lo);
  _mm_storeu_si128((__m128i *)(dst + 1 * pitchDst), r2_Lo);
  _mm_storeu_si128((__m128i *)(dst + 2 * pitchDst), r1_Lo);
  _mm_storeu_si128((__m128i *)(dst + 3 * pitchDst), r3_Lo);

  r0 = _mm_unpackhi_epi16(r0, r1);
  r2 = _mm_unpackhi_epi16(r2, r3);
  r4 = _mm_unpackhi_epi16(r4, r5);
  r6 = _mm_unpackhi_epi16(r6, r7);

  r1 = r0;
  r0 = _mm_unpacklo_epi32(r0, r2);
  r1 = _mm_unpackhi_epi32(r1, r2);
  r5 = r4;
  r4 = _mm_unpacklo_epi32(r4, r6);
  r5 = _mm_unpackhi_epi32(r5, r6);
  r2 = r0;
  r0 = _mm_unpacklo_epi64(r0, r4);
  r2 = _mm_unpackhi_epi64(r2, r4);
  r3 = r1;
  r1 = _mm_unpacklo_epi64(r1, r5);
  r3 = _mm_unpackhi_epi64(r3, r5);

  _mm_storeu_si128((__m128i *)(dst + 4 * pitchDst), r0);
  _mm_storeu_si128((__m128i *)(dst + 5 * pitchDst), r2);
  _mm_storeu_si128((__m128i *)(dst + 6 * pitchDst), r1);
  _mm_storeu_si128((__m128i *)(dst + 7 * pitchDst), r3);
}

void transpose_TX_4X8(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                      ptrdiff_t pitchDst) {
  __m128i r0, r1, r2, r3, r0_Lo, r1_Lo, r2_Lo, r3_Lo;
  r0 = _mm_load_si128(
      (__m128i *)(src + 0 * pitchSrc));  // 07,06,05,04,03,02,01,00
  r1 = _mm_load_si128(
      (__m128i *)(src + 1 * pitchSrc));  // 17,16,15,14,13,12,11,10
  r2 = _mm_load_si128(
      (__m128i *)(src + 2 * pitchSrc));  // 27,26,25,24,23,22,21,20
  r3 = _mm_load_si128(
      (__m128i *)(src + 3 * pitchSrc));  // 37,36,35,34,33,32,31,30

  r0_Lo = _mm_unpacklo_epi16(r0, r1);
  r2_Lo = _mm_unpacklo_epi16(r2, r3);

  r1_Lo = r0_Lo;
  r0_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
  r1_Lo = _mm_unpackhi_epi32(r1_Lo, r2_Lo);
  r2_Lo = _mm_srli_si128(r0_Lo, 8);
  r3_Lo = _mm_srli_si128(r1_Lo, 8);

  _mm_storel_epi64((__m128i *)(dst + 0 * pitchDst), r0_Lo);
  _mm_storel_epi64((__m128i *)(dst + 1 * pitchDst), r2_Lo);
  _mm_storel_epi64((__m128i *)(dst + 2 * pitchDst), r1_Lo);
  _mm_storel_epi64((__m128i *)(dst + 3 * pitchDst), r3_Lo);

  r0 = _mm_unpackhi_epi16(r0, r1);
  r2 = _mm_unpackhi_epi16(r2, r3);

  r1 = r0;
  r0 = _mm_unpacklo_epi32(r0, r2);
  r1 = _mm_unpackhi_epi32(r1, r2);
  r2 = _mm_srli_si128(r0, 8);
  r3 = _mm_srli_si128(r1, 8);

  _mm_storel_epi64((__m128i *)(dst + 4 * pitchDst), r0);
  _mm_storel_epi64((__m128i *)(dst + 5 * pitchDst), r2);
  _mm_storel_epi64((__m128i *)(dst + 6 * pitchDst), r1);
  _mm_storel_epi64((__m128i *)(dst + 7 * pitchDst), r3);
}
void transpose(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
               ptrdiff_t pitchDst, int width, int height) {
  for (int j = 0; j < height; j += 8)
    for (int i = 0; i < width; i += 8)
      transpose_TX_8X8(src + i * pitchSrc + j, pitchSrc, dst + j * pitchDst + i,
                       pitchDst);
}

void predict_64x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[128];
  int_least32_t leftByDiff[128];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

  for (c = 0; c < 64; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
  }
}
void predict_64x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[128];
  int_least32_t leftByDiff[128];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

  for (c = 0; c < 32; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
  }
}
void predict_64x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[128];
  int_least32_t leftByDiff[128];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

  for (c = 0; c < 16; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
  }
}
void predict_32x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[64];
  int_least32_t leftByDiff[64];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

  for (c = 0; c < 64; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
  }
}
void predict_32x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[64];
  int_least32_t leftByDiff[64];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

  for (c = 0; c < 32; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
  }
}
void predict_32x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[64];
  int_least32_t leftByDiff[64];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

  for (c = 0; c < 16; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
  }
}
void predict_32x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  uint32_t leftBy32[64];
  int_least32_t leftByDiff[64];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

  for (c = 0; c < 8; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, res3, res4, resLo, resHi;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
    b = _mm256_mullo_epi32(b, shift);
    res3 = _mm256_add_epi32(a, b);
    res3 = _mm256_srli_epi32(res3, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
    b = _mm256_mullo_epi32(b, shift);
    res4 = _mm256_add_epi32(a, b);
    res4 = _mm256_srli_epi32(res4, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
  }
}
void predict_16x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[32];
  int_least32_t leftByDiff[32];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

  for (c = 0; c < 16; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, resLo;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
  }
}
void predict_16x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[32];
  int_least32_t leftByDiff[32];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

  for (c = 0; c < 64; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, resLo;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
  }
}
void predict_16x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                        int d, int upsample) {
  uint32_t leftBy32[32];
  int_least32_t leftByDiff[32];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

  for (c = 0; c < 32; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, resLo;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
  }
}
void predict_16x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  uint32_t leftBy32[32];
  int_least32_t leftByDiff[32];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

  for (c = 0; c < 8; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, resLo;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
  }
}
void predict_16x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  uint32_t leftBy32[32];
  int_least32_t leftByDiff[32];
  __m256i a0, a1, diff, a32, a16;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

  a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
  a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
  diff = _mm256_sub_epi32(a1, a0);
  a32 = _mm256_slli_epi32(a0, 5);
  a32 = _mm256_add_epi32(a32, a16);
  _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

  for (c = 0; c < 4; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    __m256i a, b, res1, res2, resLo;
    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
    b = _mm256_mullo_epi32(b, shift);
    res2 = _mm256_add_epi32(a, b);
    res2 = _mm256_srli_epi32(res2, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
  }
}
void predict_8x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  __m256i a0, a1, diff, a32, a16;
  __m256i a, b, res1, res2, resLo;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  res2 = _mm256_setzero_si256();
  uint32_t leftBy32[16];
  int_least32_t leftByDiff[16];

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

  for (c = 0; c < 32; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm_storeu_si128((__m128i *)(dst + c * pitch),
                     _mm256_castsi256_si128(resLo));
  }
}
void predict_8x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  __m256i a0, a1, diff, a32, a16;
  __m256i a, b, res1, res2, resLo;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  res2 = _mm256_setzero_si256();
  uint32_t leftBy32[16];
  int_least32_t leftByDiff[16];

  a16 = _mm256_set1_epi32(16);
  a0 =
      _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm256_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
  a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
  a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
  _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
  _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

  for (c = 0; c < 16; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m256i inc = _mm256_set1_epi32(y);
    __m256i shift =
        _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                           _mm256_set1_epi32(0x3f)),
                          1);

    base = base + 1;
    a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
    b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
    b = _mm256_mullo_epi32(b, shift);
    res1 = _mm256_add_epi32(a, b);
    res1 = _mm256_srli_epi32(res1, 5);

    resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                     PERM4x64(0, 2, 1, 3));
    _mm_storeu_si128((__m128i *)(dst + c * pitch),
                     _mm256_castsi256_si128(resLo));
  }
}
void predict_8x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                      int d, int upsample) {
  __m256i a0, a1, diff, a32, a16;
  __m256i a, b, res1, res2, resLo;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  if (upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    leftPtr = refPel - 2;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 8; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 2;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);
      res1 = _mm256_permutevar8x32_epi32(
          res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
      b = _mm256_mullo_epi32(b, shift);
      res2 = _mm256_add_epi32(a, b);
      res2 = _mm256_srli_epi32(res2, 5);
      res2 = _mm256_permutevar8x32_epi32(
          res2, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));
      resLo = _mm256_permute2x128_si256(res1, res2, 0x20);
      res1 = _mm256_permute2x128_si256(resLo, resLo, 0x01);
      res1 = _mm256_packus_epi32(resLo, res1);
      _mm_storeu_si128((__m128i *)(dst + c * pitch),
                       _mm256_castsi256_si128(res1));
    }
  } else {
    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 8; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 1;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);

      resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                       PERM4x64(0, 2, 1, 3));
      _mm_storeu_si128((__m128i *)(dst + c * pitch),
                       _mm256_castsi256_si128(resLo));
    }
  }
}
void predict_8x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                      int d, int upsample) {
  __m256i a0, a1, diff, a32, a16;
  __m256i a, b, res1, res2, resLo;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  if (upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    leftPtr = refPel - 2;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 4; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 2;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);
      res1 = _mm256_permutevar8x32_epi32(
          res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
      b = _mm256_mullo_epi32(b, shift);
      res2 = _mm256_add_epi32(a, b);
      res2 = _mm256_srli_epi32(res2, 5);
      res2 = _mm256_permutevar8x32_epi32(
          res2, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));
      resLo = _mm256_permute2x128_si256(res1, res2, 0x20);
      res1 = _mm256_permute2x128_si256(resLo, resLo, 0x01);
      res1 = _mm256_packus_epi32(resLo, res1);
      _mm_storeu_si128((__m128i *)(dst + c * pitch),
                       _mm256_castsi256_si128(res1));
    }
  } else {
    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 4; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 1;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);

      resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                                       PERM4x64(0, 2, 1, 3));
      _mm_storeu_si128((__m128i *)(dst + c * pitch),
                       _mm256_castsi256_si128(resLo));
    }
  }
}
void predict_4x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                      int d, int upsample) {
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  if (upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1;
    __m128i res;
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];
    leftPtr = refPel - 2;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 4; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 2;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);
      res1 = _mm256_permutevar8x32_epi32(
          res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

      res = _mm256_extractf128_si256(res1, 0);
      _mm_storel_epi64((__m128i *)(dst + c * pitch),
                       _mm_packus_epi32(res, res));
    }
  } else {
    uint32_t leftBy32[8];
    int_least32_t leftByDiff[8];
    __m128i a0, a1, diff, a32, a16;
    __m128i a, b, res1, res2, resLo;
    res2 = _mm_setzero_si128();

    a16 = _mm_set1_epi32(16);
    a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
    a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
    a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
    _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

    for (c = 0; c < 4; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m128i inc = _mm_set1_epi32(y);
      __m128i shift = _mm_srli_epi32(
          _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)),
          1);

      base = base + 1;
      a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
      b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
      b = _mm_mullo_epi32(b, shift);
      res1 = _mm_add_epi32(a, b);
      res1 = _mm_srli_epi32(res1, 5);

      resLo = _mm_packus_epi32(res1, res2);
      _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
    }
  }
}
void predict_4x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                      int d, int upsample) {
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  if (upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1;
    __m128i res;
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];
    leftPtr = refPel - 2;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 8; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m256i inc = _mm256_set1_epi32(y);
      __m256i shift =
          _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                                             _mm256_set1_epi32(0x3f)),
                            1);

      base = base + 2;
      a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
      b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
      b = _mm256_mullo_epi32(b, shift);
      res1 = _mm256_add_epi32(a, b);
      res1 = _mm256_srli_epi32(res1, 5);
      res1 = _mm256_permutevar8x32_epi32(
          res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

      res = _mm256_extractf128_si256(res1, 0);
      _mm_storel_epi64((__m128i *)(dst + c * pitch),
                       _mm_packus_epi32(res, res));
    }
  } else {
    uint32_t leftBy32[8];
    int_least32_t leftByDiff[8];
    __m128i a0, a1, diff, a32, a16;
    __m128i a, b, res1, res2, resLo;
    res2 = _mm_setzero_si128();

    a16 = _mm_set1_epi32(16);
    a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
    a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
    a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
    _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

    for (c = 0; c < 8; ++c) {
      x = c + 1;
      y = -x * d;
      base = y >> frac_bits;

      __m128i inc = _mm_set1_epi32(y);
      __m128i shift = _mm_srli_epi32(
          _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)),
          1);

      base = base + 1;
      a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
      b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
      b = _mm_mullo_epi32(b, shift);
      res1 = _mm_add_epi32(a, b);
      res1 = _mm_srli_epi32(res1, 5);

      resLo = _mm_packus_epi32(res1, res2);
      _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
    }
  }
}
void predict_4x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
                       int d, int upsample) {
  uint32_t leftBy32[8];
  int_least32_t leftByDiff[8];
  __m128i a0, a1, diff, a32, a16;
  __m128i a, b, res1, res2, resLo;
  int c, x, y, base;
  const uint16_t *leftPtr = refPel - 1;
  const int frac_bits = 6 - upsample;

  res2 = _mm_setzero_si128();

  a16 = _mm_set1_epi32(16);
  a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
  a1 = _mm_cvtepu16_epi32(
      _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
  diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
  a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
  a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
  _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
  _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

  for (c = 0; c < 16; ++c) {
    x = c + 1;
    y = -x * d;
    base = y >> frac_bits;

    __m128i inc = _mm_set1_epi32(y);
    __m128i shift = _mm_srli_epi32(
        _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)), 1);

    base = base + 1;
    a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
    b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
    b = _mm_mullo_epi32(b, shift);
    res1 = _mm_add_epi32(a, b);
    res1 = _mm_srli_epi32(res1, 5);

    resLo = _mm_packus_epi32(res1, res2);
    _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
  }
}
void transpose_TX_8X4(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                      ptrdiff_t pitchDst) {
  __m128i r0, r1, r2, r3, r4, r5, r6, r7, r0_Lo, r1_Lo, r2_Lo, r3_Lo, r4_Lo,
      r5_Lo, r6_Lo;
  r0 = _mm_loadu_si128((__m128i *)(src + 0 * pitchSrc));
  r1 = _mm_srli_si128(r0, 8);
  r2 = _mm_loadu_si128((__m128i *)(src + 2 * pitchSrc));
  r3 = _mm_srli_si128(r2, 8);
  r4 = _mm_loadu_si128((__m128i *)(src + 4 * pitchSrc));
  r5 = _mm_srli_si128(r4, 8);
  r6 = _mm_loadu_si128((__m128i *)(src + 6 * pitchSrc));
  r7 = _mm_srli_si128(r6, 8);

  r0_Lo = _mm_unpacklo_epi16(r0, r1);
  r2_Lo = _mm_unpacklo_epi16(r2, r3);
  r4_Lo = _mm_unpacklo_epi16(r4, r5);
  r6_Lo = _mm_unpacklo_epi16(r6, r7);

  r1_Lo = r0_Lo;
  r0_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
  r1_Lo = _mm_unpackhi_epi32(r1_Lo, r2_Lo);
  r5_Lo = r4_Lo;
  r4_Lo = _mm_unpacklo_epi32(r4_Lo, r6_Lo);
  r5_Lo = _mm_unpackhi_epi32(r5_Lo, r6_Lo);
  r2_Lo = r0_Lo;
  r0_Lo = _mm_unpacklo_epi64(r0_Lo, r4_Lo);
  r2_Lo = _mm_unpackhi_epi64(r2_Lo, r4_Lo);
  r3_Lo = r1_Lo;
  r1_Lo = _mm_unpacklo_epi64(r1_Lo, r5_Lo);
  r3_Lo = _mm_unpackhi_epi64(r3_Lo, r5_Lo);

  _mm_storeu_si128((__m128i *)(dst + 0 * pitchDst), r0_Lo);
  _mm_storeu_si128((__m128i *)(dst + 1 * pitchDst), r2_Lo);
  _mm_storeu_si128((__m128i *)(dst + 2 * pitchDst), r1_Lo);
  _mm_storeu_si128((__m128i *)(dst + 3 * pitchDst), r3_Lo);
}
void transpose_TX_4X4(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                      ptrdiff_t pitchDst) {
  assert(pitchSrc == 4);
  (void)pitchSrc;
  if (pitchDst == 4) {
    __m128i s = _mm_loadu_si128((__m128i *)src);
    __m128i r1 = _mm_srli_si128(s, 8);
    __m128i r2 = _mm_loadu_si128((__m128i *)(src + 8));
    __m128i r3 = _mm_srli_si128(r2, 8);

    __m128i r0_Lo = _mm_unpacklo_epi16(s, r1);
    __m128i r2_Lo = _mm_unpacklo_epi16(r2, r3);
    __m128i r1_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
    r0_Lo = _mm_unpackhi_epi32(r0_Lo, r2_Lo);

    _mm_storeu_si128((__m128i *)(dst + 0 * pitchDst), r1_Lo);
    _mm_storeu_si128((__m128i *)(dst + 2 * pitchDst), r0_Lo);
  } else {
    __m128i s = _mm_loadu_si128((__m128i *)src);
    __m128i r1 = _mm_srli_si128(s, 8);
    __m128i r2 = _mm_loadu_si128((__m128i *)(src + 8));
    __m128i r3 = _mm_srli_si128(r2, 8);

    __m128i r0_Lo = _mm_unpacklo_epi16(s, r1);
    __m128i r2_Lo = _mm_unpacklo_epi16(r2, r3);
    __m128i r1_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
    r0_Lo = _mm_unpackhi_epi32(r0_Lo, r2_Lo);

    _mm_storel_epi64((__m128i *)(dst + 0 * pitchDst), r1_Lo);
    _mm_storel_epi64((__m128i *)(dst + 1 * pitchDst), _mm_srli_si128(r1_Lo, 8));
    _mm_storel_epi64((__m128i *)(dst + 2 * pitchDst), r0_Lo);
    _mm_storel_epi64((__m128i *)(dst + 3 * pitchDst), _mm_srli_si128(r0_Lo, 8));
  }
}
void transposeMul4x8(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                     ptrdiff_t pitchDst, int width, int height) {
  for (int j = 0; j < height; j += 8)
    for (int i = 0; i < width; i += 4)
      transpose_TX_4X8(src + i * pitchSrc + j, pitchSrc, dst + j * pitchDst + i,
                       pitchDst);
}
void transposeMul8x4(const uint16_t *src, ptrdiff_t pitchSrc, uint16_t *dst,
                     ptrdiff_t pitchDst, int width, int height) {
  for (int j = 0; j < height; j += 4)
    for (int i = 0; i < width; i += 8)
      transpose_TX_8X4(src + i * pitchSrc + j, pitchSrc, dst + j * pitchDst + i,
                       pitchDst);
}

void av1_highbd_dr_prediction_z2_4x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left,
                                          int upsample_above, int upsample_left,
                                          int dx, int dy, int bd) {
  int r, x, base, base1, c_pos;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  const int frac_bits_x = 6 - upsample_above;

  uint16_t dstTmp[4 * 4] = { 0 };
  uint16_t dstTmpTransp[4 * 4] = { 0 };
  int strideTmp = 4;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

 predict_4x4_sse4(above, dst, stride, dx, upsample_above);
 predict_4x4_sse4(left, dstTmp, strideTmp, dy, upsample_left);
 transpose_TX_4X4(dstTmp, strideTmp, dstTmpTransp, strideTmp);

  for (r = 0; r < 4; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storel_epi64((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_8x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left,
                                          int upsample_above, int upsample_left,
                                          int dx, int dy, int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[8 * 8] = { 0 };
  uint16_t dstTmpTransp[8 * 8] = { 0 };
  int strideTmp = 8;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_8x8_sse4(above, dst, stride, dx, upsample_above);
  predict_8x8_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmp, 8, 8);

  for (r = 0; r < 8; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storeu_si128((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_4x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left,
                                          int upsample_above, int upsample_left,
                                          int dx, int dy, int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[8 * 4] = { 0 };
  uint16_t dstTmpTransp[4 * 8] = { 0 };
  int strideTmp = 8;
  int strideTmpTransp = 4;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_4x8_sse4(above, dst, stride, dx, upsample_above);
  predict_8x4_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose_TX_4X8(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp);

  for (r = 0; r < 8; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storel_epi64((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_8x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *above,
                                          const uint16_t *left,
                                          int upsample_above, int upsample_left,
                                          int dx, int dy, int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[4 * 8] = { 0 };
  uint16_t dstTmpTransp[8 * 4] = { 0 };
  int strideTmp = 4;
  int strideTmpTransp = 8;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_8x4_sse4(above, dst, stride, dx, upsample_above);
  predict_4x8_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose_TX_8X4(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp);

  for (r = 0; r < 4; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storeu_si128((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_8x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[16 * 8] = { 0 };
  uint16_t dstTmpTransp[8 * 16] = { 0 };
  int strideTmp = 16;
  int strideTmpTransp = 8;
  const int bw = tx_size_wide[TX_8X16];
  const int bh = tx_size_high[TX_8X16];
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_8x16_sse4(above, dst, stride, dx, upsample_above);
  predict_16x8_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, bw, bh);

  for (r = 0; r < 16; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storeu_si128((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_4x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[16 * 4] = { 0 };
  uint16_t dstTmpTransp[4 * 16] = { 0 };
  int strideTmp = 16;
  int strideTmpTransp = 4;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_4x16_sse4(above, dst, stride, dx, upsample_above);
  predict_16x4_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transposeMul4x8(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 4, 16);

  for (r = 0; r < 16; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storel_epi64((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_16x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[4 * 16] = { 0 };
  uint16_t dstTmpTransp[16 * 4] = { 0 };
  int strideTmp = 4;
  int strideTmpTransp = 16;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_16x4_sse4(above, dst, stride, dx, upsample_above);
  predict_4x16_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transposeMul8x4(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 16, 4);

  for (r = 0; r < 4; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }
    ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base));
    ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base));
    a = _mm256_loadu_si256((__m256i *)(dst + r * stride));
    l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm256_and_si256(a, ma);
    leftReg = _mm256_and_si256(l, ml);
    result = _mm256_or_si256(aboveReg, leftReg);
    result = _mm256_min_epu16(clip_bd, result);
    _mm256_storeu_si256((__m256i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_16x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[8 * 16] = { 0 };
  uint16_t dstTmpTransp[16 * 8] = { 0 };
  int strideTmp = 8;
  int strideTmpTransp = 16;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_16x8_sse4(above, dst, stride, dx, upsample_above);
  predict_8x16_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 16, 8);

  for (r = 0; r < 8; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base));
    ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base));
    a = _mm256_loadu_si256((__m256i *)(dst + r * stride));
    l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm256_and_si256(a, ma);
    leftReg = _mm256_and_si256(l, ml);
    result = _mm256_or_si256(aboveReg, leftReg);
    result = _mm256_min_epu16(clip_bd, result);
    _mm256_storeu_si256((__m256i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_8x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m128i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[32 * 8] = { 0 };
  uint16_t dstTmpTransp[8 * 32] = { 0 };
  int strideTmp = 32;
  int strideTmpTransp = 8;
  __m128i clip_bd =
      (bd == 8) ? _mm_set1_epi16(255)
                : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

  predict_8x32_sse4(above, dst, stride, dx, upsample_above);
  predict_32x8_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 8, 32);

  for (r = 0; r < 32; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - (base1 >> upsample_above);
    base = 0;
    if ((c_pos >= (0 - upsample_above)) && (c_pos < 64)) {
      base = 64 - c_pos;
    }

    ma = _mm_loadu_si128((__m128i *)(z2BlendMaskabove + base));
    ml = _mm_loadu_si128((__m128i *)(z2BlendMaskleft + base));
    a = _mm_loadu_si128((__m128i *)(dst + r * stride));
    l = _mm_loadu_si128((__m128i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm_and_si128(a, ma);
    leftReg = _mm_and_si128(l, ml);
    result = _mm_or_si128(aboveReg, leftReg);
    result = _mm_min_epu16(clip_bd, result);
    _mm_storeu_si128((__m128i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_32x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *above,
                                           const uint16_t *left,
                                           int upsample_above,
                                           int upsample_left, int dx, int dy,
                                           int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[8 * 32] = { 0 };
  uint16_t dstTmpTransp[32 * 8] = { 0 };
  int strideTmp = 8;
  int strideTmpTransp = 32;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_32x8_sse4(above, dst, stride, dx, upsample_above);
  predict_8x32_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 32, 8);

  for (r = 0; r < 8; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 32; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256(
          (__m256i *)(dstTmpTransp + r * strideTmpTransp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_16x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[16 * 16] = { 0 };
  uint16_t dstTmpTransp[16 * 16] = { 0 };
  int strideTmp = 16;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_16x16_sse4(above, dst, stride, dx, upsample_above);
  predict_16x16_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmp, 16, 16);

  for (r = 0; r < 16; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base));
    ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base));
    a = _mm256_loadu_si256((__m256i *)(dst + r * stride));
    l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmp));
    aboveReg = _mm256_and_si256(a, ma);
    leftReg = _mm256_and_si256(l, ml);
    result = _mm256_or_si256(aboveReg, leftReg);
    result = _mm256_min_epu16(clip_bd, result);
    _mm256_storeu_si256((__m256i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_32x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[32 * 32] = { 0 };
  uint16_t dstTmpTransp[32 * 32] = { 0 };
  int strideTmp = 32;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_32x32_sse4(above, dst, stride, dx, upsample_above);
  predict_32x32_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmp, 32, 32);

  for (r = 0; r < 32; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 32; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_64x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[64 * 64] = { 0 };
  uint16_t dstTmpTransp[64 * 64] = { 0 };
  int strideTmp = 64;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_64x64_sse4(above, dst, stride, dx, upsample_above);
  predict_64x64_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmp, 64, 64);

  for (r = 0; r < 64; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 64; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_16x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[32 * 16] = { 0 };
  uint16_t dstTmpTransp[16 * 32] = { 0 };
  int strideTmp = 32;
  int strideTmpTransp = 16;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_16x32_sse4(above, dst, stride, dx, upsample_above);
  predict_32x16_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 16, 32);

  for (r = 0; r < 32; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base));
    ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base));
    a = _mm256_loadu_si256((__m256i *)(dst + r * stride));
    l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm256_and_si256(a, ma);
    leftReg = _mm256_and_si256(l, ml);
    result = _mm256_or_si256(aboveReg, leftReg);
    result = _mm256_min_epu16(clip_bd, result);
    _mm256_storeu_si256((__m256i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_32x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[16 * 32] = { 0 };
  uint16_t dstTmpTransp[32 * 16] = { 0 };
  int strideTmp = 16;
  int strideTmpTransp = 32;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_32x16_sse4(above, dst, stride, dx, upsample_above);
  predict_16x32_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 32, 16);

  for (r = 0; r < 16; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 32; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256(
          (__m256i *)(dstTmpTransp + r * strideTmpTransp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_32x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[64 * 32] = { 0 };
  uint16_t dstTmpTransp[32 * 64] = { 0 };
  int strideTmp = 64;
  int strideTmpTransp = 32;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_32x64_sse4(above, dst, stride, dx, upsample_above);
  predict_64x32_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 32, 64);

  for (r = 0; r < 64; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 32; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256(
          (__m256i *)(dstTmpTransp + r * strideTmpTransp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_64x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[32 * 64] = { 0 };
  uint16_t dstTmpTransp[64 * 32] = { 0 };
  int strideTmp = 32;
  int strideTmpTransp = 64;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_64x32_sse4(above, dst, stride, dx, upsample_above);
  predict_32x64_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 64, 32);

  for (r = 0; r < 32; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 64; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256(
          (__m256i *)(dstTmpTransp + r * strideTmpTransp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
void av1_highbd_dr_prediction_z2_16x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[64 * 16] = { 0 };
  uint16_t dstTmpTransp[16 * 64] = { 0 };
  int strideTmp = 64;
  int strideTmpTransp = 16;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_16x64_sse4(above, dst, stride, dx, upsample_above);
  predict_64x16_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 16, 64);

  for (r = 0; r < 64; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base));
    ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base));
    a = _mm256_loadu_si256((__m256i *)(dst + r * stride));
    l = _mm256_loadu_si256((__m256i *)(dstTmpTransp + r * strideTmpTransp));
    aboveReg = _mm256_and_si256(a, ma);
    leftReg = _mm256_and_si256(l, ml);
    result = _mm256_or_si256(aboveReg, leftReg);
    result = _mm256_min_epu16(clip_bd, result);
    _mm256_storeu_si256((__m256i *)(dst + r * stride), result);
  }
}
void av1_highbd_dr_prediction_z2_64x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *above,
                                            const uint16_t *left,
                                            int upsample_above,
                                            int upsample_left, int dx, int dy,
                                            int bd) {
  int r, c, x, base, base1, c_pos;
  const int frac_bits_x = 6 - upsample_above;
  __m256i ma, ml, a, l, aboveReg, leftReg, result;
  uint16_t dstTmp[16 * 64] = { 0 };
  uint16_t dstTmpTransp[64 * 16] = { 0 };
  int strideTmp = 16;
  int strideTmpTransp = 64;
  __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255)
                              : (bd == 10) ? _mm256_set1_epi16(1023)
                                           : _mm256_set1_epi16(4095);

  predict_64x16_sse4(above, dst, stride, dx, upsample_above);
  predict_16x64_sse4(left, dstTmp, strideTmp, dy, upsample_left);
  transpose(dstTmp, strideTmp, dstTmpTransp, strideTmpTransp, 64, 16);

  for (r = 0; r < 16; r++) {
    x = -dx * (r + 1);
    base1 = x >> frac_bits_x;
    c_pos = -1 - base1;
    base = 0;
    if (c_pos >= 0 && c_pos < 64) {
      base = 64 - c_pos;
    }

    for (c = 0; c < 64; c += 16) {
      ma = _mm256_loadu_si256((__m256i *)(z2BlendMaskabove + base + c));
      ml = _mm256_loadu_si256((__m256i *)(z2BlendMaskleft + base + c));
      a = _mm256_loadu_si256((__m256i *)(dst + r * stride + c));
      l = _mm256_loadu_si256(
          (__m256i *)(dstTmpTransp + r * strideTmpTransp + c));
      aboveReg = _mm256_and_si256(a, ma);
      leftReg = _mm256_and_si256(l, ml);
      result = _mm256_or_si256(aboveReg, leftReg);
      result = _mm256_min_epu16(clip_bd, result);
      _mm256_storeu_si256((__m256i *)(dst + r * stride + c), result);
    }
  }
}
// Directional prediction, zone 2: 90 < angle < 180
void av1_highbd_dr_prediction_z2_avx2(uint16_t *dst, ptrdiff_t stride, int bw,
                                      int bh, const uint16_t *above,
                                      const uint16_t *left, int upsample_above,
                                      int upsample_left, int dx, int dy,
                                      int bd) {
  assert(dx > 0);
  assert(dy > 0);

  if (bw == bh) {
    switch (bw) {
      case 4:
        av1_highbd_dr_prediction_z2_4x4_avx2(dst, stride, above, left,
                                             upsample_above, upsample_left, dx,
                                             dy, bd);
        break;
      case 8:
        av1_highbd_dr_prediction_z2_8x8_avx2(dst, stride, above, left,
                                             upsample_above, upsample_left, dx,
                                             dy, bd);
        break;
      case 16:
        av1_highbd_dr_prediction_z2_16x16_avx2(dst, stride, above, left,
                                               upsample_above, upsample_left,
                                               dx, dy, bd);
        break;
      case 32:
        av1_highbd_dr_prediction_z2_32x32_avx2(dst, stride, above, left,
                                               upsample_above, upsample_left,
                                               dx, dy, bd);
        break;
      case 64:
        av1_highbd_dr_prediction_z2_64x64_avx2(dst, stride, above, left,
                                               upsample_above, upsample_left,
                                               dx, dy, bd);
        break;
    }
  } else {
    if (bw < bh) {
      if (bw + bw == bh) {
        switch (bw) {
          case 4:
            av1_highbd_dr_prediction_z2_4x8_avx2(dst, stride, above, left,
                                                 upsample_above, upsample_left,
                                                 dx, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z2_8x16_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z2_16x32_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
          case 32:
            av1_highbd_dr_prediction_z2_32x64_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
        }
      } else {
        switch (bw) {
          case 4:
            av1_highbd_dr_prediction_z2_4x16_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z2_8x32_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z2_16x64_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
        }
      }
    } else {
      if (bh + bh == bw) {
        switch (bh) {
          case 4:
            av1_highbd_dr_prediction_z2_8x4_avx2(dst, stride, above, left,
                                                 upsample_above, upsample_left,
                                                 dx, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z2_16x8_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z2_32x16_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
          case 32:
            av1_highbd_dr_prediction_z2_64x32_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
        }
      } else {
        switch (bh) {
          case 4:
            av1_highbd_dr_prediction_z2_16x4_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z2_32x8_avx2(dst, stride, above, left,
                                                  upsample_above, upsample_left,
                                                  dx, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z2_64x16_avx2(dst, stride, above, left,
                                                   upsample_above,
                                                   upsample_left, dx, dy, bd);
            break;
        }
      }
    }
  }
}


void av1_highbd_dr_prediction_z3_4x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *left,
                                          int upsample_left, int dy, int bd) {
  __m128i dstvec[4], d[4];
  (void)bd;
  av1_highbd_dr_prediction_z1_4xN_internal_avx2(4, dstvec, left, upsample_left,
                                                dy);
  highbd_transpose4x8_8x4_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2],
                                   &dstvec[3], &d[0], &d[1], &d[2], &d[3]);
  _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
  _mm_storel_epi64((__m128i *)(dst + 1 * stride), d[1]);
  _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[2]);
  _mm_storel_epi64((__m128i *)(dst + 3 * stride), d[3]);
  return;
}

void av1_highbd_dr_prediction_z3_8x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *left,
                                          int upsample_left, int dy, int bd) {
  __m128i dstvec[8], d[8];
  (void)bd;
  av1_highbd_dr_prediction_z1_8xN_internal_avx2(8, dstvec, left, upsample_left,
                                                dy);
  highbd_transpose8x8_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
                           &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7],
                           &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6],
                           &d[7]);
  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
  }
}

void av1_highbd_dr_prediction_z3_4x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *left,
                                          int upsample_left, int dy, int bd) {
  __m128i dstvec[4], d[8];
  (void)bd;
  av1_highbd_dr_prediction_z1_8xN_internal_avx2(4, dstvec, left, upsample_left,
                                                dy);
  highbd_transpose4x8_8x4_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
                               &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6],
                               &d[7]);
  for (int i = 0; i < 8; i++) {
    _mm_storel_epi64((__m128i *)(dst + i * stride), d[i]);
  }
}

void av1_highbd_dr_prediction_z3_8x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                          const uint16_t *left,
                                          int upsample_left, int dy, int bd) {
  __m128i dstvec[8], d[4];
  (void)bd;
  av1_highbd_dr_prediction_z1_4xN_internal_avx2(8, dstvec, left, upsample_left,
                                                dy);
  highbd_transpose8x8_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
                               &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7],
                               &d[0], &d[1], &d[2], &d[3]);
  _mm_storeu_si128((__m128i *)(dst + 0 * stride), d[0]);
  _mm_storeu_si128((__m128i *)(dst + 1 * stride), d[1]);
  _mm_storeu_si128((__m128i *)(dst + 2 * stride), d[2]);
  _mm_storeu_si128((__m128i *)(dst + 3 * stride), d[3]);
}

void av1_highbd_dr_prediction_z3_8x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m256i dstvec[8], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(8, dstvec, left, upsample_left,
                                                 dy);
  highbd_transpose8x16_16x8_avx2(dstvec, d);
  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride),
                     _mm256_castsi256_si128(d[i]));
  }
  for (int i = 8; i < 16; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride),
                     _mm256_extracti128_si256(d[i - 8], 1));
  }
}

void av1_highbd_dr_prediction_z3_16x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m128i dstvec[16], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_8xN_internal_avx2(16, dstvec, left, upsample_left,
                                                dy);
  for (int i = 0; i < 16; i += 8) {
    highbd_transpose8x8_sse2(&dstvec[0 + i], &dstvec[1 + i], &dstvec[2 + i],
                             &dstvec[3 + i], &dstvec[4 + i], &dstvec[5 + i],
                             &dstvec[6 + i], &dstvec[7 + i], &d[0 + i],
                             &d[1 + i], &d[2 + i], &d[3 + i], &d[4 + i],
                             &d[5 + i], &d[6 + i], &d[7 + i]);
  }
  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    _mm_storeu_si128((__m128i *)(dst + i * stride + 8), d[i + 8]);
  }
}

void av1_highbd_dr_prediction_z3_4x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m256i dstvec[4], d[4], d1;
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(4, dstvec, left, upsample_left,
                                                 dy);
  highbd_transpose4x16_avx2(dstvec, d);
  for (int i = 0; i < 4; i++) {
    _mm_storel_epi64((__m128i *)(dst + i * stride),
                     _mm256_castsi256_si128(d[i]));
    d1 = _mm256_bsrli_epi128(d[i], 8);
    _mm_storel_epi64((__m128i *)(dst + (i + 4) * stride),
                     _mm256_castsi256_si128(d1));
    _mm_storel_epi64((__m128i *)(dst + (i + 8) * stride),
                     _mm256_extracti128_si256(d[i], 1));
    _mm_storel_epi64((__m128i *)(dst + (i + 12) * stride),
                     _mm256_extracti128_si256(d1, 1));
  }
}

void av1_highbd_dr_prediction_z3_16x4_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m128i dstvec[16], d[8];
  (void)bd;
  av1_highbd_dr_prediction_z1_4xN_internal_avx2(16, dstvec, left, upsample_left,
                                                dy);
  highbd_transpose16x4_8x8_sse2(dstvec, d);

  _mm_storeu_si128((__m128i *)(dst + 0 * stride), d[0]);
  _mm_storeu_si128((__m128i *)(dst + 0 * stride + 8), d[1]);
  _mm_storeu_si128((__m128i *)(dst + 1 * stride), d[2]);
  _mm_storeu_si128((__m128i *)(dst + 1 * stride + 8), d[3]);
  _mm_storeu_si128((__m128i *)(dst + 2 * stride), d[4]);
  _mm_storeu_si128((__m128i *)(dst + 2 * stride + 8), d[5]);
  _mm_storeu_si128((__m128i *)(dst + 3 * stride), d[6]);
  _mm_storeu_si128((__m128i *)(dst + 3 * stride + 8), d[7]);
}

void av1_highbd_dr_prediction_z3_8x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m256i dstvec[16], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_32xN_internal_avx2(8, dstvec, left, upsample_left,
                                                 dy);
  for (int i = 0; i < 16; i += 8) {
    highbd_transpose8x16_16x8_avx2(dstvec + i, d + i);
  }

  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride),
                     _mm256_castsi256_si128(d[i]));
  }
  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + (i + 8) * stride),
                     _mm256_extracti128_si256(d[i], 1));
  }
  for (int i = 8; i < 16; i++) {
    _mm_storeu_si128((__m128i *)(dst + (i + 8) * stride),
                     _mm256_castsi256_si128(d[i]));
  }
  for (int i = 8; i < 16; i++) {
    _mm_storeu_si128((__m128i *)(dst + (i + 16) * stride),
                     _mm256_extracti128_si256(d[i], 1));
  }
}

void av1_highbd_dr_prediction_z3_32x8_avx2(uint16_t *dst, ptrdiff_t stride,
                                           const uint16_t *left,
                                           int upsample_left, int dy, int bd) {
  __m128i dstvec[32], d[32];
  (void)bd;
  av1_highbd_dr_prediction_z1_8xN_internal_avx2(32, dstvec, left, upsample_left,
                                                dy);
  for (int i = 0; i < 32; i += 8) {
    highbd_transpose8x8_sse2(&dstvec[0 + i], &dstvec[1 + i], &dstvec[2 + i],
                             &dstvec[3 + i], &dstvec[4 + i], &dstvec[5 + i],
                             &dstvec[6 + i], &dstvec[7 + i], &d[0 + i],
                             &d[1 + i], &d[2 + i], &d[3 + i], &d[4 + i],
                             &d[5 + i], &d[6 + i], &d[7 + i]);
  }
  for (int i = 0; i < 8; i++) {
    _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    _mm_storeu_si128((__m128i *)(dst + i * stride + 8), d[i + 8]);
    _mm_storeu_si128((__m128i *)(dst + i * stride + 16), d[i + 16]);
    _mm_storeu_si128((__m128i *)(dst + i * stride + 24), d[i + 24]);
  }
}

void av1_highbd_dr_prediction_z3_16x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  __m256i dstvec[16], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(16, dstvec, left,
                                                 upsample_left, dy);
  highbd_transpose16x16_avx2(dstvec, d);

  for (int i = 0; i < 16; i++) {
    _mm256_storeu_si256((__m256i *)(dst + i * stride), d[i]);
  }
}

void av1_highbd_dr_prediction_z3_32x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  __m256i dstvec[64], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_32xN_internal_avx2(32, dstvec, left,
                                                 upsample_left, dy);

  highbd_transpose16x16_avx2(dstvec, d);
  for (int j = 0; j < 16; j++) {
    _mm256_storeu_si256((__m256i *)(dst + j * stride), d[j]);
  }
  highbd_transpose16x16_avx2(dstvec + 16, d);
  for (int j = 0; j < 16; j++) {
    _mm256_storeu_si256((__m256i *)(dst + j * stride + 16), d[j]);
  }
  highbd_transpose16x16_avx2(dstvec + 32, d);
  for (int j = 0; j < 16; j++) {
    _mm256_storeu_si256((__m256i *)(dst + (j + 16) * stride), d[j]);
  }
  highbd_transpose16x16_avx2(dstvec + 48, d);
  for (int j = 0; j < 16; j++) {
    _mm256_storeu_si256((__m256i *)(dst + (j + 16) * stride + 16), d[j]);
  }
}

void av1_highbd_dr_prediction_z3_64x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  uint16_t dstT[64 * 64];
  av1_highbd_dr_prediction_z1_64xN_avx2(64, dstT, 64, left, upsample_left, dy,
                                        bd);
  transpose(dstT, 64, dst, stride, 64, 64);
}

void av1_highbd_dr_prediction_z3_16x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  __m256i dstvec[32], d[32];
  (void)bd;
  av1_highbd_dr_prediction_z1_32xN_internal_avx2(16, dstvec, left,
                                                 upsample_left, dy);
  for (int i = 0; i < 32; i += 8) {
    highbd_transpose8x16_16x8_avx2(dstvec + i, d + i);
  }
  // store
  for (int j = 0; j < 32; j += 16) {
    for (int i = 0; i < 8; i++) {
      _mm_storeu_si128((__m128i *)(dst + (i + j) * stride),
                       _mm256_castsi256_si128(d[(i + j)]));
    }
    for (int i = 0; i < 8; i++) {
      _mm_storeu_si128((__m128i *)(dst + (i + j) * stride + 8),
                       _mm256_castsi256_si128(d[(i + j) + 8]));
    }
    for (int i = 8; i < 16; i++) {
      _mm256_storeu_si256(
          (__m256i *)(dst + (i + j) * stride),
          _mm256_inserti128_si256(
              d[(i + j)], _mm256_extracti128_si256(d[(i + j) - 8], 1), 0));
    }
  }
}

void av1_highbd_dr_prediction_z3_32x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  __m256i dstvec[32], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(32, dstvec, left,
                                                 upsample_left, dy);
  for (int i = 0; i < 32; i += 16) {
    highbd_transpose16x16_avx2((dstvec + i), d);
    for (int j = 0; j < 16; j++) {
      _mm256_storeu_si256((__m256i *)(dst + j * stride + i), d[j]);
    }
  }
}

void av1_highbd_dr_prediction_z3_32x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  uint16_t dstT[64 * 32];
  av1_highbd_dr_prediction_z1_64xN_avx2(32, dstT, 64, left, upsample_left, dy,
                                        bd);
  transpose(dstT, 64, dst, stride, 32, 64);
}

void av1_highbd_dr_prediction_z3_64x32_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  uint16_t dstT[32 * 64];
  av1_highbd_dr_prediction_z1_32xN_avx2(64, dstT, 32, left, upsample_left, dy,
                                        bd);
  transpose(dstT, 32, dst, stride, 64, 32);
  return;
}

void av1_highbd_dr_prediction_z3_16x64_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  uint16_t dstT[64 * 16];
  av1_highbd_dr_prediction_z1_64xN_avx2(16, dstT, 64, left, upsample_left, dy,
                                        bd);
  transpose(dstT, 64, dst, stride, 16, 64);
}

void av1_highbd_dr_prediction_z3_64x16_avx2(uint16_t *dst, ptrdiff_t stride,
                                            const uint16_t *left,
                                            int upsample_left, int dy, int bd) {
  __m256i dstvec[64], d[16];
  (void)bd;
  av1_highbd_dr_prediction_z1_16xN_internal_avx2(64, dstvec, left,
                                                 upsample_left, dy);
  for (int i = 0; i < 64; i += 16) {
    highbd_transpose16x16_avx2((dstvec + i), d);
    for (int j = 0; j < 16; j++) {
      _mm256_storeu_si256((__m256i *)(dst + j * stride + i), d[j]);
    }
  }
}

void av1_highbd_dr_prediction_z3_avx2(uint16_t *dst, ptrdiff_t stride, int bw,
                                      int bh, const uint16_t *above,
                                      const uint16_t *left, int upsample_left,
                                      int dx, int dy, int bd) {
  (void)above;
  (void)dx;
  assert(dx == 1);
  assert(dy > 0);

  if (bw == bh) {
    switch (bw) {
      case 4:
        av1_highbd_dr_prediction_z3_4x4_avx2(dst, stride, left, upsample_left,
                                             dy, bd);
        break;
      case 8:
        av1_highbd_dr_prediction_z3_8x8_avx2(dst, stride, left, upsample_left,
                                             dy, bd);
        break;
      case 16:
        av1_highbd_dr_prediction_z3_16x16_avx2(dst, stride, left, upsample_left,
                                               dy, bd);
        break;
      case 32:
        av1_highbd_dr_prediction_z3_32x32_avx2(dst, stride, left, upsample_left,
                                               dy, bd);
        break;
      case 64:
        av1_highbd_dr_prediction_z3_64x64_avx2(dst, stride, left, upsample_left,
                                               dy, bd);
        break;
    }
  } else {
    if (bw < bh) {
      if (bw + bw == bh) {
        switch (bw) {
          case 4:
            av1_highbd_dr_prediction_z3_4x8_avx2(dst, stride, left,
                                                 upsample_left, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z3_8x16_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z3_16x32_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
          case 32:
            av1_highbd_dr_prediction_z3_32x64_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
        }
      } else {
        switch (bw) {
          case 4:
            av1_highbd_dr_prediction_z3_4x16_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z3_8x32_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z3_16x64_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
        }
      }
    } else {
      if (bh + bh == bw) {
        switch (bh) {
          case 4:
            av1_highbd_dr_prediction_z3_8x4_avx2(dst, stride, left,
                                                 upsample_left, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z3_16x8_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z3_32x16_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
          case 32:
            av1_highbd_dr_prediction_z3_64x32_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
        }
      } else {
        switch (bh) {
          case 4:
            av1_highbd_dr_prediction_z3_16x4_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 8:
            av1_highbd_dr_prediction_z3_32x8_avx2(dst, stride, left,
                                                  upsample_left, dy, bd);
            break;
          case 16:
            av1_highbd_dr_prediction_z3_64x16_avx2(dst, stride, left,
                                                   upsample_left, dy, bd);
            break;
        }
      }
    }
  }
  return;
}
