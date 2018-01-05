/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <tmmintrin.h>

#include "./aom_config.h"
#include "./av1_rtcd.h"
#include "av1/common/filter.h"

#define WIDTH_BOUND (16)
#define HEIGHT_BOUND (16)

static INLINE void transpose_4x8(const __m128i *in, __m128i *out) {
  __m128i t0, t1;

  t0 = _mm_unpacklo_epi16(in[0], in[1]);
  t1 = _mm_unpacklo_epi16(in[2], in[3]);

  out[0] = _mm_unpacklo_epi32(t0, t1);
  out[1] = _mm_srli_si128(out[0], 8);
  out[2] = _mm_unpackhi_epi32(t0, t1);
  out[3] = _mm_srli_si128(out[2], 8);

  t0 = _mm_unpackhi_epi16(in[0], in[1]);
  t1 = _mm_unpackhi_epi16(in[2], in[3]);

  out[4] = _mm_unpacklo_epi32(t0, t1);
  out[5] = _mm_srli_si128(out[4], 8);
  // Note: We ignore out[6] and out[7] because
  // they're zero vectors.
}

typedef void (*store_pixel_t)(const __m128i *x, uint8_t *dst);

static INLINE __m128i accumulate_store(const __m128i *x, uint8_t *src) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);
  __m128i y = _mm_loadl_epi64((__m128i const *)src);
  y = _mm_unpacklo_epi8(y, zero);
  y = _mm_add_epi16(*x, y);
  y = _mm_add_epi16(y, one);
  y = _mm_srai_epi16(y, 1);
  y = _mm_packus_epi16(y, y);
  return y;
}

static INLINE void store_2_pixel_only(const __m128i *x, uint8_t *dst) {
  uint32_t temp;
  __m128i u = _mm_packus_epi16(*x, *x);
  temp = _mm_cvtsi128_si32(u);
  *(uint16_t *)dst = (uint16_t)temp;
}

static INLINE void accumulate_store_2_pixel(const __m128i *x, uint8_t *dst) {
  uint32_t temp;
  __m128i y = accumulate_store(x, dst);
  temp = _mm_cvtsi128_si32(y);
  *(uint16_t *)dst = (uint16_t)temp;
}

static INLINE void store_4_pixel_only(const __m128i *x, uint8_t *dst) {
  __m128i u = _mm_packus_epi16(*x, *x);
  *(int *)dst = _mm_cvtsi128_si32(u);
}

static INLINE void accumulate_store_4_pixel(const __m128i *x, uint8_t *dst) {
  __m128i y = accumulate_store(x, dst);
  *(int *)dst = _mm_cvtsi128_si32(y);
}

// Vertical 8-pixel parallel
typedef void (*transpose_to_dst_t)(const uint16_t *src, int src_stride,
                                   uint8_t *dst, int dst_stride);

static INLINE void transpose8x8_direct_to_dst(const uint16_t *src,
                                              int src_stride, uint8_t *dst,
                                              int dst_stride) {
  const __m128i k_256 = _mm_set1_epi16(1 << 8);
  __m128i v0, v1, v2, v3;

  __m128i u0 = _mm_loadu_si128((__m128i const *)(src + 0 * src_stride));
  __m128i u1 = _mm_loadu_si128((__m128i const *)(src + 1 * src_stride));
  __m128i u2 = _mm_loadu_si128((__m128i const *)(src + 2 * src_stride));
  __m128i u3 = _mm_loadu_si128((__m128i const *)(src + 3 * src_stride));
  __m128i u4 = _mm_loadu_si128((__m128i const *)(src + 4 * src_stride));
  __m128i u5 = _mm_loadu_si128((__m128i const *)(src + 5 * src_stride));
  __m128i u6 = _mm_loadu_si128((__m128i const *)(src + 6 * src_stride));
  __m128i u7 = _mm_loadu_si128((__m128i const *)(src + 7 * src_stride));

  u0 = _mm_mulhrs_epi16(u0, k_256);
  u1 = _mm_mulhrs_epi16(u1, k_256);
  u2 = _mm_mulhrs_epi16(u2, k_256);
  u3 = _mm_mulhrs_epi16(u3, k_256);
  u4 = _mm_mulhrs_epi16(u4, k_256);
  u5 = _mm_mulhrs_epi16(u5, k_256);
  u6 = _mm_mulhrs_epi16(u6, k_256);
  u7 = _mm_mulhrs_epi16(u7, k_256);

  v0 = _mm_packus_epi16(u0, u1);
  v1 = _mm_packus_epi16(u2, u3);
  v2 = _mm_packus_epi16(u4, u5);
  v3 = _mm_packus_epi16(u6, u7);

  u0 = _mm_unpacklo_epi8(v0, v1);
  u1 = _mm_unpackhi_epi8(v0, v1);
  u2 = _mm_unpacklo_epi8(v2, v3);
  u3 = _mm_unpackhi_epi8(v2, v3);

  u4 = _mm_unpacklo_epi8(u0, u1);
  u5 = _mm_unpacklo_epi8(u2, u3);
  u6 = _mm_unpackhi_epi8(u0, u1);
  u7 = _mm_unpackhi_epi8(u2, u3);

  u0 = _mm_unpacklo_epi32(u4, u5);
  u1 = _mm_unpackhi_epi32(u4, u5);
  u2 = _mm_unpacklo_epi32(u6, u7);
  u3 = _mm_unpackhi_epi32(u6, u7);

  u4 = _mm_srli_si128(u0, 8);
  u5 = _mm_srli_si128(u1, 8);
  u6 = _mm_srli_si128(u2, 8);
  u7 = _mm_srli_si128(u3, 8);

  _mm_storel_epi64((__m128i *)dst, u0);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 1), u4);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 2), u1);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 3), u5);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 4), u2);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 5), u6);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 6), u3);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 7), u7);
}

static INLINE void transpose8x8_accumu_to_dst(const uint16_t *src,
                                              int src_stride, uint8_t *dst,
                                              int dst_stride) {
  const __m128i k_256 = _mm_set1_epi16(1 << 8);
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);
  __m128i v0, v1, v2, v3, v4, v5, v6, v7;

  __m128i u0 = _mm_loadu_si128((__m128i const *)(src + 0 * src_stride));
  __m128i u1 = _mm_loadu_si128((__m128i const *)(src + 1 * src_stride));
  __m128i u2 = _mm_loadu_si128((__m128i const *)(src + 2 * src_stride));
  __m128i u3 = _mm_loadu_si128((__m128i const *)(src + 3 * src_stride));
  __m128i u4 = _mm_loadu_si128((__m128i const *)(src + 4 * src_stride));
  __m128i u5 = _mm_loadu_si128((__m128i const *)(src + 5 * src_stride));
  __m128i u6 = _mm_loadu_si128((__m128i const *)(src + 6 * src_stride));
  __m128i u7 = _mm_loadu_si128((__m128i const *)(src + 7 * src_stride));

  u0 = _mm_mulhrs_epi16(u0, k_256);
  u1 = _mm_mulhrs_epi16(u1, k_256);
  u2 = _mm_mulhrs_epi16(u2, k_256);
  u3 = _mm_mulhrs_epi16(u3, k_256);
  u4 = _mm_mulhrs_epi16(u4, k_256);
  u5 = _mm_mulhrs_epi16(u5, k_256);
  u6 = _mm_mulhrs_epi16(u6, k_256);
  u7 = _mm_mulhrs_epi16(u7, k_256);

  v0 = _mm_packus_epi16(u0, u1);
  v1 = _mm_packus_epi16(u2, u3);
  v2 = _mm_packus_epi16(u4, u5);
  v3 = _mm_packus_epi16(u6, u7);

  u0 = _mm_unpacklo_epi8(v0, v1);
  u1 = _mm_unpackhi_epi8(v0, v1);
  u2 = _mm_unpacklo_epi8(v2, v3);
  u3 = _mm_unpackhi_epi8(v2, v3);

  u4 = _mm_unpacklo_epi8(u0, u1);
  u5 = _mm_unpacklo_epi8(u2, u3);
  u6 = _mm_unpackhi_epi8(u0, u1);
  u7 = _mm_unpackhi_epi8(u2, u3);

  u0 = _mm_unpacklo_epi32(u4, u5);
  u1 = _mm_unpackhi_epi32(u4, u5);
  u2 = _mm_unpacklo_epi32(u6, u7);
  u3 = _mm_unpackhi_epi32(u6, u7);

  u4 = _mm_srli_si128(u0, 8);
  u5 = _mm_srli_si128(u1, 8);
  u6 = _mm_srli_si128(u2, 8);
  u7 = _mm_srli_si128(u3, 8);

  v0 = _mm_loadl_epi64((__m128i const *)(dst + 0 * dst_stride));
  v1 = _mm_loadl_epi64((__m128i const *)(dst + 1 * dst_stride));
  v2 = _mm_loadl_epi64((__m128i const *)(dst + 2 * dst_stride));
  v3 = _mm_loadl_epi64((__m128i const *)(dst + 3 * dst_stride));
  v4 = _mm_loadl_epi64((__m128i const *)(dst + 4 * dst_stride));
  v5 = _mm_loadl_epi64((__m128i const *)(dst + 5 * dst_stride));
  v6 = _mm_loadl_epi64((__m128i const *)(dst + 6 * dst_stride));
  v7 = _mm_loadl_epi64((__m128i const *)(dst + 7 * dst_stride));

  u0 = _mm_unpacklo_epi8(u0, zero);
  u1 = _mm_unpacklo_epi8(u1, zero);
  u2 = _mm_unpacklo_epi8(u2, zero);
  u3 = _mm_unpacklo_epi8(u3, zero);
  u4 = _mm_unpacklo_epi8(u4, zero);
  u5 = _mm_unpacklo_epi8(u5, zero);
  u6 = _mm_unpacklo_epi8(u6, zero);
  u7 = _mm_unpacklo_epi8(u7, zero);

  v0 = _mm_unpacklo_epi8(v0, zero);
  v1 = _mm_unpacklo_epi8(v1, zero);
  v2 = _mm_unpacklo_epi8(v2, zero);
  v3 = _mm_unpacklo_epi8(v3, zero);
  v4 = _mm_unpacklo_epi8(v4, zero);
  v5 = _mm_unpacklo_epi8(v5, zero);
  v6 = _mm_unpacklo_epi8(v6, zero);
  v7 = _mm_unpacklo_epi8(v7, zero);

  v0 = _mm_adds_epi16(u0, v0);
  v1 = _mm_adds_epi16(u4, v1);
  v2 = _mm_adds_epi16(u1, v2);
  v3 = _mm_adds_epi16(u5, v3);
  v4 = _mm_adds_epi16(u2, v4);
  v5 = _mm_adds_epi16(u6, v5);
  v6 = _mm_adds_epi16(u3, v6);
  v7 = _mm_adds_epi16(u7, v7);

  v0 = _mm_adds_epi16(v0, one);
  v1 = _mm_adds_epi16(v1, one);
  v2 = _mm_adds_epi16(v2, one);
  v3 = _mm_adds_epi16(v3, one);
  v4 = _mm_adds_epi16(v4, one);
  v5 = _mm_adds_epi16(v5, one);
  v6 = _mm_adds_epi16(v6, one);
  v7 = _mm_adds_epi16(v7, one);

  v0 = _mm_srai_epi16(v0, 1);
  v1 = _mm_srai_epi16(v1, 1);
  v2 = _mm_srai_epi16(v2, 1);
  v3 = _mm_srai_epi16(v3, 1);
  v4 = _mm_srai_epi16(v4, 1);
  v5 = _mm_srai_epi16(v5, 1);
  v6 = _mm_srai_epi16(v6, 1);
  v7 = _mm_srai_epi16(v7, 1);

  u0 = _mm_packus_epi16(v0, v1);
  u1 = _mm_packus_epi16(v2, v3);
  u2 = _mm_packus_epi16(v4, v5);
  u3 = _mm_packus_epi16(v6, v7);

  u4 = _mm_srli_si128(u0, 8);
  u5 = _mm_srli_si128(u1, 8);
  u6 = _mm_srli_si128(u2, 8);
  u7 = _mm_srli_si128(u3, 8);

  _mm_storel_epi64((__m128i *)dst, u0);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 1), u4);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 2), u1);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 3), u5);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 4), u2);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 5), u6);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 6), u3);
  _mm_storel_epi64((__m128i *)(dst + dst_stride * 7), u7);
}

static INLINE void transpose_8x16(const __m128i *in, __m128i *out) {
  __m128i t0, t1, t2, t3, u0, u1;

  t0 = _mm_unpacklo_epi16(in[0], in[1]);
  t1 = _mm_unpacklo_epi16(in[2], in[3]);
  t2 = _mm_unpacklo_epi16(in[4], in[5]);
  t3 = _mm_unpacklo_epi16(in[6], in[7]);

  u0 = _mm_unpacklo_epi32(t0, t1);
  u1 = _mm_unpacklo_epi32(t2, t3);

  out[0] = _mm_unpacklo_epi64(u0, u1);
  out[1] = _mm_unpackhi_epi64(u0, u1);

  u0 = _mm_unpackhi_epi32(t0, t1);
  u1 = _mm_unpackhi_epi32(t2, t3);

  out[2] = _mm_unpacklo_epi64(u0, u1);
  out[3] = _mm_unpackhi_epi64(u0, u1);

  t0 = _mm_unpackhi_epi16(in[0], in[1]);
  t1 = _mm_unpackhi_epi16(in[2], in[3]);
  t2 = _mm_unpackhi_epi16(in[4], in[5]);
  t3 = _mm_unpackhi_epi16(in[6], in[7]);

  u0 = _mm_unpacklo_epi32(t0, t1);
  u1 = _mm_unpacklo_epi32(t2, t3);

  out[4] = _mm_unpacklo_epi64(u0, u1);
  out[5] = _mm_unpackhi_epi64(u0, u1);

  // Ignore out[6] and out[7]
  // they're zero vectors.
}

// Vertical 4-pixel parallel
static INLINE void transpose4x4_direct_to_dst(const uint16_t *src,
                                              int src_stride, uint8_t *dst,
                                              int dst_stride) {
  const __m128i k_256 = _mm_set1_epi16(1 << 8);
  __m128i v0, v1, v2, v3;

  // TODO(luoyi): two loads, 8 elements per load (two bytes per element)
  __m128i u0 = _mm_loadl_epi64((__m128i const *)(src + 0 * src_stride));
  __m128i u1 = _mm_loadl_epi64((__m128i const *)(src + 1 * src_stride));
  __m128i u2 = _mm_loadl_epi64((__m128i const *)(src + 2 * src_stride));
  __m128i u3 = _mm_loadl_epi64((__m128i const *)(src + 3 * src_stride));

  v0 = _mm_unpacklo_epi16(u0, u1);
  v1 = _mm_unpacklo_epi16(u2, u3);

  v2 = _mm_unpacklo_epi32(v0, v1);
  v3 = _mm_unpackhi_epi32(v0, v1);

  u0 = _mm_mulhrs_epi16(v2, k_256);
  u1 = _mm_mulhrs_epi16(v3, k_256);

  u0 = _mm_packus_epi16(u0, u1);
  u1 = _mm_srli_si128(u0, 4);
  u2 = _mm_srli_si128(u0, 8);
  u3 = _mm_srli_si128(u0, 12);

  *(int *)(dst) = _mm_cvtsi128_si32(u0);
  *(int *)(dst + dst_stride) = _mm_cvtsi128_si32(u1);
  *(int *)(dst + dst_stride * 2) = _mm_cvtsi128_si32(u2);
  *(int *)(dst + dst_stride * 3) = _mm_cvtsi128_si32(u3);
}

static INLINE void transpose4x4_accumu_to_dst(const uint16_t *src,
                                              int src_stride, uint8_t *dst,
                                              int dst_stride) {
  const __m128i k_256 = _mm_set1_epi16(1 << 8);
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);

  __m128i v0, v1, v2, v3;

  __m128i u0 = _mm_loadl_epi64((__m128i const *)(src));
  __m128i u1 = _mm_loadl_epi64((__m128i const *)(src + src_stride));
  __m128i u2 = _mm_loadl_epi64((__m128i const *)(src + 2 * src_stride));
  __m128i u3 = _mm_loadl_epi64((__m128i const *)(src + 3 * src_stride));

  v0 = _mm_unpacklo_epi16(u0, u1);
  v1 = _mm_unpacklo_epi16(u2, u3);

  v2 = _mm_unpacklo_epi32(v0, v1);
  v3 = _mm_unpackhi_epi32(v0, v1);

  u0 = _mm_mulhrs_epi16(v2, k_256);
  u1 = _mm_mulhrs_epi16(v3, k_256);

  u2 = _mm_packus_epi16(u0, u1);
  u0 = _mm_unpacklo_epi8(u2, zero);
  u1 = _mm_unpackhi_epi8(u2, zero);

  // load pixel values
  v0 = _mm_loadl_epi64((__m128i const *)(dst));
  v1 = _mm_loadl_epi64((__m128i const *)(dst + dst_stride));
  v2 = _mm_loadl_epi64((__m128i const *)(dst + 2 * dst_stride));
  v3 = _mm_loadl_epi64((__m128i const *)(dst + 3 * dst_stride));

  v0 = _mm_unpacklo_epi8(v0, zero);
  v1 = _mm_unpacklo_epi8(v1, zero);
  v2 = _mm_unpacklo_epi8(v2, zero);
  v3 = _mm_unpacklo_epi8(v3, zero);

  v0 = _mm_unpacklo_epi64(v0, v1);
  v1 = _mm_unpacklo_epi64(v2, v3);

  u0 = _mm_adds_epi16(u0, v0);
  u1 = _mm_adds_epi16(u1, v1);

  u0 = _mm_adds_epi16(u0, one);
  u1 = _mm_adds_epi16(u1, one);

  u0 = _mm_srai_epi16(u0, 1);
  u1 = _mm_srai_epi16(u1, 1);

  // saturation and pack to pixels
  u0 = _mm_packus_epi16(u0, u1);
  u1 = _mm_srli_si128(u0, 4);
  u2 = _mm_srli_si128(u0, 8);
  u3 = _mm_srli_si128(u0, 12);

  *(int *)(dst) = _mm_cvtsi128_si32(u0);
  *(int *)(dst + dst_stride) = _mm_cvtsi128_si32(u1);
  *(int *)(dst + dst_stride * 2) = _mm_cvtsi128_si32(u2);
  *(int *)(dst + dst_stride * 3) = _mm_cvtsi128_si32(u3);
}

// Note:
//  This function assumes:
// (1) 10/12-taps filters
// (2) x_step_q4 = 16 then filter is fixed at the call

void av1_convolve_horiz_ssse3(const uint8_t *src, int src_stride, uint8_t *dst,
                              int dst_stride, int w, int h,
                              const InterpFilterParams filter_params,
                              const int subpel_x_q4, int x_step_q4,
                              ConvolveParams *conv_params) {
  assert(conv_params->do_average == 0 || conv_params->do_average == 1);
  (void)x_step_q4;

  if (0 == subpel_x_q4 || 16 != x_step_q4) {
    av1_convolve_horiz_c(src, src_stride, dst, dst_stride, w, h, filter_params,
                         subpel_x_q4, x_step_q4, conv_params);
    return;
  }

  av1_convolve_horiz_c(src, src_stride, dst, dst_stride, w, h, filter_params,
                       subpel_x_q4, x_step_q4, conv_params);
}

// Vertical convolution filtering
static INLINE void store_8_pixel_only(const __m128i *x, uint8_t *dst) {
  __m128i u = _mm_packus_epi16(*x, *x);
  _mm_storel_epi64((__m128i *)dst, u);
}

static INLINE void accumulate_store_8_pixel(const __m128i *x, uint8_t *dst) {
  __m128i y = accumulate_store(x, dst);
  _mm_storel_epi64((__m128i *)dst, y);
}

void av1_convolve_vert_ssse3(const uint8_t *src, int src_stride, uint8_t *dst,
                             int dst_stride, int w, int h,
                             const InterpFilterParams filter_params,
                             const int subpel_y_q4, int y_step_q4,
                             ConvolveParams *conv_params) {
  assert(conv_params->do_average == 0 || conv_params->do_average == 1);

  if (0 == subpel_y_q4 || 16 != y_step_q4) {
    av1_convolve_vert_c(src, src_stride, dst, dst_stride, w, h, filter_params,
                        subpel_y_q4, y_step_q4, conv_params);
    return;
  }

  av1_convolve_vert_c(src, src_stride, dst, dst_stride, w, h, filter_params,
                      subpel_y_q4, y_step_q4, conv_params);
}

typedef struct SimdFilter {
  InterpFilter interp_filter;
  int8_t (*simd_horiz_filter)[2][16];
  int8_t (*simd_vert_filter)[6][16];
} SimdFilter;
