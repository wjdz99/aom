/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <emmintrin.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "aom_dsp/aom_dsp_common.h"

static INLINE void highbd_sad8x4_core_sse2(__m128i *s, __m128i *r,
                                           __m128i *sad_acc) {
  const __m128i zero = _mm_setzero_si128();
  int i;
  for (i = 0; i < 4; i++) {
    s[i] = _mm_or_si128(_mm_subs_epu16(s[i], r[i]), _mm_subs_epu16(r[i], s[i]));
  }

  s[0] = _mm_add_epi16(s[0], s[1]);
  s[0] = _mm_add_epi16(s[0], s[2]);
  s[0] = _mm_add_epi16(s[0], s[3]);

  r[0] = _mm_unpacklo_epi16(s[0], zero);
  r[1] = _mm_unpackhi_epi16(s[0], zero);

  r[0] = _mm_add_epi32(r[0], r[1]);
  *sad_acc = _mm_add_epi32(*sad_acc, r[0]);
}

static INLINE void sad4x4(const uint16_t *src_ptr, int src_stride,
                          const uint16_t *ref_ptr, int ref_stride,
                          __m128i *sad_acc) {
  __m128i s[4], r[4];
  s[0] = _mm_loadl_epi64((const __m128i *)src_ptr);
  s[1] = _mm_loadl_epi64((const __m128i *)(src_ptr + src_stride));
  s[2] = _mm_loadl_epi64((const __m128i *)(src_ptr + 2 * src_stride));
  s[3] = _mm_loadl_epi64((const __m128i *)(src_ptr + 3 * src_stride));

  r[0] = _mm_loadl_epi64((const __m128i *)ref_ptr);
  r[1] = _mm_loadl_epi64((const __m128i *)(ref_ptr + ref_stride));
  r[2] = _mm_loadl_epi64((const __m128i *)(ref_ptr + 2 * ref_stride));
  r[3] = _mm_loadl_epi64((const __m128i *)(ref_ptr + 3 * ref_stride));

  highbd_sad8x4_core_sse2(s, r, sad_acc);
}

static INLINE void sad8x4(const uint16_t *src_ptr, int src_stride,
                          const uint16_t *ref_ptr, int ref_stride,
                          __m128i *sad_acc) {
  __m128i s[4], r[4];
  s[0] = _mm_loadu_si128((const __m128i *)src_ptr);
  s[1] = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride));
  s[2] = _mm_loadu_si128((const __m128i *)(src_ptr + 2 * src_stride));
  s[3] = _mm_loadu_si128((const __m128i *)(src_ptr + 3 * src_stride));

  r[0] = _mm_loadu_si128((const __m128i *)ref_ptr);
  r[1] = _mm_loadu_si128((const __m128i *)(ref_ptr + ref_stride));
  r[2] = _mm_loadu_si128((const __m128i *)(ref_ptr + 2 * ref_stride));
  r[3] = _mm_loadu_si128((const __m128i *)(ref_ptr + 3 * ref_stride));

  highbd_sad8x4_core_sse2(s, r, sad_acc);
}

static void sad16x2(const uint16_t *src_ptr, int src_stride,
                    const uint16_t *ref_ptr, int ref_stride, __m128i *sad_acc) {
  __m128i s[4], r[4];
  int row_sections = 0;

  while (row_sections < 2) {
    s[0] = _mm_loadu_si128((const __m128i *)src_ptr);
    s[1] = _mm_loadu_si128((const __m128i *)(src_ptr + 8));
    s[2] = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride));
    s[3] = _mm_loadu_si128((const __m128i *)(src_ptr + src_stride + 8));

    r[0] = _mm_loadu_si128((const __m128i *)ref_ptr);
    r[1] = _mm_loadu_si128((const __m128i *)(ref_ptr + 8));
    r[2] = _mm_loadu_si128((const __m128i *)(ref_ptr + ref_stride));
    r[3] = _mm_loadu_si128((const __m128i *)(ref_ptr + ref_stride + 8));

    highbd_sad8x4_core_sse2(s, r, sad_acc);
    row_sections += 1;
    src_ptr += src_stride << 1;
    ref_ptr += ref_stride << 1;
  }
}

static void sad32x2(const uint16_t *src_ptr, int src_stride,
                    const uint16_t *ref_ptr, int ref_stride, __m128i *sad_acc) {
  __m128i s[4], r[4];
  int i;
  for (i = 0; i < 2; i++) {
    s[0] = _mm_loadu_si128((const __m128i *)src_ptr);
    s[1] = _mm_loadu_si128((const __m128i *)(src_ptr + 8));
    s[2] = _mm_loadu_si128((const __m128i *)(src_ptr + 16));
    s[3] = _mm_loadu_si128((const __m128i *)(src_ptr + 24));

    r[0] = _mm_loadu_si128((const __m128i *)ref_ptr);
    r[1] = _mm_loadu_si128((const __m128i *)(ref_ptr + 8));
    r[2] = _mm_loadu_si128((const __m128i *)(ref_ptr + 16));
    r[3] = _mm_loadu_si128((const __m128i *)(ref_ptr + 24));

    highbd_sad8x4_core_sse2(s, r, sad_acc);
    src_ptr += src_stride;
    ref_ptr += ref_stride;
  }
}

static void sad64x1(const uint16_t *src_ptr, const uint16_t *ref_ptr,
                    __m128i *sad_acc) {
  __m128i s[4], r[4];
  int i;
  for (i = 0; i < 2; i++) {
    s[0] = _mm_loadu_si128((const __m128i *)src_ptr);
    s[1] = _mm_loadu_si128((const __m128i *)(src_ptr + 8));
    s[2] = _mm_loadu_si128((const __m128i *)(src_ptr + 16));
    s[3] = _mm_loadu_si128((const __m128i *)(src_ptr + 24));
    r[0] = _mm_loadu_si128((const __m128i *)ref_ptr);
    r[1] = _mm_loadu_si128((const __m128i *)(ref_ptr + 8));
    r[2] = _mm_loadu_si128((const __m128i *)(ref_ptr + 16));
    r[3] = _mm_loadu_si128((const __m128i *)(ref_ptr + 24));

    highbd_sad8x4_core_sse2(s, r, sad_acc);
    src_ptr += 32;
    ref_ptr += 32;
  }
}

// Combine 4 __m256i vectors to uint32_t result[4]
static INLINE void get_4d_sad_from_mm128_epi32(const __m128i *v,
                                               int32_t *sad_err,
                                               uint32_t *res) {
  __m128i sad, err;
  __m128 tmp3, tmp2, tmp1, tmp0;
  __m128i u0, u1, u2, u3;
  // transpose 4x4 uint32_t arrays and sum the results
  tmp0 = _mm_unpacklo_ps(_mm_castsi128_ps(v[0]), _mm_castsi128_ps(v[1]));
  tmp2 = _mm_unpacklo_ps(_mm_castsi128_ps(v[2]), _mm_castsi128_ps(v[3]));
  tmp1 = _mm_unpackhi_ps(_mm_castsi128_ps(v[0]), _mm_castsi128_ps(v[1]));
  tmp3 = _mm_unpackhi_ps(_mm_castsi128_ps(v[2]), _mm_castsi128_ps(v[3]));
  u0 = _mm_castps_si128(_mm_movelh_ps(tmp0, tmp2));
  u1 = _mm_castps_si128(_mm_movehl_ps(tmp2, tmp0));
  u2 = _mm_castps_si128(_mm_movelh_ps(tmp1, tmp3));
  u3 = _mm_castps_si128(_mm_movehl_ps(tmp3, tmp1));

  u0 = _mm_add_epi32(u0, u1);
  u0 = _mm_add_epi32(u0, u2);
  u0 = _mm_add_epi32(u0, u3);

  err = _mm_load_si128((__m128i *)sad_err);
  sad = _mm_add_epi32(u0, err);

  _mm_storeu_si128((__m128i *)res, sad);
}

static void convert_pointers(const uint8_t *const ref8[],
                             const uint16_t *ref[]) {
  ref[0] = CONVERT_TO_SHORTPTR(ref8[0]);
  ref[1] = CONVERT_TO_SHORTPTR(ref8[1]);
  ref[2] = CONVERT_TO_SHORTPTR(ref8[2]);
  ref[3] = CONVERT_TO_SHORTPTR(ref8[3]);
}

static void init_sad(__m128i *s) {
  s[0] = _mm_setzero_si128();
  s[1] = _mm_setzero_si128();
  s[2] = _mm_setzero_si128();
  s[3] = _mm_setzero_si128();
}

static AOM_FORCE_INLINE void aom_highbd_sad4xNx4d_sse2(
    int N, const uint8_t *src, int src_stride, const uint8_t *ref_array[],
    int ref_stride, uint32_t *sad_array, int32_t *err) {
  __m128i sad_vec[4];
  const uint16_t *refp[4];
  const uint16_t *keep = CONVERT_TO_SHORTPTR(src);
  const uint16_t *srcp;
  const int shift_for_4_rows = 2;
  int i, j;

  init_sad(sad_vec);
  convert_pointers(ref_array, refp);

  for (i = 0; i < 4; ++i) {
    srcp = keep;
    for (j = 0; j < N; j += 4) {
      sad4x4(srcp, src_stride, refp[i], ref_stride, &sad_vec[i]);
      srcp += src_stride << shift_for_4_rows;
      refp[i] += ref_stride << shift_for_4_rows;
    }
  }
  get_4d_sad_from_mm128_epi32(sad_vec, err, sad_array);
}

static AOM_FORCE_INLINE void aom_highbd_sad8xNx4d_sse2(
    int N, const uint8_t *src, int src_stride, const uint8_t *ref_array[],
    int ref_stride, uint32_t *sad_array, int32_t *err) {
  __m128i sad_vec[4];
  const uint16_t *refp[4];
  const uint16_t *keep = CONVERT_TO_SHORTPTR(src);
  const uint16_t *srcp;
  const int shift_for_4_rows = 2;
  int i, j;

  init_sad(sad_vec);
  convert_pointers(ref_array, refp);

  for (i = 0; i < 4; ++i) {
    srcp = keep;
    for (j = 0; j < N; j += 4) {
      sad8x4(srcp, src_stride, refp[i], ref_stride, &sad_vec[i]);
      srcp += src_stride << shift_for_4_rows;
      refp[i] += ref_stride << shift_for_4_rows;
    }
  }
  get_4d_sad_from_mm128_epi32(sad_vec, err, sad_array);
}

static AOM_FORCE_INLINE void aom_highbd_sad16xNx4d_sse2(
    int N, const uint8_t *src, int src_stride, const uint8_t *ref_array[],
    int ref_stride, uint32_t *sad_array, int32_t *err) {
  __m128i sad_vec[4];
  const uint16_t *refp[4];
  const uint16_t *keep = CONVERT_TO_SHORTPTR(src);
  const uint16_t *srcp;
  const int shift_for_4_rows = 2;
  int i, r;

  init_sad(sad_vec);
  convert_pointers(ref_array, refp);

  for (i = 0; i < 4; ++i) {
    srcp = keep;
    for (r = 0; r < N; r += 4) {
      sad16x2(srcp, src_stride, refp[i], ref_stride, &sad_vec[i]);
      srcp += src_stride << shift_for_4_rows;
      refp[i] += ref_stride << shift_for_4_rows;
    }
  }
  get_4d_sad_from_mm128_epi32(sad_vec, err, sad_array);
}

static AOM_FORCE_INLINE void aom_highbd_sad32xNx4d_sse2(
    int N, const uint8_t *src, int src_stride, const uint8_t *ref_array[],
    int ref_stride, uint32_t *sad_array, int32_t *err) {
  __m128i sad_vec[4];
  const uint16_t *refp[4];
  const uint16_t *keep = CONVERT_TO_SHORTPTR(src);
  const uint16_t *srcp;
  const int shift_for_rows = 1;
  int i, r;

  init_sad(sad_vec);
  convert_pointers(ref_array, refp);

  for (i = 0; i < 4; ++i) {
    srcp = keep;
    for (r = 0; r < N; r += 2) {
      sad32x2(srcp, src_stride, refp[i], ref_stride, &sad_vec[i]);
      srcp += src_stride << shift_for_rows;
      refp[i] += ref_stride << shift_for_rows;
    }
  }
  get_4d_sad_from_mm128_epi32(sad_vec, err, sad_array);
}

static AOM_FORCE_INLINE void aom_highbd_sad64xNx4d_sse2(
    int N, const uint8_t *src, int src_stride, const uint8_t *ref_array[],
    int ref_stride, uint32_t *sad_array, int32_t *err) {
  __m128i sad_vec[4];
  const uint16_t *refp[4];
  const uint16_t *keep = CONVERT_TO_SHORTPTR(src);
  const uint16_t *srcp;
  int i, r;

  init_sad(sad_vec);
  convert_pointers(ref_array, refp);

  for (i = 0; i < 4; ++i) {
    srcp = keep;
    for (r = 0; r < N; r++) {
      sad64x1(srcp, refp[i], &sad_vec[i]);
      srcp += src_stride;
      refp[i] += ref_stride;
    }
  }
  get_4d_sad_from_mm128_epi32(sad_vec, err, sad_array);
}

#define highbd_sadMxNx4d_sse2(m, n)                                          \
  void aom_highbd_sad##m##x##n##x4d_sse2(                                    \
      const uint8_t *src, int src_stride, const uint8_t *ref_array[],        \
      int ref_stride, int32_t *err, uint32_t *min_value, int32_t *min_pos) { \
    DECLARE_ALIGNED(16, uint32_t, sad_array[4]);                             \
    aom_highbd_sad##m##xNx4d_sse2(n, src, src_stride, ref_array, ref_stride, \
                                  sad_array, err);                           \
    get_min_value_pos(sad_array, 4, min_value, min_pos);                     \
  }

highbd_sadMxNx4d_sse2(4, 4);
highbd_sadMxNx4d_sse2(4, 8);
highbd_sadMxNx4d_sse2(4, 16);

highbd_sadMxNx4d_sse2(8, 4);
highbd_sadMxNx4d_sse2(8, 8);
highbd_sadMxNx4d_sse2(8, 16);
highbd_sadMxNx4d_sse2(8, 32);

highbd_sadMxNx4d_sse2(16, 4);
highbd_sadMxNx4d_sse2(16, 8);
highbd_sadMxNx4d_sse2(16, 16);
highbd_sadMxNx4d_sse2(16, 32);
highbd_sadMxNx4d_sse2(16, 64);

highbd_sadMxNx4d_sse2(32, 8);
highbd_sadMxNx4d_sse2(32, 16);
highbd_sadMxNx4d_sse2(32, 32);
highbd_sadMxNx4d_sse2(32, 64);

highbd_sadMxNx4d_sse2(64, 16);
highbd_sadMxNx4d_sse2(64, 32);
highbd_sadMxNx4d_sse2(64, 64);
