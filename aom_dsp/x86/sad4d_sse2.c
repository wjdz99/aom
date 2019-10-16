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

#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"

AOM_FORCE_INLINE void aom_sad4xNx4d_sse2(int N, const uint8_t *src,
                                         int src_stride,
                                         const uint8_t *const ref[4],
                                         int ref_stride, uint32_t res[4]) {
  __m128i src_reg, ref0_reg, ref3_reg, ref1_reg, ref2_reg;
  __m128i ref02_reg, ref13_reg;
  __m128i sum_ref0, sum_ref2;
  int i;
  const uint8_t *ref0, *ref1, *ref2, *ref3;

  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];
  sum_ref0 = _mm_setzero_si128();
  sum_ref2 = _mm_setzero_si128();
  for (i = 0; i < N; i++) {
    // load src and all refs
    src_reg = _mm_cvtsi32_si128(*((uint32_t *)src));
    src_reg = _mm_unpacklo_epi64(src_reg, src_reg);
    ref0_reg = _mm_cvtsi32_si128(*((uint32_t *)ref0));
    ref1_reg = _mm_cvtsi32_si128(*((uint32_t *)ref1));
    ref2_reg = _mm_cvtsi32_si128(*((uint32_t *)ref2));
    ref3_reg = _mm_cvtsi32_si128(*((uint32_t *)ref3));
    ref02_reg = _mm_unpacklo_epi64(ref0_reg, ref2_reg);
    ref13_reg = _mm_unpacklo_epi64(ref1_reg, ref3_reg);

    // sum of the absolute differences between every ref-i to src
    ref02_reg = _mm_sad_epu8(ref02_reg, src_reg);
    ref13_reg = _mm_sad_epu8(ref13_reg, src_reg);

    // sum every ref-i
    sum_ref0 = _mm_add_epi32(sum_ref0, ref02_reg);
    sum_ref2 = _mm_add_epi32(sum_ref2, ref13_reg);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }
  {
    sum_ref2 = _mm_slli_si128(sum_ref2, 4);
    // merge sum_ref0 and sum_ref2
    sum_ref0 = _mm_or_si128(sum_ref0, sum_ref2);
    _mm_storeu_si128((__m128i *)(res), sum_ref0);
  }
}

AOM_FORCE_INLINE void aom_sad8xNx4d_sse2(int N, const uint8_t *src,
                                         int src_stride,
                                         const uint8_t *const ref[4],
                                         int ref_stride, uint32_t res[4]) {
  __m128i src_reg, ref0_reg, ref3_reg, ref1_reg, ref2_reg;
  __m128i ref02_reg, ref13_reg;
  __m128i sum_ref0, sum_ref2;
  int i;
  const uint8_t *ref0, *ref1, *ref2, *ref3;

  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];
  sum_ref0 = _mm_setzero_si128();
  sum_ref2 = _mm_setzero_si128();
  for (i = 0; i < N; i++) {
    // load src and all refs
    src_reg = _mm_loadl_epi64((const __m128i *)src);
    src_reg = _mm_unpacklo_epi64(src_reg, src_reg);
    ref0_reg = _mm_loadl_epi64((const __m128i *)ref0);
    ref1_reg = _mm_loadl_epi64((const __m128i *)ref1);
    ref2_reg = _mm_loadl_epi64((const __m128i *)ref2);
    ref3_reg = _mm_loadl_epi64((const __m128i *)ref3);
    ref02_reg = _mm_unpacklo_epi64(ref0_reg, ref2_reg);
    ref13_reg = _mm_unpacklo_epi64(ref1_reg, ref3_reg);

    // sum of the absolute differences between every ref-i to src
    ref02_reg = _mm_sad_epu8(ref02_reg, src_reg);
    ref13_reg = _mm_sad_epu8(ref13_reg, src_reg);

    // sum every ref-i
    sum_ref0 = _mm_add_epi32(sum_ref0, ref02_reg);
    sum_ref2 = _mm_add_epi32(sum_ref2, ref13_reg);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }
  {
    sum_ref2 = _mm_slli_si128(sum_ref2, 4);
    // merge sum_ref0 and sum_ref2
    sum_ref0 = _mm_or_si128(sum_ref0, sum_ref2);
    _mm_storeu_si128((__m128i *)(res), sum_ref0);
  }
}

AOM_FORCE_INLINE void aom_sadMxNx4d_sse2(int M, int N, const uint8_t *src,
                                         int src_stride,
                                         const uint8_t *const ref[4],
                                         int ref_stride, uint32_t res[4]) {
  __m128i src_reg, ref0_reg, ref3_reg, ref1_reg, ref2_reg;
  __m128i sum_ref0, sum_ref2, sum_ref1, sum_ref3;
  int i, j;
  const uint8_t *ref0, *ref1, *ref2, *ref3;

  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  sum_ref0 = _mm_setzero_si128();
  sum_ref2 = _mm_setzero_si128();
  sum_ref1 = _mm_setzero_si128();
  sum_ref3 = _mm_setzero_si128();
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j += 16) {
    // load src and all refs
      src_reg = _mm_loadu_si128((const __m128i *)(src + j));
      ref0_reg = _mm_loadu_si128((const __m128i *)(ref0 + j));
      ref1_reg = _mm_loadu_si128((const __m128i *)(ref1 + j));
      ref2_reg = _mm_loadu_si128((const __m128i *)(ref2 + j));
      ref3_reg = _mm_loadu_si128((const __m128i *)(ref3 + j));

      // sum of the absolute differences between every ref-i to src
      ref0_reg = _mm_sad_epu8(ref0_reg, src_reg);
      ref1_reg = _mm_sad_epu8(ref1_reg, src_reg);
      ref2_reg = _mm_sad_epu8(ref2_reg, src_reg);
      ref3_reg = _mm_sad_epu8(ref3_reg, src_reg);

      // sum every ref-i
      sum_ref0 = _mm_add_epi32(sum_ref0, ref0_reg);
      sum_ref1 = _mm_add_epi32(sum_ref1, ref1_reg);
      sum_ref2 = _mm_add_epi32(sum_ref2, ref2_reg);
      sum_ref3 = _mm_add_epi32(sum_ref3, ref3_reg);
    }
    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }
  {
    __m128i sum_mlow, sum_mhigh;
    sum_ref1 = _mm_slli_si128(sum_ref1, 4);
    sum_ref3 = _mm_slli_si128(sum_ref3, 4);
    // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
    sum_ref0 = _mm_or_si128(sum_ref0, sum_ref1);
    sum_ref2 = _mm_or_si128(sum_ref2, sum_ref3);

    // merge every 64 bit from each sum_ref-i
    sum_mlow = _mm_unpacklo_epi64(sum_ref0, sum_ref2);
    sum_mhigh = _mm_unpackhi_epi64(sum_ref0, sum_ref2);

    sum_mlow = _mm_add_epi32(sum_mlow, sum_mhigh);
    _mm_storeu_si128((__m128i *)(res), sum_mlow);
  }
}

#define sad4xN_sse2(n)                                                     \
  void aom_sad4x##n##x4d_sse2(const uint8_t *src, int src_stride,          \
                              const uint8_t *const ref[4], int ref_stride, \
                              uint32_t res[4]) {                           \
    aom_sad4xNx4d_sse2(n, src, src_stride, ref, ref_stride, res);          \
  }

#define sad8xN_sse2(n)                                                     \
  void aom_sad8x##n##x4d_sse2(const uint8_t *src, int src_stride,          \
                              const uint8_t *const ref[4], int ref_stride, \
                              uint32_t res[4]) {                           \
    aom_sad8xNx4d_sse2(n, src, src_stride, ref, ref_stride, res);          \
  }

#define sadMxN_sse2(m, n)                                                      \
  void aom_sad##m##x##n##x4d_sse2(const uint8_t *src, int src_stride,          \
                                  const uint8_t *const ref[4], int ref_stride, \
                                  uint32_t res[4]) {                           \
    aom_sadMxNx4d_sse2(m, n, src, src_stride, ref, ref_stride, res);           \
  }

sad4xN_sse2(4);
sad4xN_sse2(8);
sad4xN_sse2(16);

sad8xN_sse2(4);
sad8xN_sse2(8);
sad8xN_sse2(16);
sad8xN_sse2(32);

sadMxN_sse2(16, 4);
sadMxN_sse2(16, 8);
sadMxN_sse2(16, 16);
sadMxN_sse2(16, 32);
sadMxN_sse2(16, 64);

sadMxN_sse2(32, 8);
sadMxN_sse2(32, 16);
sadMxN_sse2(32, 32);
sadMxN_sse2(32, 64);

sadMxN_sse2(64, 16);
sadMxN_sse2(64, 32);
sadMxN_sse2(64, 64);
sadMxN_sse2(64, 128);

sadMxN_sse2(128, 64);
sadMxN_sse2(128, 128);
