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
#include <immintrin.h>  // AVX2

#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"

static AOM_INLINE void aom_sadx4d_get_min(uint32_t src[4], uint32_t *min_value,
                                          int32_t *min_pos) {
  int t;
  for (t = 0; t < 4; t++) {
    if (src[t] < *min_value) {
      *min_value = src[t];
      *min_pos = t;
    }
  }
}

static AOM_INLINE void aom_sadx4d_get_min16(__m128i src, uint32_t *min_value,
                                            int32_t *min_pos) {
  DECLARE_ALIGNED(16, uint16_t, mask[8]) = { 0, 0xffff, 0, 0xffff,
                                             0, 0xffff, 0, 0xffff };
  __m128i src_32 = _mm_or_si128(src, *(__m128i *)mask);
  __m128i minpos = _mm_minpos_epu16(src_32);
  int res = _mm_cvtsi128_si32(minpos);
  (*min_value) = (res & 0xffff);
  (*min_pos) = res >> 17;
}

AOM_FORCE_INLINE void aom_sad16xNx4d_avx2(int N, const uint8_t *src,
                                          int src_stride, const uint8_t *ref[4],
                                          int ref_stride, uint32_t res[4],
                                          uint32_t err[4], uint32_t *min_value,
                                          int32_t *min_pos) {
  __m128i src_reg, ref0_reg, ref3_reg, ref1_reg, ref2_reg;
  __m128i sum_ref0, sum_ref2, sum_ref1, sum_ref3;
  int i;

  sum_ref0 = _mm_setzero_si128();
  sum_ref2 = _mm_setzero_si128();
  sum_ref1 = _mm_setzero_si128();
  sum_ref3 = _mm_setzero_si128();
  for (i = 0; i < N; i++) {
    // load src and all refs
    src_reg = _mm_load_si128((const __m128i *)(src));
    ref0_reg = _mm_loadu_si128((const __m128i *)(ref[0]));
    ref1_reg = _mm_loadu_si128((const __m128i *)(ref[1]));
    ref2_reg = _mm_loadu_si128((const __m128i *)(ref[2]));
    ref3_reg = _mm_loadu_si128((const __m128i *)(ref[3]));

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

    src += src_stride;
    ref[0] += ref_stride;
    ref[1] += ref_stride;
    ref[2] += ref_stride;
    ref[3] += ref_stride;
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

    __m128i err_reg = _mm_load_si128((const __m128i *)err);
    sum_mlow = _mm_add_epi32(sum_mlow, err_reg);

    _mm_store_si128((__m128i *)(res), sum_mlow);

    if (N < 32) {
      aom_sadx4d_get_min16(sum_mlow, min_value, min_pos);
    } else {
      aom_sadx4d_get_min(res, min_value, min_pos);
    }
  }
}

AOM_FORCE_INLINE void aom_sadMxNx4d_avx2(int M, int N, const uint8_t *src,
                                         int src_stride, const uint8_t *ref[4],
                                         int ref_stride, uint32_t res[4],
                                         uint32_t err[4], uint32_t *min_value,
                                         int32_t *min_pos) {
  __m256i src_reg, ref0_reg, ref1_reg, ref2_reg, ref3_reg;
  __m256i sum_ref0, sum_ref1, sum_ref2, sum_ref3;
  int i, j;

  sum_ref0 = _mm256_setzero_si256();
  sum_ref2 = _mm256_setzero_si256();
  sum_ref1 = _mm256_setzero_si256();
  sum_ref3 = _mm256_setzero_si256();

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j += 32) {
      // load src and all refs
      src_reg = _mm256_load_si256((const __m256i *)(src + j));
      ref0_reg = _mm256_loadu_si256((const __m256i *)(ref[0] + j));
      ref1_reg = _mm256_loadu_si256((const __m256i *)(ref[1] + j));
      ref2_reg = _mm256_loadu_si256((const __m256i *)(ref[2] + j));
      ref3_reg = _mm256_loadu_si256((const __m256i *)(ref[3] + j));

      // sum of the absolute differences between every ref-i to src
      ref0_reg = _mm256_sad_epu8(ref0_reg, src_reg);
      ref1_reg = _mm256_sad_epu8(ref1_reg, src_reg);
      ref2_reg = _mm256_sad_epu8(ref2_reg, src_reg);
      ref3_reg = _mm256_sad_epu8(ref3_reg, src_reg);
      // sum every ref-i
      sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
      sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
      sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
      sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);
    }
    src += src_stride;
    ref[0] += ref_stride;
    ref[1] += ref_stride;
    ref[2] += ref_stride;
    ref[3] += ref_stride;
  }
  {
    __m128i sum;
    __m256i sum_mlow, sum_mhigh;
    // in sum_ref-i the result is saved in the first 4 bytes
    // the other 4 bytes are zeroed.
    // sum_ref1 and sum_ref3 are shifted left by 4 bytes
    sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
    sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

    // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
    sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
    sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

    // merge every 64 bit from each sum_ref-i
    sum_mlow = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
    sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

    // add the low 64 bit to the high 64 bit
    sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

    // add the low 128 bit to the high 128 bit
    sum = _mm_add_epi32(_mm256_castsi256_si128(sum_mlow),
                        _mm256_extractf128_si256(sum_mlow, 1));

    __m128i err_reg = _mm_load_si128((const __m128i *)err);
    sum = _mm_add_epi32(sum, err_reg);

    _mm_store_si128((__m128i *)(res), sum);
    if (N < 32) {
      aom_sadx4d_get_min16(sum, min_value, min_pos);
    } else {
      aom_sadx4d_get_min(res, min_value, min_pos);
    }
  }
}

#define sad16xNx4d_avx2(n)                                              \
  void aom_sad16x##n##x4d_avx2(const uint8_t *src, int src_stride,      \
                               const uint8_t *ref[4], int ref_stride,   \
                               uint32_t res[4], uint32_t err[4],        \
                               uint32_t *min_value, int32_t *min_pos) { \
    aom_sad16xNx4d_avx2(n, src, src_stride, ref, ref_stride, res, err,  \
                        min_value, min_pos);                            \
  }

#define sadMxNx4d_avx2(m, n)                                               \
  void aom_sad##m##x##n##x4d_avx2(const uint8_t *src, int src_stride,      \
                                  const uint8_t *ref[4], int ref_stride,   \
                                  uint32_t res[4], uint32_t err[4],        \
                                  uint32_t *min_value, int32_t *min_pos) { \
    aom_sadMxNx4d_avx2(m, n, src, src_stride, ref, ref_stride, res, err,   \
                       min_value, min_pos);                                \
  }

sad16xNx4d_avx2(4);
sad16xNx4d_avx2(8);
sad16xNx4d_avx2(16);
sad16xNx4d_avx2(32);
sad16xNx4d_avx2(64);

sadMxNx4d_avx2(32, 8);
sadMxNx4d_avx2(32, 16);
sadMxNx4d_avx2(32, 32);
sadMxNx4d_avx2(32, 64);

sadMxNx4d_avx2(64, 16);
sadMxNx4d_avx2(64, 32);
sadMxNx4d_avx2(64, 64);
sadMxNx4d_avx2(64, 128);

sadMxNx4d_avx2(128, 64);
sadMxNx4d_avx2(128, 128);
