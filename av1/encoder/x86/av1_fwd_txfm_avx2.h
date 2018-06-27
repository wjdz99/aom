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

#ifndef AV1_FWD_TXFM_AVX2_H_
#define AV1_FWD_TXFM_AVX2_H_
#include <immintrin.h>

static INLINE __m256i av1_round_shift_32_avx2(__m256i vec, int bit) {
  __m256i tmp, round;
  round = _mm256_set1_epi32(1 << (bit - 1));
  tmp = _mm256_add_epi32(vec, round);
  return _mm256_srai_epi32(tmp, bit);
}

// out0 = in0*w0 + in1*w1
// out1 = -in1*w0 + in0*w1
#define btf_32_avx2_type0(w0, w1, in0, in1, out0, out1)   \
  do {                                                    \
    __m256i _in0 = in0;                                   \
    __m256i _in1 = in1;                                   \
    const __m256i ww0 = _mm256_set1_epi32(w0);            \
    const __m256i ww1 = _mm256_set1_epi32(w1);            \
    const __m256i in0_w0 = _mm256_mullo_epi32(_in0, ww0); \
    const __m256i in1_w1 = _mm256_mullo_epi32(_in1, ww1); \
    out0 = _mm256_add_epi32(in0_w0, in1_w1);              \
    out0 = av1_round_shift_32_avx2(out0, cos_bit);        \
    const __m256i in0_w1 = _mm256_mullo_epi32(_in0, ww1); \
    const __m256i in1_w0 = _mm256_mullo_epi32(_in1, ww0); \
    out1 = _mm256_sub_epi32(in0_w1, in1_w0);              \
    out1 = av1_round_shift_32_avx2(out1, cos_bit);        \
  } while (0)

#define btf_32_avx2_type1(w0, w1, in0, in1, out0, out1) \
  do {                                                  \
    btf_32_avx2_type0(w1, w0, in1, in0, out0, out1);    \
  } while (0)

// out0 = in0*w0 + in1*w1
// out1 = -in1*w0 + in0*w1
#define btf_32_type0_avx2_new(ww0, ww1, in0, in1, out0, out1) \
  do {                                                        \
    __m256i _in0 = in0;                                       \
    __m256i _in1 = in1;                                       \
    const __m256i in0_w0 = _mm256_mullo_epi32(_in0, ww0);     \
    const __m256i in1_w1 = _mm256_mullo_epi32(_in1, ww1);     \
    out0 = _mm256_add_epi32(in0_w0, in1_w1);                  \
    out0 = _mm256_add_epi32(out0, __rounding);                \
    out0 = _mm256_srai_epi32(out0, cos_bit);                  \
    const __m256i in0_w1 = _mm256_mullo_epi32(_in0, ww1);     \
    const __m256i in1_w0 = _mm256_mullo_epi32(_in1, ww0);     \
    out1 = _mm256_sub_epi32(in0_w1, in1_w0);                  \
    out1 = _mm256_add_epi32(out1, __rounding);                \
    out1 = _mm256_srai_epi32(out1, cos_bit);                  \
  } while (0)

// out0 = in0*w0 + in1*w1
// out1 = in1*w0 - in0*w1
#define btf_32_type1_avx2_new(ww0, ww1, in0, in1, out0, out1) \
  do {                                                        \
    btf_32_type0_avx2_new(ww1, ww0, in1, in0, out0, out1);    \
  } while (0)

#endif  // AV1_FWD_TXFM_AVX2_H_
