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

#ifndef DAALA_TX_KERNELS_AVX2_H_
#define DAALA_TX_KERNELS_AVX2_H_

#if CONFIG_DAALA_TX

static INLINE __m128i od_unbiased_rshift1_epi16(__m128i a) {
  return _mm_srai_epi16(_mm_add_epi16(_mm_srli_epi16(a, 15), a), 1);
}

static INLINE __m256i od_mm256_unbiased_rshift1_epi16(__m256i a) {
  return _mm256_srai_epi16(_mm256_add_epi16(_mm256_srli_epi16(a, 15), a), 1);
}

static INLINE __m256i od_mm256_unbiased_rshift1_epi32(__m256i a) {
  return _mm256_srai_epi32(_mm256_add_epi32(_mm256_srli_epi32(a, 31), a), 1);
}

static INLINE __m128i od_avg_epi16(__m128i a, __m128i b) {
  __m128i sign_bit;
  /*x86 only provides an unsigned PAVGW with a bias (ARM is better here).
    We emulate a signed one by adding an offset to convert to unsigned and
    back. We use XOR instead of addition/subtraction because it dispatches
    better on older processors.*/
  sign_bit = _mm_set1_epi16(0x8000);
  return _mm_xor_si128(
      _mm_avg_epu16(_mm_xor_si128(a, sign_bit), _mm_xor_si128(b, sign_bit)),
      sign_bit);
}

static INLINE __m256i od_mm256_avg_epi16(__m256i a, __m256i b) {
  __m256i sign_bit;
  sign_bit = _mm256_set1_epi16(0x8000);
  return _mm256_xor_si256(_mm256_avg_epu16(_mm256_xor_si256(a, sign_bit),
                                           _mm256_xor_si256(b, sign_bit)),
                          sign_bit);
}

static INLINE __m256i od_mm256_avg_epi32(__m256i a, __m256i b) {
  __m256i neg1;
  /* It's cheaper to generate -1's than 1's. */
  neg1 = _mm256_set1_epi64x(-1);
  /* There is no corresponding PAVGD, but we are not in danger of overflowing
     a 32-bit register. */
  return _mm256_srai_epi32(_mm256_add_epi32(a, _mm256_sub_epi32(b, neg1)), 1);
}

/*Like the above, but does (a - b + 1) >> 1 instead.*/
static INLINE __m128i od_hrsub_epi16(__m128i a, __m128i b) {
  __m128i sign_bit;
  sign_bit = _mm_set1_epi16(0x8000);
  return _mm_xor_si128(
      _mm_avg_epu16(_mm_xor_si128(a, sign_bit), _mm_sub_epi16(sign_bit, b)),
      sign_bit);
}

static INLINE __m256i od_mm256_hrsub_epi16(__m256i a, __m256i b) {
  __m256i sign_bit;
  sign_bit = _mm256_set1_epi16(0x8000);
  return _mm256_xor_si256(_mm256_avg_epu16(_mm256_xor_si256(a, sign_bit),
                                           _mm256_sub_epi16(sign_bit, b)),
                          sign_bit);
}

static INLINE __m256i od_mm256_hrsub_epi32(__m256i a, __m256i b) {
  __m256i neg1;
  /* It's cheaper to generate -1's than 1's. */
  neg1 = _mm256_set1_epi64x(-1);
  /* There is no corresponding PAVGD, but we are not in danger of overflowing
     a 32-bit register. */
  return _mm256_srai_epi32(_mm256_sub_epi32(a, _mm256_add_epi32(b, neg1)), 1);
}

static INLINE void od_swap_si128(__m128i *q0, __m128i *q1) {
  __m128i t;
  t = *q0;
  *q0 = *q1;
  *q1 = t;
}

static INLINE void od_mm256_swap_si256(__m256i *q0, __m256i *q1) {
  __m256i t;
  t = *q0;
  *q0 = *q1;
  *q1 = t;
}

static INLINE __m128i od_mulhrs_epi16(__m128i a, int16_t b) {
  return _mm_mulhrs_epi16(a, _mm_set1_epi16(b));
}

static INLINE __m128i od_mul_epi16(__m128i a, int32_t b, int r) {
  int32_t b_q15;
  b_q15 = b << (15 - r);
  /* b and r are in all cases compile-time constants, so these branches
     disappear when this function gets inlined. */
  if (b_q15 > 32767) {
    return _mm_add_epi16(a, od_mulhrs_epi16(a, (int16_t)(b_q15 - 32768)));
  } else if (b_q15 < -32767) {
    return _mm_sub_epi16(od_mulhrs_epi16(a, (int16_t)(32768 + b_q15)), a);
  } else {
    return od_mulhrs_epi16(a, b_q15);
  }
}

static INLINE __m256i od_mm256_mulhrs_epi16(__m256i a, int16_t b) {
  return _mm256_mulhrs_epi16(a, _mm256_set1_epi16(b));
}

static INLINE __m256i od_mm256_mul_epi16(__m256i a, int32_t b, int r) {
  int32_t b_q15;
  b_q15 = b << (15 - r);
  /* b and r are in all cases compile-time constants, so these branches
     disappear when this function gets inlined. */
  if (b_q15 > 32767) {
    return _mm256_add_epi16(a,
                            od_mm256_mulhrs_epi16(a, (int16_t)(b_q15 - 32768)));
  } else if (b_q15 < -32767) {
    return _mm256_sub_epi16(od_mm256_mulhrs_epi16(a, (int16_t)(32768 + b_q15)),
                            a);
  } else {
    return od_mm256_mulhrs_epi16(a, b_q15);
  }
}

static INLINE __m256i od_mm256_mul_epi32(__m256i a, int32_t b, int r) {
  __m256i neg1;
  /* It's cheaper to generate -1's than 1's. */
  neg1 = _mm256_set1_epi64x(-1);
  /* There's no 32-bit version of PMULHRSW on x86 like there is on ARM .*/
  a = _mm256_mullo_epi32(a, _mm256_set1_epi32(b));
  a = _mm256_srai_epi32(a, r - 1);
  a = _mm256_sub_epi32(a, neg1);
  return _mm256_srai_epi32(a, 1);
}

#undef OD_KERNEL
#undef OD_WORD
#undef OD_REG
#undef OD_ADD
#undef OD_SUB
#undef OD_RSHIFT1
#undef OD_AVG
#undef OD_HRSUB
#undef OD_MUL
#undef OD_SWAP

/* Define 8-wide 16-bit SSSE3 kernels. */

#define OD_KERNEL kernel8
#define OD_WORD epi16
#define OD_REG __m128i
#define OD_ADD _mm_add_epi16
#define OD_SUB _mm_sub_epi16
#define OD_RSHIFT1 od_unbiased_rshift1_epi16
#define OD_AVG od_avg_epi16
#define OD_HRSUB od_hrsub_epi16
#define OD_MUL od_mul_epi16
#define OD_SWAP od_swap_si128

#include "av1/common/x86/daala_tx_kernels_impl.h"

#undef OD_KERNEL
#undef OD_REG
#undef OD_ADD
#undef OD_SUB
#undef OD_RSHIFT1
#undef OD_AVG
#undef OD_HRSUB
#undef OD_MUL
#undef OD_SWAP

/* Define 16-wide 16-bit AVX2 kernels. */

#define OD_KERNEL kernel16
#define OD_REG __m256i
#define OD_ADD _mm256_add_epi16
#define OD_SUB _mm256_sub_epi16
#define OD_RSHIFT1 od_mm256_unbiased_rshift1_epi16
#define OD_AVG od_mm256_avg_epi16
#define OD_HRSUB od_mm256_hrsub_epi16
#define OD_MUL od_mm256_mul_epi16
#define OD_SWAP od_mm256_swap_si256

#include "av1/common/x86/daala_tx_kernels_impl.h"  // NOLINT

/* Define 8-wide 32-bit AVX2 kernels. */

#undef OD_KERNEL
#undef OD_WORD
#undef OD_ADD
#undef OD_SUB
#undef OD_RSHIFT1
#undef OD_AVG
#undef OD_HRSUB
#undef OD_MUL

#define OD_KERNEL kernel8
#define OD_WORD epi32
#define OD_ADD _mm256_add_epi32
#define OD_SUB _mm256_sub_epi32
#define OD_RSHIFT1 od_mm256_unbiased_rshift1_epi32
#define OD_AVG od_mm256_avg_epi32
#define OD_HRSUB od_mm256_hrsub_epi32
#define OD_MUL od_mm256_mul_epi32

#include "av1/common/x86/daala_tx_kernels_impl.h"  // NOLINT

#undef OD_KERNEL
#undef OD_WORD
#undef OD_REG
#undef OD_ADD
#undef OD_SUB
#undef OD_RSHIFT1
#undef OD_AVG
#undef OD_HRSUB
#undef OD_MUL
#undef OD_SWAP

#endif

#endif // AV1_TXMF1D_SSE2_H_
