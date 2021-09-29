/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <emmintrin.h>  // SSE2

#include "config/av1_rtcd.h"

#include "aom/aom_integer.h"

static AOM_INLINE __m128i reduce_sum_epi64(__m128i reg) {
  __m128i reg_hi = _mm_srli_si128(reg, 8);
  reg = _mm_add_epi64(reg, reg_hi);

  return reg;
}

int64_t av1_block_error_lp_sse2(const int16_t *coeff, const int16_t *dqcoeff,
                                intptr_t block_size) {
  assert(block_size % 8 == 0);
  assert(block_size >= 8);

  const __m128i zero = _mm_setzero_si128();
  __m128i accum = zero;

  for (int i = 0; i < block_size; i += 8) {
    // Load 8 elements for coeff and dqcoeff.
    const __m128i _coeff = _mm_loadu_si128((const __m128i *)coeff);
    const __m128i _dqcoeff = _mm_loadu_si128((const __m128i *)dqcoeff);
    // Compute the diff
    const __m128i diff = _mm_sub_epi16(_dqcoeff, _coeff);
    // Compute the error
    const __m128i error = _mm_madd_epi16(diff, diff);
    const __m128i error_lo = _mm_unpacklo_epi32(error, zero);
    const __m128i error_hi = _mm_unpackhi_epi32(error, zero);

    // Accumulate
    accum = _mm_add_epi64(accum, error_lo);
    accum = _mm_add_epi64(accum, error_hi);

    // Advance
    coeff += 8;
    dqcoeff += 8;
  }

  // Reduce sum the register
  accum = reduce_sum_epi64(accum);

  // Store the results.
#if ARCH_X86_64
  return _mm_cvtsi128_si64(accum);
#else
  int64_t result;
  _mm_storeu_si64(&result, accum);
  return result;
#endif  // ARCH_X86_64
}
