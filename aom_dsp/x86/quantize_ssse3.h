/*
 *  Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_X86_QUANTIZE_SSSE3_H_
#define AOM_AOM_DSP_X86_QUANTIZE_SSSE3_H_

#include <emmintrin.h>

#include "aom/aom_integer.h"
#include "config/aom_config.h"

static INLINE void calculate_dqcoeff_and_store_32x32(const __m128i qcoeff,
                                                     const __m128i dequant,
                                                     const __m128i zero,
                                                     tran_low_t *dqcoeff) {
  // Un-sign to bias rounding like C.
  const __m128i coeff = _mm_abs_epi16(qcoeff);

  const __m128i sign_0 = _mm_unpacklo_epi16(zero, qcoeff);
  const __m128i sign_1 = _mm_unpackhi_epi16(zero, qcoeff);

  const __m128i low = _mm_mullo_epi16(coeff, dequant);
  const __m128i high = _mm_mulhi_epi16(coeff, dequant);
  __m128i dqcoeff32_0 = _mm_unpacklo_epi16(low, high);
  __m128i dqcoeff32_1 = _mm_unpackhi_epi16(low, high);

  // "Divide" by 2.
  dqcoeff32_0 = _mm_srli_epi32(dqcoeff32_0, 1);
  dqcoeff32_1 = _mm_srli_epi32(dqcoeff32_1, 1);

  dqcoeff32_0 = _mm_sign_epi32(dqcoeff32_0, sign_0);
  dqcoeff32_1 = _mm_sign_epi32(dqcoeff32_1, sign_1);

  _mm_store_si128((__m128i *)(dqcoeff), dqcoeff32_0);
  _mm_store_si128((__m128i *)(dqcoeff + 4), dqcoeff32_1);
}

#endif  // AOM_AOM_DSP_X86_QUANTIZE_SSSE3_H_
