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

#ifndef AOM_DSP_X86_INV_TXFM_COMMON_AVX2_H
#define AOM_DSP_X86_INV_TXFM_COMMON_AVX2_H

#include <immintrin.h>

#include "aom_dsp/txfm_common.h"
#include "aom_dsp/x86/txfm_common_avx2.h"

static INLINE void unpack_butter_fly(const __m256i *a0, const __m256i *a1,
                                     const __m256i *c0, const __m256i *c1,
                                     __m256i *b0, __m256i *b1) {
  __m256i x0, x1;
  x0 = _mm256_unpacklo_epi16(*a0, *a1);
  x1 = _mm256_unpackhi_epi16(*a0, *a1);
  *b0 = butter_fly(&x0, &x1, c0);
  *b1 = butter_fly(&x0, &x1, c1);
}

#endif  // AOM_DSP_X86_INV_TXFM_COMMON_AVX2_H
