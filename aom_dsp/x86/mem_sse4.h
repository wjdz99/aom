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

#ifndef AOM_DSP_X86_MEM_SSE4_H_
#define AOM_DSP_X86_MEM_SSE4_H_

#include <smmintrin.h>  // SSE4.1

#include "./aom_config.h"
#include "aom/aom_integer.h"

static INLINE __m128i load_8bit_4x2_to_1_sse4_1(const uint8_t *const src,
                                                const ptrdiff_t stride) {
  const __m128i s = _mm_cvtsi32_si128(*(const int *)(src + 0 * stride));
  return _mm_insert_epi32(s, *(const int *)(src + 1 * stride), 1);
}

#endif  // AOM_DSP_X86_MEM_SSE4_H_
