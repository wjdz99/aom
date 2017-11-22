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

#ifndef AV1_DECODER_X86_DECODETXB_SSE2_H_
#define AV1_DECODER_X86_DECODETXB_SSE2_H_

#include <emmintrin.h>  // SSE2

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/x86/mem_sse2.h"
#include "av1/common/enums.h"

static INLINE void construct_dqvs_sse2(const int16_t *const dequant,
                                       __m128i *const dqvs /*dqvs[2]*/) {
  dqvs[1] = _mm_set1_epi32(dequant[1]);
  dqvs[0] = _mm_slli_si128(dqvs[1], 4);
  dqvs[0] = _mm_or_si128(dqvs[0], _mm_cvtsi32_si128(dequant[0]));
}

#endif  // AV1_DECODER_X86_DECODETXB_SSE2_H_
