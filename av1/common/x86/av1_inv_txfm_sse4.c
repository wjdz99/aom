/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <smmintrin.h> /* SSE4.1 */

#include "config/aom_config.h"
#include "config/av1_rtcd.h"

#include "av1/common/x86/av1_txfm_sse2.h"

static INLINE __m128i highbd_clamp_epi16(__m128i u, int bd) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);
  const __m128i max = _mm_sub_epi16(_mm_slli_epi16(one, bd), one);
  __m128i clamped, mask;

  mask = _mm_cmpgt_epi16(u, max);
  clamped = _mm_andnot_si128(mask, u);
  mask = _mm_and_si128(mask, max);
  clamped = _mm_or_si128(mask, clamped);
  mask = _mm_cmpgt_epi16(clamped, zero);
  clamped = _mm_and_si128(clamped, mask);

  return clamped;
}

static INLINE void load_buffer_4x4(const int32_t *coeff, __m128i *in) {
  in[0] = _mm_load_si128((const __m128i *)(coeff + 0));
  in[1] = _mm_load_si128((const __m128i *)(coeff + 4));
  in[2] = _mm_load_si128((const __m128i *)(coeff + 8));
  in[3] = _mm_load_si128((const __m128i *)(coeff + 12));
}

// Currently, this is used to perform low bitdepth iwht in av1/common/idct.c
void av1_highbd_iwht4x4_16_add_sse4_1(const tran_low_t *input, uint8_t *dest8,
                                      int stride, int bd) {
  /* 4-point reversible, orthonormal inverse Walsh-Hadamard in 3.5 adds,
     0.5 shifts per pixel. */
  __m128i op[4];
  uint16_t *dest = CONVERT_TO_SHORTPTR(dest8);

  load_buffer_4x4(input, op);

  // Shift before-hand.
  op[0] = _mm_srai_epi32(op[0], UNIT_QUANT_SHIFT);
  op[1] = _mm_srai_epi32(op[1], UNIT_QUANT_SHIFT);
  op[2] = _mm_srai_epi32(op[2], UNIT_QUANT_SHIFT);
  op[3] = _mm_srai_epi32(op[3], UNIT_QUANT_SHIFT);

  for (int i = 0; i < 2; ++i) {
    __m128i a1 = op[0];
    __m128i c1 = op[1];
    __m128i d1 = op[2];
    __m128i b1 = op[3];
    a1 = _mm_add_epi32(a1, c1);          // a1 += c1
    d1 = _mm_sub_epi32(d1, b1);          // d1 -= b1
    __m128i e1 = _mm_sub_epi32(a1, d1);  // e1 = (a1 - d1) >> 1
    e1 = _mm_srai_epi32(e1, 1);
    b1 = _mm_sub_epi32(e1, b1);  // b1 = e1 - b1
    c1 = _mm_sub_epi32(e1, c1);  // c1 = e1 - c1
    a1 = _mm_sub_epi32(a1, b1);  // a1 -= b1
    d1 = _mm_add_epi32(d1, c1);  // d1 += c1

    op[0] = a1;
    op[1] = b1;
    op[2] = c1;
    op[3] = d1;
    if (i == 0) {
      transpose_32bit_4x4(op, op);
    }
  }

  // Convert to int16_t. The C code checks that we are in range.
  op[0] = _mm_packs_epi32(op[0], op[1]);
  op[1] = _mm_packs_epi32(op[2], op[3]);

  // Load uint16_t.
  __m128i dst[2];
  __m128i tmp[4];
  tmp[0] = _mm_loadl_epi64((const __m128i *)(dest + 0 * stride));
  tmp[1] = _mm_loadl_epi64((const __m128i *)(dest + 1 * stride));
  dst[0] = _mm_unpacklo_epi64(tmp[0], tmp[1]);
  tmp[2] = _mm_loadl_epi64((const __m128i *)(dest + 2 * stride));
  tmp[3] = _mm_loadl_epi64((const __m128i *)(dest + 3 * stride));
  dst[1] = _mm_unpacklo_epi64(tmp[2], tmp[3]);

  // Add to the previous results.
  dst[0] = _mm_add_epi16(dst[0], op[0]);
  dst[1] = _mm_add_epi16(dst[1], op[1]);

  // Clamp.
  dst[0] = highbd_clamp_epi16(dst[0], bd);
  dst[1] = highbd_clamp_epi16(dst[1], bd);

  // Store.
  _mm_storel_epi64((__m128i *)(dest + 0 * stride), dst[0]);
  dst[0] = _mm_srli_si128(dst[0], 8);
  _mm_storel_epi64((__m128i *)(dest + 1 * stride), dst[0]);
  _mm_storel_epi64((__m128i *)(dest + 2 * stride), dst[1]);
  dst[1] = _mm_srli_si128(dst[1], 8);
  _mm_storel_epi64((__m128i *)(dest + 3 * stride), dst[1]);
}
