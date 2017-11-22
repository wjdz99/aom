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

#include <smmintrin.h>  // SSE4.1

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/x86/mem_sse4.h"
#include "av1/common/enums.h"
#include "av1/decoder/x86/decodetxb_sse2.h"

static INLINE void dequant_txb_kernel_sse4_1(
    const __m128i level, const __m128i mask,
    const __m128i *const dqvs /*dqvs[2]*/, const int shift,
    tran_low_t **const tcoeffs) {
  const __m128i zero = _mm_setzero_si128();
  __m128i ls[2], ts[2], masks[2];

  masks[0] = _mm_unpacklo_epi16(mask, mask);
  masks[1] = _mm_unpackhi_epi16(mask, mask);
  ts[0] = _mm_load_si128((__m128i *)(*tcoeffs + 0));
  ts[1] = _mm_load_si128((__m128i *)(*tcoeffs + 4));
  ls[0] = _mm_unpacklo_epi16(level, zero);
  ls[1] = _mm_unpackhi_epi16(level, zero);
  ts[0] = _mm_add_epi32(ls[0], ts[0]);
  ts[1] = _mm_add_epi32(ls[1], ts[1]);
  ts[0] = _mm_mullo_epi32(ts[0], dqvs[0]);
  ts[1] = _mm_mullo_epi32(ts[1], dqvs[1]);
  ts[0] = _mm_sra_epi32(ts[0], _mm_cvtsi32_si128(shift));
  ts[1] = _mm_sra_epi32(ts[1], _mm_cvtsi32_si128(shift));
  ts[0] = _mm_xor_si128(ts[0], masks[0]);
  ts[1] = _mm_xor_si128(ts[1], masks[1]);
  ts[0] = _mm_sub_epi32(ts[0], masks[0]);
  ts[1] = _mm_sub_epi32(ts[1], masks[1]);
  _mm_store_si128((__m128i *)(*tcoeffs + 0), ts[0]);
  _mm_store_si128((__m128i *)(*tcoeffs + 4), ts[1]);
  *tcoeffs += 8;
}

static INLINE void dequant_txb_4_sse4_1(const uint8_t *levels,
                                        const int8_t *signs,
                                        const int16_t *const dequant,
                                        const int height, const int shift,
                                        tran_low_t *tcoeffs) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i const1 = _mm_set1_epi8(1);
  int row = height;
  __m128i dqvs[2];

  construct_dqvs_sse2(dequant, dqvs);

  do {
    const __m128i level8 = load_8bit_4x2_to_1_sse4_1(levels, 4 + TX_PAD_HOR);
    const __m128i sign8 = _mm_loadl_epi64((__m128i *)signs);
    const __m128i level16 = _mm_unpacklo_epi8(level8, zero);
    const __m128i mask8 = _mm_cmpeq_epi8(sign8, const1);
    const __m128i mask16 = _mm_unpacklo_epi8(mask8, mask8);

    dequant_txb_kernel_sse4_1(level16, mask16, dqvs, shift, &tcoeffs);
    dqvs[0] = dqvs[1];
    levels += 2 * (4 + TX_PAD_HOR);
    signs += 8;
  } while (row -= 2);
}

static INLINE void dequant_txb_8_sse4_1(const uint8_t *levels,
                                        const int8_t *signs,
                                        const int16_t *const dequant,
                                        const int height, const int shift,
                                        tran_low_t *tcoeffs) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i const1 = _mm_set1_epi8(1);
  int row = height;
  __m128i dqvs[2];

  construct_dqvs_sse2(dequant, dqvs);

  do {
    const __m128i level8 = _mm_loadl_epi64((__m128i *)levels);
    const __m128i sign8 = _mm_loadl_epi64((__m128i *)signs);
    const __m128i level16 = _mm_unpacklo_epi8(level8, zero);
    const __m128i mask8 = _mm_cmpeq_epi8(sign8, const1);
    const __m128i mask16 = _mm_unpacklo_epi8(mask8, mask8);

    dequant_txb_kernel_sse4_1(level16, mask16, dqvs, shift, &tcoeffs);
    dqvs[0] = dqvs[1];
    levels += 8 + TX_PAD_HOR;
    signs += 8;
  } while (--row);
}

static INLINE void dequant_txb_16n_sse4_1(
    const uint8_t *levels, const int8_t *signs, const int16_t *const dequant,
    const int width, const int height, const int shift, tran_low_t *tcoeffs) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i const1 = _mm_set1_epi8(1);
  int row = height;
  __m128i dqvs[2], level16[2], mask16[2];

  construct_dqvs_sse2(dequant, dqvs);

  do {
    int col = width;
    do {
      const __m128i level8 = _mm_loadu_si128((__m128i *)levels);
      const __m128i sign8 = _mm_load_si128((__m128i *)signs);
      const __m128i mask8 = _mm_cmpeq_epi8(sign8, const1);

      level16[0] = _mm_unpacklo_epi8(level8, zero);
      level16[1] = _mm_unpackhi_epi8(level8, zero);
      mask16[0] = _mm_unpacklo_epi8(mask8, mask8);
      mask16[1] = _mm_unpackhi_epi8(mask8, mask8);

      dequant_txb_kernel_sse4_1(level16[0], mask16[0], dqvs, shift, &tcoeffs);
      dqvs[0] = dqvs[1];
      dequant_txb_kernel_sse4_1(level16[1], mask16[1], dqvs, shift, &tcoeffs);
      levels += 16;
      signs += 16;
    } while (col -= 16);
    levels += TX_PAD_HOR;
  } while (--row);
}

void av1_dequant_txb_sse4_1(const uint8_t *const levels,
                            const int8_t *const signs,
                            const int16_t *const dequant,
                            const int16_t *const scan, const int bwl,
                            const int height, const int eob, const int shift,
                            tran_low_t *const tcoeffs) {
  const int width = 1 << bwl;
  (void)scan;
  (void)eob;

  if (width == 4) {
    dequant_txb_4_sse4_1(levels, signs, dequant, height, shift, tcoeffs);
  } else if (width == 8) {
    dequant_txb_8_sse4_1(levels, signs, dequant, height, shift, tcoeffs);
  } else {
    dequant_txb_16n_sse4_1(levels, signs, dequant, width, height, shift,
                           tcoeffs);
  }
}
