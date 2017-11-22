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

#include <emmintrin.h>  // SSE2

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/x86/mem_sse2.h"
#include "av1/common/enums.h"
#include "av1/decoder/x86/decodetxb_sse2.h"

static INLINE void dequant_txb_kernel_sse2(
    const __m128i level, const __m128i mask,
    const __m128i *const dqvs /*dqvs[2]*/, const int shift,
    tran_low_t **const tcoeffs) {
  const __m128i zero = _mm_setzero_si128();
  // TODO(linfengz): Unit test is not ready to test realistic inputs. So there
  // is temporary right shifting code to handle faked saturations in unit test
  // by now.
  const __m128i left_shift = _mm_cvtsi32_si128(16 - shift);
  const __m128i tcoeff = _mm_load_si128((__m128i *)*tcoeffs);
  __m128i t16, t32[2];

  t16 = _mm_add_epi16(level, tcoeff);
  t32[0] = _mm_unpacklo_epi16(t16, zero);
  t32[1] = _mm_unpackhi_epi16(t16, zero);
  t32[0] = _mm_madd_epi16(t32[0], dqvs[0]);
  t32[1] = _mm_madd_epi16(t32[1], dqvs[1]);
  t32[0] = _mm_sll_epi32(t32[0], left_shift);  // temporary code
  t32[1] = _mm_sll_epi32(t32[1], left_shift);  // temporary code
  t32[0] = _mm_srai_epi32(t32[0], 16);
  t32[1] = _mm_srai_epi32(t32[1], 16);
  t16 = _mm_packs_epi32(t32[0], t32[1]);
  t16 = _mm_xor_si128(t16, mask);
  t16 = _mm_sub_epi16(t16, mask);
  _mm_store_si128((__m128i *)*tcoeffs, t16);
  *tcoeffs += 8;
}

static INLINE void dequant_txb_4_sse2(const uint8_t *levels,
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
    const __m128i level8 = load_8bit_4x2_to_1_sse2(levels, 4 + TX_PAD_HOR);
    const __m128i sign8 = _mm_loadl_epi64((__m128i *)signs);
    const __m128i level16 = _mm_unpacklo_epi8(level8, zero);
    const __m128i mask8 = _mm_cmpeq_epi8(sign8, const1);
    const __m128i mask16 = _mm_unpacklo_epi8(mask8, mask8);

    dequant_txb_kernel_sse2(level16, mask16, dqvs, shift, &tcoeffs);
    dqvs[0] = dqvs[1];
    levels += 2 * (4 + TX_PAD_HOR);
    signs += 8;
  } while (row -= 2);
}

static INLINE void dequant_txb_8_sse2(const uint8_t *levels,
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

    dequant_txb_kernel_sse2(level16, mask16, dqvs, shift, &tcoeffs);
    dqvs[0] = dqvs[1];
    levels += 8 + TX_PAD_HOR;
    signs += 8;
  } while (--row);
}

static INLINE void dequant_txb_16n_sse2(const uint8_t *levels,
                                        const int8_t *signs,
                                        const int16_t *const dequant,
                                        const int width, const int height,
                                        const int shift, tran_low_t *tcoeffs) {
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

      dequant_txb_kernel_sse2(level16[0], mask16[0], dqvs, shift, &tcoeffs);
      dqvs[0] = dqvs[1];
      dequant_txb_kernel_sse2(level16[1], mask16[1], dqvs, shift, &tcoeffs);
      levels += 16;
      signs += 16;
    } while (col -= 16);
    levels += TX_PAD_HOR;
  } while (--row);
}

void av1_dequant_txb_sse2(const uint8_t *const levels,
                          const int8_t *const signs,
                          const int16_t *const dequant,
                          const int16_t *const scan, const int bwl,
                          const int height, const int eob, const int shift,
                          tran_low_t *const tcoeffs) {
  const int width = 1 << bwl;
  (void)scan;
  (void)eob;

  if (width == 4) {
    dequant_txb_4_sse2(levels, signs, dequant, height, shift, tcoeffs);
  } else if (width == 8) {
    dequant_txb_8_sse2(levels, signs, dequant, height, shift, tcoeffs);
  } else {
    dequant_txb_16n_sse2(levels, signs, dequant, width, height, shift, tcoeffs);
  }
}
