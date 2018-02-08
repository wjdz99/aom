/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <smmintrin.h>
#include "./av1_rtcd.h"
#include "av1/common/enums.h"
#include "av1/common/reconintra.h"
#include "aom_dsp/x86/synonyms.h"

void filter_intra_predictor_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                   TX_SIZE tx_size, const uint8_t *above,
                                   const uint8_t *left, int mode) {
  int r, c;
  uint8_t buffer[33][33];
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];

  assert(bw <= 32 && bh <= 32);

  // The initialization is just for silencing Jenkins static analysis warnings
  for (r = 0; r < bh + 1; ++r)
    memset(buffer[r], 0, (bw + 1) * sizeof(buffer[0][0]));

  for (r = 0; r < bh; ++r) buffer[r + 1][0] = left[r];
  memcpy(buffer[0], &above[-1], (bw + 1) * sizeof(uint8_t));

  const __m128i f0 = xx_load_128(filter_intra_taps[mode][0]);
  const __m128i f1 = xx_load_128(filter_intra_taps[mode][1]);
  const __m128i f2 = xx_load_128(filter_intra_taps[mode][2]);
  const __m128i f3 = xx_load_128(filter_intra_taps[mode][3]);
  const __m128i f4 = xx_load_128(filter_intra_taps[mode][4]);
  const __m128i f5 = xx_load_128(filter_intra_taps[mode][5]);
  const __m128i f6 = xx_load_128(filter_intra_taps[mode][6]);
  const __m128i f7 = xx_load_128(filter_intra_taps[mode][7]);

  for (r = 1; r < bh + 1; r += 2) {
    for (c = 1; c < bw + 1; c += 4) {
      DECLARE_ALIGNED(16, uint8_t, p[8]);
      memcpy(p, &buffer[r - 1][c - 1], 5 * sizeof(uint8_t));
      p[5] = buffer[r][c - 1];
      p[6] = buffer[r + 1][c - 1];
      p[7] = 0;
      const __m128i in = xx_loadl_64(p);
      const __m128i in_w = _mm_cvtepu8_epi16(in);

      __m128i out0_w = _mm_mullo_epi16(in_w, f0);
      const __m128i out1_w = _mm_mullo_epi16(in_w, f1);
      out0_w = _mm_hadd_epi16(out0_w, out1_w);

      __m128i out2_w = _mm_mullo_epi16(in_w, f2);
      const __m128i out3_w = _mm_mullo_epi16(in_w, f3);
      out2_w = _mm_hadd_epi16(out2_w, out3_w);

      out0_w = _mm_hadd_epi16(out0_w, out2_w);

      __m128i out4_w = _mm_mullo_epi16(in_w, f4);
      const __m128i out5_w = _mm_mullo_epi16(in_w, f5);
      out4_w = _mm_hadd_epi16(out4_w, out5_w);

      __m128i out6_w = _mm_mullo_epi16(in_w, f6);
      const __m128i out7_w = _mm_mullo_epi16(in_w, f7);
      out6_w = _mm_hadd_epi16(out6_w, out7_w);

      out4_w = _mm_hadd_epi16(out4_w, out6_w);

      out0_w = _mm_hadd_epi16(out0_w, out4_w);

      // Rounding
      out0_w = xx_roundn_epi16(out0_w, FILTER_INTRA_SCALE_BITS);
      out4_w = _mm_unpackhi_epi64(out0_w, out0_w);
      out0_w = _mm_unpacklo_epi64(out0_w, out0_w);

      // Clipping
      const __m128i out0123 = _mm_packus_epi16(out0_w, out0_w);
      const __m128i out4567 = _mm_packus_epi16(out4_w, out4_w);

      // Storing
      xx_storel_32(&buffer[r][c], out0123);
      xx_storel_32(&buffer[r + 1][c], out4567);
    }
  }

  for (r = 0; r < bh; ++r) {
    memcpy(dst, &buffer[r + 1][1], bw * sizeof(uint8_t));
    dst += stride;
  }
}
