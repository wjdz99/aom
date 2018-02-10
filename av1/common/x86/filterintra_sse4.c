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

  // Load weights of each neighbor at 8 pixels into a vector
  const __m128i f0 = xx_load_128(av1_filter_intra_taps[mode][0]);
  const __m128i f1 = xx_load_128(av1_filter_intra_taps[mode][1]);
  const __m128i f2 = xx_load_128(av1_filter_intra_taps[mode][2]);
  const __m128i f3 = xx_load_128(av1_filter_intra_taps[mode][3]);
  const __m128i f4 = xx_load_128(av1_filter_intra_taps[mode][4]);
  const __m128i f5 = xx_load_128(av1_filter_intra_taps[mode][5]);
  const __m128i f6 = xx_load_128(av1_filter_intra_taps[mode][6]);

  for (r = 1; r < bh + 1; r += 2) {
    for (c = 1; c < bw + 1; c += 4) {
      __m128i in = _mm_set1_epi16((int16_t)buffer[r - 1][c - 1]);
      __m128i out = _mm_mullo_epi16(in, f0);

      in = _mm_set1_epi16((int16_t)buffer[r - 1][c]);
      __m128i tmpout = _mm_mullo_epi16(in, f1);
      out = _mm_add_epi16(out, tmpout);

      in = _mm_set1_epi16((int16_t)buffer[r - 1][c + 1]);
      tmpout = _mm_mullo_epi16(in, f2);
      out = _mm_add_epi16(out, tmpout);

      in = _mm_set1_epi16((int16_t)buffer[r - 1][c + 2]);
      tmpout = _mm_mullo_epi16(in, f3);
      out = _mm_add_epi16(out, tmpout);

      in = _mm_set1_epi16((int16_t)buffer[r - 1][c + 3]);
      tmpout = _mm_mullo_epi16(in, f4);
      out = _mm_add_epi16(out, tmpout);

      in = _mm_set1_epi16((int16_t)buffer[r][c - 1]);
      tmpout = _mm_mullo_epi16(in, f5);
      out = _mm_add_epi16(out, tmpout);

      in = _mm_set1_epi16((int16_t)buffer[r + 1][c - 1]);
      tmpout = _mm_mullo_epi16(in, f6);
      out = _mm_add_epi16(out, tmpout);

      // Rounding
      out = xx_roundn_epi16(out, FILTER_INTRA_SCALE_BITS);
      tmpout = _mm_unpackhi_epi64(out, out);
      out = _mm_unpacklo_epi64(out, out);

      // Clipping
      const __m128i out0123 = _mm_packus_epi16(out, out);
      const __m128i out4567 = _mm_packus_epi16(tmpout, tmpout);

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
