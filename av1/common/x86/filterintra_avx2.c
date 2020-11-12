/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <smmintrin.h>
#include <immintrin.h>

#include "config/av1_rtcd.h"

#include "aom_dsp/x86/synonyms.h"
#include "av1/common/enums.h"
#include "av1/common/reconintra.h"

// Load source from buffer.
#define LOAD_SRC(src, i, r, c)                           \
  {                                                      \
    memcpy(&src[i], &buffer[r][c], 5 * sizeof(uint8_t)); \
    src[i + 5] = buffer[r + 1][c];                       \
    src[i + 6] = buffer[r + 2][c];                       \
    src[i + 7] = 0;                                      \
  }

// Do filter intra for one 4x2 block using sse4_1 instruction set.
#define COMPUTE_FILTER_INTRA_SSE(p, filter_intra_scale)              \
  {                                                                  \
    const __m128i p_b = xx_loadl_64(p);                              \
    const __m128i in = _mm_unpacklo_epi64(p_b, p_b);                 \
    const __m128i out_01 = _mm_maddubs_epi16(in, f1f0);              \
    const __m128i out_23 = _mm_maddubs_epi16(in, f3f2);              \
    const __m128i out_45 = _mm_maddubs_epi16(in, f5f4);              \
    const __m128i out_67 = _mm_maddubs_epi16(in, f7f6);              \
    const __m128i out_0123 = _mm_hadd_epi16(out_01, out_23);         \
    const __m128i out_4567 = _mm_hadd_epi16(out_45, out_67);         \
    const __m128i out_01234567 = _mm_hadd_epi16(out_0123, out_4567); \
    const __m128i round_w = _mm_mulhrs_epi16(                        \
        out_01234567, _mm256_castsi256_si128(filter_intra_scale));   \
    out_r0 = _mm_packus_epi16(round_w, round_w);                     \
    out_r1 = _mm_srli_si128(out_r0, 4);                              \
  }

// Do filter intra for 2 4x2 block's using avx2 instruction set.
#define COMPUTE_FILTER_INTRA_AVX2(p, filter_intra_scale)                \
  {                                                                     \
    const __m128i p_b = _mm_loadu_si128((__m128i *)p);                  \
    const __m256i in =                                                  \
        _mm256_permute4x64_epi64(_mm256_castsi128_si256(p_b), 0x50);    \
    const __m256i out_01 = _mm256_maddubs_epi16(in, f1f0_256b);         \
    const __m256i out_23 = _mm256_maddubs_epi16(in, f3f2_256b);         \
    const __m256i out_45 = _mm256_maddubs_epi16(in, f5f4_256b);         \
    const __m256i out_67 = _mm256_maddubs_epi16(in, f7f6_256b);         \
    const __m256i out_0123 = _mm256_hadd_epi16(out_01, out_23);         \
    const __m256i out_4567 = _mm256_hadd_epi16(out_45, out_67);         \
    const __m256i out_01234567 = _mm256_hadd_epi16(out_0123, out_4567); \
    const __m256i round_w =                                             \
        _mm256_mulhrs_epi16(out_01234567, filter_intra_scale);          \
    out_r0_256b = _mm256_packus_epi16(round_w, round_w);                \
    out_r1_256b = _mm256_srli_si256(out_r0_256b, 4);                    \
  }

// AVX2 implementation for filter_intra_predictor
void av1_filter_intra_predictor_avx2(uint8_t *dst, ptrdiff_t stride,
                                     TX_SIZE tx_size, const uint8_t *above,
                                     const uint8_t *left, int mode) {
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];
  assert(bw <= 32 && bh <= 32);

  if (bw <= 4) {
    av1_filter_intra_predictor_sse4_1(dst, stride, tx_size, above, left, mode);
    return;
  } else {
    int r = 0, c = 0;
    __m128i out_r0, out_r1;
    __m256i out_r0_256b, out_r1_256b;
    uint8_t buffer[33][33];
    // Removed initialization of buffer, not needed

    // Fill temp buffer with top data
    for (r = 0; r < bh; ++r) buffer[r + 1][0] = left[r];
    memcpy(buffer[0], &above[-1], (bw + 1) * sizeof(uint8_t));

    // Load filter weights to 128b register
    const __m128i f1f0 = xx_load_128(av1_filter_intra_taps[mode][0]);
    const __m128i f3f2 = xx_load_128(av1_filter_intra_taps[mode][2]);
    const __m128i f5f4 = xx_load_128(av1_filter_intra_taps[mode][4]);
    const __m128i f7f6 = xx_load_128(av1_filter_intra_taps[mode][6]);
    const __m256i filter_intra_scale_bits_256 =
        _mm256_set1_epi16(1 << (15 - FILTER_INTRA_SCALE_BITS));
    // Filter weights in 256b register
    const __m256i f1f0_256b =
        _mm256_insertf128_si256(_mm256_castsi128_si256(f1f0), (f1f0), 0x1);
    const __m256i f3f2_256b =
        _mm256_insertf128_si256(_mm256_castsi128_si256(f3f2), (f3f2), 0x1);
    const __m256i f5f4_256b =
        _mm256_insertf128_si256(_mm256_castsi128_si256(f5f4), (f5f4), 0x1);
    const __m256i f7f6_256b =
        _mm256_insertf128_si256(_mm256_castsi128_si256(f7f6), (f7f6), 0x01);
    DECLARE_ALIGNED(16, uint8_t, p[16]);

    // The process of filter intra prediction happens at 4x2 block level.
    // Processing of current 4x2 block is dependent on it's left, top 4x2
    // blocks which have to be already processed.
    // So, AVX2 is implemented like mentioned below.
    // 1) Process one 4x2 block of first row.
    // 2) After processing the first 4x2 block, process 2 4x2 block's for
    // which dependency is met (E.g. After 1st 4x2 block, current row next
    // 4x2 block and next rows first 4x2 block can be processed
    // simultaneously).
    // 3) Follow the 2nd process until last 4x2 block remaining.
    // 4) Process the last 4x2 block at the end.

    // Process first 4x2 block
    {
      // Load source data
      LOAD_SRC(p, 0, 0, 0);
      COMPUTE_FILTER_INTRA_SSE(p, filter_intra_scale_bits_256);
      // Store the output to buffer
      xx_storel_32(&buffer[1][1], out_r0);
      xx_storel_32(&buffer[2][1], out_r1);
    }
    // Process 2 4x2 blocks which are diagonal
    for (r = 1; r + 1 < bh; r += 4) {
      for (c = 1; c + 3 < bw; c += 4) {
        // Load source data
        LOAD_SRC(p, 0, r - 1, c + 3);
        LOAD_SRC(p, 8, r + 1, c - 1);
        // Process 2 4x2 blocks
        COMPUTE_FILTER_INTRA_AVX2(p, filter_intra_scale_bits_256);
        // Store the output to buffer
        xx_storel_32(&buffer[r][c + 4], _mm256_castsi256_si128(out_r0_256b));
        xx_storel_32(&buffer[r + 1][c + 4],
                     _mm256_castsi256_si128(out_r1_256b));
        xx_storel_32(&buffer[r + 2][c],
                     _mm256_extracti128_si256(out_r0_256b, 1));
        xx_storel_32(&buffer[r + 3][c],
                     _mm256_extracti128_si256(out_r1_256b, 1));
      }
      // After processing whole width, checks for rows completion. If
      // enough rows are present process current row last 4x2 block with
      // next row starting 4x2 block.
      if (r + 3 < bh) {
        // Load source data
        LOAD_SRC(p, 0, r + 3, 0);
        LOAD_SRC(p, 8, r + 1, c - 1);
        // Process 2 4x2 blocks
        COMPUTE_FILTER_INTRA_AVX2(p, filter_intra_scale_bits_256);
        // Store the output to buffer
        xx_storel_32(&buffer[r + 4][1], _mm256_castsi256_si128(out_r0_256b));
        xx_storel_32(&buffer[r + 5][1], _mm256_castsi256_si128(out_r1_256b));
        xx_storel_32(&buffer[r + 2][c],
                     _mm256_extracti128_si256(out_r0_256b, 1));
        xx_storel_32(&buffer[r + 3][c],
                     _mm256_extracti128_si256(out_r1_256b, 1));
      }
    }
    // Process last 4x2 block at the end
    {
      // Load source data
      LOAD_SRC(p, 0, r - 3, c - 1);
      // Compute inta 4x2 single block
      COMPUTE_FILTER_INTRA_SSE(p, filter_intra_scale_bits_256);
      // Store the output to buffer
      xx_storel_32(&buffer[r - 2][c], out_r0);
      xx_storel_32(&buffer[r - 1][c], out_r1);
    }
    // Copy output to dst from temp buffer
    for (r = 0; r < bh; ++r) {
      memcpy(dst, &buffer[r + 1][1], bw * sizeof(uint8_t));
      dst += stride;
    }
  }
}
