/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <immintrin.h>

#include "config/av1_rtcd.h"

#include "av1/common/resize.h"

#include "aom_dsp/x86/synonyms.h"

#define PROCESS_RESIZE_Y_WD8                                           \
  /* ah0 ah1 ... ah7 */                                                \
  const __m128i AH = _mm_add_epi16(l0, l7);                            \
  /* bg0 bg1 ... bh7 */                                                \
  const __m128i BG = _mm_add_epi16(l1, l6);                            \
  /* cf0 cf1 ... cf7 */                                                \
  const __m128i CF = _mm_add_epi16(l2, l5);                            \
  /* de0 de1 ... de7 */                                                \
  const __m128i DE = _mm_add_epi16(l3, l4);                            \
                                                                       \
  /* ah0 bg0 ... ah3 bg3 */                                            \
  const __m128i AHBG_low = _mm_unpacklo_epi16(AH, BG);                 \
  /*cf0 de0 ... cf2 de2 */                                             \
  const __m128i CFDE_low = _mm_unpacklo_epi16(CF, DE);                 \
                                                                       \
  /* ah4 bg4... ah7 bg7 */                                             \
  const __m128i AHBG_hi = _mm_unpackhi_epi16(AH, BG);                  \
  /* cf4 de4... cf7 de7 */                                             \
  const __m128i CFDE_hi = _mm_unpackhi_epi16(CF, DE);                  \
                                                                       \
  /* r00 r01 r02 r03 */                                                \
  const __m128i r00 = _mm_madd_epi16(AHBG_low, coeffs_y[0]);           \
  const __m128i r01 = _mm_madd_epi16(CFDE_low, coeffs_y[1]);           \
  __m128i r0 = _mm_add_epi32(r00, r01);                                \
  /* r04 r05 r06 r07 */                                                \
  const __m128i r10 = _mm_madd_epi16(AHBG_hi, coeffs_y[0]);            \
  const __m128i r11 = _mm_madd_epi16(CFDE_hi, coeffs_y[1]);            \
  __m128i r1 = _mm_add_epi32(r10, r11);                                \
                                                                       \
  r0 = _mm_add_epi32(r0, round_const_bits);                            \
  r1 = _mm_add_epi32(r1, round_const_bits);                            \
  r0 = _mm_sra_epi32(r0, round_shift_bits);                            \
  r1 = _mm_sra_epi32(r1, round_shift_bits);                            \
                                                                       \
  /* r00 ... r07 (8 values of each 16bit) */                           \
  const __m128i res_16b = _mm_packs_epi32(r0, r1);                     \
  /* r00 ... r07 | r00 ... r07 (16 values of each 8bit) */             \
  const __m128i res_8b0 = _mm_packus_epi16(res_16b, res_16b);          \
                                                                       \
  __m128i res = _mm_min_epu8(res_8b0, clip_pixel);                     \
  res = _mm_max_epu8(res, zero);                                       \
  _mm_storel_epi64((__m128i *)&output[(i / 2) * out_stride + j], res); \
                                                                       \
  l0 = l2;                                                             \
  l1 = l3;                                                             \
  l2 = l4;                                                             \
  l3 = l5;                                                             \
  l4 = l6;                                                             \
  l5 = l7;                                                             \
  data += 2 * stride;

static INLINE void prepare_filter_coeffs(const int16_t *filter,
                                         __m128i *const coeffs /* [2] */) {
  // f0 f1 f2 f3 x x x x
  const __m128i sym_even_filter = _mm_loadl_epi64((__m128i *)filter);

  // f1 f0 f3 f2 x x x x
  const __m128i tmp1 = _mm_shufflelo_epi16(sym_even_filter, 0xb1);

  // f3 f2 f3 f2 ...
  coeffs[0] = _mm_shuffle_epi32(tmp1, 0x55);
  // f1 f0 f1 f0 ...
  coeffs[1] = _mm_shuffle_epi32(tmp1, 0x00);
}

bool av1_resize_vert_dir_sse2(uint8_t *intbuf, uint8_t *output, int out_stride,
                              int height, int height2, int stride,
                              int start_col) {
  // For the GM tool, the input layer height or width is assured to be an even
  // number. Hence the function 'down2_symodd()' is not invoked and SIMD
  // optimization of the same is not implemented.
  // When the input height is less than 8 and even, the potential input
  // heights are limited to 2, 4, or 6. These scenarios require seperate
  // handling due to padding requirements. Invoking the C function here will
  // eliminate the need for conditional statements within the subsequent SIMD
  // code to manage these cases.
  if (height & 1 || height < 8) {
    return av1_resize_vert_dir_c(intbuf, output, out_stride, height, height2,
                                 stride, start_col);
  }

  __m128i coeffs_y[2];
  const int bits = FILTER_BITS;
  const __m128i round_const_bits = _mm_set1_epi32((1 << bits) >> 1);
  const __m128i round_shift_bits = _mm_cvtsi32_si128(bits);
  const uint8_t max_pixel = 255;
  const __m128i clip_pixel = _mm_set1_epi8(max_pixel);
  const __m128i zero = _mm_setzero_si128();
  prepare_filter_coeffs(av1_down2_symeven_half_filter, coeffs_y);

  const int remain_col = stride % 8;

  for (int j = start_col; j < stride - remain_col; j += 8) {
    uint8_t *data = &intbuf[j];
    // d0 ... d7
    const __m128i l8_3 = _mm_loadl_epi64((__m128i *)(data + 0 * stride));
    // Padding top 3 rows with the last available row at the top.
    // a0 ... a7
    const __m128i l8_0 = l8_3;
    // b0 ... b7
    const __m128i l8_1 = l8_3;
    // c0 ... c7
    const __m128i l8_2 = l8_3;
    // e0 ... e7
    const __m128i l8_4 = _mm_loadl_epi64((__m128i *)(data + 1 * stride));
    // f0 ... f7
    const __m128i l8_5 = _mm_loadl_epi64((__m128i *)(data + 2 * stride));

    // Convert to 16bit as addition of 2 source pixel crosses 8 bit.
    __m128i l0 = _mm_unpacklo_epi8(l8_0, zero);  // A(128bit) = a0 - a7(16 bit)
    __m128i l1 = _mm_unpacklo_epi8(l8_1, zero);  // B(128bit) = b0 - b7(16 bit)
    __m128i l2 = _mm_unpacklo_epi8(l8_2, zero);  // C(128bit) = c0 - c7(16 bit)
    __m128i l3 = _mm_unpacklo_epi8(l8_3, zero);  // D(128bit) = d0 - d7(16 bit)
    __m128i l4 = _mm_unpacklo_epi8(l8_4, zero);  // E(128bit) = e0 - e7(16 bit)
    __m128i l5 = _mm_unpacklo_epi8(l8_5, zero);  // F(128bit) = f0 - f7(16 bit)

    // Increment the pointer such that the loading starts from row G.
    data = data + 3 * stride;
    // The core vertical SIMD processes 2 input rows simultaneously to generate
    // output corresponding to 1 row. To streamline the core loop and eliminate
    // the need for conditional checks, the remaining rows 4 are processed
    // separately.
    for (int i = 0; i < height - 4; i += 2) {
      // g0 ... g7
      __m128i l8_6 = _mm_loadl_epi64((__m128i *)(data));
      // h0 ... h7
      __m128i l8_7 = _mm_loadl_epi64((__m128i *)(data + stride));
      __m128i l6 = _mm_unpacklo_epi8(l8_6, zero);  // G(128bit):g0-g7(16b)
      __m128i l7 = _mm_unpacklo_epi8(l8_7, zero);  // H(128bit):h0-h7(16b)

      PROCESS_RESIZE_Y_WD8
    }

    __m128i l8_6 = _mm_loadl_epi64((__m128i *)(data));
    __m128i l6 = _mm_unpacklo_epi8(l8_6, zero);
    // Process the last 4 input rows here.
    for (int i = height - 4; i < height; i += 2) {
      __m128i l7 = l6;
      PROCESS_RESIZE_Y_WD8
    }
  }

  if (remain_col)
    return av1_resize_vert_dir_c(intbuf, output, out_stride, height, height2,
                                 stride, stride - remain_col);

  return true;
}

// Blends a and b using mask and returns the result.
__m128i blend(__m128i a, __m128i b, __m128i mask) {
  const __m128i masked_b = _mm_and_si128(mask, b);
  const __m128i masked_a = _mm_andnot_si128(mask, a);
  return (_mm_or_si128(masked_a, masked_b));
}

// Masks used for width 32 and 8 pixels, with left and right padding
// requirements.

static const uint8_t wd8_left_padding_mask[16] = {
  255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

static const uint8_t wd8_right_padding_mask[16] = { 0,   0,   0,   0,  0,   0,
                                                    0,   0,   0,   0,  255, 255,
                                                    255, 255, 255, 255 };

static const uint8_t mask_even[16] = {
  255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0,
};

void av1_resize_horz_dir_sse2(const uint8_t *const input, int in_stride,
                              uint8_t *intbuf, int height, int filtered_length,
                              int width2) {
  assert(height % 2 == 0);
  if (filtered_length < 8) {
    av1_resize_horz_dir_c(input, in_stride, intbuf, height, filtered_length,
                          width2);
  }

  __m128i coeffs_x[2];
  const int bits = FILTER_BITS;
  const int dst_stride = width2;
  const int remain_col = filtered_length % 8;
  const __m128i round_const_bits = _mm_set1_epi32((1 << bits) >> 1);
  const __m128i round_shift_bits = _mm_cvtsi32_si128(bits);

  const uint8_t max_pixel = 255;
  const __m128i clip_pixel = _mm_set1_epi8((char)max_pixel);
  const __m128i zero = _mm_setzero_si128();

  const __m128i wd8_start_pad_mask =
      _mm_loadu_si128((__m128i *)wd8_left_padding_mask);
  const __m128i wd8_end_pad_mask =
      _mm_loadu_si128((__m128i *)wd8_right_padding_mask);
  const __m128i wd_8_mask_even = _mm_loadu_si128((__m128i *)mask_even);
  prepare_filter_coeffs(av1_down2_symeven_half_filter, coeffs_x);

  for (int i = 0; i < height; i += 2) {
    int filter_offset = 0;
    for (int j = 0; j + 8 <= filtered_length; j += 8) {
      const int in_idx = i * in_stride + j - filter_offset;
      const int out_idx = i * dst_stride + j / 2;

      // a0 a1 a2 .... a15
      __m128i row0 = _mm_loadu_si128((__m128i *)&input[in_idx]);
      // b0 b1 b2 .... b15
      __m128i row1 = _mm_loadu_si128((__m128i *)&input[in_idx + in_stride]);
      filter_offset = 3;

      const __m128i start_pixel_row0 = _mm_set1_epi8(input[i * in_stride]);
      const __m128i start_pixel_row1 =
          _mm_set1_epi8(input[(i + 1) * in_stride]);
      const __m128i end_pixel_row0 =
          _mm_set1_epi8(input[i * in_stride + filtered_length - 1]);
      const __m128i end_pixel_row1 =
          _mm_set1_epi8(input[(i + 1) * in_stride + filtered_length - 1]);

      // Pad start pixels to the left, while processing the first pixels in the
      // row.
      if (j == 0) {
        row0 = blend(_mm_slli_si128(row0, 3), start_pixel_row0,
                     wd8_start_pad_mask);
        row1 = blend(_mm_slli_si128(row1, 3), start_pixel_row1,
                     wd8_start_pad_mask);
      }

      // Pad end pixels to the right, while processing the last pixels in the
      // row.
      const int is_last_cols8 = (j + 8 == filtered_length);
      if (is_last_cols8) {
        row0 = blend(row0, end_pixel_row0, wd8_end_pad_mask);
        row1 = blend(row1, end_pixel_row1, wd8_end_pad_mask);
      }

      // a0 a2 a4 a6 a8 a10 a12 a14
      const __m128i row0_even = _mm_and_si128(row0, wd_8_mask_even);
      // b0 b2 b4 b6 b8 b10 b12 b14
      const __m128i row1_even = _mm_and_si128(row1, wd_8_mask_even);
      // a1 a3 a5 a7 a9 a11 a13 a15
      const __m128i row0_odd =
          _mm_and_si128(_mm_srli_epi16(row0, 8), wd_8_mask_even);
      // b1 b3 b5 b7 b9 b11 b13 b15
      const __m128i row1_odd =
          _mm_and_si128(_mm_srli_epi16(row1, 8), wd_8_mask_even);

      // a0 b0 a2 b2 a4 b4 a6 b6 (each 16 bit)
      const __m128i s0 = _mm_unpacklo_epi16(row0_even, row1_even);
      // a1 b1 a3 b3 a5 b5 a7 b7
      const __m128i s1 = _mm_unpacklo_epi16(row0_odd, row1_odd);
      // a2 b2 a4 b4 a6 b6 a8 b8
      const __m128i s2 = _mm_unpacklo_epi16(_mm_srli_si128(row0_even, 2),
                                            _mm_srli_si128(row1_even, 2));
      // a3 b3 a5 b5 a7 b7 a9 b9
      const __m128i s3 = _mm_unpacklo_epi16(_mm_srli_si128(row0_odd, 2),
                                            _mm_srli_si128(row1_odd, 2));
      // a4 b4 a6 b6 a8 b8 a10 b10
      const __m128i s4 = _mm_unpacklo_epi16(_mm_srli_si128(row0_even, 4),
                                            _mm_srli_si128(row1_even, 4));
      // a5 b5 a7 b7 a9 b9 a11 b11
      const __m128i s5 = _mm_unpacklo_epi16(_mm_srli_si128(row0_odd, 4),
                                            _mm_srli_si128(row1_odd, 4));
      // a6 b6 a8 b8 a10 b10 a12 b12
      const __m128i s6 = _mm_unpacklo_epi16(_mm_srli_si128(row0_even, 6),
                                            _mm_srli_si128(row1_even, 6));
      // a7 b7 a9 b9 a11 b11 a13 b13
      const __m128i s7 = _mm_unpacklo_epi16(_mm_srli_si128(row0_odd, 6),
                                            _mm_srli_si128(row1_odd, 6));

      // a0a7 b0b7 .... a6a13 b6b13
      const __m128i s07 = _mm_add_epi16(s0, s7);
      // a1a6 b1b6  .... a7a12 b7b12
      const __m128i s16 = _mm_add_epi16(s1, s6);
      // a2a5 b2b5  .... a8a11 b8b11
      const __m128i s25 = _mm_add_epi16(s2, s5);
      // a3a4 b3b4  .... a9a10 b9b10
      const __m128i s34 = _mm_add_epi16(s3, s4);

      // a0a7 a1a6 b0b7 b1b6 a2a9 a3a8 b2b9 b3b8
      const __m128i s1607_low = _mm_unpacklo_epi16(s07, s16);
      // a2a5 a3a4 b2b5 b3b4 a4a7 a5a6 b4b7 b5b6
      const __m128i s3425_low = _mm_unpacklo_epi16(s25, s34);

      // a4a11 a5a10 b4b11 b5b10 a6a13 a7a12 b6b13 b7b12
      const __m128i s1607_high = _mm_unpackhi_epi16(s07, s16);
      // a6a9 a7a8 b6b9 b7b8 a8a11 a9a10 b8b11 b9b10
      const __m128i s3425_high = _mm_unpackhi_epi16(s25, s34);

      const __m128i r01_0 = _mm_madd_epi16(s3425_low, coeffs_x[1]);
      const __m128i r01_1 = _mm_madd_epi16(s1607_low, coeffs_x[0]);
      const __m128i r01_2 = _mm_madd_epi16(s3425_high, coeffs_x[1]);
      const __m128i r01_3 = _mm_madd_epi16(s1607_high, coeffs_x[0]);

      // Result of first 4 pixels of row0 and row1 (a0 to a3 , b0 to b3).
      // r0_0 r1_0 r0_1 r1_1
      __m128i r0 = _mm_add_epi32(r01_0, r01_1);
      r0 = _mm_add_epi32(r0, round_const_bits);
      r0 = _mm_sra_epi32(r0, round_shift_bits);

      // Result of next 4 pixels of row0 and row1 (a4 to a7 , b4 to b7).
      // r0_2 r1_2 r0_3 r1_3
      __m128i r1 = _mm_add_epi32(r01_2, r01_3);
      r1 = _mm_add_epi32(r1, round_const_bits);
      r1 = _mm_sra_epi32(r1, round_shift_bits);

      // Shuffle the pixels in order and convert the result from 32bit to 8bit.
      // r0_0 r0_1 r1_0 r1_1
      r0 = _mm_shuffle_epi32(r0, 0xd8);
      // r0_2 r0_3 r1_2 r1_3
      r1 = _mm_shuffle_epi32(r1, 0xd8);
      // r0_0 r0_1 r1_0 r1_1 r0_2 r0_3 r1_2 r1_3
      __m128i res_16 = _mm_packs_epi32(r0, r1);
      // r0_0 r0_1 r0_2 r0_3 r1_0 r1_1 r1_2 r1_3
      res_16 = _mm_shuffle_epi32(res_16, 0xd8);
      const __m128i res_8 = _mm_packus_epi16(res_16, res_16);
      __m128i res = _mm_min_epu8(res_8, clip_pixel);
      res = _mm_max_epu8(res, zero);
      const __m128i res_row1 = _mm_shuffle_epi32(res, 0x11);

      // r0_0 r0_1 r0_2 r0_3 r1_0 r1_1 r1_2 r1_3
      xx_storel_32(intbuf + out_idx, res);
      // r1_0 r1_1 r1_2 r1_3 r0_0 r0_1 r0_2 r0_3
      xx_storel_32(intbuf + out_idx + dst_stride, res_row1);
    }

    int wd_processed = filtered_length - remain_col;
    // When the remaining width is 2, the above code would not have taken
    // care of padding required for (filtered_length - 4)th pixel. Hence,
    // process that pixel again with the C code.
    wd_processed = (remain_col == 2) ? wd_processed - 2 : wd_processed;
    if (remain_col) {
      const int in_idx = (in_stride * i);
      const int out_idx = (wd_processed / 2) + width2 * i;

      down2_symeven(input + in_idx, filtered_length, intbuf + out_idx,
                    wd_processed);
      down2_symeven(input + in_idx + in_stride, filtered_length,
                    intbuf + out_idx + width2, wd_processed);
    }
  }
}
