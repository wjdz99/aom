/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <smmintrin.h>
#include <string.h>

#include "config/av1_rtcd.h"

#include "aom_dsp/x86/synonyms.h"
#include "av1/common/enums.h"
#include "av1/common/reconintra.h"

//------------------------------------------------------------------------------
// Load functions.

static inline __m128i Load4(const void *src) {
  // With new compilers such as clang 8.0.0 we can use the new _mm_loadu_si32
  // intrinsic. Both _mm_loadu_si32(src) and the code here are compiled into a
  // movss instruction.
  //
  // Until compiler support of _mm_loadu_si32 is widespread, use of
  // _mm_loadu_si32 is banned.
  int val;
  memcpy(&val, src, sizeof(val));
  return _mm_cvtsi32_si128(val);
}

static inline __m128i LoadLo8(const void *a) {
  return _mm_loadl_epi64((const __m128i *)a);
}

static inline __m128i LoadAligned16(const void *a) {
  assert(((uintptr_t)a & 0xf) == 0);
  return _mm_load_si128((const __m128i *)a);
}

//------------------------------------------------------------------------------
// Store functions.

static inline void Store4(void *dst, const __m128i x) {
  const int val = _mm_cvtsi128_si32(x);
  memcpy(dst, &val, sizeof(val));
}

//------------------------------------------------------------------------------
// Arithmetic utilities.

static inline __m128i RightShiftWithRounding_S16(const __m128i v_val_d,
                                                 int bits) {
  assert(bits <= 16);
  const __m128i v_bias_d = _mm_set1_epi16((int16_t)((1 << bits) >> 1));
  const __m128i v_tmp_d = _mm_add_epi16(v_val_d, v_bias_d);
  return _mm_srai_epi16(v_tmp_d, bits);
}

//------------------------------------------------------------------------------
// FilterIntraPredictor_SSE4_1

// This shuffle mask selects 32-bit blocks in the order 0, 1, 0, 1, which
// duplicates the first 8 bytes of a 128-bit vector into the second 8 bytes.
static const int kDuplicateFirstHalf = 0x44;

// Apply all filter taps to the given 7 packed 16-bit values, keeping the 8th
// at zero to preserve the sum.
static inline void Filter4x2_SSE4_1(uint8_t *dst, const ptrdiff_t stride,
                                    const __m128i *pixels,
                                    const __m128i *taps_0_1,
                                    const __m128i *taps_2_3,
                                    const __m128i *taps_4_5,
                                    const __m128i *taps_6_7) {
  const __m128i mul_0_01 = _mm_maddubs_epi16(*pixels, *taps_0_1);
  const __m128i mul_0_23 = _mm_maddubs_epi16(*pixels, *taps_2_3);
  // |output_half| contains 8 partial sums.
  __m128i output_half = _mm_hadd_epi16(mul_0_01, mul_0_23);
  __m128i output = _mm_hadd_epi16(output_half, output_half);
  const __m128i output_row0 =
      _mm_packus_epi16(RightShiftWithRounding_S16(output, 4),
                       /* arbitrary pack arg */ output);
  Store4(dst, output_row0);
  const __m128i mul_1_01 = _mm_maddubs_epi16(*pixels, *taps_4_5);
  const __m128i mul_1_23 = _mm_maddubs_epi16(*pixels, *taps_6_7);
  output_half = _mm_hadd_epi16(mul_1_01, mul_1_23);
  output = _mm_hadd_epi16(output_half, output_half);
  const __m128i output_row1 =
      _mm_packus_epi16(RightShiftWithRounding_S16(output, 4),
                       /* arbitrary pack arg */ output);
  Store4(dst + stride, output_row1);
}

// 4xH transform sizes are given special treatment because LoadLo8 goes out
// of bounds and every block involves the left column. This implementation
// loads TL from the top row for the first block, so it is not
static inline void Filter4xH(uint8_t *dest, ptrdiff_t stride,
                             const uint8_t *const top_ptr,
                             const uint8_t *const left_ptr, int mode,
                             const int height) {
  const __m128i taps_0_1 = LoadAligned16(av1_filter_intra_taps[mode][0]);
  const __m128i taps_2_3 = LoadAligned16(av1_filter_intra_taps[mode][2]);
  const __m128i taps_4_5 = LoadAligned16(av1_filter_intra_taps[mode][4]);
  const __m128i taps_6_7 = LoadAligned16(av1_filter_intra_taps[mode][6]);
  __m128i top = Load4(top_ptr - 1);
  __m128i pixels = _mm_insert_epi8(top, top_ptr[3], 4);
  __m128i left = (height == 4 ? Load4(left_ptr) : LoadLo8(left_ptr));
  left = _mm_slli_si128(left, 5);

  // Relative pixels: top[-1], top[0], top[1], top[2], top[3], left[0], left[1],
  // left[2], left[3], left[4], left[5], left[6], left[7]
  pixels = _mm_or_si128(left, pixels);

  // Duplicate first 8 bytes.
  pixels = _mm_shuffle_epi32(pixels, kDuplicateFirstHalf);
  Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                   &taps_6_7);
  dest += stride;  // Move to y = 1.
  pixels = Load4(dest);

  // Relative pixels: top[0], top[1], top[2], top[3], empty, left[-2], left[-1],
  // left[0], left[1], ...
  pixels = _mm_or_si128(left, pixels);

  // This mask rearranges bytes in the order: 6, 0, 1, 2, 3, 7, 8, 15. The last
  // byte is an unused value, which shall be multiplied by 0 when we apply the
  // filter.
  const int64_t kInsertTopLeftFirstMask = 0x0F08070302010006;

  // Insert left[-1] in front as TL and put left[0] and left[1] at the end.
  const __m128i pixel_order1 = _mm_set1_epi64x(kInsertTopLeftFirstMask);
  pixels = _mm_shuffle_epi8(pixels, pixel_order1);
  dest += stride;  // Move to y = 2.
  Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                   &taps_6_7);
  dest += stride;  // Move to y = 3.

  // Compute the middle 8 rows before using common code for the final 4 rows.
  // Because the common code below this block assumes that
  if (height == 16) {
    // This shift allows us to use pixel_order2 twice after shifting by 2 later.
    left = _mm_slli_si128(left, 1);
    pixels = Load4(dest);

    // Relative pixels: top[0], top[1], top[2], top[3], empty, empty, left[-4],
    // left[-3], left[-2], left[-1], left[0], left[1], left[2], left[3]
    pixels = _mm_or_si128(left, pixels);

    // This mask rearranges bytes in the order: 9, 0, 1, 2, 3, 7, 8, 15. The
    // last byte is an unused value, as above. The top-left was shifted to
    // position nine to keep two empty spaces after the top pixels.
    const int64_t kInsertTopLeftSecondMask = 0x0F0B0A0302010009;

    // Insert (relative) left[-1] in front as TL and put left[0] and left[1] at
    // the end.
    const __m128i pixel_order2 = _mm_set1_epi64x(kInsertTopLeftSecondMask);
    pixels = _mm_shuffle_epi8(pixels, pixel_order2);
    dest += stride;  // Move to y = 4.

    // First 4x2 in the if body.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);

    // Clear all but final pixel in the first 8 of left column.
    __m128i keep_top_left = _mm_srli_si128(left, 13);
    dest += stride;  // Move to y = 5.
    pixels = Load4(dest);
    left = _mm_srli_si128(left, 2);

    // Relative pixels: top[0], top[1], top[2], top[3], left[-6],
    // left[-5], left[-4], left[-3], left[-2], left[-1], left[0], left[1]
    pixels = _mm_or_si128(left, pixels);
    left = LoadLo8(left_ptr + 8);

    pixels = _mm_shuffle_epi8(pixels, pixel_order2);
    dest += stride;  // Move to y = 6.

    // Second 4x2 in the if body.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);

    // Position TL value so we can use pixel_order1.
    keep_top_left = _mm_slli_si128(keep_top_left, 6);
    dest += stride;  // Move to y = 7.
    pixels = Load4(dest);
    left = _mm_slli_si128(left, 7);
    left = _mm_or_si128(left, keep_top_left);

    // Relative pixels: top[0], top[1], top[2], top[3], empty, empty,
    // left[-1], left[0], left[1], left[2], left[3], ...
    pixels = _mm_or_si128(left, pixels);
    pixels = _mm_shuffle_epi8(pixels, pixel_order1);
    dest += stride;  // Move to y = 8.

    // Third 4x2 in the if body.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);
    dest += stride;  // Move to y = 9.

    // Prepare final inputs.
    pixels = Load4(dest);
    left = _mm_srli_si128(left, 2);

    // Relative pixels: top[0], top[1], top[2], top[3], left[-3], left[-2]
    // left[-1], left[0], left[1], left[2], left[3], ...
    pixels = _mm_or_si128(left, pixels);
    pixels = _mm_shuffle_epi8(pixels, pixel_order1);
    dest += stride;  // Move to y = 10.

    // Fourth 4x2 in the if body.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);
    dest += stride;  // Move to y = 11.
  }

  // In both the 8 and 16 case, we assume that the left vector has the next TL
  // at position 8.
  if (height > 4) {
    // Erase prior left pixels by shifting TL to position 0.
    left = _mm_srli_si128(left, 8);
    left = _mm_slli_si128(left, 6);
    pixels = Load4(dest);

    // Relative pixels: top[0], top[1], top[2], top[3], empty, empty,
    // left[-1], left[0], left[1], left[2], left[3], ...
    pixels = _mm_or_si128(left, pixels);
    pixels = _mm_shuffle_epi8(pixels, pixel_order1);
    dest += stride;  // Move to y = 12 or 4.

    // First of final two 4x2 blocks.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);
    dest += stride;  // Move to y = 13 or 5.
    pixels = Load4(dest);
    left = _mm_srli_si128(left, 2);

    // Relative pixels: top[0], top[1], top[2], top[3], left[-3], left[-2]
    // left[-1], left[0], left[1], left[2], left[3], ...
    pixels = _mm_or_si128(left, pixels);
    pixels = _mm_shuffle_epi8(pixels, pixel_order1);
    dest += stride;  // Move to y = 14 or 6.

    // Last of final two 4x2 blocks.
    Filter4x2_SSE4_1(dest, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);
  }
}

static void FilterIntraPredictor_SSE4_1(void *const dest, ptrdiff_t stride,
                                        const void *const top_row,
                                        const void *const left_column, int mode,
                                        const int width, const int height) {
  const uint8_t *const top_ptr = (const uint8_t *)top_row;
  const uint8_t *const left_ptr = (const uint8_t *)left_column;
  uint8_t *dst = (uint8_t *)dest;
  if (width == 4) {
    Filter4xH(dst, stride, top_ptr, left_ptr, mode, height);
    return;
  }

  // There is one set of 7 taps for each of the 4x2 output pixels.
  const __m128i taps_0_1 = LoadAligned16(av1_filter_intra_taps[mode][0]);
  const __m128i taps_2_3 = LoadAligned16(av1_filter_intra_taps[mode][2]);
  const __m128i taps_4_5 = LoadAligned16(av1_filter_intra_taps[mode][4]);
  const __m128i taps_6_7 = LoadAligned16(av1_filter_intra_taps[mode][6]);

  // This mask rearranges bytes in the order: 0, 1, 2, 3, 4, 8, 9, 15. The 15 at
  // the end is an unused value, which shall be multiplied by 0 when we apply
  // the filter.
  const int64_t kCondenseLeftMask = 0x0F09080403020100;

  // Takes the "left section" and puts it right after p0-p4.
  const __m128i pixel_order1 = _mm_set1_epi64x(kCondenseLeftMask);

  // This mask rearranges bytes in the order: 8, 0, 1, 2, 3, 9, 10, 15. The last
  // byte is unused as above.
  const int64_t kInsertTopLeftMask = 0x0F0A090302010008;

  // Shuffles the "top left" from the left section, to the front. Used when
  // grabbing data from left_column and not top_row.
  const __m128i pixel_order2 = _mm_set1_epi64x(kInsertTopLeftMask);

  // This first pass takes care of the cases where the top left pixel comes from
  // top_row.
  __m128i pixels = LoadLo8(top_ptr - 1);
  __m128i left = _mm_slli_si128(Load4(left_column), 8);
  pixels = _mm_or_si128(pixels, left);

  // Two sets of the same pixels to multiply with two sets of taps.
  pixels = _mm_shuffle_epi8(pixels, pixel_order1);
  Filter4x2_SSE4_1(dst, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                   &taps_6_7);
  left = _mm_srli_si128(left, 1);

  // Load
  pixels = Load4(dst + stride);

  // Because of the above shift, this OR 'invades' the final of the first 8
  // bytes of |pixels|. This is acceptable because the 8th filter tap is always
  // a padded 0.
  pixels = _mm_or_si128(pixels, left);
  pixels = _mm_shuffle_epi8(pixels, pixel_order2);
  const ptrdiff_t stride2 = stride << 1;
  const ptrdiff_t stride4 = stride << 2;
  Filter4x2_SSE4_1(dst + stride2, stride, &pixels, &taps_0_1, &taps_2_3,
                   &taps_4_5, &taps_6_7);
  dst += 4;
  for (int x = 3; x < width - 4; x += 4) {
    pixels = Load4(top_ptr + x);
    pixels = _mm_insert_epi8(pixels, top_ptr[x + 4], 4);
    pixels = _mm_insert_epi8(pixels, dst[-1], 5);
    pixels = _mm_insert_epi8(pixels, dst[stride - 1], 6);

    // Duplicate bottom half into upper half.
    pixels = _mm_shuffle_epi32(pixels, kDuplicateFirstHalf);
    Filter4x2_SSE4_1(dst, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);
    pixels = Load4(dst + stride - 1);
    pixels = _mm_insert_epi8(pixels, dst[stride + 3], 4);
    pixels = _mm_insert_epi8(pixels, dst[stride2 - 1], 5);
    pixels = _mm_insert_epi8(pixels, dst[stride + stride2 - 1], 6);

    // Duplicate bottom half into upper half.
    pixels = _mm_shuffle_epi32(pixels, kDuplicateFirstHalf);
    Filter4x2_SSE4_1(dst + stride2, stride, &pixels, &taps_0_1, &taps_2_3,
                     &taps_4_5, &taps_6_7);
    dst += 4;
  }

  // Now we handle heights that reference previous blocks rather than top_row.
  for (int y = 4; y < height; y += 4) {
    // Leftmost 4x4 block for this height.
    dst -= width;
    dst += stride4;

    // Top Left is not available by offset in these leftmost blocks.
    pixels = Load4(dst - stride);
    left = _mm_slli_si128(Load4(left_ptr + y - 1), 8);
    left = _mm_insert_epi8(left, left_ptr[y + 3], 12);
    pixels = _mm_or_si128(pixels, left);
    pixels = _mm_shuffle_epi8(pixels, pixel_order2);
    Filter4x2_SSE4_1(dst, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                     &taps_6_7);

    // The bytes shifted into positions 6 and 7 will be ignored by the shuffle.
    left = _mm_srli_si128(left, 2);
    pixels = Load4(dst + stride);
    pixels = _mm_or_si128(pixels, left);
    pixels = _mm_shuffle_epi8(pixels, pixel_order2);
    Filter4x2_SSE4_1(dst + stride2, stride, &pixels, &taps_0_1, &taps_2_3,
                     &taps_4_5, &taps_6_7);

    dst += 4;

    // Remaining 4x4 blocks for this height.
    for (int x = 4; x < width; x += 4) {
      pixels = Load4(dst - stride - 1);
      pixels = _mm_insert_epi8(pixels, dst[-stride + 3], 4);
      pixels = _mm_insert_epi8(pixels, dst[-1], 5);
      pixels = _mm_insert_epi8(pixels, dst[stride - 1], 6);

      // Duplicate bottom half into upper half.
      pixels = _mm_shuffle_epi32(pixels, kDuplicateFirstHalf);
      Filter4x2_SSE4_1(dst, stride, &pixels, &taps_0_1, &taps_2_3, &taps_4_5,
                       &taps_6_7);
      pixels = Load4(dst + stride - 1);
      pixels = _mm_insert_epi8(pixels, dst[stride + 3], 4);
      pixels = _mm_insert_epi8(pixels, dst[stride2 - 1], 5);
      pixels = _mm_insert_epi8(pixels, dst[stride2 + stride - 1], 6);

      // Duplicate bottom half into upper half.
      pixels = _mm_shuffle_epi32(pixels, kDuplicateFirstHalf);
      Filter4x2_SSE4_1(dst + stride2, stride, &pixels, &taps_0_1, &taps_2_3,
                       &taps_4_5, &taps_6_7);
      dst += 4;
    }
  }
}

void av1_filter_intra_predictor_sse4_1(uint8_t *dst, ptrdiff_t stride,
                                       TX_SIZE tx_size, const uint8_t *above,
                                       const uint8_t *left, int mode) {
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];
  FilterIntraPredictor_SSE4_1(dst, stride, above, left, mode, bw, bh);
}
