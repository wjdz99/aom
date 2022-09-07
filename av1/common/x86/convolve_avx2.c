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

#include <immintrin.h>

#include "config/av1_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/x86/convolve_avx2.h"
#include "aom_dsp/x86/convolve_common_intrin.h"
#include "aom_dsp/x86/synonyms.h"

static AOM_INLINE void av1_convolve_y_sr_avx2_general(
    const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn) {
  // right shift is F-1 because we are already dividing
  // filter co-efficients by 2
  const int right_shift_bits = (FILTER_BITS - 1);
  __m128i right_shift = _mm_cvtsi32_si128(right_shift_bits);
  __m256i right_shift_const = _mm256_set1_epi16((1 << right_shift_bits) >> 1);

  __m256i coeffs[6], s[12];
  __m128i d[10];

  int i, j, vert_tap = get_filter_tap(filter_params_y, subpel_y_qn);

  if (vert_tap == 6)
    prepare_coeffs_6t_lowbd(filter_params_y, subpel_y_qn, coeffs);
  else if (vert_tap == 12) {
    prepare_coeffs_12taps(filter_params_y, subpel_y_qn, coeffs);
  } else {
    prepare_coeffs_lowbd(filter_params_y, subpel_y_qn, coeffs);
  }

  // vert_filt as 4 tap
  if (vert_tap == 4) {
    const int fo_vert = 1;
    const uint8_t *const src_ptr = src - fo_vert * src_stride;
    for (j = 0; j < w; j += 16) {
      const uint8_t *data = &src_ptr[j];
      d[0] = _mm_loadu_si128((__m128i *)(data + 0 * src_stride));
      d[1] = _mm_loadu_si128((__m128i *)(data + 1 * src_stride));
      d[2] = _mm_loadu_si128((__m128i *)(data + 2 * src_stride));
      d[3] = _mm_loadu_si128((__m128i *)(data + 3 * src_stride));
      d[4] = _mm_loadu_si128((__m128i *)(data + 4 * src_stride));

      // Load lines a and b. Line a to lower 128, line b to upper 128
      const __m256i src_01a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[0]), _mm256_castsi128_si256(d[1]), 0x20);

      const __m256i src_12a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[1]), _mm256_castsi128_si256(d[2]), 0x20);

      const __m256i src_23a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[2]), _mm256_castsi128_si256(d[3]), 0x20);

      const __m256i src_34a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[3]), _mm256_castsi128_si256(d[4]), 0x20);

      s[0] = _mm256_unpacklo_epi8(src_01a, src_12a);
      s[1] = _mm256_unpacklo_epi8(src_23a, src_34a);

      s[3] = _mm256_unpackhi_epi8(src_01a, src_12a);
      s[4] = _mm256_unpackhi_epi8(src_23a, src_34a);

      for (i = 0; i < h; i += 2) {
        data = &src_ptr[i * src_stride + j];
        d[5] = _mm_loadu_si128((__m128i *)(data + 5 * src_stride));
        const __m256i src_45a = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(d[4]), _mm256_castsi128_si256(d[5]), 0x20);

        d[4] = _mm_loadu_si128((__m128i *)(data + 6 * src_stride));
        const __m256i src_56a = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(d[5]), _mm256_castsi128_si256(d[4]), 0x20);

        s[2] = _mm256_unpacklo_epi8(src_45a, src_56a);
        s[5] = _mm256_unpackhi_epi8(src_45a, src_56a);

        const __m256i res_lo = convolve_lowbd_4tap(s, coeffs + 1);
        /* rounding code */
        // shift by F - 1
        const __m256i res_16b_lo = _mm256_sra_epi16(
            _mm256_add_epi16(res_lo, right_shift_const), right_shift);
        // 8 bit conversion and saturation to uint8
        __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);

        if (w - j > 8) {
          const __m256i res_hi = convolve_lowbd_4tap(s + 3, coeffs + 1);

          /* rounding code */
          // shift by F - 1
          const __m256i res_16b_hi = _mm256_sra_epi16(
              _mm256_add_epi16(res_hi, right_shift_const), right_shift);
          // 8 bit conversion and saturation to uint8
          __m256i res_8b_hi = _mm256_packus_epi16(res_16b_hi, res_16b_hi);

          __m256i res_a = _mm256_unpacklo_epi64(res_8b_lo, res_8b_hi);

          const __m128i res_0 = _mm256_castsi256_si128(res_a);
          const __m128i res_1 = _mm256_extracti128_si256(res_a, 1);

          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res_0);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j + dst_stride],
                           res_1);
        } else {
          const __m128i res_0 = _mm256_castsi256_si128(res_8b_lo);
          const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
          if (w - j > 4) {
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j], res_0);
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j + dst_stride],
                             res_1);
          } else if (w - j > 2) {
            xx_storel_32(&dst[i * dst_stride + j], res_0);
            xx_storel_32(&dst[i * dst_stride + j + dst_stride], res_1);
          } else {
            __m128i *const p_0 = (__m128i *)&dst[i * dst_stride + j];
            __m128i *const p_1 =
                (__m128i *)&dst[i * dst_stride + j + dst_stride];
            *(uint16_t *)p_0 = (uint16_t)_mm_cvtsi128_si32(res_0);
            *(uint16_t *)p_1 = (uint16_t)_mm_cvtsi128_si32(res_1);
          }
        }
        s[0] = s[1];
        s[1] = s[2];

        s[3] = s[4];
        s[4] = s[5];
      }
    }
  } else if (vert_tap == 6) {
    const int fo_vert = vert_tap / 2 - 1;
    const uint8_t *const src_ptr = src - fo_vert * src_stride;

    for (j = 0; j < w; j += 16) {
      const uint8_t *data = &src_ptr[j];
      __m256i src6;

      d[0] = _mm_loadu_si128((__m128i *)(data + 0 * src_stride));
      d[1] = _mm_loadu_si128((__m128i *)(data + 1 * src_stride));
      d[2] = _mm_loadu_si128((__m128i *)(data + 2 * src_stride));
      d[3] = _mm_loadu_si128((__m128i *)(data + 3 * src_stride));
      // Load lines a and b. Line a to lower 128, line b to upper 128
      const __m256i src_01a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[0]), _mm256_castsi128_si256(d[1]), 0x20);

      const __m256i src_12a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[1]), _mm256_castsi128_si256(d[2]), 0x20);

      const __m256i src_23a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[2]), _mm256_castsi128_si256(d[3]), 0x20);

      src6 = _mm256_castsi128_si256(
          _mm_loadu_si128((__m128i *)(data + 4 * src_stride)));
      const __m256i src_34a =
          _mm256_permute2x128_si256(_mm256_castsi128_si256(d[3]), src6, 0x20);

      s[0] = _mm256_unpacklo_epi8(src_01a, src_12a);
      s[1] = _mm256_unpacklo_epi8(src_23a, src_34a);

      s[3] = _mm256_unpackhi_epi8(src_01a, src_12a);
      s[4] = _mm256_unpackhi_epi8(src_23a, src_34a);

      for (i = 0; i < h; i += 2) {
        data = &src_ptr[i * src_stride + j];
        const __m256i src_45a = _mm256_permute2x128_si256(
            src6,
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(data + 5 * src_stride))),
            0x20);

        src6 = _mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(data + 6 * src_stride)));
        const __m256i src_56a = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(data + 5 * src_stride))),
            src6, 0x20);

        s[2] = _mm256_unpacklo_epi8(src_45a, src_56a);
        s[5] = _mm256_unpackhi_epi8(src_45a, src_56a);

        const __m256i res_lo = convolve_lowbd_6tap(s, coeffs);

        /* rounding code */
        // shift by F - 1
        const __m256i res_16b_lo = _mm256_sra_epi16(
            _mm256_add_epi16(res_lo, right_shift_const), right_shift);
        // 8 bit conversion and saturation to uint8
        __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);

        if (w - j > 8) {
          const __m256i res_hi = convolve_lowbd_6tap(s + 3, coeffs);

          /* rounding code */
          // shift by F - 1
          const __m256i res_16b_hi = _mm256_sra_epi16(
              _mm256_add_epi16(res_hi, right_shift_const), right_shift);
          // 8 bit conversion and saturation to uint8
          __m256i res_8b_hi = _mm256_packus_epi16(res_16b_hi, res_16b_hi);

          __m256i res_a = _mm256_unpacklo_epi64(res_8b_lo, res_8b_hi);

          const __m128i res_0 = _mm256_castsi256_si128(res_a);
          const __m128i res_1 = _mm256_extracti128_si256(res_a, 1);

          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res_0);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j + dst_stride],
                           res_1);
        } else {
          const __m128i res_0 = _mm256_castsi256_si128(res_8b_lo);
          const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
          if (w - j > 4) {
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j], res_0);
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j + dst_stride],
                             res_1);
          } else if (w - j > 2) {
            xx_storel_32(&dst[i * dst_stride + j], res_0);
            xx_storel_32(&dst[i * dst_stride + j + dst_stride], res_1);
          } else {
            __m128i *const p_0 = (__m128i *)&dst[i * dst_stride + j];
            __m128i *const p_1 =
                (__m128i *)&dst[i * dst_stride + j + dst_stride];
            *(uint16_t *)p_0 = (uint16_t)_mm_cvtsi128_si32(res_0);
            *(uint16_t *)p_1 = (uint16_t)_mm_cvtsi128_si32(res_1);
          }
        }
        s[0] = s[1];
        s[1] = s[2];
        s[3] = s[4];
        s[4] = s[5];
      }
    }
  } else if (vert_tap == 12) {  // vert_tap == 12
    const int fo_vert = filter_params_y->taps / 2 - 1;
    const uint8_t *const src_ptr = src - fo_vert * src_stride;
    const __m256i v_zero = _mm256_setzero_si256();
    right_shift = _mm_cvtsi32_si128(FILTER_BITS);
    right_shift_const = _mm256_set1_epi32((1 << FILTER_BITS) >> 1);

    for (j = 0; j < w; j += 8) {
      const uint8_t *data = &src_ptr[j];
      __m256i src10;

      d[0] = _mm_loadl_epi64((__m128i *)(data + 0 * src_stride));
      d[1] = _mm_loadl_epi64((__m128i *)(data + 1 * src_stride));
      d[2] = _mm_loadl_epi64((__m128i *)(data + 2 * src_stride));
      d[3] = _mm_loadl_epi64((__m128i *)(data + 3 * src_stride));
      d[4] = _mm_loadl_epi64((__m128i *)(data + 4 * src_stride));
      d[5] = _mm_loadl_epi64((__m128i *)(data + 5 * src_stride));
      d[6] = _mm_loadl_epi64((__m128i *)(data + 6 * src_stride));
      d[7] = _mm_loadl_epi64((__m128i *)(data + 7 * src_stride));
      d[8] = _mm_loadl_epi64((__m128i *)(data + 8 * src_stride));
      d[9] = _mm_loadl_epi64((__m128i *)(data + 9 * src_stride));
      // Load lines a and b. Line a to lower 128, line b to upper 128
      const __m256i src_01a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[0]), _mm256_castsi128_si256(d[1]), 0x20);

      const __m256i src_12a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[1]), _mm256_castsi128_si256(d[2]), 0x20);

      const __m256i src_23a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[2]), _mm256_castsi128_si256(d[3]), 0x20);

      const __m256i src_34a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[3]), _mm256_castsi128_si256(d[4]), 0x20);

      const __m256i src_45a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[4]), _mm256_castsi128_si256(d[5]), 0x20);

      const __m256i src_56a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[5]), _mm256_castsi128_si256(d[6]), 0x20);

      const __m256i src_67a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[6]), _mm256_castsi128_si256(d[7]), 0x20);

      const __m256i src_78a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[7]), _mm256_castsi128_si256(d[8]), 0x20);

      const __m256i src_89a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[8]), _mm256_castsi128_si256(d[9]), 0x20);

      src10 = _mm256_castsi128_si256(
          _mm_loadl_epi64((__m128i *)(data + 10 * src_stride)));
      const __m256i src_910a =
          _mm256_permute2x128_si256(_mm256_castsi128_si256(d[9]), src10, 0x20);

      const __m256i src_01 = _mm256_unpacklo_epi8(src_01a, v_zero);
      const __m256i src_12 = _mm256_unpacklo_epi8(src_12a, v_zero);
      const __m256i src_23 = _mm256_unpacklo_epi8(src_23a, v_zero);
      const __m256i src_34 = _mm256_unpacklo_epi8(src_34a, v_zero);
      const __m256i src_45 = _mm256_unpacklo_epi8(src_45a, v_zero);
      const __m256i src_56 = _mm256_unpacklo_epi8(src_56a, v_zero);
      const __m256i src_67 = _mm256_unpacklo_epi8(src_67a, v_zero);
      const __m256i src_78 = _mm256_unpacklo_epi8(src_78a, v_zero);
      const __m256i src_89 = _mm256_unpacklo_epi8(src_89a, v_zero);
      const __m256i src_910 = _mm256_unpacklo_epi8(src_910a, v_zero);

      s[0] = _mm256_unpacklo_epi16(src_01, src_12);
      s[1] = _mm256_unpacklo_epi16(src_23, src_34);
      s[2] = _mm256_unpacklo_epi16(src_45, src_56);
      s[3] = _mm256_unpacklo_epi16(src_67, src_78);
      s[4] = _mm256_unpacklo_epi16(src_89, src_910);

      s[6] = _mm256_unpackhi_epi16(src_01, src_12);
      s[7] = _mm256_unpackhi_epi16(src_23, src_34);
      s[8] = _mm256_unpackhi_epi16(src_45, src_56);
      s[9] = _mm256_unpackhi_epi16(src_67, src_78);
      s[10] = _mm256_unpackhi_epi16(src_89, src_910);

      for (i = 0; i < h; i += 2) {
        data = &src_ptr[i * src_stride + j];
        const __m256i src_1011a = _mm256_permute2x128_si256(
            src10,
            _mm256_castsi128_si256(
                _mm_loadl_epi64((__m128i *)(data + 11 * src_stride))),
            0x20);

        src10 = _mm256_castsi128_si256(
            _mm_loadl_epi64((__m128i *)(data + 12 * src_stride)));

        const __m256i src_1112a = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadl_epi64((__m128i *)(data + 11 * src_stride))),
            src10, 0x20);

        const __m256i src_1011 = _mm256_unpacklo_epi8(src_1011a, v_zero);
        const __m256i src_1112 = _mm256_unpacklo_epi8(src_1112a, v_zero);

        s[5] = _mm256_unpacklo_epi16(src_1011, src_1112);
        s[11] = _mm256_unpackhi_epi16(src_1011, src_1112);

        const __m256i res_lo = convolve_12taps(s, coeffs);

        const __m256i res_32b_lo = _mm256_sra_epi32(
            _mm256_add_epi32(res_lo, right_shift_const), right_shift);
        // 8 bit conversion and saturation to uint8
        __m256i res_16b_lo = _mm256_packs_epi32(res_32b_lo, res_32b_lo);
        __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);

        if (w - j > 4) {
          const __m256i res_hi = convolve_12taps(s + 6, coeffs);

          const __m256i res_32b_hi = _mm256_sra_epi32(
              _mm256_add_epi32(res_hi, right_shift_const), right_shift);
          __m256i res_16b_hi = _mm256_packs_epi32(res_32b_hi, res_32b_hi);
          // 8 bit conversion and saturation to uint8
          __m256i res_8b_hi = _mm256_packus_epi16(res_16b_hi, res_16b_hi);

          __m256i res_a = _mm256_unpacklo_epi32(res_8b_lo, res_8b_hi);

          const __m128i res_0 = _mm256_extracti128_si256(res_a, 0);
          const __m128i res_1 = _mm256_extracti128_si256(res_a, 1);

          _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j], res_0);
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j + dst_stride],
                           res_1);
        } else {
          const __m128i res_0 = _mm256_extracti128_si256(res_8b_lo, 0);
          const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
          if (w - j > 2) {
            *(int *)&dst[i * dst_stride + j] = _mm_cvtsi128_si32(res_0);
            *(int *)&dst[i * dst_stride + j + dst_stride] =
                _mm_cvtsi128_si32(res_1);
          } else {
            *(uint16_t *)&dst[i * dst_stride + j] =
                (uint16_t)_mm_cvtsi128_si32(res_0);
            *(uint16_t *)&dst[i * dst_stride + j + dst_stride] =
                (uint16_t)_mm_cvtsi128_si32(res_1);
          }
        }
        s[0] = s[1];
        s[1] = s[2];
        s[2] = s[3];
        s[3] = s[4];
        s[4] = s[5];

        s[6] = s[7];
        s[7] = s[8];
        s[8] = s[9];
        s[9] = s[10];
        s[10] = s[11];
      }
    }
  } else {
    const int fo_vert = filter_params_y->taps / 2 - 1;
    const uint8_t *const src_ptr = src - fo_vert * src_stride;

    for (j = 0; j < w; j += 16) {
      const uint8_t *data = &src_ptr[j];
      __m256i src6;

      d[0] = _mm_loadu_si128((__m128i *)(data + 0 * src_stride));
      d[1] = _mm_loadu_si128((__m128i *)(data + 1 * src_stride));
      d[2] = _mm_loadu_si128((__m128i *)(data + 2 * src_stride));
      d[3] = _mm_loadu_si128((__m128i *)(data + 3 * src_stride));
      d[4] = _mm_loadu_si128((__m128i *)(data + 4 * src_stride));
      d[5] = _mm_loadu_si128((__m128i *)(data + 5 * src_stride));
      // Load lines a and b. Line a to lower 128, line b to upper 128
      const __m256i src_01a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[0]), _mm256_castsi128_si256(d[1]), 0x20);

      const __m256i src_12a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[1]), _mm256_castsi128_si256(d[2]), 0x20);

      const __m256i src_23a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[2]), _mm256_castsi128_si256(d[3]), 0x20);

      const __m256i src_34a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[3]), _mm256_castsi128_si256(d[4]), 0x20);

      const __m256i src_45a = _mm256_permute2x128_si256(
          _mm256_castsi128_si256(d[4]), _mm256_castsi128_si256(d[5]), 0x20);

      src6 = _mm256_castsi128_si256(
          _mm_loadu_si128((__m128i *)(data + 6 * src_stride)));
      const __m256i src_56a =
          _mm256_permute2x128_si256(_mm256_castsi128_si256(d[5]), src6, 0x20);

      s[0] = _mm256_unpacklo_epi8(src_01a, src_12a);
      s[1] = _mm256_unpacklo_epi8(src_23a, src_34a);
      s[2] = _mm256_unpacklo_epi8(src_45a, src_56a);

      s[4] = _mm256_unpackhi_epi8(src_01a, src_12a);
      s[5] = _mm256_unpackhi_epi8(src_23a, src_34a);
      s[6] = _mm256_unpackhi_epi8(src_45a, src_56a);

      for (i = 0; i < h; i += 2) {
        data = &src_ptr[i * src_stride + j];
        const __m256i src_67a = _mm256_permute2x128_si256(
            src6,
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(data + 7 * src_stride))),
            0x20);

        src6 = _mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(data + 8 * src_stride)));
        const __m256i src_78a = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(data + 7 * src_stride))),
            src6, 0x20);

        s[3] = _mm256_unpacklo_epi8(src_67a, src_78a);
        s[7] = _mm256_unpackhi_epi8(src_67a, src_78a);

        const __m256i res_lo = convolve_lowbd(s, coeffs);

        /* rounding code */
        // shift by F - 1
        const __m256i res_16b_lo = _mm256_sra_epi16(
            _mm256_add_epi16(res_lo, right_shift_const), right_shift);
        // 8 bit conversion and saturation to uint8
        __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);

        if (w - j > 8) {
          const __m256i res_hi = convolve_lowbd(s + 4, coeffs);

          /* rounding code */
          // shift by F - 1
          const __m256i res_16b_hi = _mm256_sra_epi16(
              _mm256_add_epi16(res_hi, right_shift_const), right_shift);
          // 8 bit conversion and saturation to uint8
          __m256i res_8b_hi = _mm256_packus_epi16(res_16b_hi, res_16b_hi);

          __m256i res_a = _mm256_unpacklo_epi64(res_8b_lo, res_8b_hi);

          const __m128i res_0 = _mm256_castsi256_si128(res_a);
          const __m128i res_1 = _mm256_extracti128_si256(res_a, 1);

          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res_0);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j + dst_stride],
                           res_1);
        } else {
          const __m128i res_0 = _mm256_castsi256_si128(res_8b_lo);
          const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
          if (w - j > 4) {
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j], res_0);
            _mm_storel_epi64((__m128i *)&dst[i * dst_stride + j + dst_stride],
                             res_1);
          } else if (w - j > 2) {
            xx_storel_32(&dst[i * dst_stride + j], res_0);
            xx_storel_32(&dst[i * dst_stride + j + dst_stride], res_1);
          } else {
            __m128i *const p_0 = (__m128i *)&dst[i * dst_stride + j];
            __m128i *const p_1 =
                (__m128i *)&dst[i * dst_stride + j + dst_stride];
            *(uint16_t *)p_0 = (uint16_t)_mm_cvtsi128_si32(res_0);
            *(uint16_t *)p_1 = (uint16_t)_mm_cvtsi128_si32(res_1);
          }
        }
        s[0] = s[1];
        s[1] = s[2];
        s[2] = s[3];

        s[4] = s[5];
        s[5] = s[6];
        s[6] = s[7];
      }
    }
  }
}

// Utilities for convolve y AVX2
static INLINE void sr_y_round_store_32_avx2(const __m256i res[2],
                                            uint8_t *const dst) {
  __m256i r[2];

  r[0] = sr_y_round_avx2(res[0]);
  r[1] = sr_y_round_avx2(res[1]);
  convolve_store_32_avx2(r[0], r[1], dst);
}

static INLINE void sr_y_round_store_32x2_avx2(const __m256i res[4],
                                              uint8_t *const dst,
                                              const int32_t dst_stride) {
  sr_y_round_store_32_avx2(res, dst);
  sr_y_round_store_32_avx2(res + 2, dst + dst_stride);
}

static INLINE void sr_y_2tap_32_avx2(const uint8_t *const src,
                                     const __m256i coeffs[1], const __m256i s0,
                                     __m256i *const s1, uint8_t *const dst) {
  __m256i r[2];
  y_convolve_2tap_32_avx2(src, coeffs, s0, s1, r);
  sr_y_round_store_32_avx2(r, dst);
}

void av1_convolve_y_sr_avx2(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t w,
                            int32_t h,
                            const InterpFilterParams *filter_params_y,
                            const int32_t subpel_y_q4) {
  int32_t x, y;
  __m128i coeffs_128[4];
  __m256i coeffs_256[4];

  int vert_tap = get_filter_tap(filter_params_y, subpel_y_q4);

  if (vert_tap == 2) {
    // vert_filt as 2 tap
    const uint8_t *src_ptr = src;

    y = h;

    if (subpel_y_q4 != 8) {
      if (w <= 8) {
        prepare_half_coeffs_2tap_ssse3(filter_params_y, subpel_y_q4,
                                       coeffs_128);

        if (w == 2) {
          __m128i s_16[2];

          s_16[0] = _mm_cvtsi32_si128(*(int16_t *)src_ptr);

          do {
            const __m128i res = y_convolve_2tap_2x2_ssse3(src_ptr, src_stride,
                                                          coeffs_128, s_16);
            const __m128i r = sr_y_round_sse2(res);
            pack_store_2x2_sse2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 4) {
          __m128i s_32[2];

          s_32[0] = _mm_cvtsi32_si128(*(int32_t *)src_ptr);

          do {
            const __m128i res = y_convolve_2tap_4x2_ssse3(src_ptr, src_stride,
                                                          coeffs_128, s_32);
            const __m128i r = sr_y_round_sse2(res);
            pack_store_4x2_sse2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else {
          __m128i s_64[2], s_128[2];

          assert(w == 8);

          s_64[0] = _mm_loadl_epi64((__m128i *)src_ptr);

          do {
            // Note: Faster than binding to AVX2 registers.
            s_64[1] = _mm_loadl_epi64((__m128i *)(src_ptr + src_stride));
            s_128[0] = _mm_unpacklo_epi64(s_64[0], s_64[1]);
            s_64[0] = _mm_loadl_epi64((__m128i *)(src_ptr + 2 * src_stride));
            s_128[1] = _mm_unpacklo_epi64(s_64[1], s_64[0]);
            const __m128i ss0 = _mm_unpacklo_epi8(s_128[0], s_128[1]);
            const __m128i ss1 = _mm_unpackhi_epi8(s_128[0], s_128[1]);
            const __m128i res0 = convolve_2tap_ssse3(&ss0, coeffs_128);
            const __m128i res1 = convolve_2tap_ssse3(&ss1, coeffs_128);
            const __m128i r0 = sr_y_round_sse2(res0);
            const __m128i r1 = sr_y_round_sse2(res1);
            const __m128i d = _mm_packus_epi16(r0, r1);
            _mm_storel_epi64((__m128i *)dst, d);
            _mm_storeh_epi64((__m128i *)(dst + dst_stride), d);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        }
      } else {
        prepare_half_coeffs_2tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

        if (w == 16) {
          __m128i s_128[2];

          s_128[0] = _mm_loadu_si128((__m128i *)src_ptr);

          do {
            __m256i r[2];

            y_convolve_2tap_16x2_avx2(src_ptr, src_stride, coeffs_256, s_128,
                                      r);
            sr_y_round_store_16x2_avx2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 32) {
          __m256i s_256[2];

          s_256[0] = _mm256_loadu_si256((__m256i *)src_ptr);

          do {
            sr_y_2tap_32_avx2(src_ptr + src_stride, coeffs_256, s_256[0],
                              &s_256[1], dst);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride, coeffs_256, s_256[1],
                              &s_256[0], dst + dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 64) {
          __m256i s_256[2][2];

          s_256[0][0] = _mm256_loadu_si256((__m256i *)(src_ptr + 0 * 32));
          s_256[0][1] = _mm256_loadu_si256((__m256i *)(src_ptr + 1 * 32));

          do {
            sr_y_2tap_32_avx2(src_ptr + src_stride, coeffs_256, s_256[0][0],
                              &s_256[1][0], dst);
            sr_y_2tap_32_avx2(src_ptr + src_stride + 32, coeffs_256,
                              s_256[0][1], &s_256[1][1], dst + 32);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride, coeffs_256, s_256[1][0],
                              &s_256[0][0], dst + dst_stride);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride + 32, coeffs_256,
                              s_256[1][1], &s_256[0][1], dst + dst_stride + 32);

            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else {
          __m256i s_256[2][4];

          assert(w == 128);

          s_256[0][0] = _mm256_loadu_si256((__m256i *)(src_ptr + 0 * 32));
          s_256[0][1] = _mm256_loadu_si256((__m256i *)(src_ptr + 1 * 32));
          s_256[0][2] = _mm256_loadu_si256((__m256i *)(src_ptr + 2 * 32));
          s_256[0][3] = _mm256_loadu_si256((__m256i *)(src_ptr + 3 * 32));

          do {
            sr_y_2tap_32_avx2(src_ptr + src_stride, coeffs_256, s_256[0][0],
                              &s_256[1][0], dst);
            sr_y_2tap_32_avx2(src_ptr + src_stride + 1 * 32, coeffs_256,
                              s_256[0][1], &s_256[1][1], dst + 1 * 32);
            sr_y_2tap_32_avx2(src_ptr + src_stride + 2 * 32, coeffs_256,
                              s_256[0][2], &s_256[1][2], dst + 2 * 32);
            sr_y_2tap_32_avx2(src_ptr + src_stride + 3 * 32, coeffs_256,
                              s_256[0][3], &s_256[1][3], dst + 3 * 32);

            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride, coeffs_256, s_256[1][0],
                              &s_256[0][0], dst + dst_stride);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride + 1 * 32, coeffs_256,
                              s_256[1][1], &s_256[0][1],
                              dst + dst_stride + 1 * 32);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride + 2 * 32, coeffs_256,
                              s_256[1][2], &s_256[0][2],
                              dst + dst_stride + 2 * 32);
            sr_y_2tap_32_avx2(src_ptr + 2 * src_stride + 3 * 32, coeffs_256,
                              s_256[1][3], &s_256[0][3],
                              dst + dst_stride + 3 * 32);

            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        }
      }
    } else {
      // average to get half pel
      if (w <= 8) {
        if (w == 2) {
          __m128i s_16[2];

          s_16[0] = _mm_cvtsi32_si128(*(int16_t *)src_ptr);

          do {
            s_16[1] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + src_stride));
            const __m128i d0 = _mm_avg_epu8(s_16[0], s_16[1]);
            *(int16_t *)dst = (int16_t)_mm_cvtsi128_si32(d0);
            s_16[0] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 2 * src_stride));
            const __m128i d1 = _mm_avg_epu8(s_16[1], s_16[0]);
            *(int16_t *)(dst + dst_stride) = (int16_t)_mm_cvtsi128_si32(d1);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 4) {
          __m128i s_32[2];

          s_32[0] = _mm_cvtsi32_si128(*(int32_t *)src_ptr);

          do {
            s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + src_stride));
            const __m128i d0 = _mm_avg_epu8(s_32[0], s_32[1]);
            xx_storel_32(dst, d0);
            *(uint32_t *)dst = _mm_cvtsi128_si32(d0);
            s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 2 * src_stride));
            const __m128i d1 = _mm_avg_epu8(s_32[1], s_32[0]);
            xx_storel_32(dst + dst_stride, d1);
            *(uint32_t *)(dst + dst_stride) = _mm_cvtsi128_si32(d1);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else {
          __m128i s_64[2];

          assert(w == 8);

          s_64[0] = _mm_loadl_epi64((__m128i *)src_ptr);

          do {
            // Note: Faster than binding to AVX2 registers.
            s_64[1] = _mm_loadl_epi64((__m128i *)(src_ptr + src_stride));
            const __m128i d0 = _mm_avg_epu8(s_64[0], s_64[1]);
            _mm_storel_epi64((__m128i *)dst, d0);
            s_64[0] = _mm_loadl_epi64((__m128i *)(src_ptr + 2 * src_stride));
            const __m128i d1 = _mm_avg_epu8(s_64[1], s_64[0]);
            _mm_storel_epi64((__m128i *)(dst + dst_stride), d1);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        }
      } else if (w == 16) {
        __m128i s_128[2];

        s_128[0] = _mm_loadu_si128((__m128i *)src_ptr);

        do {
          s_128[1] = _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
          const __m128i d0 = _mm_avg_epu8(s_128[0], s_128[1]);
          _mm_storeu_si128((__m128i *)dst, d0);
          s_128[0] = _mm_loadu_si128((__m128i *)(src_ptr + 2 * src_stride));
          const __m128i d1 = _mm_avg_epu8(s_128[1], s_128[0]);
          _mm_storeu_si128((__m128i *)(dst + dst_stride), d1);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 32) {
        __m256i s_256[2];

        s_256[0] = _mm256_loadu_si256((__m256i *)src_ptr);

        do {
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride, s_256[0], &s_256[1], dst);
          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride, s_256[1], &s_256[0],
                                dst + dst_stride);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 64) {
        __m256i s_256[2][2];

        s_256[0][0] = _mm256_loadu_si256((__m256i *)(src_ptr + 0 * 32));
        s_256[0][1] = _mm256_loadu_si256((__m256i *)(src_ptr + 1 * 32));

        do {
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride, s_256[0][0], &s_256[1][0],
                                dst);
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride + 32, s_256[0][1],
                                &s_256[1][1], dst + 32);

          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride, s_256[1][0],
                                &s_256[0][0], dst + dst_stride);
          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride + 32, s_256[1][1],
                                &s_256[0][1], dst + dst_stride + 32);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m256i s_256[2][4];

        assert(w == 128);

        s_256[0][0] = _mm256_loadu_si256((__m256i *)(src_ptr + 0 * 32));
        s_256[0][1] = _mm256_loadu_si256((__m256i *)(src_ptr + 1 * 32));
        s_256[0][2] = _mm256_loadu_si256((__m256i *)(src_ptr + 2 * 32));
        s_256[0][3] = _mm256_loadu_si256((__m256i *)(src_ptr + 3 * 32));

        do {
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride, s_256[0][0], &s_256[1][0],
                                dst);
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride + 1 * 32, s_256[0][1],
                                &s_256[1][1], dst + 1 * 32);
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride + 2 * 32, s_256[0][2],
                                &s_256[1][2], dst + 2 * 32);
          sr_y_2tap_32_avg_avx2(src_ptr + src_stride + 3 * 32, s_256[0][3],
                                &s_256[1][3], dst + 3 * 32);

          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride, s_256[1][0],
                                &s_256[0][0], dst + dst_stride);
          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride + 1 * 32, s_256[1][1],
                                &s_256[0][1], dst + dst_stride + 1 * 32);
          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride + 2 * 32, s_256[1][2],
                                &s_256[0][2], dst + dst_stride + 2 * 32);
          sr_y_2tap_32_avg_avx2(src_ptr + 2 * src_stride + 3 * 32, s_256[1][3],
                                &s_256[0][3], dst + dst_stride + 3 * 32);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      }
    }
  } else if (vert_tap == 4) {
    // vert_filt as 4 tap
    const uint8_t *src_ptr = src - src_stride;

    y = h;

    if (w <= 4) {
      prepare_half_coeffs_4tap_ssse3(filter_params_y, subpel_y_q4, coeffs_128);

      if (w == 2) {
        __m128i s_16[4], ss_128[2];

        s_16[0] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 0 * src_stride));
        s_16[1] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 1 * src_stride));
        s_16[2] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 2 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi16(s_16[0], s_16[1]);
        const __m128i src12 = _mm_unpacklo_epi16(s_16[1], s_16[2]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);

        do {
          src_ptr += 2 * src_stride;
          const __m128i res = y_convolve_4tap_2x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_16, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_2x2_sse2(r, dst, dst_stride);

          ss_128[0] = ss_128[1];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m128i s_32[4], ss_128[2];

        assert(w == 4);

        s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 0 * src_stride));
        s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 1 * src_stride));
        s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 2 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
        const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);

        do {
          src_ptr += 2 * src_stride;
          const __m128i res = y_convolve_4tap_4x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_32, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_4x2_sse2(r, dst, dst_stride);

          ss_128[0] = ss_128[1];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      }
    } else {
      prepare_half_coeffs_4tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

      if (w == 8) {
        __m128i s_64[4];
        __m256i ss_256[2];

        s_64[0] = _mm_loadl_epi64((__m128i *)(src_ptr + 0 * src_stride));
        s_64[1] = _mm_loadl_epi64((__m128i *)(src_ptr + 1 * src_stride));
        s_64[2] = _mm_loadl_epi64((__m128i *)(src_ptr + 2 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_64[0], s_64[1]);
        const __m256i src12 = _mm256_setr_m128i(s_64[1], s_64[2]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);

        do {
          src_ptr += 2 * src_stride;
          const __m256i res = y_convolve_4tap_8x2_avx2(
              src_ptr, src_stride, coeffs_256, s_64, ss_256);
          sr_y_round_store_8x2_avx2(res, dst, dst_stride);

          ss_256[0] = ss_256[1];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        __m128i s_128[4];
        __m256i ss_256[4], r[2];

        s_128[0] = _mm_loadu_si128((__m128i *)(src_ptr + 0 * src_stride));
        s_128[1] = _mm_loadu_si128((__m128i *)(src_ptr + 1 * src_stride));
        s_128[2] = _mm_loadu_si128((__m128i *)(src_ptr + 2 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_128[0], s_128[1]);
        const __m256i src12 = _mm256_setr_m128i(s_128[1], s_128[2]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);
        ss_256[2] = _mm256_unpackhi_epi8(src01, src12);

        do {
          src_ptr += 2 * src_stride;
          y_convolve_4tap_16x2_avx2(src_ptr, src_stride, coeffs_256, s_128,
                                    ss_256, r);
          sr_y_round_store_16x2_avx2(r, dst, dst_stride);

          ss_256[0] = ss_256[1];
          ss_256[2] = ss_256[3];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 32) {
        // AV1 standard won't have 32x4 case.
        // This only favors some optimization feature which
        // subsamples 32x8 to 32x4 and triggers 4-tap filter.

        __m256i s_256[4], ss_256[4], tt_256[4], r[4];

        s_256[0] = _mm256_loadu_si256((__m256i *)(src_ptr + 0 * src_stride));
        s_256[1] = _mm256_loadu_si256((__m256i *)(src_ptr + 1 * src_stride));
        s_256[2] = _mm256_loadu_si256((__m256i *)(src_ptr + 2 * src_stride));

        ss_256[0] = _mm256_unpacklo_epi8(s_256[0], s_256[1]);
        ss_256[2] = _mm256_unpackhi_epi8(s_256[0], s_256[1]);

        tt_256[0] = _mm256_unpacklo_epi8(s_256[1], s_256[2]);
        tt_256[2] = _mm256_unpackhi_epi8(s_256[1], s_256[2]);

        do {
          src_ptr += 2 * src_stride;
          y_convolve_4tap_32x2_avx2(src_ptr, src_stride, coeffs_256, s_256,
                                    ss_256, tt_256, r);
          sr_y_round_store_32x2_avx2(r, dst, dst_stride);

          ss_256[0] = ss_256[1];
          ss_256[2] = ss_256[3];

          tt_256[0] = tt_256[1];
          tt_256[2] = tt_256[3];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        assert(!(w % 32));

        __m256i s_256[4], ss_256[4], tt_256[4], r[4];
        x = 0;
        do {
          const uint8_t *s = src_ptr + x;
          uint8_t *d = dst + x;
          s_256[0] = _mm256_loadu_si256((__m256i *)(s + 0 * src_stride));
          s_256[1] = _mm256_loadu_si256((__m256i *)(s + 1 * src_stride));
          s_256[2] = _mm256_loadu_si256((__m256i *)(s + 2 * src_stride));

          ss_256[0] = _mm256_unpacklo_epi8(s_256[0], s_256[1]);
          ss_256[2] = _mm256_unpackhi_epi8(s_256[0], s_256[1]);

          tt_256[0] = _mm256_unpacklo_epi8(s_256[1], s_256[2]);
          tt_256[2] = _mm256_unpackhi_epi8(s_256[1], s_256[2]);

          y = h;
          do {
            s += 2 * src_stride;
            y_convolve_4tap_32x2_avx2(s, src_stride, coeffs_256, s_256, ss_256,
                                      tt_256, r);
            sr_y_round_store_32x2_avx2(r, d, dst_stride);

            ss_256[0] = ss_256[1];
            ss_256[2] = ss_256[3];

            tt_256[0] = tt_256[1];
            tt_256[2] = tt_256[3];
            d += 2 * dst_stride;
            y -= 2;
          } while (y);
          x += 32;
        } while (x < w);
      }
    }
  } else if (vert_tap == 6) {
    // vert_filt as 6 tap
    const uint8_t *src_ptr = src - 2 * src_stride;

    if (w <= 4) {
      prepare_half_coeffs_6tap_ssse3(filter_params_y, subpel_y_q4, coeffs_128);

      y = h;

      if (w == 2) {
        __m128i s_16[6], ss_128[3];

        s_16[0] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 0 * src_stride));
        s_16[1] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 1 * src_stride));
        s_16[2] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 2 * src_stride));
        s_16[3] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 3 * src_stride));
        s_16[4] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 4 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi16(s_16[0], s_16[1]);
        const __m128i src12 = _mm_unpacklo_epi16(s_16[1], s_16[2]);
        const __m128i src23 = _mm_unpacklo_epi16(s_16[2], s_16[3]);
        const __m128i src34 = _mm_unpacklo_epi16(s_16[3], s_16[4]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);
        ss_128[1] = _mm_unpacklo_epi8(src23, src34);

        do {
          src_ptr += 2 * src_stride;
          const __m128i res = y_convolve_6tap_2x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_16, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_2x2_sse2(r, dst, dst_stride);

          ss_128[0] = ss_128[1];
          ss_128[1] = ss_128[2];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m128i s_32[6], ss_128[3];

        assert(w == 4);

        s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 0 * src_stride));
        s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 1 * src_stride));
        s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 2 * src_stride));
        s_32[3] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 3 * src_stride));
        s_32[4] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 4 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
        const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);
        const __m128i src23 = _mm_unpacklo_epi32(s_32[2], s_32[3]);
        const __m128i src34 = _mm_unpacklo_epi32(s_32[3], s_32[4]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);
        ss_128[1] = _mm_unpacklo_epi8(src23, src34);

        do {
          src_ptr += 2 * src_stride;
          const __m128i res = y_convolve_6tap_4x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_32, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_4x2_sse2(r, dst, dst_stride);

          ss_128[0] = ss_128[1];
          ss_128[1] = ss_128[2];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      }
    } else {
      prepare_half_coeffs_6tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

      if (w == 8) {
        __m128i s_64[6];
        __m256i ss_256[3];

        s_64[0] = _mm_loadl_epi64((__m128i *)(src_ptr + 0 * src_stride));
        s_64[1] = _mm_loadl_epi64((__m128i *)(src_ptr + 1 * src_stride));
        s_64[2] = _mm_loadl_epi64((__m128i *)(src_ptr + 2 * src_stride));
        s_64[3] = _mm_loadl_epi64((__m128i *)(src_ptr + 3 * src_stride));
        s_64[4] = _mm_loadl_epi64((__m128i *)(src_ptr + 4 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_64[0], s_64[1]);
        const __m256i src12 = _mm256_setr_m128i(s_64[1], s_64[2]);
        const __m256i src23 = _mm256_setr_m128i(s_64[2], s_64[3]);
        const __m256i src34 = _mm256_setr_m128i(s_64[3], s_64[4]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);
        ss_256[1] = _mm256_unpacklo_epi8(src23, src34);

        y = h;
        do {
          src_ptr += 2 * src_stride;
          const __m256i res = y_convolve_6tap_8x2_avx2(
              src_ptr, src_stride, coeffs_256, s_64, ss_256);
          sr_y_round_store_8x2_avx2(res, dst, dst_stride);

          ss_256[0] = ss_256[1];
          ss_256[1] = ss_256[2];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        __m128i s_128[6];
        __m256i ss_256[6], r[2];

        s_128[0] = _mm_loadu_si128((__m128i *)(src_ptr + 0 * src_stride));
        s_128[1] = _mm_loadu_si128((__m128i *)(src_ptr + 1 * src_stride));
        s_128[2] = _mm_loadu_si128((__m128i *)(src_ptr + 2 * src_stride));
        s_128[3] = _mm_loadu_si128((__m128i *)(src_ptr + 3 * src_stride));
        s_128[4] = _mm_loadu_si128((__m128i *)(src_ptr + 4 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_128[0], s_128[1]);
        const __m256i src12 = _mm256_setr_m128i(s_128[1], s_128[2]);
        const __m256i src23 = _mm256_setr_m128i(s_128[2], s_128[3]);
        const __m256i src34 = _mm256_setr_m128i(s_128[3], s_128[4]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);
        ss_256[1] = _mm256_unpacklo_epi8(src23, src34);

        ss_256[3] = _mm256_unpackhi_epi8(src01, src12);
        ss_256[4] = _mm256_unpackhi_epi8(src23, src34);

        y = h;
        do {
          src_ptr += 2 * src_stride;
          y_convolve_6tap_16x2_avx2(src_ptr, src_stride, coeffs_256, s_128,
                                    ss_256, r);
          sr_y_round_store_16x2_avx2(r, dst, dst_stride);

          ss_256[0] = ss_256[1];
          ss_256[1] = ss_256[2];

          ss_256[3] = ss_256[4];
          ss_256[4] = ss_256[5];
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m256i s_256[6], ss_256[6], tt_256[6], r[4];

        assert(!(w % 32));

        x = 0;
        do {
          const uint8_t *s = src_ptr + x;
          uint8_t *d = dst + x;

          s_256[0] = _mm256_loadu_si256((__m256i *)(s + 0 * src_stride));
          s_256[1] = _mm256_loadu_si256((__m256i *)(s + 1 * src_stride));
          s_256[2] = _mm256_loadu_si256((__m256i *)(s + 2 * src_stride));
          s_256[3] = _mm256_loadu_si256((__m256i *)(s + 3 * src_stride));
          s_256[4] = _mm256_loadu_si256((__m256i *)(s + 4 * src_stride));

          ss_256[0] = _mm256_unpacklo_epi8(s_256[0], s_256[1]);
          ss_256[1] = _mm256_unpacklo_epi8(s_256[2], s_256[3]);
          ss_256[3] = _mm256_unpackhi_epi8(s_256[0], s_256[1]);
          ss_256[4] = _mm256_unpackhi_epi8(s_256[2], s_256[3]);

          tt_256[0] = _mm256_unpacklo_epi8(s_256[1], s_256[2]);
          tt_256[1] = _mm256_unpacklo_epi8(s_256[3], s_256[4]);
          tt_256[3] = _mm256_unpackhi_epi8(s_256[1], s_256[2]);
          tt_256[4] = _mm256_unpackhi_epi8(s_256[3], s_256[4]);

          y = h;
          do {
            s += 2 * src_stride;
            y_convolve_6tap_32x2_avx2(s, src_stride, coeffs_256, s_256, ss_256,
                                      tt_256, r);
            sr_y_round_store_32x2_avx2(r, d, dst_stride);

            ss_256[0] = ss_256[1];
            ss_256[1] = ss_256[2];
            ss_256[3] = ss_256[4];
            ss_256[4] = ss_256[5];

            tt_256[0] = tt_256[1];
            tt_256[1] = tt_256[2];
            tt_256[3] = tt_256[4];
            tt_256[4] = tt_256[5];
            d += 2 * dst_stride;
            y -= 2;
          } while (y);

          x += 32;
        } while (x < w);
      }
    }
  } else if (vert_tap == 8) {
    // vert_filt as 8 tap
    const uint8_t *src_ptr = src - 3 * src_stride;

    if (w <= 4) {
      prepare_half_coeffs_8tap_ssse3(filter_params_y, subpel_y_q4, coeffs_128);

      y = h;

      if (w == 2) {
        __m128i s_16[8], ss_128[4];

        s_16[0] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 0 * src_stride));
        s_16[1] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 1 * src_stride));
        s_16[2] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 2 * src_stride));
        s_16[3] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 3 * src_stride));
        s_16[4] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 4 * src_stride));
        s_16[5] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 5 * src_stride));
        s_16[6] = _mm_cvtsi32_si128(*(int16_t *)(src_ptr + 6 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi16(s_16[0], s_16[1]);
        const __m128i src12 = _mm_unpacklo_epi16(s_16[1], s_16[2]);
        const __m128i src23 = _mm_unpacklo_epi16(s_16[2], s_16[3]);
        const __m128i src34 = _mm_unpacklo_epi16(s_16[3], s_16[4]);
        const __m128i src45 = _mm_unpacklo_epi16(s_16[4], s_16[5]);
        const __m128i src56 = _mm_unpacklo_epi16(s_16[5], s_16[6]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);
        ss_128[1] = _mm_unpacklo_epi8(src23, src34);
        ss_128[2] = _mm_unpacklo_epi8(src45, src56);

        do {
          const __m128i res = y_convolve_8tap_2x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_16, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_2x2_sse2(r, dst, dst_stride);
          ss_128[0] = ss_128[1];
          ss_128[1] = ss_128[2];
          ss_128[2] = ss_128[3];
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m128i s_32[8], ss_128[4];

        assert(w == 4);

        s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 0 * src_stride));
        s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 1 * src_stride));
        s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 2 * src_stride));
        s_32[3] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 3 * src_stride));
        s_32[4] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 4 * src_stride));
        s_32[5] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 5 * src_stride));
        s_32[6] = _mm_cvtsi32_si128(*(int32_t *)(src_ptr + 6 * src_stride));

        const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
        const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);
        const __m128i src23 = _mm_unpacklo_epi32(s_32[2], s_32[3]);
        const __m128i src34 = _mm_unpacklo_epi32(s_32[3], s_32[4]);
        const __m128i src45 = _mm_unpacklo_epi32(s_32[4], s_32[5]);
        const __m128i src56 = _mm_unpacklo_epi32(s_32[5], s_32[6]);

        ss_128[0] = _mm_unpacklo_epi8(src01, src12);
        ss_128[1] = _mm_unpacklo_epi8(src23, src34);
        ss_128[2] = _mm_unpacklo_epi8(src45, src56);

        do {
          const __m128i res = y_convolve_8tap_4x2_ssse3(
              src_ptr, src_stride, coeffs_128, s_32, ss_128);
          const __m128i r = sr_y_round_sse2(res);
          pack_store_4x2_sse2(r, dst, dst_stride);
          ss_128[0] = ss_128[1];
          ss_128[1] = ss_128[2];
          ss_128[2] = ss_128[3];
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      }
    } else {
      prepare_half_coeffs_8tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

      if (w == 8) {
        __m128i s_64[8];
        __m256i ss_256[4];

        s_64[0] = _mm_loadl_epi64((__m128i *)(src_ptr + 0 * src_stride));
        s_64[1] = _mm_loadl_epi64((__m128i *)(src_ptr + 1 * src_stride));
        s_64[2] = _mm_loadl_epi64((__m128i *)(src_ptr + 2 * src_stride));
        s_64[3] = _mm_loadl_epi64((__m128i *)(src_ptr + 3 * src_stride));
        s_64[4] = _mm_loadl_epi64((__m128i *)(src_ptr + 4 * src_stride));
        s_64[5] = _mm_loadl_epi64((__m128i *)(src_ptr + 5 * src_stride));
        s_64[6] = _mm_loadl_epi64((__m128i *)(src_ptr + 6 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_64[0], s_64[1]);
        const __m256i src12 = _mm256_setr_m128i(s_64[1], s_64[2]);
        const __m256i src23 = _mm256_setr_m128i(s_64[2], s_64[3]);
        const __m256i src34 = _mm256_setr_m128i(s_64[3], s_64[4]);
        const __m256i src45 = _mm256_setr_m128i(s_64[4], s_64[5]);
        const __m256i src56 = _mm256_setr_m128i(s_64[5], s_64[6]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);
        ss_256[1] = _mm256_unpacklo_epi8(src23, src34);
        ss_256[2] = _mm256_unpacklo_epi8(src45, src56);

        y = h;
        do {
          const __m256i res = y_convolve_8tap_8x2_avx2(
              src_ptr, src_stride, coeffs_256, s_64, ss_256);
          sr_y_round_store_8x2_avx2(res, dst, dst_stride);
          ss_256[0] = ss_256[1];
          ss_256[1] = ss_256[2];
          ss_256[2] = ss_256[3];
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        __m128i s_128[8];
        __m256i ss_256[8], r[2];

        s_128[0] = _mm_loadu_si128((__m128i *)(src_ptr + 0 * src_stride));
        s_128[1] = _mm_loadu_si128((__m128i *)(src_ptr + 1 * src_stride));
        s_128[2] = _mm_loadu_si128((__m128i *)(src_ptr + 2 * src_stride));
        s_128[3] = _mm_loadu_si128((__m128i *)(src_ptr + 3 * src_stride));
        s_128[4] = _mm_loadu_si128((__m128i *)(src_ptr + 4 * src_stride));
        s_128[5] = _mm_loadu_si128((__m128i *)(src_ptr + 5 * src_stride));
        s_128[6] = _mm_loadu_si128((__m128i *)(src_ptr + 6 * src_stride));

        // Load lines a and b. Line a to lower 128, line b to upper 128
        const __m256i src01 = _mm256_setr_m128i(s_128[0], s_128[1]);
        const __m256i src12 = _mm256_setr_m128i(s_128[1], s_128[2]);
        const __m256i src23 = _mm256_setr_m128i(s_128[2], s_128[3]);
        const __m256i src34 = _mm256_setr_m128i(s_128[3], s_128[4]);
        const __m256i src45 = _mm256_setr_m128i(s_128[4], s_128[5]);
        const __m256i src56 = _mm256_setr_m128i(s_128[5], s_128[6]);

        ss_256[0] = _mm256_unpacklo_epi8(src01, src12);
        ss_256[1] = _mm256_unpacklo_epi8(src23, src34);
        ss_256[2] = _mm256_unpacklo_epi8(src45, src56);

        ss_256[4] = _mm256_unpackhi_epi8(src01, src12);
        ss_256[5] = _mm256_unpackhi_epi8(src23, src34);
        ss_256[6] = _mm256_unpackhi_epi8(src45, src56);

        y = h;
        do {
          y_convolve_8tap_16x2_avx2(src_ptr, src_stride, coeffs_256, s_128,
                                    ss_256, r);
          sr_y_round_store_16x2_avx2(r, dst, dst_stride);

          ss_256[0] = ss_256[1];
          ss_256[1] = ss_256[2];
          ss_256[2] = ss_256[3];

          ss_256[4] = ss_256[5];
          ss_256[5] = ss_256[6];
          ss_256[6] = ss_256[7];
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else {
        __m256i s_256[8], ss_256[8], tt_256[8], r[4];

        assert(!(w % 32));

        x = 0;
        do {
          const uint8_t *s = src_ptr + x;
          uint8_t *d = dst + x;

          s_256[0] = _mm256_loadu_si256((__m256i *)(s + 0 * src_stride));
          s_256[1] = _mm256_loadu_si256((__m256i *)(s + 1 * src_stride));
          s_256[2] = _mm256_loadu_si256((__m256i *)(s + 2 * src_stride));
          s_256[3] = _mm256_loadu_si256((__m256i *)(s + 3 * src_stride));
          s_256[4] = _mm256_loadu_si256((__m256i *)(s + 4 * src_stride));
          s_256[5] = _mm256_loadu_si256((__m256i *)(s + 5 * src_stride));
          s_256[6] = _mm256_loadu_si256((__m256i *)(s + 6 * src_stride));

          ss_256[0] = _mm256_unpacklo_epi8(s_256[0], s_256[1]);
          ss_256[1] = _mm256_unpacklo_epi8(s_256[2], s_256[3]);
          ss_256[2] = _mm256_unpacklo_epi8(s_256[4], s_256[5]);
          ss_256[4] = _mm256_unpackhi_epi8(s_256[0], s_256[1]);
          ss_256[5] = _mm256_unpackhi_epi8(s_256[2], s_256[3]);
          ss_256[6] = _mm256_unpackhi_epi8(s_256[4], s_256[5]);

          tt_256[0] = _mm256_unpacklo_epi8(s_256[1], s_256[2]);
          tt_256[1] = _mm256_unpacklo_epi8(s_256[3], s_256[4]);
          tt_256[2] = _mm256_unpacklo_epi8(s_256[5], s_256[6]);
          tt_256[4] = _mm256_unpackhi_epi8(s_256[1], s_256[2]);
          tt_256[5] = _mm256_unpackhi_epi8(s_256[3], s_256[4]);
          tt_256[6] = _mm256_unpackhi_epi8(s_256[5], s_256[6]);

          y = h;
          do {
            y_convolve_8tap_32x2_avx2(s, src_stride, coeffs_256, s_256, ss_256,
                                      tt_256, r);
            sr_y_round_store_32x2_avx2(r, d, dst_stride);

            ss_256[0] = ss_256[1];
            ss_256[1] = ss_256[2];
            ss_256[2] = ss_256[3];
            ss_256[4] = ss_256[5];
            ss_256[5] = ss_256[6];
            ss_256[6] = ss_256[7];

            tt_256[0] = tt_256[1];
            tt_256[1] = tt_256[2];
            tt_256[2] = tt_256[3];
            tt_256[4] = tt_256[5];
            tt_256[5] = tt_256[6];
            tt_256[6] = tt_256[7];
            s += 2 * src_stride;
            d += 2 * dst_stride;
            y -= 2;
          } while (y);

          x += 32;
        } while (x < w);
      }
    }
  } else {
    assert(vert_tap == 12);
    av1_convolve_y_sr_avx2_general(src, src_stride, dst, dst_stride, w, h,
                                   filter_params_y, subpel_y_q4);
  }
}

static AOM_INLINE void av1_convolve_x_sr_avx2_general(
    const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn,
    ConvolveParams *conv_params) {
  const int bits = FILTER_BITS - conv_params->round_0;
  const __m128i round_shift = _mm_cvtsi32_si128(bits);
  __m256i round_0_const =
      _mm256_set1_epi16((1 << (conv_params->round_0 - 1)) >> 1);
  __m128i round_0_shift = _mm_cvtsi32_si128(conv_params->round_0 - 1);
  __m256i round_const = _mm256_set1_epi16((1 << bits) >> 1);
  int i, horiz_tap = get_filter_tap(filter_params_x, subpel_x_qn);

  assert(bits >= 0);
  assert((FILTER_BITS - conv_params->round_1) >= 0 ||
         ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));
  assert(conv_params->round_0 > 0);

  __m256i coeffs[6], filt[4];
  filt[0] = _mm256_load_si256((__m256i const *)(filt_global_avx2));
  filt[1] = _mm256_load_si256((__m256i const *)(filt_global_avx2 + 32));

  if (horiz_tap == 6)
    prepare_coeffs_6t_lowbd(filter_params_x, subpel_x_qn, coeffs);
  else if (horiz_tap == 12) {
    prepare_coeffs_12taps(filter_params_x, subpel_x_qn, coeffs);
  } else {
    prepare_coeffs_lowbd(filter_params_x, subpel_x_qn, coeffs);
  }

  // horz_filt as 4 tap
  if (horiz_tap == 4) {
    const int fo_horiz = 1;
    const uint8_t *const src_ptr = src - fo_horiz;
    if (w <= 8) {
      for (i = 0; i < h; i += 2) {
        const __m256i data = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride]))),
            _mm256_castsi128_si256(_mm_loadu_si128(
                (__m128i *)(&src_ptr[i * src_stride + src_stride]))),
            0x20);

        __m256i res_16b = convolve_lowbd_x_4tap(data, coeffs + 1, filt);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                   round_0_shift);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                   round_shift);

        /* rounding code */
        // 8 bit conversion and saturation to uint8
        __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

        const __m128i res_0 = _mm256_castsi256_si128(res_8b);
        const __m128i res_1 = _mm256_extracti128_si256(res_8b, 1);

        if (w > 4) {
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride], res_0);
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride + dst_stride], res_1);
        } else if (w > 2) {
          xx_storel_32(&dst[i * dst_stride], res_0);
          xx_storel_32(&dst[i * dst_stride + dst_stride], res_1);
        } else {
          __m128i *const p_0 = (__m128i *)&dst[i * dst_stride];
          __m128i *const p_1 = (__m128i *)&dst[i * dst_stride + dst_stride];
          *(uint16_t *)p_0 = (uint16_t)_mm_cvtsi128_si32(res_0);
          *(uint16_t *)p_1 = (uint16_t)_mm_cvtsi128_si32(res_1);
        }
      }
    } else {
      for (i = 0; i < h; ++i) {
        for (int j = 0; j < w; j += 16) {
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 8 9 10 11 12 13 14 15 16 17
          // 18 19 20 21 22 23
          const __m256i data = _mm256_inserti128_si256(
              _mm256_loadu_si256((__m256i *)&src_ptr[(i * src_stride) + j]),
              _mm_loadu_si128((__m128i *)&src_ptr[(i * src_stride) + (j + 8)]),
              1);

          __m256i res_16b = convolve_lowbd_x_4tap(data, coeffs + 1, filt);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                     round_0_shift);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                     round_shift);

          /* rounding code */
          // 8 bit conversion and saturation to uint8
          __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

          // Store values into the destination buffer
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
          res_8b = _mm256_permute4x64_epi64(res_8b, 216);
          __m128i res = _mm256_castsi256_si128(res_8b);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res);
        }
      }
    }
  } else if (horiz_tap == 6) {
    const int fo_horiz = horiz_tap / 2 - 1;
    const uint8_t *const src_ptr = src - fo_horiz;
    filt[2] = _mm256_load_si256((__m256i const *)(filt_global_avx2 + 32 * 2));
    filt[3] = _mm256_load_si256((__m256i const *)(filt_global_avx2 + 32 * 3));

    if (w <= 8) {
      for (i = 0; i < h; i += 2) {
        const __m256i data = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride]))),
            _mm256_castsi128_si256(_mm_loadu_si128(
                (__m128i *)(&src_ptr[i * src_stride + src_stride]))),
            0x20);

        __m256i res_16b = convolve_lowbd_x_6tap(data, coeffs, filt);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                   round_0_shift);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                   round_shift);

        /* rounding code */
        // 8 bit conversion and saturation to uint8
        __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

        const __m128i res_0 = _mm256_castsi256_si128(res_8b);
        const __m128i res_1 = _mm256_extracti128_si256(res_8b, 1);
        if (w > 4) {
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride], res_0);
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride + dst_stride], res_1);
        } else if (w > 2) {
          xx_storel_32(&dst[i * dst_stride], res_0);
          xx_storel_32(&dst[i * dst_stride + dst_stride], res_1);
        } else {
          __m128i *const p_0 = (__m128i *)&dst[i * dst_stride];
          __m128i *const p_1 = (__m128i *)&dst[i * dst_stride + dst_stride];
          *(uint16_t *)p_0 = _mm_cvtsi128_si32(res_0);
          *(uint16_t *)p_1 = _mm_cvtsi128_si32(res_1);
        }
      }
    } else {
      for (i = 0; i < h; ++i) {
        for (int j = 0; j < w; j += 16) {
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 8 9 10 11 12 13 14 15 16 17
          // 18 19 20 21 22 23
          const __m256i data = _mm256_inserti128_si256(
              _mm256_loadu_si256((__m256i *)&src_ptr[(i * src_stride) + j]),
              _mm_loadu_si128((__m128i *)&src_ptr[(i * src_stride) + (j + 8)]),
              1);

          __m256i res_16b = convolve_lowbd_x_6tap(data, coeffs, filt);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                     round_0_shift);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                     round_shift);

          /* rounding code */
          // 8 bit conversion and saturation to uint8
          __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

          // Store values into the destination buffer
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
          res_8b = _mm256_permute4x64_epi64(res_8b, 216);
          __m128i res = _mm256_castsi256_si128(res_8b);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res);
        }
      }
    }
  } else if (horiz_tap == 12) {  // horiz_tap == 12
    const int fo_horiz = filter_params_x->taps / 2 - 1;
    const uint8_t *const src_ptr = src - fo_horiz;
    const __m256i v_zero = _mm256_setzero_si256();
    round_0_const = _mm256_set1_epi32((1 << (conv_params->round_0)) >> 1);
    round_const = _mm256_set1_epi32((1 << bits) >> 1);
    round_0_shift = _mm_cvtsi32_si128(conv_params->round_0);
    __m256i s[6];

    if (w <= 4) {
      for (i = 0; i < h; i += 2) {
        const __m256i data = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride]))),
            _mm256_castsi128_si256(_mm_loadu_si128(
                (__m128i *)(&src_ptr[i * src_stride + src_stride]))),
            0x20);
        // row0 0..7 row1 0..7
        const __m256i s_16l = _mm256_unpacklo_epi8(data, v_zero);
        // row0 8..F row1 8..F
        const __m256i s_16h = _mm256_unpackhi_epi8(data, v_zero);

        // row0 00 00 01 01 .. 03 03 row1 00 00 01 01 .. 03 03
        const __m256i s_ll = _mm256_unpacklo_epi16(s_16l, s_16l);
        // row0 04 04 .. 07 07 row1 04 04 .. 07 07
        const __m256i s_lh = _mm256_unpackhi_epi16(s_16l, s_16l);

        // row0 08 08 09 09 .. 0B 0B row1 08 08 09 09 .. 0B 0B
        const __m256i s_hl = _mm256_unpacklo_epi16(s_16h, s_16h);
        // row0 0C 0C .. 0F 0F row1 0C 0C .. 0F 0F
        const __m256i s_hh = _mm256_unpackhi_epi16(s_16h, s_16h);

        // 00 01 01 02 02 03 03 04 10 11 11 12 12 13 13 14
        s[0] = _mm256_alignr_epi8(s_lh, s_ll, 2);
        // 02 03 03 04 04 05 05 06 12 13 13 14 14 15 15 16
        s[1] = _mm256_alignr_epi8(s_lh, s_ll, 10);
        // 04 05 05 06 06 07 07 08 14 15 15 16 16 17 17 18
        s[2] = _mm256_alignr_epi8(s_hl, s_lh, 2);
        // 06 07 07 08 08 09 09 0A 16 17 17 18 18 19 19 1A
        s[3] = _mm256_alignr_epi8(s_hl, s_lh, 10);
        // 08 09 09 0A 0A 0B 0B 0C 18 19 19 1A 1A 1B 1B 1C
        s[4] = _mm256_alignr_epi8(s_hh, s_hl, 2);
        // 0A 0B 0B 0C 0C 0D 0D 0E 1A 1B 1B 1C 1C 1D 1D 1E
        s[5] = _mm256_alignr_epi8(s_hh, s_hl, 10);

        const __m256i res_lo = convolve_12taps(s, coeffs);

        __m256i res_32b_lo = _mm256_sra_epi32(
            _mm256_add_epi32(res_lo, round_0_const), round_0_shift);

        // 00 01 02 03 10 12 13 14
        res_32b_lo = _mm256_sra_epi32(_mm256_add_epi32(res_32b_lo, round_const),
                                      round_shift);
        // 8 bit conversion and saturation to uint8
        // 00 01 02 03 00 01 02 03 10 11 12 13 10 11 12 13
        __m256i res_16b_lo = _mm256_packs_epi32(res_32b_lo, res_32b_lo);
        // 00 01 02 03 00 01 02 03 00 01 02 03 00 01 02 03
        // 10 11 12 13 10 11 12 13 10 11 12 13 10 11 12 13
        __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);

        // 00 01 02 03 00 01 02 03 00 01 02 03 00 01 02 03
        const __m128i res_0 = _mm256_extracti128_si256(res_8b_lo, 0);
        // 10 11 12 13 10 11 12 13 10 11 12 13 10 11 12 13
        const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
        if (w > 2) {
          // 00 01 02 03
          *(int *)&dst[i * dst_stride] = _mm_cvtsi128_si32(res_0);
          // 10 11 12 13
          *(int *)&dst[i * dst_stride + dst_stride] = _mm_cvtsi128_si32(res_1);
        } else {
          // 00 01
          *(uint16_t *)&dst[i * dst_stride] =
              (uint16_t)_mm_cvtsi128_si32(res_0);
          // 10 11
          *(uint16_t *)&dst[i * dst_stride + dst_stride] =
              (uint16_t)_mm_cvtsi128_si32(res_1);
        }
      }
    } else {
      for (i = 0; i < h; i++) {
        for (int j = 0; j < w; j += 8) {
          const __m256i data = _mm256_permute2x128_si256(
              _mm256_castsi128_si256(
                  _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride + j]))),
              _mm256_castsi128_si256(_mm_loadu_si128(
                  (__m128i *)(&src_ptr[i * src_stride + j + 4]))),
              0x20);
          // row0 0..7 4..B
          const __m256i s_16l = _mm256_unpacklo_epi8(data, v_zero);
          // row0 8..F C..13
          const __m256i s_16h = _mm256_unpackhi_epi8(data, v_zero);

          // row0 00 00 01 01 .. 03 03 04 04 05 05 .. 07 07
          const __m256i s_ll = _mm256_unpacklo_epi16(s_16l, s_16l);
          // row0 04 04 .. 07 07 08 08 .. 0B 0B
          const __m256i s_lh = _mm256_unpackhi_epi16(s_16l, s_16l);

          // row0 08 08 09 09 .. 0B 0B 0C 0C 0D 0D .. 0F 0F
          const __m256i s_hl = _mm256_unpacklo_epi16(s_16h, s_16h);
          // row0 0C 0C 0D 0D .. 0F 0F 10 10 11 11 .. 13 13
          const __m256i s_hh = _mm256_unpackhi_epi16(s_16h, s_16h);

          s[0] = _mm256_alignr_epi8(s_lh, s_ll, 2);
          s[1] = _mm256_alignr_epi8(s_lh, s_ll, 10);
          s[2] = _mm256_alignr_epi8(s_hl, s_lh, 2);
          s[3] = _mm256_alignr_epi8(s_hl, s_lh, 10);
          s[4] = _mm256_alignr_epi8(s_hh, s_hl, 2);
          s[5] = _mm256_alignr_epi8(s_hh, s_hl, 10);

          const __m256i res_lo = convolve_12taps(s, coeffs);

          __m256i res_32b_lo = _mm256_sra_epi32(
              _mm256_add_epi32(res_lo, round_0_const), round_0_shift);

          res_32b_lo = _mm256_sra_epi32(
              _mm256_add_epi32(res_32b_lo, round_const), round_shift);
          // 8 bit conversion and saturation to uint8
          __m256i res_16b_lo = _mm256_packs_epi32(res_32b_lo, res_32b_lo);
          __m256i res_8b_lo = _mm256_packus_epi16(res_16b_lo, res_16b_lo);
          const __m128i res_0 = _mm256_extracti128_si256(res_8b_lo, 0);
          const __m128i res_1 = _mm256_extracti128_si256(res_8b_lo, 1);
          *(int *)&dst[i * dst_stride + j] = _mm_cvtsi128_si32(res_0);
          *(int *)&dst[i * dst_stride + j + 4] = _mm_cvtsi128_si32(res_1);
        }
      }
    }
  } else {
    const int fo_horiz = filter_params_x->taps / 2 - 1;
    const uint8_t *const src_ptr = src - fo_horiz;
    filt[2] = _mm256_load_si256((__m256i const *)(filt_global_avx2 + 32 * 2));
    filt[3] = _mm256_load_si256((__m256i const *)(filt_global_avx2 + 32 * 3));

    if (w <= 8) {
      for (i = 0; i < h; i += 2) {
        const __m256i data = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride]))),
            _mm256_castsi128_si256(_mm_loadu_si128(
                (__m128i *)(&src_ptr[i * src_stride + src_stride]))),
            0x20);

        __m256i res_16b = convolve_lowbd_x(data, coeffs, filt);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                   round_0_shift);

        res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                   round_shift);

        /* rounding code */
        // 8 bit conversion and saturation to uint8
        __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

        const __m128i res_0 = _mm256_castsi256_si128(res_8b);
        const __m128i res_1 = _mm256_extracti128_si256(res_8b, 1);
        if (w > 4) {
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride], res_0);
          _mm_storel_epi64((__m128i *)&dst[i * dst_stride + dst_stride], res_1);
        } else if (w > 2) {
          xx_storel_32(&dst[i * dst_stride], res_0);
          xx_storel_32(&dst[i * dst_stride + dst_stride], res_1);
        } else {
          __m128i *const p_0 = (__m128i *)&dst[i * dst_stride];
          __m128i *const p_1 = (__m128i *)&dst[i * dst_stride + dst_stride];
          *(uint16_t *)p_0 = (uint16_t)_mm_cvtsi128_si32(res_0);
          *(uint16_t *)p_1 = (uint16_t)_mm_cvtsi128_si32(res_1);
        }
      }
    } else {
      for (i = 0; i < h; ++i) {
        for (int j = 0; j < w; j += 16) {
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 8 9 10 11 12 13 14 15 16 17
          // 18 19 20 21 22 23
          const __m256i data = _mm256_inserti128_si256(
              _mm256_loadu_si256((__m256i *)&src_ptr[(i * src_stride) + j]),
              _mm_loadu_si128((__m128i *)&src_ptr[(i * src_stride) + (j + 8)]),
              1);

          __m256i res_16b = convolve_lowbd_x(data, coeffs, filt);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_0_const),
                                     round_0_shift);

          res_16b = _mm256_sra_epi16(_mm256_add_epi16(res_16b, round_const),
                                     round_shift);

          /* rounding code */
          // 8 bit conversion and saturation to uint8
          __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

          // Store values into the destination buffer
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
          res_8b = _mm256_permute4x64_epi64(res_8b, 216);
          __m128i res = _mm256_castsi256_si128(res_8b);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res);
        }
      }
    }
  }
}

static INLINE void sr_x_2tap_32_avx2(const uint8_t *const src,
                                     const __m256i coeffs[1],
                                     uint8_t *const dst) {
  __m256i r[2];

  x_convolve_2tap_32_avx2(src, coeffs, r);
  sr_x_round_store_32_avx2(r, dst);
}

static INLINE void sr_x_6tap_32_avx2(const uint8_t *const src,
                                     const __m256i coeffs[3],
                                     const __m256i filt[3],
                                     uint8_t *const dst) {
  __m256i r[2];

  x_convolve_6tap_32_avx2(src, coeffs, filt, r);
  sr_x_round_store_32_avx2(r, dst);
}

static AOM_FORCE_INLINE void sr_x_8tap_32_avx2(const uint8_t *const src,
                                               const __m256i coeffs[4],
                                               const __m256i filt[4],
                                               uint8_t *const dst) {
  __m256i r[2];

  x_convolve_8tap_32_avx2(src, coeffs, filt, r);
  sr_x_round_store_32_avx2(r, dst);
}

void av1_convolve_x_sr_avx2(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t w,
                            int32_t h,
                            const InterpFilterParams *filter_params_x,
                            const int32_t subpel_x_q4,
                            ConvolveParams *conv_params) {
  int32_t y = h;
  __m128i coeffs_128[4];
  __m256i coeffs_256[4];

  assert(conv_params->round_0 == 3);
  assert((FILTER_BITS - conv_params->round_1) >= 0 ||
         ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));
  (void)conv_params;

  const int horz_tap = get_filter_tap(filter_params_x, subpel_x_q4);

  if (horz_tap == 2) {
    // horz_filt as 2 tap
    const uint8_t *src_ptr = src;

    if (subpel_x_q4 != 8) {
      if (w <= 8) {
        prepare_half_coeffs_2tap_ssse3(filter_params_x, subpel_x_q4,
                                       coeffs_128);

        if (w == 2) {
          do {
            const __m128i res =
                x_convolve_2tap_2x2_sse4_1(src_ptr, src_stride, coeffs_128);
            const __m128i r = sr_x_round_sse2(res);
            pack_store_2x2_sse2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 4) {
          do {
            const __m128i res =
                x_convolve_2tap_4x2_ssse3(src_ptr, src_stride, coeffs_128);
            const __m128i r = sr_x_round_sse2(res);
            pack_store_4x2_sse2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else {
          assert(w == 8);

          do {
            __m128i res[2];

            x_convolve_2tap_8x2_ssse3(src_ptr, src_stride, coeffs_128, res);
            res[0] = sr_x_round_sse2(res[0]);
            res[1] = sr_x_round_sse2(res[1]);
            const __m128i d = _mm_packus_epi16(res[0], res[1]);
            _mm_storel_epi64((__m128i *)dst, d);
            _mm_storeh_epi64((__m128i *)(dst + dst_stride), d);

            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        }
      } else {
        prepare_half_coeffs_2tap_avx2(filter_params_x, subpel_x_q4, coeffs_256);

        if (w == 16) {
          do {
            __m256i r[2];

            x_convolve_2tap_16x2_avx2(src_ptr, src_stride, coeffs_256, r);
            sr_x_round_store_16x2_avx2(r, dst, dst_stride);
            src_ptr += 2 * src_stride;
            dst += 2 * dst_stride;
            y -= 2;
          } while (y);
        } else if (w == 32) {
          do {
            sr_x_2tap_32_avx2(src_ptr, coeffs_256, dst);
            src_ptr += src_stride;
            dst += dst_stride;
          } while (--y);
        } else if (w == 64) {
          do {
            sr_x_2tap_32_avx2(src_ptr + 0 * 32, coeffs_256, dst + 0 * 32);
            sr_x_2tap_32_avx2(src_ptr + 1 * 32, coeffs_256, dst + 1 * 32);
            src_ptr += src_stride;
            dst += dst_stride;
          } while (--y);
        } else {
          assert(w == 128);

          do {
            sr_x_2tap_32_avx2(src_ptr + 0 * 32, coeffs_256, dst + 0 * 32);
            sr_x_2tap_32_avx2(src_ptr + 1 * 32, coeffs_256, dst + 1 * 32);
            sr_x_2tap_32_avx2(src_ptr + 2 * 32, coeffs_256, dst + 2 * 32);
            sr_x_2tap_32_avx2(src_ptr + 3 * 32, coeffs_256, dst + 3 * 32);
            src_ptr += src_stride;
            dst += dst_stride;
          } while (--y);
        }
      }
    } else {
      // average to get half pel
      if (w == 2) {
        do {
          __m128i s_128;

          s_128 = load_u8_4x2_sse4_1(src_ptr, src_stride);
          const __m128i s1 = _mm_srli_si128(s_128, 1);
          const __m128i d = _mm_avg_epu8(s_128, s1);
          *(uint16_t *)dst = (uint16_t)_mm_cvtsi128_si32(d);
          *(uint16_t *)(dst + dst_stride) = _mm_extract_epi16(d, 2);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 4) {
        do {
          __m128i s_128;

          s_128 = load_u8_8x2_sse2(src_ptr, src_stride);
          const __m128i s1 = _mm_srli_si128(s_128, 1);
          const __m128i d = _mm_avg_epu8(s_128, s1);
          xx_storel_32(dst, d);
          *(int32_t *)(dst + dst_stride) = _mm_extract_epi32(d, 2);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 8) {
        do {
          const __m128i s00 = _mm_loadu_si128((__m128i *)src_ptr);
          const __m128i s10 =
              _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
          const __m128i s01 = _mm_srli_si128(s00, 1);
          const __m128i s11 = _mm_srli_si128(s10, 1);
          const __m128i d0 = _mm_avg_epu8(s00, s01);
          const __m128i d1 = _mm_avg_epu8(s10, s11);
          _mm_storel_epi64((__m128i *)dst, d0);
          _mm_storel_epi64((__m128i *)(dst + dst_stride), d1);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        do {
          const __m128i s00 = _mm_loadu_si128((__m128i *)src_ptr);
          const __m128i s01 = _mm_loadu_si128((__m128i *)(src_ptr + 1));
          const __m128i s10 =
              _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
          const __m128i s11 =
              _mm_loadu_si128((__m128i *)(src_ptr + src_stride + 1));
          const __m128i d0 = _mm_avg_epu8(s00, s01);
          const __m128i d1 = _mm_avg_epu8(s10, s11);
          _mm_storeu_si128((__m128i *)dst, d0);
          _mm_storeu_si128((__m128i *)(dst + dst_stride), d1);

          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 32) {
        do {
          sr_x_2tap_32_avg_avx2(src_ptr, dst);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else if (w == 64) {
        do {
          sr_x_2tap_32_avg_avx2(src_ptr + 0 * 32, dst + 0 * 32);
          sr_x_2tap_32_avg_avx2(src_ptr + 1 * 32, dst + 1 * 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else {
        assert(w == 128);

        do {
          sr_x_2tap_32_avg_avx2(src_ptr + 0 * 32, dst + 0 * 32);
          sr_x_2tap_32_avg_avx2(src_ptr + 1 * 32, dst + 1 * 32);
          sr_x_2tap_32_avg_avx2(src_ptr + 2 * 32, dst + 2 * 32);
          sr_x_2tap_32_avg_avx2(src_ptr + 3 * 32, dst + 3 * 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      }
    }
  } else if (horz_tap == 4) {
    // horz_filt as 4 tap
    const uint8_t *src_ptr = src - 1;

    prepare_half_coeffs_4tap_ssse3(filter_params_x, subpel_x_q4, coeffs_128);

    if (w == 2) {
      do {
        const __m128i res =
            x_convolve_4tap_2x2_ssse3(src_ptr, src_stride, coeffs_128);
        const __m128i r = sr_x_round_sse2(res);
        pack_store_2x2_sse2(r, dst, dst_stride);
        src_ptr += 2 * src_stride;
        dst += 2 * dst_stride;
        y -= 2;
      } while (y);
    } else if (w == 4) {
      do {
        const __m128i res =
            x_convolve_4tap_4x2_ssse3(src_ptr, src_stride, coeffs_128);
        const __m128i r = sr_x_round_sse2(res);
        pack_store_4x2_sse2(r, dst, dst_stride);
        src_ptr += 2 * src_stride;
        dst += 2 * dst_stride;
        y -= 2;
      } while (y);
    } else if (w == 8) {
      // TODO(chiyotsai@google.com): Reuse the old SIMD code here. Need to
      // rewrite this for better performance later.
      __m256i filt_256[2];
      prepare_coeffs_lowbd(filter_params_x, subpel_x_q4, coeffs_256);

      filt_256[0] = _mm256_loadu_si256((__m256i const *)filt1_global_avx2);
      filt_256[1] = _mm256_loadu_si256((__m256i const *)filt2_global_avx2);
      for (int i = 0; i < h; i += 2) {
        const __m256i data = _mm256_permute2x128_si256(
            _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)(&src_ptr[i * src_stride]))),
            _mm256_castsi128_si256(_mm_loadu_si128(
                (__m128i *)(&src_ptr[i * src_stride + src_stride]))),
            0x20);

        __m256i res_16b = convolve_lowbd_x_4tap(data, coeffs_256 + 1, filt_256);
        res_16b = sr_x_round_avx2(res_16b);

        __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

        const __m128i res_0 = _mm256_castsi256_si128(res_8b);
        const __m128i res_1 = _mm256_extracti128_si256(res_8b, 1);

        _mm_storel_epi64((__m128i *)&dst[i * dst_stride], res_0);
        _mm_storel_epi64((__m128i *)&dst[i * dst_stride + dst_stride], res_1);
      }
    } else {
      assert(!(w % 16));
      // TODO(chiyotsai@google.com): Reuse the old SIMD code here. Need to
      // rewrite this for better performance later.
      __m256i filt_256[2];
      prepare_coeffs_lowbd(filter_params_x, subpel_x_q4, coeffs_256);
      filt_256[0] = _mm256_loadu_si256((__m256i const *)filt1_global_avx2);
      filt_256[1] = _mm256_loadu_si256((__m256i const *)filt2_global_avx2);

      for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; j += 16) {
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 8 9 10 11 12 13 14 15 16 17
          // 18 19 20 21 22 23
          const __m256i data = _mm256_inserti128_si256(
              _mm256_loadu_si256((__m256i *)&src_ptr[(i * src_stride) + j]),
              _mm_loadu_si128((__m128i *)&src_ptr[(i * src_stride) + (j + 8)]),
              1);

          __m256i res_16b =
              convolve_lowbd_x_4tap(data, coeffs_256 + 1, filt_256);
          res_16b = sr_x_round_avx2(res_16b);

          /* rounding code */
          // 8 bit conversion and saturation to uint8
          __m256i res_8b = _mm256_packus_epi16(res_16b, res_16b);

          // Store values into the destination buffer
          // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
          res_8b = _mm256_permute4x64_epi64(res_8b, 216);
          __m128i res = _mm256_castsi256_si128(res_8b);
          _mm_storeu_si128((__m128i *)&dst[i * dst_stride + j], res);
        }
      }
    }
  } else {
    __m256i filt_256[4];

    filt_256[0] = _mm256_loadu_si256((__m256i const *)filt1_global_avx2);
    filt_256[1] = _mm256_loadu_si256((__m256i const *)filt2_global_avx2);
    filt_256[2] = _mm256_loadu_si256((__m256i const *)filt3_global_avx2);

    if (horz_tap == 6) {
      // horz_filt as 6 tap
      const uint8_t *src_ptr = src - 2;

      prepare_half_coeffs_6tap_avx2(filter_params_x, subpel_x_q4, coeffs_256);

      if (w == 8) {
        do {
          const __m256i res = x_convolve_6tap_8x2_avx2(src_ptr, src_stride,
                                                       coeffs_256, filt_256);
          sr_x_round_store_8x2_avx2(res, dst, dst_stride);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        do {
          __m256i r[2];

          x_convolve_6tap_16x2_avx2(src_ptr, src_stride, coeffs_256, filt_256,
                                    r);
          sr_x_round_store_16x2_avx2(r, dst, dst_stride);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 32) {
        do {
          sr_x_6tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else if (w == 64) {
        do {
          sr_x_6tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          sr_x_6tap_32_avx2(src_ptr + 32, coeffs_256, filt_256, dst + 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else {
        assert(w == 128);

        do {
          sr_x_6tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          sr_x_6tap_32_avx2(src_ptr + 1 * 32, coeffs_256, filt_256,
                            dst + 1 * 32);
          sr_x_6tap_32_avx2(src_ptr + 2 * 32, coeffs_256, filt_256,
                            dst + 2 * 32);
          sr_x_6tap_32_avx2(src_ptr + 3 * 32, coeffs_256, filt_256,
                            dst + 3 * 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      }
    } else if (horz_tap == 8) {
      // horz_filt as 8 tap
      const uint8_t *src_ptr = src - 3;

      filt_256[3] = _mm256_loadu_si256((__m256i const *)filt4_global_avx2);

      prepare_half_coeffs_8tap_avx2(filter_params_x, subpel_x_q4, coeffs_256);

      if (w == 8) {
        do {
          const __m256i res = x_convolve_8tap_8x2_avx2(src_ptr, src_stride,
                                                       coeffs_256, filt_256);
          sr_x_round_store_8x2_avx2(res, dst, dst_stride);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 16) {
        do {
          __m256i r[2];

          x_convolve_8tap_16x2_avx2(src_ptr, src_stride, coeffs_256, filt_256,
                                    r);
          sr_x_round_store_16x2_avx2(r, dst, dst_stride);
          src_ptr += 2 * src_stride;
          dst += 2 * dst_stride;
          y -= 2;
        } while (y);
      } else if (w == 32) {
        do {
          sr_x_8tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else if (w == 64) {
        do {
          sr_x_8tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          sr_x_8tap_32_avx2(src_ptr + 32, coeffs_256, filt_256, dst + 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      } else {
        assert(w == 128);

        do {
          sr_x_8tap_32_avx2(src_ptr, coeffs_256, filt_256, dst);
          sr_x_8tap_32_avx2(src_ptr + 1 * 32, coeffs_256, filt_256,
                            dst + 1 * 32);
          sr_x_8tap_32_avx2(src_ptr + 2 * 32, coeffs_256, filt_256,
                            dst + 2 * 32);
          sr_x_8tap_32_avx2(src_ptr + 3 * 32, coeffs_256, filt_256,
                            dst + 3 * 32);
          src_ptr += src_stride;
          dst += dst_stride;
        } while (--y);
      }
    } else if (horz_tap == 12) {
      av1_convolve_x_sr_avx2_general(src, src_stride, dst, dst_stride, w, h,
                                     filter_params_x, subpel_x_q4, conv_params);
    }
  }
}
