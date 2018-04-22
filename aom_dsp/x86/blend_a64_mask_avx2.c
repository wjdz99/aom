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

#include <immintrin.h>  // AVX2

#include <assert.h>

#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/blend.h"

#include "./aom_dsp_rtcd.h"

void aom_highbd_blend_a64_d16_mask_avx2(
    uint8_t *dst_8, uint32_t dst_stride, const CONV_BUF_TYPE *src0,
    uint32_t src0_stride, const CONV_BUF_TYPE *src1, uint32_t src1_stride,
    const uint8_t *mask, uint32_t mask_stride, int h, int w, int subh, int subw,
    ConvolveParams *conv_params, const int bd) {
  const int offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
  const int round_offset = (1 << (offset_bits - conv_params->round_1)) +
                           (1 << (offset_bits - conv_params->round_1 - 1));
  const int round_bits =
      2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);

  assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
  assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

  assert(h >= 1);
  assert(w >= 1);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  // excerpt from clip_pixel_highbd()
  // set saturation_value to (1 << bd) - 1
  unsigned int saturation_value;
  switch (bd) {
    case 8:
    default: saturation_value = 255; break;
    case 10: saturation_value = 1023; break;
    case 12: saturation_value = 4095; break;
  }
  const int n16 = w / 16;
  const __m256i vround_offset = _mm256_set1_epi32(round_offset);
  const __m256i vsaturation_value = _mm256_set1_epi32(saturation_value);
  const __m256i vround_half = _mm256_set1_epi32(1 << (round_bits - 1));
  const __m256i v0 = _mm256_setzero_si256();
  const __m256i v64 = _mm256_set1_epi16(64);
  const __m128i vround_bits = _mm_set1_epi64x(round_bits);

  if (subw == 0 && subh == 0) {
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < n16; ++j) {
        __m256i m0 = _mm256_cvtepi8_epi16(
            _mm_loadu_si128((const __m128i *)&mask[j * 16]));
        __m256i s0 = _mm256_loadu_si256((const __m256i *)&src0[j * 16]);
        __m256i s1 = _mm256_loadu_si256((const __m256i *)&src1[j * 16]);
        __m256i m1 = _mm256_subs_epu8(v64, m0);
        __m256i slo = _mm256_unpacklo_epi16(s0, s1);
        __m256i shi = _mm256_unpackhi_epi16(s0, s1);
        __m256i mlo = _mm256_unpacklo_epi16(m0, m1);
        __m256i mhi = _mm256_unpackhi_epi16(m0, m1);
        s0 = _mm256_permute2f128_si256(slo, shi, 0x20);
        s1 = _mm256_permute2f128_si256(slo, shi, 0x31);
        m0 = _mm256_permute2f128_si256(mlo, mhi, 0x20);
        m1 = _mm256_permute2f128_si256(mlo, mhi, 0x31);
        s0 = _mm256_madd_epi16(s0, m0);
        s1 = _mm256_madd_epi16(s1, m1);
        s0 = _mm256_srai_epi32(s0, AOM_BLEND_A64_ROUND_BITS);
        s1 = _mm256_srai_epi32(s1, AOM_BLEND_A64_ROUND_BITS);
        s0 = _mm256_sub_epi32(s0, vround_offset);
        s1 = _mm256_sub_epi32(s1, vround_offset);
        s0 = _mm256_add_epi32(s0, vround_half);
        s1 = _mm256_add_epi32(s1, vround_half);
        s0 = _mm256_max_epi32(s0, v0);
        s1 = _mm256_max_epi32(s1, v0);
        s0 = _mm256_srl_epi32(s0, vround_bits);
        s1 = _mm256_srl_epi32(s1, vround_bits);
        s0 = _mm256_min_epi32(s0, vsaturation_value);
        s1 = _mm256_min_epi32(s1, vsaturation_value);
        s0 = _mm256_packus_epi32(s0, s1);
        s0 = _mm256_permute4x64_epi64(s0, _MM_SHUFFLE(3, 1, 2, 0));
        _mm256_storeu_si256((__m256i *)&dst[j * 16], s0);
      }
      for (int j = n16 * 16; j < w; ++j) {
        int32_t res;
        const int m = mask[j];
        res = ((m * src0[j] + (AOM_BLEND_A64_MAX_ALPHA - m) * src1[j]) >>
               AOM_BLEND_A64_ROUND_BITS);
        res -= round_offset;
        unsigned int v = negative_to_zero(ROUND_POWER_OF_TWO(res, round_bits));
        dst[j] = AOMMIN(v, saturation_value);
      }
      mask += mask_stride;
      src0 += src0_stride;
      src1 += src1_stride;
      dst += dst_stride;
    }
  } else if (subw == 1 && subh == 1) {
    const __m256i v1 = _mm256_set1_epi8(1);
    const __m256i v2 = _mm256_set1_epi16(2);
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < n16; ++j) {
        __m256i m0 = _mm256_adds_epu8(
            _mm256_loadu_si256((const __m256i *)&mask[j * 32]),
            _mm256_loadu_si256((const __m256i *)&mask[mask_stride + j * 32]));
        m0 = _mm256_maddubs_epi16(m0, v1);
        m0 = _mm256_add_epi16(m0, v2);
        m0 = _mm256_srai_epi16(m0, 2);
        __m256i s0 = _mm256_loadu_si256((const __m256i *)&src0[j * 16]);
        __m256i s1 = _mm256_loadu_si256((const __m256i *)&src1[j * 16]);
        __m256i m1 = _mm256_subs_epu8(v64, m0);
        __m256i slo = _mm256_unpacklo_epi16(s0, s1);
        __m256i shi = _mm256_unpackhi_epi16(s0, s1);
        __m256i mlo = _mm256_unpacklo_epi16(m0, m1);
        __m256i mhi = _mm256_unpackhi_epi16(m0, m1);
        s0 = _mm256_permute2f128_si256(slo, shi, 0x20);
        s1 = _mm256_permute2f128_si256(slo, shi, 0x31);
        m0 = _mm256_permute2f128_si256(mlo, mhi, 0x20);
        m1 = _mm256_permute2f128_si256(mlo, mhi, 0x31);
        s0 = _mm256_madd_epi16(s0, m0);
        s1 = _mm256_madd_epi16(s1, m1);
        s0 = _mm256_srai_epi32(s0, AOM_BLEND_A64_ROUND_BITS);
        s1 = _mm256_srai_epi32(s1, AOM_BLEND_A64_ROUND_BITS);
        s0 = _mm256_sub_epi32(s0, vround_offset);
        s1 = _mm256_sub_epi32(s1, vround_offset);
        s0 = _mm256_add_epi32(s0, vround_half);
        s1 = _mm256_add_epi32(s1, vround_half);
        s0 = _mm256_max_epi32(s0, v0);
        s1 = _mm256_max_epi32(s1, v0);
        s0 = _mm256_srl_epi32(s0, vround_bits);
        s1 = _mm256_srl_epi32(s1, vround_bits);
        s0 = _mm256_min_epi32(s0, vsaturation_value);
        s1 = _mm256_min_epi32(s1, vsaturation_value);
        s0 = _mm256_packus_epi32(s0, s1);
        s0 = _mm256_permute4x64_epi64(s0, _MM_SHUFFLE(3, 1, 2, 0));
        _mm256_storeu_si256((__m256i *)&dst[j * 16], s0);
      }
      for (int j = n16 * 16; j < w; ++j) {
        int32_t res;
        const int m = ROUND_POWER_OF_TWO(
            mask[2 * j] + mask[mask_stride + 2 * j] + mask[2 * j + 1] +
                mask[mask_stride + 2 * j + 1],
            2);
        res = (m * src0[j] + (AOM_BLEND_A64_MAX_ALPHA - m) * src1[j]) >>
              AOM_BLEND_A64_ROUND_BITS;
        res -= round_offset;
        unsigned int v = negative_to_zero(ROUND_POWER_OF_TWO(res, round_bits));
        dst[j] = AOMMIN(v, saturation_value);
      }
      mask += 2 * mask_stride;
      src0 += src0_stride;
      src1 += src1_stride;
      dst += dst_stride;
    }
  } else if (subw == 1 && subh == 0) {
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < w; ++j) {
        int32_t res;
        const int m = AOM_BLEND_AVG(mask[2 * j], mask[2 * j + 1]);
        res = (m * src0[j] + (AOM_BLEND_A64_MAX_ALPHA - m) * src1[j]) >>
              AOM_BLEND_A64_ROUND_BITS;
        res -= round_offset;
        unsigned int v = negative_to_zero(ROUND_POWER_OF_TWO(res, round_bits));
        dst[j] = AOMMIN(v, saturation_value);
      }
      mask += mask_stride;
      src0 += src0_stride;
      src1 += src1_stride;
      dst += dst_stride;
    }
  } else {
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < w; ++j) {
        int32_t res;
        const int m = AOM_BLEND_AVG(mask[j], mask[mask_stride + j]);
        res = (m * src0[j] + (AOM_BLEND_A64_MAX_ALPHA - m) * src1[j]) >>
              AOM_BLEND_A64_ROUND_BITS;
        res -= round_offset;
        unsigned int v = negative_to_zero(ROUND_POWER_OF_TWO(res, round_bits));
        dst[j] = AOMMIN(v, saturation_value);
      }
      mask += 2 * mask_stride;
      src0 += src0_stride;
      src1 += src1_stride;
      dst += dst_stride;
    }
  }
}
