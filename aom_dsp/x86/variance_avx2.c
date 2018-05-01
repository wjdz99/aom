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

#include <immintrin.h>
#include "./aom_dsp_rtcd.h"
#include "aom_dsp/x86/masked_variance_intrin_ssse3.h"

typedef void (*get_var_avx2)(const uint8_t *src, int src_stride,
                             const uint8_t *ref, int ref_stride,
                             unsigned int *sse, int *sum);

void aom_get32x32var_avx2(const uint8_t *src, int src_stride,
                          const uint8_t *ref, int ref_stride, unsigned int *sse,
                          int *sum);

static void variance_avx2(const uint8_t *src, int src_stride,
                          const uint8_t *ref, int ref_stride, int w, int h,
                          unsigned int *sse, int *sum, get_var_avx2 var_fn,
                          int block_size) {
  int i, j;

  *sse = 0;
  *sum = 0;

  for (i = 0; i < h; i += 16) {
    for (j = 0; j < w; j += block_size) {
      unsigned int sse0;
      int sum0;
      var_fn(&src[src_stride * i + j], src_stride, &ref[ref_stride * i + j],
             ref_stride, &sse0, &sum0);
      *sse += sse0;
      *sum += sum0;
    }
  }
}

unsigned int aom_variance16x16_avx2(const uint8_t *src, int src_stride,
                                    const uint8_t *ref, int ref_stride,
                                    unsigned int *sse) {
  int sum;
  unsigned int variance;
  aom_get16x16var_avx2(src, src_stride, ref, ref_stride, sse, &sum);
  variance = *sse - (((uint32_t)((int64_t)sum * sum)) >> 8);
  _mm256_zeroupper();
  return variance;
}

unsigned int aom_mse16x16_avx2(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride,
                               unsigned int *sse) {
  int sum;
  aom_get16x16var_avx2(src, src_stride, ref, ref_stride, sse, &sum);
  _mm256_zeroupper();
  return *sse;
}

#define AOM_VAR_AVX2(bw, bh, w, bits)                                         \
  unsigned int aom_variance##bw##x##bh##_avx2(                                \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, \
      unsigned int *sse) {                                                    \
    int sum;                                                                  \
    unsigned int variance;                                                    \
    variance_avx2(src, src_stride, ref, ref_stride, bw, bh, sse, &sum,        \
                  aom_get##w##x##w##var_avx2, w);                             \
    variance = *sse - (uint32_t)(((int64_t)sum * sum) >> bits);               \
    _mm256_zeroupper();                                                       \
    return variance;                                                          \
  }

AOM_VAR_AVX2(32, 16, 16, 9);
AOM_VAR_AVX2(32, 32, 32, 10);
AOM_VAR_AVX2(64, 64, 32, 12);
AOM_VAR_AVX2(64, 32, 32, 11);
AOM_VAR_AVX2(128, 128, 32, 14);
AOM_VAR_AVX2(128, 64, 32, 13);
AOM_VAR_AVX2(64, 128, 32, 13);

unsigned int aom_sub_pixel_variance32xh_avx2(const uint8_t *src, int src_stride,
                                             int x_offset, int y_offset,
                                             const uint8_t *dst, int dst_stride,
                                             int height, unsigned int *sse);

unsigned int aom_sub_pixel_avg_variance32xh_avx2(
    const uint8_t *src, int src_stride, int x_offset, int y_offset,
    const uint8_t *dst, int dst_stride, const uint8_t *sec, int sec_stride,
    int height, unsigned int *sseptr);

unsigned int aom_sub_pixel_variance64x64_avx2(const uint8_t *src,
                                              int src_stride, int x_offset,
                                              int y_offset, const uint8_t *dst,
                                              int dst_stride,
                                              unsigned int *sse) {
  unsigned int sse1;
  const int se1 = aom_sub_pixel_variance32xh_avx2(
      src, src_stride, x_offset, y_offset, dst, dst_stride, 64, &sse1);
  unsigned int sse2;
  const int se2 =
      aom_sub_pixel_variance32xh_avx2(src + 32, src_stride, x_offset, y_offset,
                                      dst + 32, dst_stride, 64, &sse2);
  const int se = se1 + se2;
  unsigned int variance;
  *sse = sse1 + sse2;

  variance = *sse - (uint32_t)(((int64_t)se * se) >> 12);
  _mm256_zeroupper();
  return variance;
}

unsigned int aom_sub_pixel_variance32x32_avx2(const uint8_t *src,
                                              int src_stride, int x_offset,
                                              int y_offset, const uint8_t *dst,
                                              int dst_stride,
                                              unsigned int *sse) {
  const int se = aom_sub_pixel_variance32xh_avx2(
      src, src_stride, x_offset, y_offset, dst, dst_stride, 32, sse);

  const unsigned int variance = *sse - (uint32_t)(((int64_t)se * se) >> 10);
  _mm256_zeroupper();
  return variance;
}

unsigned int aom_sub_pixel_avg_variance64x64_avx2(
    const uint8_t *src, int src_stride, int x_offset, int y_offset,
    const uint8_t *dst, int dst_stride, unsigned int *sse, const uint8_t *sec) {
  unsigned int sse1;
  const int se1 = aom_sub_pixel_avg_variance32xh_avx2(
      src, src_stride, x_offset, y_offset, dst, dst_stride, sec, 64, 64, &sse1);
  unsigned int sse2;
  const int se2 = aom_sub_pixel_avg_variance32xh_avx2(
      src + 32, src_stride, x_offset, y_offset, dst + 32, dst_stride, sec + 32,
      64, 64, &sse2);
  const int se = se1 + se2;
  unsigned int variance;

  *sse = sse1 + sse2;

  variance = *sse - (uint32_t)(((int64_t)se * se) >> 12);
  _mm256_zeroupper();
  return variance;
}

unsigned int aom_sub_pixel_avg_variance32x32_avx2(
    const uint8_t *src, int src_stride, int x_offset, int y_offset,
    const uint8_t *dst, int dst_stride, unsigned int *sse, const uint8_t *sec) {
  // Process 32 elements in parallel.
  const int se = aom_sub_pixel_avg_variance32xh_avx2(
      src, src_stride, x_offset, y_offset, dst, dst_stride, sec, 32, 32, sse);

  const unsigned int variance = *sse - (uint32_t)(((int64_t)se * se) >> 10);
  _mm256_zeroupper();
  return variance;
}

static INLINE __m256i mm256_loadu2(const uint8_t *p0, const uint8_t *p1) {
  const __m256i d =
      _mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)p1));
  return _mm256_insertf128_si256(d, _mm_loadu_si128((const __m128i *)p0), 1);
}

static INLINE void comp_mask_pred_line_avx2(const __m256i s0, const __m256i s1,
                                            const __m256i a,
                                            uint8_t *comp_pred) {
  const __m256i alpha_max = _mm256_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
  const int16_t round_bits = 15 - AOM_BLEND_A64_ROUND_BITS;
  const __m256i round_offset = _mm256_set1_epi16(1 << (round_bits));

  const __m256i ma = _mm256_sub_epi8(alpha_max, a);

  const __m256i ssAL = _mm256_unpacklo_epi8(s0, s1);
  const __m256i aaAL = _mm256_unpacklo_epi8(a, ma);
  const __m256i ssAH = _mm256_unpackhi_epi8(s0, s1);
  const __m256i aaAH = _mm256_unpackhi_epi8(a, ma);

  const __m256i blendAL = _mm256_maddubs_epi16(ssAL, aaAL);
  const __m256i blendAH = _mm256_maddubs_epi16(ssAH, aaAH);
  const __m256i roundAL = _mm256_mulhrs_epi16(blendAL, round_offset);
  const __m256i roundAH = _mm256_mulhrs_epi16(blendAH, round_offset);

  const __m256i roundA = _mm256_packus_epi16(roundAL, roundAH);
  _mm256_storeu_si256((__m256i *)(comp_pred), roundA);
}

void aom_comp_mask_pred_avx2(uint8_t *comp_pred, const uint8_t *pred, int width,
                             int height, const uint8_t *ref, int ref_stride,
                             const uint8_t *mask, int mask_stride,
                             int invert_mask) {
  int i = 0;
  const uint8_t *src0 = invert_mask ? pred : ref;
  const uint8_t *src1 = invert_mask ? ref : pred;
  const int stride0 = invert_mask ? width : ref_stride;
  const int stride1 = invert_mask ? ref_stride : width;
  if (width == 8) {
    comp_mask_pred_8_ssse3(comp_pred, height, src0, stride0, src1, stride1,
                           mask, mask_stride);
  } else if (width == 16) {
    do {
      const __m256i sA0 = mm256_loadu2(src0 + stride0, src0);
      const __m256i sA1 = mm256_loadu2(src1 + stride1, src1);
      const __m256i aA = mm256_loadu2(mask + mask_stride, mask);
      src0 += (stride0 << 1);
      src1 += (stride1 << 1);
      mask += (mask_stride << 1);
      const __m256i sB0 = mm256_loadu2(src0 + stride0, src0);
      const __m256i sB1 = mm256_loadu2(src1 + stride1, src1);
      const __m256i aB = mm256_loadu2(mask + mask_stride, mask);
      src0 += (stride0 << 1);
      src1 += (stride1 << 1);
      mask += (mask_stride << 1);
      // comp_pred's stride == width == 16
      comp_mask_pred_line_avx2(sA0, sA1, aA, comp_pred);
      comp_mask_pred_line_avx2(sB0, sB1, aB, comp_pred + 32);
      comp_pred += (16 << 2);
      i += 4;
    } while (i < height);
  } else {  // for width == 32
    do {
      const __m256i sA0 = _mm256_lddqu_si256((const __m256i *)(src0));
      const __m256i sA1 = _mm256_lddqu_si256((const __m256i *)(src1));
      const __m256i aA = _mm256_lddqu_si256((const __m256i *)(mask));

      const __m256i sB0 = _mm256_lddqu_si256((const __m256i *)(src0 + stride0));
      const __m256i sB1 = _mm256_lddqu_si256((const __m256i *)(src1 + stride1));
      const __m256i aB =
          _mm256_lddqu_si256((const __m256i *)(mask + mask_stride));

      comp_mask_pred_line_avx2(sA0, sA1, aA, comp_pred);
      comp_mask_pred_line_avx2(sB0, sB1, aB, comp_pred + 32);
      comp_pred += (32 << 1);

      src0 += (stride0 << 1);
      src1 += (stride1 << 1);
      mask += (mask_stride << 1);
      i += 2;
    } while (i < height);
  }
}

void aom_highbd_comp_mask_upsampled_pred_avx2(
    uint16_t *comp_pred, const uint8_t *pred8, int width, int height,
    int subpel_x_q3, int subpel_y_q3, const uint8_t *ref8, int ref_stride,
    const uint8_t *mask, int mask_stride, int invert_mask, int bd) {
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
  aom_highbd_upsampled_pred(comp_pred, width, height, subpel_x_q3, subpel_y_q3,
                            ref8, ref_stride, bd);
  const __m256i v64 = _mm256_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);
  const __m256i vhalf = _mm256_set1_epi32(1 << (AOM_BLEND_A64_ROUND_BITS - 1));
  uint16_t *pv0;
  uint16_t *pv1;
  if (invert_mask) {
    pv0 = pred;
    pv1 = comp_pred;
  } else {
    pv0 = comp_pred;
    pv1 = pred;
  }
  switch (width) {
    case 128:
    case 64:
    case 32: {
      const int n32 = width / 32;
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < n32; ++j) {
          __m256i vmaskA = _mm256_cvtepu8_epi16(
              _mm_loadu_si128((const __m128i *)(mask + j * 32 + 0)));
          __m256i vmaskB = _mm256_cvtepu8_epi16(
              _mm_loadu_si128((const __m128i *)(mask + j * 32 + 16)));
          __m256i v0A = _mm256_loadu_si256((const __m256i *)(pv0 + j * 32 + 0));
          __m256i v0B =
              _mm256_loadu_si256((const __m256i *)(pv0 + j * 32 + 16));
          __m256i v1A = _mm256_loadu_si256((const __m256i *)(pv1 + j * 32 + 0));
          __m256i v1B =
              _mm256_loadu_si256((const __m256i *)(pv1 + j * 32 + 16));
          __m256i vmaskA2 = _mm256_sub_epi16(v64, vmaskA);
          __m256i vmaskB2 = _mm256_sub_epi16(v64, vmaskB);
          __m256i vmaskAlo = _mm256_unpacklo_epi16(vmaskA, vmaskA2);
          __m256i vmaskAhi = _mm256_unpackhi_epi16(vmaskA, vmaskA2);
          __m256i vmaskBlo = _mm256_unpacklo_epi16(vmaskB, vmaskB2);
          __m256i vmaskBhi = _mm256_unpackhi_epi16(vmaskB, vmaskB2);
          __m256i v01Alo = _mm256_unpacklo_epi16(v0A, v1A);
          __m256i v01Ahi = _mm256_unpackhi_epi16(v0A, v1A);
          __m256i v01Blo = _mm256_unpacklo_epi16(v0B, v1B);
          __m256i v01Bhi = _mm256_unpackhi_epi16(v0B, v1B);
          vmaskAlo = _mm256_madd_epi16(vmaskAlo, v01Alo);
          vmaskAhi = _mm256_madd_epi16(vmaskAhi, v01Ahi);
          vmaskBlo = _mm256_madd_epi16(vmaskBlo, v01Blo);
          vmaskBhi = _mm256_madd_epi16(vmaskBhi, v01Bhi);
          vmaskAlo = _mm256_add_epi32(vmaskAlo, vhalf);
          vmaskAhi = _mm256_add_epi32(vmaskAhi, vhalf);
          vmaskBlo = _mm256_add_epi32(vmaskBlo, vhalf);
          vmaskBhi = _mm256_add_epi32(vmaskBhi, vhalf);
          vmaskAlo = _mm256_srai_epi32(vmaskAlo, AOM_BLEND_A64_ROUND_BITS);
          vmaskAhi = _mm256_srai_epi32(vmaskAhi, AOM_BLEND_A64_ROUND_BITS);
          vmaskBlo = _mm256_srai_epi32(vmaskBlo, AOM_BLEND_A64_ROUND_BITS);
          vmaskBhi = _mm256_srai_epi32(vmaskBhi, AOM_BLEND_A64_ROUND_BITS);
          vmaskAlo = _mm256_packs_epi32(vmaskAlo, vmaskAhi);
          vmaskBlo = _mm256_packs_epi32(vmaskBlo, vmaskBhi);
          _mm256_storeu_si256((__m256i *)(comp_pred + j * 32 + 0), vmaskAlo);
          _mm256_storeu_si256((__m256i *)(comp_pred + j * 32 + 16), vmaskBlo);
        }
        comp_pred += width;
        pv0 += width;
        pv1 += width;
        mask += mask_stride;
      }
    } break;
    case 16: {
      for (int i = 0; i < height; ++i) {
        __m256i vmask =
            _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)mask));
        __m256i v0 = _mm256_loadu_si256((const __m256i *)pv0);
        __m256i v1 = _mm256_loadu_si256((const __m256i *)pv1);
        __m256i vmask2 = _mm256_sub_epi16(v64, vmask);
        __m256i vmasklo = _mm256_unpacklo_epi16(vmask, vmask2);
        __m256i vmaskhi = _mm256_unpackhi_epi16(vmask, vmask2);
        __m256i v01lo = _mm256_unpacklo_epi16(v0, v1);
        __m256i v01hi = _mm256_unpackhi_epi16(v0, v1);
        vmasklo = _mm256_madd_epi16(vmasklo, v01lo);
        vmaskhi = _mm256_madd_epi16(vmaskhi, v01hi);
        vmasklo = _mm256_add_epi32(vmasklo, vhalf);
        vmaskhi = _mm256_add_epi32(vmaskhi, vhalf);
        vmasklo = _mm256_srai_epi32(vmasklo, AOM_BLEND_A64_ROUND_BITS);
        vmaskhi = _mm256_srai_epi32(vmaskhi, AOM_BLEND_A64_ROUND_BITS);
        vmasklo = _mm256_packs_epi32(vmasklo, vmaskhi);
        _mm256_storeu_si256((__m256i *)comp_pred, vmasklo);
        comp_pred += width;
        pv0 += width;
        pv1 += width;
        mask += mask_stride;
      }
    } break;
    case 8:
      for (int i = 0; i < height; ++i) {
        __m128i vmask =
            _mm_cvtepu8_epi16(_mm_loadl_epi64((const __m128i *)mask));
        __m128i v0 = _mm_loadu_si128((const __m128i *)(pv0));
        __m128i v1 = _mm_loadu_si128((const __m128i *)(pv1));
        __m128i vmask2 = _mm_sub_epi16(_mm256_castsi256_si128(v64), vmask);
        __m128i vmasklo = _mm_unpacklo_epi16(vmask, vmask2);
        __m128i vmaskhi = _mm_unpackhi_epi16(vmask, vmask2);
        __m128i v01lo = _mm_unpacklo_epi16(v0, v1);
        __m128i v01hi = _mm_unpackhi_epi16(v0, v1);
        vmasklo = _mm_madd_epi16(vmasklo, v01lo);
        vmaskhi = _mm_madd_epi16(vmaskhi, v01hi);
        vmasklo = _mm_add_epi32(vmasklo, _mm256_castsi256_si128(vhalf));
        vmaskhi = _mm_add_epi32(vmaskhi, _mm256_castsi256_si128(vhalf));
        vmasklo = _mm_srai_epi32(vmasklo, AOM_BLEND_A64_ROUND_BITS);
        vmaskhi = _mm_srai_epi32(vmaskhi, AOM_BLEND_A64_ROUND_BITS);
        vmasklo = _mm_packs_epi32(vmasklo, vmaskhi);
        _mm_storeu_si128((__m128i *)comp_pred, vmasklo);
        comp_pred += width;
        pv0 += width;
        pv1 += width;
        mask += mask_stride;
      }
      break;
    case 4:
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
          comp_pred[j] = AOM_BLEND_A64(mask[j], pv0[j], pv1[j]);
        }
        comp_pred += width;
        pv0 += width;
        pv1 += width;
        mask += mask_stride;
      }
      break;
  }
}
