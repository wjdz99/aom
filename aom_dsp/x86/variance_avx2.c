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

void aom_highbd_variance64_avx2(const uint8_t *a8, int a_stride,
                                const uint8_t *b8, int b_stride, int w, int h,
                                uint64_t *sse, int64_t *sum) {
  const uint16_t *a = CONVERT_TO_SHORTPTR(a8);
  const uint16_t *b = CONVERT_TO_SHORTPTR(b8);

  switch (w) {
    case 128:
    case 64:
    case 32:
    case 16: {
      const int n16 = w / 16;
      __m256i vsum = _mm256_setzero_si256();
      __m256i vsse = _mm256_setzero_si256();
      for (int i = 0; i < h; ++i) {
        __m256i vsum0;
        __m256i vsum1;
        __m256i vsse0;
        __m256i vsse1;
        {
          __m256i al =
              _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(a + 0)));
          __m256i ah =
              _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(a + 8)));
          __m256i bl =
              _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(b + 0)));
          __m256i bh =
              _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(b + 8)));
          vsum0 = _mm256_sub_epi32(al, bl);
          vsum1 = _mm256_sub_epi32(ah, bh);
          vsse0 = _mm256_mullo_epi32(vsum0, vsum0);
          vsse1 = _mm256_mullo_epi32(vsum1, vsum1);
        }
        for (int j = 1; j < n16; ++j) {
          __m256i al = _mm256_cvtepu16_epi32(
              _mm_loadu_si128((const __m128i *)(a + j * 16 + 0)));
          __m256i ah = _mm256_cvtepu16_epi32(
              _mm_loadu_si128((const __m128i *)(a + j * 16 + 8)));
          __m256i bl = _mm256_cvtepu16_epi32(
              _mm_loadu_si128((const __m128i *)(b + j * 16 + 0)));
          __m256i bh = _mm256_cvtepu16_epi32(
              _mm_loadu_si128((const __m128i *)(b + j * 16 + 8)));
          __m256i dl = _mm256_sub_epi32(al, bl);
          __m256i dh = _mm256_sub_epi32(ah, bh);
          __m256i dl2 = _mm256_mullo_epi32(dl, dl);
          __m256i dh2 = _mm256_mullo_epi32(dh, dh);
          vsum0 = _mm256_add_epi32(vsum0, dl);
          vsum1 = _mm256_add_epi32(vsum1, dh);
          vsse0 = _mm256_add_epi32(vsse0, dl2);
          vsse1 = _mm256_add_epi32(vsse1, dh2);
        }
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_castsi256_si128(vsum0)));
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_extracti128_si256(vsum0, 1)));
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_castsi256_si128(vsum1)));
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_extracti128_si256(vsum1, 1)));
        vsse = _mm256_add_epi64(
            vsse, _mm256_cvtepu32_epi64(_mm256_castsi256_si128(vsse0)));
        vsse = _mm256_add_epi64(
            vsse, _mm256_cvtepu32_epi64(_mm256_extracti128_si256(vsse0, 1)));
        vsse = _mm256_add_epi64(
            vsse, _mm256_cvtepu32_epi64(_mm256_castsi256_si128(vsse1)));
        vsse = _mm256_add_epi64(
            vsse, _mm256_cvtepu32_epi64(_mm256_extracti128_si256(vsse1, 1)));
        a += a_stride;
        b += b_stride;
      }
      __m128i xsum = _mm_add_epi64(_mm256_castsi256_si128(vsum),
                                   _mm256_extracti128_si256(vsum, 1));
      __m128i xsse = _mm_add_epi64(_mm256_castsi256_si128(vsse),
                                   _mm256_extracti128_si256(vsse, 1));
      xsum = _mm_add_epi64(xsum, _mm_srli_si128(xsum, 8));
      xsse = _mm_add_epi64(xsse, _mm_srli_si128(xsse, 8));
      _mm_storel_epi64((__m128i *)sum, xsum);
      _mm_storel_epi64((__m128i *)sse, xsse);
    } break;
    case 8: {
      __m256i vsum = _mm256_setzero_si256();
      __m256i vsse = _mm256_setzero_si256();
      for (int i = 0; i < h; ++i) {
        __m256i va = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)a));
        __m256i vb = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)b));
        __m256i vd = _mm256_sub_epi32(va, vb);
        __m256i vd2 = _mm256_mullo_epi32(vd, vd);
        __m256i vdl = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(vd));
        __m256i vdh = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(vd, 1));
        __m256i vd2l = _mm256_cvtepu32_epi64(_mm256_castsi256_si128(vd2));
        __m256i vd2h = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(vd2, 1));
        vsum = _mm256_add_epi64(vsum, vdl);
        vsum = _mm256_add_epi64(vsum, vdh);
        vsse = _mm256_add_epi64(vsse, vd2l);
        vsse = _mm256_add_epi64(vsse, vd2h);
        a += a_stride;
        b += b_stride;
      }
      __m128i vsum2 = _mm_add_epi64(_mm256_castsi256_si128(vsum),
                                    _mm256_extracti128_si256(vsum, 1));
      __m128i vsse2 = _mm_add_epi64(_mm256_castsi256_si128(vsse),
                                    _mm256_extracti128_si256(vsse, 1));
      vsum2 = _mm_add_epi64(vsum2, _mm_srli_si128(vsum2, 8));
      vsse2 = _mm_add_epi64(vsse2, _mm_srli_si128(vsse2, 8));
      _mm_storel_epi64((__m128i *)sum, vsum2);
      _mm_storel_epi64((__m128i *)sse, vsse2);
    } break;
    case 4: {
      __m256i vsum = _mm256_setzero_si256();
      __m256i vsse = _mm256_setzero_si256();
      for (int i = 0; i < h; ++i) {
        __m128i va = _mm_cvtepu16_epi32(_mm_loadl_epi64((const __m128i *)a));
        __m128i vb = _mm_cvtepu16_epi32(_mm_loadl_epi64((const __m128i *)b));
        __m128i vd = _mm_sub_epi32(va, vb);
        __m128i vd2 = _mm_mullo_epi32(vd, vd);
        vsum = _mm256_add_epi64(vsum, _mm256_cvtepi32_epi64(vd));
        vsse = _mm256_add_epi64(vsse, _mm256_cvtepi32_epi64(vd2));
        a += a_stride;
        b += b_stride;
      }
      __m128i vsum2 = _mm_add_epi64(_mm256_castsi256_si128(vsum),
                                    _mm256_extracti128_si256(vsum, 1));
      __m128i vsse2 = _mm_add_epi64(_mm256_castsi256_si128(vsse),
                                    _mm256_extracti128_si256(vsse, 1));
      vsum2 = _mm_add_epi64(vsum2, _mm_srli_si128(vsum2, 8));
      vsse2 = _mm_add_epi64(vsse2, _mm_srli_si128(vsse2, 8));
      _mm_storel_epi64((__m128i *)sum, vsum2);
      _mm_storel_epi64((__m128i *)sse, vsse2);
    } break;
    case 2: {
      int64_t tsum = 0;
      uint64_t tsse = 0;
      for (int i = 0; i < h; ++i) {
        for (int j = 0; j < 2; ++j) {
          const int diff0 = a[0] - b[0];
          const int diff1 = a[1] - b[1];
          tsum += diff0;
          tsum += diff1;
          tsse += (uint32_t)(diff0 * diff0);
          tsse += (uint32_t)(diff1 * diff1);
        }
        a += a_stride;
        b += b_stride;
      }
      *sum = tsum;
      *sse = tsse;
    } break;
    default: {
      int64_t tsum = 0;
      uint64_t tsse = 0;
      for (int i = 0; i < h; ++i) {
        int32_t lsum = 0;
        for (int j = 0; j < w; ++j) {
          const int diff = a[j] - b[j];
          lsum += diff;
          tsse += (uint32_t)(diff * diff);
        }
        tsum += lsum;
        a += a_stride;
        b += b_stride;
      }
      *sum = tsum;
      *sse = tsse;
    } break;
  }
}
