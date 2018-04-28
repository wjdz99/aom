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
    case 32: {
      __m256i vsum = _mm256_setzero_si256();
      __m256i vsse = _mm256_setzero_si256();
      const __m256i maskFFFF = _mm256_set1_epi32(0xFFFF);
      const __m256i maskFFFFFFFF = _mm256_set1_epi64x(0xFFFFFFFF);
      const int n32 = w / 32;
      for (int i = 0; i < h; ++i) {
        __m256i vsum32 = _mm256_setzero_si256();
        const __m256i *va = (const __m256i *)a;
        const __m256i *vb = (const __m256i *)b;
        for (int j = 0; j < n32; ++j) {
          __m256i a0 = _mm256_loadu_si256(&va[2 * j + 0]);
          __m256i a1 = _mm256_loadu_si256(&va[2 * j + 1]);
          __m256i b0 = _mm256_loadu_si256(&vb[2 * j + 0]);
          __m256i b1 = _mm256_loadu_si256(&vb[2 * j + 1]);
          __m256i a0l = _mm256_and_si256(a0, maskFFFF);
          __m256i a1l = _mm256_and_si256(a1, maskFFFF);
          __m256i b0l = _mm256_and_si256(b0, maskFFFF);
          __m256i b1l = _mm256_and_si256(b1, maskFFFF);
          __m256i a0h = _mm256_srli_epi32(a0, 16);
          __m256i a1h = _mm256_srli_epi32(a1, 16);
          __m256i b0h = _mm256_srli_epi32(b0, 16);
          __m256i b1h = _mm256_srli_epi32(b1, 16);

          a0l = _mm256_sub_epi32(a0l, b0l);
          a1l = _mm256_sub_epi32(a1l, b1l);
          a0h = _mm256_sub_epi32(a0h, b0h);
          a1h = _mm256_sub_epi32(a1h, b1h);

          b0l = _mm256_add_epi32(a0l, a1l);
          b1l = _mm256_add_epi32(a0h, a1h);
          b0l = _mm256_add_epi32(b0l, b1l);
          vsum32 = _mm256_add_epi32(vsum32, b0l);

          a0l = _mm256_mullo_epi32(a0l, a0l);
          a1l = _mm256_mullo_epi32(a1l, a1l);
          a0h = _mm256_mullo_epi32(a0h, a0h);
          a1h = _mm256_mullo_epi32(a1h, a1h);

          b0l = _mm256_and_si256(a0l, maskFFFFFFFF);
          b1l = _mm256_and_si256(a1l, maskFFFFFFFF);
          b0h = _mm256_and_si256(a0h, maskFFFFFFFF);
          b1h = _mm256_and_si256(a1h, maskFFFFFFFF);

          a0l = _mm256_srli_epi64(a0l, 32);
          a1l = _mm256_srli_epi64(a1l, 32);
          a0h = _mm256_srli_epi64(a0h, 32);
          a1h = _mm256_srli_epi64(a1h, 32);

          b0l = _mm256_add_epi64(b0l, a0l);
          b1l = _mm256_add_epi64(b1l, a1l);
          b0h = _mm256_add_epi64(b0h, a0h);
          b1h = _mm256_add_epi64(b1h, a1h);

          b0l = _mm256_add_epi64(b0l, b1l);
          b0h = _mm256_add_epi64(b0h, b1h);
          b0l = _mm256_add_epi64(b0l, b0h);
          vsse = _mm256_add_epi64(vsse, b0l);
        }
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_castsi256_si128(vsum32)));
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_extracti128_si256(vsum32, 1)));
        a += a_stride;
        b += b_stride;
      }
      __m128i xsum = _mm_add_epi64(_mm256_castsi256_si128(vsum),
                                   _mm256_extracti128_si256(vsum, 1));
      __m128i xsse = _mm_add_epi64(_mm256_castsi256_si128(vsse),
                                   _mm256_extracti128_si256(vsse, 1));
      xsum = _mm_add_epi64(xsum, _mm_bsrli_si128(xsum, 8));
      xsse = _mm_add_epi64(xsse, _mm_bsrli_si128(xsse, 8));
      _mm_storel_epi64((__m128i *)sum, xsum);
      _mm_storel_epi64((__m128i *)sse, xsse);
    } break;
    case 16: {
      __m256i vsum = _mm256_setzero_si256();
      __m256i vsse = _mm256_setzero_si256();
      __m256i maskFFFF = _mm256_set1_epi32(0xFFFF);
      __m256i maskFFFFFFFF = _mm256_set1_epi64x(0xFFFFFFFF);
      for (int i = 0; i < h; ++i) {
        __m256i vsum32 = _mm256_setzero_si256();

        __m256i a0 = _mm256_loadu_si256((const __m256i *)a);
        __m256i b0 = _mm256_loadu_si256((const __m256i *)b);
        __m256i a0l = _mm256_and_si256(a0, maskFFFF);
        __m256i b0l = _mm256_and_si256(b0, maskFFFF);
        __m256i a0h = _mm256_srli_epi32(a0, 16);
        __m256i b0h = _mm256_srli_epi32(b0, 16);

        a0l = _mm256_sub_epi32(a0l, b0l);
        a0h = _mm256_sub_epi32(a0h, b0h);

        b0l = _mm256_add_epi32(a0l, a0h);
        vsum32 = _mm256_add_epi32(vsum32, b0l);

        a0l = _mm256_mullo_epi32(a0l, a0l);
        a0h = _mm256_mullo_epi32(a0h, a0h);

        b0l = _mm256_and_si256(a0l, maskFFFFFFFF);
        b0h = _mm256_and_si256(a0h, maskFFFFFFFF);

        a0l = _mm256_srli_epi64(a0l, 32);
        a0h = _mm256_srli_epi64(a0h, 32);

        b0l = _mm256_add_epi64(b0l, a0l);
        b0h = _mm256_add_epi64(b0h, a0h);

        b0l = _mm256_add_epi64(b0l, b0h);
        vsse = _mm256_add_epi64(vsse, b0l);

        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_castsi256_si128(vsum32)));
        vsum = _mm256_add_epi64(
            vsum, _mm256_cvtepi32_epi64(_mm256_extracti128_si256(vsum32, 1)));
        a += a_stride;
        b += b_stride;
      }
      __m128i xsum = _mm_add_epi64(_mm256_castsi256_si128(vsum),
                                   _mm256_extracti128_si256(vsum, 1));
      __m128i xsse = _mm_add_epi64(_mm256_castsi256_si128(vsse),
                                   _mm256_extracti128_si256(vsse, 1));
      xsum = _mm_add_epi64(xsum, _mm_bsrli_si128(xsum, 8));
      xsse = _mm_add_epi64(xsse, _mm_bsrli_si128(xsse, 8));
      _mm_storel_epi64((__m128i *)sum, xsum);
      _mm_storel_epi64((__m128i *)sse, xsse);
    } break;
    case 8: {
      __m128i vsum = _mm_setzero_si128();
      __m128i vsse = _mm_setzero_si128();
      const __m128i maskFFFF = _mm_set1_epi32(0xFFFF);
      const __m128i maskFFFFFFFF = _mm_set1_epi64x(0xFFFFFFFF);
      for (int i = 0; i < h; ++i) {
        __m128i vsum32 = _mm_setzero_si128();

        __m128i a0 = _mm_loadu_si128((const __m128i *)a);
        __m128i b0 = _mm_loadu_si128((const __m128i *)b);
        __m128i a0l = _mm_and_si128(a0, maskFFFF);
        __m128i b0l = _mm_and_si128(b0, maskFFFF);
        __m128i a0h = _mm_srli_epi32(a0, 16);
        __m128i b0h = _mm_srli_epi32(b0, 16);

        a0l = _mm_sub_epi32(a0l, b0l);
        a0h = _mm_sub_epi32(a0h, b0h);

        b0l = _mm_add_epi32(a0l, a0h);
        vsum32 = _mm_add_epi32(vsum32, b0l);

        a0l = _mm_mullo_epi32(a0l, a0l);
        a0h = _mm_mullo_epi32(a0h, a0h);

        b0l = _mm_and_si128(a0l, maskFFFFFFFF);
        b0h = _mm_and_si128(a0h, maskFFFFFFFF);

        a0l = _mm_srli_epi64(a0l, 32);
        a0h = _mm_srli_epi64(a0h, 32);

        b0l = _mm_add_epi64(b0l, a0l);
        b0h = _mm_add_epi64(b0h, a0h);

        b0l = _mm_add_epi64(b0l, b0h);
        vsse = _mm_add_epi64(vsse, b0l);

        vsum = _mm_add_epi64(vsum, _mm_cvtepi32_epi64(vsum32));
        a += a_stride;
        b += b_stride;
      }
      vsum = _mm_add_epi64(vsum, _mm_bsrli_si128(vsum, 8));
      vsse = _mm_add_epi64(vsse, _mm_bsrli_si128(vsse, 8));
      _mm_storel_epi64((__m128i *)sum, vsum);
      _mm_storel_epi64((__m128i *)sse, vsse);
    } break;
    case 4: {
      __m128i vsum = _mm_setzero_si128();
      __m128i vsse = _mm_setzero_si128();
      const __m128i maskFFFF = _mm_set1_epi32(0xFFFF);
      const __m128i maskFFFFFFFF = _mm_set1_epi64x(0xFFFFFFFF);
      for (int i = 0; i < h; ++i) {
        __m128i vsum32 = _mm_setzero_si128();

        __m128i a0 = _mm_loadu_si64((const __m128i *)a);
        __m128i b0 = _mm_loadu_si64((const __m128i *)b);
        __m128i a0l = _mm_and_si128(a0, maskFFFF);
        __m128i b0l = _mm_and_si128(b0, maskFFFF);
        __m128i a0h = _mm_srli_epi32(a0, 16);
        __m128i b0h = _mm_srli_epi32(b0, 16);

        a0l = _mm_sub_epi32(a0l, b0l);
        a0h = _mm_sub_epi32(a0h, b0h);

        b0l = _mm_add_epi32(a0l, a0h);
        vsum32 = _mm_add_epi32(vsum32, b0l);

        a0l = _mm_mullo_epi32(a0l, a0l);
        a0h = _mm_mullo_epi32(a0h, a0h);

        b0l = _mm_and_si128(a0l, maskFFFFFFFF);
        b0h = _mm_and_si128(a0h, maskFFFFFFFF);

        a0l = _mm_srli_epi64(a0l, 32);
        a0h = _mm_srli_epi64(a0h, 32);

        b0l = _mm_add_epi64(b0l, a0l);
        b0h = _mm_add_epi64(b0h, a0h);

        b0l = _mm_add_epi64(b0l, b0h);
        vsse = _mm_add_epi64(vsse, b0l);

        vsum = _mm_add_epi64(vsum, _mm_cvtepi32_epi64(vsum32));
        a += a_stride;
        b += b_stride;
      }
      _mm_storel_epi64((__m128i *)sum, vsum);
      _mm_storel_epi64((__m128i *)sse, vsse);
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
