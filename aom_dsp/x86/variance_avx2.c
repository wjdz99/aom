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

static INLINE void get16xn_var_avx2(const unsigned char *src_ptr,
                                    int source_stride,
                                    const unsigned char *ref_ptr,
                                    int recon_stride, int height,
                                    unsigned int *SSE, int *Sum) {
  __m256i src, src_expand_low, src_expand_high, ref, ref_expand_low;
  __m256i ref_expand_high, madd_low, madd_high;
  unsigned int src_2strides, ref_2strides;
  __m256i zero_reg = _mm256_set1_epi16(0);
  __m256i sum_ref_src = _mm256_set1_epi16(0);
  __m256i madd_ref_src = _mm256_set1_epi16(0);

  // processing two strides in a 256 bit register reducing the number
  // of loop stride by half (comparing to the sse2 code)
  src_2strides = source_stride << 1;
  ref_2strides = recon_stride << 1;
  height = (height >> 1);
  for (int i = 0; i < height; i++) {
    src = _mm256_castsi128_si256(_mm_loadu_si128((__m128i const *)(src_ptr)));
    src = _mm256_inserti128_si256(
        src, _mm_loadu_si128((__m128i const *)(src_ptr + source_stride)), 1);

    ref = _mm256_castsi128_si256(_mm_loadu_si128((__m128i const *)(ref_ptr)));
    ref = _mm256_inserti128_si256(
        ref, _mm_loadu_si128((__m128i const *)(ref_ptr + recon_stride)), 1);

    // expanding to 16 bit each lane
    src_expand_low = _mm256_unpacklo_epi8(src, zero_reg);
    src_expand_high = _mm256_unpackhi_epi8(src, zero_reg);

    ref_expand_low = _mm256_unpacklo_epi8(ref, zero_reg);
    ref_expand_high = _mm256_unpackhi_epi8(ref, zero_reg);

    // src-ref
    src_expand_low = _mm256_sub_epi16(src_expand_low, ref_expand_low);
    src_expand_high = _mm256_sub_epi16(src_expand_high, ref_expand_high);

    // madd low (src - ref)
    madd_low = _mm256_madd_epi16(src_expand_low, src_expand_low);

    // add high to low
    src_expand_low = _mm256_add_epi16(src_expand_low, src_expand_high);

    // madd high (src - ref)
    madd_high = _mm256_madd_epi16(src_expand_high, src_expand_high);

    sum_ref_src = _mm256_add_epi16(sum_ref_src, src_expand_low);

    // add high to low
    madd_ref_src =
        _mm256_add_epi32(madd_ref_src, _mm256_add_epi32(madd_low, madd_high));

    src_ptr += src_2strides;
    ref_ptr += ref_2strides;
  }

  {
    __m128i sum_res, madd_res;
    __m128i expand_sum_low, expand_sum_high, expand_sum;
    __m128i expand_madd_low, expand_madd_high, expand_madd;
    __m128i ex_expand_sum_low, ex_expand_sum_high, ex_expand_sum;

    // extract the low lane and add it to the high lane
    sum_res = _mm_add_epi16(_mm256_castsi256_si128(sum_ref_src),
                            _mm256_extractf128_si256(sum_ref_src, 1));

    madd_res = _mm_add_epi32(_mm256_castsi256_si128(madd_ref_src),
                             _mm256_extractf128_si256(madd_ref_src, 1));

    // padding each 2 bytes with another 2 zeroed bytes
    expand_sum_low =
        _mm_unpacklo_epi16(_mm256_castsi256_si128(zero_reg), sum_res);
    expand_sum_high =
        _mm_unpackhi_epi16(_mm256_castsi256_si128(zero_reg), sum_res);

    // shifting the sign 16 bits right
    expand_sum_low = _mm_srai_epi32(expand_sum_low, 16);
    expand_sum_high = _mm_srai_epi32(expand_sum_high, 16);

    expand_sum = _mm_add_epi32(expand_sum_low, expand_sum_high);

    // expand each 32 bits of the madd result to 64 bits
    expand_madd_low =
        _mm_unpacklo_epi32(madd_res, _mm256_castsi256_si128(zero_reg));
    expand_madd_high =
        _mm_unpackhi_epi32(madd_res, _mm256_castsi256_si128(zero_reg));

    expand_madd = _mm_add_epi32(expand_madd_low, expand_madd_high);

    ex_expand_sum_low =
        _mm_unpacklo_epi32(expand_sum, _mm256_castsi256_si128(zero_reg));
    ex_expand_sum_high =
        _mm_unpackhi_epi32(expand_sum, _mm256_castsi256_si128(zero_reg));

    ex_expand_sum = _mm_add_epi32(ex_expand_sum_low, ex_expand_sum_high);

    // shift 8 bytes eight
    madd_res = _mm_srli_si128(expand_madd, 8);
    sum_res = _mm_srli_si128(ex_expand_sum, 8);

    madd_res = _mm_add_epi32(madd_res, expand_madd);
    sum_res = _mm_add_epi32(sum_res, ex_expand_sum);

    *((int *)SSE) = _mm_cvtsi128_si32(madd_res);

    *((int *)Sum) = _mm_cvtsi128_si32(sum_res);
  }
  _mm256_zeroupper();
}

static INLINE void get32xn_var_avx2(const unsigned char *src_ptr,
                                    int source_stride,
                                    const unsigned char *ref_ptr,
                                    int recon_stride, int height,
                                    unsigned int *SSE, int *Sum) {
  __m256i src, src_expand_low, src_expand_high, ref, ref_expand_low;
  __m256i ref_expand_high, madd_low, madd_high;
  __m256i zero_reg = _mm256_set1_epi16(0);
  __m256i sum_ref_src = _mm256_set1_epi16(0);
  __m256i madd_ref_src = _mm256_set1_epi16(0);

  // processing 32 elements in parallel
  for (int i = 0; i < height; i++) {
    src = _mm256_loadu_si256((__m256i const *)(src_ptr));

    ref = _mm256_loadu_si256((__m256i const *)(ref_ptr));

    // expanding to 16 bit each lane
    src_expand_low = _mm256_unpacklo_epi8(src, zero_reg);
    src_expand_high = _mm256_unpackhi_epi8(src, zero_reg);

    ref_expand_low = _mm256_unpacklo_epi8(ref, zero_reg);
    ref_expand_high = _mm256_unpackhi_epi8(ref, zero_reg);

    // src-ref
    src_expand_low = _mm256_sub_epi16(src_expand_low, ref_expand_low);
    src_expand_high = _mm256_sub_epi16(src_expand_high, ref_expand_high);

    // madd low (src - ref)
    madd_low = _mm256_madd_epi16(src_expand_low, src_expand_low);

    // add high to low
    src_expand_low = _mm256_add_epi16(src_expand_low, src_expand_high);

    // madd high (src - ref)
    madd_high = _mm256_madd_epi16(src_expand_high, src_expand_high);

    sum_ref_src = _mm256_add_epi16(sum_ref_src, src_expand_low);

    // add high to low
    madd_ref_src =
        _mm256_add_epi32(madd_ref_src, _mm256_add_epi32(madd_low, madd_high));

    src_ptr += source_stride;
    ref_ptr += recon_stride;
  }

  {
    __m256i expand_sum_low, expand_sum_high, expand_sum;
    __m256i expand_madd_low, expand_madd_high, expand_madd;
    __m256i ex_expand_sum_low, ex_expand_sum_high, ex_expand_sum;

    // padding each 2 bytes with another 2 zeroed bytes
    expand_sum_low = _mm256_unpacklo_epi16(zero_reg, sum_ref_src);
    expand_sum_high = _mm256_unpackhi_epi16(zero_reg, sum_ref_src);

    // shifting the sign 16 bits right
    expand_sum_low = _mm256_srai_epi32(expand_sum_low, 16);
    expand_sum_high = _mm256_srai_epi32(expand_sum_high, 16);

    expand_sum = _mm256_add_epi32(expand_sum_low, expand_sum_high);

    // expand each 32 bits of the madd result to 64 bits
    expand_madd_low = _mm256_unpacklo_epi32(madd_ref_src, zero_reg);
    expand_madd_high = _mm256_unpackhi_epi32(madd_ref_src, zero_reg);

    expand_madd = _mm256_add_epi32(expand_madd_low, expand_madd_high);

    ex_expand_sum_low = _mm256_unpacklo_epi32(expand_sum, zero_reg);
    ex_expand_sum_high = _mm256_unpackhi_epi32(expand_sum, zero_reg);

    ex_expand_sum = _mm256_add_epi32(ex_expand_sum_low, ex_expand_sum_high);

    // shift 8 bytes eight
    madd_ref_src = _mm256_srli_si256(expand_madd, 8);
    sum_ref_src = _mm256_srli_si256(ex_expand_sum, 8);

    madd_ref_src = _mm256_add_epi32(madd_ref_src, expand_madd);
    sum_ref_src = _mm256_add_epi32(sum_ref_src, ex_expand_sum);

    // extract the low lane and the high lane and add the results
    *((int *)SSE) =
        _mm_cvtsi128_si32(_mm256_castsi256_si128(madd_ref_src)) +
        _mm_cvtsi128_si32(_mm256_extractf128_si256(madd_ref_src, 1));

    *((int *)Sum) = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_ref_src)) +
                    _mm_cvtsi128_si32(_mm256_extractf128_si256(sum_ref_src, 1));
  }
  _mm256_zeroupper();
}

#define AOM_VAR_AVX2(bw, bh, uw, uh, bits)                                     \
  unsigned int aom_variance##bw##x##bh##_avx2(                                 \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,  \
      unsigned int *sse) {                                                     \
    int sum = 0;                                                               \
    unsigned int variance;                                                     \
    *sse = 0;                                                                  \
    for (int i = 0; i < bh; i += uh) {                                         \
      const uint8_t *src_i = src + i * src_stride;                             \
      const uint8_t *ref_i = ref + i * ref_stride;                             \
      for (int j = 0; j < bw; j += uw) {                                       \
        unsigned int sse0 = 0;                                                 \
        int sum0 = 0;                                                          \
        get##uw##xn_var_avx2(src_i + j, src_stride, ref_i + j, ref_stride, uh, \
                             &sse0, &sum0);                                    \
        *sse += sse0;                                                          \
        sum += sum0;                                                           \
      }                                                                        \
    }                                                                          \
    variance = *sse - (uint32_t)(((int64_t)sum * sum) >> bits);                \
    _mm256_zeroupper();                                                        \
    return variance;                                                           \
  }

AOM_VAR_AVX2(16, 4, 16, 4, 6);
AOM_VAR_AVX2(16, 8, 16, 8, 7);
AOM_VAR_AVX2(16, 16, 16, 16, 8);
AOM_VAR_AVX2(16, 32, 16, 16, 9);
AOM_VAR_AVX2(16, 64, 16, 16, 10);

AOM_VAR_AVX2(32, 8, 32, 8, 8);
AOM_VAR_AVX2(32, 16, 32, 16, 9);
AOM_VAR_AVX2(32, 32, 32, 16, 10);
AOM_VAR_AVX2(32, 64, 32, 16, 11);

AOM_VAR_AVX2(64, 16, 32, 16, 10);
AOM_VAR_AVX2(64, 32, 32, 16, 11);
AOM_VAR_AVX2(64, 64, 32, 16, 12);
AOM_VAR_AVX2(64, 128, 32, 16, 13);

AOM_VAR_AVX2(128, 128, 32, 16, 14);
AOM_VAR_AVX2(128, 64, 32, 16, 13);

unsigned int aom_mse16x16_avx2(const uint8_t *src, int src_stride,
                               const uint8_t *ref, int ref_stride,
                               unsigned int *sse) {
  aom_variance16x16_avx2(src, src_stride, ref, ref_stride, sse);
  return *sse;
}

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
