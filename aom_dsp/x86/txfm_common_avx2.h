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

#ifndef AOM_AOM_DSP_X86_TXFM_COMMON_AVX2_H_
#define AOM_AOM_DSP_X86_TXFM_COMMON_AVX2_H_

#include <emmintrin.h>
#include "aom/aom_integer.h"
#include "aom_dsp/x86/synonyms.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE __m256i pair_set_w16_epi16(int16_t a, int16_t b) {
  return _mm256_set1_epi32(
      (int32_t)(((uint16_t)(a)) | (((uint32_t)(b)) << 16)));
}

static INLINE void btf_16_w16_avx2(const __m256i w0, const __m256i w1,
                                   __m256i *in0, __m256i *in1, const __m256i _r,
                                   const int32_t cos_bit) {
  __m256i t0 = _mm256_unpacklo_epi16(*in0, *in1);
  __m256i t1 = _mm256_unpackhi_epi16(*in0, *in1);
  __m256i u0 = _mm256_madd_epi16(t0, w0);
  __m256i u1 = _mm256_madd_epi16(t1, w0);
  __m256i v0 = _mm256_madd_epi16(t0, w1);
  __m256i v1 = _mm256_madd_epi16(t1, w1);

  __m256i a0 = _mm256_add_epi32(u0, _r);
  __m256i a1 = _mm256_add_epi32(u1, _r);
  __m256i b0 = _mm256_add_epi32(v0, _r);
  __m256i b1 = _mm256_add_epi32(v1, _r);

  __m256i c0 = _mm256_srai_epi32(a0, cos_bit);
  __m256i c1 = _mm256_srai_epi32(a1, cos_bit);
  __m256i d0 = _mm256_srai_epi32(b0, cos_bit);
  __m256i d1 = _mm256_srai_epi32(b1, cos_bit);

  *in0 = _mm256_packs_epi32(c0, c1);
  *in1 = _mm256_packs_epi32(d0, d1);
}

static INLINE void btf_16_adds_subs_avx2(__m256i *in0, __m256i *in1) {
  const __m256i _in0 = *in0;
  const __m256i _in1 = *in1;
  *in0 = _mm256_adds_epi16(_in0, _in1);
  *in1 = _mm256_subs_epi16(_in0, _in1);
}

static INLINE void btf_32_add_sub_avx2(__m256i *in0, __m256i *in1) {
  const __m256i _in0 = *in0;
  const __m256i _in1 = *in1;
  *in0 = _mm256_add_epi32(_in0, _in1);
  *in1 = _mm256_sub_epi32(_in0, _in1);
}

static INLINE void btf_16_adds_subs_out_avx2(__m256i *out0, __m256i *out1,
                                             __m256i in0, __m256i in1) {
  const __m256i _in0 = in0;
  const __m256i _in1 = in1;
  *out0 = _mm256_adds_epi16(_in0, _in1);
  *out1 = _mm256_subs_epi16(_in0, _in1);
}

static INLINE void btf_32_add_sub_out_avx2(__m256i *out0, __m256i *out1,
                                           __m256i in0, __m256i in1) {
  const __m256i _in0 = in0;
  const __m256i _in1 = in1;
  *out0 = _mm256_add_epi32(_in0, _in1);
  *out1 = _mm256_sub_epi32(_in0, _in1);
}

static INLINE __m256i load_16bit_to_16bit_avx2(const int16_t *a) {
  return _mm256_load_si256((const __m256i *)a);
}

static INLINE void load_buffer_16bit_to_16bit_avx2(const int16_t *in,
                                                   int stride, __m256i *out,
                                                   int out_size) {
  for (int i = 0; i < out_size; ++i) {
    out[i] = load_16bit_to_16bit_avx2(in + i * stride);
  }
}

static INLINE void load_buffer_16bit_to_16bit_flip_avx2(const int16_t *in,
                                                        int stride,
                                                        __m256i *out,
                                                        int out_size) {
  for (int i = 0; i < out_size; ++i) {
    out[out_size - i - 1] = load_16bit_to_16bit_avx2(in + i * stride);
  }
}

static INLINE __m256i load_32bit_to_16bit_w16_avx2(const int32_t *a) {
  const __m256i a_low = _mm256_lddqu_si256((const __m256i *)a);
  const __m256i b = _mm256_packs_epi32(a_low, *(const __m256i *)(a + 8));
  return _mm256_permute4x64_epi64(b, 0xD8);
}

static INLINE void load_buffer_32bit_to_16bit_w16_avx2(const int32_t *in,
                                                       int stride, __m256i *out,
                                                       int out_size) {
  for (int i = 0; i < out_size; ++i) {
    out[i] = load_32bit_to_16bit_w16_avx2(in + i * stride);
  }
}

static INLINE void transpose_16bit_16x16_avx2(const __m256i *const in,
                                              __m256i *const out) {
  // Unpack 16 bit elements. Goes from:
  // in[0]: 00 01 02 03  08 09 0a 0b  04 05 06 07  0c 0d 0e 0f
  // in[1]: 10 11 12 13  18 19 1a 1b  14 15 16 17  1c 1d 1e 1f
  // in[2]: 20 21 22 23  28 29 2a 2b  24 25 26 27  2c 2d 2e 2f
  // in[3]: 30 31 32 33  38 39 3a 3b  34 35 36 37  3c 3d 3e 3f
  // in[4]: 40 41 42 43  48 49 4a 4b  44 45 46 47  4c 4d 4e 4f
  // in[5]: 50 51 52 53  58 59 5a 5b  54 55 56 57  5c 5d 5e 5f
  // in[6]: 60 61 62 63  68 69 6a 6b  64 65 66 67  6c 6d 6e 6f
  // in[7]: 70 71 72 73  78 79 7a 7b  74 75 76 77  7c 7d 7e 7f
  // in[8]: 80 81 82 83  88 89 8a 8b  84 85 86 87  8c 8d 8e 8f
  // to:
  // a0:    00 10 01 11  02 12 03 13  04 14 05 15  06 16 07 17
  // a1:    20 30 21 31  22 32 23 33  24 34 25 35  26 36 27 37
  // a2:    40 50 41 51  42 52 43 53  44 54 45 55  46 56 47 57
  // a3:    60 70 61 71  62 72 63 73  64 74 65 75  66 76 67 77
  // ...
  __m256i a[16];
  for (int i = 0; i < 16; i += 2) {
    a[i / 2 + 0] = _mm256_unpacklo_epi16(in[i], in[i + 1]);
    a[i / 2 + 8] = _mm256_unpackhi_epi16(in[i], in[i + 1]);
  }
  __m256i b[16];
  for (int i = 0; i < 16; i += 2) {
    b[i / 2 + 0] = _mm256_unpacklo_epi32(a[i], a[i + 1]);
    b[i / 2 + 8] = _mm256_unpackhi_epi32(a[i], a[i + 1]);
  }
  __m256i c[16];
  for (int i = 0; i < 16; i += 2) {
    c[i / 2 + 0] = _mm256_unpacklo_epi64(b[i], b[i + 1]);
    c[i / 2 + 8] = _mm256_unpackhi_epi64(b[i], b[i + 1]);
  }
  out[0 + 0] = _mm256_permute2x128_si256(c[0], c[1], 0x20);
  out[1 + 0] = _mm256_permute2x128_si256(c[8], c[9], 0x20);
  out[2 + 0] = _mm256_permute2x128_si256(c[4], c[5], 0x20);
  out[3 + 0] = _mm256_permute2x128_si256(c[12], c[13], 0x20);

  out[0 + 8] = _mm256_permute2x128_si256(c[0], c[1], 0x31);
  out[1 + 8] = _mm256_permute2x128_si256(c[8], c[9], 0x31);
  out[2 + 8] = _mm256_permute2x128_si256(c[4], c[5], 0x31);
  out[3 + 8] = _mm256_permute2x128_si256(c[12], c[13], 0x31);

  out[4 + 0] = _mm256_permute2x128_si256(c[0 + 2], c[1 + 2], 0x20);
  out[5 + 0] = _mm256_permute2x128_si256(c[8 + 2], c[9 + 2], 0x20);
  out[6 + 0] = _mm256_permute2x128_si256(c[4 + 2], c[5 + 2], 0x20);
  out[7 + 0] = _mm256_permute2x128_si256(c[12 + 2], c[13 + 2], 0x20);

  out[4 + 8] = _mm256_permute2x128_si256(c[0 + 2], c[1 + 2], 0x31);
  out[5 + 8] = _mm256_permute2x128_si256(c[8 + 2], c[9 + 2], 0x31);
  out[6 + 8] = _mm256_permute2x128_si256(c[4 + 2], c[5 + 2], 0x31);
  out[7 + 8] = _mm256_permute2x128_si256(c[12 + 2], c[13 + 2], 0x31);
}
static INLINE void transpose_16bit_16x8_avx2(const __m256i *const in,
                                             __m256i *const out) {
  const __m256i a0 = _mm256_unpacklo_epi16(in[0], in[1]);
  const __m256i a1 = _mm256_unpacklo_epi16(in[2], in[3]);
  const __m256i a2 = _mm256_unpacklo_epi16(in[4], in[5]);
  const __m256i a3 = _mm256_unpacklo_epi16(in[6], in[7]);
  const __m256i a4 = _mm256_unpackhi_epi16(in[0], in[1]);
  const __m256i a5 = _mm256_unpackhi_epi16(in[2], in[3]);
  const __m256i a6 = _mm256_unpackhi_epi16(in[4], in[5]);
  const __m256i a7 = _mm256_unpackhi_epi16(in[6], in[7]);

  const __m256i b0 = _mm256_unpacklo_epi32(a0, a1);
  const __m256i b1 = _mm256_unpacklo_epi32(a2, a3);
  const __m256i b2 = _mm256_unpacklo_epi32(a4, a5);
  const __m256i b3 = _mm256_unpacklo_epi32(a6, a7);
  const __m256i b4 = _mm256_unpackhi_epi32(a0, a1);
  const __m256i b5 = _mm256_unpackhi_epi32(a2, a3);
  const __m256i b6 = _mm256_unpackhi_epi32(a4, a5);
  const __m256i b7 = _mm256_unpackhi_epi32(a6, a7);

  out[0] = _mm256_unpacklo_epi64(b0, b1);
  out[1] = _mm256_unpackhi_epi64(b0, b1);
  out[2] = _mm256_unpacklo_epi64(b4, b5);
  out[3] = _mm256_unpackhi_epi64(b4, b5);
  out[4] = _mm256_unpacklo_epi64(b2, b3);
  out[5] = _mm256_unpackhi_epi64(b2, b3);
  out[6] = _mm256_unpacklo_epi64(b6, b7);
  out[7] = _mm256_unpackhi_epi64(b6, b7);
}
static INLINE void flip_buf_avx2(__m256i *in, __m256i *out, int size) {
  for (int i = 0; i < size; ++i) {
    out[size - i - 1] = in[i];
  }
}

static INLINE void round_shift_16bit_w16_avx2(__m256i *in, int size, int bit) {
  if (bit < 0) {
    bit = -bit;
    __m256i round = _mm256_set1_epi16(1 << (bit - 1));
    for (int i = 0; i < size; ++i) {
      in[i] = _mm256_adds_epi16(in[i], round);
      in[i] = _mm256_srai_epi16(in[i], bit);
    }
  } else if (bit > 0) {
    for (int i = 0; i < size; ++i) {
      in[i] = _mm256_slli_epi16(in[i], bit);
    }
  }
}

static INLINE __m256i av1_round_shift_32_avx2(__m256i vec, int bit) {
  __m256i tmp, round;
  round = _mm256_set1_epi32(1 << (bit - 1));
  tmp = _mm256_add_epi32(vec, round);
  return _mm256_srai_epi32(tmp, bit);
}

static INLINE void av1_round_shift_array_32_avx2(__m256i *input,
                                                 __m256i *output,
                                                 const int size,
                                                 const int bit) {
  if (bit > 0) {
    int i;
    for (i = 0; i < size; i++) {
      output[i] = av1_round_shift_32_avx2(input[i], bit);
    }
  } else {
    int i;
    for (i = 0; i < size; i++) {
      output[i] = _mm256_slli_epi32(input[i], -bit);
    }
  }
}

static INLINE void av1_round_shift_rect_array_32_avx2(__m256i *input,
                                                      __m256i *output,
                                                      const int size,
                                                      const int bit,
                                                      const int val) {
  const __m256i sqrt2 = _mm256_set1_epi32(val);
  if (bit > 0) {
    int i;
    for (i = 0; i < size; i++) {
      const __m256i r0 = av1_round_shift_32_avx2(input[i], bit);
      const __m256i r1 = _mm256_mullo_epi32(sqrt2, r0);
      output[i] = av1_round_shift_32_avx2(r1, NewSqrt2Bits);
    }
  } else {
    int i;
    for (i = 0; i < size; i++) {
      const __m256i r0 = _mm256_slli_epi32(input[i], -bit);
      const __m256i r1 = _mm256_mullo_epi32(sqrt2, r0);
      output[i] = av1_round_shift_32_avx2(r1, NewSqrt2Bits);
    }
  }
}
static INLINE __m256i scale_round_avx2(const __m256i a, const int scale) {
  const __m256i scale_rounding =
      pair_set_w16_epi16(scale, 1 << (NewSqrt2Bits - 1));
  const __m256i b = _mm256_madd_epi16(a, scale_rounding);
  return _mm256_srai_epi32(b, NewSqrt2Bits);
}
static INLINE void store_rect_16bit_to_32bit_w8_avx2(const __m256i a,
                                                     int32_t *const b) {
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i a_lo = _mm256_unpacklo_epi16(a, one);
  const __m256i a_hi = _mm256_unpackhi_epi16(a, one);
  const __m256i b_lo = scale_round_avx2(a_lo, NewSqrt2);
  const __m256i b_hi = scale_round_avx2(a_hi, NewSqrt2);
  const __m256i temp = _mm256_permute2f128_si256(b_lo, b_hi, 0x31);
  _mm_store_si128((__m128i *)b, _mm256_castsi256_si128(b_lo));
  _mm_store_si128((__m128i *)(b + 4), _mm256_castsi256_si128(b_hi));
  _mm256_store_si256((__m256i *)(b + 64), temp);
}
static INLINE void store_rect_buffer_16bit_to_32bit_w8_avx2(
    const __m256i *const in, int32_t *const out, const int stride,
    const int out_size) {
  for (int i = 0; i < out_size; ++i) {
    store_rect_16bit_to_32bit_w8_avx2(in[i], out + i * stride);
  }
}
static INLINE void pack_reg(const __m128i *in1, const __m128i *in2,
                            __m256i *out) {
  out[0] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[0]), in2[0], 0x1);
  out[1] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[1]), in2[1], 0x1);
  out[2] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[2]), in2[2], 0x1);
  out[3] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[3]), in2[3], 0x1);
  out[4] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[4]), in2[4], 0x1);
  out[5] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[5]), in2[5], 0x1);
  out[6] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[6]), in2[6], 0x1);
  out[7] = _mm256_insertf128_si256(_mm256_castsi128_si256(in1[7]), in2[7], 0x1);
}
static INLINE void extract_reg(const __m256i *in, __m128i *out1) {
  out1[0] = _mm256_castsi256_si128(in[0]);
  out1[1] = _mm256_castsi256_si128(in[1]);
  out1[2] = _mm256_castsi256_si128(in[2]);
  out1[3] = _mm256_castsi256_si128(in[3]);
  out1[4] = _mm256_castsi256_si128(in[4]);
  out1[5] = _mm256_castsi256_si128(in[5]);
  out1[6] = _mm256_castsi256_si128(in[6]);
  out1[7] = _mm256_castsi256_si128(in[7]);

  out1[8] = _mm256_extractf128_si256(in[0], 0x01);
  out1[9] = _mm256_extractf128_si256(in[1], 0x01);
  out1[10] = _mm256_extractf128_si256(in[2], 0x01);
  out1[11] = _mm256_extractf128_si256(in[3], 0x01);
  out1[12] = _mm256_extractf128_si256(in[4], 0x01);
  out1[13] = _mm256_extractf128_si256(in[5], 0x01);
  out1[14] = _mm256_extractf128_si256(in[6], 0x01);
  out1[15] = _mm256_extractf128_si256(in[7], 0x01);
}
#ifdef __cplusplus
}
#endif

#endif  // AOM_AOM_DSP_X86_TXFM_COMMON_AVX2_H_
