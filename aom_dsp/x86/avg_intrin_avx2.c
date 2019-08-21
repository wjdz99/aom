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

#include "config/aom_dsp_rtcd.h"
#include "aom/aom_integer.h"
#include "aom_dsp/x86/bitdepth_conversion_avx2.h"
#include "aom_ports/mem.h"

static INLINE __m256i load_tran_low(const tran_low_t *a) {
  const __m256i a_low = _mm256_loadu_si256((const __m256i *)a);
  const __m256i a_high = _mm256_loadu_si256((const __m256i *)(a + 8));
  return _mm256_packs_epi32(a_low, a_high);
}

static INLINE void store_tran_low(__m256i a, tran_low_t *b) {
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i a_hi = _mm256_mulhi_epi16(a, one);
  const __m256i a_lo = _mm256_mullo_epi16(a, one);
  const __m256i a_1 = _mm256_unpacklo_epi16(a_lo, a_hi);
  const __m256i a_2 = _mm256_unpackhi_epi16(a_lo, a_hi);
  _mm256_storeu_si256((__m256i *)b, a_1);
  _mm256_storeu_si256((__m256i *)(b + 8), a_2);
}
static void hadamard_col8x2_avx2(__m256i *in, int iter) {
  __m256i a0 = in[0];
  __m256i a1 = in[1];
  __m256i a2 = in[2];
  __m256i a3 = in[3];
  __m256i a4 = in[4];
  __m256i a5 = in[5];
  __m256i a6 = in[6];
  __m256i a7 = in[7];

  __m256i b0 = _mm256_add_epi16(a0, a1);
  __m256i b1 = _mm256_sub_epi16(a0, a1);
  __m256i b2 = _mm256_add_epi16(a2, a3);
  __m256i b3 = _mm256_sub_epi16(a2, a3);
  __m256i b4 = _mm256_add_epi16(a4, a5);
  __m256i b5 = _mm256_sub_epi16(a4, a5);
  __m256i b6 = _mm256_add_epi16(a6, a7);
  __m256i b7 = _mm256_sub_epi16(a6, a7);

  a0 = _mm256_add_epi16(b0, b2);
  a1 = _mm256_add_epi16(b1, b3);
  a2 = _mm256_sub_epi16(b0, b2);
  a3 = _mm256_sub_epi16(b1, b3);
  a4 = _mm256_add_epi16(b4, b6);
  a5 = _mm256_add_epi16(b5, b7);
  a6 = _mm256_sub_epi16(b4, b6);
  a7 = _mm256_sub_epi16(b5, b7);

  if (iter == 0) {
    b0 = _mm256_add_epi16(a0, a4);
    b7 = _mm256_add_epi16(a1, a5);
    b3 = _mm256_add_epi16(a2, a6);
    b4 = _mm256_add_epi16(a3, a7);
    b2 = _mm256_sub_epi16(a0, a4);
    b6 = _mm256_sub_epi16(a1, a5);
    b1 = _mm256_sub_epi16(a2, a6);
    b5 = _mm256_sub_epi16(a3, a7);

    a0 = _mm256_unpacklo_epi16(b0, b1);
    a1 = _mm256_unpacklo_epi16(b2, b3);
    a2 = _mm256_unpackhi_epi16(b0, b1);
    a3 = _mm256_unpackhi_epi16(b2, b3);
    a4 = _mm256_unpacklo_epi16(b4, b5);
    a5 = _mm256_unpacklo_epi16(b6, b7);
    a6 = _mm256_unpackhi_epi16(b4, b5);
    a7 = _mm256_unpackhi_epi16(b6, b7);

    b0 = _mm256_unpacklo_epi32(a0, a1);
    b1 = _mm256_unpacklo_epi32(a4, a5);
    b2 = _mm256_unpackhi_epi32(a0, a1);
    b3 = _mm256_unpackhi_epi32(a4, a5);
    b4 = _mm256_unpacklo_epi32(a2, a3);
    b5 = _mm256_unpacklo_epi32(a6, a7);
    b6 = _mm256_unpackhi_epi32(a2, a3);
    b7 = _mm256_unpackhi_epi32(a6, a7);

    in[0] = _mm256_unpacklo_epi64(b0, b1);
    in[1] = _mm256_unpackhi_epi64(b0, b1);
    in[2] = _mm256_unpacklo_epi64(b2, b3);
    in[3] = _mm256_unpackhi_epi64(b2, b3);
    in[4] = _mm256_unpacklo_epi64(b4, b5);
    in[5] = _mm256_unpackhi_epi64(b4, b5);
    in[6] = _mm256_unpacklo_epi64(b6, b7);
    in[7] = _mm256_unpackhi_epi64(b6, b7);
  } else {
    in[0] = _mm256_add_epi16(a0, a4);
    in[7] = _mm256_add_epi16(a1, a5);
    in[3] = _mm256_add_epi16(a2, a6);
    in[4] = _mm256_add_epi16(a3, a7);
    in[2] = _mm256_sub_epi16(a0, a4);
    in[6] = _mm256_sub_epi16(a1, a5);
    in[1] = _mm256_sub_epi16(a2, a6);
    in[5] = _mm256_sub_epi16(a3, a7);
  }
}

static void hadamard_8x8x2_avx2(const int16_t *src_diff, ptrdiff_t src_stride,
                                int16_t *coeff) {
  __m256i src[8];
  src[0] = _mm256_loadu_si256((const __m256i *)src_diff);
  src[1] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[2] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[3] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[4] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[5] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[6] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[7] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));

  hadamard_col8x2_avx2(src, 0);
  hadamard_col8x2_avx2(src, 1);

  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[0], src[1], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[2], src[3], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[4], src[5], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[6], src[7], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[0], src[1], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[2], src[3], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[4], src[5], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[6], src[7], 0x31));
}

static INLINE void hadamard_16x16_avx2(const int16_t *src_diff,
                                       ptrdiff_t src_stride, tran_low_t *coeff,
                                       int is_final) {
  DECLARE_ALIGNED(32, int16_t, temp_coeff[16 * 16]);
  int16_t *t_coeff = temp_coeff;
  int16_t *coeff16 = (int16_t *)coeff;
  int idx;
  for (idx = 0; idx < 2; ++idx) {
    const int16_t *src_ptr = src_diff + idx * 8 * src_stride;
    hadamard_8x8x2_avx2(src_ptr, src_stride, t_coeff + (idx * 64 * 2));
  }

  for (idx = 0; idx < 64; idx += 16) {
    const __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    const __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 64));
    const __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 128));
    const __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 192));

    __m256i b0 = _mm256_add_epi16(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi16(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi16(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi16(coeff2, coeff3);

    b0 = _mm256_srai_epi16(b0, 1);
    b1 = _mm256_srai_epi16(b1, 1);
    b2 = _mm256_srai_epi16(b2, 1);
    b3 = _mm256_srai_epi16(b3, 1);
    if (is_final) {
      store_tran_low(_mm256_add_epi16(b0, b2), coeff);
      store_tran_low(_mm256_add_epi16(b1, b3), coeff + 64);
      store_tran_low(_mm256_sub_epi16(b0, b2), coeff + 128);
      store_tran_low(_mm256_sub_epi16(b1, b3), coeff + 192);
      coeff += 16;
    } else {
      _mm256_storeu_si256((__m256i *)coeff16, _mm256_add_epi16(b0, b2));
      _mm256_storeu_si256((__m256i *)(coeff16 + 64), _mm256_add_epi16(b1, b3));
      _mm256_storeu_si256((__m256i *)(coeff16 + 128), _mm256_sub_epi16(b0, b2));
      _mm256_storeu_si256((__m256i *)(coeff16 + 192), _mm256_sub_epi16(b1, b3));
      coeff16 += 16;
    }
    t_coeff += 16;
  }
}

void aom_hadamard_16x16_avx2(const int16_t *src_diff, ptrdiff_t src_stride,
                             tran_low_t *coeff) {
  hadamard_16x16_avx2(src_diff, src_stride, coeff, 1);
}

void aom_hadamard_32x32_avx2(const int16_t *src_diff, ptrdiff_t src_stride,
                             tran_low_t *coeff) {
  // For high bitdepths, it is unnecessary to store_tran_low
  // (mult/unpack/store), then load_tran_low (load/pack) the same memory in the
  // next stage.  Output to an intermediate buffer first, then store_tran_low()
  // in the final stage.
  DECLARE_ALIGNED(32, int16_t, temp_coeff[32 * 32]);
  int16_t *t_coeff = temp_coeff;
  int idx;
  for (idx = 0; idx < 4; ++idx) {
    // src_diff: 9 bit, dynamic range [-255, 255]
    const int16_t *src_ptr =
        src_diff + (idx >> 1) * 16 * src_stride + (idx & 0x01) * 16;
    hadamard_16x16_avx2(src_ptr, src_stride,
                        (tran_low_t *)(t_coeff + idx * 256), 0);
  }

  for (idx = 0; idx < 256; idx += 16) {
    const __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    const __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 256));
    const __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 512));
    const __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 768));

    __m256i b0 = _mm256_add_epi16(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi16(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi16(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi16(coeff2, coeff3);

    b0 = _mm256_srai_epi16(b0, 2);
    b1 = _mm256_srai_epi16(b1, 2);
    b2 = _mm256_srai_epi16(b2, 2);
    b3 = _mm256_srai_epi16(b3, 2);

    store_tran_low(_mm256_add_epi16(b0, b2), coeff);
    store_tran_low(_mm256_add_epi16(b1, b3), coeff + 256);
    store_tran_low(_mm256_sub_epi16(b0, b2), coeff + 512);
    store_tran_low(_mm256_sub_epi16(b1, b3), coeff + 768);

    coeff += 16;
    t_coeff += 16;
  }
}

static void highbd_hadamard_col8_avx2(__m256i *in, int iter) {
  __m256i a0 = in[0];
  __m256i a1 = in[1];
  __m256i a2 = in[2];
  __m256i a3 = in[3];
  __m256i a4 = in[4];
  __m256i a5 = in[5];
  __m256i a6 = in[6];
  __m256i a7 = in[7];

  __m256i b0 = _mm256_add_epi32(a0, a1);
  __m256i b1 = _mm256_sub_epi32(a0, a1);
  __m256i b2 = _mm256_add_epi32(a2, a3);
  __m256i b3 = _mm256_sub_epi32(a2, a3);
  __m256i b4 = _mm256_add_epi32(a4, a5);
  __m256i b5 = _mm256_sub_epi32(a4, a5);
  __m256i b6 = _mm256_add_epi32(a6, a7);
  __m256i b7 = _mm256_sub_epi32(a6, a7);

  a0 = _mm256_add_epi32(b0, b2);
  a1 = _mm256_add_epi32(b1, b3);
  a2 = _mm256_sub_epi32(b0, b2);
  a3 = _mm256_sub_epi32(b1, b3);
  a4 = _mm256_add_epi32(b4, b6);
  a5 = _mm256_add_epi32(b5, b7);
  a6 = _mm256_sub_epi32(b4, b6);
  a7 = _mm256_sub_epi32(b5, b7);

  if (iter == 0) {
    b0 = _mm256_add_epi32(a0, a4);
    b7 = _mm256_add_epi32(a1, a5);
    b3 = _mm256_add_epi32(a2, a6);
    b4 = _mm256_add_epi32(a3, a7);
    b2 = _mm256_sub_epi32(a0, a4);
    b6 = _mm256_sub_epi32(a1, a5);
    b1 = _mm256_sub_epi32(a2, a6);
    b5 = _mm256_sub_epi32(a3, a7);

    a0 = _mm256_unpacklo_epi32(b0, b1);
    a1 = _mm256_unpacklo_epi32(b2, b3);
    a2 = _mm256_unpackhi_epi32(b0, b1);
    a3 = _mm256_unpackhi_epi32(b2, b3);
    a4 = _mm256_unpacklo_epi32(b4, b5);
    a5 = _mm256_unpacklo_epi32(b6, b7);
    a6 = _mm256_unpackhi_epi32(b4, b5);
    a7 = _mm256_unpackhi_epi32(b6, b7);

    b0 = _mm256_unpacklo_epi64(a0, a1);
    b1 = _mm256_unpacklo_epi64(a4, a5);
    b2 = _mm256_unpackhi_epi64(a0, a1);
    b3 = _mm256_unpackhi_epi64(a4, a5);
    b4 = _mm256_unpacklo_epi64(a2, a3);
    b5 = _mm256_unpacklo_epi64(a6, a7);
    b6 = _mm256_unpackhi_epi64(a2, a3);
    b7 = _mm256_unpackhi_epi64(a6, a7);

    in[0] = _mm256_permute2x128_si256(b0, b1, 0x20);
    in[1] = _mm256_permute2x128_si256(b0, b1, 0x31);
    in[2] = _mm256_permute2x128_si256(b2, b3, 0x20);
    in[3] = _mm256_permute2x128_si256(b2, b3, 0x31);
    in[4] = _mm256_permute2x128_si256(b4, b5, 0x20);
    in[5] = _mm256_permute2x128_si256(b4, b5, 0x31);
    in[6] = _mm256_permute2x128_si256(b6, b7, 0x20);
    in[7] = _mm256_permute2x128_si256(b6, b7, 0x31);
  } else {
    in[0] = _mm256_add_epi32(a0, a4);
    in[7] = _mm256_add_epi32(a1, a5);
    in[3] = _mm256_add_epi32(a2, a6);
    in[4] = _mm256_add_epi32(a3, a7);
    in[2] = _mm256_sub_epi32(a0, a4);
    in[6] = _mm256_sub_epi32(a1, a5);
    in[1] = _mm256_sub_epi32(a2, a6);
    in[5] = _mm256_sub_epi32(a3, a7);
  }
}

void aom_highbd_hadamard_8x8_avx2(const int16_t *src_diff, ptrdiff_t src_stride,
                                  tran_low_t *coeff) {
  __m128i src16[8];
  __m256i src32[8];

  src16[0] = _mm_loadu_si128((const __m128i *)src_diff);
  src16[1] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[2] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[3] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[4] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[5] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[6] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));
  src16[7] = _mm_loadu_si128((const __m128i *)(src_diff += src_stride));

  src32[0] = _mm256_cvtepi16_epi32(src16[0]);
  src32[1] = _mm256_cvtepi16_epi32(src16[1]);
  src32[2] = _mm256_cvtepi16_epi32(src16[2]);
  src32[3] = _mm256_cvtepi16_epi32(src16[3]);
  src32[4] = _mm256_cvtepi16_epi32(src16[4]);
  src32[5] = _mm256_cvtepi16_epi32(src16[5]);
  src32[6] = _mm256_cvtepi16_epi32(src16[6]);
  src32[7] = _mm256_cvtepi16_epi32(src16[7]);

  highbd_hadamard_col8_avx2(src32, 0);
  highbd_hadamard_col8_avx2(src32, 1);

  _mm256_storeu_si256((__m256i *)coeff, src32[0]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[1]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[2]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[3]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[4]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[5]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[6]);
  coeff += 8;
  _mm256_storeu_si256((__m256i *)coeff, src32[7]);
}

void aom_highbd_hadamard_16x16_avx2(const int16_t *src_diff,
                                    ptrdiff_t src_stride, tran_low_t *coeff) {
  int idx;
  tran_low_t *t_coeff = coeff;
  for (idx = 0; idx < 4; ++idx) {
    const int16_t *src_ptr =
        src_diff + (idx >> 1) * 8 * src_stride + (idx & 0x01) * 8;
    aom_highbd_hadamard_8x8_avx2(src_ptr, src_stride, t_coeff + idx * 64);
  }

  for (idx = 0; idx < 64; idx += 8) {
    __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 64));
    __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 128));
    __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 192));

    __m256i b0 = _mm256_add_epi32(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi32(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi32(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi32(coeff2, coeff3);

    b0 = _mm256_srai_epi32(b0, 1);
    b1 = _mm256_srai_epi32(b1, 1);
    b2 = _mm256_srai_epi32(b2, 1);
    b3 = _mm256_srai_epi32(b3, 1);

    coeff0 = _mm256_add_epi32(b0, b2);
    coeff1 = _mm256_add_epi32(b1, b3);
    coeff2 = _mm256_sub_epi32(b0, b2);
    coeff3 = _mm256_sub_epi32(b1, b3);

    _mm256_storeu_si256((__m256i *)coeff, coeff0);
    _mm256_storeu_si256((__m256i *)(coeff + 64), coeff1);
    _mm256_storeu_si256((__m256i *)(coeff + 128), coeff2);
    _mm256_storeu_si256((__m256i *)(coeff + 192), coeff3);

    coeff += 8;
    t_coeff += 8;
  }
}

void aom_highbd_hadamard_32x32_avx2(const int16_t *src_diff,
                                    ptrdiff_t src_stride, tran_low_t *coeff) {
  int idx;
  tran_low_t *t_coeff = coeff;
  for (idx = 0; idx < 4; ++idx) {
    const int16_t *src_ptr =
        src_diff + (idx >> 1) * 16 * src_stride + (idx & 0x01) * 16;
    aom_highbd_hadamard_16x16_avx2(src_ptr, src_stride, t_coeff + idx * 256);
  }

  for (idx = 0; idx < 256; idx += 8) {
    __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 256));
    __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 512));
    __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 768));

    __m256i b0 = _mm256_add_epi32(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi32(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi32(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi32(coeff2, coeff3);

    b0 = _mm256_srai_epi32(b0, 2);
    b1 = _mm256_srai_epi32(b1, 2);
    b2 = _mm256_srai_epi32(b2, 2);
    b3 = _mm256_srai_epi32(b3, 2);

    coeff0 = _mm256_add_epi32(b0, b2);
    coeff1 = _mm256_add_epi32(b1, b3);
    coeff2 = _mm256_sub_epi32(b0, b2);
    coeff3 = _mm256_sub_epi32(b1, b3);

    _mm256_storeu_si256((__m256i *)coeff, coeff0);
    _mm256_storeu_si256((__m256i *)(coeff + 256), coeff1);
    _mm256_storeu_si256((__m256i *)(coeff + 512), coeff2);
    _mm256_storeu_si256((__m256i *)(coeff + 768), coeff3);

    coeff += 8;
    t_coeff += 8;
  }
}

int aom_satd_avx2(const tran_low_t *coeff, int length) {
  __m256i accum = _mm256_setzero_si256();
  int i;

  for (i = 0; i < length; i += 8, coeff += 8) {
    const __m256i src_line = _mm256_loadu_si256((const __m256i *)coeff);
    const __m256i abs = _mm256_abs_epi32(src_line);
    accum = _mm256_add_epi32(accum, abs);
  }

  {  // 32 bit horizontal add
    const __m256i a = _mm256_srli_si256(accum, 8);
    const __m256i b = _mm256_add_epi32(accum, a);
    const __m256i c = _mm256_srli_epi64(b, 32);
    const __m256i d = _mm256_add_epi32(b, c);
    const __m128i accum_128 = _mm_add_epi32(_mm256_castsi256_si128(d),
                                            _mm256_extractf128_si256(d, 1));
    return _mm_cvtsi128_si32(accum_128);
  }
}

static void hadamard_row16_1d(const int16_t *src_diff, ptrdiff_t src_stride,
                              int16_t *coeff, int is_final) {
  (void)src_stride;
  int16_t b0 = src_diff[0] + src_diff[1];     // 0+1
  int16_t b1 = src_diff[0] - src_diff[1];     // 0-1
  int16_t b2 = src_diff[2] + src_diff[3];     // 2+3
  int16_t b3 = src_diff[2] - src_diff[3];     // 2-3
  int16_t b4 = src_diff[4] + src_diff[5];     // 4+5
  int16_t b5 = src_diff[4] - src_diff[5];     // 4-5
  int16_t b6 = src_diff[6] + src_diff[7];     // 6+7
  int16_t b7 = src_diff[6] - src_diff[7];     // 6-7
  int16_t b8 = src_diff[8] + src_diff[9];     // 8+9
  int16_t b9 = src_diff[8] - src_diff[9];     // 8-9
  int16_t b10 = src_diff[10] + src_diff[11];  // 10+11
  int16_t b11 = src_diff[10] - src_diff[11];  // 10-11
  int16_t b12 = src_diff[12] + src_diff[13];  // 12+13
  int16_t b13 = src_diff[12] - src_diff[13];  // 12-13
  int16_t b14 = src_diff[14] + src_diff[15];  // 14+15
  int16_t b15 = src_diff[14] - src_diff[15];  // 14-15

  int16_t c0 = b0 + b2;     // 0+1+2+3
  int16_t c1 = b1 + b3;     // 0-1+2-3
  int16_t c2 = b0 - b2;     // 0+1-2-3
  int16_t c3 = b1 - b3;     // 0-1-2+3
  int16_t c4 = b4 + b6;     // 4+5+6+7
  int16_t c5 = b5 + b7;     // 4-5+6-7
  int16_t c6 = b4 - b6;     // 4+5-6-7
  int16_t c7 = b5 - b7;     // 4-5-6+7
  int16_t c8 = b8 + b10;    // 8+9+10+11
  int16_t c9 = b9 + b11;    // 8-9+10-11
  int16_t c10 = b8 - b10;   // 8+9-10-11
  int16_t c11 = b9 - b11;   // 8-9-10+11
  int16_t c12 = b12 + b14;  // 12+13+14+15
  int16_t c13 = b13 + b15;  // 12-13+14-15
  int16_t c14 = b12 - b14;  // 12+13-14-15
  int16_t c15 = b13 - b15;  // 12-13-14+15

  int16_t d0 = c0 + c4;     // 0+1+2+3+4+5+6+7		coeff[0] = c0 + c4;
  int16_t d1 = c1 + c5;     // 0-1+2-3+4-5+6-7		coeff[7] = c1 + c5;
  int16_t d2 = c2 + c6;     // 0+1-2-3+4+5-6-7		coeff[3] = c2 + c6;
  int16_t d3 = c3 + c7;     // 0-1-2+3+4-5-6+7		coeff[4] = c3 + c7;
  int16_t d4 = c0 - c4;     // 0+1-2-3-4-5-6-7		coeff[2] = c0 - c4;
  int16_t d5 = c1 - c5;     // 0-1-2+3-4+5-6+7		coeff[6] = c1 - c5;
  int16_t d6 = c2 - c6;     // 0+1-2-3-4-5+6+7		coeff[1] = c2 - c6;
  int16_t d7 = c3 - c7;     // 0-1-2+3-4+5+6-7		coeff[5] = c3 - c7;
  int16_t d8 = c8 + c12;    // 8+9+10+11+12+13+14+15
  int16_t d9 = c9 + c13;    // 8-9+10-11+12-13+14-15
  int16_t d10 = c10 + c14;  // 8+9-10-11+12+13-14-15
  int16_t d11 = c11 + c15;  // 8-9-10+11+12-13-14+15
  int16_t d12 = c8 - c12;   // 8+9-10-11-12-13-14-15
  int16_t d13 = c9 - c13;   // 8-9-10+11-12+13-14+15
  int16_t d14 = c10 - c14;  // 8+9-10-11-12-13+14+15
  int16_t d15 = c11 - c15;  // 8-9-10+11-12+13+14-15

  if (is_final) {
    coeff[0] = d0 + d8;
    coeff[2] = d1 + d9;
    coeff[4] = d2 + d10;
    coeff[6] = d3 + d11;
    coeff[8] = d4 + d12;
    coeff[10] = d5 + d13;
    coeff[12] = d6 + d14;
    coeff[14] = d7 + d15;
    coeff[1] = d0 - d8;
    coeff[3] = d1 - d9;
    coeff[5] = d2 - d10;
    coeff[7] = d3 - d11;
    coeff[9] = d4 - d12;
    coeff[11] = d5 - d13;
    coeff[13] = d6 - d14;
    coeff[15] = d7 - d15;
  } else {
    coeff[0 * 16] = d0 + d8;
    coeff[2 * 16] = d1 + d9;
    coeff[4 * 16] = d2 + d10;
    coeff[6 * 16] = d3 + d11;
    coeff[8 * 16] = d4 + d12;
    coeff[10 * 16] = d5 + d13;
    coeff[12 * 16] = d6 + d14;
    coeff[14 * 16] = d7 + d15;
    coeff[1 * 16] = d0 - d8;
    coeff[3 * 16] = d1 - d9;
    coeff[5 * 16] = d2 - d10;
    coeff[7 * 16] = d3 - d11;
    coeff[9 * 16] = d4 - d12;
    coeff[11 * 16] = d5 - d13;
    coeff[13 * 16] = d6 - d14;
    coeff[15 * 16] = d7 - d15;
  }
}

static void hadamard_col16_1d(const int16_t *src_diff, ptrdiff_t src_stride,
                              int16_t *coeff, int is_final) {
  (void)is_final;
  int16_t b0 = src_diff[0 * src_stride] + src_diff[1 * src_stride];     // 0+1
  int16_t b1 = src_diff[0 * src_stride] - src_diff[1 * src_stride];     // 0-1
  int16_t b2 = src_diff[2 * src_stride] + src_diff[3 * src_stride];     // 2+3
  int16_t b3 = src_diff[2 * src_stride] - src_diff[3 * src_stride];     // 2-3
  int16_t b4 = src_diff[4 * src_stride] + src_diff[5 * src_stride];     // 4+5
  int16_t b5 = src_diff[4 * src_stride] - src_diff[5 * src_stride];     // 4-5
  int16_t b6 = src_diff[6 * src_stride] + src_diff[7 * src_stride];     // 6+7
  int16_t b7 = src_diff[6 * src_stride] - src_diff[7 * src_stride];     // 6-7
  int16_t b8 = src_diff[8 * src_stride] + src_diff[9 * src_stride];     // 8+9
  int16_t b9 = src_diff[8 * src_stride] - src_diff[9 * src_stride];     // 8-9
  int16_t b10 = src_diff[10 * src_stride] + src_diff[11 * src_stride];  // 10+11
  int16_t b11 = src_diff[10 * src_stride] - src_diff[11 * src_stride];  // 10-11
  int16_t b12 = src_diff[12 * src_stride] + src_diff[13 * src_stride];  // 12+13
  int16_t b13 = src_diff[12 * src_stride] - src_diff[13 * src_stride];  // 12-13
  int16_t b14 = src_diff[14 * src_stride] + src_diff[15 * src_stride];  // 14+15
  int16_t b15 = src_diff[14 * src_stride] - src_diff[15 * src_stride];  // 14-15

  int16_t c0 = b0 + b2;     // 0+1+2+3
  int16_t c1 = b1 + b3;     // 0-1+2-3
  int16_t c2 = b0 - b2;     // 0+1-2-3
  int16_t c3 = b1 - b3;     // 0-1-2+3
  int16_t c4 = b4 + b6;     // 4+5+6+7
  int16_t c5 = b5 + b7;     // 4-5+6-7
  int16_t c6 = b4 - b6;     // 4+5-6-7
  int16_t c7 = b5 - b7;     // 4-5-6+7
  int16_t c8 = b8 + b10;    // 8+9+10+11
  int16_t c9 = b9 + b11;    // 8-9+10-11
  int16_t c10 = b8 - b10;   // 8+9-10-11
  int16_t c11 = b9 - b11;   // 8-9-10+11
  int16_t c12 = b12 + b14;  // 12+13+14+15
  int16_t c13 = b13 + b15;  // 12-13+14-15
  int16_t c14 = b12 - b14;  // 12+13-14-15
  int16_t c15 = b13 - b15;  // 12-13-14+15

  int16_t d0 = c0 + c4;     // 0+1+2+3+4+5+6+7		coeff[0] = c0 + c4;
  int16_t d1 = c1 + c5;     // 0-1+2-3+4-5+6-7		coeff[7] = c1 + c5;
  int16_t d2 = c2 + c6;     // 0+1-2-3+4+5-6-7		coeff[3] = c2 + c6;
  int16_t d3 = c3 + c7;     // 0-1-2+3+4-5-6+7		coeff[4] = c3 + c7;
  int16_t d4 = c0 - c4;     // 0+1-2-3-4-5-6-7		coeff[2] = c0 - c4;
  int16_t d5 = c1 - c5;     // 0-1-2+3-4+5-6+7		coeff[6] = c1 - c5;
  int16_t d6 = c2 - c6;     // 0+1-2-3-4-5+6+7		coeff[1] = c2 - c6;
  int16_t d7 = c3 - c7;     // 0-1-2+3-4+5+6-7		coeff[5] = c3 - c7;
  int16_t d8 = c8 + c12;    // 8+9+10+11+12+13+14+15
  int16_t d9 = c9 + c13;    // 8-9+10-11+12-13+14-15
  int16_t d10 = c10 + c14;  // 8+9-10-11+12+13-14-15
  int16_t d11 = c11 + c15;  // 8-9-10+11+12-13-14+15
  int16_t d12 = c8 - c12;   // 8+9-10-11-12-13-14-15
  int16_t d13 = c9 - c13;   // 8-9-10+11-12+13-14+15
  int16_t d14 = c10 - c14;  // 8+9-10-11-12-13+14+15
  int16_t d15 = c11 - c15;  // 8-9-10+11-12+13+14-15

  coeff[0] = d0 + d8;
  coeff[2] = d1 + d9;
  coeff[4] = d2 + d10;
  coeff[6] = d3 + d11;
  coeff[8] = d4 + d12;
  coeff[10] = d5 + d13;
  coeff[12] = d6 + d14;
  coeff[14] = d7 + d15;
  coeff[1] = d0 - d8;
  coeff[3] = d1 - d9;
  coeff[5] = d2 - d10;
  coeff[7] = d3 - d11;
  coeff[9] = d4 - d12;
  coeff[11] = d5 - d13;
  coeff[13] = d6 - d14;
  coeff[15] = d7 - d15;
}

// In place 16x16 2D Hadamard transform
void aom_hadamard_16x16_using_1d(const int16_t *src_diff, ptrdiff_t src_stride,
                                 tran_low_t *coeff, int *had_satd) {
  int idx;
  int16_t buffer[256];
  int16_t buffer2[256];
  int16_t *tmp_buf = &buffer[0];
  const int16_t *tmp_src = src_diff;
  // vert/col transform
  tmp_buf = &buffer2[0];
  for (idx = 0; idx < 16; ++idx) {
    hadamard_col16_1d(tmp_src, src_stride, tmp_buf,
                      0);  // src_diff: 9 bit
                           // dynamic range [-255, 255]
    tmp_buf += 16;
    tmp_src++;
  }

  tmp_buf = &buffer2[0];
  for (int i = 0; i < 256; i++) had_satd[1] += abs(tmp_buf[i]);

  tmp_buf = &buffer[0];
  // horz/row transform
  for (idx = 0; idx < 16; ++idx) {
    hadamard_row16_1d(src_diff, src_stride, tmp_buf,
                      0);  // src_diff: 9 bit
                           // dynamic range [-255, 255]
    tmp_buf++;
    src_diff += src_stride;
  }
  tmp_buf = &buffer[0];
  for (int i = 0; i < 256; i++) had_satd[0] += abs(tmp_buf[i]);

  for (idx = 0; idx < 16; ++idx) {
    hadamard_row16_1d(tmp_buf, 16, buffer2 + 16 * idx,
                      1);  // tmp_buf: 12 bit
                           // dynamic range [-2040, 2040]
                           // buffer2: 15 bit
                           // dynamic range [-16320, 16320]
    tmp_buf += 16;
  }

  for (idx = 0; idx < 256; ++idx) {
    coeff[idx] = ((tran_low_t)buffer2[idx]) >> 1;
    had_satd[2] += abs(coeff[idx]);
  }

  had_satd[0] = (had_satd[0] >> 2);  // multiplied by (1/4)
  had_satd[1] = (had_satd[1] >> 2);
  had_satd[2] = (had_satd[2] >> 3);
}
static void hadamard_row8(const int16_t *src_diff, ptrdiff_t src_stride,
                          int16_t *coeff) {
  (void)src_stride;
  int16_t b0 = src_diff[0] + src_diff[1];
  int16_t b1 = src_diff[0] - src_diff[1];
  int16_t b2 = src_diff[2] + src_diff[3];
  int16_t b3 = src_diff[2] - src_diff[3];
  int16_t b4 = src_diff[4] + src_diff[5];
  int16_t b5 = src_diff[4] - src_diff[5];
  int16_t b6 = src_diff[6] + src_diff[7];
  int16_t b7 = src_diff[6] - src_diff[7];

  int16_t c0 = b0 + b2;
  int16_t c1 = b1 + b3;
  int16_t c2 = b0 - b2;
  int16_t c3 = b1 - b3;
  int16_t c4 = b4 + b6;
  int16_t c5 = b5 + b7;
  int16_t c6 = b4 - b6;
  int16_t c7 = b5 - b7;

  coeff[0] = c0 + c4;
  coeff[7] = c1 + c5;
  coeff[3] = c2 + c6;
  coeff[4] = c3 + c7;
  coeff[2] = c0 - c4;
  coeff[6] = c1 - c5;
  coeff[1] = c2 - c6;
  coeff[5] = c3 - c7;
}
static void hadamard_col8(const int16_t *src_diff, ptrdiff_t src_stride,
                          int16_t *coeff) {
  int16_t b0 = src_diff[0 * src_stride] + src_diff[1 * src_stride];
  int16_t b1 = src_diff[0 * src_stride] - src_diff[1 * src_stride];
  int16_t b2 = src_diff[2 * src_stride] + src_diff[3 * src_stride];
  int16_t b3 = src_diff[2 * src_stride] - src_diff[3 * src_stride];
  int16_t b4 = src_diff[4 * src_stride] + src_diff[5 * src_stride];
  int16_t b5 = src_diff[4 * src_stride] - src_diff[5 * src_stride];
  int16_t b6 = src_diff[6 * src_stride] + src_diff[7 * src_stride];
  int16_t b7 = src_diff[6 * src_stride] - src_diff[7 * src_stride];

  int16_t c0 = b0 + b2;
  int16_t c1 = b1 + b3;
  int16_t c2 = b0 - b2;
  int16_t c3 = b1 - b3;
  int16_t c4 = b4 + b6;
  int16_t c5 = b5 + b7;
  int16_t c6 = b4 - b6;
  int16_t c7 = b5 - b7;

  coeff[0] = c0 + c4;
  coeff[7] = c1 + c5;
  coeff[3] = c2 + c6;
  coeff[4] = c3 + c7;
  coeff[2] = c0 - c4;
  coeff[6] = c1 - c5;
  coeff[1] = c2 - c6;
  coeff[5] = c3 - c7;
}

void aom_hadamard_8x8_using_1d(const int16_t *src_diff, ptrdiff_t src_stride,
                               tran_low_t *coeff, int *had_satd) {
  int idx;
  int16_t buffer[64];
  int16_t buffer2[64];
  int16_t *tmp_buf = &buffer[0];
  const int16_t *tmp_src = src_diff;

  // horz/row 1d-txfm
  for (idx = 0; idx < 8; ++idx) {
    hadamard_row8(tmp_src, src_stride, tmp_buf);  // src_diff: 9 bit
                                                  // dynamic range [-255, 255]
    tmp_buf += 8;
    tmp_src += src_stride;
  }
  tmp_buf = &buffer[0];
  for (int i = 0; i < 64; i++) had_satd[0] += abs(tmp_buf[i]);

  // col/vert 1d-txfm
  for (idx = 0; idx < 8; ++idx) {
    hadamard_col8(src_diff, src_stride, tmp_buf);  // src_diff: 9 bit
                                                   // dynamic range [-255, 255]
    tmp_buf += 8;
    ++src_diff;
  }
  tmp_buf = &buffer[0];
  for (int i = 0; i < 64; i++) had_satd[1] += abs(tmp_buf[i]);

  for (idx = 0; idx < 8; ++idx) {
    hadamard_col8(tmp_buf, 8,
                  buffer2 + 8 * idx);  // tmp_buf: 12 bit
                                       // dynamic range [-2040, 2040]
                                       // buffer2: 15 bit
                                       // dynamic range [-16320, 16320]
    ++tmp_buf;
  }

  for (idx = 0; idx < 64; ++idx) {
    coeff[idx] = (tran_low_t)buffer2[idx];
    had_satd[2] += abs(coeff[idx]);
  }

  had_satd[0] = (int)(had_satd[0] * 0.3535);  // multiplied by 1/(2 * sqrt(2))
  had_satd[1] = (int)(had_satd[1] * 0.3535);
  had_satd[2] = (had_satd[2] >> 3);
}
static INLINE void store_tran_low_tx_size(__m256i a, tran_low_t *b) {
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i a_hi = _mm256_mulhi_epi16(a, one);
  const __m256i a_lo = _mm256_mullo_epi16(a, one);
  const __m256i a_1 = _mm256_unpacklo_epi16(a_lo, a_hi);
  const __m256i a_2 = _mm256_unpackhi_epi16(a_lo, a_hi);
  _mm256_storeu_si256((__m256i *)b, a_1);
  _mm256_storeu_si256((__m256i *)(b + 8), a_2);
}
static void hadamard_col8x2_avx2_tx_size(__m256i *in, int iter) {
  __m256i a0 = in[0];
  __m256i a1 = in[1];
  __m256i a2 = in[2];
  __m256i a3 = in[3];
  __m256i a4 = in[4];
  __m256i a5 = in[5];
  __m256i a6 = in[6];
  __m256i a7 = in[7];

  __m256i b0 = _mm256_add_epi16(a0, a1);
  __m256i b1 = _mm256_sub_epi16(a0, a1);
  __m256i b2 = _mm256_add_epi16(a2, a3);
  __m256i b3 = _mm256_sub_epi16(a2, a3);
  __m256i b4 = _mm256_add_epi16(a4, a5);
  __m256i b5 = _mm256_sub_epi16(a4, a5);
  __m256i b6 = _mm256_add_epi16(a6, a7);
  __m256i b7 = _mm256_sub_epi16(a6, a7);

  a0 = _mm256_add_epi16(b0, b2);
  a1 = _mm256_add_epi16(b1, b3);
  a2 = _mm256_sub_epi16(b0, b2);
  a3 = _mm256_sub_epi16(b1, b3);
  a4 = _mm256_add_epi16(b4, b6);
  a5 = _mm256_add_epi16(b5, b7);
  a6 = _mm256_sub_epi16(b4, b6);
  a7 = _mm256_sub_epi16(b5, b7);

  if (iter == 0) {
    b0 = _mm256_add_epi16(a0, a4);
    b7 = _mm256_add_epi16(a1, a5);
    b3 = _mm256_add_epi16(a2, a6);
    b4 = _mm256_add_epi16(a3, a7);
    b2 = _mm256_sub_epi16(a0, a4);
    b6 = _mm256_sub_epi16(a1, a5);
    b1 = _mm256_sub_epi16(a2, a6);
    b5 = _mm256_sub_epi16(a3, a7);

    a0 = _mm256_unpacklo_epi16(b0, b1);
    a1 = _mm256_unpacklo_epi16(b2, b3);
    a2 = _mm256_unpackhi_epi16(b0, b1);
    a3 = _mm256_unpackhi_epi16(b2, b3);
    a4 = _mm256_unpacklo_epi16(b4, b5);
    a5 = _mm256_unpacklo_epi16(b6, b7);
    a6 = _mm256_unpackhi_epi16(b4, b5);
    a7 = _mm256_unpackhi_epi16(b6, b7);

    b0 = _mm256_unpacklo_epi32(a0, a1);
    b1 = _mm256_unpacklo_epi32(a4, a5);
    b2 = _mm256_unpackhi_epi32(a0, a1);
    b3 = _mm256_unpackhi_epi32(a4, a5);
    b4 = _mm256_unpacklo_epi32(a2, a3);
    b5 = _mm256_unpacklo_epi32(a6, a7);
    b6 = _mm256_unpackhi_epi32(a2, a3);
    b7 = _mm256_unpackhi_epi32(a6, a7);

    in[0] = _mm256_unpacklo_epi64(b0, b1);
    in[1] = _mm256_unpackhi_epi64(b0, b1);
    in[2] = _mm256_unpacklo_epi64(b2, b3);
    in[3] = _mm256_unpackhi_epi64(b2, b3);
    in[4] = _mm256_unpacklo_epi64(b4, b5);
    in[5] = _mm256_unpackhi_epi64(b4, b5);
    in[6] = _mm256_unpacklo_epi64(b6, b7);
    in[7] = _mm256_unpackhi_epi64(b6, b7);
  } else {
    in[0] = _mm256_add_epi16(a0, a4);
    in[7] = _mm256_add_epi16(a1, a5);
    in[3] = _mm256_add_epi16(a2, a6);
    in[4] = _mm256_add_epi16(a3, a7);
    in[2] = _mm256_sub_epi16(a0, a4);
    in[6] = _mm256_sub_epi16(a1, a5);
    in[1] = _mm256_sub_epi16(a2, a6);
    in[5] = _mm256_sub_epi16(a3, a7);

    in[0] = _mm256_srai_epi16(in[0], 3);
    in[7] = _mm256_srai_epi16(in[7], 3);
    in[3] = _mm256_srai_epi16(in[3], 3);
    in[4] = _mm256_srai_epi16(in[4], 3);
    in[2] = _mm256_srai_epi16(in[2], 3);
    in[6] = _mm256_srai_epi16(in[6], 3);
    in[1] = _mm256_srai_epi16(in[1], 3);
    in[5] = _mm256_srai_epi16(in[5], 3);
  }
}
static void hadamard_8x8x2_avx2_tx_size(const int16_t *src_diff,
                                        ptrdiff_t src_stride, int16_t *coeff) {
  __m256i src[8];
  src[0] = _mm256_loadu_si256((const __m256i *)src_diff);
  src[1] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[2] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[3] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[4] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[5] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[6] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));
  src[7] = _mm256_loadu_si256((const __m256i *)(src_diff += src_stride));

  hadamard_col8x2_avx2_tx_size(src, 0);
  hadamard_col8x2_avx2_tx_size(src, 1);

  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[0], src[1], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[2], src[3], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[4], src[5], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[6], src[7], 0x20));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[0], src[1], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[2], src[3], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[4], src[5], 0x31));
  coeff += 16;
  _mm256_storeu_si256((__m256i *)coeff,
                      _mm256_permute2x128_si256(src[6], src[7], 0x31));
}

static INLINE void hadamard_16x16_avx2_tx_size(const int16_t *src_diff,
                                               ptrdiff_t src_stride,
                                               tran_low_t *coeff,
                                               int is_final) {
  DECLARE_ALIGNED(32, int16_t, temp_coeff[16 * 16]);
  int16_t *t_coeff = temp_coeff;
  int16_t *coeff16 = (int16_t *)coeff;
  int idx;
  for (idx = 0; idx < 2; ++idx) {
    const int16_t *src_ptr = src_diff + idx * 8 * src_stride;
    hadamard_8x8x2_avx2_tx_size(src_ptr, src_stride, t_coeff + (idx * 64 * 2));
  }

  for (idx = 0; idx < 64; idx += 16) {
    const __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    const __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 64));
    const __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 128));
    const __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 192));

    __m256i b0 = _mm256_add_epi16(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi16(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi16(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi16(coeff2, coeff3);

    b0 = _mm256_srai_epi16(b0, 1);
    b1 = _mm256_srai_epi16(b1, 1);
    b2 = _mm256_srai_epi16(b2, 1);
    b3 = _mm256_srai_epi16(b3, 1);
    if (is_final) {
      store_tran_low_tx_size(_mm256_add_epi16(b0, b2), coeff);
      store_tran_low_tx_size(_mm256_add_epi16(b1, b3), coeff + 64);
      store_tran_low_tx_size(_mm256_sub_epi16(b0, b2), coeff + 128);
      store_tran_low_tx_size(_mm256_sub_epi16(b1, b3), coeff + 192);
      coeff += 16;
    } else {
      _mm256_storeu_si256((__m256i *)coeff16, _mm256_add_epi16(b0, b2));
      _mm256_storeu_si256((__m256i *)(coeff16 + 64), _mm256_add_epi16(b1, b3));
      _mm256_storeu_si256((__m256i *)(coeff16 + 128), _mm256_sub_epi16(b0, b2));
      _mm256_storeu_si256((__m256i *)(coeff16 + 192), _mm256_sub_epi16(b1, b3));
      coeff16 += 16;
    }
    t_coeff += 16;
  }
}

void aom_hadamard_32x32_avx2_tx_size(const int16_t *src_diff,
                                     ptrdiff_t src_stride, tran_low_t *coeff) {
  // For high bitdepths, it is unnecessary to store_tran_low
  // (mult/unpack/store), then load_tran_low (load/pack) the same memory in the
  // next stage.  Output to an intermediate buffer first, then store_tran_low()
  // in the final stage.
  DECLARE_ALIGNED(32, int16_t, temp_coeff[32 * 32]);
  int16_t *t_coeff = temp_coeff;
  int idx;
  for (idx = 0; idx < 4; ++idx) {
    // src_diff: 9 bit, dynamic range [-255, 255]
    const int16_t *src_ptr =
        src_diff + (idx >> 1) * 16 * src_stride + (idx & 0x01) * 16;
    hadamard_16x16_avx2_tx_size(src_ptr, src_stride,
                                (tran_low_t *)(t_coeff + idx * 256), 0);
  }

  for (idx = 0; idx < 256; idx += 16) {
    const __m256i coeff0 = _mm256_loadu_si256((const __m256i *)t_coeff);
    const __m256i coeff1 = _mm256_loadu_si256((const __m256i *)(t_coeff + 256));
    const __m256i coeff2 = _mm256_loadu_si256((const __m256i *)(t_coeff + 512));
    const __m256i coeff3 = _mm256_loadu_si256((const __m256i *)(t_coeff + 768));

    __m256i b0 = _mm256_add_epi16(coeff0, coeff1);
    __m256i b1 = _mm256_sub_epi16(coeff0, coeff1);
    __m256i b2 = _mm256_add_epi16(coeff2, coeff3);
    __m256i b3 = _mm256_sub_epi16(coeff2, coeff3);

    b0 = _mm256_srai_epi16(b0, 1);
    b1 = _mm256_srai_epi16(b1, 1);
    b2 = _mm256_srai_epi16(b2, 1);
    b3 = _mm256_srai_epi16(b3, 1);

    store_tran_low_tx_size(_mm256_add_epi16(b0, b2), coeff);
    store_tran_low_tx_size(_mm256_add_epi16(b1, b3), coeff + 256);
    store_tran_low_tx_size(_mm256_sub_epi16(b0, b2), coeff + 512);
    store_tran_low_tx_size(_mm256_sub_epi16(b1, b3), coeff + 768);

    coeff += 16;
    t_coeff += 16;
  }
}
