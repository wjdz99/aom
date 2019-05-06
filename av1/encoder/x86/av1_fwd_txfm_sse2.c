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

#include "av1/common/x86/av1_txfm_sse2.h"
#include "av1/encoder/av1_fwd_txfm1d_cfg.h"
#include "av1/encoder/x86/av1_fwd_txfm_sse2.h"

// TODO(linfengz): refine fdct4x8 and fadst4x8 optimization (if possible).

static void fdct4x4_new_sse2(const __m128i *input, __m128i *output,
                             int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  const __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  const __m128i cospi_p16_p48 = pair_set_epi16(cospi[16], cospi[48]);
  const __m128i cospi_p48_m16 = pair_set_epi16(cospi[48], -cospi[16]);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));
  __m128i u[4], v[4];

  u[0] = _mm_unpacklo_epi16(input[0], input[1]);
  u[1] = _mm_unpacklo_epi16(input[3], input[2]);

  v[0] = _mm_add_epi16(u[0], u[1]);
  v[1] = _mm_sub_epi16(u[0], u[1]);

  u[0] = _mm_madd_epi16(v[0], cospi_p32_p32);  // 0
  u[1] = _mm_madd_epi16(v[0], cospi_p32_m32);  // 2
  u[2] = _mm_madd_epi16(v[1], cospi_p16_p48);  // 1
  u[3] = _mm_madd_epi16(v[1], cospi_p48_m16);  // 3

  v[0] = _mm_add_epi32(u[0], __rounding);
  v[1] = _mm_add_epi32(u[1], __rounding);
  v[2] = _mm_add_epi32(u[2], __rounding);
  v[3] = _mm_add_epi32(u[3], __rounding);
  u[0] = _mm_srai_epi32(v[0], cos_bit);
  u[1] = _mm_srai_epi32(v[1], cos_bit);
  u[2] = _mm_srai_epi32(v[2], cos_bit);
  u[3] = _mm_srai_epi32(v[3], cos_bit);

  output[0] = _mm_packs_epi32(u[0], u[1]);
  output[1] = _mm_packs_epi32(u[2], u[3]);
  output[2] = _mm_srli_si128(output[0], 8);
  output[3] = _mm_srli_si128(output[1], 8);
}

static void fdct8x4_new_sse2(const __m128i *input, __m128i *output,
                             int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p48_p16 = pair_set_epi16(cospi[48], cospi[16]);
  __m128i cospi_m16_p48 = pair_set_epi16(-cospi[16], cospi[48]);

  // stage 1
  __m128i x1[4];
  x1[0] = _mm_adds_epi16(input[0], input[3]);
  x1[3] = _mm_subs_epi16(input[0], input[3]);
  x1[1] = _mm_adds_epi16(input[1], input[2]);
  x1[2] = _mm_subs_epi16(input[1], input[2]);

  // stage 2
  __m128i x2[4];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[0], x1[1], x2[0], x2[1]);
  btf_16_sse2(cospi_p48_p16, cospi_m16_p48, x1[2], x1[3], x2[2], x2[3]);

  // stage 3
  output[0] = x2[0];
  output[1] = x2[2];
  output[2] = x2[1];
  output[3] = x2[3];
}

static void fdct4x8_new_sse2(const __m128i *input, __m128i *output,
                             int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_m32_p32 = pair_set_epi16(-cospi[32], cospi[32]);
  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p48_p16 = pair_set_epi16(cospi[48], cospi[16]);
  __m128i cospi_m16_p48 = pair_set_epi16(-cospi[16], cospi[48]);
  __m128i cospi_p56_p08 = pair_set_epi16(cospi[56], cospi[8]);
  __m128i cospi_m08_p56 = pair_set_epi16(-cospi[8], cospi[56]);
  __m128i cospi_p24_p40 = pair_set_epi16(cospi[24], cospi[40]);
  __m128i cospi_m40_p24 = pair_set_epi16(-cospi[40], cospi[24]);

  // stage 1
  __m128i x1[8];
  x1[0] = _mm_adds_epi16(input[0], input[7]);
  x1[7] = _mm_subs_epi16(input[0], input[7]);
  x1[1] = _mm_adds_epi16(input[1], input[6]);
  x1[6] = _mm_subs_epi16(input[1], input[6]);
  x1[2] = _mm_adds_epi16(input[2], input[5]);
  x1[5] = _mm_subs_epi16(input[2], input[5]);
  x1[3] = _mm_adds_epi16(input[3], input[4]);
  x1[4] = _mm_subs_epi16(input[3], input[4]);

  // stage 2
  __m128i x2[8];
  x2[0] = _mm_adds_epi16(x1[0], x1[3]);
  x2[3] = _mm_subs_epi16(x1[0], x1[3]);
  x2[1] = _mm_adds_epi16(x1[1], x1[2]);
  x2[2] = _mm_subs_epi16(x1[1], x1[2]);
  x2[4] = x1[4];
  btf_16_w4_sse2(&cospi_m32_p32, &cospi_p32_p32, __rounding, cos_bit, &x1[5],
                 &x1[6], &x2[5], &x2[6]);
  x2[7] = x1[7];

  // stage 3
  __m128i x3[8];
  btf_16_w4_sse2(&cospi_p32_p32, &cospi_p32_m32, __rounding, cos_bit, &x2[0],
                 &x2[1], &x3[0], &x3[1]);
  btf_16_w4_sse2(&cospi_p48_p16, &cospi_m16_p48, __rounding, cos_bit, &x2[2],
                 &x2[3], &x3[2], &x3[3]);
  x3[4] = _mm_adds_epi16(x2[4], x2[5]);
  x3[5] = _mm_subs_epi16(x2[4], x2[5]);
  x3[6] = _mm_subs_epi16(x2[7], x2[6]);
  x3[7] = _mm_adds_epi16(x2[7], x2[6]);

  // stage 4
  __m128i x4[8];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_w4_sse2(&cospi_p56_p08, &cospi_m08_p56, __rounding, cos_bit, &x3[4],
                 &x3[7], &x4[4], &x4[7]);
  btf_16_w4_sse2(&cospi_p24_p40, &cospi_m40_p24, __rounding, cos_bit, &x3[5],
                 &x3[6], &x4[5], &x4[6]);

  // stage 5
  output[0] = x4[0];
  output[1] = x4[4];
  output[2] = x4[2];
  output[3] = x4[6];
  output[4] = x4[1];
  output[5] = x4[5];
  output[6] = x4[3];
  output[7] = x4[7];
}

static void fdct8x8_new_sse2(const __m128i *input, __m128i *output,
                             int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_m32_p32 = pair_set_epi16(-cospi[32], cospi[32]);
  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p48_p16 = pair_set_epi16(cospi[48], cospi[16]);
  __m128i cospi_m16_p48 = pair_set_epi16(-cospi[16], cospi[48]);
  __m128i cospi_p56_p08 = pair_set_epi16(cospi[56], cospi[8]);
  __m128i cospi_m08_p56 = pair_set_epi16(-cospi[8], cospi[56]);
  __m128i cospi_p24_p40 = pair_set_epi16(cospi[24], cospi[40]);
  __m128i cospi_m40_p24 = pair_set_epi16(-cospi[40], cospi[24]);

  // stage 1
  __m128i x1[8];
  x1[0] = _mm_adds_epi16(input[0], input[7]);
  x1[7] = _mm_subs_epi16(input[0], input[7]);
  x1[1] = _mm_adds_epi16(input[1], input[6]);
  x1[6] = _mm_subs_epi16(input[1], input[6]);
  x1[2] = _mm_adds_epi16(input[2], input[5]);
  x1[5] = _mm_subs_epi16(input[2], input[5]);
  x1[3] = _mm_adds_epi16(input[3], input[4]);
  x1[4] = _mm_subs_epi16(input[3], input[4]);

  // stage 2
  __m128i x2[8];
  x2[0] = _mm_adds_epi16(x1[0], x1[3]);
  x2[3] = _mm_subs_epi16(x1[0], x1[3]);
  x2[1] = _mm_adds_epi16(x1[1], x1[2]);
  x2[2] = _mm_subs_epi16(x1[1], x1[2]);
  x2[4] = x1[4];
  btf_16_sse2(cospi_m32_p32, cospi_p32_p32, x1[5], x1[6], x2[5], x2[6]);
  x2[7] = x1[7];

  // stage 3
  __m128i x3[8];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x2[0], x2[1], x3[0], x3[1]);
  btf_16_sse2(cospi_p48_p16, cospi_m16_p48, x2[2], x2[3], x3[2], x3[3]);
  x3[4] = _mm_adds_epi16(x2[4], x2[5]);
  x3[5] = _mm_subs_epi16(x2[4], x2[5]);
  x3[6] = _mm_subs_epi16(x2[7], x2[6]);
  x3[7] = _mm_adds_epi16(x2[7], x2[6]);

  // stage 4
  __m128i x4[8];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_sse2(cospi_p56_p08, cospi_m08_p56, x3[4], x3[7], x4[4], x4[7]);
  btf_16_sse2(cospi_p24_p40, cospi_m40_p24, x3[5], x3[6], x4[5], x4[6]);

  // stage 5
  output[0] = x4[0];
  output[1] = x4[4];
  output[2] = x4[2];
  output[3] = x4[6];
  output[4] = x4[1];
  output[5] = x4[5];
  output[6] = x4[3];
  output[7] = x4[7];
}

static void fdct8x16_new_sse2(const __m128i *input, __m128i *output,
                              int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_m32_p32 = pair_set_epi16(-cospi[32], cospi[32]);
  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p48_p16 = pair_set_epi16(cospi[48], cospi[16]);
  __m128i cospi_m16_p48 = pair_set_epi16(-cospi[16], cospi[48]);
  __m128i cospi_m48_m16 = pair_set_epi16(-cospi[48], -cospi[16]);
  __m128i cospi_p56_p08 = pair_set_epi16(cospi[56], cospi[8]);
  __m128i cospi_m08_p56 = pair_set_epi16(-cospi[8], cospi[56]);
  __m128i cospi_p24_p40 = pair_set_epi16(cospi[24], cospi[40]);
  __m128i cospi_m40_p24 = pair_set_epi16(-cospi[40], cospi[24]);
  __m128i cospi_p60_p04 = pair_set_epi16(cospi[60], cospi[4]);
  __m128i cospi_m04_p60 = pair_set_epi16(-cospi[4], cospi[60]);
  __m128i cospi_p28_p36 = pair_set_epi16(cospi[28], cospi[36]);
  __m128i cospi_m36_p28 = pair_set_epi16(-cospi[36], cospi[28]);
  __m128i cospi_p44_p20 = pair_set_epi16(cospi[44], cospi[20]);
  __m128i cospi_m20_p44 = pair_set_epi16(-cospi[20], cospi[44]);
  __m128i cospi_p12_p52 = pair_set_epi16(cospi[12], cospi[52]);
  __m128i cospi_m52_p12 = pair_set_epi16(-cospi[52], cospi[12]);

  // stage 1
  __m128i x1[16];
  x1[0] = _mm_adds_epi16(input[0], input[15]);
  x1[15] = _mm_subs_epi16(input[0], input[15]);
  x1[1] = _mm_adds_epi16(input[1], input[14]);
  x1[14] = _mm_subs_epi16(input[1], input[14]);
  x1[2] = _mm_adds_epi16(input[2], input[13]);
  x1[13] = _mm_subs_epi16(input[2], input[13]);
  x1[3] = _mm_adds_epi16(input[3], input[12]);
  x1[12] = _mm_subs_epi16(input[3], input[12]);
  x1[4] = _mm_adds_epi16(input[4], input[11]);
  x1[11] = _mm_subs_epi16(input[4], input[11]);
  x1[5] = _mm_adds_epi16(input[5], input[10]);
  x1[10] = _mm_subs_epi16(input[5], input[10]);
  x1[6] = _mm_adds_epi16(input[6], input[9]);
  x1[9] = _mm_subs_epi16(input[6], input[9]);
  x1[7] = _mm_adds_epi16(input[7], input[8]);
  x1[8] = _mm_subs_epi16(input[7], input[8]);

  // stage 2
  __m128i x2[16];
  x2[0] = _mm_adds_epi16(x1[0], x1[7]);
  x2[7] = _mm_subs_epi16(x1[0], x1[7]);
  x2[1] = _mm_adds_epi16(x1[1], x1[6]);
  x2[6] = _mm_subs_epi16(x1[1], x1[6]);
  x2[2] = _mm_adds_epi16(x1[2], x1[5]);
  x2[5] = _mm_subs_epi16(x1[2], x1[5]);
  x2[3] = _mm_adds_epi16(x1[3], x1[4]);
  x2[4] = _mm_subs_epi16(x1[3], x1[4]);
  x2[8] = x1[8];
  x2[9] = x1[9];
  btf_16_sse2(cospi_m32_p32, cospi_p32_p32, x1[10], x1[13], x2[10], x2[13]);
  btf_16_sse2(cospi_m32_p32, cospi_p32_p32, x1[11], x1[12], x2[11], x2[12]);
  x2[14] = x1[14];
  x2[15] = x1[15];

  // stage 3
  __m128i x3[16];
  x3[0] = _mm_adds_epi16(x2[0], x2[3]);
  x3[3] = _mm_subs_epi16(x2[0], x2[3]);
  x3[1] = _mm_adds_epi16(x2[1], x2[2]);
  x3[2] = _mm_subs_epi16(x2[1], x2[2]);
  x3[4] = x2[4];
  btf_16_sse2(cospi_m32_p32, cospi_p32_p32, x2[5], x2[6], x3[5], x3[6]);
  x3[7] = x2[7];
  x3[8] = _mm_adds_epi16(x2[8], x2[11]);
  x3[11] = _mm_subs_epi16(x2[8], x2[11]);
  x3[9] = _mm_adds_epi16(x2[9], x2[10]);
  x3[10] = _mm_subs_epi16(x2[9], x2[10]);
  x3[12] = _mm_subs_epi16(x2[15], x2[12]);
  x3[15] = _mm_adds_epi16(x2[15], x2[12]);
  x3[13] = _mm_subs_epi16(x2[14], x2[13]);
  x3[14] = _mm_adds_epi16(x2[14], x2[13]);

  // stage 4
  __m128i x4[16];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x3[0], x3[1], x4[0], x4[1]);
  btf_16_sse2(cospi_p48_p16, cospi_m16_p48, x3[2], x3[3], x4[2], x4[3]);
  x4[4] = _mm_adds_epi16(x3[4], x3[5]);
  x4[5] = _mm_subs_epi16(x3[4], x3[5]);
  x4[6] = _mm_subs_epi16(x3[7], x3[6]);
  x4[7] = _mm_adds_epi16(x3[7], x3[6]);
  x4[8] = x3[8];
  btf_16_sse2(cospi_m16_p48, cospi_p48_p16, x3[9], x3[14], x4[9], x4[14]);
  btf_16_sse2(cospi_m48_m16, cospi_m16_p48, x3[10], x3[13], x4[10], x4[13]);
  x4[11] = x3[11];
  x4[12] = x3[12];
  x4[15] = x3[15];

  // stage 5
  __m128i x5[16];
  x5[0] = x4[0];
  x5[1] = x4[1];
  x5[2] = x4[2];
  x5[3] = x4[3];
  btf_16_sse2(cospi_p56_p08, cospi_m08_p56, x4[4], x4[7], x5[4], x5[7]);
  btf_16_sse2(cospi_p24_p40, cospi_m40_p24, x4[5], x4[6], x5[5], x5[6]);
  x5[8] = _mm_adds_epi16(x4[8], x4[9]);
  x5[9] = _mm_subs_epi16(x4[8], x4[9]);
  x5[10] = _mm_subs_epi16(x4[11], x4[10]);
  x5[11] = _mm_adds_epi16(x4[11], x4[10]);
  x5[12] = _mm_adds_epi16(x4[12], x4[13]);
  x5[13] = _mm_subs_epi16(x4[12], x4[13]);
  x5[14] = _mm_subs_epi16(x4[15], x4[14]);
  x5[15] = _mm_adds_epi16(x4[15], x4[14]);

  // stage 6
  __m128i x6[16];
  x6[0] = x5[0];
  x6[1] = x5[1];
  x6[2] = x5[2];
  x6[3] = x5[3];
  x6[4] = x5[4];
  x6[5] = x5[5];
  x6[6] = x5[6];
  x6[7] = x5[7];
  btf_16_sse2(cospi_p60_p04, cospi_m04_p60, x5[8], x5[15], x6[8], x6[15]);
  btf_16_sse2(cospi_p28_p36, cospi_m36_p28, x5[9], x5[14], x6[9], x6[14]);
  btf_16_sse2(cospi_p44_p20, cospi_m20_p44, x5[10], x5[13], x6[10], x6[13]);
  btf_16_sse2(cospi_p12_p52, cospi_m52_p12, x5[11], x5[12], x6[11], x6[12]);

  // stage 7
  output[0] = x6[0];
  output[1] = x6[8];
  output[2] = x6[4];
  output[3] = x6[12];
  output[4] = x6[2];
  output[5] = x6[10];
  output[6] = x6[6];
  output[7] = x6[14];
  output[8] = x6[1];
  output[9] = x6[9];
  output[10] = x6[5];
  output[11] = x6[13];
  output[12] = x6[3];
  output[13] = x6[11];
  output[14] = x6[7];
  output[15] = x6[15];
}

static void fadst4x4_new_sse2(const __m128i *input, __m128i *output,
                              int8_t cos_bit) {
  const int32_t *sinpi = sinpi_arr(cos_bit);
  const __m128i sinpi_p01_p02 = pair_set_epi16(sinpi[1], sinpi[2]);
  const __m128i sinpi_p04_m01 = pair_set_epi16(sinpi[4], -sinpi[1]);
  const __m128i sinpi_p03_p04 = pair_set_epi16(sinpi[3], sinpi[4]);
  const __m128i sinpi_m03_p02 = pair_set_epi16(-sinpi[3], sinpi[2]);
  const __m128i sinpi_p03_p03 = _mm_set1_epi16((int16_t)sinpi[3]);
  const __m128i __zero = _mm_set1_epi16(0);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));
  const __m128i in7 = _mm_add_epi16(input[0], input[1]);
  __m128i u[8], v[8];

  u[0] = _mm_unpacklo_epi16(input[0], input[1]);
  u[1] = _mm_unpacklo_epi16(input[2], input[3]);
  u[2] = _mm_unpacklo_epi16(in7, __zero);
  u[3] = _mm_unpacklo_epi16(input[2], __zero);
  u[4] = _mm_unpacklo_epi16(input[3], __zero);

  v[0] = _mm_madd_epi16(u[0], sinpi_p01_p02);  // s0 + s2
  v[1] = _mm_madd_epi16(u[1], sinpi_p03_p04);  // s4 + s5
  v[2] = _mm_madd_epi16(u[2], sinpi_p03_p03);  // x1
  v[3] = _mm_madd_epi16(u[0], sinpi_p04_m01);  // s1 - s3
  v[4] = _mm_madd_epi16(u[1], sinpi_m03_p02);  // -s4 + s6
  v[5] = _mm_madd_epi16(u[3], sinpi_p03_p03);  // s4
  v[6] = _mm_madd_epi16(u[4], sinpi_p03_p03);

  u[0] = _mm_add_epi32(v[0], v[1]);
  u[1] = _mm_sub_epi32(v[2], v[6]);
  u[2] = _mm_add_epi32(v[3], v[4]);
  u[3] = _mm_sub_epi32(u[2], u[0]);
  u[4] = _mm_slli_epi32(v[5], 2);
  u[5] = _mm_sub_epi32(u[4], v[5]);
  u[6] = _mm_add_epi32(u[3], u[5]);

  v[0] = _mm_add_epi32(u[0], __rounding);
  v[1] = _mm_add_epi32(u[1], __rounding);
  v[2] = _mm_add_epi32(u[2], __rounding);
  v[3] = _mm_add_epi32(u[6], __rounding);

  u[0] = _mm_srai_epi32(v[0], cos_bit);
  u[1] = _mm_srai_epi32(v[1], cos_bit);
  u[2] = _mm_srai_epi32(v[2], cos_bit);
  u[3] = _mm_srai_epi32(v[3], cos_bit);

  output[0] = _mm_packs_epi32(u[0], u[2]);
  output[1] = _mm_packs_epi32(u[1], u[3]);
  output[2] = _mm_srli_si128(output[0], 8);
  output[3] = _mm_srli_si128(output[1], 8);
}

static void fadst4x8_new_sse2(const __m128i *input, __m128i *output,
                              int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __zero = _mm_setzero_si128();
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p16_p48 = pair_set_epi16(cospi[16], cospi[48]);
  __m128i cospi_p48_m16 = pair_set_epi16(cospi[48], -cospi[16]);
  __m128i cospi_m48_p16 = pair_set_epi16(-cospi[48], cospi[16]);
  __m128i cospi_p04_p60 = pair_set_epi16(cospi[4], cospi[60]);
  __m128i cospi_p60_m04 = pair_set_epi16(cospi[60], -cospi[4]);
  __m128i cospi_p20_p44 = pair_set_epi16(cospi[20], cospi[44]);
  __m128i cospi_p44_m20 = pair_set_epi16(cospi[44], -cospi[20]);
  __m128i cospi_p36_p28 = pair_set_epi16(cospi[36], cospi[28]);
  __m128i cospi_p28_m36 = pair_set_epi16(cospi[28], -cospi[36]);
  __m128i cospi_p52_p12 = pair_set_epi16(cospi[52], cospi[12]);
  __m128i cospi_p12_m52 = pair_set_epi16(cospi[12], -cospi[52]);

  // stage 1
  __m128i x1[8];
  x1[0] = input[0];
  x1[1] = _mm_subs_epi16(__zero, input[7]);
  x1[2] = _mm_subs_epi16(__zero, input[3]);
  x1[3] = input[4];
  x1[4] = _mm_subs_epi16(__zero, input[1]);
  x1[5] = input[6];
  x1[6] = input[2];
  x1[7] = _mm_subs_epi16(__zero, input[5]);

  // stage 2
  __m128i x2[8];
  x2[0] = x1[0];
  x2[1] = x1[1];
  btf_16_w4_sse2(&cospi_p32_p32, &cospi_p32_m32, __rounding, cos_bit, &x1[2],
                 &x1[3], &x2[2], &x2[3]);
  x2[4] = x1[4];
  x2[5] = x1[5];
  btf_16_w4_sse2(&cospi_p32_p32, &cospi_p32_m32, __rounding, cos_bit, &x1[6],
                 &x1[7], &x2[6], &x2[7]);

  // stage 3
  __m128i x3[8];
  x3[0] = _mm_adds_epi16(x2[0], x2[2]);
  x3[2] = _mm_subs_epi16(x2[0], x2[2]);
  x3[1] = _mm_adds_epi16(x2[1], x2[3]);
  x3[3] = _mm_subs_epi16(x2[1], x2[3]);
  x3[4] = _mm_adds_epi16(x2[4], x2[6]);
  x3[6] = _mm_subs_epi16(x2[4], x2[6]);
  x3[5] = _mm_adds_epi16(x2[5], x2[7]);
  x3[7] = _mm_subs_epi16(x2[5], x2[7]);

  // stage 4
  __m128i x4[8];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_w4_sse2(&cospi_p16_p48, &cospi_p48_m16, __rounding, cos_bit, &x3[4],
                 &x3[5], &x4[4], &x4[5]);
  btf_16_w4_sse2(&cospi_m48_p16, &cospi_p16_p48, __rounding, cos_bit, &x3[6],
                 &x3[7], &x4[6], &x4[7]);

  // stage 5
  __m128i x5[8];
  x5[0] = _mm_adds_epi16(x4[0], x4[4]);
  x5[4] = _mm_subs_epi16(x4[0], x4[4]);
  x5[1] = _mm_adds_epi16(x4[1], x4[5]);
  x5[5] = _mm_subs_epi16(x4[1], x4[5]);
  x5[2] = _mm_adds_epi16(x4[2], x4[6]);
  x5[6] = _mm_subs_epi16(x4[2], x4[6]);
  x5[3] = _mm_adds_epi16(x4[3], x4[7]);
  x5[7] = _mm_subs_epi16(x4[3], x4[7]);

  // stage 6
  __m128i x6[8];
  btf_16_w4_sse2(&cospi_p04_p60, &cospi_p60_m04, __rounding, cos_bit, &x5[0],
                 &x5[1], &x6[0], &x6[1]);
  btf_16_w4_sse2(&cospi_p20_p44, &cospi_p44_m20, __rounding, cos_bit, &x5[2],
                 &x5[3], &x6[2], &x6[3]);
  btf_16_w4_sse2(&cospi_p36_p28, &cospi_p28_m36, __rounding, cos_bit, &x5[4],
                 &x5[5], &x6[4], &x6[5]);
  btf_16_w4_sse2(&cospi_p52_p12, &cospi_p12_m52, __rounding, cos_bit, &x5[6],
                 &x5[7], &x6[6], &x6[7]);

  // stage 7
  output[0] = x6[1];
  output[1] = x6[6];
  output[2] = x6[3];
  output[3] = x6[4];
  output[4] = x6[5];
  output[5] = x6[2];
  output[6] = x6[7];
  output[7] = x6[0];
}

static void fadst8x4_new_sse2(const __m128i *input, __m128i *output,
                              int8_t cos_bit) {
  const int32_t *sinpi = sinpi_arr(cos_bit);
  const __m128i sinpi_p01_p02 = pair_set_epi16(sinpi[1], sinpi[2]);
  const __m128i sinpi_p04_m01 = pair_set_epi16(sinpi[4], -sinpi[1]);
  const __m128i sinpi_p03_p04 = pair_set_epi16(sinpi[3], sinpi[4]);
  const __m128i sinpi_m03_p02 = pair_set_epi16(-sinpi[3], sinpi[2]);
  const __m128i sinpi_p03_p03 = _mm_set1_epi16((int16_t)sinpi[3]);
  const __m128i __zero = _mm_set1_epi16(0);
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));
  const __m128i in7 = _mm_add_epi16(input[0], input[1]);
  __m128i u_lo[8], u_hi[8], v_lo[8], v_hi[8];

  u_lo[0] = _mm_unpacklo_epi16(input[0], input[1]);
  u_hi[0] = _mm_unpackhi_epi16(input[0], input[1]);
  u_lo[1] = _mm_unpacklo_epi16(input[2], input[3]);
  u_hi[1] = _mm_unpackhi_epi16(input[2], input[3]);
  u_lo[2] = _mm_unpacklo_epi16(in7, __zero);
  u_hi[2] = _mm_unpackhi_epi16(in7, __zero);
  u_lo[3] = _mm_unpacklo_epi16(input[2], __zero);
  u_hi[3] = _mm_unpackhi_epi16(input[2], __zero);
  u_lo[4] = _mm_unpacklo_epi16(input[3], __zero);
  u_hi[4] = _mm_unpackhi_epi16(input[3], __zero);

  v_lo[0] = _mm_madd_epi16(u_lo[0], sinpi_p01_p02);  // s0 + s2
  v_hi[0] = _mm_madd_epi16(u_hi[0], sinpi_p01_p02);  // s0 + s2
  v_lo[1] = _mm_madd_epi16(u_lo[1], sinpi_p03_p04);  // s4 + s5
  v_hi[1] = _mm_madd_epi16(u_hi[1], sinpi_p03_p04);  // s4 + s5
  v_lo[2] = _mm_madd_epi16(u_lo[2], sinpi_p03_p03);  // x1
  v_hi[2] = _mm_madd_epi16(u_hi[2], sinpi_p03_p03);  // x1
  v_lo[3] = _mm_madd_epi16(u_lo[0], sinpi_p04_m01);  // s1 - s3
  v_hi[3] = _mm_madd_epi16(u_hi[0], sinpi_p04_m01);  // s1 - s3
  v_lo[4] = _mm_madd_epi16(u_lo[1], sinpi_m03_p02);  // -s4 + s6
  v_hi[4] = _mm_madd_epi16(u_hi[1], sinpi_m03_p02);  // -s4 + s6
  v_lo[5] = _mm_madd_epi16(u_lo[3], sinpi_p03_p03);  // s4
  v_hi[5] = _mm_madd_epi16(u_hi[3], sinpi_p03_p03);  // s4
  v_lo[6] = _mm_madd_epi16(u_lo[4], sinpi_p03_p03);
  v_hi[6] = _mm_madd_epi16(u_hi[4], sinpi_p03_p03);

  u_lo[0] = _mm_add_epi32(v_lo[0], v_lo[1]);
  u_hi[0] = _mm_add_epi32(v_hi[0], v_hi[1]);
  u_lo[1] = _mm_sub_epi32(v_lo[2], v_lo[6]);
  u_hi[1] = _mm_sub_epi32(v_hi[2], v_hi[6]);
  u_lo[2] = _mm_add_epi32(v_lo[3], v_lo[4]);
  u_hi[2] = _mm_add_epi32(v_hi[3], v_hi[4]);
  u_lo[3] = _mm_sub_epi32(u_lo[2], u_lo[0]);
  u_hi[3] = _mm_sub_epi32(u_hi[2], u_hi[0]);
  u_lo[4] = _mm_slli_epi32(v_lo[5], 2);
  u_hi[4] = _mm_slli_epi32(v_hi[5], 2);
  u_lo[5] = _mm_sub_epi32(u_lo[4], v_lo[5]);
  u_hi[5] = _mm_sub_epi32(u_hi[4], v_hi[5]);
  u_lo[6] = _mm_add_epi32(u_lo[3], u_lo[5]);
  u_hi[6] = _mm_add_epi32(u_hi[3], u_hi[5]);

  v_lo[0] = _mm_add_epi32(u_lo[0], __rounding);
  v_hi[0] = _mm_add_epi32(u_hi[0], __rounding);
  v_lo[1] = _mm_add_epi32(u_lo[1], __rounding);
  v_hi[1] = _mm_add_epi32(u_hi[1], __rounding);
  v_lo[2] = _mm_add_epi32(u_lo[2], __rounding);
  v_hi[2] = _mm_add_epi32(u_hi[2], __rounding);
  v_lo[3] = _mm_add_epi32(u_lo[6], __rounding);
  v_hi[3] = _mm_add_epi32(u_hi[6], __rounding);

  u_lo[0] = _mm_srai_epi32(v_lo[0], cos_bit);
  u_hi[0] = _mm_srai_epi32(v_hi[0], cos_bit);
  u_lo[1] = _mm_srai_epi32(v_lo[1], cos_bit);
  u_hi[1] = _mm_srai_epi32(v_hi[1], cos_bit);
  u_lo[2] = _mm_srai_epi32(v_lo[2], cos_bit);
  u_hi[2] = _mm_srai_epi32(v_hi[2], cos_bit);
  u_lo[3] = _mm_srai_epi32(v_lo[3], cos_bit);
  u_hi[3] = _mm_srai_epi32(v_hi[3], cos_bit);

  output[0] = _mm_packs_epi32(u_lo[0], u_hi[0]);
  output[1] = _mm_packs_epi32(u_lo[1], u_hi[1]);
  output[2] = _mm_packs_epi32(u_lo[2], u_hi[2]);
  output[3] = _mm_packs_epi32(u_lo[3], u_hi[3]);
}

static void fadst8x8_new_sse2(const __m128i *input, __m128i *output,
                              int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __zero = _mm_setzero_si128();
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p16_p48 = pair_set_epi16(cospi[16], cospi[48]);
  __m128i cospi_p48_m16 = pair_set_epi16(cospi[48], -cospi[16]);
  __m128i cospi_m48_p16 = pair_set_epi16(-cospi[48], cospi[16]);
  __m128i cospi_p04_p60 = pair_set_epi16(cospi[4], cospi[60]);
  __m128i cospi_p60_m04 = pair_set_epi16(cospi[60], -cospi[4]);
  __m128i cospi_p20_p44 = pair_set_epi16(cospi[20], cospi[44]);
  __m128i cospi_p44_m20 = pair_set_epi16(cospi[44], -cospi[20]);
  __m128i cospi_p36_p28 = pair_set_epi16(cospi[36], cospi[28]);
  __m128i cospi_p28_m36 = pair_set_epi16(cospi[28], -cospi[36]);
  __m128i cospi_p52_p12 = pair_set_epi16(cospi[52], cospi[12]);
  __m128i cospi_p12_m52 = pair_set_epi16(cospi[12], -cospi[52]);

  // stage 1
  __m128i x1[8];
  x1[0] = input[0];
  x1[1] = _mm_subs_epi16(__zero, input[7]);
  x1[2] = _mm_subs_epi16(__zero, input[3]);
  x1[3] = input[4];
  x1[4] = _mm_subs_epi16(__zero, input[1]);
  x1[5] = input[6];
  x1[6] = input[2];
  x1[7] = _mm_subs_epi16(__zero, input[5]);

  // stage 2
  __m128i x2[8];
  x2[0] = x1[0];
  x2[1] = x1[1];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[2], x1[3], x2[2], x2[3]);
  x2[4] = x1[4];
  x2[5] = x1[5];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[6], x1[7], x2[6], x2[7]);

  // stage 3
  __m128i x3[8];
  x3[0] = _mm_adds_epi16(x2[0], x2[2]);
  x3[2] = _mm_subs_epi16(x2[0], x2[2]);
  x3[1] = _mm_adds_epi16(x2[1], x2[3]);
  x3[3] = _mm_subs_epi16(x2[1], x2[3]);
  x3[4] = _mm_adds_epi16(x2[4], x2[6]);
  x3[6] = _mm_subs_epi16(x2[4], x2[6]);
  x3[5] = _mm_adds_epi16(x2[5], x2[7]);
  x3[7] = _mm_subs_epi16(x2[5], x2[7]);

  // stage 4
  __m128i x4[8];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_sse2(cospi_p16_p48, cospi_p48_m16, x3[4], x3[5], x4[4], x4[5]);
  btf_16_sse2(cospi_m48_p16, cospi_p16_p48, x3[6], x3[7], x4[6], x4[7]);

  // stage 5
  __m128i x5[8];
  x5[0] = _mm_adds_epi16(x4[0], x4[4]);
  x5[4] = _mm_subs_epi16(x4[0], x4[4]);
  x5[1] = _mm_adds_epi16(x4[1], x4[5]);
  x5[5] = _mm_subs_epi16(x4[1], x4[5]);
  x5[2] = _mm_adds_epi16(x4[2], x4[6]);
  x5[6] = _mm_subs_epi16(x4[2], x4[6]);
  x5[3] = _mm_adds_epi16(x4[3], x4[7]);
  x5[7] = _mm_subs_epi16(x4[3], x4[7]);

  // stage 6
  __m128i x6[8];
  btf_16_sse2(cospi_p04_p60, cospi_p60_m04, x5[0], x5[1], x6[0], x6[1]);
  btf_16_sse2(cospi_p20_p44, cospi_p44_m20, x5[2], x5[3], x6[2], x6[3]);
  btf_16_sse2(cospi_p36_p28, cospi_p28_m36, x5[4], x5[5], x6[4], x6[5]);
  btf_16_sse2(cospi_p52_p12, cospi_p12_m52, x5[6], x5[7], x6[6], x6[7]);

  // stage 7
  output[0] = x6[1];
  output[1] = x6[6];
  output[2] = x6[3];
  output[3] = x6[4];
  output[4] = x6[5];
  output[5] = x6[2];
  output[6] = x6[7];
  output[7] = x6[0];
}

static void fadst8x16_new_sse2(const __m128i *input, __m128i *output,
                               int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m128i __zero = _mm_setzero_si128();
  const __m128i __rounding = _mm_set1_epi32(1 << (cos_bit - 1));

  __m128i cospi_p32_p32 = pair_set_epi16(cospi[32], cospi[32]);
  __m128i cospi_p32_m32 = pair_set_epi16(cospi[32], -cospi[32]);
  __m128i cospi_p16_p48 = pair_set_epi16(cospi[16], cospi[48]);
  __m128i cospi_p48_m16 = pair_set_epi16(cospi[48], -cospi[16]);
  __m128i cospi_m48_p16 = pair_set_epi16(-cospi[48], cospi[16]);
  __m128i cospi_p08_p56 = pair_set_epi16(cospi[8], cospi[56]);
  __m128i cospi_p56_m08 = pair_set_epi16(cospi[56], -cospi[8]);
  __m128i cospi_p40_p24 = pair_set_epi16(cospi[40], cospi[24]);
  __m128i cospi_p24_m40 = pair_set_epi16(cospi[24], -cospi[40]);
  __m128i cospi_m56_p08 = pair_set_epi16(-cospi[56], cospi[8]);
  __m128i cospi_m24_p40 = pair_set_epi16(-cospi[24], cospi[40]);
  __m128i cospi_p02_p62 = pair_set_epi16(cospi[2], cospi[62]);
  __m128i cospi_p62_m02 = pair_set_epi16(cospi[62], -cospi[2]);
  __m128i cospi_p10_p54 = pair_set_epi16(cospi[10], cospi[54]);
  __m128i cospi_p54_m10 = pair_set_epi16(cospi[54], -cospi[10]);
  __m128i cospi_p18_p46 = pair_set_epi16(cospi[18], cospi[46]);
  __m128i cospi_p46_m18 = pair_set_epi16(cospi[46], -cospi[18]);
  __m128i cospi_p26_p38 = pair_set_epi16(cospi[26], cospi[38]);
  __m128i cospi_p38_m26 = pair_set_epi16(cospi[38], -cospi[26]);
  __m128i cospi_p34_p30 = pair_set_epi16(cospi[34], cospi[30]);
  __m128i cospi_p30_m34 = pair_set_epi16(cospi[30], -cospi[34]);
  __m128i cospi_p42_p22 = pair_set_epi16(cospi[42], cospi[22]);
  __m128i cospi_p22_m42 = pair_set_epi16(cospi[22], -cospi[42]);
  __m128i cospi_p50_p14 = pair_set_epi16(cospi[50], cospi[14]);
  __m128i cospi_p14_m50 = pair_set_epi16(cospi[14], -cospi[50]);
  __m128i cospi_p58_p06 = pair_set_epi16(cospi[58], cospi[6]);
  __m128i cospi_p06_m58 = pair_set_epi16(cospi[6], -cospi[58]);

  // stage 1
  __m128i x1[16];
  x1[0] = input[0];
  x1[1] = _mm_subs_epi16(__zero, input[15]);
  x1[2] = _mm_subs_epi16(__zero, input[7]);
  x1[3] = input[8];
  x1[4] = _mm_subs_epi16(__zero, input[3]);
  x1[5] = input[12];
  x1[6] = input[4];
  x1[7] = _mm_subs_epi16(__zero, input[11]);
  x1[8] = _mm_subs_epi16(__zero, input[1]);
  x1[9] = input[14];
  x1[10] = input[6];
  x1[11] = _mm_subs_epi16(__zero, input[9]);
  x1[12] = input[2];
  x1[13] = _mm_subs_epi16(__zero, input[13]);
  x1[14] = _mm_subs_epi16(__zero, input[5]);
  x1[15] = input[10];

  // stage 2
  __m128i x2[16];
  x2[0] = x1[0];
  x2[1] = x1[1];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[2], x1[3], x2[2], x2[3]);
  x2[4] = x1[4];
  x2[5] = x1[5];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[6], x1[7], x2[6], x2[7]);
  x2[8] = x1[8];
  x2[9] = x1[9];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[10], x1[11], x2[10], x2[11]);
  x2[12] = x1[12];
  x2[13] = x1[13];
  btf_16_sse2(cospi_p32_p32, cospi_p32_m32, x1[14], x1[15], x2[14], x2[15]);

  // stage 3
  __m128i x3[16];
  x3[0] = _mm_adds_epi16(x2[0], x2[2]);
  x3[2] = _mm_subs_epi16(x2[0], x2[2]);
  x3[1] = _mm_adds_epi16(x2[1], x2[3]);
  x3[3] = _mm_subs_epi16(x2[1], x2[3]);
  x3[4] = _mm_adds_epi16(x2[4], x2[6]);
  x3[6] = _mm_subs_epi16(x2[4], x2[6]);
  x3[5] = _mm_adds_epi16(x2[5], x2[7]);
  x3[7] = _mm_subs_epi16(x2[5], x2[7]);
  x3[8] = _mm_adds_epi16(x2[8], x2[10]);
  x3[10] = _mm_subs_epi16(x2[8], x2[10]);
  x3[9] = _mm_adds_epi16(x2[9], x2[11]);
  x3[11] = _mm_subs_epi16(x2[9], x2[11]);
  x3[12] = _mm_adds_epi16(x2[12], x2[14]);
  x3[14] = _mm_subs_epi16(x2[12], x2[14]);
  x3[13] = _mm_adds_epi16(x2[13], x2[15]);
  x3[15] = _mm_subs_epi16(x2[13], x2[15]);

  // stage 4
  __m128i x4[16];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_sse2(cospi_p16_p48, cospi_p48_m16, x3[4], x3[5], x4[4], x4[5]);
  btf_16_sse2(cospi_m48_p16, cospi_p16_p48, x3[6], x3[7], x4[6], x4[7]);
  x4[8] = x3[8];
  x4[9] = x3[9];
  x4[10] = x3[10];
  x4[11] = x3[11];
  btf_16_sse2(cospi_p16_p48, cospi_p48_m16, x3[12], x3[13], x4[12], x4[13]);
  btf_16_sse2(cospi_m48_p16, cospi_p16_p48, x3[14], x3[15], x4[14], x4[15]);

  // stage 5
  __m128i x5[16];
  x5[0] = _mm_adds_epi16(x4[0], x4[4]);
  x5[4] = _mm_subs_epi16(x4[0], x4[4]);
  x5[1] = _mm_adds_epi16(x4[1], x4[5]);
  x5[5] = _mm_subs_epi16(x4[1], x4[5]);
  x5[2] = _mm_adds_epi16(x4[2], x4[6]);
  x5[6] = _mm_subs_epi16(x4[2], x4[6]);
  x5[3] = _mm_adds_epi16(x4[3], x4[7]);
  x5[7] = _mm_subs_epi16(x4[3], x4[7]);
  x5[8] = _mm_adds_epi16(x4[8], x4[12]);
  x5[12] = _mm_subs_epi16(x4[8], x4[12]);
  x5[9] = _mm_adds_epi16(x4[9], x4[13]);
  x5[13] = _mm_subs_epi16(x4[9], x4[13]);
  x5[10] = _mm_adds_epi16(x4[10], x4[14]);
  x5[14] = _mm_subs_epi16(x4[10], x4[14]);
  x5[11] = _mm_adds_epi16(x4[11], x4[15]);
  x5[15] = _mm_subs_epi16(x4[11], x4[15]);

  // stage 6
  __m128i x6[16];
  x6[0] = x5[0];
  x6[1] = x5[1];
  x6[2] = x5[2];
  x6[3] = x5[3];
  x6[4] = x5[4];
  x6[5] = x5[5];
  x6[6] = x5[6];
  x6[7] = x5[7];
  btf_16_sse2(cospi_p08_p56, cospi_p56_m08, x5[8], x5[9], x6[8], x6[9]);
  btf_16_sse2(cospi_p40_p24, cospi_p24_m40, x5[10], x5[11], x6[10], x6[11]);
  btf_16_sse2(cospi_m56_p08, cospi_p08_p56, x5[12], x5[13], x6[12], x6[13]);
  btf_16_sse2(cospi_m24_p40, cospi_p40_p24, x5[14], x5[15], x6[14], x6[15]);

  // stage 7
  __m128i x7[16];
  x7[0] = _mm_adds_epi16(x6[0], x6[8]);
  x7[8] = _mm_subs_epi16(x6[0], x6[8]);
  x7[1] = _mm_adds_epi16(x6[1], x6[9]);
  x7[9] = _mm_subs_epi16(x6[1], x6[9]);
  x7[2] = _mm_adds_epi16(x6[2], x6[10]);
  x7[10] = _mm_subs_epi16(x6[2], x6[10]);
  x7[3] = _mm_adds_epi16(x6[3], x6[11]);
  x7[11] = _mm_subs_epi16(x6[3], x6[11]);
  x7[4] = _mm_adds_epi16(x6[4], x6[12]);
  x7[12] = _mm_subs_epi16(x6[4], x6[12]);
  x7[5] = _mm_adds_epi16(x6[5], x6[13]);
  x7[13] = _mm_subs_epi16(x6[5], x6[13]);
  x7[6] = _mm_adds_epi16(x6[6], x6[14]);
  x7[14] = _mm_subs_epi16(x6[6], x6[14]);
  x7[7] = _mm_adds_epi16(x6[7], x6[15]);
  x7[15] = _mm_subs_epi16(x6[7], x6[15]);

  // stage 8
  __m128i x8[16];
  btf_16_sse2(cospi_p02_p62, cospi_p62_m02, x7[0], x7[1], x8[0], x8[1]);
  btf_16_sse2(cospi_p10_p54, cospi_p54_m10, x7[2], x7[3], x8[2], x8[3]);
  btf_16_sse2(cospi_p18_p46, cospi_p46_m18, x7[4], x7[5], x8[4], x8[5]);
  btf_16_sse2(cospi_p26_p38, cospi_p38_m26, x7[6], x7[7], x8[6], x8[7]);
  btf_16_sse2(cospi_p34_p30, cospi_p30_m34, x7[8], x7[9], x8[8], x8[9]);
  btf_16_sse2(cospi_p42_p22, cospi_p22_m42, x7[10], x7[11], x8[10], x8[11]);
  btf_16_sse2(cospi_p50_p14, cospi_p14_m50, x7[12], x7[13], x8[12], x8[13]);
  btf_16_sse2(cospi_p58_p06, cospi_p06_m58, x7[14], x7[15], x8[14], x8[15]);

  // stage 9
  output[0] = x8[1];
  output[1] = x8[14];
  output[2] = x8[3];
  output[3] = x8[12];
  output[4] = x8[5];
  output[5] = x8[10];
  output[6] = x8[7];
  output[7] = x8[8];
  output[8] = x8[9];
  output[9] = x8[6];
  output[10] = x8[11];
  output[11] = x8[4];
  output[12] = x8[13];
  output[13] = x8[2];
  output[14] = x8[15];
  output[15] = x8[0];
}

static const transform_1d_sse2 col_txfm4x4_arr[TX_TYPES] = {
  fdct4x4_new_sse2,       // DCT_DCT
  fadst4x4_new_sse2,      // ADST_DCT
  fdct4x4_new_sse2,       // DCT_ADST
  fadst4x4_new_sse2,      // ADST_ADST
  fadst4x4_new_sse2,      // FLIPADST_DCT
  fdct4x4_new_sse2,       // DCT_FLIPADST
  fadst4x4_new_sse2,      // FLIPADST_FLIPADST
  fadst4x4_new_sse2,      // ADST_FLIPADST
  fadst4x4_new_sse2,      // FLIPADST_ADST
  fidentity4x4_new_sse2,  // IDTX
  fdct4x4_new_sse2,       // V_DCT
  fidentity4x4_new_sse2,  // H_DCT
  fadst4x4_new_sse2,      // V_ADST
  fidentity4x4_new_sse2,  // H_ADST
  fadst4x4_new_sse2,      // V_FLIPADST
  fidentity4x4_new_sse2   // H_FLIPADST
};

static const transform_1d_sse2 row_txfm4x4_arr[TX_TYPES] = {
  fdct4x4_new_sse2,       // DCT_DCT
  fdct4x4_new_sse2,       // ADST_DCT
  fadst4x4_new_sse2,      // DCT_ADST
  fadst4x4_new_sse2,      // ADST_ADST
  fdct4x4_new_sse2,       // FLIPADST_DCT
  fadst4x4_new_sse2,      // DCT_FLIPADST
  fadst4x4_new_sse2,      // FLIPADST_FLIPADST
  fadst4x4_new_sse2,      // ADST_FLIPADST
  fadst4x4_new_sse2,      // FLIPADST_ADST
  fidentity4x4_new_sse2,  // IDTX
  fidentity4x4_new_sse2,  // V_DCT
  fdct4x4_new_sse2,       // H_DCT
  fidentity4x4_new_sse2,  // V_ADST
  fadst4x4_new_sse2,      // H_ADST
  fidentity4x4_new_sse2,  // V_FLIPADST
  fadst4x4_new_sse2       // H_FLIPADST
};

static const transform_1d_sse2 col_txfm4x8_arr[TX_TYPES] = {
  fdct4x8_new_sse2,       // DCT_DCT
  fadst4x8_new_sse2,      // ADST_DCT
  fdct4x8_new_sse2,       // DCT_ADST
  fadst4x8_new_sse2,      // ADST_ADST
  fadst4x8_new_sse2,      // FLIPADST_DCT
  fdct4x8_new_sse2,       // DCT_FLIPADST
  fadst4x8_new_sse2,      // FLIPADST_FLIPADST
  fadst4x8_new_sse2,      // ADST_FLIPADST
  fadst4x8_new_sse2,      // FLIPADST_ADST
  fidentity8x8_new_sse2,  // IDTX
  fdct4x8_new_sse2,       // V_DCT
  fidentity8x8_new_sse2,  // H_DCT
  fadst4x8_new_sse2,      // V_ADST
  fidentity8x8_new_sse2,  // H_ADST
  fadst4x8_new_sse2,      // V_FLIPADST
  fidentity8x8_new_sse2   // H_FLIPADST
};

static const transform_1d_sse2 row_txfm8x4_arr[TX_TYPES] = {
  fdct8x4_new_sse2,       // DCT_DCT
  fdct8x4_new_sse2,       // ADST_DCT
  fadst8x4_new_sse2,      // DCT_ADST
  fadst8x4_new_sse2,      // ADST_ADST
  fdct8x4_new_sse2,       // FLIPADST_DCT
  fadst8x4_new_sse2,      // DCT_FLIPADST
  fadst8x4_new_sse2,      // FLIPADST_FLIPADST
  fadst8x4_new_sse2,      // ADST_FLIPADST
  fadst8x4_new_sse2,      // FLIPADST_ADST
  fidentity8x4_new_sse2,  // IDTX
  fidentity8x4_new_sse2,  // V_DCT
  fdct8x4_new_sse2,       // H_DCT
  fidentity8x4_new_sse2,  // V_ADST
  fadst8x4_new_sse2,      // H_ADST
  fidentity8x4_new_sse2,  // V_FLIPADST
  fadst8x4_new_sse2       // H_FLIPADST
};

static const transform_1d_sse2 col_txfm8x4_arr[TX_TYPES] = {
  fdct8x4_new_sse2,       // DCT_DCT
  fadst8x4_new_sse2,      // ADST_DCT
  fdct8x4_new_sse2,       // DCT_ADST
  fadst8x4_new_sse2,      // ADST_ADST
  fadst8x4_new_sse2,      // FLIPADST_DCT
  fdct8x4_new_sse2,       // DCT_FLIPADST
  fadst8x4_new_sse2,      // FLIPADST_FLIPADST
  fadst8x4_new_sse2,      // ADST_FLIPADST
  fadst8x4_new_sse2,      // FLIPADST_ADST
  fidentity8x4_new_sse2,  // IDTX
  fdct8x4_new_sse2,       // V_DCT
  fidentity8x4_new_sse2,  // H_DCT
  fadst8x4_new_sse2,      // V_ADST
  fidentity8x4_new_sse2,  // H_ADST
  fadst8x4_new_sse2,      // V_FLIPADST
  fidentity8x4_new_sse2   // H_FLIPADST
};

static const transform_1d_sse2 row_txfm4x8_arr[TX_TYPES] = {
  fdct4x8_new_sse2,       // DCT_DCT
  fdct4x8_new_sse2,       // ADST_DCT
  fadst4x8_new_sse2,      // DCT_ADST
  fadst4x8_new_sse2,      // ADST_ADST
  fdct4x8_new_sse2,       // FLIPADST_DCT
  fadst4x8_new_sse2,      // DCT_FLIPADST
  fadst4x8_new_sse2,      // FLIPADST_FLIPADST
  fadst4x8_new_sse2,      // ADST_FLIPADST
  fadst4x8_new_sse2,      // FLIPADST_ADST
  fidentity8x8_new_sse2,  // IDTX
  fidentity8x8_new_sse2,  // V_DCT
  fdct4x8_new_sse2,       // H_DCT
  fidentity8x8_new_sse2,  // V_ADST
  fadst4x8_new_sse2,      // H_ADST
  fidentity8x8_new_sse2,  // V_FLIPADST
  fadst4x8_new_sse2       // H_FLIPADST
};

static const transform_1d_sse2 col_txfm8x8_arr[TX_TYPES] = {
  fdct8x8_new_sse2,       // DCT_DCT
  fadst8x8_new_sse2,      // ADST_DCT
  fdct8x8_new_sse2,       // DCT_ADST
  fadst8x8_new_sse2,      // ADST_ADST
  fadst8x8_new_sse2,      // FLIPADST_DCT
  fdct8x8_new_sse2,       // DCT_FLIPADST
  fadst8x8_new_sse2,      // FLIPADST_FLIPADST
  fadst8x8_new_sse2,      // ADST_FLIPADST
  fadst8x8_new_sse2,      // FLIPADST_ADST
  fidentity8x8_new_sse2,  // IDTX
  fdct8x8_new_sse2,       // V_DCT
  fidentity8x8_new_sse2,  // H_DCT
  fadst8x8_new_sse2,      // V_ADST
  fidentity8x8_new_sse2,  // H_ADST
  fadst8x8_new_sse2,      // V_FLIPADST
  fidentity8x8_new_sse2,  // H_FLIPADST
};

static const transform_1d_sse2 row_txfm8x8_arr[TX_TYPES] = {
  fdct8x8_new_sse2,       // DCT_DCT
  fdct8x8_new_sse2,       // ADST_DCT
  fadst8x8_new_sse2,      // DCT_ADST
  fadst8x8_new_sse2,      // ADST_ADST
  fdct8x8_new_sse2,       // FLIPADST_DCT
  fadst8x8_new_sse2,      // DCT_FLIPADST
  fadst8x8_new_sse2,      // FLIPADST_FLIPADST
  fadst8x8_new_sse2,      // ADST_FLIPADST
  fadst8x8_new_sse2,      // FLIPADST_ADST
  fidentity8x8_new_sse2,  // IDTX
  fidentity8x8_new_sse2,  // V_DCT
  fdct8x8_new_sse2,       // H_DCT
  fidentity8x8_new_sse2,  // V_ADST
  fadst8x8_new_sse2,      // H_ADST
  fidentity8x8_new_sse2,  // V_FLIPADST
  fadst8x8_new_sse2       // H_FLIPADST
};

static const transform_1d_sse2 col_txfm8x16_arr[TX_TYPES] = {
  fdct8x16_new_sse2,       // DCT_DCT
  fadst8x16_new_sse2,      // ADST_DCT
  fdct8x16_new_sse2,       // DCT_ADST
  fadst8x16_new_sse2,      // ADST_ADST
  fadst8x16_new_sse2,      // FLIPADST_DCT
  fdct8x16_new_sse2,       // DCT_FLIPADST
  fadst8x16_new_sse2,      // FLIPADST_FLIPADST
  fadst8x16_new_sse2,      // ADST_FLIPADST
  fadst8x16_new_sse2,      // FLIPADST_ADST
  fidentity8x16_new_sse2,  // IDTX
  fdct8x16_new_sse2,       // V_DCT
  fidentity8x16_new_sse2,  // H_DCT
  fadst8x16_new_sse2,      // V_ADST
  fidentity8x16_new_sse2,  // H_ADST
  fadst8x16_new_sse2,      // V_FLIPADST
  fidentity8x16_new_sse2   // H_FLIPADST
};

static const transform_1d_sse2 row_txfm8x16_arr[TX_TYPES] = {
  fdct8x16_new_sse2,       // DCT_DCT
  fdct8x16_new_sse2,       // ADST_DCT
  fadst8x16_new_sse2,      // DCT_ADST
  fadst8x16_new_sse2,      // ADST_ADST
  fdct8x16_new_sse2,       // FLIPADST_DCT
  fadst8x16_new_sse2,      // DCT_FLIPADST
  fadst8x16_new_sse2,      // FLIPADST_FLIPADST
  fadst8x16_new_sse2,      // ADST_FLIPADST
  fadst8x16_new_sse2,      // FLIPADST_ADST
  fidentity8x16_new_sse2,  // IDTX
  fidentity8x16_new_sse2,  // V_DCT
  fdct8x16_new_sse2,       // H_DCT
  fidentity8x16_new_sse2,  // V_ADST
  fadst8x16_new_sse2,      // H_ADST
  fidentity8x16_new_sse2,  // V_FLIPADST
  fadst8x16_new_sse2       // H_FLIPADST
};

static const transform_1d_sse2 row_txfm8x32_arr[TX_TYPES] = {
  fdct8x32_new_sse2,       // DCT_DCT
  NULL,                    // ADST_DCT
  NULL,                    // DCT_ADST
  NULL,                    // ADST_ADST
  NULL,                    // FLIPADST_DCT
  NULL,                    // DCT_FLIPADST
  NULL,                    // FLIPADST_FLIPADST
  NULL,                    // ADST_FLIPADST
  NULL,                    // FLIPADST_ADST
  fidentity8x32_new_sse2,  // IDTX
  fidentity8x32_new_sse2,  // V_DCT
  fdct8x32_new_sse2,       // H_DCT
  NULL,                    // V_ADST
  NULL,                    // H_ADST
  NULL,                    // V_FLIPADST
  NULL                     // H_FLIPADST
};

void av1_lowbd_fwd_txfm2d_4x4_sse2(const int16_t *input, int32_t *output,
                                   int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[4], buf1[4], *buf;
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_4X4];
  const int txw_idx = get_txw_idx(TX_4X4);
  const int txh_idx = get_txh_idx(TX_4X4);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 4;
  const int height = 4;
  const transform_1d_sse2 col_txfm = col_txfm4x4_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm4x4_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip) {
    load_buffer_16bit_to_16bit_w4_flip(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit_w4(input, stride, buf0, height);
  }
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_4x4(buf0, buf1);

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_4x4(buf, buf);
  store_buffer_16bit_to_32bit_w4(buf, output, width, height);
}

void av1_lowbd_fwd_txfm2d_4x8_sse2(const int16_t *input, int32_t *output,
                                   int stride, TX_TYPE tx_type, int bd) {
  (void)stride;
  (void)bd;
  __m128i buf0[8], buf1[8], *buf;
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_4X8];
  const int txw_idx = get_txw_idx(TX_4X8);
  const int txh_idx = get_txh_idx(TX_4X8);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 4;
  const int height = 8;
  const transform_1d_sse2 col_txfm = col_txfm4x8_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x4_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip) {
    load_buffer_16bit_to_16bit_w4_flip(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit_w4(input, stride, buf0, height);
  }
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_4x8(buf0, buf1);

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_8x4(buf, buf);
  store_rect_buffer_16bit_to_32bit_w4(buf, output, width, height);
}

void av1_lowbd_fwd_txfm2d_4x16_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[16], buf1[16];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_4X16];
  const int txw_idx = get_txw_idx(TX_4X16);
  const int txh_idx = get_txh_idx(TX_4X16);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 4;
  const int height = 16;
  const transform_1d_sse2 col_txfm = col_txfm8x16_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x4_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip) {
    load_buffer_16bit_to_16bit_w4_flip(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit_w4(input, stride, buf0, height);
  }
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_4x8(buf0, buf1);
  transpose_16bit_4x8(buf0 + 8, buf1 + 8);

  for (int i = 0; i < 2; i++) {
    __m128i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_sse2(buf1 + 8 * i, buf, width);
    } else {
      buf = buf1 + 8 * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    transpose_16bit_8x4(buf, buf);
    store_buffer_16bit_to_32bit_w4(buf, output + 8 * width * i, width, 8);
  }
}

void av1_lowbd_fwd_txfm2d_8x4_sse2(const int16_t *input, int32_t *output,
                                   int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[8], buf1[8], *buf;
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_8X4];
  const int txw_idx = get_txw_idx(TX_8X4);
  const int txh_idx = get_txh_idx(TX_8X4);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 8;
  const int height = 4;
  const transform_1d_sse2 col_txfm = col_txfm8x4_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm4x8_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip)
    load_buffer_16bit_to_16bit_flip(input, stride, buf0, height);
  else
    load_buffer_16bit_to_16bit(input, stride, buf0, height);
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_8x8(buf0, buf1);

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_8x8(buf, buf);
  store_rect_buffer_16bit_to_32bit_w8(buf, output, width, height);
}

void av1_lowbd_fwd_txfm2d_8x8_sse2(const int16_t *input, int32_t *output,
                                   int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[8], buf1[8], *buf;
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_8X8];
  const int txw_idx = get_txw_idx(TX_8X8);
  const int txh_idx = get_txh_idx(TX_8X8);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 8;
  const int height = 8;
  const transform_1d_sse2 col_txfm = col_txfm8x8_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x8_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip)
    load_buffer_16bit_to_16bit_flip(input, stride, buf0, height);
  else
    load_buffer_16bit_to_16bit(input, stride, buf0, height);
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_8x8(buf0, buf1);

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_8x8(buf, buf);
  store_buffer_16bit_to_32bit_w8(buf, output, width, height);
}

void av1_lowbd_fwd_txfm2d_8x16_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[16], buf1[16];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_8X16];
  const int txw_idx = get_txw_idx(TX_8X16);
  const int txh_idx = get_txh_idx(TX_8X16);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 8;
  const int height = 16;
  const transform_1d_sse2 col_txfm = col_txfm8x16_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x8_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip) {
    load_buffer_16bit_to_16bit_flip(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit(input, stride, buf0, height);
  }
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_8x8(buf0, buf1);
  transpose_16bit_8x8(buf0 + 8, buf1 + 8);

  for (int i = 0; i < 2; i++) {
    __m128i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_sse2(buf1 + width * i, buf, width);
    } else {
      buf = buf1 + width * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    transpose_16bit_8x8(buf, buf);
    store_rect_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width, 8);
  }
}

void av1_lowbd_fwd_txfm2d_8x32_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[32], buf1[32];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_8X32];
  const int txw_idx = get_txw_idx(TX_8X32);
  const int txh_idx = get_txh_idx(TX_8X32);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 8;
  const int height = 32;
  const transform_1d_sse2 col_txfm = col_txfm8x32_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x8_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  if (ud_flip) {
    load_buffer_16bit_to_16bit_flip(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit(input, stride, buf0, height);
  }
  round_shift_16bit(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit(buf0, height, shift[1]);
  transpose_16bit_8x8(buf0, buf1);
  transpose_16bit_8x8(buf0 + 8, buf1 + 8);
  transpose_16bit_8x8(buf0 + 16, buf1 + 16);
  transpose_16bit_8x8(buf0 + 24, buf1 + 24);

  for (int i = 0; i < 4; i++) {
    __m128i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_sse2(buf1 + width * i, buf, width);
    } else {
      buf = buf1 + width * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    transpose_16bit_8x8(buf, buf);
    store_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width, 8);
  }
}

void av1_lowbd_fwd_txfm2d_16x4_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[16], buf1[16];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_16X4];
  const int txw_idx = get_txw_idx(TX_16X4);
  const int txh_idx = get_txh_idx(TX_16X4);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 16;
  const int height = 4;
  const transform_1d_sse2 col_txfm = col_txfm8x4_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x16_arr[tx_type];
  __m128i *buf;
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  for (int i = 0; i < 2; i++) {
    if (ud_flip) {
      load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
    } else {
      load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
    }
    round_shift_16bit(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit(buf0, height, shift[1]);
    transpose_16bit_8x4(buf0, buf1 + 8 * i);
  }

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_4x8(buf, buf);
  store_buffer_16bit_to_32bit_w8(buf, output, width, height);
  transpose_16bit_4x8(buf + 8, buf + 8);
  store_buffer_16bit_to_32bit_w8(buf + 8, output + 8, width, height);
}

void av1_lowbd_fwd_txfm2d_16x8_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[16], buf1[16];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_16X8];
  const int txw_idx = get_txw_idx(TX_16X8);
  const int txh_idx = get_txh_idx(TX_16X8);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 16;
  const int height = 8;
  const transform_1d_sse2 col_txfm = col_txfm8x8_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x16_arr[tx_type];
  __m128i *buf;
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  for (int i = 0; i < 2; i++) {
    if (ud_flip) {
      load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
    } else {
      load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
    }
    round_shift_16bit(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit(buf0, height, shift[1]);
    transpose_16bit_8x8(buf0, buf1 + 8 * i);
  }

  if (lr_flip) {
    buf = buf0;
    flip_buf_sse2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit(buf, width, shift[2]);
  transpose_16bit_8x8(buf, buf);
  store_rect_buffer_16bit_to_32bit_w8(buf, output, width, height);
  transpose_16bit_8x8(buf + 8, buf + 8);
  store_rect_buffer_16bit_to_32bit_w8(buf + 8, output + 8, width, height);
}

void av1_lowbd_fwd_txfm2d_16x16_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[16], buf1[32];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_16X16];
  const int txw_idx = get_txw_idx(TX_16X16);
  const int txh_idx = get_txh_idx(TX_16X16);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 16;
  const int height = 16;
  const transform_1d_sse2 col_txfm = col_txfm8x16_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x16_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  for (int i = 0; i < 2; i++) {
    if (ud_flip) {
      load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
    } else {
      load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
    }
    round_shift_16bit(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit(buf0, height, shift[1]);
    transpose_16bit_8x8(buf0, buf1 + 0 * width + 8 * i);
    transpose_16bit_8x8(buf0 + 8, buf1 + 1 * width + 8 * i);
  }

  for (int i = 0; i < 2; i++) {
    __m128i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_sse2(buf1 + width * i, buf, width);
    } else {
      buf = buf1 + width * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    transpose_16bit_8x8(buf, buf);
    store_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width, 8);
    transpose_16bit_8x8(buf + 8, buf + 8);
    store_buffer_16bit_to_32bit_w8(buf + 8, output + 8 * width * i + 8, width,
                                   8);
  }
}

void av1_lowbd_fwd_txfm2d_16x32_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[32], buf1[64];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_16X32];
  const int txw_idx = get_txw_idx(TX_16X32);
  const int txh_idx = get_txh_idx(TX_16X32);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 16;
  const int height = 32;
  const transform_1d_sse2 col_txfm = col_txfm8x32_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x16_arr[tx_type];

  if (col_txfm != NULL && row_txfm != NULL) {
    int ud_flip, lr_flip;
    get_flip_cfg(tx_type, &ud_flip, &lr_flip);

    for (int i = 0; i < 2; i++) {
      if (ud_flip) {
        load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
      } else {
        load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
      }
      round_shift_16bit(buf0, height, shift[0]);
      col_txfm(buf0, buf0, cos_bit_col);
      round_shift_16bit(buf0, height, shift[1]);
      transpose_16bit_8x8(buf0 + 0 * 8, buf1 + 0 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 1 * 8, buf1 + 1 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 2 * 8, buf1 + 2 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 3 * 8, buf1 + 3 * width + 8 * i);
    }

    for (int i = 0; i < 4; i++) {
      __m128i *buf;
      if (lr_flip) {
        buf = buf0;
        flip_buf_sse2(buf1 + width * i, buf, width);
      } else {
        buf = buf1 + width * i;
      }
      row_txfm(buf, buf, cos_bit_row);
      round_shift_16bit(buf, width, shift[2]);
      transpose_16bit_8x8(buf, buf);
      store_rect_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width,
                                          8);
      transpose_16bit_8x8(buf + 8, buf + 8);
      store_rect_buffer_16bit_to_32bit_w8(buf + 8, output + 8 * width * i + 8,
                                          width, 8);
    }
  } else {
    av1_fwd_txfm2d_16x32_c(input, output, stride, tx_type, bd);
  }
}

void av1_lowbd_fwd_txfm2d_32x8_sse2(const int16_t *input, int32_t *output,
                                    int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[32], buf1[32];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_32X8];
  const int txw_idx = get_txw_idx(TX_32X8);
  const int txh_idx = get_txh_idx(TX_32X8);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 32;
  const int height = 8;
  const transform_1d_sse2 col_txfm = col_txfm8x8_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x32_arr[tx_type];

  if (col_txfm != NULL && row_txfm != NULL) {
    int ud_flip, lr_flip;
    get_flip_cfg(tx_type, &ud_flip, &lr_flip);

    for (int i = 0; i < 4; i++) {
      if (ud_flip) {
        load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
      } else {
        load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
      }
      round_shift_16bit(buf0, height, shift[0]);
      col_txfm(buf0, buf0, cos_bit_col);
      round_shift_16bit(buf0, height, shift[1]);
      transpose_16bit_8x8(buf0, buf1 + 0 * width + 8 * i);
    }

    for (int i = 0; i < 1; i++) {
      __m128i *buf;
      if (lr_flip) {
        buf = buf0;
        flip_buf_sse2(buf1 + width * i, buf, width);
      } else {
        buf = buf1 + width * i;
      }
      row_txfm(buf, buf, cos_bit_row);
      round_shift_16bit(buf, width, shift[2]);
      transpose_16bit_8x8(buf, buf);
      store_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width,
                                     height);
      transpose_16bit_8x8(buf + 8, buf + 8);
      store_buffer_16bit_to_32bit_w8(buf + 8, output + 8 * width * i + 8, width,
                                     height);
      transpose_16bit_8x8(buf + 16, buf + 16);
      store_buffer_16bit_to_32bit_w8(buf + 16, output + 8 * width * i + 16,
                                     width, height);
      transpose_16bit_8x8(buf + 24, buf + 24);
      store_buffer_16bit_to_32bit_w8(buf + 24, output + 8 * width * i + 24,
                                     width, height);
    }
  } else {
    av1_fwd_txfm2d_32x16_c(input, output, stride, tx_type, bd);
  }
}

void av1_lowbd_fwd_txfm2d_32x16_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[32], buf1[64];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_32X16];
  const int txw_idx = get_txw_idx(TX_32X16);
  const int txh_idx = get_txh_idx(TX_32X16);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 32;
  const int height = 16;
  const transform_1d_sse2 col_txfm = col_txfm8x16_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x32_arr[tx_type];

  if (col_txfm != NULL && row_txfm != NULL) {
    int ud_flip, lr_flip;
    get_flip_cfg(tx_type, &ud_flip, &lr_flip);

    for (int i = 0; i < 4; i++) {
      if (ud_flip) {
        load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
      } else {
        load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
      }
      round_shift_16bit(buf0, height, shift[0]);
      col_txfm(buf0, buf0, cos_bit_col);
      round_shift_16bit(buf0, height, shift[1]);
      transpose_16bit_8x8(buf0, buf1 + 0 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 8, buf1 + 1 * width + 8 * i);
    }

    for (int i = 0; i < 2; i++) {
      __m128i *buf;
      if (lr_flip) {
        buf = buf0;
        flip_buf_sse2(buf1 + width * i, buf, width);
      } else {
        buf = buf1 + width * i;
      }
      row_txfm(buf, buf, cos_bit_row);
      round_shift_16bit(buf, width, shift[2]);
      transpose_16bit_8x8(buf, buf);
      store_rect_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width,
                                          8);
      transpose_16bit_8x8(buf + 8, buf + 8);
      store_rect_buffer_16bit_to_32bit_w8(buf + 8, output + 8 * width * i + 8,
                                          width, 8);
      transpose_16bit_8x8(buf + 16, buf + 16);
      store_rect_buffer_16bit_to_32bit_w8(buf + 16, output + 8 * width * i + 16,
                                          width, 8);
      transpose_16bit_8x8(buf + 24, buf + 24);
      store_rect_buffer_16bit_to_32bit_w8(buf + 24, output + 8 * width * i + 24,
                                          width, 8);
    }
  } else {
    av1_fwd_txfm2d_32x16_c(input, output, stride, tx_type, bd);
  }
}

void av1_lowbd_fwd_txfm2d_32x32_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m128i buf0[32], buf1[128];
  const int8_t *shift = av1_fwd_txfm_shift_ls[TX_32X32];
  const int txw_idx = get_txw_idx(TX_32X32);
  const int txh_idx = get_txh_idx(TX_32X32);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 32;
  const int height = 32;
  const transform_1d_sse2 col_txfm = col_txfm8x32_arr[tx_type];
  const transform_1d_sse2 row_txfm = row_txfm8x32_arr[tx_type];

  if (col_txfm != NULL && row_txfm != NULL) {
    int ud_flip, lr_flip;
    get_flip_cfg(tx_type, &ud_flip, &lr_flip);

    for (int i = 0; i < 4; i++) {
      if (ud_flip) {
        load_buffer_16bit_to_16bit_flip(input + 8 * i, stride, buf0, height);
      } else {
        load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
      }
      round_shift_16bit(buf0, height, shift[0]);
      col_txfm(buf0, buf0, cos_bit_col);
      round_shift_16bit(buf0, height, shift[1]);
      transpose_16bit_8x8(buf0 + 0 * 8, buf1 + 0 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 1 * 8, buf1 + 1 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 2 * 8, buf1 + 2 * width + 8 * i);
      transpose_16bit_8x8(buf0 + 3 * 8, buf1 + 3 * width + 8 * i);
    }

    for (int i = 0; i < 4; i++) {
      __m128i *buf;
      if (lr_flip) {
        buf = buf0;
        flip_buf_sse2(buf1 + width * i, buf, width);
      } else {
        buf = buf1 + width * i;
      }
      row_txfm(buf, buf, cos_bit_row);
      round_shift_16bit(buf, width, shift[2]);
      transpose_16bit_8x8(buf, buf);
      store_buffer_16bit_to_32bit_w8(buf, output + 8 * width * i, width, 8);
      transpose_16bit_8x8(buf + 8, buf + 8);
      store_buffer_16bit_to_32bit_w8(buf + 8, output + 8 * width * i + 8, width,
                                     8);
      transpose_16bit_8x8(buf + 16, buf + 16);
      store_buffer_16bit_to_32bit_w8(buf + 16, output + 8 * width * i + 16,
                                     width, 8);
      transpose_16bit_8x8(buf + 24, buf + 24);
      store_buffer_16bit_to_32bit_w8(buf + 24, output + 8 * width * i + 24,
                                     width, 8);
    }
  } else {
    av1_fwd_txfm2d_32x32_c(input, output, stride, tx_type, bd);
  }
}

void av1_lowbd_fwd_txfm2d_64x16_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_64X16;
  __m128i buf0[64], buf1[128];
  const int8_t *shift = av1_fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_sse2 col_txfm = fdct8x16_new_sse2;
  const transform_1d_sse2 row_txfm = fdct8x64_new_sse2;
  const int width_div8 = (width >> 3);
  const int height_div8 = (height >> 3);

  for (int i = 0; i < width_div8; i++) {
    load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
    round_shift_16bit(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit(buf0, height, shift[1]);
    for (int j = 0; j < height_div8; ++j) {
      transpose_16bit_8x8(buf0 + j * 8, buf1 + j * width + 8 * i);
    }
  }

  for (int i = 0; i < height_div8; i++) {
    __m128i *buf = buf1 + width * i;
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    int32_t *output8 = output + 8 * 32 * i;
    for (int j = 0; j < 4; ++j) {
      __m128i *buf8 = buf + 8 * j;
      transpose_16bit_8x8(buf8, buf8);
      store_buffer_16bit_to_32bit_w8(buf8, output8 + 8 * j, 32, 8);
    }
  }
}

void av1_lowbd_fwd_txfm2d_16x64_sse2(const int16_t *input, int32_t *output,
                                     int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_16X64;
  __m128i buf0[64], buf1[128];
  const int8_t *shift = av1_fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = av1_fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = av1_fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_sse2 col_txfm = fdct8x64_new_sse2;
  const transform_1d_sse2 row_txfm = fdct8x16_new_sse2;
  const int width_div8 = (width >> 3);
  const int height_div8 = (height >> 3);

  for (int i = 0; i < width_div8; i++) {
    load_buffer_16bit_to_16bit(input + 8 * i, stride, buf0, height);
    round_shift_16bit(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit(buf0, height, shift[1]);
    for (int j = 0; j < height_div8; ++j) {
      transpose_16bit_8x8(buf0 + j * 8, buf1 + j * width + 8 * i);
    }
  }

  for (int i = 0; i < AOMMIN(4, height_div8); i++) {
    __m128i *buf = buf1 + width * i;
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit(buf, width, shift[2]);
    int32_t *output8 = output + 8 * width * i;
    for (int j = 0; j < width_div8; ++j) {
      __m128i *buf8 = buf + 8 * j;
      transpose_16bit_8x8(buf8, buf8);
      store_buffer_16bit_to_32bit_w8(buf8, output8 + 8 * j, width, 8);
    }
  }
  // Zero out the bottom 16x32 area.
  memset(output + 16 * 32, 0, 16 * 32 * sizeof(*output));
}

static FwdTxfm2dFunc fwd_txfm2d_func_ls[TX_SIZES_ALL] = {
  av1_lowbd_fwd_txfm2d_4x4_sse2,    // 4x4 transform
  av1_lowbd_fwd_txfm2d_8x8_sse2,    // 8x8 transform
  av1_lowbd_fwd_txfm2d_16x16_sse2,  // 16x16 transform
  av1_lowbd_fwd_txfm2d_32x32_sse2,  // 32x32 transform
  NULL,                             // 64x64 transform
  av1_lowbd_fwd_txfm2d_4x8_sse2,    // 4x8 transform
  av1_lowbd_fwd_txfm2d_8x4_sse2,    // 8x4 transform
  av1_lowbd_fwd_txfm2d_8x16_sse2,   // 8x16 transform
  av1_lowbd_fwd_txfm2d_16x8_sse2,   // 16x8 transform
  av1_lowbd_fwd_txfm2d_16x32_sse2,  // 16x32 transform
  av1_lowbd_fwd_txfm2d_32x16_sse2,  // 32x16 transform
  NULL,                             // 32x64 transform
  NULL,                             // 64x32 transform
  av1_lowbd_fwd_txfm2d_4x16_sse2,   // 4x16 transform
  av1_lowbd_fwd_txfm2d_16x4_sse2,   // 16x4 transform
  av1_lowbd_fwd_txfm2d_8x32_sse2,   // 8x32 transform
  av1_lowbd_fwd_txfm2d_32x8_sse2,   // 32x8 transform
  av1_lowbd_fwd_txfm2d_16x64_sse2,  // 16x64 transform
  av1_lowbd_fwd_txfm2d_64x16_sse2,  // 64x16 transform
};

void av1_lowbd_fwd_txfm_sse2(const int16_t *src_diff, tran_low_t *coeff,
                             int diff_stride, TxfmParam *txfm_param) {
  FwdTxfm2dFunc fwd_txfm2d_func = fwd_txfm2d_func_ls[txfm_param->tx_size];

  if ((fwd_txfm2d_func == NULL) ||
      (txfm_param->lossless && txfm_param->tx_size == TX_4X4))
    av1_lowbd_fwd_txfm_c(src_diff, coeff, diff_stride, txfm_param);
  else
    fwd_txfm2d_func(src_diff, coeff, diff_stride, txfm_param->tx_type,
                    txfm_param->bd);
}
