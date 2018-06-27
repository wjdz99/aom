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

#include "config/av1_rtcd.h"

#include "av1/common/enums.h"
#include "av1/common/av1_txfm.h"
#include "av1/encoder/x86/av1_fwd_txfm_avx2.h"
#include "av1/common/x86/av1_txfm_sse2.h"
#include "av1/encoder/av1_fwd_txfm1d_cfg.h"
#include "av1/encoder/x86/av1_txfm1d_sse4.h"
#include "av1/encoder/x86/av1_fwd_txfm_sse2.h"
#include "aom_dsp/x86/txfm_common_avx2.h"

static INLINE void fdct16x16_new_avx2(const __m256i *input, __m256i *output,
                                      int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m256i __rounding = _mm256_set1_epi32(1 << (cos_bit - 1));

  __m256i cospi_m32_p32 = pair_set_w16_epi16(-cospi[32], cospi[32]);
  __m256i cospi_p32_p32 = pair_set_w16_epi16(cospi[32], cospi[32]);
  __m256i cospi_p32_m32 = pair_set_w16_epi16(cospi[32], -cospi[32]);
  __m256i cospi_p48_p16 = pair_set_w16_epi16(cospi[48], cospi[16]);
  __m256i cospi_m16_p48 = pair_set_w16_epi16(-cospi[16], cospi[48]);
  __m256i cospi_m48_m16 = pair_set_w16_epi16(-cospi[48], -cospi[16]);
  __m256i cospi_p56_p08 = pair_set_w16_epi16(cospi[56], cospi[8]);
  __m256i cospi_m08_p56 = pair_set_w16_epi16(-cospi[8], cospi[56]);
  __m256i cospi_p24_p40 = pair_set_w16_epi16(cospi[24], cospi[40]);
  __m256i cospi_m40_p24 = pair_set_w16_epi16(-cospi[40], cospi[24]);
  __m256i cospi_p60_p04 = pair_set_w16_epi16(cospi[60], cospi[4]);
  __m256i cospi_m04_p60 = pair_set_w16_epi16(-cospi[4], cospi[60]);
  __m256i cospi_p28_p36 = pair_set_w16_epi16(cospi[28], cospi[36]);
  __m256i cospi_m36_p28 = pair_set_w16_epi16(-cospi[36], cospi[28]);
  __m256i cospi_p44_p20 = pair_set_w16_epi16(cospi[44], cospi[20]);
  __m256i cospi_m20_p44 = pair_set_w16_epi16(-cospi[20], cospi[44]);
  __m256i cospi_p12_p52 = pair_set_w16_epi16(cospi[12], cospi[52]);
  __m256i cospi_m52_p12 = pair_set_w16_epi16(-cospi[52], cospi[12]);

  // stage 1
  __m256i x1[16];
  x1[0] = _mm256_adds_epi16(input[0], input[15]);
  x1[15] = _mm256_subs_epi16(input[0], input[15]);
  x1[1] = _mm256_adds_epi16(input[1], input[14]);
  x1[14] = _mm256_subs_epi16(input[1], input[14]);
  x1[2] = _mm256_adds_epi16(input[2], input[13]);
  x1[13] = _mm256_subs_epi16(input[2], input[13]);
  x1[3] = _mm256_adds_epi16(input[3], input[12]);
  x1[12] = _mm256_subs_epi16(input[3], input[12]);
  x1[4] = _mm256_adds_epi16(input[4], input[11]);
  x1[11] = _mm256_subs_epi16(input[4], input[11]);
  x1[5] = _mm256_adds_epi16(input[5], input[10]);
  x1[10] = _mm256_subs_epi16(input[5], input[10]);
  x1[6] = _mm256_adds_epi16(input[6], input[9]);
  x1[9] = _mm256_subs_epi16(input[6], input[9]);
  x1[7] = _mm256_adds_epi16(input[7], input[8]);
  x1[8] = _mm256_subs_epi16(input[7], input[8]);

  // stage 2
  __m256i x2[16];
  x2[0] = _mm256_adds_epi16(x1[0], x1[7]);
  x2[7] = _mm256_subs_epi16(x1[0], x1[7]);
  x2[1] = _mm256_adds_epi16(x1[1], x1[6]);
  x2[6] = _mm256_subs_epi16(x1[1], x1[6]);
  x2[2] = _mm256_adds_epi16(x1[2], x1[5]);
  x2[5] = _mm256_subs_epi16(x1[2], x1[5]);
  x2[3] = _mm256_adds_epi16(x1[3], x1[4]);
  x2[4] = _mm256_subs_epi16(x1[3], x1[4]);
  x2[8] = x1[8];
  x2[9] = x1[9];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[10], x1[13], x2[10], x2[13]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[11], x1[12], x2[11], x2[12]);
  x2[14] = x1[14];
  x2[15] = x1[15];

  // stage 3
  __m256i x3[16];
  x3[0] = _mm256_adds_epi16(x2[0], x2[3]);
  x3[3] = _mm256_subs_epi16(x2[0], x2[3]);
  x3[1] = _mm256_adds_epi16(x2[1], x2[2]);
  x3[2] = _mm256_subs_epi16(x2[1], x2[2]);
  x3[4] = x2[4];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[5], x2[6], x3[5], x3[6]);
  x3[7] = x2[7];
  x3[8] = _mm256_adds_epi16(x2[8], x2[11]);
  x3[11] = _mm256_subs_epi16(x2[8], x2[11]);
  x3[9] = _mm256_adds_epi16(x2[9], x2[10]);
  x3[10] = _mm256_subs_epi16(x2[9], x2[10]);
  x3[12] = _mm256_subs_epi16(x2[15], x2[12]);
  x3[15] = _mm256_adds_epi16(x2[15], x2[12]);
  x3[13] = _mm256_subs_epi16(x2[14], x2[13]);
  x3[14] = _mm256_adds_epi16(x2[14], x2[13]);

  // stage 4
  __m256i x4[16];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x3[0], x3[1], x4[0], x4[1]);
  btf_16_w16_avx2(cospi_p48_p16, cospi_m16_p48, x3[2], x3[3], x4[2], x4[3]);
  x4[4] = _mm256_adds_epi16(x3[4], x3[5]);
  x4[5] = _mm256_subs_epi16(x3[4], x3[5]);
  x4[6] = _mm256_subs_epi16(x3[7], x3[6]);
  x4[7] = _mm256_adds_epi16(x3[7], x3[6]);
  x4[8] = x3[8];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[9], x3[14], x4[9], x4[14]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[10], x3[13], x4[10], x4[13]);
  x4[11] = x3[11];
  x4[12] = x3[12];
  x4[15] = x3[15];

  // stage 5
  __m256i x5[16];
  x5[0] = x4[0];
  x5[1] = x4[1];
  x5[2] = x4[2];
  x5[3] = x4[3];
  btf_16_w16_avx2(cospi_p56_p08, cospi_m08_p56, x4[4], x4[7], x5[4], x5[7]);
  btf_16_w16_avx2(cospi_p24_p40, cospi_m40_p24, x4[5], x4[6], x5[5], x5[6]);
  x5[8] = _mm256_adds_epi16(x4[8], x4[9]);
  x5[9] = _mm256_subs_epi16(x4[8], x4[9]);
  x5[10] = _mm256_subs_epi16(x4[11], x4[10]);
  x5[11] = _mm256_adds_epi16(x4[11], x4[10]);
  x5[12] = _mm256_adds_epi16(x4[12], x4[13]);
  x5[13] = _mm256_subs_epi16(x4[12], x4[13]);
  x5[14] = _mm256_subs_epi16(x4[15], x4[14]);
  x5[15] = _mm256_adds_epi16(x4[15], x4[14]);

  // stage 6
  __m256i x6[16];
  x6[0] = x5[0];
  x6[1] = x5[1];
  x6[2] = x5[2];
  x6[3] = x5[3];
  x6[4] = x5[4];
  x6[5] = x5[5];
  x6[6] = x5[6];
  x6[7] = x5[7];
  btf_16_w16_avx2(cospi_p60_p04, cospi_m04_p60, x5[8], x5[15], x6[8], x6[15]);
  btf_16_w16_avx2(cospi_p28_p36, cospi_m36_p28, x5[9], x5[14], x6[9], x6[14]);
  btf_16_w16_avx2(cospi_p44_p20, cospi_m20_p44, x5[10], x5[13], x6[10], x6[13]);
  btf_16_w16_avx2(cospi_p12_p52, cospi_m52_p12, x5[11], x5[12], x6[11], x6[12]);

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

static INLINE void fdct16x32_new_avx2(const __m256i *input, __m256i *output,
                                      int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m256i __rounding = _mm256_set1_epi32(1 << (cos_bit - 1));

  __m256i cospi_m32_p32 = pair_set_w16_epi16(-cospi[32], cospi[32]);
  __m256i cospi_p32_p32 = pair_set_w16_epi16(cospi[32], cospi[32]);
  __m256i cospi_m16_p48 = pair_set_w16_epi16(-cospi[16], cospi[48]);
  __m256i cospi_p48_p16 = pair_set_w16_epi16(cospi[48], cospi[16]);
  __m256i cospi_m48_m16 = pair_set_w16_epi16(-cospi[48], -cospi[16]);
  __m256i cospi_p32_m32 = pair_set_w16_epi16(cospi[32], -cospi[32]);
  __m256i cospi_p56_p08 = pair_set_w16_epi16(cospi[56], cospi[8]);
  __m256i cospi_m08_p56 = pair_set_w16_epi16(-cospi[8], cospi[56]);
  __m256i cospi_p24_p40 = pair_set_w16_epi16(cospi[24], cospi[40]);
  __m256i cospi_m40_p24 = pair_set_w16_epi16(-cospi[40], cospi[24]);
  __m256i cospi_m56_m08 = pair_set_w16_epi16(-cospi[56], -cospi[8]);
  __m256i cospi_m24_m40 = pair_set_w16_epi16(-cospi[24], -cospi[40]);
  __m256i cospi_p60_p04 = pair_set_w16_epi16(cospi[60], cospi[4]);
  __m256i cospi_m04_p60 = pair_set_w16_epi16(-cospi[4], cospi[60]);
  __m256i cospi_p28_p36 = pair_set_w16_epi16(cospi[28], cospi[36]);
  __m256i cospi_m36_p28 = pair_set_w16_epi16(-cospi[36], cospi[28]);
  __m256i cospi_p44_p20 = pair_set_w16_epi16(cospi[44], cospi[20]);
  __m256i cospi_m20_p44 = pair_set_w16_epi16(-cospi[20], cospi[44]);
  __m256i cospi_p12_p52 = pair_set_w16_epi16(cospi[12], cospi[52]);
  __m256i cospi_m52_p12 = pair_set_w16_epi16(-cospi[52], cospi[12]);
  __m256i cospi_p62_p02 = pair_set_w16_epi16(cospi[62], cospi[2]);
  __m256i cospi_m02_p62 = pair_set_w16_epi16(-cospi[2], cospi[62]);
  __m256i cospi_p30_p34 = pair_set_w16_epi16(cospi[30], cospi[34]);
  __m256i cospi_m34_p30 = pair_set_w16_epi16(-cospi[34], cospi[30]);
  __m256i cospi_p46_p18 = pair_set_w16_epi16(cospi[46], cospi[18]);
  __m256i cospi_m18_p46 = pair_set_w16_epi16(-cospi[18], cospi[46]);
  __m256i cospi_p14_p50 = pair_set_w16_epi16(cospi[14], cospi[50]);
  __m256i cospi_m50_p14 = pair_set_w16_epi16(-cospi[50], cospi[14]);
  __m256i cospi_p54_p10 = pair_set_w16_epi16(cospi[54], cospi[10]);
  __m256i cospi_m10_p54 = pair_set_w16_epi16(-cospi[10], cospi[54]);
  __m256i cospi_p22_p42 = pair_set_w16_epi16(cospi[22], cospi[42]);
  __m256i cospi_m42_p22 = pair_set_w16_epi16(-cospi[42], cospi[22]);
  __m256i cospi_p38_p26 = pair_set_w16_epi16(cospi[38], cospi[26]);
  __m256i cospi_m26_p38 = pair_set_w16_epi16(-cospi[26], cospi[38]);
  __m256i cospi_p06_p58 = pair_set_w16_epi16(cospi[6], cospi[58]);
  __m256i cospi_m58_p06 = pair_set_w16_epi16(-cospi[58], cospi[6]);

  // stage 1
  __m256i x1[32];
  x1[0] = _mm256_adds_epi16(input[0], input[31]);
  x1[31] = _mm256_subs_epi16(input[0], input[31]);
  x1[1] = _mm256_adds_epi16(input[1], input[30]);
  x1[30] = _mm256_subs_epi16(input[1], input[30]);
  x1[2] = _mm256_adds_epi16(input[2], input[29]);
  x1[29] = _mm256_subs_epi16(input[2], input[29]);
  x1[3] = _mm256_adds_epi16(input[3], input[28]);
  x1[28] = _mm256_subs_epi16(input[3], input[28]);
  x1[4] = _mm256_adds_epi16(input[4], input[27]);
  x1[27] = _mm256_subs_epi16(input[4], input[27]);
  x1[5] = _mm256_adds_epi16(input[5], input[26]);
  x1[26] = _mm256_subs_epi16(input[5], input[26]);
  x1[6] = _mm256_adds_epi16(input[6], input[25]);
  x1[25] = _mm256_subs_epi16(input[6], input[25]);
  x1[7] = _mm256_adds_epi16(input[7], input[24]);
  x1[24] = _mm256_subs_epi16(input[7], input[24]);
  x1[8] = _mm256_adds_epi16(input[8], input[23]);
  x1[23] = _mm256_subs_epi16(input[8], input[23]);
  x1[9] = _mm256_adds_epi16(input[9], input[22]);
  x1[22] = _mm256_subs_epi16(input[9], input[22]);
  x1[10] = _mm256_adds_epi16(input[10], input[21]);
  x1[21] = _mm256_subs_epi16(input[10], input[21]);
  x1[11] = _mm256_adds_epi16(input[11], input[20]);
  x1[20] = _mm256_subs_epi16(input[11], input[20]);
  x1[12] = _mm256_adds_epi16(input[12], input[19]);
  x1[19] = _mm256_subs_epi16(input[12], input[19]);
  x1[13] = _mm256_adds_epi16(input[13], input[18]);
  x1[18] = _mm256_subs_epi16(input[13], input[18]);
  x1[14] = _mm256_adds_epi16(input[14], input[17]);
  x1[17] = _mm256_subs_epi16(input[14], input[17]);
  x1[15] = _mm256_adds_epi16(input[15], input[16]);
  x1[16] = _mm256_subs_epi16(input[15], input[16]);

  // stage 2
  __m256i x2[32];
  x2[0] = _mm256_adds_epi16(x1[0], x1[15]);
  x2[15] = _mm256_subs_epi16(x1[0], x1[15]);
  x2[1] = _mm256_adds_epi16(x1[1], x1[14]);
  x2[14] = _mm256_subs_epi16(x1[1], x1[14]);
  x2[2] = _mm256_adds_epi16(x1[2], x1[13]);
  x2[13] = _mm256_subs_epi16(x1[2], x1[13]);
  x2[3] = _mm256_adds_epi16(x1[3], x1[12]);
  x2[12] = _mm256_subs_epi16(x1[3], x1[12]);
  x2[4] = _mm256_adds_epi16(x1[4], x1[11]);
  x2[11] = _mm256_subs_epi16(x1[4], x1[11]);
  x2[5] = _mm256_adds_epi16(x1[5], x1[10]);
  x2[10] = _mm256_subs_epi16(x1[5], x1[10]);
  x2[6] = _mm256_adds_epi16(x1[6], x1[9]);
  x2[9] = _mm256_subs_epi16(x1[6], x1[9]);
  x2[7] = _mm256_adds_epi16(x1[7], x1[8]);
  x2[8] = _mm256_subs_epi16(x1[7], x1[8]);
  x2[16] = x1[16];
  x2[17] = x1[17];
  x2[18] = x1[18];
  x2[19] = x1[19];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[20], x1[27], x2[20], x2[27]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[21], x1[26], x2[21], x2[26]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[22], x1[25], x2[22], x2[25]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[23], x1[24], x2[23], x2[24]);
  x2[28] = x1[28];
  x2[29] = x1[29];
  x2[30] = x1[30];
  x2[31] = x1[31];

  // stage 3
  __m256i x3[32];
  x3[0] = _mm256_adds_epi16(x2[0], x2[7]);
  x3[7] = _mm256_subs_epi16(x2[0], x2[7]);
  x3[1] = _mm256_adds_epi16(x2[1], x2[6]);
  x3[6] = _mm256_subs_epi16(x2[1], x2[6]);
  x3[2] = _mm256_adds_epi16(x2[2], x2[5]);
  x3[5] = _mm256_subs_epi16(x2[2], x2[5]);
  x3[3] = _mm256_adds_epi16(x2[3], x2[4]);
  x3[4] = _mm256_subs_epi16(x2[3], x2[4]);
  x3[8] = x2[8];
  x3[9] = x2[9];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[10], x2[13], x3[10], x3[13]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[11], x2[12], x3[11], x3[12]);
  x3[14] = x2[14];
  x3[15] = x2[15];
  x3[16] = _mm256_adds_epi16(x2[16], x2[23]);
  x3[23] = _mm256_subs_epi16(x2[16], x2[23]);
  x3[17] = _mm256_adds_epi16(x2[17], x2[22]);
  x3[22] = _mm256_subs_epi16(x2[17], x2[22]);
  x3[18] = _mm256_adds_epi16(x2[18], x2[21]);
  x3[21] = _mm256_subs_epi16(x2[18], x2[21]);
  x3[19] = _mm256_adds_epi16(x2[19], x2[20]);
  x3[20] = _mm256_subs_epi16(x2[19], x2[20]);
  x3[24] = _mm256_subs_epi16(x2[31], x2[24]);
  x3[31] = _mm256_adds_epi16(x2[31], x2[24]);
  x3[25] = _mm256_subs_epi16(x2[30], x2[25]);
  x3[30] = _mm256_adds_epi16(x2[30], x2[25]);
  x3[26] = _mm256_subs_epi16(x2[29], x2[26]);
  x3[29] = _mm256_adds_epi16(x2[29], x2[26]);
  x3[27] = _mm256_subs_epi16(x2[28], x2[27]);
  x3[28] = _mm256_adds_epi16(x2[28], x2[27]);

  // stage 4
  __m256i x4[32];
  x4[0] = _mm256_adds_epi16(x3[0], x3[3]);
  x4[3] = _mm256_subs_epi16(x3[0], x3[3]);
  x4[1] = _mm256_adds_epi16(x3[1], x3[2]);
  x4[2] = _mm256_subs_epi16(x3[1], x3[2]);
  x4[4] = x3[4];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x3[5], x3[6], x4[5], x4[6]);
  x4[7] = x3[7];
  x4[8] = _mm256_adds_epi16(x3[8], x3[11]);
  x4[11] = _mm256_subs_epi16(x3[8], x3[11]);
  x4[9] = _mm256_adds_epi16(x3[9], x3[10]);
  x4[10] = _mm256_subs_epi16(x3[9], x3[10]);
  x4[12] = _mm256_subs_epi16(x3[15], x3[12]);
  x4[15] = _mm256_adds_epi16(x3[15], x3[12]);
  x4[13] = _mm256_subs_epi16(x3[14], x3[13]);
  x4[14] = _mm256_adds_epi16(x3[14], x3[13]);
  x4[16] = x3[16];
  x4[17] = x3[17];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[18], x3[29], x4[18], x4[29]);
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[19], x3[28], x4[19], x4[28]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[20], x3[27], x4[20], x4[27]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[21], x3[26], x4[21], x4[26]);
  x4[22] = x3[22];
  x4[23] = x3[23];
  x4[24] = x3[24];
  x4[25] = x3[25];
  x4[30] = x3[30];
  x4[31] = x3[31];

  // stage 5
  __m256i x5[32];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x4[0], x4[1], x5[0], x5[1]);
  btf_16_w16_avx2(cospi_p48_p16, cospi_m16_p48, x4[2], x4[3], x5[2], x5[3]);
  x5[4] = _mm256_adds_epi16(x4[4], x4[5]);
  x5[5] = _mm256_subs_epi16(x4[4], x4[5]);
  x5[6] = _mm256_subs_epi16(x4[7], x4[6]);
  x5[7] = _mm256_adds_epi16(x4[7], x4[6]);
  x5[8] = x4[8];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x4[9], x4[14], x5[9], x5[14]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x4[10], x4[13], x5[10], x5[13]);
  x5[11] = x4[11];
  x5[12] = x4[12];
  x5[15] = x4[15];
  x5[16] = _mm256_adds_epi16(x4[16], x4[19]);
  x5[19] = _mm256_subs_epi16(x4[16], x4[19]);
  x5[17] = _mm256_adds_epi16(x4[17], x4[18]);
  x5[18] = _mm256_subs_epi16(x4[17], x4[18]);
  x5[20] = _mm256_subs_epi16(x4[23], x4[20]);
  x5[23] = _mm256_adds_epi16(x4[23], x4[20]);
  x5[21] = _mm256_subs_epi16(x4[22], x4[21]);
  x5[22] = _mm256_adds_epi16(x4[22], x4[21]);
  x5[24] = _mm256_adds_epi16(x4[24], x4[27]);
  x5[27] = _mm256_subs_epi16(x4[24], x4[27]);
  x5[25] = _mm256_adds_epi16(x4[25], x4[26]);
  x5[26] = _mm256_subs_epi16(x4[25], x4[26]);
  x5[28] = _mm256_subs_epi16(x4[31], x4[28]);
  x5[31] = _mm256_adds_epi16(x4[31], x4[28]);
  x5[29] = _mm256_subs_epi16(x4[30], x4[29]);
  x5[30] = _mm256_adds_epi16(x4[30], x4[29]);

  // stage 6
  __m256i x6[32];
  x6[0] = x5[0];
  x6[1] = x5[1];
  x6[2] = x5[2];
  x6[3] = x5[3];
  btf_16_w16_avx2(cospi_p56_p08, cospi_m08_p56, x5[4], x5[7], x6[4], x6[7]);
  btf_16_w16_avx2(cospi_p24_p40, cospi_m40_p24, x5[5], x5[6], x6[5], x6[6]);
  x6[8] = _mm256_adds_epi16(x5[8], x5[9]);
  x6[9] = _mm256_subs_epi16(x5[8], x5[9]);
  x6[10] = _mm256_subs_epi16(x5[11], x5[10]);
  x6[11] = _mm256_adds_epi16(x5[11], x5[10]);
  x6[12] = _mm256_adds_epi16(x5[12], x5[13]);
  x6[13] = _mm256_subs_epi16(x5[12], x5[13]);
  x6[14] = _mm256_subs_epi16(x5[15], x5[14]);
  x6[15] = _mm256_adds_epi16(x5[15], x5[14]);
  x6[16] = x5[16];
  btf_16_w16_avx2(cospi_m08_p56, cospi_p56_p08, x5[17], x5[30], x6[17], x6[30]);
  btf_16_w16_avx2(cospi_m56_m08, cospi_m08_p56, x5[18], x5[29], x6[18], x6[29]);
  x6[19] = x5[19];
  x6[20] = x5[20];
  btf_16_w16_avx2(cospi_m40_p24, cospi_p24_p40, x5[21], x5[26], x6[21], x6[26]);
  btf_16_w16_avx2(cospi_m24_m40, cospi_m40_p24, x5[22], x5[25], x6[22], x6[25]);
  x6[23] = x5[23];
  x6[24] = x5[24];
  x6[27] = x5[27];
  x6[28] = x5[28];
  x6[31] = x5[31];

  // stage 7
  __m256i x7[32];
  x7[0] = x6[0];
  x7[1] = x6[1];
  x7[2] = x6[2];
  x7[3] = x6[3];
  x7[4] = x6[4];
  x7[5] = x6[5];
  x7[6] = x6[6];
  x7[7] = x6[7];
  btf_16_w16_avx2(cospi_p60_p04, cospi_m04_p60, x6[8], x6[15], x7[8], x7[15]);
  btf_16_w16_avx2(cospi_p28_p36, cospi_m36_p28, x6[9], x6[14], x7[9], x7[14]);
  btf_16_w16_avx2(cospi_p44_p20, cospi_m20_p44, x6[10], x6[13], x7[10], x7[13]);
  btf_16_w16_avx2(cospi_p12_p52, cospi_m52_p12, x6[11], x6[12], x7[11], x7[12]);
  x7[16] = _mm256_adds_epi16(x6[16], x6[17]);
  x7[17] = _mm256_subs_epi16(x6[16], x6[17]);
  x7[18] = _mm256_subs_epi16(x6[19], x6[18]);
  x7[19] = _mm256_adds_epi16(x6[19], x6[18]);
  x7[20] = _mm256_adds_epi16(x6[20], x6[21]);
  x7[21] = _mm256_subs_epi16(x6[20], x6[21]);
  x7[22] = _mm256_subs_epi16(x6[23], x6[22]);
  x7[23] = _mm256_adds_epi16(x6[23], x6[22]);
  x7[24] = _mm256_adds_epi16(x6[24], x6[25]);
  x7[25] = _mm256_subs_epi16(x6[24], x6[25]);
  x7[26] = _mm256_subs_epi16(x6[27], x6[26]);
  x7[27] = _mm256_adds_epi16(x6[27], x6[26]);
  x7[28] = _mm256_adds_epi16(x6[28], x6[29]);
  x7[29] = _mm256_subs_epi16(x6[28], x6[29]);
  x7[30] = _mm256_subs_epi16(x6[31], x6[30]);
  x7[31] = _mm256_adds_epi16(x6[31], x6[30]);

  // stage 8
  __m256i x8[32];
  x8[0] = x7[0];
  x8[1] = x7[1];
  x8[2] = x7[2];
  x8[3] = x7[3];
  x8[4] = x7[4];
  x8[5] = x7[5];
  x8[6] = x7[6];
  x8[7] = x7[7];
  x8[8] = x7[8];
  x8[9] = x7[9];
  x8[10] = x7[10];
  x8[11] = x7[11];
  x8[12] = x7[12];
  x8[13] = x7[13];
  x8[14] = x7[14];
  x8[15] = x7[15];
  btf_16_w16_avx2(cospi_p62_p02, cospi_m02_p62, x7[16], x7[31], x8[16], x8[31]);
  btf_16_w16_avx2(cospi_p30_p34, cospi_m34_p30, x7[17], x7[30], x8[17], x8[30]);
  btf_16_w16_avx2(cospi_p46_p18, cospi_m18_p46, x7[18], x7[29], x8[18], x8[29]);
  btf_16_w16_avx2(cospi_p14_p50, cospi_m50_p14, x7[19], x7[28], x8[19], x8[28]);
  btf_16_w16_avx2(cospi_p54_p10, cospi_m10_p54, x7[20], x7[27], x8[20], x8[27]);
  btf_16_w16_avx2(cospi_p22_p42, cospi_m42_p22, x7[21], x7[26], x8[21], x8[26]);
  btf_16_w16_avx2(cospi_p38_p26, cospi_m26_p38, x7[22], x7[25], x8[22], x8[25]);
  btf_16_w16_avx2(cospi_p06_p58, cospi_m58_p06, x7[23], x7[24], x8[23], x8[24]);

  // stage 9
  output[0] = x8[0];
  output[1] = x8[16];
  output[2] = x8[8];
  output[3] = x8[24];
  output[4] = x8[4];
  output[5] = x8[20];
  output[6] = x8[12];
  output[7] = x8[28];
  output[8] = x8[2];
  output[9] = x8[18];
  output[10] = x8[10];
  output[11] = x8[26];
  output[12] = x8[6];
  output[13] = x8[22];
  output[14] = x8[14];
  output[15] = x8[30];
  output[16] = x8[1];
  output[17] = x8[17];
  output[18] = x8[9];
  output[19] = x8[25];
  output[20] = x8[5];
  output[21] = x8[21];
  output[22] = x8[13];
  output[23] = x8[29];
  output[24] = x8[3];
  output[25] = x8[19];
  output[26] = x8[11];
  output[27] = x8[27];
  output[28] = x8[7];
  output[29] = x8[23];
  output[30] = x8[15];
  output[31] = x8[31];
}

static INLINE void fdct16x64_new_avx2(const __m256i *input, __m256i *output,
                                      int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m256i __rounding = _mm256_set1_epi32(1 << (cos_bit - 1));

  __m256i cospi_m32_p32 = pair_set_w16_epi16(-cospi[32], cospi[32]);
  __m256i cospi_p32_p32 = pair_set_w16_epi16(cospi[32], cospi[32]);
  __m256i cospi_m16_p48 = pair_set_w16_epi16(-cospi[16], cospi[48]);
  __m256i cospi_p48_p16 = pair_set_w16_epi16(cospi[48], cospi[16]);
  __m256i cospi_m48_m16 = pair_set_w16_epi16(-cospi[48], -cospi[16]);
  __m256i cospi_p32_m32 = pair_set_w16_epi16(cospi[32], -cospi[32]);
  __m256i cospi_m08_p56 = pair_set_w16_epi16(-cospi[8], cospi[56]);
  __m256i cospi_p56_p08 = pair_set_w16_epi16(cospi[56], cospi[8]);
  __m256i cospi_m56_m08 = pair_set_w16_epi16(-cospi[56], -cospi[8]);
  __m256i cospi_m40_p24 = pair_set_w16_epi16(-cospi[40], cospi[24]);
  __m256i cospi_p24_p40 = pair_set_w16_epi16(cospi[24], cospi[40]);
  __m256i cospi_m24_m40 = pair_set_w16_epi16(-cospi[24], -cospi[40]);
  __m256i cospi_p60_p04 = pair_set_w16_epi16(cospi[60], cospi[4]);
  __m256i cospi_m04_p60 = pair_set_w16_epi16(-cospi[4], cospi[60]);
  __m256i cospi_p28_p36 = pair_set_w16_epi16(cospi[28], cospi[36]);
  __m256i cospi_m36_p28 = pair_set_w16_epi16(-cospi[36], cospi[28]);
  __m256i cospi_p44_p20 = pair_set_w16_epi16(cospi[44], cospi[20]);
  __m256i cospi_m20_p44 = pair_set_w16_epi16(-cospi[20], cospi[44]);
  __m256i cospi_p12_p52 = pair_set_w16_epi16(cospi[12], cospi[52]);
  __m256i cospi_m52_p12 = pair_set_w16_epi16(-cospi[52], cospi[12]);
  __m256i cospi_m60_m04 = pair_set_w16_epi16(-cospi[60], -cospi[4]);
  __m256i cospi_m28_m36 = pair_set_w16_epi16(-cospi[28], -cospi[36]);
  __m256i cospi_m44_m20 = pair_set_w16_epi16(-cospi[44], -cospi[20]);
  __m256i cospi_m12_m52 = pair_set_w16_epi16(-cospi[12], -cospi[52]);
  __m256i cospi_p62_p02 = pair_set_w16_epi16(cospi[62], cospi[2]);
  __m256i cospi_m02_p62 = pair_set_w16_epi16(-cospi[2], cospi[62]);
  __m256i cospi_p30_p34 = pair_set_w16_epi16(cospi[30], cospi[34]);
  __m256i cospi_m34_p30 = pair_set_w16_epi16(-cospi[34], cospi[30]);
  __m256i cospi_p46_p18 = pair_set_w16_epi16(cospi[46], cospi[18]);
  __m256i cospi_m18_p46 = pair_set_w16_epi16(-cospi[18], cospi[46]);
  __m256i cospi_p14_p50 = pair_set_w16_epi16(cospi[14], cospi[50]);
  __m256i cospi_m50_p14 = pair_set_w16_epi16(-cospi[50], cospi[14]);
  __m256i cospi_p54_p10 = pair_set_w16_epi16(cospi[54], cospi[10]);
  __m256i cospi_m10_p54 = pair_set_w16_epi16(-cospi[10], cospi[54]);
  __m256i cospi_p22_p42 = pair_set_w16_epi16(cospi[22], cospi[42]);
  __m256i cospi_m42_p22 = pair_set_w16_epi16(-cospi[42], cospi[22]);
  __m256i cospi_p38_p26 = pair_set_w16_epi16(cospi[38], cospi[26]);
  __m256i cospi_m26_p38 = pair_set_w16_epi16(-cospi[26], cospi[38]);
  __m256i cospi_p06_p58 = pair_set_w16_epi16(cospi[6], cospi[58]);
  __m256i cospi_m58_p06 = pair_set_w16_epi16(-cospi[58], cospi[6]);
  __m256i cospi_p63_p01 = pair_set_w16_epi16(cospi[63], cospi[1]);
  __m256i cospi_m01_p63 = pair_set_w16_epi16(-cospi[1], cospi[63]);
  __m256i cospi_p31_p33 = pair_set_w16_epi16(cospi[31], cospi[33]);
  __m256i cospi_m33_p31 = pair_set_w16_epi16(-cospi[33], cospi[31]);
  __m256i cospi_p47_p17 = pair_set_w16_epi16(cospi[47], cospi[17]);
  __m256i cospi_m17_p47 = pair_set_w16_epi16(-cospi[17], cospi[47]);
  __m256i cospi_p15_p49 = pair_set_w16_epi16(cospi[15], cospi[49]);
  __m256i cospi_m49_p15 = pair_set_w16_epi16(-cospi[49], cospi[15]);
  __m256i cospi_p55_p09 = pair_set_w16_epi16(cospi[55], cospi[9]);
  __m256i cospi_m09_p55 = pair_set_w16_epi16(-cospi[9], cospi[55]);
  __m256i cospi_p23_p41 = pair_set_w16_epi16(cospi[23], cospi[41]);
  __m256i cospi_m41_p23 = pair_set_w16_epi16(-cospi[41], cospi[23]);
  __m256i cospi_p39_p25 = pair_set_w16_epi16(cospi[39], cospi[25]);
  __m256i cospi_m25_p39 = pair_set_w16_epi16(-cospi[25], cospi[39]);
  __m256i cospi_p07_p57 = pair_set_w16_epi16(cospi[7], cospi[57]);
  __m256i cospi_m57_p07 = pair_set_w16_epi16(-cospi[57], cospi[7]);
  __m256i cospi_p59_p05 = pair_set_w16_epi16(cospi[59], cospi[5]);
  __m256i cospi_m05_p59 = pair_set_w16_epi16(-cospi[5], cospi[59]);
  __m256i cospi_p27_p37 = pair_set_w16_epi16(cospi[27], cospi[37]);
  __m256i cospi_m37_p27 = pair_set_w16_epi16(-cospi[37], cospi[27]);
  __m256i cospi_p43_p21 = pair_set_w16_epi16(cospi[43], cospi[21]);
  __m256i cospi_m21_p43 = pair_set_w16_epi16(-cospi[21], cospi[43]);
  __m256i cospi_p11_p53 = pair_set_w16_epi16(cospi[11], cospi[53]);
  __m256i cospi_m53_p11 = pair_set_w16_epi16(-cospi[53], cospi[11]);
  __m256i cospi_p51_p13 = pair_set_w16_epi16(cospi[51], cospi[13]);
  __m256i cospi_m13_p51 = pair_set_w16_epi16(-cospi[13], cospi[51]);
  __m256i cospi_p19_p45 = pair_set_w16_epi16(cospi[19], cospi[45]);
  __m256i cospi_m45_p19 = pair_set_w16_epi16(-cospi[45], cospi[19]);
  __m256i cospi_p35_p29 = pair_set_w16_epi16(cospi[35], cospi[29]);
  __m256i cospi_m29_p35 = pair_set_w16_epi16(-cospi[29], cospi[35]);
  __m256i cospi_p03_p61 = pair_set_w16_epi16(cospi[3], cospi[61]);
  __m256i cospi_m61_p03 = pair_set_w16_epi16(-cospi[61], cospi[3]);

  // stage 1
  __m256i x1[64];
  x1[0] = _mm256_adds_epi16(input[0], input[63]);
  x1[63] = _mm256_subs_epi16(input[0], input[63]);
  x1[1] = _mm256_adds_epi16(input[1], input[62]);
  x1[62] = _mm256_subs_epi16(input[1], input[62]);
  x1[2] = _mm256_adds_epi16(input[2], input[61]);
  x1[61] = _mm256_subs_epi16(input[2], input[61]);
  x1[3] = _mm256_adds_epi16(input[3], input[60]);
  x1[60] = _mm256_subs_epi16(input[3], input[60]);
  x1[4] = _mm256_adds_epi16(input[4], input[59]);
  x1[59] = _mm256_subs_epi16(input[4], input[59]);
  x1[5] = _mm256_adds_epi16(input[5], input[58]);
  x1[58] = _mm256_subs_epi16(input[5], input[58]);
  x1[6] = _mm256_adds_epi16(input[6], input[57]);
  x1[57] = _mm256_subs_epi16(input[6], input[57]);
  x1[7] = _mm256_adds_epi16(input[7], input[56]);
  x1[56] = _mm256_subs_epi16(input[7], input[56]);
  x1[8] = _mm256_adds_epi16(input[8], input[55]);
  x1[55] = _mm256_subs_epi16(input[8], input[55]);
  x1[9] = _mm256_adds_epi16(input[9], input[54]);
  x1[54] = _mm256_subs_epi16(input[9], input[54]);
  x1[10] = _mm256_adds_epi16(input[10], input[53]);
  x1[53] = _mm256_subs_epi16(input[10], input[53]);
  x1[11] = _mm256_adds_epi16(input[11], input[52]);
  x1[52] = _mm256_subs_epi16(input[11], input[52]);
  x1[12] = _mm256_adds_epi16(input[12], input[51]);
  x1[51] = _mm256_subs_epi16(input[12], input[51]);
  x1[13] = _mm256_adds_epi16(input[13], input[50]);
  x1[50] = _mm256_subs_epi16(input[13], input[50]);
  x1[14] = _mm256_adds_epi16(input[14], input[49]);
  x1[49] = _mm256_subs_epi16(input[14], input[49]);
  x1[15] = _mm256_adds_epi16(input[15], input[48]);
  x1[48] = _mm256_subs_epi16(input[15], input[48]);
  x1[16] = _mm256_adds_epi16(input[16], input[47]);
  x1[47] = _mm256_subs_epi16(input[16], input[47]);
  x1[17] = _mm256_adds_epi16(input[17], input[46]);
  x1[46] = _mm256_subs_epi16(input[17], input[46]);
  x1[18] = _mm256_adds_epi16(input[18], input[45]);
  x1[45] = _mm256_subs_epi16(input[18], input[45]);
  x1[19] = _mm256_adds_epi16(input[19], input[44]);
  x1[44] = _mm256_subs_epi16(input[19], input[44]);
  x1[20] = _mm256_adds_epi16(input[20], input[43]);
  x1[43] = _mm256_subs_epi16(input[20], input[43]);
  x1[21] = _mm256_adds_epi16(input[21], input[42]);
  x1[42] = _mm256_subs_epi16(input[21], input[42]);
  x1[22] = _mm256_adds_epi16(input[22], input[41]);
  x1[41] = _mm256_subs_epi16(input[22], input[41]);
  x1[23] = _mm256_adds_epi16(input[23], input[40]);
  x1[40] = _mm256_subs_epi16(input[23], input[40]);
  x1[24] = _mm256_adds_epi16(input[24], input[39]);
  x1[39] = _mm256_subs_epi16(input[24], input[39]);
  x1[25] = _mm256_adds_epi16(input[25], input[38]);
  x1[38] = _mm256_subs_epi16(input[25], input[38]);
  x1[26] = _mm256_adds_epi16(input[26], input[37]);
  x1[37] = _mm256_subs_epi16(input[26], input[37]);
  x1[27] = _mm256_adds_epi16(input[27], input[36]);
  x1[36] = _mm256_subs_epi16(input[27], input[36]);
  x1[28] = _mm256_adds_epi16(input[28], input[35]);
  x1[35] = _mm256_subs_epi16(input[28], input[35]);
  x1[29] = _mm256_adds_epi16(input[29], input[34]);
  x1[34] = _mm256_subs_epi16(input[29], input[34]);
  x1[30] = _mm256_adds_epi16(input[30], input[33]);
  x1[33] = _mm256_subs_epi16(input[30], input[33]);
  x1[31] = _mm256_adds_epi16(input[31], input[32]);
  x1[32] = _mm256_subs_epi16(input[31], input[32]);

  // stage 2
  __m256i x2[64];
  x2[0] = _mm256_adds_epi16(x1[0], x1[31]);
  x2[31] = _mm256_subs_epi16(x1[0], x1[31]);
  x2[1] = _mm256_adds_epi16(x1[1], x1[30]);
  x2[30] = _mm256_subs_epi16(x1[1], x1[30]);
  x2[2] = _mm256_adds_epi16(x1[2], x1[29]);
  x2[29] = _mm256_subs_epi16(x1[2], x1[29]);
  x2[3] = _mm256_adds_epi16(x1[3], x1[28]);
  x2[28] = _mm256_subs_epi16(x1[3], x1[28]);
  x2[4] = _mm256_adds_epi16(x1[4], x1[27]);
  x2[27] = _mm256_subs_epi16(x1[4], x1[27]);
  x2[5] = _mm256_adds_epi16(x1[5], x1[26]);
  x2[26] = _mm256_subs_epi16(x1[5], x1[26]);
  x2[6] = _mm256_adds_epi16(x1[6], x1[25]);
  x2[25] = _mm256_subs_epi16(x1[6], x1[25]);
  x2[7] = _mm256_adds_epi16(x1[7], x1[24]);
  x2[24] = _mm256_subs_epi16(x1[7], x1[24]);
  x2[8] = _mm256_adds_epi16(x1[8], x1[23]);
  x2[23] = _mm256_subs_epi16(x1[8], x1[23]);
  x2[9] = _mm256_adds_epi16(x1[9], x1[22]);
  x2[22] = _mm256_subs_epi16(x1[9], x1[22]);
  x2[10] = _mm256_adds_epi16(x1[10], x1[21]);
  x2[21] = _mm256_subs_epi16(x1[10], x1[21]);
  x2[11] = _mm256_adds_epi16(x1[11], x1[20]);
  x2[20] = _mm256_subs_epi16(x1[11], x1[20]);
  x2[12] = _mm256_adds_epi16(x1[12], x1[19]);
  x2[19] = _mm256_subs_epi16(x1[12], x1[19]);
  x2[13] = _mm256_adds_epi16(x1[13], x1[18]);
  x2[18] = _mm256_subs_epi16(x1[13], x1[18]);
  x2[14] = _mm256_adds_epi16(x1[14], x1[17]);
  x2[17] = _mm256_subs_epi16(x1[14], x1[17]);
  x2[15] = _mm256_adds_epi16(x1[15], x1[16]);
  x2[16] = _mm256_subs_epi16(x1[15], x1[16]);
  x2[32] = x1[32];
  x2[33] = x1[33];
  x2[34] = x1[34];
  x2[35] = x1[35];
  x2[36] = x1[36];
  x2[37] = x1[37];
  x2[38] = x1[38];
  x2[39] = x1[39];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[40], x1[55], x2[40], x2[55]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[41], x1[54], x2[41], x2[54]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[42], x1[53], x2[42], x2[53]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[43], x1[52], x2[43], x2[52]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[44], x1[51], x2[44], x2[51]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[45], x1[50], x2[45], x2[50]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[46], x1[49], x2[46], x2[49]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x1[47], x1[48], x2[47], x2[48]);
  x2[56] = x1[56];
  x2[57] = x1[57];
  x2[58] = x1[58];
  x2[59] = x1[59];
  x2[60] = x1[60];
  x2[61] = x1[61];
  x2[62] = x1[62];
  x2[63] = x1[63];

  // stage 3
  __m256i x3[64];
  x3[0] = _mm256_adds_epi16(x2[0], x2[15]);
  x3[15] = _mm256_subs_epi16(x2[0], x2[15]);
  x3[1] = _mm256_adds_epi16(x2[1], x2[14]);
  x3[14] = _mm256_subs_epi16(x2[1], x2[14]);
  x3[2] = _mm256_adds_epi16(x2[2], x2[13]);
  x3[13] = _mm256_subs_epi16(x2[2], x2[13]);
  x3[3] = _mm256_adds_epi16(x2[3], x2[12]);
  x3[12] = _mm256_subs_epi16(x2[3], x2[12]);
  x3[4] = _mm256_adds_epi16(x2[4], x2[11]);
  x3[11] = _mm256_subs_epi16(x2[4], x2[11]);
  x3[5] = _mm256_adds_epi16(x2[5], x2[10]);
  x3[10] = _mm256_subs_epi16(x2[5], x2[10]);
  x3[6] = _mm256_adds_epi16(x2[6], x2[9]);
  x3[9] = _mm256_subs_epi16(x2[6], x2[9]);
  x3[7] = _mm256_adds_epi16(x2[7], x2[8]);
  x3[8] = _mm256_subs_epi16(x2[7], x2[8]);
  x3[16] = x2[16];
  x3[17] = x2[17];
  x3[18] = x2[18];
  x3[19] = x2[19];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[20], x2[27], x3[20], x3[27]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[21], x2[26], x3[21], x3[26]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[22], x2[25], x3[22], x3[25]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x2[23], x2[24], x3[23], x3[24]);
  x3[28] = x2[28];
  x3[29] = x2[29];
  x3[30] = x2[30];
  x3[31] = x2[31];
  x3[32] = _mm256_adds_epi16(x2[32], x2[47]);
  x3[47] = _mm256_subs_epi16(x2[32], x2[47]);
  x3[33] = _mm256_adds_epi16(x2[33], x2[46]);
  x3[46] = _mm256_subs_epi16(x2[33], x2[46]);
  x3[34] = _mm256_adds_epi16(x2[34], x2[45]);
  x3[45] = _mm256_subs_epi16(x2[34], x2[45]);
  x3[35] = _mm256_adds_epi16(x2[35], x2[44]);
  x3[44] = _mm256_subs_epi16(x2[35], x2[44]);
  x3[36] = _mm256_adds_epi16(x2[36], x2[43]);
  x3[43] = _mm256_subs_epi16(x2[36], x2[43]);
  x3[37] = _mm256_adds_epi16(x2[37], x2[42]);
  x3[42] = _mm256_subs_epi16(x2[37], x2[42]);
  x3[38] = _mm256_adds_epi16(x2[38], x2[41]);
  x3[41] = _mm256_subs_epi16(x2[38], x2[41]);
  x3[39] = _mm256_adds_epi16(x2[39], x2[40]);
  x3[40] = _mm256_subs_epi16(x2[39], x2[40]);
  x3[48] = _mm256_subs_epi16(x2[63], x2[48]);
  x3[63] = _mm256_adds_epi16(x2[63], x2[48]);
  x3[49] = _mm256_subs_epi16(x2[62], x2[49]);
  x3[62] = _mm256_adds_epi16(x2[62], x2[49]);
  x3[50] = _mm256_subs_epi16(x2[61], x2[50]);
  x3[61] = _mm256_adds_epi16(x2[61], x2[50]);
  x3[51] = _mm256_subs_epi16(x2[60], x2[51]);
  x3[60] = _mm256_adds_epi16(x2[60], x2[51]);
  x3[52] = _mm256_subs_epi16(x2[59], x2[52]);
  x3[59] = _mm256_adds_epi16(x2[59], x2[52]);
  x3[53] = _mm256_subs_epi16(x2[58], x2[53]);
  x3[58] = _mm256_adds_epi16(x2[58], x2[53]);
  x3[54] = _mm256_subs_epi16(x2[57], x2[54]);
  x3[57] = _mm256_adds_epi16(x2[57], x2[54]);
  x3[55] = _mm256_subs_epi16(x2[56], x2[55]);
  x3[56] = _mm256_adds_epi16(x2[56], x2[55]);

  // stage 4
  __m256i x4[64];
  x4[0] = _mm256_adds_epi16(x3[0], x3[7]);
  x4[7] = _mm256_subs_epi16(x3[0], x3[7]);
  x4[1] = _mm256_adds_epi16(x3[1], x3[6]);
  x4[6] = _mm256_subs_epi16(x3[1], x3[6]);
  x4[2] = _mm256_adds_epi16(x3[2], x3[5]);
  x4[5] = _mm256_subs_epi16(x3[2], x3[5]);
  x4[3] = _mm256_adds_epi16(x3[3], x3[4]);
  x4[4] = _mm256_subs_epi16(x3[3], x3[4]);
  x4[8] = x3[8];
  x4[9] = x3[9];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x3[10], x3[13], x4[10], x4[13]);
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x3[11], x3[12], x4[11], x4[12]);
  x4[14] = x3[14];
  x4[15] = x3[15];
  x4[16] = _mm256_adds_epi16(x3[16], x3[23]);
  x4[23] = _mm256_subs_epi16(x3[16], x3[23]);
  x4[17] = _mm256_adds_epi16(x3[17], x3[22]);
  x4[22] = _mm256_subs_epi16(x3[17], x3[22]);
  x4[18] = _mm256_adds_epi16(x3[18], x3[21]);
  x4[21] = _mm256_subs_epi16(x3[18], x3[21]);
  x4[19] = _mm256_adds_epi16(x3[19], x3[20]);
  x4[20] = _mm256_subs_epi16(x3[19], x3[20]);
  x4[24] = _mm256_subs_epi16(x3[31], x3[24]);
  x4[31] = _mm256_adds_epi16(x3[31], x3[24]);
  x4[25] = _mm256_subs_epi16(x3[30], x3[25]);
  x4[30] = _mm256_adds_epi16(x3[30], x3[25]);
  x4[26] = _mm256_subs_epi16(x3[29], x3[26]);
  x4[29] = _mm256_adds_epi16(x3[29], x3[26]);
  x4[27] = _mm256_subs_epi16(x3[28], x3[27]);
  x4[28] = _mm256_adds_epi16(x3[28], x3[27]);
  x4[32] = x3[32];
  x4[33] = x3[33];
  x4[34] = x3[34];
  x4[35] = x3[35];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[36], x3[59], x4[36], x4[59]);
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[37], x3[58], x4[37], x4[58]);
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[38], x3[57], x4[38], x4[57]);
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x3[39], x3[56], x4[39], x4[56]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[40], x3[55], x4[40], x4[55]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[41], x3[54], x4[41], x4[54]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[42], x3[53], x4[42], x4[53]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x3[43], x3[52], x4[43], x4[52]);
  x4[44] = x3[44];
  x4[45] = x3[45];
  x4[46] = x3[46];
  x4[47] = x3[47];
  x4[48] = x3[48];
  x4[49] = x3[49];
  x4[50] = x3[50];
  x4[51] = x3[51];
  x4[60] = x3[60];
  x4[61] = x3[61];
  x4[62] = x3[62];
  x4[63] = x3[63];

  // stage 5
  __m256i x5[64];
  x5[0] = _mm256_adds_epi16(x4[0], x4[3]);
  x5[3] = _mm256_subs_epi16(x4[0], x4[3]);
  x5[1] = _mm256_adds_epi16(x4[1], x4[2]);
  x5[2] = _mm256_subs_epi16(x4[1], x4[2]);
  x5[4] = x4[4];
  btf_16_w16_avx2(cospi_m32_p32, cospi_p32_p32, x4[5], x4[6], x5[5], x5[6]);
  x5[7] = x4[7];
  x5[8] = _mm256_adds_epi16(x4[8], x4[11]);
  x5[11] = _mm256_subs_epi16(x4[8], x4[11]);
  x5[9] = _mm256_adds_epi16(x4[9], x4[10]);
  x5[10] = _mm256_subs_epi16(x4[9], x4[10]);
  x5[12] = _mm256_subs_epi16(x4[15], x4[12]);
  x5[15] = _mm256_adds_epi16(x4[15], x4[12]);
  x5[13] = _mm256_subs_epi16(x4[14], x4[13]);
  x5[14] = _mm256_adds_epi16(x4[14], x4[13]);
  x5[16] = x4[16];
  x5[17] = x4[17];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x4[18], x4[29], x5[18], x5[29]);
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x4[19], x4[28], x5[19], x5[28]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x4[20], x4[27], x5[20], x5[27]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x4[21], x4[26], x5[21], x5[26]);
  x5[22] = x4[22];
  x5[23] = x4[23];
  x5[24] = x4[24];
  x5[25] = x4[25];
  x5[30] = x4[30];
  x5[31] = x4[31];
  x5[32] = _mm256_adds_epi16(x4[32], x4[39]);
  x5[39] = _mm256_subs_epi16(x4[32], x4[39]);
  x5[33] = _mm256_adds_epi16(x4[33], x4[38]);
  x5[38] = _mm256_subs_epi16(x4[33], x4[38]);
  x5[34] = _mm256_adds_epi16(x4[34], x4[37]);
  x5[37] = _mm256_subs_epi16(x4[34], x4[37]);
  x5[35] = _mm256_adds_epi16(x4[35], x4[36]);
  x5[36] = _mm256_subs_epi16(x4[35], x4[36]);
  x5[40] = _mm256_subs_epi16(x4[47], x4[40]);
  x5[47] = _mm256_adds_epi16(x4[47], x4[40]);
  x5[41] = _mm256_subs_epi16(x4[46], x4[41]);
  x5[46] = _mm256_adds_epi16(x4[46], x4[41]);
  x5[42] = _mm256_subs_epi16(x4[45], x4[42]);
  x5[45] = _mm256_adds_epi16(x4[45], x4[42]);
  x5[43] = _mm256_subs_epi16(x4[44], x4[43]);
  x5[44] = _mm256_adds_epi16(x4[44], x4[43]);
  x5[48] = _mm256_adds_epi16(x4[48], x4[55]);
  x5[55] = _mm256_subs_epi16(x4[48], x4[55]);
  x5[49] = _mm256_adds_epi16(x4[49], x4[54]);
  x5[54] = _mm256_subs_epi16(x4[49], x4[54]);
  x5[50] = _mm256_adds_epi16(x4[50], x4[53]);
  x5[53] = _mm256_subs_epi16(x4[50], x4[53]);
  x5[51] = _mm256_adds_epi16(x4[51], x4[52]);
  x5[52] = _mm256_subs_epi16(x4[51], x4[52]);
  x5[56] = _mm256_subs_epi16(x4[63], x4[56]);
  x5[63] = _mm256_adds_epi16(x4[63], x4[56]);
  x5[57] = _mm256_subs_epi16(x4[62], x4[57]);
  x5[62] = _mm256_adds_epi16(x4[62], x4[57]);
  x5[58] = _mm256_subs_epi16(x4[61], x4[58]);
  x5[61] = _mm256_adds_epi16(x4[61], x4[58]);
  x5[59] = _mm256_subs_epi16(x4[60], x4[59]);
  x5[60] = _mm256_adds_epi16(x4[60], x4[59]);

  // stage 6
  __m256i x6[64];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x5[0], x5[1], x6[0], x6[1]);
  btf_16_w16_avx2(cospi_p48_p16, cospi_m16_p48, x5[2], x5[3], x6[2], x6[3]);
  x6[4] = _mm256_adds_epi16(x5[4], x5[5]);
  x6[5] = _mm256_subs_epi16(x5[4], x5[5]);
  x6[6] = _mm256_subs_epi16(x5[7], x5[6]);
  x6[7] = _mm256_adds_epi16(x5[7], x5[6]);
  x6[8] = x5[8];
  btf_16_w16_avx2(cospi_m16_p48, cospi_p48_p16, x5[9], x5[14], x6[9], x6[14]);
  btf_16_w16_avx2(cospi_m48_m16, cospi_m16_p48, x5[10], x5[13], x6[10], x6[13]);
  x6[11] = x5[11];
  x6[12] = x5[12];
  x6[15] = x5[15];
  x6[16] = _mm256_adds_epi16(x5[16], x5[19]);
  x6[19] = _mm256_subs_epi16(x5[16], x5[19]);
  x6[17] = _mm256_adds_epi16(x5[17], x5[18]);
  x6[18] = _mm256_subs_epi16(x5[17], x5[18]);
  x6[20] = _mm256_subs_epi16(x5[23], x5[20]);
  x6[23] = _mm256_adds_epi16(x5[23], x5[20]);
  x6[21] = _mm256_subs_epi16(x5[22], x5[21]);
  x6[22] = _mm256_adds_epi16(x5[22], x5[21]);
  x6[24] = _mm256_adds_epi16(x5[24], x5[27]);
  x6[27] = _mm256_subs_epi16(x5[24], x5[27]);
  x6[25] = _mm256_adds_epi16(x5[25], x5[26]);
  x6[26] = _mm256_subs_epi16(x5[25], x5[26]);
  x6[28] = _mm256_subs_epi16(x5[31], x5[28]);
  x6[31] = _mm256_adds_epi16(x5[31], x5[28]);
  x6[29] = _mm256_subs_epi16(x5[30], x5[29]);
  x6[30] = _mm256_adds_epi16(x5[30], x5[29]);
  x6[32] = x5[32];
  x6[33] = x5[33];
  btf_16_w16_avx2(cospi_m08_p56, cospi_p56_p08, x5[34], x5[61], x6[34], x6[61]);
  btf_16_w16_avx2(cospi_m08_p56, cospi_p56_p08, x5[35], x5[60], x6[35], x6[60]);
  btf_16_w16_avx2(cospi_m56_m08, cospi_m08_p56, x5[36], x5[59], x6[36], x6[59]);
  btf_16_w16_avx2(cospi_m56_m08, cospi_m08_p56, x5[37], x5[58], x6[37], x6[58]);
  x6[38] = x5[38];
  x6[39] = x5[39];
  x6[40] = x5[40];
  x6[41] = x5[41];
  btf_16_w16_avx2(cospi_m40_p24, cospi_p24_p40, x5[42], x5[53], x6[42], x6[53]);
  btf_16_w16_avx2(cospi_m40_p24, cospi_p24_p40, x5[43], x5[52], x6[43], x6[52]);
  btf_16_w16_avx2(cospi_m24_m40, cospi_m40_p24, x5[44], x5[51], x6[44], x6[51]);
  btf_16_w16_avx2(cospi_m24_m40, cospi_m40_p24, x5[45], x5[50], x6[45], x6[50]);
  x6[46] = x5[46];
  x6[47] = x5[47];
  x6[48] = x5[48];
  x6[49] = x5[49];
  x6[54] = x5[54];
  x6[55] = x5[55];
  x6[56] = x5[56];
  x6[57] = x5[57];
  x6[62] = x5[62];
  x6[63] = x5[63];

  // stage 7
  __m256i x7[64];
  x7[0] = x6[0];
  x7[1] = x6[1];
  x7[2] = x6[2];
  x7[3] = x6[3];
  btf_16_w16_avx2(cospi_p56_p08, cospi_m08_p56, x6[4], x6[7], x7[4], x7[7]);
  btf_16_w16_avx2(cospi_p24_p40, cospi_m40_p24, x6[5], x6[6], x7[5], x7[6]);
  x7[8] = _mm256_adds_epi16(x6[8], x6[9]);
  x7[9] = _mm256_subs_epi16(x6[8], x6[9]);
  x7[10] = _mm256_subs_epi16(x6[11], x6[10]);
  x7[11] = _mm256_adds_epi16(x6[11], x6[10]);
  x7[12] = _mm256_adds_epi16(x6[12], x6[13]);
  x7[13] = _mm256_subs_epi16(x6[12], x6[13]);
  x7[14] = _mm256_subs_epi16(x6[15], x6[14]);
  x7[15] = _mm256_adds_epi16(x6[15], x6[14]);
  x7[16] = x6[16];
  btf_16_w16_avx2(cospi_m08_p56, cospi_p56_p08, x6[17], x6[30], x7[17], x7[30]);
  btf_16_w16_avx2(cospi_m56_m08, cospi_m08_p56, x6[18], x6[29], x7[18], x7[29]);
  x7[19] = x6[19];
  x7[20] = x6[20];
  btf_16_w16_avx2(cospi_m40_p24, cospi_p24_p40, x6[21], x6[26], x7[21], x7[26]);
  btf_16_w16_avx2(cospi_m24_m40, cospi_m40_p24, x6[22], x6[25], x7[22], x7[25]);
  x7[23] = x6[23];
  x7[24] = x6[24];
  x7[27] = x6[27];
  x7[28] = x6[28];
  x7[31] = x6[31];
  x7[32] = _mm256_adds_epi16(x6[32], x6[35]);
  x7[35] = _mm256_subs_epi16(x6[32], x6[35]);
  x7[33] = _mm256_adds_epi16(x6[33], x6[34]);
  x7[34] = _mm256_subs_epi16(x6[33], x6[34]);
  x7[36] = _mm256_subs_epi16(x6[39], x6[36]);
  x7[39] = _mm256_adds_epi16(x6[39], x6[36]);
  x7[37] = _mm256_subs_epi16(x6[38], x6[37]);
  x7[38] = _mm256_adds_epi16(x6[38], x6[37]);
  x7[40] = _mm256_adds_epi16(x6[40], x6[43]);
  x7[43] = _mm256_subs_epi16(x6[40], x6[43]);
  x7[41] = _mm256_adds_epi16(x6[41], x6[42]);
  x7[42] = _mm256_subs_epi16(x6[41], x6[42]);
  x7[44] = _mm256_subs_epi16(x6[47], x6[44]);
  x7[47] = _mm256_adds_epi16(x6[47], x6[44]);
  x7[45] = _mm256_subs_epi16(x6[46], x6[45]);
  x7[46] = _mm256_adds_epi16(x6[46], x6[45]);
  x7[48] = _mm256_adds_epi16(x6[48], x6[51]);
  x7[51] = _mm256_subs_epi16(x6[48], x6[51]);
  x7[49] = _mm256_adds_epi16(x6[49], x6[50]);
  x7[50] = _mm256_subs_epi16(x6[49], x6[50]);
  x7[52] = _mm256_subs_epi16(x6[55], x6[52]);
  x7[55] = _mm256_adds_epi16(x6[55], x6[52]);
  x7[53] = _mm256_subs_epi16(x6[54], x6[53]);
  x7[54] = _mm256_adds_epi16(x6[54], x6[53]);
  x7[56] = _mm256_adds_epi16(x6[56], x6[59]);
  x7[59] = _mm256_subs_epi16(x6[56], x6[59]);
  x7[57] = _mm256_adds_epi16(x6[57], x6[58]);
  x7[58] = _mm256_subs_epi16(x6[57], x6[58]);
  x7[60] = _mm256_subs_epi16(x6[63], x6[60]);
  x7[63] = _mm256_adds_epi16(x6[63], x6[60]);
  x7[61] = _mm256_subs_epi16(x6[62], x6[61]);
  x7[62] = _mm256_adds_epi16(x6[62], x6[61]);

  // stage 8
  __m256i x8[64];
  x8[0] = x7[0];
  x8[1] = x7[1];
  x8[2] = x7[2];
  x8[3] = x7[3];
  x8[4] = x7[4];
  x8[5] = x7[5];
  x8[6] = x7[6];
  x8[7] = x7[7];
  btf_16_w16_avx2(cospi_p60_p04, cospi_m04_p60, x7[8], x7[15], x8[8], x8[15]);
  btf_16_w16_avx2(cospi_p28_p36, cospi_m36_p28, x7[9], x7[14], x8[9], x8[14]);
  btf_16_w16_avx2(cospi_p44_p20, cospi_m20_p44, x7[10], x7[13], x8[10], x8[13]);
  btf_16_w16_avx2(cospi_p12_p52, cospi_m52_p12, x7[11], x7[12], x8[11], x8[12]);
  x8[16] = _mm256_adds_epi16(x7[16], x7[17]);
  x8[17] = _mm256_subs_epi16(x7[16], x7[17]);
  x8[18] = _mm256_subs_epi16(x7[19], x7[18]);
  x8[19] = _mm256_adds_epi16(x7[19], x7[18]);
  x8[20] = _mm256_adds_epi16(x7[20], x7[21]);
  x8[21] = _mm256_subs_epi16(x7[20], x7[21]);
  x8[22] = _mm256_subs_epi16(x7[23], x7[22]);
  x8[23] = _mm256_adds_epi16(x7[23], x7[22]);
  x8[24] = _mm256_adds_epi16(x7[24], x7[25]);
  x8[25] = _mm256_subs_epi16(x7[24], x7[25]);
  x8[26] = _mm256_subs_epi16(x7[27], x7[26]);
  x8[27] = _mm256_adds_epi16(x7[27], x7[26]);
  x8[28] = _mm256_adds_epi16(x7[28], x7[29]);
  x8[29] = _mm256_subs_epi16(x7[28], x7[29]);
  x8[30] = _mm256_subs_epi16(x7[31], x7[30]);
  x8[31] = _mm256_adds_epi16(x7[31], x7[30]);
  x8[32] = x7[32];
  btf_16_w16_avx2(cospi_m04_p60, cospi_p60_p04, x7[33], x7[62], x8[33], x8[62]);
  btf_16_w16_avx2(cospi_m60_m04, cospi_m04_p60, x7[34], x7[61], x8[34], x8[61]);
  x8[35] = x7[35];
  x8[36] = x7[36];
  btf_16_w16_avx2(cospi_m36_p28, cospi_p28_p36, x7[37], x7[58], x8[37], x8[58]);
  btf_16_w16_avx2(cospi_m28_m36, cospi_m36_p28, x7[38], x7[57], x8[38], x8[57]);
  x8[39] = x7[39];
  x8[40] = x7[40];
  btf_16_w16_avx2(cospi_m20_p44, cospi_p44_p20, x7[41], x7[54], x8[41], x8[54]);
  btf_16_w16_avx2(cospi_m44_m20, cospi_m20_p44, x7[42], x7[53], x8[42], x8[53]);
  x8[43] = x7[43];
  x8[44] = x7[44];
  btf_16_w16_avx2(cospi_m52_p12, cospi_p12_p52, x7[45], x7[50], x8[45], x8[50]);
  btf_16_w16_avx2(cospi_m12_m52, cospi_m52_p12, x7[46], x7[49], x8[46], x8[49]);
  x8[47] = x7[47];
  x8[48] = x7[48];
  x8[51] = x7[51];
  x8[52] = x7[52];
  x8[55] = x7[55];
  x8[56] = x7[56];
  x8[59] = x7[59];
  x8[60] = x7[60];
  x8[63] = x7[63];

  // stage 9
  __m256i x9[64];
  x9[0] = x8[0];
  x9[1] = x8[1];
  x9[2] = x8[2];
  x9[3] = x8[3];
  x9[4] = x8[4];
  x9[5] = x8[5];
  x9[6] = x8[6];
  x9[7] = x8[7];
  x9[8] = x8[8];
  x9[9] = x8[9];
  x9[10] = x8[10];
  x9[11] = x8[11];
  x9[12] = x8[12];
  x9[13] = x8[13];
  x9[14] = x8[14];
  x9[15] = x8[15];
  btf_16_w16_avx2(cospi_p62_p02, cospi_m02_p62, x8[16], x8[31], x9[16], x9[31]);
  btf_16_w16_avx2(cospi_p30_p34, cospi_m34_p30, x8[17], x8[30], x9[17], x9[30]);
  btf_16_w16_avx2(cospi_p46_p18, cospi_m18_p46, x8[18], x8[29], x9[18], x9[29]);
  btf_16_w16_avx2(cospi_p14_p50, cospi_m50_p14, x8[19], x8[28], x9[19], x9[28]);
  btf_16_w16_avx2(cospi_p54_p10, cospi_m10_p54, x8[20], x8[27], x9[20], x9[27]);
  btf_16_w16_avx2(cospi_p22_p42, cospi_m42_p22, x8[21], x8[26], x9[21], x9[26]);
  btf_16_w16_avx2(cospi_p38_p26, cospi_m26_p38, x8[22], x8[25], x9[22], x9[25]);
  btf_16_w16_avx2(cospi_p06_p58, cospi_m58_p06, x8[23], x8[24], x9[23], x9[24]);
  x9[32] = _mm256_adds_epi16(x8[32], x8[33]);
  x9[33] = _mm256_subs_epi16(x8[32], x8[33]);
  x9[34] = _mm256_subs_epi16(x8[35], x8[34]);
  x9[35] = _mm256_adds_epi16(x8[35], x8[34]);
  x9[36] = _mm256_adds_epi16(x8[36], x8[37]);
  x9[37] = _mm256_subs_epi16(x8[36], x8[37]);
  x9[38] = _mm256_subs_epi16(x8[39], x8[38]);
  x9[39] = _mm256_adds_epi16(x8[39], x8[38]);
  x9[40] = _mm256_adds_epi16(x8[40], x8[41]);
  x9[41] = _mm256_subs_epi16(x8[40], x8[41]);
  x9[42] = _mm256_subs_epi16(x8[43], x8[42]);
  x9[43] = _mm256_adds_epi16(x8[43], x8[42]);
  x9[44] = _mm256_adds_epi16(x8[44], x8[45]);
  x9[45] = _mm256_subs_epi16(x8[44], x8[45]);
  x9[46] = _mm256_subs_epi16(x8[47], x8[46]);
  x9[47] = _mm256_adds_epi16(x8[47], x8[46]);
  x9[48] = _mm256_adds_epi16(x8[48], x8[49]);
  x9[49] = _mm256_subs_epi16(x8[48], x8[49]);
  x9[50] = _mm256_subs_epi16(x8[51], x8[50]);
  x9[51] = _mm256_adds_epi16(x8[51], x8[50]);
  x9[52] = _mm256_adds_epi16(x8[52], x8[53]);
  x9[53] = _mm256_subs_epi16(x8[52], x8[53]);
  x9[54] = _mm256_subs_epi16(x8[55], x8[54]);
  x9[55] = _mm256_adds_epi16(x8[55], x8[54]);
  x9[56] = _mm256_adds_epi16(x8[56], x8[57]);
  x9[57] = _mm256_subs_epi16(x8[56], x8[57]);
  x9[58] = _mm256_subs_epi16(x8[59], x8[58]);
  x9[59] = _mm256_adds_epi16(x8[59], x8[58]);
  x9[60] = _mm256_adds_epi16(x8[60], x8[61]);
  x9[61] = _mm256_subs_epi16(x8[60], x8[61]);
  x9[62] = _mm256_subs_epi16(x8[63], x8[62]);
  x9[63] = _mm256_adds_epi16(x8[63], x8[62]);

  // stage 10
  __m256i x10[64];
  x10[0] = x9[0];
  x10[1] = x9[1];
  x10[2] = x9[2];
  x10[3] = x9[3];
  x10[4] = x9[4];
  x10[5] = x9[5];
  x10[6] = x9[6];
  x10[7] = x9[7];
  x10[8] = x9[8];
  x10[9] = x9[9];
  x10[10] = x9[10];
  x10[11] = x9[11];
  x10[12] = x9[12];
  x10[13] = x9[13];
  x10[14] = x9[14];
  x10[15] = x9[15];
  x10[16] = x9[16];
  x10[17] = x9[17];
  x10[18] = x9[18];
  x10[19] = x9[19];
  x10[20] = x9[20];
  x10[21] = x9[21];
  x10[22] = x9[22];
  x10[23] = x9[23];
  x10[24] = x9[24];
  x10[25] = x9[25];
  x10[26] = x9[26];
  x10[27] = x9[27];
  x10[28] = x9[28];
  x10[29] = x9[29];
  x10[30] = x9[30];
  x10[31] = x9[31];
  btf_16_w16_avx2(cospi_p63_p01, cospi_m01_p63, x9[32], x9[63], x10[32],
                  x10[63]);
  btf_16_w16_avx2(cospi_p31_p33, cospi_m33_p31, x9[33], x9[62], x10[33],
                  x10[62]);
  btf_16_w16_avx2(cospi_p47_p17, cospi_m17_p47, x9[34], x9[61], x10[34],
                  x10[61]);
  btf_16_w16_avx2(cospi_p15_p49, cospi_m49_p15, x9[35], x9[60], x10[35],
                  x10[60]);
  btf_16_w16_avx2(cospi_p55_p09, cospi_m09_p55, x9[36], x9[59], x10[36],
                  x10[59]);
  btf_16_w16_avx2(cospi_p23_p41, cospi_m41_p23, x9[37], x9[58], x10[37],
                  x10[58]);
  btf_16_w16_avx2(cospi_p39_p25, cospi_m25_p39, x9[38], x9[57], x10[38],
                  x10[57]);
  btf_16_w16_avx2(cospi_p07_p57, cospi_m57_p07, x9[39], x9[56], x10[39],
                  x10[56]);
  btf_16_w16_avx2(cospi_p59_p05, cospi_m05_p59, x9[40], x9[55], x10[40],
                  x10[55]);
  btf_16_w16_avx2(cospi_p27_p37, cospi_m37_p27, x9[41], x9[54], x10[41],
                  x10[54]);
  btf_16_w16_avx2(cospi_p43_p21, cospi_m21_p43, x9[42], x9[53], x10[42],
                  x10[53]);
  btf_16_w16_avx2(cospi_p11_p53, cospi_m53_p11, x9[43], x9[52], x10[43],
                  x10[52]);
  btf_16_w16_avx2(cospi_p51_p13, cospi_m13_p51, x9[44], x9[51], x10[44],
                  x10[51]);
  btf_16_w16_avx2(cospi_p19_p45, cospi_m45_p19, x9[45], x9[50], x10[45],
                  x10[50]);
  btf_16_w16_avx2(cospi_p35_p29, cospi_m29_p35, x9[46], x9[49], x10[46],
                  x10[49]);
  btf_16_w16_avx2(cospi_p03_p61, cospi_m61_p03, x9[47], x9[48], x10[47],
                  x10[48]);

  // stage 11
  output[0] = x10[0];
  output[1] = x10[32];
  output[2] = x10[16];
  output[3] = x10[48];
  output[4] = x10[8];
  output[5] = x10[40];
  output[6] = x10[24];
  output[7] = x10[56];
  output[8] = x10[4];
  output[9] = x10[36];
  output[10] = x10[20];
  output[11] = x10[52];
  output[12] = x10[12];
  output[13] = x10[44];
  output[14] = x10[28];
  output[15] = x10[60];
  output[16] = x10[2];
  output[17] = x10[34];
  output[18] = x10[18];
  output[19] = x10[50];
  output[20] = x10[10];
  output[21] = x10[42];
  output[22] = x10[26];
  output[23] = x10[58];
  output[24] = x10[6];
  output[25] = x10[38];
  output[26] = x10[22];
  output[27] = x10[54];
  output[28] = x10[14];
  output[29] = x10[46];
  output[30] = x10[30];
  output[31] = x10[62];
  output[32] = x10[1];
  output[33] = x10[33];
  output[34] = x10[17];
  output[35] = x10[49];
  output[36] = x10[9];
  output[37] = x10[41];
  output[38] = x10[25];
  output[39] = x10[57];
  output[40] = x10[5];
  output[41] = x10[37];
  output[42] = x10[21];
  output[43] = x10[53];
  output[44] = x10[13];
  output[45] = x10[45];
  output[46] = x10[29];
  output[47] = x10[61];
  output[48] = x10[3];
  output[49] = x10[35];
  output[50] = x10[19];
  output[51] = x10[51];
  output[52] = x10[11];
  output[53] = x10[43];
  output[54] = x10[27];
  output[55] = x10[59];
  output[56] = x10[7];
  output[57] = x10[39];
  output[58] = x10[23];
  output[59] = x10[55];
  output[60] = x10[15];
  output[61] = x10[47];
  output[62] = x10[31];
  output[63] = x10[63];
}

static INLINE void av1_fdct32_new_avx2(const __m256i *input, __m256i *output,
                                       int8_t cos_bit) {
  __m256i buf0[32];
  __m256i buf1[32];
  const int32_t *cospi;
  // stage 0
  // stage 1
  buf1[0] = _mm256_add_epi32(input[0], input[31]);
  buf1[31] = _mm256_sub_epi32(input[0], input[31]);
  buf1[1] = _mm256_add_epi32(input[1], input[30]);
  buf1[30] = _mm256_sub_epi32(input[1], input[30]);
  buf1[2] = _mm256_add_epi32(input[2], input[29]);
  buf1[29] = _mm256_sub_epi32(input[2], input[29]);
  buf1[3] = _mm256_add_epi32(input[3], input[28]);
  buf1[28] = _mm256_sub_epi32(input[3], input[28]);
  buf1[4] = _mm256_add_epi32(input[4], input[27]);
  buf1[27] = _mm256_sub_epi32(input[4], input[27]);
  buf1[5] = _mm256_add_epi32(input[5], input[26]);
  buf1[26] = _mm256_sub_epi32(input[5], input[26]);
  buf1[6] = _mm256_add_epi32(input[6], input[25]);
  buf1[25] = _mm256_sub_epi32(input[6], input[25]);
  buf1[7] = _mm256_add_epi32(input[7], input[24]);
  buf1[24] = _mm256_sub_epi32(input[7], input[24]);
  buf1[8] = _mm256_add_epi32(input[8], input[23]);
  buf1[23] = _mm256_sub_epi32(input[8], input[23]);
  buf1[9] = _mm256_add_epi32(input[9], input[22]);
  buf1[22] = _mm256_sub_epi32(input[9], input[22]);
  buf1[10] = _mm256_add_epi32(input[10], input[21]);
  buf1[21] = _mm256_sub_epi32(input[10], input[21]);
  buf1[11] = _mm256_add_epi32(input[11], input[20]);
  buf1[20] = _mm256_sub_epi32(input[11], input[20]);
  buf1[12] = _mm256_add_epi32(input[12], input[19]);
  buf1[19] = _mm256_sub_epi32(input[12], input[19]);
  buf1[13] = _mm256_add_epi32(input[13], input[18]);
  buf1[18] = _mm256_sub_epi32(input[13], input[18]);
  buf1[14] = _mm256_add_epi32(input[14], input[17]);
  buf1[17] = _mm256_sub_epi32(input[14], input[17]);
  buf1[15] = _mm256_add_epi32(input[15], input[16]);
  buf1[16] = _mm256_sub_epi32(input[15], input[16]);

  // stage 2
  cospi = cospi_arr(cos_bit);
  buf0[0] = _mm256_add_epi32(buf1[0], buf1[15]);
  buf0[15] = _mm256_sub_epi32(buf1[0], buf1[15]);
  buf0[1] = _mm256_add_epi32(buf1[1], buf1[14]);
  buf0[14] = _mm256_sub_epi32(buf1[1], buf1[14]);
  buf0[2] = _mm256_add_epi32(buf1[2], buf1[13]);
  buf0[13] = _mm256_sub_epi32(buf1[2], buf1[13]);
  buf0[3] = _mm256_add_epi32(buf1[3], buf1[12]);
  buf0[12] = _mm256_sub_epi32(buf1[3], buf1[12]);
  buf0[4] = _mm256_add_epi32(buf1[4], buf1[11]);
  buf0[11] = _mm256_sub_epi32(buf1[4], buf1[11]);
  buf0[5] = _mm256_add_epi32(buf1[5], buf1[10]);
  buf0[10] = _mm256_sub_epi32(buf1[5], buf1[10]);
  buf0[6] = _mm256_add_epi32(buf1[6], buf1[9]);
  buf0[9] = _mm256_sub_epi32(buf1[6], buf1[9]);
  buf0[7] = _mm256_add_epi32(buf1[7], buf1[8]);
  buf0[8] = _mm256_sub_epi32(buf1[7], buf1[8]);
  buf0[16] = buf1[16];
  buf0[17] = buf1[17];
  buf0[18] = buf1[18];
  buf0[19] = buf1[19];
  btf_32_avx2_type0(-cospi[32], cospi[32], buf1[20], buf1[27], buf0[20],
                    buf0[27], cos_bit);
  btf_32_avx2_type0(-cospi[32], cospi[32], buf1[21], buf1[26], buf0[21],
                    buf0[26], cos_bit);
  btf_32_avx2_type0(-cospi[32], cospi[32], buf1[22], buf1[25], buf0[22],
                    buf0[25], cos_bit);
  btf_32_avx2_type0(-cospi[32], cospi[32], buf1[23], buf1[24], buf0[23],
                    buf0[24], cos_bit);
  buf0[28] = buf1[28];
  buf0[29] = buf1[29];
  buf0[30] = buf1[30];
  buf0[31] = buf1[31];

  // stage 3
  cospi = cospi_arr(cos_bit);
  buf1[0] = _mm256_add_epi32(buf0[0], buf0[7]);
  buf1[7] = _mm256_sub_epi32(buf0[0], buf0[7]);
  buf1[1] = _mm256_add_epi32(buf0[1], buf0[6]);
  buf1[6] = _mm256_sub_epi32(buf0[1], buf0[6]);
  buf1[2] = _mm256_add_epi32(buf0[2], buf0[5]);
  buf1[5] = _mm256_sub_epi32(buf0[2], buf0[5]);
  buf1[3] = _mm256_add_epi32(buf0[3], buf0[4]);
  buf1[4] = _mm256_sub_epi32(buf0[3], buf0[4]);
  buf1[8] = buf0[8];
  buf1[9] = buf0[9];
  btf_32_avx2_type0(-cospi[32], cospi[32], buf0[10], buf0[13], buf1[10],
                    buf1[13], cos_bit);
  btf_32_avx2_type0(-cospi[32], cospi[32], buf0[11], buf0[12], buf1[11],
                    buf1[12], cos_bit);
  buf1[14] = buf0[14];
  buf1[15] = buf0[15];
  buf1[16] = _mm256_add_epi32(buf0[16], buf0[23]);
  buf1[23] = _mm256_sub_epi32(buf0[16], buf0[23]);
  buf1[17] = _mm256_add_epi32(buf0[17], buf0[22]);
  buf1[22] = _mm256_sub_epi32(buf0[17], buf0[22]);
  buf1[18] = _mm256_add_epi32(buf0[18], buf0[21]);
  buf1[21] = _mm256_sub_epi32(buf0[18], buf0[21]);
  buf1[19] = _mm256_add_epi32(buf0[19], buf0[20]);
  buf1[20] = _mm256_sub_epi32(buf0[19], buf0[20]);
  buf1[24] = _mm256_sub_epi32(buf0[31], buf0[24]);
  buf1[31] = _mm256_add_epi32(buf0[31], buf0[24]);
  buf1[25] = _mm256_sub_epi32(buf0[30], buf0[25]);
  buf1[30] = _mm256_add_epi32(buf0[30], buf0[25]);
  buf1[26] = _mm256_sub_epi32(buf0[29], buf0[26]);
  buf1[29] = _mm256_add_epi32(buf0[29], buf0[26]);
  buf1[27] = _mm256_sub_epi32(buf0[28], buf0[27]);
  buf1[28] = _mm256_add_epi32(buf0[28], buf0[27]);

  // stage 4
  cospi = cospi_arr(cos_bit);
  buf0[0] = _mm256_add_epi32(buf1[0], buf1[3]);
  buf0[3] = _mm256_sub_epi32(buf1[0], buf1[3]);
  buf0[1] = _mm256_add_epi32(buf1[1], buf1[2]);
  buf0[2] = _mm256_sub_epi32(buf1[1], buf1[2]);
  buf0[4] = buf1[4];
  btf_32_avx2_type0(-cospi[32], cospi[32], buf1[5], buf1[6], buf0[5], buf0[6],
                    cos_bit);
  buf0[7] = buf1[7];
  buf0[8] = _mm256_add_epi32(buf1[8], buf1[11]);
  buf0[11] = _mm256_sub_epi32(buf1[8], buf1[11]);
  buf0[9] = _mm256_add_epi32(buf1[9], buf1[10]);
  buf0[10] = _mm256_sub_epi32(buf1[9], buf1[10]);
  buf0[12] = _mm256_sub_epi32(buf1[15], buf1[12]);
  buf0[15] = _mm256_add_epi32(buf1[15], buf1[12]);
  buf0[13] = _mm256_sub_epi32(buf1[14], buf1[13]);
  buf0[14] = _mm256_add_epi32(buf1[14], buf1[13]);
  buf0[16] = buf1[16];
  buf0[17] = buf1[17];
  btf_32_avx2_type0(-cospi[16], cospi[48], buf1[18], buf1[29], buf0[18],
                    buf0[29], cos_bit);
  btf_32_avx2_type0(-cospi[16], cospi[48], buf1[19], buf1[28], buf0[19],
                    buf0[28], cos_bit);
  btf_32_avx2_type0(-cospi[48], -cospi[16], buf1[20], buf1[27], buf0[20],
                    buf0[27], cos_bit);
  btf_32_avx2_type0(-cospi[48], -cospi[16], buf1[21], buf1[26], buf0[21],
                    buf0[26], cos_bit);
  buf0[22] = buf1[22];
  buf0[23] = buf1[23];
  buf0[24] = buf1[24];
  buf0[25] = buf1[25];
  buf0[30] = buf1[30];
  buf0[31] = buf1[31];

  // stage 5
  cospi = cospi_arr(cos_bit);
  btf_32_avx2_type0(cospi[32], cospi[32], buf0[0], buf0[1], buf1[0], buf1[1],
                    cos_bit);
  btf_32_avx2_type1(cospi[48], cospi[16], buf0[2], buf0[3], buf1[2], buf1[3],
                    cos_bit);
  buf1[4] = _mm256_add_epi32(buf0[4], buf0[5]);
  buf1[5] = _mm256_sub_epi32(buf0[4], buf0[5]);
  buf1[6] = _mm256_sub_epi32(buf0[7], buf0[6]);
  buf1[7] = _mm256_add_epi32(buf0[7], buf0[6]);
  buf1[8] = buf0[8];
  btf_32_avx2_type0(-cospi[16], cospi[48], buf0[9], buf0[14], buf1[9], buf1[14],
                    cos_bit);
  btf_32_avx2_type0(-cospi[48], -cospi[16], buf0[10], buf0[13], buf1[10],
                    buf1[13], cos_bit);
  buf1[11] = buf0[11];
  buf1[12] = buf0[12];
  buf1[15] = buf0[15];
  buf1[16] = _mm256_add_epi32(buf0[16], buf0[19]);
  buf1[19] = _mm256_sub_epi32(buf0[16], buf0[19]);
  buf1[17] = _mm256_add_epi32(buf0[17], buf0[18]);
  buf1[18] = _mm256_sub_epi32(buf0[17], buf0[18]);
  buf1[20] = _mm256_sub_epi32(buf0[23], buf0[20]);
  buf1[23] = _mm256_add_epi32(buf0[23], buf0[20]);
  buf1[21] = _mm256_sub_epi32(buf0[22], buf0[21]);
  buf1[22] = _mm256_add_epi32(buf0[22], buf0[21]);
  buf1[24] = _mm256_add_epi32(buf0[24], buf0[27]);
  buf1[27] = _mm256_sub_epi32(buf0[24], buf0[27]);
  buf1[25] = _mm256_add_epi32(buf0[25], buf0[26]);
  buf1[26] = _mm256_sub_epi32(buf0[25], buf0[26]);
  buf1[28] = _mm256_sub_epi32(buf0[31], buf0[28]);
  buf1[31] = _mm256_add_epi32(buf0[31], buf0[28]);
  buf1[29] = _mm256_sub_epi32(buf0[30], buf0[29]);
  buf1[30] = _mm256_add_epi32(buf0[30], buf0[29]);

  // stage 6
  cospi = cospi_arr(cos_bit);
  buf0[0] = buf1[0];
  buf0[1] = buf1[1];
  buf0[2] = buf1[2];
  buf0[3] = buf1[3];
  btf_32_avx2_type1(cospi[56], cospi[8], buf1[4], buf1[7], buf0[4], buf0[7],
                    cos_bit);
  btf_32_avx2_type1(cospi[24], cospi[40], buf1[5], buf1[6], buf0[5], buf0[6],
                    cos_bit);
  buf0[8] = _mm256_add_epi32(buf1[8], buf1[9]);
  buf0[9] = _mm256_sub_epi32(buf1[8], buf1[9]);
  buf0[10] = _mm256_sub_epi32(buf1[11], buf1[10]);
  buf0[11] = _mm256_add_epi32(buf1[11], buf1[10]);
  buf0[12] = _mm256_add_epi32(buf1[12], buf1[13]);
  buf0[13] = _mm256_sub_epi32(buf1[12], buf1[13]);
  buf0[14] = _mm256_sub_epi32(buf1[15], buf1[14]);
  buf0[15] = _mm256_add_epi32(buf1[15], buf1[14]);
  buf0[16] = buf1[16];
  btf_32_avx2_type0(-cospi[8], cospi[56], buf1[17], buf1[30], buf0[17],
                    buf0[30], cos_bit);
  btf_32_avx2_type0(-cospi[56], -cospi[8], buf1[18], buf1[29], buf0[18],
                    buf0[29], cos_bit);
  buf0[19] = buf1[19];
  buf0[20] = buf1[20];
  btf_32_avx2_type0(-cospi[40], cospi[24], buf1[21], buf1[26], buf0[21],
                    buf0[26], cos_bit);
  btf_32_avx2_type0(-cospi[24], -cospi[40], buf1[22], buf1[25], buf0[22],
                    buf0[25], cos_bit);
  buf0[23] = buf1[23];
  buf0[24] = buf1[24];
  buf0[27] = buf1[27];
  buf0[28] = buf1[28];
  buf0[31] = buf1[31];

  // stage 7
  cospi = cospi_arr(cos_bit);
  buf1[0] = buf0[0];
  buf1[1] = buf0[1];
  buf1[2] = buf0[2];
  buf1[3] = buf0[3];
  buf1[4] = buf0[4];
  buf1[5] = buf0[5];
  buf1[6] = buf0[6];
  buf1[7] = buf0[7];
  btf_32_avx2_type1(cospi[60], cospi[4], buf0[8], buf0[15], buf1[8], buf1[15],
                    cos_bit);
  btf_32_avx2_type1(cospi[28], cospi[36], buf0[9], buf0[14], buf1[9], buf1[14],
                    cos_bit);
  btf_32_avx2_type1(cospi[44], cospi[20], buf0[10], buf0[13], buf1[10],
                    buf1[13], cos_bit);
  btf_32_avx2_type1(cospi[12], cospi[52], buf0[11], buf0[12], buf1[11],
                    buf1[12], cos_bit);
  buf1[16] = _mm256_add_epi32(buf0[16], buf0[17]);
  buf1[17] = _mm256_sub_epi32(buf0[16], buf0[17]);
  buf1[18] = _mm256_sub_epi32(buf0[19], buf0[18]);
  buf1[19] = _mm256_add_epi32(buf0[19], buf0[18]);
  buf1[20] = _mm256_add_epi32(buf0[20], buf0[21]);
  buf1[21] = _mm256_sub_epi32(buf0[20], buf0[21]);
  buf1[22] = _mm256_sub_epi32(buf0[23], buf0[22]);
  buf1[23] = _mm256_add_epi32(buf0[23], buf0[22]);
  buf1[24] = _mm256_add_epi32(buf0[24], buf0[25]);
  buf1[25] = _mm256_sub_epi32(buf0[24], buf0[25]);
  buf1[26] = _mm256_sub_epi32(buf0[27], buf0[26]);
  buf1[27] = _mm256_add_epi32(buf0[27], buf0[26]);
  buf1[28] = _mm256_add_epi32(buf0[28], buf0[29]);
  buf1[29] = _mm256_sub_epi32(buf0[28], buf0[29]);
  buf1[30] = _mm256_sub_epi32(buf0[31], buf0[30]);
  buf1[31] = _mm256_add_epi32(buf0[31], buf0[30]);

  // stage 8
  cospi = cospi_arr(cos_bit);
  buf0[0] = buf1[0];
  buf0[1] = buf1[1];
  buf0[2] = buf1[2];
  buf0[3] = buf1[3];
  buf0[4] = buf1[4];
  buf0[5] = buf1[5];
  buf0[6] = buf1[6];
  buf0[7] = buf1[7];
  buf0[8] = buf1[8];
  buf0[9] = buf1[9];
  buf0[10] = buf1[10];
  buf0[11] = buf1[11];
  buf0[12] = buf1[12];
  buf0[13] = buf1[13];
  buf0[14] = buf1[14];
  buf0[15] = buf1[15];
  btf_32_avx2_type1(cospi[62], cospi[2], buf1[16], buf1[31], buf0[16], buf0[31],
                    cos_bit);
  btf_32_avx2_type1(cospi[30], cospi[34], buf1[17], buf1[30], buf0[17],
                    buf0[30], cos_bit);
  btf_32_avx2_type1(cospi[46], cospi[18], buf1[18], buf1[29], buf0[18],
                    buf0[29], cos_bit);
  btf_32_avx2_type1(cospi[14], cospi[50], buf1[19], buf1[28], buf0[19],
                    buf0[28], cos_bit);
  btf_32_avx2_type1(cospi[54], cospi[10], buf1[20], buf1[27], buf0[20],
                    buf0[27], cos_bit);
  btf_32_avx2_type1(cospi[22], cospi[42], buf1[21], buf1[26], buf0[21],
                    buf0[26], cos_bit);
  btf_32_avx2_type1(cospi[38], cospi[26], buf1[22], buf1[25], buf0[22],
                    buf0[25], cos_bit);
  btf_32_avx2_type1(cospi[6], cospi[58], buf1[23], buf1[24], buf0[23], buf0[24],
                    cos_bit);

  // stage 9
  output[0] = buf0[0];
  output[1] = buf0[16];
  output[2] = buf0[8];
  output[3] = buf0[24];
  output[4] = buf0[4];
  output[5] = buf0[20];
  output[6] = buf0[12];
  output[7] = buf0[28];
  output[8] = buf0[2];
  output[9] = buf0[18];
  output[10] = buf0[10];
  output[11] = buf0[26];
  output[12] = buf0[6];
  output[13] = buf0[22];
  output[14] = buf0[14];
  output[15] = buf0[30];
  output[16] = buf0[1];
  output[17] = buf0[17];
  output[18] = buf0[9];
  output[19] = buf0[25];
  output[20] = buf0[5];
  output[21] = buf0[21];
  output[22] = buf0[13];
  output[23] = buf0[29];
  output[24] = buf0[3];
  output[25] = buf0[19];
  output[26] = buf0[11];
  output[27] = buf0[27];
  output[28] = buf0[7];
  output[29] = buf0[23];
  output[30] = buf0[15];
  output[31] = buf0[31];
}

static INLINE void av1_fdct64_new_avx2(const __m256i *input, __m256i *output,
                                       int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m256i __rounding = _mm256_set1_epi32(1 << (cos_bit - 1));

  __m256i cospi_m32 = _mm256_set1_epi32(-cospi[32]);
  __m256i cospi_p32 = _mm256_set1_epi32(cospi[32]);
  __m256i cospi_m16 = _mm256_set1_epi32(-cospi[16]);
  __m256i cospi_p48 = _mm256_set1_epi32(cospi[48]);
  __m256i cospi_m48 = _mm256_set1_epi32(-cospi[48]);
  __m256i cospi_p16 = _mm256_set1_epi32(cospi[16]);
  __m256i cospi_m08 = _mm256_set1_epi32(-cospi[8]);
  __m256i cospi_p56 = _mm256_set1_epi32(cospi[56]);
  __m256i cospi_m56 = _mm256_set1_epi32(-cospi[56]);
  __m256i cospi_m40 = _mm256_set1_epi32(-cospi[40]);
  __m256i cospi_p24 = _mm256_set1_epi32(cospi[24]);
  __m256i cospi_m24 = _mm256_set1_epi32(-cospi[24]);
  __m256i cospi_p08 = _mm256_set1_epi32(cospi[8]);
  __m256i cospi_p40 = _mm256_set1_epi32(cospi[40]);
  __m256i cospi_p60 = _mm256_set1_epi32(cospi[60]);
  __m256i cospi_p04 = _mm256_set1_epi32(cospi[4]);
  __m256i cospi_p28 = _mm256_set1_epi32(cospi[28]);
  __m256i cospi_p36 = _mm256_set1_epi32(cospi[36]);
  __m256i cospi_p44 = _mm256_set1_epi32(cospi[44]);
  __m256i cospi_p20 = _mm256_set1_epi32(cospi[20]);
  __m256i cospi_p12 = _mm256_set1_epi32(cospi[12]);
  __m256i cospi_p52 = _mm256_set1_epi32(cospi[52]);
  __m256i cospi_m04 = _mm256_set1_epi32(-cospi[4]);
  __m256i cospi_m60 = _mm256_set1_epi32(-cospi[60]);
  __m256i cospi_m36 = _mm256_set1_epi32(-cospi[36]);
  __m256i cospi_m28 = _mm256_set1_epi32(-cospi[28]);
  __m256i cospi_m20 = _mm256_set1_epi32(-cospi[20]);
  __m256i cospi_m44 = _mm256_set1_epi32(-cospi[44]);
  __m256i cospi_m52 = _mm256_set1_epi32(-cospi[52]);
  __m256i cospi_m12 = _mm256_set1_epi32(-cospi[12]);
  __m256i cospi_p62 = _mm256_set1_epi32(cospi[62]);
  __m256i cospi_p02 = _mm256_set1_epi32(cospi[2]);
  __m256i cospi_p30 = _mm256_set1_epi32(cospi[30]);
  __m256i cospi_p34 = _mm256_set1_epi32(cospi[34]);
  __m256i cospi_p46 = _mm256_set1_epi32(cospi[46]);
  __m256i cospi_p18 = _mm256_set1_epi32(cospi[18]);
  __m256i cospi_p14 = _mm256_set1_epi32(cospi[14]);
  __m256i cospi_p50 = _mm256_set1_epi32(cospi[50]);
  __m256i cospi_p54 = _mm256_set1_epi32(cospi[54]);
  __m256i cospi_p10 = _mm256_set1_epi32(cospi[10]);
  __m256i cospi_p22 = _mm256_set1_epi32(cospi[22]);
  __m256i cospi_p42 = _mm256_set1_epi32(cospi[42]);
  __m256i cospi_p38 = _mm256_set1_epi32(cospi[38]);
  __m256i cospi_p26 = _mm256_set1_epi32(cospi[26]);
  __m256i cospi_p06 = _mm256_set1_epi32(cospi[6]);
  __m256i cospi_p58 = _mm256_set1_epi32(cospi[58]);
  __m256i cospi_p63 = _mm256_set1_epi32(cospi[63]);
  __m256i cospi_p01 = _mm256_set1_epi32(cospi[1]);
  __m256i cospi_p31 = _mm256_set1_epi32(cospi[31]);
  __m256i cospi_p33 = _mm256_set1_epi32(cospi[33]);
  __m256i cospi_p47 = _mm256_set1_epi32(cospi[47]);
  __m256i cospi_p17 = _mm256_set1_epi32(cospi[17]);
  __m256i cospi_p15 = _mm256_set1_epi32(cospi[15]);
  __m256i cospi_p49 = _mm256_set1_epi32(cospi[49]);
  __m256i cospi_p55 = _mm256_set1_epi32(cospi[55]);
  __m256i cospi_p09 = _mm256_set1_epi32(cospi[9]);
  __m256i cospi_p23 = _mm256_set1_epi32(cospi[23]);
  __m256i cospi_p41 = _mm256_set1_epi32(cospi[41]);
  __m256i cospi_p39 = _mm256_set1_epi32(cospi[39]);
  __m256i cospi_p25 = _mm256_set1_epi32(cospi[25]);
  __m256i cospi_p07 = _mm256_set1_epi32(cospi[7]);
  __m256i cospi_p57 = _mm256_set1_epi32(cospi[57]);
  __m256i cospi_p59 = _mm256_set1_epi32(cospi[59]);
  __m256i cospi_p05 = _mm256_set1_epi32(cospi[5]);
  __m256i cospi_p27 = _mm256_set1_epi32(cospi[27]);
  __m256i cospi_p37 = _mm256_set1_epi32(cospi[37]);
  __m256i cospi_p43 = _mm256_set1_epi32(cospi[43]);
  __m256i cospi_p21 = _mm256_set1_epi32(cospi[21]);
  __m256i cospi_p11 = _mm256_set1_epi32(cospi[11]);
  __m256i cospi_p53 = _mm256_set1_epi32(cospi[53]);
  __m256i cospi_p51 = _mm256_set1_epi32(cospi[51]);
  __m256i cospi_p13 = _mm256_set1_epi32(cospi[13]);
  __m256i cospi_p19 = _mm256_set1_epi32(cospi[19]);
  __m256i cospi_p45 = _mm256_set1_epi32(cospi[45]);
  __m256i cospi_p35 = _mm256_set1_epi32(cospi[35]);
  __m256i cospi_p29 = _mm256_set1_epi32(cospi[29]);
  __m256i cospi_p03 = _mm256_set1_epi32(cospi[3]);
  __m256i cospi_p61 = _mm256_set1_epi32(cospi[61]);

  // stage 1
  __m256i x1[64];
  x1[0] = _mm256_add_epi32(input[0], input[63]);
  x1[63] = _mm256_sub_epi32(input[0], input[63]);
  x1[1] = _mm256_add_epi32(input[1], input[62]);
  x1[62] = _mm256_sub_epi32(input[1], input[62]);
  x1[2] = _mm256_add_epi32(input[2], input[61]);
  x1[61] = _mm256_sub_epi32(input[2], input[61]);
  x1[3] = _mm256_add_epi32(input[3], input[60]);
  x1[60] = _mm256_sub_epi32(input[3], input[60]);
  x1[4] = _mm256_add_epi32(input[4], input[59]);
  x1[59] = _mm256_sub_epi32(input[4], input[59]);
  x1[5] = _mm256_add_epi32(input[5], input[58]);
  x1[58] = _mm256_sub_epi32(input[5], input[58]);
  x1[6] = _mm256_add_epi32(input[6], input[57]);
  x1[57] = _mm256_sub_epi32(input[6], input[57]);
  x1[7] = _mm256_add_epi32(input[7], input[56]);
  x1[56] = _mm256_sub_epi32(input[7], input[56]);
  x1[8] = _mm256_add_epi32(input[8], input[55]);
  x1[55] = _mm256_sub_epi32(input[8], input[55]);
  x1[9] = _mm256_add_epi32(input[9], input[54]);
  x1[54] = _mm256_sub_epi32(input[9], input[54]);
  x1[10] = _mm256_add_epi32(input[10], input[53]);
  x1[53] = _mm256_sub_epi32(input[10], input[53]);
  x1[11] = _mm256_add_epi32(input[11], input[52]);
  x1[52] = _mm256_sub_epi32(input[11], input[52]);
  x1[12] = _mm256_add_epi32(input[12], input[51]);
  x1[51] = _mm256_sub_epi32(input[12], input[51]);
  x1[13] = _mm256_add_epi32(input[13], input[50]);
  x1[50] = _mm256_sub_epi32(input[13], input[50]);
  x1[14] = _mm256_add_epi32(input[14], input[49]);
  x1[49] = _mm256_sub_epi32(input[14], input[49]);
  x1[15] = _mm256_add_epi32(input[15], input[48]);
  x1[48] = _mm256_sub_epi32(input[15], input[48]);
  x1[16] = _mm256_add_epi32(input[16], input[47]);
  x1[47] = _mm256_sub_epi32(input[16], input[47]);
  x1[17] = _mm256_add_epi32(input[17], input[46]);
  x1[46] = _mm256_sub_epi32(input[17], input[46]);
  x1[18] = _mm256_add_epi32(input[18], input[45]);
  x1[45] = _mm256_sub_epi32(input[18], input[45]);
  x1[19] = _mm256_add_epi32(input[19], input[44]);
  x1[44] = _mm256_sub_epi32(input[19], input[44]);
  x1[20] = _mm256_add_epi32(input[20], input[43]);
  x1[43] = _mm256_sub_epi32(input[20], input[43]);
  x1[21] = _mm256_add_epi32(input[21], input[42]);
  x1[42] = _mm256_sub_epi32(input[21], input[42]);
  x1[22] = _mm256_add_epi32(input[22], input[41]);
  x1[41] = _mm256_sub_epi32(input[22], input[41]);
  x1[23] = _mm256_add_epi32(input[23], input[40]);
  x1[40] = _mm256_sub_epi32(input[23], input[40]);
  x1[24] = _mm256_add_epi32(input[24], input[39]);
  x1[39] = _mm256_sub_epi32(input[24], input[39]);
  x1[25] = _mm256_add_epi32(input[25], input[38]);
  x1[38] = _mm256_sub_epi32(input[25], input[38]);
  x1[26] = _mm256_add_epi32(input[26], input[37]);
  x1[37] = _mm256_sub_epi32(input[26], input[37]);
  x1[27] = _mm256_add_epi32(input[27], input[36]);
  x1[36] = _mm256_sub_epi32(input[27], input[36]);
  x1[28] = _mm256_add_epi32(input[28], input[35]);
  x1[35] = _mm256_sub_epi32(input[28], input[35]);
  x1[29] = _mm256_add_epi32(input[29], input[34]);
  x1[34] = _mm256_sub_epi32(input[29], input[34]);
  x1[30] = _mm256_add_epi32(input[30], input[33]);
  x1[33] = _mm256_sub_epi32(input[30], input[33]);
  x1[31] = _mm256_add_epi32(input[31], input[32]);
  x1[32] = _mm256_sub_epi32(input[31], input[32]);

  // stage 2
  __m256i x2[64];
  x2[0] = _mm256_add_epi32(x1[0], x1[31]);
  x2[31] = _mm256_sub_epi32(x1[0], x1[31]);
  x2[1] = _mm256_add_epi32(x1[1], x1[30]);
  x2[30] = _mm256_sub_epi32(x1[1], x1[30]);
  x2[2] = _mm256_add_epi32(x1[2], x1[29]);
  x2[29] = _mm256_sub_epi32(x1[2], x1[29]);
  x2[3] = _mm256_add_epi32(x1[3], x1[28]);
  x2[28] = _mm256_sub_epi32(x1[3], x1[28]);
  x2[4] = _mm256_add_epi32(x1[4], x1[27]);
  x2[27] = _mm256_sub_epi32(x1[4], x1[27]);
  x2[5] = _mm256_add_epi32(x1[5], x1[26]);
  x2[26] = _mm256_sub_epi32(x1[5], x1[26]);
  x2[6] = _mm256_add_epi32(x1[6], x1[25]);
  x2[25] = _mm256_sub_epi32(x1[6], x1[25]);
  x2[7] = _mm256_add_epi32(x1[7], x1[24]);
  x2[24] = _mm256_sub_epi32(x1[7], x1[24]);
  x2[8] = _mm256_add_epi32(x1[8], x1[23]);
  x2[23] = _mm256_sub_epi32(x1[8], x1[23]);
  x2[9] = _mm256_add_epi32(x1[9], x1[22]);
  x2[22] = _mm256_sub_epi32(x1[9], x1[22]);
  x2[10] = _mm256_add_epi32(x1[10], x1[21]);
  x2[21] = _mm256_sub_epi32(x1[10], x1[21]);
  x2[11] = _mm256_add_epi32(x1[11], x1[20]);
  x2[20] = _mm256_sub_epi32(x1[11], x1[20]);
  x2[12] = _mm256_add_epi32(x1[12], x1[19]);
  x2[19] = _mm256_sub_epi32(x1[12], x1[19]);
  x2[13] = _mm256_add_epi32(x1[13], x1[18]);
  x2[18] = _mm256_sub_epi32(x1[13], x1[18]);
  x2[14] = _mm256_add_epi32(x1[14], x1[17]);
  x2[17] = _mm256_sub_epi32(x1[14], x1[17]);
  x2[15] = _mm256_add_epi32(x1[15], x1[16]);
  x2[16] = _mm256_sub_epi32(x1[15], x1[16]);
  x2[32] = x1[32];
  x2[33] = x1[33];
  x2[34] = x1[34];
  x2[35] = x1[35];
  x2[36] = x1[36];
  x2[37] = x1[37];
  x2[38] = x1[38];
  x2[39] = x1[39];
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[40], x1[55], x2[40], x2[55],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[41], x1[54], x2[41], x2[54],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[42], x1[53], x2[42], x2[53],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[43], x1[52], x2[43], x2[52],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[44], x1[51], x2[44], x2[51],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[45], x1[50], x2[45], x2[50],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[46], x1[49], x2[46], x2[49],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x1[47], x1[48], x2[47], x2[48],
                        __rounding, cos_bit);
  x2[56] = x1[56];
  x2[57] = x1[57];
  x2[58] = x1[58];
  x2[59] = x1[59];
  x2[60] = x1[60];
  x2[61] = x1[61];
  x2[62] = x1[62];
  x2[63] = x1[63];

  // stage 3
  __m256i x3[64];
  x3[0] = _mm256_add_epi32(x2[0], x2[15]);
  x3[15] = _mm256_sub_epi32(x2[0], x2[15]);
  x3[1] = _mm256_add_epi32(x2[1], x2[14]);
  x3[14] = _mm256_sub_epi32(x2[1], x2[14]);
  x3[2] = _mm256_add_epi32(x2[2], x2[13]);
  x3[13] = _mm256_sub_epi32(x2[2], x2[13]);
  x3[3] = _mm256_add_epi32(x2[3], x2[12]);
  x3[12] = _mm256_sub_epi32(x2[3], x2[12]);
  x3[4] = _mm256_add_epi32(x2[4], x2[11]);
  x3[11] = _mm256_sub_epi32(x2[4], x2[11]);
  x3[5] = _mm256_add_epi32(x2[5], x2[10]);
  x3[10] = _mm256_sub_epi32(x2[5], x2[10]);
  x3[6] = _mm256_add_epi32(x2[6], x2[9]);
  x3[9] = _mm256_sub_epi32(x2[6], x2[9]);
  x3[7] = _mm256_add_epi32(x2[7], x2[8]);
  x3[8] = _mm256_sub_epi32(x2[7], x2[8]);
  x3[16] = x2[16];
  x3[17] = x2[17];
  x3[18] = x2[18];
  x3[19] = x2[19];
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x2[20], x2[27], x3[20], x3[27],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x2[21], x2[26], x3[21], x3[26],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x2[22], x2[25], x3[22], x3[25],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x2[23], x2[24], x3[23], x3[24],
                        __rounding, cos_bit);
  x3[28] = x2[28];
  x3[29] = x2[29];
  x3[30] = x2[30];
  x3[31] = x2[31];
  x3[32] = _mm256_add_epi32(x2[32], x2[47]);
  x3[47] = _mm256_sub_epi32(x2[32], x2[47]);
  x3[33] = _mm256_add_epi32(x2[33], x2[46]);
  x3[46] = _mm256_sub_epi32(x2[33], x2[46]);
  x3[34] = _mm256_add_epi32(x2[34], x2[45]);
  x3[45] = _mm256_sub_epi32(x2[34], x2[45]);
  x3[35] = _mm256_add_epi32(x2[35], x2[44]);
  x3[44] = _mm256_sub_epi32(x2[35], x2[44]);
  x3[36] = _mm256_add_epi32(x2[36], x2[43]);
  x3[43] = _mm256_sub_epi32(x2[36], x2[43]);
  x3[37] = _mm256_add_epi32(x2[37], x2[42]);
  x3[42] = _mm256_sub_epi32(x2[37], x2[42]);
  x3[38] = _mm256_add_epi32(x2[38], x2[41]);
  x3[41] = _mm256_sub_epi32(x2[38], x2[41]);
  x3[39] = _mm256_add_epi32(x2[39], x2[40]);
  x3[40] = _mm256_sub_epi32(x2[39], x2[40]);
  x3[48] = _mm256_sub_epi32(x2[63], x2[48]);
  x3[63] = _mm256_add_epi32(x2[63], x2[48]);
  x3[49] = _mm256_sub_epi32(x2[62], x2[49]);
  x3[62] = _mm256_add_epi32(x2[62], x2[49]);
  x3[50] = _mm256_sub_epi32(x2[61], x2[50]);
  x3[61] = _mm256_add_epi32(x2[61], x2[50]);
  x3[51] = _mm256_sub_epi32(x2[60], x2[51]);
  x3[60] = _mm256_add_epi32(x2[60], x2[51]);
  x3[52] = _mm256_sub_epi32(x2[59], x2[52]);
  x3[59] = _mm256_add_epi32(x2[59], x2[52]);
  x3[53] = _mm256_sub_epi32(x2[58], x2[53]);
  x3[58] = _mm256_add_epi32(x2[58], x2[53]);
  x3[54] = _mm256_sub_epi32(x2[57], x2[54]);
  x3[57] = _mm256_add_epi32(x2[57], x2[54]);
  x3[55] = _mm256_sub_epi32(x2[56], x2[55]);
  x3[56] = _mm256_add_epi32(x2[56], x2[55]);

  // stage 4
  __m256i x4[64];
  x4[0] = _mm256_add_epi32(x3[0], x3[7]);
  x4[7] = _mm256_sub_epi32(x3[0], x3[7]);
  x4[1] = _mm256_add_epi32(x3[1], x3[6]);
  x4[6] = _mm256_sub_epi32(x3[1], x3[6]);
  x4[2] = _mm256_add_epi32(x3[2], x3[5]);
  x4[5] = _mm256_sub_epi32(x3[2], x3[5]);
  x4[3] = _mm256_add_epi32(x3[3], x3[4]);
  x4[4] = _mm256_sub_epi32(x3[3], x3[4]);
  x4[8] = x3[8];
  x4[9] = x3[9];
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x3[10], x3[13], x4[10], x4[13],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x3[11], x3[12], x4[11], x4[12],
                        __rounding, cos_bit);
  x4[14] = x3[14];
  x4[15] = x3[15];
  x4[16] = _mm256_add_epi32(x3[16], x3[23]);
  x4[23] = _mm256_sub_epi32(x3[16], x3[23]);
  x4[17] = _mm256_add_epi32(x3[17], x3[22]);
  x4[22] = _mm256_sub_epi32(x3[17], x3[22]);
  x4[18] = _mm256_add_epi32(x3[18], x3[21]);
  x4[21] = _mm256_sub_epi32(x3[18], x3[21]);
  x4[19] = _mm256_add_epi32(x3[19], x3[20]);
  x4[20] = _mm256_sub_epi32(x3[19], x3[20]);
  x4[24] = _mm256_sub_epi32(x3[31], x3[24]);
  x4[31] = _mm256_add_epi32(x3[31], x3[24]);
  x4[25] = _mm256_sub_epi32(x3[30], x3[25]);
  x4[30] = _mm256_add_epi32(x3[30], x3[25]);
  x4[26] = _mm256_sub_epi32(x3[29], x3[26]);
  x4[29] = _mm256_add_epi32(x3[29], x3[26]);
  x4[27] = _mm256_sub_epi32(x3[28], x3[27]);
  x4[28] = _mm256_add_epi32(x3[28], x3[27]);
  x4[32] = x3[32];
  x4[33] = x3[33];
  x4[34] = x3[34];
  x4[35] = x3[35];
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x3[36], x3[59], x4[36], x4[59],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x3[37], x3[58], x4[37], x4[58],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x3[38], x3[57], x4[38], x4[57],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x3[39], x3[56], x4[39], x4[56],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x3[40], x3[55], x4[40], x4[55],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x3[41], x3[54], x4[41], x4[54],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x3[42], x3[53], x4[42], x4[53],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x3[43], x3[52], x4[43], x4[52],
                        __rounding, cos_bit);
  x4[44] = x3[44];
  x4[45] = x3[45];
  x4[46] = x3[46];
  x4[47] = x3[47];
  x4[48] = x3[48];
  x4[49] = x3[49];
  x4[50] = x3[50];
  x4[51] = x3[51];
  x4[60] = x3[60];
  x4[61] = x3[61];
  x4[62] = x3[62];
  x4[63] = x3[63];

  // stage 5
  __m256i x5[64];
  x5[0] = _mm256_add_epi32(x4[0], x4[3]);
  x5[3] = _mm256_sub_epi32(x4[0], x4[3]);
  x5[1] = _mm256_add_epi32(x4[1], x4[2]);
  x5[2] = _mm256_sub_epi32(x4[1], x4[2]);
  x5[4] = x4[4];
  btf_32_type0_avx2_new(cospi_m32, cospi_p32, x4[5], x4[6], x5[5], x5[6],
                        __rounding, cos_bit);
  x5[7] = x4[7];
  x5[8] = _mm256_add_epi32(x4[8], x4[11]);
  x5[11] = _mm256_sub_epi32(x4[8], x4[11]);
  x5[9] = _mm256_add_epi32(x4[9], x4[10]);
  x5[10] = _mm256_sub_epi32(x4[9], x4[10]);
  x5[12] = _mm256_sub_epi32(x4[15], x4[12]);
  x5[15] = _mm256_add_epi32(x4[15], x4[12]);
  x5[13] = _mm256_sub_epi32(x4[14], x4[13]);
  x5[14] = _mm256_add_epi32(x4[14], x4[13]);
  x5[16] = x4[16];
  x5[17] = x4[17];
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x4[18], x4[29], x5[18], x5[29],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x4[19], x4[28], x5[19], x5[28],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x4[20], x4[27], x5[20], x5[27],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x4[21], x4[26], x5[21], x5[26],
                        __rounding, cos_bit);
  x5[22] = x4[22];
  x5[23] = x4[23];
  x5[24] = x4[24];
  x5[25] = x4[25];
  x5[30] = x4[30];
  x5[31] = x4[31];
  x5[32] = _mm256_add_epi32(x4[32], x4[39]);
  x5[39] = _mm256_sub_epi32(x4[32], x4[39]);
  x5[33] = _mm256_add_epi32(x4[33], x4[38]);
  x5[38] = _mm256_sub_epi32(x4[33], x4[38]);
  x5[34] = _mm256_add_epi32(x4[34], x4[37]);
  x5[37] = _mm256_sub_epi32(x4[34], x4[37]);
  x5[35] = _mm256_add_epi32(x4[35], x4[36]);
  x5[36] = _mm256_sub_epi32(x4[35], x4[36]);
  x5[40] = _mm256_sub_epi32(x4[47], x4[40]);
  x5[47] = _mm256_add_epi32(x4[47], x4[40]);
  x5[41] = _mm256_sub_epi32(x4[46], x4[41]);
  x5[46] = _mm256_add_epi32(x4[46], x4[41]);
  x5[42] = _mm256_sub_epi32(x4[45], x4[42]);
  x5[45] = _mm256_add_epi32(x4[45], x4[42]);
  x5[43] = _mm256_sub_epi32(x4[44], x4[43]);
  x5[44] = _mm256_add_epi32(x4[44], x4[43]);
  x5[48] = _mm256_add_epi32(x4[48], x4[55]);
  x5[55] = _mm256_sub_epi32(x4[48], x4[55]);
  x5[49] = _mm256_add_epi32(x4[49], x4[54]);
  x5[54] = _mm256_sub_epi32(x4[49], x4[54]);
  x5[50] = _mm256_add_epi32(x4[50], x4[53]);
  x5[53] = _mm256_sub_epi32(x4[50], x4[53]);
  x5[51] = _mm256_add_epi32(x4[51], x4[52]);
  x5[52] = _mm256_sub_epi32(x4[51], x4[52]);
  x5[56] = _mm256_sub_epi32(x4[63], x4[56]);
  x5[63] = _mm256_add_epi32(x4[63], x4[56]);
  x5[57] = _mm256_sub_epi32(x4[62], x4[57]);
  x5[62] = _mm256_add_epi32(x4[62], x4[57]);
  x5[58] = _mm256_sub_epi32(x4[61], x4[58]);
  x5[61] = _mm256_add_epi32(x4[61], x4[58]);
  x5[59] = _mm256_sub_epi32(x4[60], x4[59]);
  x5[60] = _mm256_add_epi32(x4[60], x4[59]);

  // stage 6
  __m256i x6[64];
  btf_32_type0_avx2_new(cospi_p32, cospi_p32, x5[0], x5[1], x6[0], x6[1],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p48, cospi_p16, x5[2], x5[3], x6[2], x6[3],
                        __rounding, cos_bit);
  x6[4] = _mm256_add_epi32(x5[4], x5[5]);
  x6[5] = _mm256_sub_epi32(x5[4], x5[5]);
  x6[6] = _mm256_sub_epi32(x5[7], x5[6]);
  x6[7] = _mm256_add_epi32(x5[7], x5[6]);
  x6[8] = x5[8];
  btf_32_type0_avx2_new(cospi_m16, cospi_p48, x5[9], x5[14], x6[9], x6[14],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m48, cospi_m16, x5[10], x5[13], x6[10], x6[13],
                        __rounding, cos_bit);
  x6[11] = x5[11];
  x6[12] = x5[12];
  x6[15] = x5[15];
  x6[16] = _mm256_add_epi32(x5[16], x5[19]);
  x6[19] = _mm256_sub_epi32(x5[16], x5[19]);
  x6[17] = _mm256_add_epi32(x5[17], x5[18]);
  x6[18] = _mm256_sub_epi32(x5[17], x5[18]);
  x6[20] = _mm256_sub_epi32(x5[23], x5[20]);
  x6[23] = _mm256_add_epi32(x5[23], x5[20]);
  x6[21] = _mm256_sub_epi32(x5[22], x5[21]);
  x6[22] = _mm256_add_epi32(x5[22], x5[21]);
  x6[24] = _mm256_add_epi32(x5[24], x5[27]);
  x6[27] = _mm256_sub_epi32(x5[24], x5[27]);
  x6[25] = _mm256_add_epi32(x5[25], x5[26]);
  x6[26] = _mm256_sub_epi32(x5[25], x5[26]);
  x6[28] = _mm256_sub_epi32(x5[31], x5[28]);
  x6[31] = _mm256_add_epi32(x5[31], x5[28]);
  x6[29] = _mm256_sub_epi32(x5[30], x5[29]);
  x6[30] = _mm256_add_epi32(x5[30], x5[29]);
  x6[32] = x5[32];
  x6[33] = x5[33];
  btf_32_type0_avx2_new(cospi_m08, cospi_p56, x5[34], x5[61], x6[34], x6[61],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m08, cospi_p56, x5[35], x5[60], x6[35], x6[60],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m56, cospi_m08, x5[36], x5[59], x6[36], x6[59],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m56, cospi_m08, x5[37], x5[58], x6[37], x6[58],
                        __rounding, cos_bit);
  x6[38] = x5[38];
  x6[39] = x5[39];
  x6[40] = x5[40];
  x6[41] = x5[41];
  btf_32_type0_avx2_new(cospi_m40, cospi_p24, x5[42], x5[53], x6[42], x6[53],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m40, cospi_p24, x5[43], x5[52], x6[43], x6[52],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m24, cospi_m40, x5[44], x5[51], x6[44], x6[51],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m24, cospi_m40, x5[45], x5[50], x6[45], x6[50],
                        __rounding, cos_bit);
  x6[46] = x5[46];
  x6[47] = x5[47];
  x6[48] = x5[48];
  x6[49] = x5[49];
  x6[54] = x5[54];
  x6[55] = x5[55];
  x6[56] = x5[56];
  x6[57] = x5[57];
  x6[62] = x5[62];
  x6[63] = x5[63];

  // stage 7
  __m256i x7[64];
  x7[0] = x6[0];
  x7[1] = x6[1];
  x7[2] = x6[2];
  x7[3] = x6[3];
  btf_32_type1_avx2_new(cospi_p56, cospi_p08, x6[4], x6[7], x7[4], x7[7],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p24, cospi_p40, x6[5], x6[6], x7[5], x7[6],
                        __rounding, cos_bit);
  x7[8] = _mm256_add_epi32(x6[8], x6[9]);
  x7[9] = _mm256_sub_epi32(x6[8], x6[9]);
  x7[10] = _mm256_sub_epi32(x6[11], x6[10]);
  x7[11] = _mm256_add_epi32(x6[11], x6[10]);
  x7[12] = _mm256_add_epi32(x6[12], x6[13]);
  x7[13] = _mm256_sub_epi32(x6[12], x6[13]);
  x7[14] = _mm256_sub_epi32(x6[15], x6[14]);
  x7[15] = _mm256_add_epi32(x6[15], x6[14]);
  x7[16] = x6[16];
  btf_32_type0_avx2_new(cospi_m08, cospi_p56, x6[17], x6[30], x7[17], x7[30],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m56, cospi_m08, x6[18], x6[29], x7[18], x7[29],
                        __rounding, cos_bit);
  x7[19] = x6[19];
  x7[20] = x6[20];
  btf_32_type0_avx2_new(cospi_m40, cospi_p24, x6[21], x6[26], x7[21], x7[26],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m24, cospi_m40, x6[22], x6[25], x7[22], x7[25],
                        __rounding, cos_bit);
  x7[23] = x6[23];
  x7[24] = x6[24];
  x7[27] = x6[27];
  x7[28] = x6[28];
  x7[31] = x6[31];
  x7[32] = _mm256_add_epi32(x6[32], x6[35]);
  x7[35] = _mm256_sub_epi32(x6[32], x6[35]);
  x7[33] = _mm256_add_epi32(x6[33], x6[34]);
  x7[34] = _mm256_sub_epi32(x6[33], x6[34]);
  x7[36] = _mm256_sub_epi32(x6[39], x6[36]);
  x7[39] = _mm256_add_epi32(x6[39], x6[36]);
  x7[37] = _mm256_sub_epi32(x6[38], x6[37]);
  x7[38] = _mm256_add_epi32(x6[38], x6[37]);
  x7[40] = _mm256_add_epi32(x6[40], x6[43]);
  x7[43] = _mm256_sub_epi32(x6[40], x6[43]);
  x7[41] = _mm256_add_epi32(x6[41], x6[42]);
  x7[42] = _mm256_sub_epi32(x6[41], x6[42]);
  x7[44] = _mm256_sub_epi32(x6[47], x6[44]);
  x7[47] = _mm256_add_epi32(x6[47], x6[44]);
  x7[45] = _mm256_sub_epi32(x6[46], x6[45]);
  x7[46] = _mm256_add_epi32(x6[46], x6[45]);
  x7[48] = _mm256_add_epi32(x6[48], x6[51]);
  x7[51] = _mm256_sub_epi32(x6[48], x6[51]);
  x7[49] = _mm256_add_epi32(x6[49], x6[50]);
  x7[50] = _mm256_sub_epi32(x6[49], x6[50]);
  x7[52] = _mm256_sub_epi32(x6[55], x6[52]);
  x7[55] = _mm256_add_epi32(x6[55], x6[52]);
  x7[53] = _mm256_sub_epi32(x6[54], x6[53]);
  x7[54] = _mm256_add_epi32(x6[54], x6[53]);
  x7[56] = _mm256_add_epi32(x6[56], x6[59]);
  x7[59] = _mm256_sub_epi32(x6[56], x6[59]);
  x7[57] = _mm256_add_epi32(x6[57], x6[58]);
  x7[58] = _mm256_sub_epi32(x6[57], x6[58]);
  x7[60] = _mm256_sub_epi32(x6[63], x6[60]);
  x7[63] = _mm256_add_epi32(x6[63], x6[60]);
  x7[61] = _mm256_sub_epi32(x6[62], x6[61]);
  x7[62] = _mm256_add_epi32(x6[62], x6[61]);

  // stage 8
  __m256i x8[64];
  x8[0] = x7[0];
  x8[1] = x7[1];
  x8[2] = x7[2];
  x8[3] = x7[3];
  x8[4] = x7[4];
  x8[5] = x7[5];
  x8[6] = x7[6];
  x8[7] = x7[7];
  btf_32_type1_avx2_new(cospi_p60, cospi_p04, x7[8], x7[15], x8[8], x8[15],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p28, cospi_p36, x7[9], x7[14], x8[9], x8[14],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p44, cospi_p20, x7[10], x7[13], x8[10], x8[13],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p12, cospi_p52, x7[11], x7[12], x8[11], x8[12],
                        __rounding, cos_bit);
  x8[16] = _mm256_add_epi32(x7[16], x7[17]);
  x8[17] = _mm256_sub_epi32(x7[16], x7[17]);
  x8[18] = _mm256_sub_epi32(x7[19], x7[18]);
  x8[19] = _mm256_add_epi32(x7[19], x7[18]);
  x8[20] = _mm256_add_epi32(x7[20], x7[21]);
  x8[21] = _mm256_sub_epi32(x7[20], x7[21]);
  x8[22] = _mm256_sub_epi32(x7[23], x7[22]);
  x8[23] = _mm256_add_epi32(x7[23], x7[22]);
  x8[24] = _mm256_add_epi32(x7[24], x7[25]);
  x8[25] = _mm256_sub_epi32(x7[24], x7[25]);
  x8[26] = _mm256_sub_epi32(x7[27], x7[26]);
  x8[27] = _mm256_add_epi32(x7[27], x7[26]);
  x8[28] = _mm256_add_epi32(x7[28], x7[29]);
  x8[29] = _mm256_sub_epi32(x7[28], x7[29]);
  x8[30] = _mm256_sub_epi32(x7[31], x7[30]);
  x8[31] = _mm256_add_epi32(x7[31], x7[30]);
  x8[32] = x7[32];
  btf_32_type0_avx2_new(cospi_m04, cospi_p60, x7[33], x7[62], x8[33], x8[62],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m60, cospi_m04, x7[34], x7[61], x8[34], x8[61],
                        __rounding, cos_bit);
  x8[35] = x7[35];
  x8[36] = x7[36];
  btf_32_type0_avx2_new(cospi_m36, cospi_p28, x7[37], x7[58], x8[37], x8[58],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m28, cospi_m36, x7[38], x7[57], x8[38], x8[57],
                        __rounding, cos_bit);
  x8[39] = x7[39];
  x8[40] = x7[40];
  btf_32_type0_avx2_new(cospi_m20, cospi_p44, x7[41], x7[54], x8[41], x8[54],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m44, cospi_m20, x7[42], x7[53], x8[42], x8[53],
                        __rounding, cos_bit);
  x8[43] = x7[43];
  x8[44] = x7[44];
  btf_32_type0_avx2_new(cospi_m52, cospi_p12, x7[45], x7[50], x8[45], x8[50],
                        __rounding, cos_bit);
  btf_32_type0_avx2_new(cospi_m12, cospi_m52, x7[46], x7[49], x8[46], x8[49],
                        __rounding, cos_bit);
  x8[47] = x7[47];
  x8[48] = x7[48];
  x8[51] = x7[51];
  x8[52] = x7[52];
  x8[55] = x7[55];
  x8[56] = x7[56];
  x8[59] = x7[59];
  x8[60] = x7[60];
  x8[63] = x7[63];

  // stage 9
  __m256i x9[64];
  x9[0] = x8[0];
  x9[1] = x8[1];
  x9[2] = x8[2];
  x9[3] = x8[3];
  x9[4] = x8[4];
  x9[5] = x8[5];
  x9[6] = x8[6];
  x9[7] = x8[7];
  x9[8] = x8[8];
  x9[9] = x8[9];
  x9[10] = x8[10];
  x9[11] = x8[11];
  x9[12] = x8[12];
  x9[13] = x8[13];
  x9[14] = x8[14];
  x9[15] = x8[15];
  btf_32_type1_avx2_new(cospi_p62, cospi_p02, x8[16], x8[31], x9[16], x9[31],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p30, cospi_p34, x8[17], x8[30], x9[17], x9[30],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p46, cospi_p18, x8[18], x8[29], x9[18], x9[29],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p14, cospi_p50, x8[19], x8[28], x9[19], x9[28],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p54, cospi_p10, x8[20], x8[27], x9[20], x9[27],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p22, cospi_p42, x8[21], x8[26], x9[21], x9[26],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p38, cospi_p26, x8[22], x8[25], x9[22], x9[25],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p06, cospi_p58, x8[23], x8[24], x9[23], x9[24],
                        __rounding, cos_bit);
  x9[32] = _mm256_add_epi32(x8[32], x8[33]);
  x9[33] = _mm256_sub_epi32(x8[32], x8[33]);
  x9[34] = _mm256_sub_epi32(x8[35], x8[34]);
  x9[35] = _mm256_add_epi32(x8[35], x8[34]);
  x9[36] = _mm256_add_epi32(x8[36], x8[37]);
  x9[37] = _mm256_sub_epi32(x8[36], x8[37]);
  x9[38] = _mm256_sub_epi32(x8[39], x8[38]);
  x9[39] = _mm256_add_epi32(x8[39], x8[38]);
  x9[40] = _mm256_add_epi32(x8[40], x8[41]);
  x9[41] = _mm256_sub_epi32(x8[40], x8[41]);
  x9[42] = _mm256_sub_epi32(x8[43], x8[42]);
  x9[43] = _mm256_add_epi32(x8[43], x8[42]);
  x9[44] = _mm256_add_epi32(x8[44], x8[45]);
  x9[45] = _mm256_sub_epi32(x8[44], x8[45]);
  x9[46] = _mm256_sub_epi32(x8[47], x8[46]);
  x9[47] = _mm256_add_epi32(x8[47], x8[46]);
  x9[48] = _mm256_add_epi32(x8[48], x8[49]);
  x9[49] = _mm256_sub_epi32(x8[48], x8[49]);
  x9[50] = _mm256_sub_epi32(x8[51], x8[50]);
  x9[51] = _mm256_add_epi32(x8[51], x8[50]);
  x9[52] = _mm256_add_epi32(x8[52], x8[53]);
  x9[53] = _mm256_sub_epi32(x8[52], x8[53]);
  x9[54] = _mm256_sub_epi32(x8[55], x8[54]);
  x9[55] = _mm256_add_epi32(x8[55], x8[54]);
  x9[56] = _mm256_add_epi32(x8[56], x8[57]);
  x9[57] = _mm256_sub_epi32(x8[56], x8[57]);
  x9[58] = _mm256_sub_epi32(x8[59], x8[58]);
  x9[59] = _mm256_add_epi32(x8[59], x8[58]);
  x9[60] = _mm256_add_epi32(x8[60], x8[61]);
  x9[61] = _mm256_sub_epi32(x8[60], x8[61]);
  x9[62] = _mm256_sub_epi32(x8[63], x8[62]);
  x9[63] = _mm256_add_epi32(x8[63], x8[62]);

  // stage 10
  __m256i x10[64];
  x10[0] = x9[0];
  x10[1] = x9[1];
  x10[2] = x9[2];
  x10[3] = x9[3];
  x10[4] = x9[4];
  x10[5] = x9[5];
  x10[6] = x9[6];
  x10[7] = x9[7];
  x10[8] = x9[8];
  x10[9] = x9[9];
  x10[10] = x9[10];
  x10[11] = x9[11];
  x10[12] = x9[12];
  x10[13] = x9[13];
  x10[14] = x9[14];
  x10[15] = x9[15];
  x10[16] = x9[16];
  x10[17] = x9[17];
  x10[18] = x9[18];
  x10[19] = x9[19];
  x10[20] = x9[20];
  x10[21] = x9[21];
  x10[22] = x9[22];
  x10[23] = x9[23];
  x10[24] = x9[24];
  x10[25] = x9[25];
  x10[26] = x9[26];
  x10[27] = x9[27];
  x10[28] = x9[28];
  x10[29] = x9[29];
  x10[30] = x9[30];
  x10[31] = x9[31];
  btf_32_type1_avx2_new(cospi_p63, cospi_p01, x9[32], x9[63], x10[32], x10[63],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p31, cospi_p33, x9[33], x9[62], x10[33], x10[62],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p47, cospi_p17, x9[34], x9[61], x10[34], x10[61],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p15, cospi_p49, x9[35], x9[60], x10[35], x10[60],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p55, cospi_p09, x9[36], x9[59], x10[36], x10[59],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p23, cospi_p41, x9[37], x9[58], x10[37], x10[58],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p39, cospi_p25, x9[38], x9[57], x10[38], x10[57],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p07, cospi_p57, x9[39], x9[56], x10[39], x10[56],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p59, cospi_p05, x9[40], x9[55], x10[40], x10[55],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p27, cospi_p37, x9[41], x9[54], x10[41], x10[54],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p43, cospi_p21, x9[42], x9[53], x10[42], x10[53],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p11, cospi_p53, x9[43], x9[52], x10[43], x10[52],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p51, cospi_p13, x9[44], x9[51], x10[44], x10[51],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p19, cospi_p45, x9[45], x9[50], x10[45], x10[50],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p35, cospi_p29, x9[46], x9[49], x10[46], x10[49],
                        __rounding, cos_bit);
  btf_32_type1_avx2_new(cospi_p03, cospi_p61, x9[47], x9[48], x10[47], x10[48],
                        __rounding, cos_bit);

  // stage 11
  output[0] = x10[0];
  output[1] = x10[32];
  output[2] = x10[16];
  output[3] = x10[48];
  output[4] = x10[8];
  output[5] = x10[40];
  output[6] = x10[24];
  output[7] = x10[56];
  output[8] = x10[4];
  output[9] = x10[36];
  output[10] = x10[20];
  output[11] = x10[52];
  output[12] = x10[12];
  output[13] = x10[44];
  output[14] = x10[28];
  output[15] = x10[60];
  output[16] = x10[2];
  output[17] = x10[34];
  output[18] = x10[18];
  output[19] = x10[50];
  output[20] = x10[10];
  output[21] = x10[42];
  output[22] = x10[26];
  output[23] = x10[58];
  output[24] = x10[6];
  output[25] = x10[38];
  output[26] = x10[22];
  output[27] = x10[54];
  output[28] = x10[14];
  output[29] = x10[46];
  output[30] = x10[30];
  output[31] = x10[62];
  output[32] = x10[1];
  output[33] = x10[33];
  output[34] = x10[17];
  output[35] = x10[49];
  output[36] = x10[9];
  output[37] = x10[41];
  output[38] = x10[25];
  output[39] = x10[57];
  output[40] = x10[5];
  output[41] = x10[37];
  output[42] = x10[21];
  output[43] = x10[53];
  output[44] = x10[13];
  output[45] = x10[45];
  output[46] = x10[29];
  output[47] = x10[61];
  output[48] = x10[3];
  output[49] = x10[35];
  output[50] = x10[19];
  output[51] = x10[51];
  output[52] = x10[11];
  output[53] = x10[43];
  output[54] = x10[27];
  output[55] = x10[59];
  output[56] = x10[7];
  output[57] = x10[39];
  output[58] = x10[23];
  output[59] = x10[55];
  output[60] = x10[15];
  output[61] = x10[47];
  output[62] = x10[31];
  output[63] = x10[63];
}

static INLINE void fadst16x16_new_avx2(const __m256i *input, __m256i *output,
                                       int8_t cos_bit) {
  const int32_t *cospi = cospi_arr(cos_bit);
  const __m256i __zero = _mm256_setzero_si256();
  const __m256i __rounding = _mm256_set1_epi32(1 << (cos_bit - 1));

  __m256i cospi_p32_p32 = pair_set_w16_epi16(cospi[32], cospi[32]);
  __m256i cospi_p32_m32 = pair_set_w16_epi16(cospi[32], -cospi[32]);
  __m256i cospi_p16_p48 = pair_set_w16_epi16(cospi[16], cospi[48]);
  __m256i cospi_p48_m16 = pair_set_w16_epi16(cospi[48], -cospi[16]);
  __m256i cospi_m48_p16 = pair_set_w16_epi16(-cospi[48], cospi[16]);
  __m256i cospi_p08_p56 = pair_set_w16_epi16(cospi[8], cospi[56]);
  __m256i cospi_p56_m08 = pair_set_w16_epi16(cospi[56], -cospi[8]);
  __m256i cospi_p40_p24 = pair_set_w16_epi16(cospi[40], cospi[24]);
  __m256i cospi_p24_m40 = pair_set_w16_epi16(cospi[24], -cospi[40]);
  __m256i cospi_m56_p08 = pair_set_w16_epi16(-cospi[56], cospi[8]);
  __m256i cospi_m24_p40 = pair_set_w16_epi16(-cospi[24], cospi[40]);
  __m256i cospi_p02_p62 = pair_set_w16_epi16(cospi[2], cospi[62]);
  __m256i cospi_p62_m02 = pair_set_w16_epi16(cospi[62], -cospi[2]);
  __m256i cospi_p10_p54 = pair_set_w16_epi16(cospi[10], cospi[54]);
  __m256i cospi_p54_m10 = pair_set_w16_epi16(cospi[54], -cospi[10]);
  __m256i cospi_p18_p46 = pair_set_w16_epi16(cospi[18], cospi[46]);
  __m256i cospi_p46_m18 = pair_set_w16_epi16(cospi[46], -cospi[18]);
  __m256i cospi_p26_p38 = pair_set_w16_epi16(cospi[26], cospi[38]);
  __m256i cospi_p38_m26 = pair_set_w16_epi16(cospi[38], -cospi[26]);
  __m256i cospi_p34_p30 = pair_set_w16_epi16(cospi[34], cospi[30]);
  __m256i cospi_p30_m34 = pair_set_w16_epi16(cospi[30], -cospi[34]);
  __m256i cospi_p42_p22 = pair_set_w16_epi16(cospi[42], cospi[22]);
  __m256i cospi_p22_m42 = pair_set_w16_epi16(cospi[22], -cospi[42]);
  __m256i cospi_p50_p14 = pair_set_w16_epi16(cospi[50], cospi[14]);
  __m256i cospi_p14_m50 = pair_set_w16_epi16(cospi[14], -cospi[50]);
  __m256i cospi_p58_p06 = pair_set_w16_epi16(cospi[58], cospi[6]);
  __m256i cospi_p06_m58 = pair_set_w16_epi16(cospi[6], -cospi[58]);

  // stage 1
  __m256i x1[16];
  x1[0] = input[0];
  x1[1] = _mm256_subs_epi16(__zero, input[15]);
  x1[2] = _mm256_subs_epi16(__zero, input[7]);
  x1[3] = input[8];
  x1[4] = _mm256_subs_epi16(__zero, input[3]);
  x1[5] = input[12];
  x1[6] = input[4];
  x1[7] = _mm256_subs_epi16(__zero, input[11]);
  x1[8] = _mm256_subs_epi16(__zero, input[1]);
  x1[9] = input[14];
  x1[10] = input[6];
  x1[11] = _mm256_subs_epi16(__zero, input[9]);
  x1[12] = input[2];
  x1[13] = _mm256_subs_epi16(__zero, input[13]);
  x1[14] = _mm256_subs_epi16(__zero, input[5]);
  x1[15] = input[10];

  // stage 2
  __m256i x2[16];
  x2[0] = x1[0];
  x2[1] = x1[1];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x1[2], x1[3], x2[2], x2[3]);
  x2[4] = x1[4];
  x2[5] = x1[5];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x1[6], x1[7], x2[6], x2[7]);
  x2[8] = x1[8];
  x2[9] = x1[9];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x1[10], x1[11], x2[10], x2[11]);
  x2[12] = x1[12];
  x2[13] = x1[13];
  btf_16_w16_avx2(cospi_p32_p32, cospi_p32_m32, x1[14], x1[15], x2[14], x2[15]);

  // stage 3
  __m256i x3[16];
  x3[0] = _mm256_adds_epi16(x2[0], x2[2]);
  x3[2] = _mm256_subs_epi16(x2[0], x2[2]);
  x3[1] = _mm256_adds_epi16(x2[1], x2[3]);
  x3[3] = _mm256_subs_epi16(x2[1], x2[3]);
  x3[4] = _mm256_adds_epi16(x2[4], x2[6]);
  x3[6] = _mm256_subs_epi16(x2[4], x2[6]);
  x3[5] = _mm256_adds_epi16(x2[5], x2[7]);
  x3[7] = _mm256_subs_epi16(x2[5], x2[7]);
  x3[8] = _mm256_adds_epi16(x2[8], x2[10]);
  x3[10] = _mm256_subs_epi16(x2[8], x2[10]);
  x3[9] = _mm256_adds_epi16(x2[9], x2[11]);
  x3[11] = _mm256_subs_epi16(x2[9], x2[11]);
  x3[12] = _mm256_adds_epi16(x2[12], x2[14]);
  x3[14] = _mm256_subs_epi16(x2[12], x2[14]);
  x3[13] = _mm256_adds_epi16(x2[13], x2[15]);
  x3[15] = _mm256_subs_epi16(x2[13], x2[15]);

  // stage 4
  __m256i x4[16];
  x4[0] = x3[0];
  x4[1] = x3[1];
  x4[2] = x3[2];
  x4[3] = x3[3];
  btf_16_w16_avx2(cospi_p16_p48, cospi_p48_m16, x3[4], x3[5], x4[4], x4[5]);
  btf_16_w16_avx2(cospi_m48_p16, cospi_p16_p48, x3[6], x3[7], x4[6], x4[7]);
  x4[8] = x3[8];
  x4[9] = x3[9];
  x4[10] = x3[10];
  x4[11] = x3[11];
  btf_16_w16_avx2(cospi_p16_p48, cospi_p48_m16, x3[12], x3[13], x4[12], x4[13]);
  btf_16_w16_avx2(cospi_m48_p16, cospi_p16_p48, x3[14], x3[15], x4[14], x4[15]);

  // stage 5
  __m256i x5[16];
  x5[0] = _mm256_adds_epi16(x4[0], x4[4]);
  x5[4] = _mm256_subs_epi16(x4[0], x4[4]);
  x5[1] = _mm256_adds_epi16(x4[1], x4[5]);
  x5[5] = _mm256_subs_epi16(x4[1], x4[5]);
  x5[2] = _mm256_adds_epi16(x4[2], x4[6]);
  x5[6] = _mm256_subs_epi16(x4[2], x4[6]);
  x5[3] = _mm256_adds_epi16(x4[3], x4[7]);
  x5[7] = _mm256_subs_epi16(x4[3], x4[7]);
  x5[8] = _mm256_adds_epi16(x4[8], x4[12]);
  x5[12] = _mm256_subs_epi16(x4[8], x4[12]);
  x5[9] = _mm256_adds_epi16(x4[9], x4[13]);
  x5[13] = _mm256_subs_epi16(x4[9], x4[13]);
  x5[10] = _mm256_adds_epi16(x4[10], x4[14]);
  x5[14] = _mm256_subs_epi16(x4[10], x4[14]);
  x5[11] = _mm256_adds_epi16(x4[11], x4[15]);
  x5[15] = _mm256_subs_epi16(x4[11], x4[15]);

  // stage 6
  __m256i x6[16];
  x6[0] = x5[0];
  x6[1] = x5[1];
  x6[2] = x5[2];
  x6[3] = x5[3];
  x6[4] = x5[4];
  x6[5] = x5[5];
  x6[6] = x5[6];
  x6[7] = x5[7];
  btf_16_w16_avx2(cospi_p08_p56, cospi_p56_m08, x5[8], x5[9], x6[8], x6[9]);
  btf_16_w16_avx2(cospi_p40_p24, cospi_p24_m40, x5[10], x5[11], x6[10], x6[11]);
  btf_16_w16_avx2(cospi_m56_p08, cospi_p08_p56, x5[12], x5[13], x6[12], x6[13]);
  btf_16_w16_avx2(cospi_m24_p40, cospi_p40_p24, x5[14], x5[15], x6[14], x6[15]);

  // stage 7
  __m256i x7[16];
  x7[0] = _mm256_adds_epi16(x6[0], x6[8]);
  x7[8] = _mm256_subs_epi16(x6[0], x6[8]);
  x7[1] = _mm256_adds_epi16(x6[1], x6[9]);
  x7[9] = _mm256_subs_epi16(x6[1], x6[9]);
  x7[2] = _mm256_adds_epi16(x6[2], x6[10]);
  x7[10] = _mm256_subs_epi16(x6[2], x6[10]);
  x7[3] = _mm256_adds_epi16(x6[3], x6[11]);
  x7[11] = _mm256_subs_epi16(x6[3], x6[11]);
  x7[4] = _mm256_adds_epi16(x6[4], x6[12]);
  x7[12] = _mm256_subs_epi16(x6[4], x6[12]);
  x7[5] = _mm256_adds_epi16(x6[5], x6[13]);
  x7[13] = _mm256_subs_epi16(x6[5], x6[13]);
  x7[6] = _mm256_adds_epi16(x6[6], x6[14]);
  x7[14] = _mm256_subs_epi16(x6[6], x6[14]);
  x7[7] = _mm256_adds_epi16(x6[7], x6[15]);
  x7[15] = _mm256_subs_epi16(x6[7], x6[15]);

  // stage 8
  __m256i x8[16];
  btf_16_w16_avx2(cospi_p02_p62, cospi_p62_m02, x7[0], x7[1], x8[0], x8[1]);
  btf_16_w16_avx2(cospi_p10_p54, cospi_p54_m10, x7[2], x7[3], x8[2], x8[3]);
  btf_16_w16_avx2(cospi_p18_p46, cospi_p46_m18, x7[4], x7[5], x8[4], x8[5]);
  btf_16_w16_avx2(cospi_p26_p38, cospi_p38_m26, x7[6], x7[7], x8[6], x8[7]);
  btf_16_w16_avx2(cospi_p34_p30, cospi_p30_m34, x7[8], x7[9], x8[8], x8[9]);
  btf_16_w16_avx2(cospi_p42_p22, cospi_p22_m42, x7[10], x7[11], x8[10], x8[11]);
  btf_16_w16_avx2(cospi_p50_p14, cospi_p14_m50, x7[12], x7[13], x8[12], x8[13]);
  btf_16_w16_avx2(cospi_p58_p06, cospi_p06_m58, x7[14], x7[15], x8[14], x8[15]);

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

static INLINE __m256i scale_round_avx2(const __m256i a, const int scale) {
  const __m256i scale_rounding =
      pair_set_w16_epi16(scale, 1 << (NewSqrt2Bits - 1));
  const __m256i b = _mm256_madd_epi16(a, scale_rounding);
  return _mm256_srai_epi32(b, NewSqrt2Bits);
}

static INLINE void fidentity16x16_new_avx2(const __m256i *input,
                                           __m256i *output, int8_t cos_bit) {
  (void)cos_bit;
  const __m256i one = _mm256_set1_epi16(1);

  for (int i = 0; i < 16; ++i) {
    const __m256i a_lo = _mm256_unpacklo_epi16(input[i], one);
    const __m256i a_hi = _mm256_unpackhi_epi16(input[i], one);
    const __m256i b_lo = scale_round_avx2(a_lo, 2 * NewSqrt2);
    const __m256i b_hi = scale_round_avx2(a_hi, 2 * NewSqrt2);
    output[i] = _mm256_packs_epi32(b_lo, b_hi);
  }
}

static INLINE void fidentity16x32_new_avx2(const __m256i *input,
                                           __m256i *output, int8_t cos_bit) {
  (void)cos_bit;
  for (int i = 0; i < 32; ++i) {
    output[i] = _mm256_slli_epi16(input[i], 2);
  }
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
                                                      const int bit) {
  const __m256i sqrt2 = _mm256_set1_epi32(NewSqrt2);
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
static INLINE void transpose_32_8x8_avx2(int stride, const __m256i *inputA,
                                         __m256i *output) {
  __m256i temp0 = _mm256_unpacklo_epi32(inputA[0], inputA[2]);
  __m256i temp1 = _mm256_unpackhi_epi32(inputA[0], inputA[2]);
  __m256i temp2 = _mm256_unpacklo_epi32(inputA[1], inputA[3]);
  __m256i temp3 = _mm256_unpackhi_epi32(inputA[1], inputA[3]);
  __m256i temp4 = _mm256_unpacklo_epi32(inputA[4], inputA[6]);
  __m256i temp5 = _mm256_unpackhi_epi32(inputA[4], inputA[6]);
  __m256i temp6 = _mm256_unpacklo_epi32(inputA[5], inputA[7]);
  __m256i temp7 = _mm256_unpackhi_epi32(inputA[5], inputA[7]);

  __m256i t0 = _mm256_unpacklo_epi32(temp0, temp2);
  __m256i t1 = _mm256_unpackhi_epi32(temp0, temp2);
  __m256i t2 = _mm256_unpacklo_epi32(temp1, temp3);
  __m256i t3 = _mm256_unpackhi_epi32(temp1, temp3);
  __m256i t4 = _mm256_unpacklo_epi32(temp4, temp6);
  __m256i t5 = _mm256_unpackhi_epi32(temp4, temp6);
  __m256i t6 = _mm256_unpacklo_epi32(temp5, temp7);
  __m256i t7 = _mm256_unpackhi_epi32(temp5, temp7);

  output[0 * stride] = _mm256_permute2x128_si256(t0, t4, 0x20);
  output[1 * stride] = _mm256_permute2x128_si256(t1, t5, 0x20);
  output[2 * stride] = _mm256_permute2x128_si256(t2, t6, 0x20);
  output[3 * stride] = _mm256_permute2x128_si256(t3, t7, 0x20);
  output[4 * stride] = _mm256_permute2x128_si256(t0, t4, 0x31);
  output[5 * stride] = _mm256_permute2x128_si256(t1, t5, 0x31);
  output[6 * stride] = _mm256_permute2x128_si256(t2, t6, 0x31);
  output[7 * stride] = _mm256_permute2x128_si256(t3, t7, 0x31);
}

// Store 8 16 bit values. Sign extend the values.
static INLINE void store_buffer_16bit_to_32bit_w16_avx2(const __m256i *const in,
                                                        int32_t *out,
                                                        const int stride,
                                                        const int out_size) {
  for (int i = 0; i < out_size; ++i) {
    _mm256_store_si256((__m256i *)(out),
                       _mm256_cvtepi16_epi32(_mm256_castsi256_si128(in[i])));
    _mm256_store_si256(
        (__m256i *)(out + 8),
        _mm256_cvtepi16_epi32(_mm256_extracti128_si256(in[i], 1)));
    out += stride;
  }
}

static INLINE void store_rect_16bit_to_32bit_avx2(const __m256i a,
                                                  int32_t *const b) {
  const __m256i one = _mm256_set1_epi16(1);
  const __m256i a_reoder = _mm256_permute4x64_epi64(a, 0xd8);
  const __m256i a_lo = _mm256_unpacklo_epi16(a_reoder, one);
  const __m256i a_hi = _mm256_unpackhi_epi16(a_reoder, one);
  const __m256i b_lo = scale_round_avx2(a_lo, NewSqrt2);
  const __m256i b_hi = scale_round_avx2(a_hi, NewSqrt2);
  _mm256_store_si256((__m256i *)b, b_lo);
  _mm256_store_si256((__m256i *)(b + 8), b_hi);
}
static INLINE void store_rect_buffer_16bit_to_32bit_w16_avx2(
    const __m256i *const in, int32_t *const out, const int stride,
    const int out_size) {
  for (int i = 0; i < out_size; ++i) {
    store_rect_16bit_to_32bit_avx2(in[i], out + i * stride);
  }
}

static const transform_1d_avx2 col_txfm16x32_arr[TX_TYPES] = {
  fdct16x32_new_avx2,       // DCT_DCT
  NULL,                     // ADST_DCT
  NULL,                     // DCT_ADST
  NULL,                     // ADST_ADST
  NULL,                     // FLIPADST_DCT
  NULL,                     // DCT_FLIPADST
  NULL,                     // FLIPADST_FLIPADST
  NULL,                     // ADST_FLIPADST
  NULL,                     // FLIPADST_ADST
  fidentity16x32_new_avx2,  // IDTX
  fdct16x32_new_avx2,       // V_DCT
  fidentity16x32_new_avx2,  // H_DCT
  NULL,                     // V_ADST
  NULL,                     // H_ADST
  NULL,                     // V_FLIPADST
  NULL                      // H_FLIPADST
};

static const transform_1d_avx2 row_txfm16x32_arr[TX_TYPES] = {
  fdct16x32_new_avx2,       // DCT_DCT
  NULL,                     // ADST_DCT
  NULL,                     // DCT_ADST
  NULL,                     // ADST_ADST
  NULL,                     // FLIPADST_DCT
  NULL,                     // DCT_FLIPADST
  NULL,                     // FLIPADST_FLIPADST
  NULL,                     // ADST_FLIPADST
  NULL,                     // FLIPADST_ADST
  fidentity16x32_new_avx2,  // IDTX
  fidentity16x32_new_avx2,  // V_DCT
  fdct16x32_new_avx2,       // H_DCT
  NULL,                     // V_ADST
  NULL,                     // H_ADST
  NULL,                     // V_FLIPADST
  NULL                      // H_FLIPADST
};

static const transform_1d_avx2 col_txfm16x16_arr[TX_TYPES] = {
  fdct16x16_new_avx2,       // DCT_DCT
  fadst16x16_new_avx2,      // ADST_DCT
  fdct16x16_new_avx2,       // DCT_ADST
  fadst16x16_new_avx2,      // ADST_ADST
  fadst16x16_new_avx2,      // FLIPADST_DCT
  fdct16x16_new_avx2,       // DCT_FLIPADST
  fadst16x16_new_avx2,      // FLIPADST_FLIPADST
  fadst16x16_new_avx2,      // ADST_FLIPADST
  fadst16x16_new_avx2,      // FLIPADST_ADST
  fidentity16x16_new_avx2,  // IDTX
  fdct16x16_new_avx2,       // V_DCT
  fidentity16x16_new_avx2,  // H_DCT
  fadst16x16_new_avx2,      // V_ADST
  fidentity16x16_new_avx2,  // H_ADST
  fadst16x16_new_avx2,      // V_FLIPADST
  fidentity16x16_new_avx2   // H_FLIPADST
};

static const transform_1d_avx2 row_txfm16x16_arr[TX_TYPES] = {
  fdct16x16_new_avx2,       // DCT_DCT
  fdct16x16_new_avx2,       // ADST_DCT
  fadst16x16_new_avx2,      // DCT_ADST
  fadst16x16_new_avx2,      // ADST_ADST
  fdct16x16_new_avx2,       // FLIPADST_DCT
  fadst16x16_new_avx2,      // DCT_FLIPADST
  fadst16x16_new_avx2,      // FLIPADST_FLIPADST
  fadst16x16_new_avx2,      // ADST_FLIPADST
  fadst16x16_new_avx2,      // FLIPADST_ADST
  fidentity16x16_new_avx2,  // IDTX
  fidentity16x16_new_avx2,  // V_DCT
  fdct16x16_new_avx2,       // H_DCT
  fidentity16x16_new_avx2,  // V_ADST
  fadst16x16_new_avx2,      // H_ADST
  fidentity16x16_new_avx2,  // V_FLIPADST
  fadst16x16_new_avx2       // H_FLIPADST
};

static void lowbd_fwd_txfm2d_16x16_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  const TX_SIZE tx_size = TX_16X16;
  __m256i buf0[16], buf1[16];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = col_txfm16x16_arr[tx_type];
  const transform_1d_avx2 row_txfm = row_txfm16x16_arr[tx_type];
  int ud_flip, lr_flip;

  get_flip_cfg(tx_type, &ud_flip, &lr_flip);
  const int32_t i = 0;
  if (ud_flip) {
    load_buffer_16bit_to_16bit_flip_avx2(input + 16 * i, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
  }
  round_shift_16bit_w16_avx2(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit_w16_avx2(buf0, height, shift[1]);
  transpose_16bit_16x16_avx2(buf0, buf1 + 0 * width + 16 * i);

  __m256i *buf;
  if (lr_flip) {
    buf = buf0;
    flip_buf_avx2(buf1 + width * i, buf, width);
  } else {
    buf = buf1 + width * i;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit_w16_avx2(buf, width, shift[2]);
  transpose_16bit_16x16_avx2(buf, buf);
  store_buffer_16bit_to_32bit_w16_avx2(buf, output + 16 * width * i, width, 16);
}

static void lowbd_fwd_txfm2d_32x32_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  const TX_SIZE tx_size = TX_32X32;
  __m256i buf0[32], buf1[128];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = col_txfm16x32_arr[tx_type];
  const transform_1d_avx2 row_txfm = row_txfm16x32_arr[tx_type];

  int ud_flip, lr_flip;
  get_flip_cfg(tx_type, &ud_flip, &lr_flip);

  for (int i = 0; i < 2; i++) {
    if (ud_flip) {
      load_buffer_16bit_to_16bit_flip_avx2(input + 16 * i, stride, buf0,
                                           height);
    } else {
      load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    }
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    transpose_16bit_16x16_avx2(buf0 + 0 * 16, buf1 + 0 * width + 16 * i);
    transpose_16bit_16x16_avx2(buf0 + 1 * 16, buf1 + 1 * width + 16 * i);
  }

  for (int i = 0; i < 2; i++) {
    __m256i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_avx2(buf1 + width * i, buf, width);
    } else {
      buf = buf1 + width * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit_w16_avx2(buf, width, shift[2]);
    transpose_16bit_16x16_avx2(buf, buf);
    store_buffer_16bit_to_32bit_w16_avx2(buf, output + 16 * width * i, width,
                                         16);
    transpose_16bit_16x16_avx2(buf + 16, buf + 16);
    store_buffer_16bit_to_32bit_w16_avx2(buf + 16, output + 16 * width * i + 16,
                                         width, 16);
  }
}

static void lowbd_fwd_txfm2d_64x64_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_64X64;
  __m256i buf0[64], buf1[256];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = fdct16x64_new_avx2;
  const int width_div16 = (width >> 4);
  const int height_div16 = (height >> 4);

  for (int i = 0; i < width_div16; i++) {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    for (int j = 0; j < AOMMIN(2, height_div16); ++j) {
      transpose_16bit_16x16_avx2(buf0 + j * 16, buf1 + j * width + 16 * i);
    }
  }

  for (int i = 0; i < AOMMIN(2, height_div16); i++) {
    __m256i bufA[64];
    __m256i bufB[64];
    __m128i *buf = (__m128i *)(buf1 + width * i);
    for (int j = 0; j < width; ++j) {
      bufA[j] = _mm256_cvtepi16_epi32(buf[j * 2]);
      bufB[j] = _mm256_cvtepi16_epi32(buf[j * 2 + 1]);
    }
    av1_fdct64_new_avx2(bufA, bufA, cos_bit_row);
    av1_fdct64_new_avx2(bufB, bufB, cos_bit_row);
    av1_round_shift_array_32_avx2(bufA, bufA, 32, -shift[2]);
    av1_round_shift_array_32_avx2(bufB, bufB, 32, -shift[2]);

    int32_t *output8 = output + 16 * 32 * i;
    for (int j = 0; j < 4; ++j) {
      __m256i *out = (__m256i *)(output8 + 8 * j);
      transpose_32_8x8_avx2(4, bufA + 8 * j, out);
      transpose_32_8x8_avx2(4, bufB + 8 * j, out + 8 * 4);
    }
  }
}

static void lowbd_fwd_txfm2d_16x32_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  const TX_SIZE tx_size = TX_16X32;
  __m256i buf0[32], buf1[32];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = col_txfm16x32_arr[tx_type];
  const transform_1d_avx2 row_txfm = row_txfm16x16_arr[tx_type];

  int ud_flip, lr_flip;
  get_flip_cfg(tx_type, &ud_flip, &lr_flip);

  if (ud_flip) {
    load_buffer_16bit_to_16bit_flip_avx2(input, stride, buf0, height);
  } else {
    load_buffer_16bit_to_16bit_avx2(input, stride, buf0, height);
  }
  round_shift_16bit_w16_avx2(buf0, height, shift[0]);
  col_txfm(buf0, buf0, cos_bit_col);
  round_shift_16bit_w16_avx2(buf0, height, shift[1]);
  transpose_16bit_16x16_avx2(buf0, buf1);
  transpose_16bit_16x16_avx2(buf0 + 16, buf1 + 16);

  for (int i = 0; i < 2; i++) {
    __m256i *buf;
    if (lr_flip) {
      buf = buf0;
      flip_buf_avx2(buf1 + width * i, buf, width);
    } else {
      buf = buf1 + width * i;
    }
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit_w16_avx2(buf, width, shift[2]);
    transpose_16bit_16x16_avx2(buf, buf);
    store_rect_buffer_16bit_to_32bit_w16_avx2(buf, output + 16 * width * i,
                                              width, 16);
  }
}

static void lowbd_fwd_txfm2d_32x16_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  __m256i buf0[32], buf1[64];
  const int8_t *shift = fwd_txfm_shift_ls[TX_32X16];
  const int txw_idx = get_txw_idx(TX_32X16);
  const int txh_idx = get_txh_idx(TX_32X16);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = 32;
  const int height = 16;
  const transform_1d_avx2 col_txfm = col_txfm16x16_arr[tx_type];
  const transform_1d_avx2 row_txfm = row_txfm16x32_arr[tx_type];

  int ud_flip, lr_flip;
  get_flip_cfg(tx_type, &ud_flip, &lr_flip);

  for (int i = 0; i < 2; i++) {
    if (ud_flip) {
      load_buffer_16bit_to_16bit_flip_avx2(input + 16 * i, stride, buf0,
                                           height);
    } else {
      load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    }
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    transpose_16bit_16x16_avx2(buf0, buf1 + 0 * width + 16 * i);
  }

  __m256i *buf;
  if (lr_flip) {
    buf = buf0;
    flip_buf_avx2(buf1, buf, width);
  } else {
    buf = buf1;
  }
  row_txfm(buf, buf, cos_bit_row);
  round_shift_16bit_w16_avx2(buf, width, shift[2]);
  transpose_16bit_16x16_avx2(buf, buf);
  store_rect_buffer_16bit_to_32bit_w16_avx2(buf, output, width, 16);

  transpose_16bit_16x16_avx2(buf + 16, buf + 16);
  store_rect_buffer_16bit_to_32bit_w16_avx2(buf + 16, output + 16, width, 16);
}

static void lowbd_fwd_txfm2d_64x32_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  const TX_SIZE tx_size = TX_64X32;
  __m256i buf0[64], buf1[256];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = col_txfm16x32_arr[tx_type];
  const int width_div16 = (width >> 4);
  const int height_div16 = (height >> 4);

  for (int i = 0; i < width_div16; i++) {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    for (int j = 0; j < AOMMIN(4, height_div16); ++j) {
      transpose_16bit_16x16_avx2(buf0 + j * 16, buf1 + j * width + 16 * i);
    }
  }
  assert(tx_type == DCT_DCT);
  for (int i = 0; i < AOMMIN(2, height_div16); i++) {
    __m256i bufA[64];
    __m256i bufB[64];
    __m128i *buf = (__m128i *)(buf1 + width * i);
    for (int j = 0; j < width; ++j) {
      bufA[j] = _mm256_cvtepi16_epi32(buf[j * 2]);
      bufB[j] = _mm256_cvtepi16_epi32(buf[j * 2 + 1]);
    }
    av1_fdct64_new_avx2(bufA, bufA, cos_bit_row);
    av1_fdct64_new_avx2(bufB, bufB, cos_bit_row);
    av1_round_shift_rect_array_32_avx2(bufA, bufA, 32, -shift[2]);
    av1_round_shift_rect_array_32_avx2(bufB, bufB, 32, -shift[2]);

    int32_t *output8 = output + 16 * 32 * i;
    for (int j = 0; j < 4; ++j) {
      __m256i *out = (__m256i *)(output8 + 8 * j);
      transpose_32_8x8_avx2(4, bufA + 8 * j, out);
      transpose_32_8x8_avx2(4, bufB + 8 * j, out + 8 * 4);
    }
  }
}

static void lowbd_fwd_txfm2d_32x64_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_32X64;
  __m256i buf0[64], buf1[256];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = fdct16x64_new_avx2;
  const int width_div16 = (width >> 4);
  const int height_div16 = (height >> 4);

  for (int i = 0; i < width_div16; i++) {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    for (int j = 0; j < AOMMIN(2, height_div16); ++j) {
      transpose_16bit_16x16_avx2(buf0 + j * 16, buf1 + j * width + 16 * i);
    }
  }

  for (int i = 0; i < AOMMIN(2, height_div16); i++) {
    __m256i bufA[32];
    __m256i bufB[32];
    __m128i *buf = (__m128i *)(buf1 + width * i);
    for (int j = 0; j < width; ++j) {
      bufA[j] = _mm256_cvtepi16_epi32(buf[j * 2]);
      bufB[j] = _mm256_cvtepi16_epi32(buf[j * 2 + 1]);
    }
    av1_fdct32_new_avx2(bufA, bufA, cos_bit_row);
    av1_fdct32_new_avx2(bufB, bufB, cos_bit_row);
    av1_round_shift_rect_array_32_avx2(bufA, bufA, 32, -shift[2]);
    av1_round_shift_rect_array_32_avx2(bufB, bufB, 32, -shift[2]);

    int32_t *output8 = output + 16 * 32 * i;
    for (int j = 0; j < 4; ++j) {
      __m256i *out = (__m256i *)(output8 + 8 * j);
      transpose_32_8x8_avx2(4, bufA + 8 * j, out);
      transpose_32_8x8_avx2(4, bufB + 8 * j, out + 8 * 4);
    }
  }
}

static void lowbd_fwd_txfm2d_16x64_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_16X64;
  __m256i buf0[64], buf1[64];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = fdct16x64_new_avx2;
  const transform_1d_avx2 row_txfm = fdct16x16_new_avx2;
  const int width_div16 = (width >> 4);
  const int height_div16 = (height >> 4);

  for (int i = 0; i < width_div16; i++) {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    for (int j = 0; j < height_div16; ++j) {
      transpose_16bit_16x16_avx2(buf0 + j * 16, buf1 + j * width + 16 * i);
    }
  }

  for (int i = 0; i < AOMMIN(4, height_div16); i++) {
    __m256i *buf = buf1 + width * i;
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit_w16_avx2(buf, width, shift[2]);
    int32_t *output16 = output + 16 * width * i;
    for (int j = 0; j < width_div16; ++j) {
      __m256i *buf16 = buf + 16 * j;
      transpose_16bit_16x16_avx2(buf16, buf16);
      store_buffer_16bit_to_32bit_w16_avx2(buf16, output16 + 16 * j, width, 16);
    }
  }
  // Zero out the bottom 16x32 area.
  memset(output + 16 * 32, 0, 16 * 32 * sizeof(*output));
}

static void lowbd_fwd_txfm2d_64x16_avx2(const int16_t *input, int32_t *output,
                                        int stride, TX_TYPE tx_type, int bd) {
  (void)bd;
  (void)tx_type;
  assert(tx_type == DCT_DCT);
  const TX_SIZE tx_size = TX_64X16;
  __m256i buf0[64], buf1[64];
  const int8_t *shift = fwd_txfm_shift_ls[tx_size];
  const int txw_idx = get_txw_idx(tx_size);
  const int txh_idx = get_txh_idx(tx_size);
  const int cos_bit_col = fwd_cos_bit_col[txw_idx][txh_idx];
  const int cos_bit_row = fwd_cos_bit_row[txw_idx][txh_idx];
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const transform_1d_avx2 col_txfm = fdct16x16_new_avx2;
  const transform_1d_avx2 row_txfm = fdct16x64_new_avx2;
  const int width_div16 = (width >> 4);
  const int height_div16 = (height >> 4);

  for (int i = 0; i < width_div16; i++) {
    load_buffer_16bit_to_16bit_avx2(input + 16 * i, stride, buf0, height);
    round_shift_16bit_w16_avx2(buf0, height, shift[0]);
    col_txfm(buf0, buf0, cos_bit_col);
    round_shift_16bit_w16_avx2(buf0, height, shift[1]);
    for (int j = 0; j < height_div16; ++j) {
      transpose_16bit_16x16_avx2(buf0 + j * 16, buf1 + j * width + 16 * i);
    }
  }

  for (int i = 0; i < height_div16; i++) {
    __m256i *buf = buf1 + width * i;
    row_txfm(buf, buf, cos_bit_row);
    round_shift_16bit_w16_avx2(buf, width, shift[2]);
    int32_t *output16 = output + 16 * 32 * i;
    for (int j = 0; j < 2; ++j) {
      __m256i *buf16 = buf + 16 * j;
      transpose_16bit_16x16_avx2(buf16, buf16);
      store_buffer_16bit_to_32bit_w16_avx2(buf16, output16 + 16 * j, 32, 16);
    }
  }
}

static FwdTxfm2dFunc fwd_txfm2d_func_ls[TX_SIZES_ALL] = {
  av1_lowbd_fwd_txfm2d_4x4_sse2,   // 4x4 transform
  av1_lowbd_fwd_txfm2d_8x8_sse2,   // 8x8 transform
  lowbd_fwd_txfm2d_16x16_avx2,     // 16x16 transform
  lowbd_fwd_txfm2d_32x32_avx2,     // 32x32 transform
  lowbd_fwd_txfm2d_64x64_avx2,     // 64x64 transform
  av1_lowbd_fwd_txfm2d_4x8_sse2,   // 4x8 transform
  av1_lowbd_fwd_txfm2d_8x4_sse2,   // 8x4 transform
  av1_lowbd_fwd_txfm2d_8x16_sse2,  // 8x16 transform
  av1_lowbd_fwd_txfm2d_16x8_sse2,  // 16x8 transform
  lowbd_fwd_txfm2d_16x32_avx2,     // 16x32 transform
  lowbd_fwd_txfm2d_32x16_avx2,     // 32x16 transform
  lowbd_fwd_txfm2d_32x64_avx2,     // 32x64 transform
  lowbd_fwd_txfm2d_64x32_avx2,     // 64x32 transform
  av1_lowbd_fwd_txfm2d_4x16_sse2,  // 4x16 transform
  av1_lowbd_fwd_txfm2d_16x4_sse2,  // 16x4 transform
  av1_lowbd_fwd_txfm2d_8x32_sse2,  // 8x32 transform
  av1_lowbd_fwd_txfm2d_32x8_sse2,  // 32x8 transform
  lowbd_fwd_txfm2d_16x64_avx2,     // 16x64 transform
  lowbd_fwd_txfm2d_64x16_avx2,     // 64x16 transform
};

void av1_lowbd_fwd_txfm_avx2(const int16_t *src_diff, tran_low_t *coeff,
                             int diff_stride, TxfmParam *txfm_param) {
  FwdTxfm2dFunc fwd_txfm2d_func = fwd_txfm2d_func_ls[txfm_param->tx_size];
  if ((fwd_txfm2d_func == NULL) ||
      (txfm_param->lossless && txfm_param->tx_size == TX_4X4)) {
    av1_lowbd_fwd_txfm_c(src_diff, coeff, diff_stride, txfm_param);
  } else {
    fwd_txfm2d_func(src_diff, coeff, diff_stride, txfm_param->tx_type,
                    txfm_param->bd);
  }
}
