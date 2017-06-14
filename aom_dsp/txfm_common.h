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

#ifndef AOM_DSP_TXFM_COMMON_H_
#define AOM_DSP_TXFM_COMMON_H_

#include "aom_dsp/aom_dsp_common.h"

// Constants and Macros used by all idct/dct functions
#define DCT_CONST_BITS 14
#define DCT_CONST_ROUNDING (1 << (DCT_CONST_BITS - 1))

#define UNIT_QUANT_SHIFT 2
#define UNIT_QUANT_FACTOR (1 << UNIT_QUANT_SHIFT)

// Constants:
//  for (int i = 1; i< 32; ++i)
//    printf("static const int cospi_%d_64 = %.0f;\n", i,
//           round(16384 * cos(i*M_PI/64)));
// Note: sin(k*Pi/64) = cos((32-k)*Pi/64)
static const tran_high_t cospi_1_64 = 16364;
static const tran_high_t cospi_2_64 = 16305;
static const tran_high_t cospi_3_64 = 16207;
static const tran_high_t cospi_4_64 = 16069;
static const tran_high_t cospi_5_64 = 15893;
static const tran_high_t cospi_6_64 = 15679;
static const tran_high_t cospi_7_64 = 15426;
static const tran_high_t cospi_8_64 = 15137;
static const tran_high_t cospi_9_64 = 14811;
static const tran_high_t cospi_10_64 = 14449;
static const tran_high_t cospi_11_64 = 14053;
static const tran_high_t cospi_12_64 = 13623;
static const tran_high_t cospi_13_64 = 13160;
static const tran_high_t cospi_14_64 = 12665;
static const tran_high_t cospi_15_64 = 12140;
static const tran_high_t cospi_16_64 = 11585;
static const tran_high_t cospi_17_64 = 11003;
static const tran_high_t cospi_18_64 = 10394;
static const tran_high_t cospi_19_64 = 9760;
static const tran_high_t cospi_20_64 = 9102;
static const tran_high_t cospi_21_64 = 8423;
static const tran_high_t cospi_22_64 = 7723;
static const tran_high_t cospi_23_64 = 7005;
static const tran_high_t cospi_24_64 = 6270;
static const tran_high_t cospi_25_64 = 5520;
static const tran_high_t cospi_26_64 = 4756;
static const tran_high_t cospi_27_64 = 3981;
static const tran_high_t cospi_28_64 = 3196;
static const tran_high_t cospi_29_64 = 2404;
static const tran_high_t cospi_30_64 = 1606;
static const tran_high_t cospi_31_64 = 804;

//  16384 * sqrt(2) * sin(kPi/9) * 2 / 3
static const tran_high_t sinpi_1_9 = 5283;
static const tran_high_t sinpi_2_9 = 9929;
static const tran_high_t sinpi_3_9 = 13377;
static const tran_high_t sinpi_4_9 = 15212;

// 16384 * sqrt(2)
static const tran_high_t Sqrt2 = 23170;

static INLINE tran_high_t fdct_round_shift(tran_high_t input) {
  tran_high_t rv = ROUND_POWER_OF_TWO(input, DCT_CONST_BITS);
  return rv;
}

#if CONFIG_LGT
// LGT4--a modified ADST4, for inter
// Inter LGT4 name: lgt4_160
// Self loops: 1.600, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000
static const tran_high_t lgtbasis4_inter[4][4] = {
  { 3809, 9358, 13567, 15834 },
  { 10673, 15348, 2189, -13513 },
  { 15057, -157, -14961, 9290 },
  { 13481, -14619, 11144, -4151 },
};

// LGT4--a modified ADST4, for intra
// Intra LGT4 name: lgt4_140
// Self loops: 1.400, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000
static const tran_high_t lgtbasis4_intra[4][4] = {
  { 4206, 9518, 13524, 15674 },
  { 11552, 14833, 1560, -13453 },
  { 15391, -1906, -14393, 9445 },
  { 12201, -14921, 12016, -4581 },
};

// LGT8--a modified ADST8, for inter
// Inter LGT8 name: lgt8_180
// Self loops: 1.800, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000
static const tran_high_t lgtbasis8_inter[8][8] = {
  { 1766, 4877, 7804, 10435, 12670, 14425, 15633, 16249 },
  { 5225, 12894, 16276, 14248, 7483, -1770, -10434, -15629 },
  { 8470, 16291, 9830, -5249, -15726, -12417, 1778, 14414 },
  { 11359, 13769, -5683, -16111, -958, 15716, 7436, -12651 },
  { 13740, 6078, -15914, -386, 16052, -5356, -14136, 10412 },
  { 15412, -4152, -10972, 15886, -6016, -9452, 16124, -7792 },
  { 15876, -12954, 5056, 4784, -12786, 15878, -12870, 4920 },
  { 12813, -14379, 14826, -14119, 12314, -9551, 6045, -2069 },
};

// LGT8--a modified ADST8, for intra
// Intra LGT8 name: lgt8_160
// Self loops: 1.600, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000
static const tran_high_t lgtbasis8_intra[8][8] = {
  { 1961, 5024, 7901, 10483, 12675, 14394, 15578, 16180 },
  { 5789, 13159, 16229, 13996, 7188, -1968, -10481, -15569 },
  { 9335, 16225, 9130, -5834, -15770, -12113, 1984, 14372 },
  { 12396, 12872, -6751, -15833, -191, 15749, 7096, -12638 },
  { 14715, 4122, -16033, 1007, 15711, -6033, -13781, 10441 },
  { 15878, -6725, -8995, 15932, -7312, -8448, 15958, -7886 },
  { 14991, -14559, 7883, 2174, -11298, 15577, -13176, 5125 },
  { 10110, -13122, 14795, -14957, 13592, -10839, 6980, -2409 },
};
#endif  // CONFIG_LGT
#endif  // AOM_DSP_TXFM_COMMON_H_
