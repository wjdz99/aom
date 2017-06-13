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

// A modified ADST4, LGT4
// LGT4 name: lgt4_orig
// Self loops: 1.000, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000
static const tran_high_t lgtbasis4[4][4] = {
  { 5283, 9929, 13377, 15212 },
  { 13377, 13377, 0, -13377 },
  { 15212, -5283, -13377, 9929 },
  { 9929, -15212, 13377, -5283 },
};

// A modified ADST8, LGT8
// LGT8 name: lgt8_orig
// Self loops: 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
// Edges: 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000
static const tran_high_t lgtbasis8[8][8] = {
  { 1606, 4756, 7723, 10394, 12665, 14449, 15679, 16305 },
  { 4756, 12665, 16305, 14449, 7723, -1606, -10394, -15679 },
  { 7723, 16305, 10394, -4756, -15679, -12665, 1606, 14449 },
  { 10394, 14449, -4756, -16305, -1606, 15679, 7723, -12665 },
  { 12665, 7723, -15679, -1606, 16305, -4756, -14449, 10394 },
  { 14449, -1606, -12665, 15679, -4756, -10394, 16305, -7723 },
  { 15679, -10394, 1606, 7723, -14449, 16305, -12665, 4756 },
  { 16305, -15679, 14449, -12665, 10394, -7723, 4756, -1606 },
};

#endif  // CONFIG_LGT

#endif  // AOM_DSP_TXFM_COMMON_H_
