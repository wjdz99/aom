/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_TXFM_CFG_H_
#define AV1_TXFM_CFG_H_
// 4x4 1D

static const int8_t tx_cos_bit_col_dct_4[4] = { 13, 13, 13, 13 };
static const int8_t tx_cos_bit_row_dct_4[4] = { 13, 13, 13, 13 };
static const int8_t tx_cos_bit_col_adst_4[6] = { 13, 13, 13, 13, 13, 13 };
static const int8_t tx_cos_bit_row_adst_4[6] = { 13, 13, 13, 13, 13, 13 };

// 8x8 1D

static const int8_t tx_cos_bit_col_dct_8[6] = { 13, 13, 13, 13, 13, 13 };
static const int8_t tx_cos_bit_row_dct_8[6] = { 13, 13, 13, 13, 13, 13 };
static const int8_t tx_cos_bit_col_adst_8[8] = {
  13, 13, 13, 13, 13, 13, 13, 13
};
static const int8_t tx_cos_bit_row_adst_8[8] = {
  13, 13, 13, 13, 13, 13, 13, 13
};

// 16x16 1D

static const int8_t tx_cos_bit_col_dct_16[8] = {
  13, 13, 13, 13, 13, 13, 13, 13
};
static const int8_t tx_cos_bit_row_dct_16[8] = {
  12, 12, 12, 12, 12, 12, 12, 12
};
static const int8_t tx_cos_bit_col_adst_16[10] = { 13, 13, 13, 13, 13,
                                                   13, 13, 13, 13, 13 };
static const int8_t tx_cos_bit_row_adst_16[10] = { 12, 12, 12, 12, 12,
                                                   12, 12, 12, 12, 12 };
// 32x32 1D

static const int8_t tx_cos_bit_col_dct_32[10] = { 12, 12, 12, 12, 12,
                                                  12, 12, 12, 12, 12 };
static const int8_t tx_cos_bit_row_dct_32[10] = { 12, 12, 12, 12, 12,
                                                  12, 12, 12, 12, 12 };
static const int8_t tx_cos_bit_col_adst_32[12] = { 12, 12, 12, 12, 12, 12,
                                                   12, 12, 12, 12, 12, 12 };
static const int8_t tx_cos_bit_row_adst_32[12] = { 12, 12, 12, 12, 12, 12,
                                                   12, 12, 12, 12, 12, 12 };
// 64x64 1D

static const int8_t tx_cos_bit_col_dct_64[12] = { 13, 13, 13, 13, 13, 13,
                                                  13, 13, 13, 13, 13, 13 };
static const int8_t tx_cos_bit_row_dct_64[12] = { 12, 12, 12, 12, 12, 12,
                                                  12, 12, 12, 12, 12, 12 };

#endif  // AV1_TXFM_CFG_H_
