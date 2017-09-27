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
enum { tx_cos_bit_col_dct_4 = 13 };
enum { tx_cos_bit_row_dct_4 = 13 };
enum { tx_cos_bit_col_adst_4 = 13 };
enum { tx_cos_bit_row_adst_4 = 13 };

// 8x8 1D
enum { tx_cos_bit_col_dct_8 = 13 };
enum { tx_cos_bit_row_dct_8 = 13 };
enum { tx_cos_bit_col_adst_8 = 13 };
enum { tx_cos_bit_row_adst_8 = 13 };

// 16x16 1D
enum { tx_cos_bit_col_dct_16 = 13 };
enum { tx_cos_bit_row_dct_16 = 12 };
enum { tx_cos_bit_col_adst_16 = 13 };
enum { tx_cos_bit_row_adst_16 = 12 };

// 32x32 1D
enum { tx_cos_bit_col_dct_32 = 12 };
enum { tx_cos_bit_row_dct_32 = 12 };
enum { tx_cos_bit_col_adst_32 = 12 };
enum { tx_cos_bit_row_adst_32 = 12 };

// 64x64 1D
enum { tx_cos_bit_col_dct_64 = 13 };
enum { tx_cos_bit_row_dct_64 = 12 };

#endif  // AV1_TXFM_CFG_H_
