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

#ifndef AV1_COMMON_IDCT_H_
#define AV1_COMMON_IDCT_H_

#include <assert.h>

#include "./aom_config.h"
#include "av1/common/blockd.h"
#include "av1/common/common.h"
#include "av1/common/enums.h"
#include "aom_dsp/inv_txfm.h"
#include "aom_dsp/txfm_common.h"
#include "aom_ports/mem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct INV_TXFM_PARAM {
#if CONFIG_ADAPT_SCAN
  const int16_t *eob_threshold;
#endif
  TX_TYPE tx_type;
  TX_SIZE tx_size;
  int eob;
  int lossless;
#if CONFIG_HIGHBITDEPTH
  int bd;
#endif
} INV_TXFM_PARAM;

typedef void (*transform_1d)(const tran_low_t *, tran_low_t *);

typedef struct {
  transform_1d cols, rows;  // vertical and horizontal
} transform_2d;

#if CONFIG_HIGHBITDEPTH
typedef void (*highbd_transform_1d)(const tran_low_t *, tran_low_t *, int bd);

typedef struct {
  highbd_transform_1d cols, rows;  // vertical and horizontal
} highbd_transform_2d;
#endif  // CONFIG_HIGHBITDEPTH

#define MAX_TX_SCALE 1
int av1_get_tx_scale(const TX_SIZE tx_size);

void av1_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob);
void av1_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                     int eob);

void av1_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type, int lossless);
void av1_inv_txfm_add_8x4(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type);
void av1_inv_txfm_add_4x8(const tran_low_t *input, uint8_t *dest, int stride,
                          int eob, TX_TYPE tx_type);
void av1_inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                      INV_TXFM_PARAM *inv_txfm_param);
void av1_inverse_transform_block(const MACROBLOCKD *xd,
                                 const tran_low_t *dqcoeff, TX_TYPE tx_type,
                                 TX_SIZE tx_size, uint8_t *dst, int stride,
                                 int eob);
void av1_inverse_transform_block_facade(MACROBLOCKD *xd, int plane, int block,
                                        int blk_row, int blk_col, int eob);
#if CONFIG_HIGHBITDEPTH
void av1_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);
void av1_highbd_idct4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);
void av1_highbd_inv_txfm_add_4x4(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type,
                                 int lossless);
void av1_highbd_inv_txfm_add_4x8(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type);
void av1_highbd_inv_txfm_add_8x4(const tran_low_t *input, uint8_t *dest,
                                 int stride, int eob, int bd, TX_TYPE tx_type);
void av1_highbd_inv_txfm_add(const tran_low_t *input, uint8_t *dest, int stride,
                             INV_TXFM_PARAM *inv_txfm_param);
#endif  // CONFIG_HIGHBITDEPTH
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_IDCT_H_
