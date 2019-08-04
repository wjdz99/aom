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

#ifndef AOM_AV1_COMMON_IDCT_H_
#define AOM_AV1_COMMON_IDCT_H_

#include "config/aom_config.h"

#include "av1/common/blockd.h"
#include "av1/common/common.h"
#include "av1/common/enums.h"
#include "aom_dsp/txfm_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*transform_1d)(const tran_low_t *, tran_low_t *);

typedef struct {
  transform_1d cols, rows;  // vertical and horizontal
} transform_2d;

#define MAX_TX_SCALE 1
int av1_get_tx_scale(const TX_SIZE tx_size);

void av1_inverse_transform_block(const MACROBLOCKD *xd,
                                 const tran_low_t *dqcoeff, int plane,
                                 TX_TYPE tx_type, TX_SIZE tx_size, uint8_t *dst,
                                 int stride, int eob, int reduced_tx_set);
void av1_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd);

static INLINE const int32_t *cast_to_int32(const tran_low_t *input) {
  assert(sizeof(int32_t) == sizeof(tran_low_t));
  return (const int32_t *)input;
}

#if CONFIG_VQ4X4
#define NUM_CODEWORDS 16
#define VQ_DEBUG 0
#define VQ_BS_DEBUG 0
// TODO(kslu): refine the entropy coding of gain
#define VQ_GAIN_BITS 8
#define VQ_CODEWORD_BITS 4

// codeword4x4[i] = The i-th unit-norm codeword * 2^8
static const int32_t codeword4x4[NUM_CODEWORDS][16] = {
  { 50, 53, 46, 43, 69, 72, 60, 55, 74, 78, 65, 60, 73, 79, 68, 64 },
  { 13, 37, 61, 66, 14, 48, 84, 92, 13, 48, 86, 97, 16, 47, 82, 96 },
  { -5, -12, -14, -8, 5, 7, 12, 22, 36, 66, 89, 92, 54, 98, 124, 119 },
  { -3, -12, 13, 78, -7, -23, 25, 140, -7, -25, 25, 149, -4, -17, 20, 119 },
  { 38, 73, 91, 82, 46, 89, 121, 118, 14, 23, 32, 40, -10, -20, -29, -24 },
  { 19, 59, 23, -22, 36, 108, 43, -39, 41, 134, 60, -43, 30, 117, 60, -32 },
  { 64, 35, -17, -14, 121, 70, -28, -26, 129, 81, -30, -30, 92, 65, -22, -26 },
  { 2, 2, -7, -19, 0, 1, -16, -24, -7, -12, 7, 62, -11, -12, 82, 231 },
  { 21, 63, 141, 190, -1, 1, 16, 50, -8, -17, -29, -30, 0, -1, -5, -8 },
  { -2, -6, -5, -2, -14, -32, -38, -30, 5, 5, -12, -23, 78, 154, 150, 93 },
  { -7, -16, -27, -26, 37, 60, 67, 59, 62, 106, 134, 133, -1, -3, 3, 14 },
  { 164, 160, 57, 4, 71, 60, 17, -3, -13, -17, -9, -4, -13, -13, -6, 0 },
  { -20, 17, 68, 1, -35, 27, 111, 4, -44, 35, 139, 9, -41, 37, 141, 13 },
  { -27, -7, 4, 2, -27, -12, 1, 4, 70, 26, -12, -3, 210, 115, -18, -20 },
  { 73, -24, -7, 11, 122, -37, -11, 17, 145, -45, -11, 22, 128, -46, -9, 22 },
  { -38, -58, -65, -55, 72, 108, 124, 117, -29, -41, -44, -35, 17, 22, 24, 24 },
};

void av1_vec_dequant(const MACROBLOCKD *xd, int plane, int blk_row, int blk_col,
                     uint8_t *dst, int stride, TX_SIZE tx_size);
#endif  // CONFIG_VQ4X4

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_IDCT_H_
