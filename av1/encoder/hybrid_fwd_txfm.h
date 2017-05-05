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

#ifndef AV1_ENCODER_HYBRID_FWD_TXFM_H_
#define AV1_ENCODER_HYBRID_FWD_TXFM_H_

#include "./aom_config.h"

typedef enum FwdTxfmOpt { FWD_TXFM_OPT_NORMAL } FwdTxfmOpt;

typedef struct FwdTxfmParam {
  TxType tx_type;
  TxSize tx_size;
  int lossless;
#if CONFIG_HIGHBITDEPTH
  int bd;
#endif  // CONFIG_HIGHBITDEPTH
} FwdTxfmParam;

#ifdef __cplusplus
extern "C" {
#endif

void av1_fwd_txfm(const int16_t *src_diff, TranLowT *coeff, int diff_stride,
                  FwdTxfmParam *fwd_txfm_param);

#if CONFIG_HIGHBITDEPTH
void av1_highbd_fwd_txfm(const int16_t *src_diff, TranLowT *coeff,
                         int diff_stride, FwdTxfmParam *fwd_txfm_param);
#endif  // CONFIG_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_HYBRID_FWD_TXFM_H_
