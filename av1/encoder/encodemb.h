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

#ifndef AV1_ENCODER_ENCODEMB_H_
#define AV1_ENCODER_ENCODEMB_H_

#include "./aom_config.h"
#include "av1/common/onyxc_int.h"
#include "av1/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

struct OptimizeCtx {
  EntropyContext ta[MAX_MB_PLANE][2 * MAX_MIB_SIZE];
  EntropyContext tl[MAX_MB_PLANE][2 * MAX_MIB_SIZE];
};

struct EncodeBArgs {
  Av1Common *cm;
  Macroblock *x;
  struct OptimizeCtx *ctx;
  int8_t *skip;
  EntropyContext *ta;
  EntropyContext *tl;
  int8_t enable_optimize_b;
};

typedef enum Av1XformQuant {
  AV1_XFORM_QUANT_FP = 0,
  AV1_XFORM_QUANT_B = 1,
  AV1_XFORM_QUANT_DC = 2,
  AV1_XFORM_QUANT_SKIP_QUANT,
  AV1_XFORM_QUANT_TYPES,
} Av1XformQuant;

void av1_encode_sb(Av1Common *cm, Macroblock *x, BlockSize bsize, int mi_row,
                   int mi_col);
#if CONFIG_SUPERTX
void av1_encode_sb_supertx(Av1Common *cm, Macroblock *x, BlockSize bsize);
#endif  // CONFIG_SUPERTX
void av1_encode_sby_pass1(Av1Common *cm, Macroblock *x, BlockSize bsize);
void av1_xform_quant(const Av1Common *cm, Macroblock *x, int plane, int block,
                     int blk_row, int blk_col, BlockSize plane_bsize,
                     TxSize tx_size, int ctx, Av1XformQuant xform_quant_idx);

int av1_optimize_b(const Av1Common *cm, Macroblock *mb, int plane, int block,
                   TxSize tx_size, int ctx);

void av1_subtract_txb(Macroblock *x, int plane, BlockSize plane_bsize,
                      int blk_col, int blk_row, TxSize tx_size);

void av1_subtract_plane(Macroblock *x, BlockSize bsize, int plane);

void av1_set_txb_context(Macroblock *x, int plane, int block, TxSize tx_size,
                         EntropyContext *a, EntropyContext *l);

void av1_encode_block_intra(int plane, int block, int blk_row, int blk_col,
                            BlockSize plane_bsize, TxSize tx_size, void *arg);

void av1_encode_intra_block_plane(Av1Common *cm, Macroblock *x, BlockSize bsize,
                                  int plane, int enable_optimize_b, int mi_row,
                                  int mi_col);

#if CONFIG_PVQ
PvqSkipType av1_pvq_encode_helper(Macroblock *x, TranLowT *const coeff,
                                  TranLowT *ref_coeff, TranLowT *const dqcoeff,
                                  uint16_t *eob, const int16_t *quant,
                                  int plane, int tx_size, TxType tx_type,
                                  int *rate, int speed, PvqInfo *pvq_info);

void av1_store_pvq_enc_info(PvqInfo *pvq_info, int *qg, int *theta, int *k,
                            OdCoeff *y, int nb_bands, const int *off, int *size,
                            int skip_rest, int skip_dir, int bs);
#endif

#if CONFIG_CFL
void av1_predict_intra_block_encoder_facade(Macroblockd *xd, int plane,
                                            int block_idx, int blk_col,
                                            int blk_row, TxSize tx_size);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_ENCODEMB_H_
