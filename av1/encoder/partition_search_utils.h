/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_PARTITION_SEARCH_UTILS_H_
#define AOM_AV1_ENCODER_PARTITION_SEARCH_UTILS_H_

#include "av1/encoder/encodemb.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodeframe_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE void store_pred_mv(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
  memcpy(ctx->pred_mv, x->pred_mv, sizeof(x->pred_mv));
}

static INLINE void load_pred_mv(MACROBLOCK *x,
                                const PICK_MODE_CONTEXT *const ctx) {
  memcpy(x->pred_mv, ctx->pred_mv, sizeof(x->pred_mv));
}

// Checks to see if a super block is on a horizontal image edge.
// In most cases this is the "real" edge unless there are formatting
// bars embedded in the stream.
int av1_active_h_edge(const AV1_COMP *cpi, int mi_row, int mi_step);

// Checks to see if a super block is on a vertical image edge.
// In most cases this is the "real" edge unless there are formatting
// bars embedded in the stream.
int av1_active_v_edge(const AV1_COMP *cpi, int mi_col, int mi_step);

// Record the ref frames that have been selected by square partition blocks.
void av1_update_picked_ref_frames_mask(MACROBLOCK *const x, int ref_type,
                                       BLOCK_SIZE bsize, int mib_size,
                                       int mi_row, int mi_col);

void av1_setup_block_rdmult(const struct AV1_COMP *cpi,
                            struct macroblock *const x, int mi_row, int mi_col,
                            BLOCK_SIZE bsize, AQ_MODE aq_mode,
                            MB_MODE_INFO *mbmi);

void av1_update_state(const AV1_COMP *const cpi, ThreadData *td,
                      const PICK_MODE_CONTEXT *const ctx, int mi_row,
                      int mi_col, BLOCK_SIZE bsize, RUN_TYPE dry_run);

void av1_encode_b(const AV1_COMP *const cpi, TileDataEnc *tile_data,
                  ThreadData *td, TOKENEXTRA **tp, int mi_row, int mi_col,
                  RUN_TYPE dry_run, BLOCK_SIZE bsize, PARTITION_TYPE partition,
                  const PICK_MODE_CONTEXT *const ctx, int *rate);

void av1_encode_sb(const AV1_COMP *const cpi, ThreadData *td,
                   TileDataEnc *tile_data, TOKENEXTRA **tp, int mi_row,
                   int mi_col, RUN_TYPE dry_run, BLOCK_SIZE bsize,
                   PC_TREE *pc_tree, PARTITION_TREE *ptree, int *rate);

void av1_encode_superblock(const AV1_COMP *const cpi, TileDataEnc *tile_data,
                           ThreadData *td, TOKENEXTRA **t, RUN_TYPE dry_run,
                           BLOCK_SIZE bsize, int *rate);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_ENCODEFRAME_UTILS_H_
