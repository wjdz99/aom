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

#ifndef VP10_ENCODER_RDOPT_H_
#define VP10_ENCODER_RDOPT_H_

#include "av1/common/blockd.h"

#include "av1/encoder/block.h"
#include "av1/encoder/context_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TileInfo;
struct VP10_COMP;
struct macroblock;
struct RD_COST;

void av1_rd_pick_intra_mode_sb(struct VP10_COMP *cpi, struct macroblock *x,
                                struct RD_COST *rd_cost, BLOCK_SIZE bsize,
                                PICK_MODE_CONTEXT *ctx, int64_t best_rd);

unsigned int av1_get_sby_perpixel_variance(VP10_COMP *cpi,
                                            const struct buf_2d *ref,
                                            BLOCK_SIZE bs);
#if CONFIG_VPX_HIGHBITDEPTH
unsigned int av1_high_get_sby_perpixel_variance(VP10_COMP *cpi,
                                                 const struct buf_2d *ref,
                                                 BLOCK_SIZE bs, int bd);
#endif

void av1_rd_pick_inter_mode_sb(struct VP10_COMP *cpi,
                                struct TileDataEnc *tile_data,
                                struct macroblock *x, int mi_row, int mi_col,
                                struct RD_COST *rd_cost, BLOCK_SIZE bsize,
                                PICK_MODE_CONTEXT *ctx, int64_t best_rd_so_far);

void av1_rd_pick_inter_mode_sb_seg_skip(
    struct VP10_COMP *cpi, struct TileDataEnc *tile_data, struct macroblock *x,
    struct RD_COST *rd_cost, BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
    int64_t best_rd_so_far);

int av1_internal_image_edge(struct VP10_COMP *cpi);
int av1_active_h_edge(struct VP10_COMP *cpi, int mi_row, int mi_step);
int av1_active_v_edge(struct VP10_COMP *cpi, int mi_col, int mi_step);
int av1_active_edge_sb(struct VP10_COMP *cpi, int mi_row, int mi_col);

void av1_rd_pick_inter_mode_sub8x8(struct VP10_COMP *cpi,
                                    struct TileDataEnc *tile_data,
                                    struct macroblock *x, int mi_row,
                                    int mi_col, struct RD_COST *rd_cost,
                                    BLOCK_SIZE bsize, PICK_MODE_CONTEXT *ctx,
                                    int64_t best_rd_so_far);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_RDOPT_H_
