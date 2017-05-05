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

#ifndef AV1_ENCODER_RDOPT_H_
#define AV1_ENCODER_RDOPT_H_

#include "av1/common/blockd.h"

#include "av1/encoder/block.h"
#include "av1/encoder/context_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TileInfo;
struct Av1Comp;
struct Macroblock;
struct RdStats;

#if CONFIG_RD_DEBUG
static INLINE void av1_update_txb_coeff_cost(RdStats *rd_stats, int plane,
                                             TxSize tx_size, int blk_row,
                                             int blk_col, int txb_coeff_cost) {
  (void)blk_row;
  (void)blk_col;
  (void)tx_size;
  rd_stats->txb_coeff_cost[plane] += txb_coeff_cost;

#if CONFIG_VAR_TX
  {
    const int txb_h = tx_size_high_unit[tx_size];
    const int txb_w = tx_size_wide_unit[tx_size];
    int idx, idy;
    for (idy = 0; idy < txb_h; ++idy)
      for (idx = 0; idx < txb_w; ++idx)
        rd_stats->txb_coeff_cost_map[plane][blk_row + idy][blk_col + idx] = 0;

    rd_stats->txb_coeff_cost_map[plane][blk_row][blk_col] = txb_coeff_cost;
  }
  assert(blk_row < TXB_COEFF_COST_MAP_SIZE);
  assert(blk_col < TXB_COEFF_COST_MAP_SIZE);
#endif
}
#endif

typedef enum OutputStatus {
  OUTPUT_HAS_PREDICTED_PIXELS,
  OUTPUT_HAS_DECODED_PIXELS
} OutputStatus;

void av1_dist_block(const Av1Comp *cpi, Macroblock *x, int plane,
                    BlockSize plane_bsize, int block, int blk_row, int blk_col,
                    TxSize tx_size, int64_t *out_dist, int64_t *out_sse,
                    OutputStatus output_status);

#if !CONFIG_PVQ || CONFIG_VAR_TX
int av1_cost_coeffs(const Av1Comp *const cpi, Macroblock *x, int plane,
                    int block, TxSize tx_size, const ScanOrder *scan_order,
                    const EntropyContext *a, const EntropyContext *l,
                    int use_fast_coef_costing);
#endif
void av1_rd_pick_intra_mode_sb(const struct Av1Comp *cpi, struct Macroblock *x,
                               struct RdStats *rd_cost, BlockSize bsize,
                               PickModeContext *ctx, int64_t best_rd);

unsigned int av1_get_sby_perpixel_variance(const Av1Comp *cpi,
                                           const struct Buf2d *ref,
                                           BlockSize bs);
#if CONFIG_HIGHBITDEPTH
unsigned int av1_high_get_sby_perpixel_variance(const Av1Comp *cpi,
                                                const struct Buf2d *ref,
                                                BlockSize bs, int bd);
#endif

void av1_rd_pick_inter_mode_sb(const struct Av1Comp *cpi,
                               struct TileDataEnc *tile_data,
                               struct Macroblock *x, int mi_row, int mi_col,
                               struct RdStats *rd_cost,
#if CONFIG_SUPERTX
                               int *returnrate_nocoef,
#endif  // CONFIG_SUPERTX
                               BlockSize bsize, PickModeContext *ctx,
                               int64_t best_rd_so_far);

void av1_rd_pick_inter_mode_sb_seg_skip(const struct Av1Comp *cpi,
                                        struct TileDataEnc *tile_data,
                                        struct Macroblock *x, int mi_row,
                                        int mi_col, struct RdStats *rd_cost,
                                        BlockSize bsize, PickModeContext *ctx,
                                        int64_t best_rd_so_far);

int av1_internal_image_edge(const struct Av1Comp *cpi);
int av1_active_h_edge(const struct Av1Comp *cpi, int mi_row, int mi_step);
int av1_active_v_edge(const struct Av1Comp *cpi, int mi_col, int mi_step);
int av1_active_edge_sb(const struct Av1Comp *cpi, int mi_row, int mi_col);

void av1_rd_pick_inter_mode_sub8x8(const struct Av1Comp *cpi,
                                   struct TileDataEnc *tile_data,
                                   struct Macroblock *x, int mi_row, int mi_col,
                                   struct RdStats *rd_cost,
#if CONFIG_SUPERTX
                                   int *returnrate_nocoef,
#endif  // CONFIG_SUPERTX
                                   BlockSize bsize, PickModeContext *ctx,
                                   int64_t best_rd_so_far);

#if CONFIG_MOTION_VAR && CONFIG_NCOBMC
void av1_check_ncobmc_rd(const struct Av1Comp *cpi, struct Macroblock *x,
                         int mi_row, int mi_col);
#endif  // CONFIG_MOTION_VAR && CONFIG_NCOBMC

#if CONFIG_SUPERTX
#if CONFIG_VAR_TX
void av1_tx_block_rd_b(const Av1Comp *cpi, Macroblock *x, TxSize tx_size,
                       int blk_row, int blk_col, int plane, int block,
                       int plane_bsize, const EntropyContext *a,
                       const EntropyContext *l, RdStats *rd_stats);
#endif

void av1_txfm_rd_in_plane_supertx(Macroblock *x, const Av1Comp *cpi, int *rate,
                                  int64_t *distortion, int *skippable,
                                  int64_t *sse, int64_t ref_best_rd, int plane,
                                  BlockSize bsize, TxSize tx_size,
                                  int use_fast_coef_casting);
#endif  // CONFIG_SUPERTX

#ifdef __cplusplus
}  // extern "C"
#endif

int av1_tx_type_cost(const Av1Comp *cpi, const Macroblockd *xd, BlockSize bsize,
                     int plane, TxSize tx_size, TxType tx_type);

#endif  // AV1_ENCODER_RDOPT_H_
