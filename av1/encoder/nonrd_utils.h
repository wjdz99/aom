/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_NONRD_UTILS_H_
#define AOM_AV1_ENCODER_NONRD_UTILS_H_

void av1_estimate_intra_mode(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                             int best_early_term, unsigned int ref_cost_intra,
                             int reuse_prediction, struct buf_2d *orig_dst,
                             PRED_BUFFER *tmp_buffers,
                             PRED_BUFFER **this_mode_pred, RD_STATS *best_rdc,
                             BEST_PICKMODE *best_pickmode,
                             PICK_MODE_CONTEXT *ctx);

void av1_model_skip_for_sb_y_large_64(AV1_COMP *cpi, BLOCK_SIZE bsize,
                                      int mi_row, int mi_col, MACROBLOCK *x,
                                      MACROBLOCKD *xd, RD_STATS *rd_stats,
                                      int *early_term, int calculate_rd,
                                      int64_t best_sse,
                                      unsigned int *var_output,
                                      unsigned int var_prune_threshold);

void av1_model_skip_for_sb_y_large(AV1_COMP *cpi, BLOCK_SIZE bsize, int mi_row,
                                   int mi_col, MACROBLOCK *x, MACROBLOCKD *xd,
                                   RD_STATS *rd_stats, int *early_term,
                                   int calculate_rd, int64_t best_sse,
                                   unsigned int *var_output,
                                   unsigned int var_prune_threshold);

void av1_model_rd_for_sb_y(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                           MACROBLOCK *x, MACROBLOCKD *xd, RD_STATS *rd_stats,
                           unsigned int *var_out, int calculate_rd,
                           int *early_term);

int64_t av1_model_rd_for_sb_uv(AV1_COMP *cpi, BLOCK_SIZE plane_bsize,
                               MACROBLOCK *x, MACROBLOCKD *xd,
                               RD_STATS *this_rdc, int start_plane,
                               int stop_plane);

void av1_set_color_sensitivity(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                               int y_sad, unsigned int source_variance,
                               struct buf_2d yv12_mb[MAX_MB_PLANE]);

#if CONFIG_INTERNAL_STATS
void av1_store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                              int mode_index);
#else
void av1_store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx);
#endif  // CONFIG_INTERNAL_STATS

void av1_block_yrd(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                   BLOCK_SIZE bsize, TX_SIZE tx_size, int is_inter_mode);

void av1_block_yrd_idtx(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                        BLOCK_SIZE bsize, TX_SIZE tx_size);

#endif  // AOM_AV1_ENCODER_NONRD_UTILS_H_
