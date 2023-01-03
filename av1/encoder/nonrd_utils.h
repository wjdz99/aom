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

static INLINE void init_mbmi_nonrd(MB_MODE_INFO *mbmi,
                                   PREDICTION_MODE pred_mode,
                                   MV_REFERENCE_FRAME ref_frame0,
                                   MV_REFERENCE_FRAME ref_frame1,
                                   const AV1_COMMON *cm) {
  PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  mbmi->ref_mv_idx = 0;
  mbmi->mode = pred_mode;
  mbmi->uv_mode = UV_DC_PRED;
  mbmi->ref_frame[0] = ref_frame0;
  mbmi->ref_frame[1] = ref_frame1;
  pmi->palette_size[PLANE_TYPE_Y] = 0;
  pmi->palette_size[PLANE_TYPE_UV] = 0;
  mbmi->filter_intra_mode_info.use_filter_intra = 0;
  mbmi->mv[0].as_int = mbmi->mv[1].as_int = 0;
  mbmi->motion_mode = SIMPLE_TRANSLATION;
  mbmi->num_proj_ref = 1;
  mbmi->interintra_mode = 0;
  set_default_interp_filters(mbmi, cm->features.interp_filter);
}

static INLINE int get_pred_buffer(PRED_BUFFER *p, int len) {
  for (int buf_idx = 0; buf_idx < len; buf_idx++) {
    if (!p[buf_idx].in_use) {
      p[buf_idx].in_use = 1;
      return buf_idx;
    }
  }
  return -1;
}

static INLINE void free_pred_buffer(PRED_BUFFER *p) {
  if (p != NULL) p->in_use = 0;
}

#if CONFIG_INTERNAL_STATS
void av1_store_coding_context_nonrd(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                                    int mode_index);
#else
void av1_store_coding_context_nonrd(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx);
#endif  // CONFIG_INTERNAL_STATS

void av1_block_yrd(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                   BLOCK_SIZE bsize, TX_SIZE tx_size, int is_inter_mode);

void av1_block_yrd_idtx(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                        BLOCK_SIZE bsize, TX_SIZE tx_size);

int64_t av1_model_rd_for_sb_uv(AV1_COMP *cpi, BLOCK_SIZE plane_bsize,
                               MACROBLOCK *x, MACROBLOCKD *xd,
                               RD_STATS *this_rdc, int start_plane,
                               int stop_plane);

void av1_estimate_block_intra(int plane, int block, int row, int col,
                              BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                              void *arg);

void av1_estimate_intra_mode(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                             int best_early_term, unsigned int ref_cost_intra,
                             int reuse_prediction, struct buf_2d *orig_dst,
                             PRED_BUFFER *tmp_buffers,
                             PRED_BUFFER **this_mode_pred, RD_STATS *best_rdc,
                             BEST_PICKMODE *best_pickmode,
                             PICK_MODE_CONTEXT *ctx);
#endif  // AOM_AV1_ENCODER_NONRD_UTILS_H_
