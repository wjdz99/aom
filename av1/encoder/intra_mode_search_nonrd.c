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

#include "config/aom_dsp_rtcd.h"

#include "av1/common/reconinter.h"

#include "av1/encoder/encodemv.h"
#include "av1/encoder/nonrd_opt.h"
#include "av1/encoder/nonrd_utils.h"

static void compute_intra_yprediction(const AV1_COMMON *cm,
                                      PREDICTION_MODE mode, BLOCK_SIZE bsize,
                                      MACROBLOCK *x, MACROBLOCKD *xd) {
  const SequenceHeader *seq_params = cm->seq_params;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  uint8_t *const src_buf_base = p->src.buf;
  uint8_t *const dst_buf_base = pd->dst.buf;
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int plane = 0;
  int row, col;
  // block and transform sizes, in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // transform size varies per plane, look it up in a common way.
  const TX_SIZE tx_size = max_txsize_lookup[bsize];
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
  // If mb_to_right_edge is < 0 we are in a situation in which
  // the current block size extends into the UMV and we won't
  // visit the sub blocks that are wholly within the UMV.
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, plane);
  const int max_blocks_high = max_block_high(xd, plane_bsize, plane);
  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (row = 0; row < max_blocks_high; row += (1 << tx_size)) {
    // Skip visiting the sub blocks that are wholly within the UMV.
    for (col = 0; col < max_blocks_wide; col += (1 << tx_size)) {
      p->src.buf = &src_buf_base[4 * (row * (int64_t)src_stride + col)];
      pd->dst.buf = &dst_buf_base[4 * (row * (int64_t)dst_stride + col)];
      av1_predict_intra_block(
          xd, seq_params->sb_size, seq_params->enable_intra_edge_filter,
          block_size_wide[bsize], block_size_high[bsize], tx_size, mode, 0, 0,
          FILTER_INTRA_MODES, pd->dst.buf, dst_stride, pd->dst.buf, dst_stride,
          0, 0, plane);
    }
  }
  p->src.buf = src_buf_base;
  pd->dst.buf = dst_buf_base;
}

// Checks whether Intra mode needs to be pruned based on
// 'intra_y_mode_bsize_mask_nrd' and 'prune_hv_pred_modes_using_blksad'
// speed features.
static INLINE bool is_prune_intra_mode(
    AV1_COMP *cpi, int mode_index, int force_intra_check, BLOCK_SIZE bsize,
    uint8_t segment_id, SOURCE_SAD source_sad_nonrd,
    uint8_t color_sensitivity[MAX_MB_PLANE - 1]) {
  const PREDICTION_MODE this_mode = intra_mode_list[mode_index];
  if (mode_index > 2 || force_intra_check == 0) {
    if (!((1 << this_mode) & cpi->sf.rt_sf.intra_y_mode_bsize_mask_nrd[bsize]))
      return true;

    if (this_mode == DC_PRED) return false;

    if (!cpi->sf.rt_sf.prune_hv_pred_modes_using_src_sad) return false;

    const bool has_color_sensitivity =
        color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] &&
        color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)];
    if (has_color_sensitivity &&
        (cpi->rc.frame_source_sad > 1.1 * cpi->rc.avg_source_sad ||
         cyclic_refresh_segment_id_boosted(segment_id) ||
         source_sad_nonrd > kMedSad))
      return false;

    return true;
  }
  return false;
}

/*!\brief Estimation of RD cost of an intra mode for Non-RD optimized case.
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Calculates RD Cost for an intra mode for a single TX block using Hadamard
 * transform.
 * \param[in]    plane          Color plane
 * \param[in]    block          Index of a TX block in a prediction block
 * \param[in]    row            Row of a current TX block
 * \param[in]    col            Column of a current TX block
 * \param[in]    plane_bsize    Block size of a current prediction block
 * \param[in]    tx_size        Transform size
 * \param[in]    arg            Pointer to a structure that holds parameters
 *                              for intra mode search
 *
 * \remark Nothing is returned. Instead, best mode and RD Cost of the best mode
 * are set in \c args->rdc and \c args->mode
 */
void av1_estimate_block_intra(int plane, int block, int row, int col,
                              BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                              void *arg) {
  struct estimate_block_intra_args *const args = arg;
  AV1_COMP *const cpi = args->cpi;
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const BLOCK_SIZE bsize_tx = txsize_to_bsize[tx_size];
  uint8_t *const src_buf_base = p->src.buf;
  uint8_t *const dst_buf_base = pd->dst.buf;
  const int64_t src_stride = p->src.stride;
  const int64_t dst_stride = pd->dst.stride;
  RD_STATS this_rdc;

  (void)block;
  (void)plane_bsize;

  av1_predict_intra_block_facade(cm, xd, plane, col, row, tx_size);
  av1_invalid_rd_stats(&this_rdc);

  p->src.buf = &src_buf_base[4 * (row * src_stride + col)];
  pd->dst.buf = &dst_buf_base[4 * (row * dst_stride + col)];

  if (plane == 0) {
    av1_block_yrd(x, &this_rdc, &args->skippable, bsize_tx,
                  AOMMIN(tx_size, TX_16X16), 0);
  } else {
    av1_model_rd_for_sb_uv(cpi, bsize_tx, x, xd, &this_rdc, plane, plane);
  }

  p->src.buf = src_buf_base;
  pd->dst.buf = dst_buf_base;
  args->rdc->rate += this_rdc.rate;
  args->rdc->dist += this_rdc.dist;
}

/*!\brief Estimates best intra mode for inter mode search
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 *
 * Using heuristics based on best inter mode, block size, and other decides
 * whether to check intra modes. If so, estimates and selects best intra mode
 * from the reduced set of intra modes (max 4 intra modes checked)
 *
 * \param[in]    cpi                      Top-level encoder structure
 * \param[in]    x                        Pointer to structure holding all the
 *                                        data for the current macroblock
 * \param[in]    bsize                    Current block size
 * \param[in]    best_early_term          Flag, indicating that TX for the
 *                                        best inter mode was skipped
 * \param[in]    ref_cost_intra           Cost of signalling intra mode
 * \param[in]    reuse_prediction         Flag, indicating prediction re-use
 * \param[in]    orig_dst                 Original destination buffer
 * \param[in]    tmp_buffers              Pointer to a temporary buffers for
 *                                        prediction re-use
 * \param[out]   this_mode_pred           Pointer to store prediction buffer
 *                                        for prediction re-use
 * \param[in]    best_rdc                 Pointer to RD cost for the best
 *                                        selected intra mode
 * \param[in]    best_pickmode            Pointer to a structure containing
 *                                        best mode picked so far
 * \param[in]    ctx                      Pointer to structure holding coding
 *                                        contexts and modes for the block
 *
 * \remark Nothing is returned. Instead, calculated RD cost is placed to
 * \c best_rdc and best selected mode is placed to \c best_pickmode
 *
 */
void av1_estimate_intra_mode(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                             int best_early_term, unsigned int ref_cost_intra,
                             int reuse_prediction, struct buf_2d *orig_dst,
                             PRED_BUFFER *tmp_buffers,
                             PRED_BUFFER **this_mode_pred, RD_STATS *best_rdc,
                             BEST_PICKMODE *best_pickmode,
                             PICK_MODE_CONTEXT *ctx) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mi = xd->mi[0];
  const TxfmSearchParams *txfm_params = &x->txfm_search_params;
  const unsigned char segment_id = mi->segment_id;
  const int *const rd_threshes = cpi->rd.threshes[segment_id][bsize];
  const int *const rd_thresh_freq_fact = x->thresh_freq_fact[bsize];
  const bool is_screen_content =
      cpi->oxcf.tune_cfg.content == AOM_CONTENT_SCREEN;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  const REAL_TIME_SPEED_FEATURES *const rt_sf = &cpi->sf.rt_sf;

  const CommonQuantParams *quant_params = &cm->quant_params;

  RD_STATS this_rdc;

  int intra_cost_penalty = av1_get_intra_cost_penalty(
      quant_params->base_qindex, quant_params->y_dc_delta_q,
      cm->seq_params->bit_depth);
  int64_t inter_mode_thresh =
      RDCOST(x->rdmult, ref_cost_intra + intra_cost_penalty, 0);
  int perform_intra_pred = rt_sf->check_intra_pred_nonrd;
  int force_intra_check = 0;
  // For spatial enhancement layer: turn off intra prediction if the
  // previous spatial layer as golden ref is not chosen as best reference.
  // only do this for temporal enhancement layer and on non-key frames.
  if (cpi->svc.spatial_layer_id > 0 &&
      best_pickmode->best_ref_frame != GOLDEN_FRAME &&
      cpi->svc.temporal_layer_id > 0 &&
      !cpi->svc.layer_context[cpi->svc.temporal_layer_id].is_key_frame)
    perform_intra_pred = 0;

  int do_early_exit_rdthresh = 1;

  uint32_t spatial_var_thresh = 50;
  int motion_thresh = 32;
  // Adjust thresholds to make intra mode likely tested if the other
  // references (golden, alt) are skipped/not checked. For now always
  // adjust for svc mode.
  if (cpi->ppi->use_svc || (rt_sf->use_nonrd_altref_frame == 0 &&
                            rt_sf->nonrd_prune_ref_frame_search > 0)) {
    spatial_var_thresh = 150;
    motion_thresh = 0;
  }

  // Some adjustments to checking intra mode based on source variance.
  if (x->source_variance < spatial_var_thresh) {
    // If the best inter mode is large motion or non-LAST ref reduce intra cost
    // penalty, so intra mode is more likely tested.
    if (best_rdc->rdcost != INT64_MAX &&
        (best_pickmode->best_ref_frame != LAST_FRAME ||
         abs(mi->mv[0].as_mv.row) >= motion_thresh ||
         abs(mi->mv[0].as_mv.col) >= motion_thresh)) {
      intra_cost_penalty = intra_cost_penalty >> 2;
      inter_mode_thresh =
          RDCOST(x->rdmult, ref_cost_intra + intra_cost_penalty, 0);
      do_early_exit_rdthresh = 0;
    }
    if ((x->source_variance < AOMMAX(50, (spatial_var_thresh >> 1)) &&
         x->content_state_sb.source_sad_nonrd >= kHighSad) ||
        (is_screen_content && x->source_variance < 50 &&
         ((bsize >= BLOCK_32X32 &&
           x->content_state_sb.source_sad_nonrd != kZeroSad) ||
          x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] == 1 ||
          x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] == 1)))
      force_intra_check = 1;
    // For big blocks worth checking intra (since only DC will be checked),
    // even if best_early_term is set.
    if (bsize >= BLOCK_32X32) best_early_term = 0;
  } else if (rt_sf->source_metrics_sb_nonrd &&
             x->content_state_sb.source_sad_nonrd <= kLowSad) {
    perform_intra_pred = 0;
  }

  if (best_rdc->skip_txfm && best_pickmode->best_mode_initial_skip_flag) {
    if (rt_sf->skip_intra_pred == 1 && best_pickmode->best_mode != NEWMV)
      perform_intra_pred = 0;
    else if (rt_sf->skip_intra_pred == 2)
      perform_intra_pred = 0;
  }

  if (!(best_rdc->rdcost == INT64_MAX || force_intra_check ||
        (perform_intra_pred && !best_early_term &&
         bsize <= cpi->sf.part_sf.max_intra_bsize))) {
    return;
  }

  // Early exit based on RD cost calculated using known rate. When
  // is_screen_content is true, more bias is given to intra modes. Hence,
  // considered conservative threshold in early exit for the same.
  const int64_t known_rd = is_screen_content
                               ? CALC_BIASED_RDCOST(inter_mode_thresh)
                               : inter_mode_thresh;
  if (known_rd > best_rdc->rdcost) return;

  struct estimate_block_intra_args args = { cpi, x, DC_PRED, 1, 0 };
  TX_SIZE intra_tx_size = AOMMIN(
      AOMMIN(max_txsize_lookup[bsize],
             tx_mode_to_biggest_tx_size[txfm_params->tx_mode_search_type]),
      TX_16X16);
  if (is_screen_content && cpi->rc.high_source_sad &&
      x->source_variance > spatial_var_thresh && bsize <= BLOCK_16X16)
    intra_tx_size = TX_4X4;

  PRED_BUFFER *const best_pred = best_pickmode->best_pred;
  if (reuse_prediction && best_pred != NULL) {
    const int bh = block_size_high[bsize];
    const int bw = block_size_wide[bsize];
    if (best_pred->data == orig_dst->buf) {
      *this_mode_pred = &tmp_buffers[get_pred_buffer(tmp_buffers, 3)];
      aom_convolve_copy(best_pred->data, best_pred->stride,
                        (*this_mode_pred)->data, (*this_mode_pred)->stride, bw,
                        bh);
      best_pickmode->best_pred = *this_mode_pred;
    }
  }
  pd->dst = *orig_dst;

  for (int midx = 0; midx < RTC_INTRA_MODES; ++midx) {
    const PREDICTION_MODE this_mode = intra_mode_list[midx];
    const THR_MODES mode_index = mode_idx[INTRA_FRAME][mode_offset(this_mode)];
    const int64_t mode_rd_thresh = rd_threshes[mode_index];

    if (is_prune_intra_mode(cpi, midx, force_intra_check, bsize, segment_id,
                            x->content_state_sb.source_sad_nonrd,
                            x->color_sensitivity))
      continue;

    if (is_screen_content && rt_sf->source_metrics_sb_nonrd) {
      // For spatially flat blocks with zero motion only check
      // DC mode.
      if (x->content_state_sb.source_sad_nonrd == kZeroSad &&
          x->source_variance == 0 && this_mode != DC_PRED)
        continue;
      // Only test Intra for big blocks if spatial_variance is small.
      else if (bsize > BLOCK_32X32 && x->source_variance > 50)
        continue;
    }

    if (rd_less_than_thresh(best_rdc->rdcost, mode_rd_thresh,
                            rd_thresh_freq_fact[mode_index]) &&
        (do_early_exit_rdthresh || this_mode == SMOOTH_PRED)) {
      continue;
    }
    const BLOCK_SIZE uv_bsize =
        get_plane_block_size(bsize, xd->plane[AOM_PLANE_U].subsampling_x,
                             xd->plane[AOM_PLANE_U].subsampling_y);

    mi->mode = this_mode;
    mi->ref_frame[0] = INTRA_FRAME;
    mi->ref_frame[1] = NONE_FRAME;

    av1_invalid_rd_stats(&this_rdc);
    args.mode = this_mode;
    args.skippable = 1;
    args.rdc = &this_rdc;
    mi->tx_size = intra_tx_size;
    compute_intra_yprediction(cm, this_mode, bsize, x, xd);
    // Look into selecting tx_size here, based on prediction residual.
    av1_block_yrd(x, &this_rdc, &args.skippable, bsize, mi->tx_size, 0);
    // TODO(kyslov@) Need to account for skippable
    if (x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)]) {
      av1_foreach_transformed_block_in_plane(xd, uv_bsize, AOM_PLANE_U,
                                             av1_estimate_block_intra, &args);
    }
    if (x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)]) {
      av1_foreach_transformed_block_in_plane(xd, uv_bsize, AOM_PLANE_V,
                                             av1_estimate_block_intra, &args);
    }

    int mode_cost = 0;
    if (av1_is_directional_mode(this_mode) && av1_use_angle_delta(bsize)) {
      mode_cost +=
          x->mode_costs.angle_delta_cost[this_mode - V_PRED]
                                        [MAX_ANGLE_DELTA +
                                         mi->angle_delta[PLANE_TYPE_Y]];
    }
    if (this_mode == DC_PRED && av1_filter_intra_allowed_bsize(cm, bsize)) {
      mode_cost += x->mode_costs.filter_intra_cost[bsize][0];
    }
    this_rdc.rate += ref_cost_intra;
    this_rdc.rate += intra_cost_penalty;
    this_rdc.rate += mode_cost;
    this_rdc.rdcost = RDCOST(x->rdmult, this_rdc.rate, this_rdc.dist);

    if (is_screen_content && rt_sf->source_metrics_sb_nonrd) {
      // For blocks with low spatial variance and color sad,
      // favor the intra-modes, only on scene/slide change.
      if (cpi->rc.high_source_sad && x->source_variance < 800 &&
          (x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] ||
           x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)]))
        this_rdc.rdcost = CALC_BIASED_RDCOST(this_rdc.rdcost);
      // Otherwise bias against intra for blocks with zero
      // motion and no color, on non-scene/slide changes.
      else if (!cpi->rc.high_source_sad && x->source_variance > 0 &&
               x->content_state_sb.source_sad_nonrd == kZeroSad &&
               x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] == 0 &&
               x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] == 0)
        this_rdc.rdcost = (3 * this_rdc.rdcost) >> 1;
    }

    if (this_rdc.rdcost < best_rdc->rdcost) {
      *best_rdc = this_rdc;
      best_pickmode->best_mode = this_mode;
      best_pickmode->best_tx_size = mi->tx_size;
      best_pickmode->best_ref_frame = INTRA_FRAME;
      best_pickmode->best_second_ref_frame = NONE;
      best_pickmode->best_mode_skip_txfm = this_rdc.skip_txfm;
      if (!this_rdc.skip_txfm) {
        memcpy(ctx->blk_skip, x->txfm_search_info.blk_skip,
               sizeof(x->txfm_search_info.blk_skip[0]) * ctx->num_4x4_blk);
      }
      mi->uv_mode = this_mode;
      mi->mv[0].as_int = INVALID_MV;
      mi->mv[1].as_int = INVALID_MV;
    }
  }
  mi->tx_size = best_pickmode->best_tx_size;
}
