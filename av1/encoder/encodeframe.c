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

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/binary_codes_writer.h"
#include "aom_ports/mem.h"
#include "aom_ports/aom_timer.h"

#if CONFIG_MISMATCH_DEBUG
#include "aom_util/debug_util.h"
#endif  // CONFIG_MISMATCH_DEBUG

#include "av1/common/cfl.h"
#include "av1/common/common.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/idct.h"
#include "av1/common/mv.h"
#include "av1/common/mvref_common.h"
#include "av1/common/pred_common.h"
#include "av1/common/quant_common.h"
#include "av1/common/reconintra.h"
#include "av1/common/reconinter.h"
#include "av1/common/seg_common.h"
#include "av1/common/tile_common.h"
#include "av1/common/warped_motion.h"

#include "av1/encoder/allintra_vis.h"
#include "av1/encoder/aq_complexity.h"
#include "av1/encoder/aq_cyclicrefresh.h"
#include "av1/encoder/aq_variance.h"
#include "av1/encoder/global_motion_facade.h"
#include "av1/encoder/encodeframe.h"
#include "av1/encoder/encodeframe_utils.h"
#include "av1/encoder/encodemb.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/ethread.h"
#include "av1/encoder/extend.h"
#include "av1/encoder/intra_mode_search_utils.h"
#include "av1/encoder/ml.h"
#include "av1/encoder/motion_search_facade.h"
#include "av1/encoder/partition_strategy.h"
#if !CONFIG_REALTIME_ONLY
#include "av1/encoder/partition_model_weights.h"
#endif
#include "av1/encoder/partition_search.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/tokenize.h"
#include "av1/encoder/tpl_model.h"
#include "av1/encoder/var_based_part.h"

#if CONFIG_TUNE_VMAF
#include "av1/encoder/tune_vmaf.h"
#endif

/*!\cond */
// This is used as a reference when computing the source variance for the
//  purposes of activity masking.
// Eventually this should be replaced by custom no-reference routines,
//  which will be faster.
const uint8_t AV1_VAR_OFFS[MAX_SB_SIZE] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128
};

static const uint16_t AV1_HIGH_VAR_OFFS_8[MAX_SB_SIZE] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128
};

static const uint16_t AV1_HIGH_VAR_OFFS_10[MAX_SB_SIZE] = {
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4
};

static const uint16_t AV1_HIGH_VAR_OFFS_12[MAX_SB_SIZE] = {
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16, 128 * 16,
  128 * 16, 128 * 16
};
/*!\endcond */

unsigned int av1_get_sby_perpixel_variance(const AV1_COMP *cpi,
                                           const struct buf_2d *ref,
                                           BLOCK_SIZE bs) {
  unsigned int sse;
  const unsigned int var =
      cpi->ppi->fn_ptr[bs].vf(ref->buf, ref->stride, AV1_VAR_OFFS, 0, &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

unsigned int av1_high_get_sby_perpixel_variance(const AV1_COMP *cpi,
                                                const struct buf_2d *ref,
                                                BLOCK_SIZE bs, int bd) {
  unsigned int var, sse;
  assert(bd == 8 || bd == 10 || bd == 12);
  const int off_index = (bd - 8) >> 1;
  const uint16_t *high_var_offs[3] = { AV1_HIGH_VAR_OFFS_8,
                                       AV1_HIGH_VAR_OFFS_10,
                                       AV1_HIGH_VAR_OFFS_12 };
  var = cpi->ppi->fn_ptr[bs].vf(ref->buf, ref->stride,
                                CONVERT_TO_BYTEPTR(high_var_offs[off_index]), 0,
                                &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

void av1_setup_src_planes(MACROBLOCK *x, const YV12_BUFFER_CONFIG *src,
                          int mi_row, int mi_col, const int num_planes,
                          BLOCK_SIZE bsize) {
  // Set current frame pointer.
  x->e_mbd.cur_buf = src;

  // We use AOMMIN(num_planes, MAX_MB_PLANE) instead of num_planes to quiet
  // the static analysis warnings.
  for (int i = 0; i < AOMMIN(num_planes, MAX_MB_PLANE); i++) {
    const int is_uv = i > 0;
    setup_pred_plane(
        &x->plane[i].src, bsize, src->buffers[i], src->crop_widths[is_uv],
        src->crop_heights[is_uv], src->strides[is_uv], mi_row, mi_col, NULL,
        x->e_mbd.plane[i].subsampling_x, x->e_mbd.plane[i].subsampling_y);
  }
}

<<<<<<< HEAD   (de3a75 Set mv_step_param based on previous frame in speed >= 2)
=======
static EdgeInfo edge_info(const struct buf_2d *ref, const BLOCK_SIZE bsize,
                          const bool high_bd, const int bd) {
  const int width = block_size_wide[bsize];
  const int height = block_size_high[bsize];
  // Implementation requires width to be a multiple of 8. It also requires
  // height to be a multiple of 4, but this is always the case.
  assert(height % 4 == 0);
  if (width % 8 != 0) {
    EdgeInfo ei = { .magnitude = 0, .x = 0, .y = 0 };
    return ei;
  }
  return av1_edge_exists(ref->buf, ref->stride, width, height, high_bd, bd);
}

static int use_pb_simple_motion_pred_sse(const AV1_COMP *const cpi) {
  // TODO(debargha, yuec): Not in use, need to implement a speed feature
  // utilizing this data point, and replace '0' by the corresponding speed
  // feature flag.
  return 0 && !frame_is_intra_only(&cpi->common);
}

// This function will copy the winner reference mode information from block
// level (x->mbmi_ext) to frame level (cpi->mbmi_ext_frame_base). This frame
// level buffer (cpi->mbmi_ext_frame_base) will be used during bitstream
// preparation.
static INLINE void copy_winner_ref_mode_from_mbmi_ext(MACROBLOCK *const x) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = xd->mi[0];
  uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
  memcpy(x->mbmi_ext_frame->ref_mv_stack,
         x->mbmi_ext->ref_mv_stack[ref_frame_type],
         sizeof(x->mbmi_ext->ref_mv_stack[USABLE_REF_MV_STACK_SIZE]));
  memcpy(x->mbmi_ext_frame->weight, x->mbmi_ext->weight[ref_frame_type],
         sizeof(x->mbmi_ext->weight[USABLE_REF_MV_STACK_SIZE]));
  x->mbmi_ext_frame->mode_context = x->mbmi_ext->mode_context[ref_frame_type];
  x->mbmi_ext_frame->ref_mv_count = x->mbmi_ext->ref_mv_count[ref_frame_type];
  memcpy(x->mbmi_ext_frame->global_mvs, x->mbmi_ext->global_mvs,
         sizeof(x->mbmi_ext->global_mvs));
}

static AOM_INLINE void pick_sb_modes(AV1_COMP *const cpi,
                                     TileDataEnc *tile_data,
                                     MACROBLOCK *const x, int mi_row,
                                     int mi_col, RD_STATS *rd_cost,
                                     PARTITION_TYPE partition, BLOCK_SIZE bsize,
                                     PICK_MODE_CONTEXT *ctx, RD_STATS best_rd,
                                     int pick_mode_type) {
  if (best_rd.rdcost < 0) {
    ctx->rd_stats.rdcost = INT64_MAX;
    ctx->rd_stats.skip = 0;
    av1_invalid_rd_stats(rd_cost);
    return;
  }

  set_offsets(cpi, &tile_data->tile_info, x, mi_row, mi_col, bsize);

  if (ctx->rd_mode_is_ready) {
    assert(ctx->mic.sb_type == bsize);
    assert(ctx->mic.partition == partition);
    rd_cost->rate = ctx->rd_stats.rate;
    rd_cost->dist = ctx->rd_stats.dist;
    rd_cost->rdcost = ctx->rd_stats.rdcost;
    return;
  }

  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;
  struct macroblock_plane *const p = x->plane;
  struct macroblockd_plane *const pd = xd->plane;
  const AQ_MODE aq_mode = cpi->oxcf.aq_mode;
  int i;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, rd_pick_sb_modes_time);
#endif

  aom_clear_system_state();

  mbmi = xd->mi[0];
  mbmi->sb_type = bsize;
  mbmi->partition = partition;

#if CONFIG_RD_DEBUG
  mbmi->mi_row = mi_row;
  mbmi->mi_col = mi_col;
#endif

  xd->tx_type_map = x->tx_type_map;
  xd->tx_type_map_stride = mi_size_wide[bsize];

  for (i = 0; i < num_planes; ++i) {
    p[i].coeff = ctx->coeff[i];
    p[i].qcoeff = ctx->qcoeff[i];
    pd[i].dqcoeff = ctx->dqcoeff[i];
    p[i].eobs = ctx->eobs[i];
    p[i].txb_entropy_ctx = ctx->txb_entropy_ctx[i];
  }

  for (i = 0; i < 2; ++i) pd[i].color_index_map = ctx->color_index_map[i];

  ctx->skippable = 0;
  // Set to zero to make sure we do not use the previous encoded frame stats
  mbmi->skip = 0;
  // Reset skip mode flag.
  mbmi->skip_mode = 0;
  x->skip_chroma_rd =
      !is_chroma_reference(mi_row, mi_col, bsize, xd->plane[1].subsampling_x,
                           xd->plane[1].subsampling_y);

  if (is_cur_buf_hbd(xd)) {
    x->source_variance = av1_high_get_sby_perpixel_variance(
        cpi, &x->plane[0].src, bsize, xd->bd);
  } else {
    x->source_variance =
        av1_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);
  }
  if (use_pb_simple_motion_pred_sse(cpi)) {
    const MV ref_mv_full = { .row = 0, .col = 0 };
    unsigned int var = 0;
    av1_simple_motion_sse_var(cpi, x, mi_row, mi_col, bsize, ref_mv_full, 0,
                              &x->simple_motion_pred_sse, &var);
  }

  int try_filter = 1;
  // If the threshold for disabling wedge search is zero, it means the feature
  // should not be used. Use a value that will always succeed in the check.
  if (!try_filter && cpi->sf.disable_wedge_search_edge_thresh == 0) {
    x->edge_strength = UINT16_MAX;
    x->edge_strength_x = UINT16_MAX;
    x->edge_strength_y = UINT16_MAX;
  } else {
    EdgeInfo ei =
        edge_info(&x->plane[0].src, bsize, is_cur_buf_hbd(xd), xd->bd);
    x->edge_strength = ei.magnitude;
    x->edge_strength_x = ei.x;
    x->edge_strength_y = ei.y;
  }

  // Initialize default mode evaluation params
  set_mode_eval_params(cpi, x, DEFAULT_EVAL);

  // Save rdmult before it might be changed, so it can be restored later.
  const int orig_rdmult = x->rdmult;
  setup_block_rdmult(cpi, x, mi_row, mi_col, bsize, aq_mode, mbmi);
  // Set error per bit for current rdmult
  set_error_per_bit(x, x->rdmult);
  av1_rd_cost_update(x->rdmult, &best_rd);

  // Find best coding mode & reconstruct the MB so it is available
  // as a predictor for MBs that follow in the SB
  if (frame_is_intra_only(cm)) {
#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, av1_rd_pick_intra_mode_sb_time);
#endif
    av1_rd_pick_intra_mode_sb(cpi, x, rd_cost, bsize, ctx, best_rd.rdcost);
#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, av1_rd_pick_intra_mode_sb_time);
#endif
  } else {
#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, av1_rd_pick_inter_mode_sb_time);
#endif
    if (segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
      av1_rd_pick_inter_mode_sb_seg_skip(cpi, tile_data, x, mi_row, mi_col,
                                         rd_cost, bsize, ctx, best_rd.rdcost);
    } else {
      // TODO(kyslov): do the same for pick_intra_mode and
      //               pick_inter_mode_sb_seg_skip
      switch (pick_mode_type) {
        case PICK_MODE_RD:
          av1_rd_pick_inter_mode_sb(cpi, tile_data, x, rd_cost, bsize, ctx,
                                    best_rd.rdcost);
          break;
        case PICK_MODE_NONRD:
          av1_nonrd_pick_inter_mode_sb(cpi, tile_data, x, rd_cost, bsize, ctx,
                                       best_rd.rdcost);
          break;
        default: assert(0 && "Unknown pick mode type.");
      }
    }
#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, av1_rd_pick_inter_mode_sb_time);
#endif
  }

  // Examine the resulting rate and for AQ mode 2 make a segment choice.
  if ((rd_cost->rate != INT_MAX) && (aq_mode == COMPLEXITY_AQ) &&
      (bsize >= BLOCK_16X16) &&
      (cm->current_frame.frame_type == KEY_FRAME ||
       cpi->refresh_alt_ref_frame || cpi->refresh_bwd_ref_frame ||
       (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref))) {
    av1_caq_select_segment(cpi, x, bsize, mi_row, mi_col, rd_cost->rate);
  }

  x->rdmult = orig_rdmult;

  // TODO(jingning) The rate-distortion optimization flow needs to be
  // refactored to provide proper exit/return handle.
  if (rd_cost->rate == INT_MAX) rd_cost->rdcost = INT64_MAX;

  ctx->rd_stats.rate = rd_cost->rate;
  ctx->rd_stats.dist = rd_cost->dist;
  ctx->rd_stats.rdcost = rd_cost->rdcost;

#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, rd_pick_sb_modes_time);
#endif
}

static AOM_INLINE void update_inter_mode_stats(FRAME_CONTEXT *fc,
                                               FRAME_COUNTS *counts,
                                               PREDICTION_MODE mode,
                                               int16_t mode_context) {
  (void)counts;

  int16_t mode_ctx = mode_context & NEWMV_CTX_MASK;
  if (mode == NEWMV) {
#if CONFIG_ENTROPY_STATS
    ++counts->newmv_mode[mode_ctx][0];
#endif
    update_cdf(fc->newmv_cdf[mode_ctx], 0, 2);
    return;
  }

#if CONFIG_ENTROPY_STATS
  ++counts->newmv_mode[mode_ctx][1];
#endif
  update_cdf(fc->newmv_cdf[mode_ctx], 1, 2);

  mode_ctx = (mode_context >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
  if (mode == GLOBALMV) {
#if CONFIG_ENTROPY_STATS
    ++counts->zeromv_mode[mode_ctx][0];
#endif
    update_cdf(fc->zeromv_cdf[mode_ctx], 0, 2);
    return;
  }

#if CONFIG_ENTROPY_STATS
  ++counts->zeromv_mode[mode_ctx][1];
#endif
  update_cdf(fc->zeromv_cdf[mode_ctx], 1, 2);

  mode_ctx = (mode_context >> REFMV_OFFSET) & REFMV_CTX_MASK;
#if CONFIG_ENTROPY_STATS
  ++counts->refmv_mode[mode_ctx][mode != NEARESTMV];
#endif
  update_cdf(fc->refmv_cdf[mode_ctx], mode != NEARESTMV, 2);
}

static AOM_INLINE void update_palette_cdf(MACROBLOCKD *xd,
                                          const MB_MODE_INFO *const mbmi,
                                          FRAME_COUNTS *counts) {
  FRAME_CONTEXT *fc = xd->tile_ctx;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  const int palette_bsize_ctx = av1_get_palette_bsize_ctx(bsize);

  (void)counts;

  if (mbmi->mode == DC_PRED) {
    const int n = pmi->palette_size[0];
    const int palette_mode_ctx = av1_get_palette_mode_ctx(xd);

#if CONFIG_ENTROPY_STATS
    ++counts->palette_y_mode[palette_bsize_ctx][palette_mode_ctx][n > 0];
#endif
    update_cdf(fc->palette_y_mode_cdf[palette_bsize_ctx][palette_mode_ctx],
               n > 0, 2);
    if (n > 0) {
#if CONFIG_ENTROPY_STATS
      ++counts->palette_y_size[palette_bsize_ctx][n - PALETTE_MIN_SIZE];
#endif
      update_cdf(fc->palette_y_size_cdf[palette_bsize_ctx],
                 n - PALETTE_MIN_SIZE, PALETTE_SIZES);
    }
  }

  if (mbmi->uv_mode == UV_DC_PRED) {
    const int n = pmi->palette_size[1];
    const int palette_uv_mode_ctx = (pmi->palette_size[0] > 0);

#if CONFIG_ENTROPY_STATS
    ++counts->palette_uv_mode[palette_uv_mode_ctx][n > 0];
#endif
    update_cdf(fc->palette_uv_mode_cdf[palette_uv_mode_ctx], n > 0, 2);

    if (n > 0) {
#if CONFIG_ENTROPY_STATS
      ++counts->palette_uv_size[palette_bsize_ctx][n - PALETTE_MIN_SIZE];
#endif
      update_cdf(fc->palette_uv_size_cdf[palette_bsize_ctx],
                 n - PALETTE_MIN_SIZE, PALETTE_SIZES);
    }
  }
}

static AOM_INLINE void sum_intra_stats(const AV1_COMMON *const cm,
                                       FRAME_COUNTS *counts, MACROBLOCKD *xd,
                                       const MB_MODE_INFO *const mbmi,
                                       const MB_MODE_INFO *above_mi,
                                       const MB_MODE_INFO *left_mi,
                                       const int intraonly, const int mi_row,
                                       const int mi_col) {
  FRAME_CONTEXT *fc = xd->tile_ctx;
  const PREDICTION_MODE y_mode = mbmi->mode;
  (void)counts;
  const BLOCK_SIZE bsize = mbmi->sb_type;

  if (intraonly) {
#if CONFIG_ENTROPY_STATS
    const PREDICTION_MODE above = av1_above_block_mode(above_mi);
    const PREDICTION_MODE left = av1_left_block_mode(left_mi);
    const int above_ctx = intra_mode_context[above];
    const int left_ctx = intra_mode_context[left];
    ++counts->kf_y_mode[above_ctx][left_ctx][y_mode];
#endif  // CONFIG_ENTROPY_STATS
    update_cdf(get_y_mode_cdf(fc, above_mi, left_mi), y_mode, INTRA_MODES);
  } else {
#if CONFIG_ENTROPY_STATS
    ++counts->y_mode[size_group_lookup[bsize]][y_mode];
#endif  // CONFIG_ENTROPY_STATS
    update_cdf(fc->y_mode_cdf[size_group_lookup[bsize]], y_mode, INTRA_MODES);
  }

  if (av1_filter_intra_allowed(cm, mbmi)) {
    const int use_filter_intra_mode =
        mbmi->filter_intra_mode_info.use_filter_intra;
#if CONFIG_ENTROPY_STATS
    ++counts->filter_intra[mbmi->sb_type][use_filter_intra_mode];
    if (use_filter_intra_mode) {
      ++counts
            ->filter_intra_mode[mbmi->filter_intra_mode_info.filter_intra_mode];
    }
#endif  // CONFIG_ENTROPY_STATS
    update_cdf(fc->filter_intra_cdfs[mbmi->sb_type], use_filter_intra_mode, 2);
    if (use_filter_intra_mode) {
      update_cdf(fc->filter_intra_mode_cdf,
                 mbmi->filter_intra_mode_info.filter_intra_mode,
                 FILTER_INTRA_MODES);
    }
  }
  if (av1_is_directional_mode(mbmi->mode) && av1_use_angle_delta(bsize)) {
#if CONFIG_ENTROPY_STATS
    ++counts->angle_delta[mbmi->mode - V_PRED]
                         [mbmi->angle_delta[PLANE_TYPE_Y] + MAX_ANGLE_DELTA];
#endif
    update_cdf(fc->angle_delta_cdf[mbmi->mode - V_PRED],
               mbmi->angle_delta[PLANE_TYPE_Y] + MAX_ANGLE_DELTA,
               2 * MAX_ANGLE_DELTA + 1);
  }

  if (!is_chroma_reference(mi_row, mi_col, bsize,
                           xd->plane[AOM_PLANE_U].subsampling_x,
                           xd->plane[AOM_PLANE_U].subsampling_y))
    return;

  const UV_PREDICTION_MODE uv_mode = mbmi->uv_mode;
  const CFL_ALLOWED_TYPE cfl_allowed = is_cfl_allowed(xd);
#if CONFIG_ENTROPY_STATS
  ++counts->uv_mode[cfl_allowed][y_mode][uv_mode];
#endif  // CONFIG_ENTROPY_STATS
  update_cdf(fc->uv_mode_cdf[cfl_allowed][y_mode], uv_mode,
             UV_INTRA_MODES - !cfl_allowed);
  if (uv_mode == UV_CFL_PRED) {
    const int8_t joint_sign = mbmi->cfl_alpha_signs;
    const uint8_t idx = mbmi->cfl_alpha_idx;

#if CONFIG_ENTROPY_STATS
    ++counts->cfl_sign[joint_sign];
#endif
    update_cdf(fc->cfl_sign_cdf, joint_sign, CFL_JOINT_SIGNS);
    if (CFL_SIGN_U(joint_sign) != CFL_SIGN_ZERO) {
      aom_cdf_prob *cdf_u = fc->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];

#if CONFIG_ENTROPY_STATS
      ++counts->cfl_alpha[CFL_CONTEXT_U(joint_sign)][CFL_IDX_U(idx)];
#endif
      update_cdf(cdf_u, CFL_IDX_U(idx), CFL_ALPHABET_SIZE);
    }
    if (CFL_SIGN_V(joint_sign) != CFL_SIGN_ZERO) {
      aom_cdf_prob *cdf_v = fc->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];

#if CONFIG_ENTROPY_STATS
      ++counts->cfl_alpha[CFL_CONTEXT_V(joint_sign)][CFL_IDX_V(idx)];
#endif
      update_cdf(cdf_v, CFL_IDX_V(idx), CFL_ALPHABET_SIZE);
    }
  }
  if (av1_is_directional_mode(get_uv_mode(uv_mode)) &&
      av1_use_angle_delta(bsize)) {
#if CONFIG_ENTROPY_STATS
    ++counts->angle_delta[uv_mode - UV_V_PRED]
                         [mbmi->angle_delta[PLANE_TYPE_UV] + MAX_ANGLE_DELTA];
#endif
    update_cdf(fc->angle_delta_cdf[uv_mode - UV_V_PRED],
               mbmi->angle_delta[PLANE_TYPE_UV] + MAX_ANGLE_DELTA,
               2 * MAX_ANGLE_DELTA + 1);
  }
  if (av1_allow_palette(cm->allow_screen_content_tools, bsize)) {
    update_palette_cdf(xd, mbmi, counts);
  }
}

static AOM_INLINE void update_stats(const AV1_COMMON *const cm, ThreadData *td,
                                    int mi_row, int mi_col) {
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const MB_MODE_INFO_EXT *const mbmi_ext = x->mbmi_ext;
  const CurrentFrame *const current_frame = &cm->current_frame;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  FRAME_CONTEXT *fc = xd->tile_ctx;
  const int seg_ref_active =
      segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_REF_FRAME);

  if (current_frame->skip_mode_info.skip_mode_flag && !seg_ref_active &&
      is_comp_ref_allowed(bsize)) {
    const int skip_mode_ctx = av1_get_skip_mode_context(xd);
#if CONFIG_ENTROPY_STATS
    td->counts->skip_mode[skip_mode_ctx][mbmi->skip_mode]++;
#endif
    update_cdf(fc->skip_mode_cdfs[skip_mode_ctx], mbmi->skip_mode, 2);
  }

  if (!mbmi->skip_mode && !seg_ref_active) {
    const int skip_ctx = av1_get_skip_context(xd);
#if CONFIG_ENTROPY_STATS
    td->counts->skip[skip_ctx][mbmi->skip]++;
#endif
    update_cdf(fc->skip_cdfs[skip_ctx], mbmi->skip, 2);
  }

#if CONFIG_ENTROPY_STATS
  // delta quant applies to both intra and inter
  const int super_block_upper_left =
      ((mi_row & (cm->seq_params.mib_size - 1)) == 0) &&
      ((mi_col & (cm->seq_params.mib_size - 1)) == 0);
  const DeltaQInfo *const delta_q_info = &cm->delta_q_info;
  if (delta_q_info->delta_q_present_flag &&
      (bsize != cm->seq_params.sb_size || !mbmi->skip) &&
      super_block_upper_left) {
    const int dq =
        (mbmi->current_qindex - xd->current_qindex) / delta_q_info->delta_q_res;
    const int absdq = abs(dq);
    for (int i = 0; i < AOMMIN(absdq, DELTA_Q_SMALL); ++i) {
      td->counts->delta_q[i][1]++;
    }
    if (absdq < DELTA_Q_SMALL) td->counts->delta_q[absdq][0]++;
    if (delta_q_info->delta_lf_present_flag) {
      if (delta_q_info->delta_lf_multi) {
        const int frame_lf_count =
            av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
        for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
          const int delta_lf = (mbmi->delta_lf[lf_id] - xd->delta_lf[lf_id]) /
                               delta_q_info->delta_lf_res;
          const int abs_delta_lf = abs(delta_lf);
          for (int i = 0; i < AOMMIN(abs_delta_lf, DELTA_LF_SMALL); ++i) {
            td->counts->delta_lf_multi[lf_id][i][1]++;
          }
          if (abs_delta_lf < DELTA_LF_SMALL)
            td->counts->delta_lf_multi[lf_id][abs_delta_lf][0]++;
        }
      } else {
        const int delta_lf =
            (mbmi->delta_lf_from_base - xd->delta_lf_from_base) /
            delta_q_info->delta_lf_res;
        const int abs_delta_lf = abs(delta_lf);
        for (int i = 0; i < AOMMIN(abs_delta_lf, DELTA_LF_SMALL); ++i) {
          td->counts->delta_lf[i][1]++;
        }
        if (abs_delta_lf < DELTA_LF_SMALL)
          td->counts->delta_lf[abs_delta_lf][0]++;
      }
    }
  }
#endif

  if (!is_inter_block(mbmi)) {
    sum_intra_stats(cm, td->counts, xd, mbmi, xd->above_mbmi, xd->left_mbmi,
                    frame_is_intra_only(cm), mi_row, mi_col);
  }

  if (av1_allow_intrabc(cm)) {
    update_cdf(fc->intrabc_cdf, is_intrabc_block(mbmi), 2);
#if CONFIG_ENTROPY_STATS
    ++td->counts->intrabc[is_intrabc_block(mbmi)];
#endif  // CONFIG_ENTROPY_STATS
  }

  if (frame_is_intra_only(cm) || mbmi->skip_mode) return;

  FRAME_COUNTS *const counts = td->counts;
  const int inter_block = is_inter_block(mbmi);

  if (!seg_ref_active) {
#if CONFIG_ENTROPY_STATS
    counts->intra_inter[av1_get_intra_inter_context(xd)][inter_block]++;
#endif
    update_cdf(fc->intra_inter_cdf[av1_get_intra_inter_context(xd)],
               inter_block, 2);
    // If the segment reference feature is enabled we have only a single
    // reference frame allowed for the segment so exclude it from
    // the reference frame counts used to work out probabilities.
    if (inter_block) {
      const MV_REFERENCE_FRAME ref0 = mbmi->ref_frame[0];
      const MV_REFERENCE_FRAME ref1 = mbmi->ref_frame[1];
      if (current_frame->reference_mode == REFERENCE_MODE_SELECT) {
        if (is_comp_ref_allowed(bsize)) {
#if CONFIG_ENTROPY_STATS
          counts->comp_inter[av1_get_reference_mode_context(xd)]
                            [has_second_ref(mbmi)]++;
#endif  // CONFIG_ENTROPY_STATS
          update_cdf(av1_get_reference_mode_cdf(xd), has_second_ref(mbmi), 2);
        }
      }

      if (has_second_ref(mbmi)) {
        const COMP_REFERENCE_TYPE comp_ref_type = has_uni_comp_refs(mbmi)
                                                      ? UNIDIR_COMP_REFERENCE
                                                      : BIDIR_COMP_REFERENCE;
        update_cdf(av1_get_comp_reference_type_cdf(xd), comp_ref_type,
                   COMP_REFERENCE_TYPES);
#if CONFIG_ENTROPY_STATS
        counts->comp_ref_type[av1_get_comp_reference_type_context(xd)]
                             [comp_ref_type]++;
#endif  // CONFIG_ENTROPY_STATS

        if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
          const int bit = (ref0 == BWDREF_FRAME);
          update_cdf(av1_get_pred_cdf_uni_comp_ref_p(xd), bit, 2);
#if CONFIG_ENTROPY_STATS
          counts
              ->uni_comp_ref[av1_get_pred_context_uni_comp_ref_p(xd)][0][bit]++;
#endif  // CONFIG_ENTROPY_STATS
          if (!bit) {
            const int bit1 = (ref1 == LAST3_FRAME || ref1 == GOLDEN_FRAME);
            update_cdf(av1_get_pred_cdf_uni_comp_ref_p1(xd), bit1, 2);
#if CONFIG_ENTROPY_STATS
            counts->uni_comp_ref[av1_get_pred_context_uni_comp_ref_p1(xd)][1]
                                [bit1]++;
#endif  // CONFIG_ENTROPY_STATS
            if (bit1) {
              update_cdf(av1_get_pred_cdf_uni_comp_ref_p2(xd),
                         ref1 == GOLDEN_FRAME, 2);
#if CONFIG_ENTROPY_STATS
              counts->uni_comp_ref[av1_get_pred_context_uni_comp_ref_p2(xd)][2]
                                  [ref1 == GOLDEN_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
            }
          }
        } else {
          const int bit = (ref0 == GOLDEN_FRAME || ref0 == LAST3_FRAME);
          update_cdf(av1_get_pred_cdf_comp_ref_p(xd), bit, 2);
#if CONFIG_ENTROPY_STATS
          counts->comp_ref[av1_get_pred_context_comp_ref_p(xd)][0][bit]++;
#endif  // CONFIG_ENTROPY_STATS
          if (!bit) {
            update_cdf(av1_get_pred_cdf_comp_ref_p1(xd), ref0 == LAST2_FRAME,
                       2);
#if CONFIG_ENTROPY_STATS
            counts->comp_ref[av1_get_pred_context_comp_ref_p1(xd)][1]
                            [ref0 == LAST2_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          } else {
            update_cdf(av1_get_pred_cdf_comp_ref_p2(xd), ref0 == GOLDEN_FRAME,
                       2);
#if CONFIG_ENTROPY_STATS
            counts->comp_ref[av1_get_pred_context_comp_ref_p2(xd)][2]
                            [ref0 == GOLDEN_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          }
          update_cdf(av1_get_pred_cdf_comp_bwdref_p(xd), ref1 == ALTREF_FRAME,
                     2);
#if CONFIG_ENTROPY_STATS
          counts->comp_bwdref[av1_get_pred_context_comp_bwdref_p(xd)][0]
                             [ref1 == ALTREF_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          if (ref1 != ALTREF_FRAME) {
            update_cdf(av1_get_pred_cdf_comp_bwdref_p1(xd),
                       ref1 == ALTREF2_FRAME, 2);
#if CONFIG_ENTROPY_STATS
            counts->comp_bwdref[av1_get_pred_context_comp_bwdref_p1(xd)][1]
                               [ref1 == ALTREF2_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          }
        }
      } else {
        const int bit = (ref0 >= BWDREF_FRAME);
        update_cdf(av1_get_pred_cdf_single_ref_p1(xd), bit, 2);
#if CONFIG_ENTROPY_STATS
        counts->single_ref[av1_get_pred_context_single_ref_p1(xd)][0][bit]++;
#endif  // CONFIG_ENTROPY_STATS
        if (bit) {
          assert(ref0 <= ALTREF_FRAME);
          update_cdf(av1_get_pred_cdf_single_ref_p2(xd), ref0 == ALTREF_FRAME,
                     2);
#if CONFIG_ENTROPY_STATS
          counts->single_ref[av1_get_pred_context_single_ref_p2(xd)][1]
                            [ref0 == ALTREF_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          if (ref0 != ALTREF_FRAME) {
            update_cdf(av1_get_pred_cdf_single_ref_p6(xd),
                       ref0 == ALTREF2_FRAME, 2);
#if CONFIG_ENTROPY_STATS
            counts->single_ref[av1_get_pred_context_single_ref_p6(xd)][5]
                              [ref0 == ALTREF2_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          }
        } else {
          const int bit1 = !(ref0 == LAST2_FRAME || ref0 == LAST_FRAME);
          update_cdf(av1_get_pred_cdf_single_ref_p3(xd), bit1, 2);
#if CONFIG_ENTROPY_STATS
          counts->single_ref[av1_get_pred_context_single_ref_p3(xd)][2][bit1]++;
#endif  // CONFIG_ENTROPY_STATS
          if (!bit1) {
            update_cdf(av1_get_pred_cdf_single_ref_p4(xd), ref0 != LAST_FRAME,
                       2);
#if CONFIG_ENTROPY_STATS
            counts->single_ref[av1_get_pred_context_single_ref_p4(xd)][3]
                              [ref0 != LAST_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          } else {
            update_cdf(av1_get_pred_cdf_single_ref_p5(xd), ref0 != LAST3_FRAME,
                       2);
#if CONFIG_ENTROPY_STATS
            counts->single_ref[av1_get_pred_context_single_ref_p5(xd)][4]
                              [ref0 != LAST3_FRAME]++;
#endif  // CONFIG_ENTROPY_STATS
          }
        }
      }

      if (cm->seq_params.enable_interintra_compound &&
          is_interintra_allowed(mbmi)) {
        const int bsize_group = size_group_lookup[bsize];
        if (mbmi->ref_frame[1] == INTRA_FRAME) {
#if CONFIG_ENTROPY_STATS
          counts->interintra[bsize_group][1]++;
#endif
          update_cdf(fc->interintra_cdf[bsize_group], 1, 2);
#if CONFIG_ENTROPY_STATS
          counts->interintra_mode[bsize_group][mbmi->interintra_mode]++;
#endif
          update_cdf(fc->interintra_mode_cdf[bsize_group],
                     mbmi->interintra_mode, INTERINTRA_MODES);
          if (is_interintra_wedge_used(bsize)) {
#if CONFIG_ENTROPY_STATS
            counts->wedge_interintra[bsize][mbmi->use_wedge_interintra]++;
#endif
            update_cdf(fc->wedge_interintra_cdf[bsize],
                       mbmi->use_wedge_interintra, 2);
            if (mbmi->use_wedge_interintra) {
#if CONFIG_ENTROPY_STATS
              counts->wedge_idx[bsize][mbmi->interintra_wedge_index]++;
#endif
              update_cdf(fc->wedge_idx_cdf[bsize], mbmi->interintra_wedge_index,
                         16);
            }
          }
        } else {
#if CONFIG_ENTROPY_STATS
          counts->interintra[bsize_group][0]++;
#endif
          update_cdf(fc->interintra_cdf[bsize_group], 0, 2);
        }
      }

      const MOTION_MODE motion_allowed =
          cm->switchable_motion_mode
              ? motion_mode_allowed(xd->global_motion, xd, mbmi,
                                    cm->allow_warped_motion)
              : SIMPLE_TRANSLATION;
      if (mbmi->ref_frame[1] != INTRA_FRAME) {
        if (motion_allowed == WARPED_CAUSAL) {
#if CONFIG_ENTROPY_STATS
          counts->motion_mode[bsize][mbmi->motion_mode]++;
#endif
          update_cdf(fc->motion_mode_cdf[bsize], mbmi->motion_mode,
                     MOTION_MODES);
        } else if (motion_allowed == OBMC_CAUSAL) {
#if CONFIG_ENTROPY_STATS
          counts->obmc[bsize][mbmi->motion_mode == OBMC_CAUSAL]++;
#endif
          update_cdf(fc->obmc_cdf[bsize], mbmi->motion_mode == OBMC_CAUSAL, 2);
        }
      }

      if (has_second_ref(mbmi)) {
        assert(current_frame->reference_mode != SINGLE_REFERENCE &&
               is_inter_compound_mode(mbmi->mode) &&
               mbmi->motion_mode == SIMPLE_TRANSLATION);

        const int masked_compound_used = is_any_masked_compound_used(bsize) &&
                                         cm->seq_params.enable_masked_compound;
        if (masked_compound_used) {
          const int comp_group_idx_ctx = get_comp_group_idx_context(xd);
#if CONFIG_ENTROPY_STATS
          ++counts->comp_group_idx[comp_group_idx_ctx][mbmi->comp_group_idx];
#endif
          update_cdf(fc->comp_group_idx_cdf[comp_group_idx_ctx],
                     mbmi->comp_group_idx, 2);
        }

        if (mbmi->comp_group_idx == 0) {
          const int comp_index_ctx = get_comp_index_context(cm, xd);
#if CONFIG_ENTROPY_STATS
          ++counts->compound_index[comp_index_ctx][mbmi->compound_idx];
#endif
          update_cdf(fc->compound_index_cdf[comp_index_ctx], mbmi->compound_idx,
                     2);
        } else {
          assert(masked_compound_used);
          if (is_interinter_compound_used(COMPOUND_WEDGE, bsize)) {
#if CONFIG_ENTROPY_STATS
            ++counts->compound_type[bsize][mbmi->interinter_comp.type -
                                           COMPOUND_WEDGE];
#endif
            update_cdf(fc->compound_type_cdf[bsize],
                       mbmi->interinter_comp.type - COMPOUND_WEDGE,
                       MASKED_COMPOUND_TYPES);
          }
        }
      }
      if (mbmi->interinter_comp.type == COMPOUND_WEDGE) {
        if (is_interinter_compound_used(COMPOUND_WEDGE, bsize)) {
#if CONFIG_ENTROPY_STATS
          counts->wedge_idx[bsize][mbmi->interinter_comp.wedge_index]++;
#endif
          update_cdf(fc->wedge_idx_cdf[bsize],
                     mbmi->interinter_comp.wedge_index, 16);
        }
      }
    }
  }

  if (inter_block && cm->interp_filter == SWITCHABLE &&
      mbmi->motion_mode != WARPED_CAUSAL &&
      !is_nontrans_global_motion(xd, mbmi)) {
    update_filter_type_cdf(xd, mbmi);
  }
  if (inter_block &&
      !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP)) {
    const PREDICTION_MODE mode = mbmi->mode;
    const int16_t mode_ctx =
        av1_mode_context_analyzer(mbmi_ext->mode_context, mbmi->ref_frame);
    if (has_second_ref(mbmi)) {
#if CONFIG_ENTROPY_STATS
      ++counts->inter_compound_mode[mode_ctx][INTER_COMPOUND_OFFSET(mode)];
#endif
      update_cdf(fc->inter_compound_mode_cdf[mode_ctx],
                 INTER_COMPOUND_OFFSET(mode), INTER_COMPOUND_MODES);
    } else {
      update_inter_mode_stats(fc, counts, mode, mode_ctx);
    }

    const int new_mv = mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV;
    if (new_mv) {
      const uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
      for (int idx = 0; idx < 2; ++idx) {
        if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
          const uint8_t drl_ctx =
              av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);
          update_cdf(fc->drl_cdf[drl_ctx], mbmi->ref_mv_idx != idx, 2);
#if CONFIG_ENTROPY_STATS
          ++counts->drl_mode[drl_ctx][mbmi->ref_mv_idx != idx];
#endif
          if (mbmi->ref_mv_idx == idx) break;
        }
      }
    }

    if (have_nearmv_in_inter_mode(mbmi->mode)) {
      const uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
      for (int idx = 1; idx < 3; ++idx) {
        if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
          const uint8_t drl_ctx =
              av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);
          update_cdf(fc->drl_cdf[drl_ctx], mbmi->ref_mv_idx != idx - 1, 2);
#if CONFIG_ENTROPY_STATS
          ++counts->drl_mode[drl_ctx][mbmi->ref_mv_idx != idx - 1];
#endif
          if (mbmi->ref_mv_idx == idx - 1) break;
        }
      }
    }
    if (have_newmv_in_inter_mode(mbmi->mode)) {
      const int allow_hp = cm->cur_frame_force_integer_mv
                               ? MV_SUBPEL_NONE
                               : cm->allow_high_precision_mv;
      if (new_mv) {
        for (int ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
          const int_mv ref_mv = av1_get_ref_mv(x, ref);
          av1_update_mv_stats(&mbmi->mv[ref].as_mv, &ref_mv.as_mv, &fc->nmvc,
                              allow_hp);
        }
      } else if (mbmi->mode == NEAREST_NEWMV || mbmi->mode == NEAR_NEWMV) {
        const int ref = 1;
        const int_mv ref_mv = av1_get_ref_mv(x, ref);
        av1_update_mv_stats(&mbmi->mv[ref].as_mv, &ref_mv.as_mv, &fc->nmvc,
                            allow_hp);
      } else if (mbmi->mode == NEW_NEARESTMV || mbmi->mode == NEW_NEARMV) {
        const int ref = 0;
        const int_mv ref_mv = av1_get_ref_mv(x, ref);
        av1_update_mv_stats(&mbmi->mv[ref].as_mv, &ref_mv.as_mv, &fc->nmvc,
                            allow_hp);
      }
    }
  }
}

typedef struct {
  ENTROPY_CONTEXT a[MAX_MIB_SIZE * MAX_MB_PLANE];
  ENTROPY_CONTEXT l[MAX_MIB_SIZE * MAX_MB_PLANE];
  PARTITION_CONTEXT sa[MAX_MIB_SIZE];
  PARTITION_CONTEXT sl[MAX_MIB_SIZE];
  TXFM_CONTEXT *p_ta;
  TXFM_CONTEXT *p_tl;
  TXFM_CONTEXT ta[MAX_MIB_SIZE];
  TXFM_CONTEXT tl[MAX_MIB_SIZE];
} RD_SEARCH_MACROBLOCK_CONTEXT;

static AOM_INLINE void restore_context(MACROBLOCK *x,
                                       const RD_SEARCH_MACROBLOCK_CONTEXT *ctx,
                                       int mi_row, int mi_col, BLOCK_SIZE bsize,
                                       const int num_planes) {
  MACROBLOCKD *xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide =
      block_size_wide[bsize] >> tx_size_wide_log2[0];
  const int num_4x4_blocks_high =
      block_size_high[bsize] >> tx_size_high_log2[0];
  int mi_width = mi_size_wide[bsize];
  int mi_height = mi_size_high[bsize];
  for (p = 0; p < num_planes; p++) {
    int tx_col = mi_col;
    int tx_row = mi_row & MAX_MIB_MASK;
    memcpy(xd->above_context[p] + (tx_col >> xd->plane[p].subsampling_x),
           ctx->a + num_4x4_blocks_wide * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(xd->left_context[p] + (tx_row >> xd->plane[p].subsampling_y),
           ctx->l + num_4x4_blocks_high * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(xd->above_seg_context + mi_col, ctx->sa,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(xd->left_seg_context + (mi_row & MAX_MIB_MASK), ctx->sl,
         sizeof(xd->left_seg_context[0]) * mi_height);
  xd->above_txfm_context = ctx->p_ta;
  xd->left_txfm_context = ctx->p_tl;
  memcpy(xd->above_txfm_context, ctx->ta,
         sizeof(*xd->above_txfm_context) * mi_width);
  memcpy(xd->left_txfm_context, ctx->tl,
         sizeof(*xd->left_txfm_context) * mi_height);
}

static AOM_INLINE void save_context(const MACROBLOCK *x,
                                    RD_SEARCH_MACROBLOCK_CONTEXT *ctx,
                                    int mi_row, int mi_col, BLOCK_SIZE bsize,
                                    const int num_planes) {
  const MACROBLOCKD *xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide =
      block_size_wide[bsize] >> tx_size_wide_log2[0];
  const int num_4x4_blocks_high =
      block_size_high[bsize] >> tx_size_high_log2[0];
  int mi_width = mi_size_wide[bsize];
  int mi_height = mi_size_high[bsize];

  // buffer the above/left context information of the block in search.
  for (p = 0; p < num_planes; ++p) {
    int tx_col = mi_col;
    int tx_row = mi_row & MAX_MIB_MASK;
    memcpy(ctx->a + num_4x4_blocks_wide * p,
           xd->above_context[p] + (tx_col >> xd->plane[p].subsampling_x),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(ctx->l + num_4x4_blocks_high * p,
           xd->left_context[p] + (tx_row >> xd->plane[p].subsampling_y),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(ctx->sa, xd->above_seg_context + mi_col,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(ctx->sl, xd->left_seg_context + (mi_row & MAX_MIB_MASK),
         sizeof(xd->left_seg_context[0]) * mi_height);
  memcpy(ctx->ta, xd->above_txfm_context,
         sizeof(*xd->above_txfm_context) * mi_width);
  memcpy(ctx->tl, xd->left_txfm_context,
         sizeof(*xd->left_txfm_context) * mi_height);
  ctx->p_ta = xd->above_txfm_context;
  ctx->p_tl = xd->left_txfm_context;
}

static AOM_INLINE void encode_b(const AV1_COMP *const cpi,
                                TileDataEnc *tile_data, ThreadData *td,
                                TOKENEXTRA **tp, int mi_row, int mi_col,
                                RUN_TYPE dry_run, BLOCK_SIZE bsize,
                                PARTITION_TYPE partition,
                                PICK_MODE_CONTEXT *const ctx, int *rate) {
  TileInfo *const tile = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;

  set_offsets_without_segment_id(cpi, tile, x, mi_row, mi_col, bsize);
  const int origin_mult = x->rdmult;
  setup_block_rdmult(cpi, x, mi_row, mi_col, bsize, NO_AQ, NULL);
  MB_MODE_INFO *mbmi = xd->mi[0];
  mbmi->partition = partition;
  update_state(cpi, td, ctx, mi_row, mi_col, bsize, dry_run);

  if (!dry_run) {
    x->mbmi_ext_frame->cb_offset = x->cb_offset;
    assert(x->cb_offset <
           (1 << num_pels_log2_lookup[cpi->common.seq_params.sb_size]));
  }

  encode_superblock(cpi, tile_data, td, tp, dry_run, bsize, rate);

  if (!dry_run) {
    const AV1_COMMON *const cm = &cpi->common;
    x->cb_offset += block_size_wide[bsize] * block_size_high[bsize];
    if (bsize == cpi->common.seq_params.sb_size && mbmi->skip == 1 &&
        cm->delta_q_info.delta_lf_present_flag) {
      const int frame_lf_count =
          av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
      for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id)
        mbmi->delta_lf[lf_id] = xd->delta_lf[lf_id];
      mbmi->delta_lf_from_base = xd->delta_lf_from_base;
    }
    if (has_second_ref(mbmi)) {
      if (mbmi->compound_idx == 0 ||
          mbmi->interinter_comp.type == COMPOUND_AVERAGE)
        mbmi->comp_group_idx = 0;
      else
        mbmi->comp_group_idx = 1;
    }

    // delta quant applies to both intra and inter
    const int super_block_upper_left =
        ((mi_row & (cm->seq_params.mib_size - 1)) == 0) &&
        ((mi_col & (cm->seq_params.mib_size - 1)) == 0);
    const DeltaQInfo *const delta_q_info = &cm->delta_q_info;
    if (delta_q_info->delta_q_present_flag &&
        (bsize != cm->seq_params.sb_size || !mbmi->skip) &&
        super_block_upper_left) {
      xd->current_qindex = mbmi->current_qindex;
      if (delta_q_info->delta_lf_present_flag) {
        if (delta_q_info->delta_lf_multi) {
          const int frame_lf_count =
              av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
          for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
            xd->delta_lf[lf_id] = mbmi->delta_lf[lf_id];
          }
        } else {
          xd->delta_lf_from_base = mbmi->delta_lf_from_base;
        }
      }
    }

    RD_COUNTS *rdc = &td->rd_counts;
    if (mbmi->skip_mode) {
      assert(!frame_is_intra_only(cm));
      rdc->skip_mode_used_flag = 1;
      if (cm->current_frame.reference_mode == REFERENCE_MODE_SELECT) {
        assert(has_second_ref(mbmi));
        rdc->compound_ref_used_flag = 1;
      }
      set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);
    } else {
      const int seg_ref_active =
          segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_REF_FRAME);
      if (!seg_ref_active) {
        // If the segment reference feature is enabled we have only a single
        // reference frame allowed for the segment so exclude it from
        // the reference frame counts used to work out probabilities.
        if (is_inter_block(mbmi)) {
          av1_collect_neighbors_ref_counts(xd);
          if (cm->current_frame.reference_mode == REFERENCE_MODE_SELECT) {
            if (has_second_ref(mbmi)) {
              // This flag is also updated for 4x4 blocks
              rdc->compound_ref_used_flag = 1;
            }
          }
          set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);
        }
      }
    }

    if (tile_data->allow_update_cdf) {
      update_stats(&cpi->common, td, mi_row, mi_col);
    }

    // Gather obmc count to update the probability.
    if (cpi->sf.prune_obmc_prob_thresh > 0) {
      const int inter_block = is_inter_block(mbmi);
      const int seg_ref_active =
          segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_REF_FRAME);
      if (!seg_ref_active && inter_block) {
        const MOTION_MODE motion_allowed =
            cm->switchable_motion_mode
                ? motion_mode_allowed(xd->global_motion, xd, mbmi,
                                      cm->allow_warped_motion)
                : SIMPLE_TRANSLATION;
        if (mbmi->ref_frame[1] != INTRA_FRAME &&
            motion_allowed == OBMC_CAUSAL) {
          td->rd_counts.obmc_used[bsize][mbmi->motion_mode == OBMC_CAUSAL]++;
        }
      }
    }
  }
  // TODO(Ravi/Remya): Move this copy function to a better logical place
  copy_winner_ref_mode_from_mbmi_ext(x);
  x->rdmult = origin_mult;
}

static AOM_INLINE void encode_sb(const AV1_COMP *const cpi, ThreadData *td,
                                 TileDataEnc *tile_data, TOKENEXTRA **tp,
                                 int mi_row, int mi_col, RUN_TYPE dry_run,
                                 BLOCK_SIZE bsize, PC_TREE *pc_tree,
                                 int *rate) {
  assert(bsize < BLOCK_SIZES_ALL);
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  assert(bsize < BLOCK_SIZES_ALL);
  const int hbs = mi_size_wide[bsize] / 2;
  const int is_partition_root = bsize >= BLOCK_8X8;
  const int ctx = is_partition_root
                      ? partition_plane_context(xd, mi_row, mi_col, bsize)
                      : -1;
  const PARTITION_TYPE partition = pc_tree->partitioning;
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, partition);
  int quarter_step = mi_size_wide[bsize] / 4;
  int i;
  BLOCK_SIZE bsize2 = get_partition_subsize(bsize, PARTITION_SPLIT);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  if (!dry_run && ctx >= 0) {
    const int has_rows = (mi_row + hbs) < cm->mi_rows;
    const int has_cols = (mi_col + hbs) < cm->mi_cols;

    if (has_rows && has_cols) {
#if CONFIG_ENTROPY_STATS
      td->counts->partition[ctx][partition]++;
#endif

      if (tile_data->allow_update_cdf) {
        FRAME_CONTEXT *fc = xd->tile_ctx;
        update_cdf(fc->partition_cdf[ctx], partition,
                   partition_cdf_length(bsize));
      }
    }
  }

  switch (partition) {
    case PARTITION_NONE:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, subsize,
               partition, &pc_tree->none, rate);
      break;
    case PARTITION_VERT:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, subsize,
               partition, &pc_tree->vertical[0], rate);
      if (mi_col + hbs < cm->mi_cols) {
        encode_b(cpi, tile_data, td, tp, mi_row, mi_col + hbs, dry_run, subsize,
                 partition, &pc_tree->vertical[1], rate);
      }
      break;
    case PARTITION_HORZ:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, subsize,
               partition, &pc_tree->horizontal[0], rate);
      if (mi_row + hbs < cm->mi_rows) {
        encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col, dry_run, subsize,
                 partition, &pc_tree->horizontal[1], rate);
      }
      break;
    case PARTITION_SPLIT:
      encode_sb(cpi, td, tile_data, tp, mi_row, mi_col, dry_run, subsize,
                pc_tree->split[0], rate);
      encode_sb(cpi, td, tile_data, tp, mi_row, mi_col + hbs, dry_run, subsize,
                pc_tree->split[1], rate);
      encode_sb(cpi, td, tile_data, tp, mi_row + hbs, mi_col, dry_run, subsize,
                pc_tree->split[2], rate);
      encode_sb(cpi, td, tile_data, tp, mi_row + hbs, mi_col + hbs, dry_run,
                subsize, pc_tree->split[3], rate);
      break;

    case PARTITION_HORZ_A:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, bsize2,
               partition, &pc_tree->horizontala[0], rate);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col + hbs, dry_run, bsize2,
               partition, &pc_tree->horizontala[1], rate);
      encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col, dry_run, subsize,
               partition, &pc_tree->horizontala[2], rate);
      break;
    case PARTITION_HORZ_B:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, subsize,
               partition, &pc_tree->horizontalb[0], rate);
      encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col, dry_run, bsize2,
               partition, &pc_tree->horizontalb[1], rate);
      encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col + hbs, dry_run,
               bsize2, partition, &pc_tree->horizontalb[2], rate);
      break;
    case PARTITION_VERT_A:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, bsize2,
               partition, &pc_tree->verticala[0], rate);
      encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col, dry_run, bsize2,
               partition, &pc_tree->verticala[1], rate);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col + hbs, dry_run, subsize,
               partition, &pc_tree->verticala[2], rate);

      break;
    case PARTITION_VERT_B:
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, dry_run, subsize,
               partition, &pc_tree->verticalb[0], rate);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col + hbs, dry_run, bsize2,
               partition, &pc_tree->verticalb[1], rate);
      encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col + hbs, dry_run,
               bsize2, partition, &pc_tree->verticalb[2], rate);
      break;
    case PARTITION_HORZ_4:
      for (i = 0; i < 4; ++i) {
        int this_mi_row = mi_row + i * quarter_step;
        if (i > 0 && this_mi_row >= cm->mi_rows) break;

        encode_b(cpi, tile_data, td, tp, this_mi_row, mi_col, dry_run, subsize,
                 partition, &pc_tree->horizontal4[i], rate);
      }
      break;
    case PARTITION_VERT_4:
      for (i = 0; i < 4; ++i) {
        int this_mi_col = mi_col + i * quarter_step;
        if (i > 0 && this_mi_col >= cm->mi_cols) break;
        encode_b(cpi, tile_data, td, tp, mi_row, this_mi_col, dry_run, subsize,
                 partition, &pc_tree->vertical4[i], rate);
      }
      break;
    default: assert(0 && "Invalid partition type."); break;
  }

  update_ext_partition_context(xd, mi_row, mi_col, subsize, bsize, partition);
}

static AOM_INLINE void set_partial_sb_partition(
    const AV1_COMMON *const cm, MB_MODE_INFO *mi, int bh_in, int bw_in,
    int mi_rows_remaining, int mi_cols_remaining, BLOCK_SIZE bsize,
    MB_MODE_INFO **mib) {
  int bh = bh_in;
  int r, c;
  for (r = 0; r < cm->seq_params.mib_size; r += bh) {
    int bw = bw_in;
    for (c = 0; c < cm->seq_params.mib_size; c += bw) {
      const int index = r * cm->mi_stride + c;
      mib[index] = mi + index;
      mib[index]->sb_type = find_partition_size(
          bsize, mi_rows_remaining - r, mi_cols_remaining - c, &bh, &bw);
    }
  }
}

// This function attempts to set all mode info entries in a given superblock
// to the same block partition size.
// However, at the bottom and right borders of the image the requested size
// may not be allowed in which case this code attempts to choose the largest
// allowable partition.
static AOM_INLINE void set_fixed_partitioning(AV1_COMP *cpi,
                                              const TileInfo *const tile,
                                              MB_MODE_INFO **mib, int mi_row,
                                              int mi_col, BLOCK_SIZE bsize) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_rows_remaining = tile->mi_row_end - mi_row;
  const int mi_cols_remaining = tile->mi_col_end - mi_col;
  int block_row, block_col;
  MB_MODE_INFO *const mi_upper_left =
      cm->mi + get_alloc_mi_idx(cm, mi_row, mi_col);
  int bh = mi_size_high[bsize];
  int bw = mi_size_wide[bsize];

  assert(bsize >= cm->mi_alloc_bsize &&
         "Attempted to use bsize < cm->mi_alloc_bsize");
  assert((mi_rows_remaining > 0) && (mi_cols_remaining > 0));

  // Apply the requested partition size to the SB if it is all "in image"
  if ((mi_cols_remaining >= cm->seq_params.mib_size) &&
      (mi_rows_remaining >= cm->seq_params.mib_size)) {
    for (block_row = 0; block_row < cm->seq_params.mib_size; block_row += bh) {
      for (block_col = 0; block_col < cm->seq_params.mib_size;
           block_col += bw) {
        const int grid_index = get_mi_grid_idx(cm, block_row, block_col);
        const int mi_index = get_alloc_mi_idx(cm, block_row, block_col);
        mib[grid_index] = mi_upper_left + mi_index;
        mib[grid_index]->sb_type = bsize;
      }
    }
  } else {
    // Else this is a partial SB.
    set_partial_sb_partition(cm, mi_upper_left, bh, bw, mi_rows_remaining,
                             mi_cols_remaining, bsize, mib);
  }
}

static AOM_INLINE void rd_use_partition(
    AV1_COMP *cpi, ThreadData *td, TileDataEnc *tile_data, MB_MODE_INFO **mib,
    TOKENEXTRA **tp, int mi_row, int mi_col, BLOCK_SIZE bsize, int *rate,
    int64_t *dist, int do_recon, PC_TREE *pc_tree) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int bs = mi_size_wide[bsize];
  const int hbs = bs / 2;
  int i;
  const int pl = (bsize >= BLOCK_8X8)
                     ? partition_plane_context(xd, mi_row, mi_col, bsize)
                     : 0;
  const PARTITION_TYPE partition =
      (bsize >= BLOCK_8X8) ? get_partition(cm, mi_row, mi_col, bsize)
                           : PARTITION_NONE;
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, partition);
  RD_SEARCH_MACROBLOCK_CONTEXT x_ctx;
  RD_STATS last_part_rdc, none_rdc, chosen_rdc, invalid_rdc;
  BLOCK_SIZE sub_subsize = BLOCK_4X4;
  int splits_below = 0;
  BLOCK_SIZE bs_type = mib[0]->sb_type;
  int do_partition_search = 1;
  PICK_MODE_CONTEXT *ctx_none = &pc_tree->none;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  assert(mi_size_wide[bsize] == mi_size_high[bsize]);

  av1_invalid_rd_stats(&last_part_rdc);
  av1_invalid_rd_stats(&none_rdc);
  av1_invalid_rd_stats(&chosen_rdc);
  av1_invalid_rd_stats(&invalid_rdc);

  pc_tree->partitioning = partition;

  xd->above_txfm_context = cm->above_txfm_context[tile_info->tile_row] + mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);
  save_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);

  if (bsize == BLOCK_16X16 && cpi->vaq_refresh) {
    set_offsets(cpi, tile_info, x, mi_row, mi_col, bsize);
    x->mb_energy = av1_log_block_var(cpi, x, bsize);
  }

  // Save rdmult before it might be changed, so it can be restored later.
  const int orig_rdmult = x->rdmult;
  setup_block_rdmult(cpi, x, mi_row, mi_col, bsize, NO_AQ, NULL);

  if (do_partition_search &&
      cpi->sf.partition_search_type == SEARCH_PARTITION &&
      cpi->sf.adjust_partitioning_from_last_frame) {
    // Check if any of the sub blocks are further split.
    if (partition == PARTITION_SPLIT && subsize > BLOCK_8X8) {
      sub_subsize = get_partition_subsize(subsize, PARTITION_SPLIT);
      splits_below = 1;
      for (i = 0; i < 4; i++) {
        int jj = i >> 1, ii = i & 0x01;
        MB_MODE_INFO *this_mi = mib[jj * hbs * cm->mi_stride + ii * hbs];
        if (this_mi && this_mi->sb_type >= sub_subsize) {
          splits_below = 0;
        }
      }
    }

    // If partition is not none try none unless each of the 4 splits are split
    // even further..
    if (partition != PARTITION_NONE && !splits_below &&
        mi_row + hbs < cm->mi_rows && mi_col + hbs < cm->mi_cols) {
      pc_tree->partitioning = PARTITION_NONE;
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &none_rdc,
                    PARTITION_NONE, bsize, ctx_none, invalid_rdc, PICK_MODE_RD);

      if (none_rdc.rate < INT_MAX) {
        none_rdc.rate += x->partition_cost[pl][PARTITION_NONE];
        none_rdc.rdcost = RDCOST(x->rdmult, none_rdc.rate, none_rdc.dist);
      }

      restore_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);
      mib[0]->sb_type = bs_type;
      pc_tree->partitioning = partition;
    }
  }

  switch (partition) {
    case PARTITION_NONE:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                    PARTITION_NONE, bsize, ctx_none, invalid_rdc, PICK_MODE_RD);
      break;
    case PARTITION_HORZ:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                    PARTITION_HORZ, subsize, &pc_tree->horizontal[0],
                    invalid_rdc, PICK_MODE_RD);
      if (last_part_rdc.rate != INT_MAX && bsize >= BLOCK_8X8 &&
          mi_row + hbs < cm->mi_rows) {
        RD_STATS tmp_rdc;
        const PICK_MODE_CONTEXT *const ctx_h = &pc_tree->horizontal[0];
        av1_init_rd_stats(&tmp_rdc);
        update_state(cpi, td, ctx_h, mi_row, mi_col, subsize, 1);
        encode_superblock(cpi, tile_data, td, tp, DRY_RUN_NORMAL, subsize,
                          NULL);
        pick_sb_modes(cpi, tile_data, x, mi_row + hbs, mi_col, &tmp_rdc,
                      PARTITION_HORZ, subsize, &pc_tree->horizontal[1],
                      invalid_rdc, PICK_MODE_RD);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          av1_invalid_rd_stats(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
      }
      break;
    case PARTITION_VERT:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &last_part_rdc,
                    PARTITION_VERT, subsize, &pc_tree->vertical[0], invalid_rdc,
                    PICK_MODE_RD);
      if (last_part_rdc.rate != INT_MAX && bsize >= BLOCK_8X8 &&
          mi_col + hbs < cm->mi_cols) {
        RD_STATS tmp_rdc;
        const PICK_MODE_CONTEXT *const ctx_v = &pc_tree->vertical[0];
        av1_init_rd_stats(&tmp_rdc);
        update_state(cpi, td, ctx_v, mi_row, mi_col, subsize, 1);
        encode_superblock(cpi, tile_data, td, tp, DRY_RUN_NORMAL, subsize,
                          NULL);
        pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + hbs, &tmp_rdc,
                      PARTITION_VERT, subsize,
                      &pc_tree->vertical[bsize > BLOCK_8X8], invalid_rdc,
                      PICK_MODE_RD);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          av1_invalid_rd_stats(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
        last_part_rdc.rdcost += tmp_rdc.rdcost;
      }
      break;
    case PARTITION_SPLIT:
      last_part_rdc.rate = 0;
      last_part_rdc.dist = 0;
      last_part_rdc.rdcost = 0;
      for (i = 0; i < 4; i++) {
        int x_idx = (i & 1) * hbs;
        int y_idx = (i >> 1) * hbs;
        int jj = i >> 1, ii = i & 0x01;
        RD_STATS tmp_rdc;
        if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
          continue;

        av1_init_rd_stats(&tmp_rdc);
        rd_use_partition(cpi, td, tile_data,
                         mib + jj * hbs * cm->mi_stride + ii * hbs, tp,
                         mi_row + y_idx, mi_col + x_idx, subsize, &tmp_rdc.rate,
                         &tmp_rdc.dist, i != 3, pc_tree->split[i]);
        if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
          av1_invalid_rd_stats(&last_part_rdc);
          break;
        }
        last_part_rdc.rate += tmp_rdc.rate;
        last_part_rdc.dist += tmp_rdc.dist;
      }
      break;
    case PARTITION_VERT_A:
    case PARTITION_VERT_B:
    case PARTITION_HORZ_A:
    case PARTITION_HORZ_B:
    case PARTITION_HORZ_4:
    case PARTITION_VERT_4:
      assert(0 && "Cannot handle extended partition types");
    default: assert(0); break;
  }

  if (last_part_rdc.rate < INT_MAX) {
    last_part_rdc.rate += x->partition_cost[pl][partition];
    last_part_rdc.rdcost =
        RDCOST(x->rdmult, last_part_rdc.rate, last_part_rdc.dist);
  }

  if (do_partition_search && cpi->sf.adjust_partitioning_from_last_frame &&
      cpi->sf.partition_search_type == SEARCH_PARTITION &&
      partition != PARTITION_SPLIT && bsize > BLOCK_8X8 &&
      (mi_row + bs < cm->mi_rows || mi_row + hbs == cm->mi_rows) &&
      (mi_col + bs < cm->mi_cols || mi_col + hbs == cm->mi_cols)) {
    BLOCK_SIZE split_subsize = get_partition_subsize(bsize, PARTITION_SPLIT);
    chosen_rdc.rate = 0;
    chosen_rdc.dist = 0;

    restore_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);
    pc_tree->partitioning = PARTITION_SPLIT;

    // Split partition.
    for (i = 0; i < 4; i++) {
      int x_idx = (i & 1) * hbs;
      int y_idx = (i >> 1) * hbs;
      RD_STATS tmp_rdc;

      if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
        continue;

      save_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);
      pc_tree->split[i]->partitioning = PARTITION_NONE;
      pick_sb_modes(cpi, tile_data, x, mi_row + y_idx, mi_col + x_idx, &tmp_rdc,
                    PARTITION_SPLIT, split_subsize, &pc_tree->split[i]->none,
                    invalid_rdc, PICK_MODE_RD);

      restore_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);
      if (tmp_rdc.rate == INT_MAX || tmp_rdc.dist == INT64_MAX) {
        av1_invalid_rd_stats(&chosen_rdc);
        break;
      }

      chosen_rdc.rate += tmp_rdc.rate;
      chosen_rdc.dist += tmp_rdc.dist;

      if (i != 3)
        encode_sb(cpi, td, tile_data, tp, mi_row + y_idx, mi_col + x_idx,
                  OUTPUT_ENABLED, split_subsize, pc_tree->split[i], NULL);

      chosen_rdc.rate += x->partition_cost[pl][PARTITION_NONE];
    }
    if (chosen_rdc.rate < INT_MAX) {
      chosen_rdc.rate += x->partition_cost[pl][PARTITION_SPLIT];
      chosen_rdc.rdcost = RDCOST(x->rdmult, chosen_rdc.rate, chosen_rdc.dist);
    }
  }

  // If last_part is better set the partitioning to that.
  if (last_part_rdc.rdcost < chosen_rdc.rdcost) {
    mib[0]->sb_type = bsize;
    if (bsize >= BLOCK_8X8) pc_tree->partitioning = partition;
    chosen_rdc = last_part_rdc;
  }
  // If none was better set the partitioning to that.
  if (none_rdc.rdcost < chosen_rdc.rdcost) {
    if (bsize >= BLOCK_8X8) pc_tree->partitioning = PARTITION_NONE;
    chosen_rdc = none_rdc;
  }

  restore_context(x, &x_ctx, mi_row, mi_col, bsize, num_planes);

  // We must have chosen a partitioning and encoding or we'll fail later on.
  // No other opportunities for success.
  if (bsize == cm->seq_params.sb_size)
    assert(chosen_rdc.rate < INT_MAX && chosen_rdc.dist < INT64_MAX);

  if (do_recon) {
    if (bsize == cm->seq_params.sb_size) {
      // NOTE: To get estimate for rate due to the tokens, use:
      // int rate_coeffs = 0;
      // encode_sb(cpi, td, tile_data, tp, mi_row, mi_col, DRY_RUN_COSTCOEFFS,
      //           bsize, pc_tree, &rate_coeffs);
      x->cb_offset = 0;
      encode_sb(cpi, td, tile_data, tp, mi_row, mi_col, OUTPUT_ENABLED, bsize,
                pc_tree, NULL);
    } else {
      encode_sb(cpi, td, tile_data, tp, mi_row, mi_col, DRY_RUN_NORMAL, bsize,
                pc_tree, NULL);
    }
  }

  *rate = chosen_rdc.rate;
  *dist = chosen_rdc.dist;
  x->rdmult = orig_rdmult;
}

static int is_leaf_split_partition(AV1_COMMON *cm, int mi_row, int mi_col,
                                   BLOCK_SIZE bsize) {
  const int bs = mi_size_wide[bsize];
  const int hbs = bs / 2;
  assert(bsize >= BLOCK_8X8);
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, PARTITION_SPLIT);
  for (int i = 0; i < 4; i++) {
    int x_idx = (i & 1) * hbs;
    int y_idx = (i >> 1) * hbs;
    if ((mi_row + y_idx >= cm->mi_rows) || (mi_col + x_idx >= cm->mi_cols))
      return 0;
    if (get_partition(cm, mi_row + y_idx, mi_col + x_idx, subsize) !=
        PARTITION_NONE)
      return 0;
  }
  return 1;
}

static AOM_INLINE void nonrd_use_partition(AV1_COMP *cpi, ThreadData *td,
                                           TileDataEnc *tile_data,
                                           MB_MODE_INFO **mib, TOKENEXTRA **tp,
                                           int mi_row, int mi_col,
                                           BLOCK_SIZE bsize, PC_TREE *pc_tree) {
  AV1_COMMON *const cm = &cpi->common;
  TileInfo *const tile_info = &tile_data->tile_info;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  // Only square blocks from 8x8 to 128x128 are supported
  assert(bsize >= BLOCK_8X8 && bsize <= BLOCK_128X128);
  const int bs = mi_size_wide[bsize];
  const int hbs = bs / 2;
  const PARTITION_TYPE partition =
      (bsize >= BLOCK_8X8) ? get_partition(cm, mi_row, mi_col, bsize)
                           : PARTITION_NONE;
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, partition);
  assert(subsize <= BLOCK_LARGEST);
  const int pl = (bsize >= BLOCK_8X8)
                     ? partition_plane_context(xd, mi_row, mi_col, bsize)
                     : 0;

  RD_STATS dummy_cost;
  av1_invalid_rd_stats(&dummy_cost);
  RD_STATS invalid_rd;
  av1_invalid_rd_stats(&invalid_rd);

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  assert(mi_size_wide[bsize] == mi_size_high[bsize]);

  pc_tree->partitioning = partition;

  xd->above_txfm_context = cm->above_txfm_context[tile_info->tile_row] + mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);

  switch (partition) {
    case PARTITION_NONE:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &dummy_cost,
                    PARTITION_NONE, bsize, &pc_tree->none, invalid_rd,
                    PICK_MODE_NONRD);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, 0, bsize, partition,
               &pc_tree->none, NULL);
      break;
    case PARTITION_VERT:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &dummy_cost,
                    PARTITION_VERT, subsize, &pc_tree->vertical[0], invalid_rd,
                    PICK_MODE_NONRD);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, 0, subsize,
               PARTITION_VERT, &pc_tree->vertical[0], NULL);
      if (mi_col + hbs < cm->mi_cols && bsize > BLOCK_8X8) {
        pick_sb_modes(cpi, tile_data, x, mi_row, mi_col + hbs, &dummy_cost,
                      PARTITION_VERT, subsize, &pc_tree->vertical[1],
                      invalid_rd, PICK_MODE_NONRD);
        encode_b(cpi, tile_data, td, tp, mi_row, mi_col + hbs, 0, subsize,
                 PARTITION_VERT, &pc_tree->vertical[1], NULL);
      }
      break;
    case PARTITION_HORZ:
      pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &dummy_cost,
                    PARTITION_HORZ, subsize, &pc_tree->horizontal[0],
                    invalid_rd, PICK_MODE_NONRD);
      encode_b(cpi, tile_data, td, tp, mi_row, mi_col, 0, subsize,
               PARTITION_HORZ, &pc_tree->horizontal[0], NULL);

      if (mi_row + hbs < cm->mi_rows && bsize > BLOCK_8X8) {
        pick_sb_modes(cpi, tile_data, x, mi_row + hbs, mi_col, &dummy_cost,
                      PARTITION_HORZ, subsize, &pc_tree->horizontal[1],
                      invalid_rd, PICK_MODE_NONRD);
        encode_b(cpi, tile_data, td, tp, mi_row + hbs, mi_col, 0, subsize,
                 PARTITION_HORZ, &pc_tree->horizontal[1], NULL);
      }
      break;
    case PARTITION_SPLIT:
      if (cpi->sf.nonrd_merge_partition &&
          is_leaf_split_partition(cm, mi_row, mi_col, bsize) &&
          !frame_is_intra_only(cm)) {
        RD_SEARCH_MACROBLOCK_CONTEXT x_ctx;
        RD_STATS split_rdc, none_rdc;
        av1_init_rd_stats(&split_rdc);
        av1_invalid_rd_stats(&none_rdc);
        save_context(x, &x_ctx, mi_row, mi_col, bsize, 3);
        xd->above_txfm_context =
            cm->above_txfm_context[tile_info->tile_row] + mi_col;
        xd->left_txfm_context =
            xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);
        pc_tree->partitioning = PARTITION_NONE;
        pick_sb_modes(cpi, tile_data, x, mi_row, mi_col, &none_rdc,
                      PARTITION_NONE, bsize, &pc_tree->none, invalid_rd,
                      PICK_MODE_NONRD);
        none_rdc.rate += x->partition_cost[pl][PARTITION_NONE];
        none_rdc.rdcost = RDCOST(x->rdmult, none_rdc.rate, none_rdc.dist);
        restore_context(x, &x_ctx, mi_row, mi_col, bsize, 3);

        for (int i = 0; i < 4; i++) {
          RD_STATS block_rdc;
          av1_invalid_rd_stats(&block_rdc);
          int x_idx = (i & 1) * hbs;
          int y_idx = (i >> 1) * hbs;
          if ((mi_row + y_idx >= cm->mi_rows) ||
              (mi_col + x_idx >= cm->mi_cols))
            continue;
          xd->above_txfm_context =
              cm->above_txfm_context[tile_info->tile_row] + mi_col + x_idx;
          xd->left_txfm_context =
              xd->left_txfm_context_buffer + ((mi_row + y_idx) & MAX_MIB_MASK);
          pc_tree->split[i]->partitioning = PARTITION_NONE;
          pick_sb_modes(cpi, tile_data, x, mi_row + y_idx, mi_col + x_idx,
                        &block_rdc, PARTITION_NONE, subsize,
                        &pc_tree->split[i]->none, invalid_rd, PICK_MODE_NONRD);
          split_rdc.rate += block_rdc.rate;
          split_rdc.dist += block_rdc.dist;

          encode_b(cpi, tile_data, td, tp, mi_row + y_idx, mi_col + x_idx, 1,
                   subsize, PARTITION_NONE, &pc_tree->split[i]->none, NULL);
        }
        restore_context(x, &x_ctx, mi_row, mi_col, bsize, 3);
        split_rdc.rate += x->partition_cost[pl][PARTITION_SPLIT];
        split_rdc.rdcost = RDCOST(x->rdmult, split_rdc.rate, split_rdc.dist);
        if (none_rdc.rdcost < split_rdc.rdcost) {
          set_offsets_without_segment_id(cpi, &tile_data->tile_info, x, mi_row,
                                         mi_col, bsize);
          mib[0]->sb_type = bsize;
          pc_tree->partitioning = PARTITION_NONE;
          encode_b(cpi, tile_data, td, tp, mi_row, mi_col, 0, bsize, partition,
                   &pc_tree->none, NULL);
        } else {
          set_offsets_without_segment_id(cpi, &tile_data->tile_info, x, mi_row,
                                         mi_col, bsize);
          mib[0]->sb_type = subsize;
          pc_tree->partitioning = PARTITION_SPLIT;
          for (int i = 0; i < 4; i++) {
            int x_idx = (i & 1) * hbs;
            int y_idx = (i >> 1) * hbs;
            if ((mi_row + y_idx >= cm->mi_rows) ||
                (mi_col + x_idx >= cm->mi_cols))
              continue;

            set_offsets_without_segment_id(cpi, &tile_data->tile_info, x,
                                           mi_row + y_idx, mi_col + x_idx,
                                           subsize);
            encode_b(cpi, tile_data, td, tp, mi_row + y_idx, mi_col + x_idx, 0,
                     subsize, PARTITION_NONE, &pc_tree->split[i]->none, NULL);
          }
        }
      } else {
        for (int i = 0; i < 4; i++) {
          int x_idx = (i & 1) * hbs;
          int y_idx = (i >> 1) * hbs;
          int jj = i >> 1, ii = i & 0x01;
          if ((mi_row + y_idx >= cm->mi_rows) ||
              (mi_col + x_idx >= cm->mi_cols))
            continue;
          nonrd_use_partition(
              cpi, td, tile_data, mib + jj * hbs * cm->mi_stride + ii * hbs, tp,
              mi_row + y_idx, mi_col + x_idx, subsize, pc_tree->split[i]);
        }
      }
      break;
    case PARTITION_VERT_A:
    case PARTITION_VERT_B:
    case PARTITION_HORZ_A:
    case PARTITION_HORZ_B:
    case PARTITION_HORZ_4:
    case PARTITION_VERT_4:
      assert(0 && "Cannot handle extended partition types");
    default: assert(0); break;
  }
}

>>>>>>> CHANGE (5861e8 WIP: interp filter)
#if !CONFIG_REALTIME_ONLY
/*!\brief Assigns different quantization parameters to each super
 * block based on its TPL weight.
 *
 * \ingroup tpl_modelling
 *
 * \param[in]     cpi         Top level encoder instance structure
 * \param[in,out] td          Thread data structure
 * \param[in,out] x           Macro block level data for this block.
 * \param[in]     tile_info   Tile infromation / identification
 * \param[in]     mi_row      Block row (in "MI_SIZE" units) index
 * \param[in]     mi_col      Block column (in "MI_SIZE" units) index
 * \param[out]    num_planes  Number of image planes (e.g. Y,U,V)
 *
 * \return No return value but updates macroblock and thread data
 * related to the q / q delta to be used.
 */
static AOM_INLINE void setup_delta_q(AV1_COMP *const cpi, ThreadData *td,
                                     MACROBLOCK *const x,
                                     const TileInfo *const tile_info,
                                     int mi_row, int mi_col, int num_planes) {
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const DeltaQInfo *const delta_q_info = &cm->delta_q_info;
  assert(delta_q_info->delta_q_present_flag);

  const BLOCK_SIZE sb_size = cm->seq_params->sb_size;
  // Delta-q modulation based on variance
  av1_setup_src_planes(x, cpi->source, mi_row, mi_col, num_planes, sb_size);

  const int delta_q_res = delta_q_info->delta_q_res;
  int current_qindex = cm->quant_params.base_qindex;
  if (cpi->oxcf.q_cfg.deltaq_mode == DELTA_Q_PERCEPTUAL) {
    if (DELTA_Q_PERCEPTUAL_MODULATION == 1) {
      const int block_wavelet_energy_level =
          av1_block_wavelet_energy_level(cpi, x, sb_size);
      x->sb_energy_level = block_wavelet_energy_level;
      current_qindex = av1_compute_q_from_energy_level_deltaq_mode(
          cpi, block_wavelet_energy_level);
    } else {
      const int block_var_level = av1_log_block_var(cpi, x, sb_size);
      x->sb_energy_level = block_var_level;
      current_qindex =
          av1_compute_q_from_energy_level_deltaq_mode(cpi, block_var_level);
    }
  } else if (cpi->oxcf.q_cfg.deltaq_mode == DELTA_Q_OBJECTIVE &&
             cpi->oxcf.algo_cfg.enable_tpl_model) {
    // Setup deltaq based on tpl stats
    current_qindex =
        av1_get_q_for_deltaq_objective(cpi, sb_size, mi_row, mi_col);
  } else if (cpi->oxcf.q_cfg.deltaq_mode == DELTA_Q_PERCEPTUAL_AI) {
    current_qindex = av1_get_sbq_perceptual_ai(cpi, sb_size, mi_row, mi_col);
  } else if (cpi->oxcf.q_cfg.deltaq_mode == DELTA_Q_USER_RATING_BASED) {
    current_qindex = av1_get_sbq_user_rating_based(cpi, mi_row, mi_col);
  } else if (cpi->oxcf.q_cfg.enable_hdr_deltaq) {
    current_qindex = av1_get_q_for_hdr(cpi, x, sb_size, mi_row, mi_col);
  }

  MACROBLOCKD *const xd = &x->e_mbd;
  current_qindex = av1_adjust_q_from_delta_q_res(
      delta_q_res, xd->current_base_qindex, current_qindex);

  x->delta_qindex = current_qindex - cm->quant_params.base_qindex;
  av1_set_offsets(cpi, tile_info, x, mi_row, mi_col, sb_size);
  xd->mi[0]->current_qindex = current_qindex;
  av1_init_plane_quantizers(cpi, x, xd->mi[0]->segment_id);

  // keep track of any non-zero delta-q used
  td->deltaq_used |= (x->delta_qindex != 0);

  if (cpi->oxcf.tool_cfg.enable_deltalf_mode) {
    const int delta_lf_res = delta_q_info->delta_lf_res;
    const int lfmask = ~(delta_lf_res - 1);
    const int delta_lf_from_base =
        ((x->delta_qindex / 4 + delta_lf_res / 2) & lfmask);
    const int8_t delta_lf =
        (int8_t)clamp(delta_lf_from_base, -MAX_LOOP_FILTER, MAX_LOOP_FILTER);
    const int frame_lf_count =
        av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
    const int mib_size = cm->seq_params->mib_size;

    // pre-set the delta lf for loop filter. Note that this value is set
    // before mi is assigned for each block in current superblock
    for (int j = 0; j < AOMMIN(mib_size, mi_params->mi_rows - mi_row); j++) {
      for (int k = 0; k < AOMMIN(mib_size, mi_params->mi_cols - mi_col); k++) {
        const int grid_idx = get_mi_grid_idx(mi_params, mi_row + j, mi_col + k);
        mi_params->mi_alloc[grid_idx].delta_lf_from_base = delta_lf;
        for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
          mi_params->mi_alloc[grid_idx].delta_lf[lf_id] = delta_lf;
        }
      }
    }
  }
}

static void init_ref_frame_space(AV1_COMP *cpi, ThreadData *td, int mi_row,
                                 int mi_col) {
  const AV1_COMMON *cm = &cpi->common;
  const GF_GROUP *const gf_group = &cpi->ppi->gf_group;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  MACROBLOCK *x = &td->mb;
  const int frame_idx = cpi->gf_frame_index;
  TplParams *const tpl_data = &cpi->ppi->tpl_data;
  const uint8_t block_mis_log2 = tpl_data->tpl_stats_block_mis_log2;

  av1_zero(x->tpl_keep_ref_frame);

  if (!av1_tpl_stats_ready(tpl_data, frame_idx)) return;
  if (!is_frame_tpl_eligible(gf_group, cpi->gf_frame_index)) return;
  if (cpi->oxcf.q_cfg.aq_mode != NO_AQ) return;

  const int is_overlay =
      cpi->ppi->gf_group.update_type[frame_idx] == OVERLAY_UPDATE;
  if (is_overlay) {
    memset(x->tpl_keep_ref_frame, 1, sizeof(x->tpl_keep_ref_frame));
    return;
  }

  TplDepFrame *tpl_frame = &tpl_data->tpl_frame[frame_idx];
  TplDepStats *tpl_stats = tpl_frame->tpl_stats_ptr;
  const int tpl_stride = tpl_frame->stride;
  int64_t inter_cost[INTER_REFS_PER_FRAME] = { 0 };
  const int step = 1 << block_mis_log2;
  const BLOCK_SIZE sb_size = cm->seq_params->sb_size;

  const int mi_row_end =
      AOMMIN(mi_size_high[sb_size] + mi_row, mi_params->mi_rows);
  const int mi_cols_sr = av1_pixels_to_mi(cm->superres_upscaled_width);
  const int mi_col_sr =
      coded_to_superres_mi(mi_col, cm->superres_scale_denominator);
  const int mi_col_end_sr =
      AOMMIN(coded_to_superres_mi(mi_col + mi_size_wide[sb_size],
                                  cm->superres_scale_denominator),
             mi_cols_sr);
  const int row_step = step;
  const int col_step_sr =
      coded_to_superres_mi(step, cm->superres_scale_denominator);
  for (int row = mi_row; row < mi_row_end; row += row_step) {
    for (int col = mi_col_sr; col < mi_col_end_sr; col += col_step_sr) {
      const TplDepStats *this_stats =
          &tpl_stats[av1_tpl_ptr_pos(row, col, tpl_stride, block_mis_log2)];
      int64_t tpl_pred_error[INTER_REFS_PER_FRAME] = { 0 };
      // Find the winner ref frame idx for the current block
      int64_t best_inter_cost = this_stats->pred_error[0];
      int best_rf_idx = 0;
      for (int idx = 1; idx < INTER_REFS_PER_FRAME; ++idx) {
        if ((this_stats->pred_error[idx] < best_inter_cost) &&
            (this_stats->pred_error[idx] != 0)) {
          best_inter_cost = this_stats->pred_error[idx];
          best_rf_idx = idx;
        }
      }
      // tpl_pred_error is the pred_error reduction of best_ref w.r.t.
      // LAST_FRAME.
      tpl_pred_error[best_rf_idx] = this_stats->pred_error[best_rf_idx] -
                                    this_stats->pred_error[LAST_FRAME - 1];

      for (int rf_idx = 1; rf_idx < INTER_REFS_PER_FRAME; ++rf_idx)
        inter_cost[rf_idx] += tpl_pred_error[rf_idx];
    }
  }

  int rank_index[INTER_REFS_PER_FRAME - 1];
  for (int idx = 0; idx < INTER_REFS_PER_FRAME - 1; ++idx) {
    rank_index[idx] = idx + 1;
    for (int i = idx; i > 0; --i) {
      if (inter_cost[rank_index[i - 1]] > inter_cost[rank_index[i]]) {
        const int tmp = rank_index[i - 1];
        rank_index[i - 1] = rank_index[i];
        rank_index[i] = tmp;
      }
    }
  }

  x->tpl_keep_ref_frame[INTRA_FRAME] = 1;
  x->tpl_keep_ref_frame[LAST_FRAME] = 1;

  int cutoff_ref = 0;
  for (int idx = 0; idx < INTER_REFS_PER_FRAME - 1; ++idx) {
    x->tpl_keep_ref_frame[rank_index[idx] + LAST_FRAME] = 1;
    if (idx > 2) {
      if (!cutoff_ref) {
        // If the predictive coding gains are smaller than the previous more
        // relevant frame over certain amount, discard this frame and all the
        // frames afterwards.
        if (llabs(inter_cost[rank_index[idx]]) <
                llabs(inter_cost[rank_index[idx - 1]]) / 8 ||
            inter_cost[rank_index[idx]] == 0)
          cutoff_ref = 1;
      }

      if (cutoff_ref) x->tpl_keep_ref_frame[rank_index[idx] + LAST_FRAME] = 0;
    }
  }
}

static AOM_INLINE void adjust_rdmult_tpl_model(AV1_COMP *cpi, MACROBLOCK *x,
                                               int mi_row, int mi_col) {
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int orig_rdmult = cpi->rd.RDMULT;

  assert(IMPLIES(cpi->ppi->gf_group.size > 0,
                 cpi->gf_frame_index < cpi->ppi->gf_group.size));
  const int gf_group_index = cpi->gf_frame_index;
  if (cpi->oxcf.algo_cfg.enable_tpl_model && cpi->oxcf.q_cfg.aq_mode == NO_AQ &&
      cpi->oxcf.q_cfg.deltaq_mode == NO_DELTA_Q && gf_group_index > 0 &&
      cpi->ppi->gf_group.update_type[gf_group_index] == ARF_UPDATE) {
    const int dr =
        av1_get_rdmult_delta(cpi, sb_size, mi_row, mi_col, orig_rdmult);
    x->rdmult = dr;
  }
}
#endif  // !CONFIG_REALTIME_ONLY

#if CONFIG_RT_ML_PARTITIONING
// Get a prediction(stored in x->est_pred) for the whole superblock.
static void get_estimated_pred(AV1_COMP *cpi, const TileInfo *const tile,
                               MACROBLOCK *x, int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int is_key_frame = frame_is_intra_only(cm);
  MACROBLOCKD *xd = &x->e_mbd;

  // TODO(kyslov) Extend to 128x128
  assert(cm->seq_params->sb_size == BLOCK_64X64);

  av1_set_offsets(cpi, tile, x, mi_row, mi_col, BLOCK_64X64);

  if (!is_key_frame) {
    MB_MODE_INFO *mi = xd->mi[0];
    const YV12_BUFFER_CONFIG *yv12 = get_ref_frame_yv12_buf(cm, LAST_FRAME);

    assert(yv12 != NULL);

    av1_setup_pre_planes(xd, 0, yv12, mi_row, mi_col,
                         get_ref_scale_factors(cm, LAST_FRAME), 1);
    mi->ref_frame[0] = LAST_FRAME;
    mi->ref_frame[1] = NONE;
    mi->bsize = BLOCK_64X64;
    mi->mv[0].as_int = 0;
    mi->interp_filters = av1_broadcast_interp_filter(BILINEAR);

    set_ref_ptrs(cm, xd, mi->ref_frame[0], mi->ref_frame[1]);

    xd->plane[0].dst.buf = x->est_pred;
    xd->plane[0].dst.stride = 64;
    av1_enc_build_inter_predictor_y(xd, mi_row, mi_col);
  } else {
#if CONFIG_AV1_HIGHBITDEPTH
    switch (xd->bd) {
      case 8: memset(x->est_pred, 128, 64 * 64 * sizeof(x->est_pred[0])); break;
      case 10:
        memset(x->est_pred, 128 * 4, 64 * 64 * sizeof(x->est_pred[0]));
        break;
      case 12:
        memset(x->est_pred, 128 * 16, 64 * 64 * sizeof(x->est_pred[0]));
        break;
    }
#else
    memset(x->est_pred, 128, 64 * 64 * sizeof(x->est_pred[0]));
#endif  // CONFIG_VP9_HIGHBITDEPTH
  }
}
#endif  // CONFIG_RT_ML_PARTITIONING

#define AVG_CDF_WEIGHT_LEFT 3
#define AVG_CDF_WEIGHT_TOP_RIGHT 1

/*!\brief Encode a superblock (minimal RD search involved)
 *
 * \ingroup partition_search
 * Encodes the superblock by a pre-determined partition pattern, only minor
 * rd-based searches are allowed to adjust the initial pattern. It is only used
 * by realtime encoding.
 */
static AOM_INLINE void encode_nonrd_sb(AV1_COMP *cpi, ThreadData *td,
                                       TileDataEnc *tile_data, TokenExtra **tp,
                                       const int mi_row, const int mi_col,
                                       const int seg_skip) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  const SPEED_FEATURES *const sf = &cpi->sf;
  const TileInfo *const tile_info = &tile_data->tile_info;
  MB_MODE_INFO **mi = cm->mi_params.mi_grid_base +
                      get_mi_grid_idx(&cm->mi_params, mi_row, mi_col);
  const BLOCK_SIZE sb_size = cm->seq_params->sb_size;

  // Grade the temporal variation of the sb, the grade will be used to decide
  // fast mode search strategy for coding blocks
  if (sf->rt_sf.source_metrics_sb_nonrd &&
      cpi->svc.number_spatial_layers <= 1 &&
      cm->current_frame.frame_type != KEY_FRAME) {
    int offset = cpi->source->y_stride * (mi_row << 2) + (mi_col << 2);
    av1_source_content_sb(cpi, x, offset);
  }
#if CONFIG_RT_ML_PARTITIONING
  if (sf->part_sf.partition_search_type == ML_BASED_PARTITION) {
    PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);
    RD_STATS dummy_rdc;
    get_estimated_pred(cpi, tile_info, x, mi_row, mi_col);
    av1_nonrd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col,
                             BLOCK_64X64, &dummy_rdc, 1, INT64_MAX, pc_root);
    av1_free_pc_tree_recursive(pc_root, av1_num_planes(cm), 0, 0);
    return;
  }
#endif
  // Set the partition
  if (sf->part_sf.partition_search_type == FIXED_PARTITION || seg_skip) {
    // set a fixed-size partition
    av1_set_offsets(cpi, tile_info, x, mi_row, mi_col, sb_size);
    const BLOCK_SIZE bsize =
        seg_skip ? sb_size : sf->part_sf.fixed_partition_size;
    av1_set_fixed_partitioning(cpi, tile_info, mi, mi_row, mi_col, bsize);
  } else if (sf->part_sf.partition_search_type == VAR_BASED_PARTITION) {
    // set a variance-based partition
    av1_set_offsets(cpi, tile_info, x, mi_row, mi_col, sb_size);
    av1_choose_var_based_partitioning(cpi, tile_info, td, x, mi_row, mi_col);
  }
  assert(sf->part_sf.partition_search_type == FIXED_PARTITION || seg_skip ||
         sf->part_sf.partition_search_type == VAR_BASED_PARTITION);
  set_cb_offsets(td->mb.cb_offset, 0, 0);

  // Adjust and encode the superblock
  PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);

  // Initialize the flag to skip cdef to 1.
  if (sf->rt_sf.skip_cdef_sb) {
    // If 128x128 block is used, we need to set the flag for all 4 64x64 sub
    // "blocks".
    const int block64_in_sb = (sb_size == BLOCK_128X128) ? 2 : 1;
    for (int r = 0; r < block64_in_sb; ++r) {
      for (int c = 0; c < block64_in_sb; ++c) {
        const int idx_in_sb =
            r * MI_SIZE_64X64 * cm->mi_params.mi_stride + c * MI_SIZE_64X64;
        if (mi[idx_in_sb]) mi[idx_in_sb]->skip_cdef_curr_sb = 1;
      }
    }
  }

  av1_nonrd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, sb_size,
                          pc_root);

  if (sf->rt_sf.skip_cdef_sb) {
    // If 128x128 block is used, we need to set the flag for all 4 64x64 sub
    // "blocks".
    const int block64_in_sb = (sb_size == BLOCK_128X128) ? 2 : 1;
    const int skip = mi[0]->skip_cdef_curr_sb;
    for (int r = 0; r < block64_in_sb; ++r) {
      for (int c = 0; c < block64_in_sb; ++c) {
        const int idx_in_sb =
            r * MI_SIZE_64X64 * cm->mi_params.mi_stride + c * MI_SIZE_64X64;
        if (mi[idx_in_sb]) mi[idx_in_sb]->skip_cdef_curr_sb = skip;
      }
    }
  }
  av1_free_pc_tree_recursive(pc_root, av1_num_planes(cm), 0, 0);
}

// This function initializes the stats for encode_rd_sb.
static INLINE void init_encode_rd_sb(AV1_COMP *cpi, ThreadData *td,
                                     const TileDataEnc *tile_data,
                                     SIMPLE_MOTION_DATA_TREE *sms_root,
                                     RD_STATS *rd_cost, int mi_row, int mi_col,
                                     int gather_tpl_data) {
  const AV1_COMMON *cm = &cpi->common;
  const TileInfo *tile_info = &tile_data->tile_info;
  MACROBLOCK *x = &td->mb;

  const SPEED_FEATURES *sf = &cpi->sf;
  const int use_simple_motion_search =
      (sf->part_sf.simple_motion_search_split ||
       sf->part_sf.simple_motion_search_prune_rect ||
       sf->part_sf.simple_motion_search_early_term_none ||
       sf->part_sf.ml_early_term_after_part_split_level) &&
      !frame_is_intra_only(cm);
  if (use_simple_motion_search) {
    av1_init_simple_motion_search_mvs_for_sb(cpi, tile_info, x, sms_root,
                                             mi_row, mi_col);
  }

#if !CONFIG_REALTIME_ONLY
  if (!(has_no_stats_stage(cpi) && cpi->oxcf.mode == REALTIME &&
        cpi->oxcf.gf_cfg.lag_in_frames == 0)) {
    init_ref_frame_space(cpi, td, mi_row, mi_col);
    x->sb_energy_level = 0;
    x->part_search_info.cnn_output_valid = 0;
    if (gather_tpl_data) {
      if (cm->delta_q_info.delta_q_present_flag) {
        const int num_planes = av1_num_planes(cm);
        const BLOCK_SIZE sb_size = cm->seq_params->sb_size;
        setup_delta_q(cpi, td, x, tile_info, mi_row, mi_col, num_planes);
        av1_tpl_rdmult_setup_sb(cpi, x, sb_size, mi_row, mi_col);
      }
      if (cpi->oxcf.algo_cfg.enable_tpl_model) {
        adjust_rdmult_tpl_model(cpi, x, mi_row, mi_col);
      }
    }
  }
#else
  (void)tile_info;
  (void)mi_row;
  (void)mi_col;
  (void)gather_tpl_data;
#endif

  reset_mb_rd_record(x->txfm_search_info.mb_rd_record);
  av1_zero(x->picked_ref_frames_mask);
  av1_invalid_rd_stats(rd_cost);
}

/*!\brief Encode a superblock (RD-search-based)
 *
 * \ingroup partition_search
 * Conducts partition search for a superblock, based on rate-distortion costs,
 * from scratch or adjusting from a pre-calculated partition pattern.
 */
static AOM_INLINE void encode_rd_sb(AV1_COMP *cpi, ThreadData *td,
                                    TileDataEnc *tile_data, TokenExtra **tp,
                                    const int mi_row, const int mi_col,
                                    const int seg_skip) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  const SPEED_FEATURES *const sf = &cpi->sf;
  const TileInfo *const tile_info = &tile_data->tile_info;
  MB_MODE_INFO **mi = cm->mi_params.mi_grid_base +
                      get_mi_grid_idx(&cm->mi_params, mi_row, mi_col);
  const BLOCK_SIZE sb_size = cm->seq_params->sb_size;
  const int num_planes = av1_num_planes(cm);
  int dummy_rate;
  int64_t dummy_dist;
  RD_STATS dummy_rdc;
  SIMPLE_MOTION_DATA_TREE *const sms_root = td->sms_root;

#if CONFIG_REALTIME_ONLY
  (void)seg_skip;
#endif  // CONFIG_REALTIME_ONLY

  init_encode_rd_sb(cpi, td, tile_data, sms_root, &dummy_rdc, mi_row, mi_col,
                    1);

  // Encode the superblock
  if (sf->part_sf.partition_search_type == VAR_BASED_PARTITION) {
#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, rd_use_partition_time);
#endif
    // partition search starting from a variance-based partition
    av1_set_offsets(cpi, tile_info, x, mi_row, mi_col, sb_size);
    av1_choose_var_based_partitioning(cpi, tile_info, td, x, mi_row, mi_col);
    PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);
    av1_rd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, sb_size,
                         &dummy_rate, &dummy_dist, 1, pc_root);
    av1_free_pc_tree_recursive(pc_root, num_planes, 0, 0);
#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, rd_use_partition_time);
#endif
  }
#if !CONFIG_REALTIME_ONLY
  else if (sf->part_sf.partition_search_type == FIXED_PARTITION || seg_skip) {
    // partition search by adjusting a fixed-size partition
    av1_set_offsets(cpi, tile_info, x, mi_row, mi_col, sb_size);
    const BLOCK_SIZE bsize =
        seg_skip ? sb_size : sf->part_sf.fixed_partition_size;
    av1_set_fixed_partitioning(cpi, tile_info, mi, mi_row, mi_col, bsize);
    PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);
    av1_rd_use_partition(cpi, td, tile_data, mi, tp, mi_row, mi_col, sb_size,
                         &dummy_rate, &dummy_dist, 1, pc_root);
    av1_free_pc_tree_recursive(pc_root, num_planes, 0, 0);
  } else {
    // The most exhaustive recursive partition search
    SuperBlockEnc *sb_enc = &x->sb_enc;
    // No stats for overlay frames. Exclude key frame.
    av1_get_tpl_stats_sb(cpi, sb_size, mi_row, mi_col, sb_enc);

    // Reset the tree for simple motion search data
    av1_reset_simple_motion_tree_partition(sms_root, sb_size);

#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, rd_pick_partition_time);
#endif

    // Estimate the maximum square partition block size, which will be used
    // as the starting block size for partitioning the sb
    set_max_min_partition_size(sb_enc, cpi, x, sf, sb_size, mi_row, mi_col);

    // The superblock can be searched only once, or twice consecutively for
    // better quality. Note that the meaning of passes here is different from
    // the general concept of 1-pass/2-pass encoders.
    const int num_passes =
        cpi->oxcf.unit_test_cfg.sb_multipass_unit_test ? 2 : 1;

    if (num_passes == 1) {
#if CONFIG_PARTITION_SEARCH_ORDER
      if (cpi->ext_part_controller.ready && !frame_is_intra_only(cm)) {
        av1_reset_part_sf(&cpi->sf.part_sf);
        av1_reset_sf_for_ext_part(cpi);
        RD_STATS this_rdc;
        av1_rd_partition_search(cpi, td, tile_data, tp, sms_root, mi_row,
                                mi_col, sb_size, &this_rdc);
      } else {
        PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);
        av1_rd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, sb_size,
                              &dummy_rdc, dummy_rdc, pc_root, sms_root, NULL,
                              SB_SINGLE_PASS, NULL);
      }
#else
      PC_TREE *const pc_root = av1_alloc_pc_tree_node(sb_size);
      av1_rd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, sb_size,
                            &dummy_rdc, dummy_rdc, pc_root, sms_root, NULL,
                            SB_SINGLE_PASS, NULL);
#endif  // CONFIG_PARTITION_SEARCH_ORDER
    } else {
      // First pass
      SB_FIRST_PASS_STATS sb_fp_stats;
      av1_backup_sb_state(&sb_fp_stats, cpi, td, tile_data, mi_row, mi_col);
      PC_TREE *const pc_root_p0 = av1_alloc_pc_tree_node(sb_size);
      av1_rd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, sb_size,
                            &dummy_rdc, dummy_rdc, pc_root_p0, sms_root, NULL,
                            SB_DRY_PASS, NULL);

      // Second pass
      init_encode_rd_sb(cpi, td, tile_data, sms_root, &dummy_rdc, mi_row,
                        mi_col, 0);
      av1_reset_mbmi(&cm->mi_params, sb_size, mi_row, mi_col);
      av1_reset_simple_motion_tree_partition(sms_root, sb_size);

      av1_restore_sb_state(&sb_fp_stats, cpi, td, tile_data, mi_row, mi_col);

      PC_TREE *const pc_root_p1 = av1_alloc_pc_tree_node(sb_size);
      av1_rd_pick_partition(cpi, td, tile_data, tp, mi_row, mi_col, sb_size,
                            &dummy_rdc, dummy_rdc, pc_root_p1, sms_root, NULL,
                            SB_WET_PASS, NULL);
    }
    // Reset to 0 so that it wouldn't be used elsewhere mistakenly.
    sb_enc->tpl_data_count = 0;
#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, rd_pick_partition_time);
#endif
  }
#endif  // !CONFIG_REALTIME_ONLY

  // Update the inter rd model
  // TODO(angiebird): Let inter_mode_rd_model_estimation support multi-tile.
  if (cpi->sf.inter_sf.inter_mode_rd_model_estimation == 1 &&
      cm->tiles.cols == 1 && cm->tiles.rows == 1) {
    av1_inter_mode_data_fit(tile_data, x->rdmult);
  }
}

// Check if the cost update of symbols mode, coeff and dv are tile or off.
static AOM_INLINE int is_mode_coeff_dv_upd_freq_tile_or_off(
    const AV1_COMP *const cpi, const CostUpdateFreq *const cost_upd_freq) {
  const INTER_MODE_SPEED_FEATURES *const inter_sf = &cpi->sf.inter_sf;

  return ((cost_upd_freq->coeff >= COST_UPD_TILE ||
           inter_sf->coeff_cost_upd_level == INTERNAL_COST_UPD_OFF) &&
          (cost_upd_freq->mode >= COST_UPD_TILE ||
           inter_sf->mode_cost_upd_level == INTERNAL_COST_UPD_OFF) &&
          (cost_upd_freq->dv >= COST_UPD_TILE ||
           cpi->sf.intra_sf.dv_cost_upd_level == INTERNAL_COST_UPD_OFF));
}

// When row-mt is enabled and cost update frequencies are set to off/tile,
// processing of current SB can start even before processing of top-right SB
// is finished. This function checks if it is sufficient to wait for top SB
// to finish processing before current SB starts processing.
static AOM_INLINE int delay_wait_for_top_right_sb(const AV1_COMP *const cpi) {
  const MODE mode = cpi->oxcf.mode;
  if (mode == GOOD) return 0;

  const CostUpdateFreq *const cost_upd_freq = &cpi->oxcf.cost_upd_freq;
  if (mode == ALLINTRA)
    return is_mode_coeff_dv_upd_freq_tile_or_off(cpi, cost_upd_freq);
  else if (mode == REALTIME)
    return (is_mode_coeff_dv_upd_freq_tile_or_off(cpi, cost_upd_freq) &&
            (cost_upd_freq->mv >= COST_UPD_TILE ||
             cpi->sf.inter_sf.mv_cost_upd_level == INTERNAL_COST_UPD_OFF));
  else
    return 0;
}

/*!\brief Encode a superblock row by breaking it into superblocks
 *
 * \ingroup partition_search
 * \callgraph
 * \callergraph
 * Do partition and mode search for an sb row: one row of superblocks filling up
 * the width of the current tile.
 */
static AOM_INLINE void encode_sb_row(AV1_COMP *cpi, ThreadData *td,
                                     TileDataEnc *tile_data, int mi_row,
                                     TokenExtra **tp) {
  AV1_COMMON *const cm = &cpi->common;
  const TileInfo *const tile_info = &tile_data->tile_info;
  MultiThreadInfo *const mt_info = &cpi->mt_info;
  AV1EncRowMultiThreadInfo *const enc_row_mt = &mt_info->enc_row_mt;
  AV1EncRowMultiThreadSync *const row_mt_sync = &tile_data->row_mt_sync;
  bool row_mt_enabled = mt_info->row_mt_enabled;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int sb_cols_in_tile = av1_get_sb_cols_in_tile(cm, tile_data->tile_info);
  const BLOCK_SIZE sb_size = cm->seq_params->sb_size;
  const int mib_size = cm->seq_params->mib_size;
  const int mib_size_log2 = cm->seq_params->mib_size_log2;
  const int sb_row = (mi_row - tile_info->mi_row_start) >> mib_size_log2;
  const int use_nonrd_mode = cpi->sf.rt_sf.use_nonrd_pick_mode;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, encode_sb_row_time);
#endif

  // Initialize the left context for the new SB row
  av1_zero_left_context(xd);

  // Reset delta for quantizer and loof filters at the beginning of every tile
  if (mi_row == tile_info->mi_row_start || row_mt_enabled) {
    if (cm->delta_q_info.delta_q_present_flag)
      xd->current_base_qindex = cm->quant_params.base_qindex;
    if (cm->delta_q_info.delta_lf_present_flag) {
      av1_reset_loop_filter_delta(xd, av1_num_planes(cm));
    }
  }

  reset_thresh_freq_fact(x);

  // Code each SB in the row
  for (int mi_col = tile_info->mi_col_start, sb_col_in_tile = 0;
       mi_col < tile_info->mi_col_end; mi_col += mib_size, sb_col_in_tile++) {
    // In realtime/allintra mode and when frequency of cost updates is off/tile,
    // wait for the top superblock to finish encoding. Otherwise, wait for the
    // top-right superblock to finish encoding.
    (*(enc_row_mt->sync_read_ptr))(
        row_mt_sync, sb_row, sb_col_in_tile - delay_wait_for_top_right_sb(cpi));
    const int update_cdf = tile_data->allow_update_cdf && row_mt_enabled;
    if (update_cdf && (tile_info->mi_row_start != mi_row)) {
      if ((tile_info->mi_col_start == mi_col)) {
        // restore frame context at the 1st column sb
        memcpy(xd->tile_ctx, x->row_ctx, sizeof(*xd->tile_ctx));
      } else {
        // update context
        int wt_left = AVG_CDF_WEIGHT_LEFT;
        int wt_tr = AVG_CDF_WEIGHT_TOP_RIGHT;
        if (tile_info->mi_col_end > (mi_col + mib_size))
          av1_avg_cdf_symbols(xd->tile_ctx, x->row_ctx + sb_col_in_tile,
                              wt_left, wt_tr);
        else
          av1_avg_cdf_symbols(xd->tile_ctx, x->row_ctx + sb_col_in_tile - 1,
                              wt_left, wt_tr);
      }
    }

    // Update the rate cost tables for some symbols
    av1_set_cost_upd_freq(cpi, td, tile_info, mi_row, mi_col);

    // Reset color coding related parameters
    x->color_sensitivity_sb[0] = 0;
    x->color_sensitivity_sb[1] = 0;
    x->color_sensitivity[0] = 0;
    x->color_sensitivity[1] = 0;
    x->content_state_sb.source_sad = kMedSad;
    x->content_state_sb.lighting_change = 0;
    x->content_state_sb.low_sumdiff = 0;

    xd->cur_frame_force_integer_mv = cm->features.cur_frame_force_integer_mv;
    x->source_variance = UINT_MAX;
    td->mb.cb_coef_buff = av1_get_cb_coeff_buffer(cpi, mi_row, mi_col);

    // Get segment id and skip flag
    const struct segmentation *const seg = &cm->seg;
    int seg_skip = 0;
    if (seg->enabled) {
      const uint8_t *const map =
          seg->update_map ? cpi->enc_seg.map : cm->last_frame_seg_map;
      const int segment_id =
          map ? get_segment_id(&cm->mi_params, map, sb_size, mi_row, mi_col)
              : 0;
      seg_skip = segfeature_active(seg, segment_id, SEG_LVL_SKIP);
    }

    produce_gradients_for_sb(cpi, x, sb_size, mi_row, mi_col);

    // encode the superblock
    if (use_nonrd_mode) {
      encode_nonrd_sb(cpi, td, tile_data, tp, mi_row, mi_col, seg_skip);
    } else {
      encode_rd_sb(cpi, td, tile_data, tp, mi_row, mi_col, seg_skip);
    }

    // Update the top-right context in row_mt coding
    if (update_cdf && (tile_info->mi_row_end > (mi_row + mib_size))) {
      if (sb_cols_in_tile == 1)
        memcpy(x->row_ctx, xd->tile_ctx, sizeof(*xd->tile_ctx));
      else if (sb_col_in_tile >= 1)
        memcpy(x->row_ctx + sb_col_in_tile - 1, xd->tile_ctx,
               sizeof(*xd->tile_ctx));
    }
    (*(enc_row_mt->sync_write_ptr))(row_mt_sync, sb_row, sb_col_in_tile,
                                    sb_cols_in_tile);
  }
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, encode_sb_row_time);
#endif
}

static AOM_INLINE void init_encode_frame_mb_context(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  // Copy data over into macro block data structures.
  av1_setup_src_planes(x, cpi->source, 0, 0, num_planes,
                       cm->seq_params->sb_size);

  av1_setup_block_planes(xd, cm->seq_params->subsampling_x,
                         cm->seq_params->subsampling_y, num_planes);
}

void av1_alloc_tile_data(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int tile_cols = cm->tiles.cols;
  const int tile_rows = cm->tiles.rows;

  if (cpi->tile_data != NULL) aom_free(cpi->tile_data);
  CHECK_MEM_ERROR(
      cm, cpi->tile_data,
      aom_memalign(32, tile_cols * tile_rows * sizeof(*cpi->tile_data)));

  cpi->allocated_tiles = tile_cols * tile_rows;
}

void av1_init_tile_data(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  const int tile_cols = cm->tiles.cols;
  const int tile_rows = cm->tiles.rows;
  int tile_col, tile_row;
  TokenInfo *const token_info = &cpi->token_info;
  TokenExtra *pre_tok = token_info->tile_tok[0][0];
  TokenList *tplist = token_info->tplist[0][0];
  unsigned int tile_tok = 0;
  int tplist_count = 0;

  if (!is_stat_generation_stage(cpi) &&
      cm->features.allow_screen_content_tools) {
    // Number of tokens for which token info needs to be allocated.
    unsigned int tokens_required =
        get_token_alloc(cm->mi_params.mb_rows, cm->mi_params.mb_cols,
                        MAX_SB_SIZE_LOG2, num_planes);
    // Allocate/reallocate memory for token related info if the number of tokens
    // required is more than the number of tokens already allocated. This could
    // occur in case of the following:
    // 1) If the memory is not yet allocated
    // 2) If the frame dimensions have changed
    const bool realloc_tokens = tokens_required > token_info->tokens_allocated;
    if (realloc_tokens) {
      free_token_info(token_info);
      alloc_token_info(cm, token_info, tokens_required);
      pre_tok = token_info->tile_tok[0][0];
      tplist = token_info->tplist[0][0];
    }
  }

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileDataEnc *const tile_data =
          &cpi->tile_data[tile_row * tile_cols + tile_col];
      TileInfo *const tile_info = &tile_data->tile_info;
      av1_tile_init(tile_info, cm, tile_row, tile_col);
      tile_data->firstpass_top_mv = kZeroMv;
      tile_data->abs_sum_level = 0;

      if (is_token_info_allocated(token_info)) {
        token_info->tile_tok[tile_row][tile_col] = pre_tok + tile_tok;
        pre_tok = token_info->tile_tok[tile_row][tile_col];
        tile_tok = allocated_tokens(
            *tile_info, cm->seq_params->mib_size_log2 + MI_SIZE_LOG2,
            num_planes);
        token_info->tplist[tile_row][tile_col] = tplist + tplist_count;
        tplist = token_info->tplist[tile_row][tile_col];
        tplist_count = av1_get_sb_rows_in_tile(cm, tile_data->tile_info);
      }
      tile_data->allow_update_cdf = !cm->tiles.large_scale;
      tile_data->allow_update_cdf = tile_data->allow_update_cdf &&
                                    !cm->features.disable_cdf_update &&
                                    !delay_wait_for_top_right_sb(cpi);
      tile_data->tctx = *cm->fc;
    }
  }
}

// Populate the start palette token info prior to encoding an SB row.
static AOM_INLINE void get_token_start(AV1_COMP *cpi, const TileInfo *tile_info,
                                       int tile_row, int tile_col, int mi_row,
                                       TokenExtra **tp) {
  const TokenInfo *token_info = &cpi->token_info;
  if (!is_token_info_allocated(token_info)) return;

  const AV1_COMMON *cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  TokenList *const tplist = cpi->token_info.tplist[tile_row][tile_col];
  const int sb_row_in_tile =
      (mi_row - tile_info->mi_row_start) >> cm->seq_params->mib_size_log2;

  get_start_tok(cpi, tile_row, tile_col, mi_row, tp,
                cm->seq_params->mib_size_log2 + MI_SIZE_LOG2, num_planes);
  assert(tplist != NULL);
  tplist[sb_row_in_tile].start = *tp;
}

// Populate the token count after encoding an SB row.
static AOM_INLINE void populate_token_count(AV1_COMP *cpi,
                                            const TileInfo *tile_info,
                                            int tile_row, int tile_col,
                                            int mi_row, TokenExtra *tok) {
  const TokenInfo *token_info = &cpi->token_info;
  if (!is_token_info_allocated(token_info)) return;

  const AV1_COMMON *cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  TokenList *const tplist = token_info->tplist[tile_row][tile_col];
  const int sb_row_in_tile =
      (mi_row - tile_info->mi_row_start) >> cm->seq_params->mib_size_log2;
  const int tile_mb_cols =
      (tile_info->mi_col_end - tile_info->mi_col_start + 2) >> 2;
  const int num_mb_rows_in_sb =
      ((1 << (cm->seq_params->mib_size_log2 + MI_SIZE_LOG2)) + 8) >> 4;
  tplist[sb_row_in_tile].count =
      (unsigned int)(tok - tplist[sb_row_in_tile].start);

  assert((unsigned int)(tok - tplist[sb_row_in_tile].start) <=
         get_token_alloc(num_mb_rows_in_sb, tile_mb_cols,
                         cm->seq_params->mib_size_log2 + MI_SIZE_LOG2,
                         num_planes));

  (void)num_planes;
  (void)tile_mb_cols;
  (void)num_mb_rows_in_sb;
}

/*!\brief Encode a superblock row
 *
 * \ingroup partition_search
 */
void av1_encode_sb_row(AV1_COMP *cpi, ThreadData *td, int tile_row,
                       int tile_col, int mi_row) {
  AV1_COMMON *const cm = &cpi->common;
  const int tile_cols = cm->tiles.cols;
  TileDataEnc *this_tile = &cpi->tile_data[tile_row * tile_cols + tile_col];
  const TileInfo *const tile_info = &this_tile->tile_info;
  TokenExtra *tok = NULL;

  get_token_start(cpi, tile_info, tile_row, tile_col, mi_row, &tok);

  encode_sb_row(cpi, td, this_tile, mi_row, &tok);

  populate_token_count(cpi, tile_info, tile_row, tile_col, mi_row, tok);
}

/*!\brief Encode a tile
 *
 * \ingroup partition_search
 */
void av1_encode_tile(AV1_COMP *cpi, ThreadData *td, int tile_row,
                     int tile_col) {
  AV1_COMMON *const cm = &cpi->common;
  TileDataEnc *const this_tile =
      &cpi->tile_data[tile_row * cm->tiles.cols + tile_col];
  const TileInfo *const tile_info = &this_tile->tile_info;

  if (!cpi->sf.rt_sf.use_nonrd_pick_mode) av1_inter_mode_data_init(this_tile);

  av1_zero_above_context(cm, &td->mb.e_mbd, tile_info->mi_col_start,
                         tile_info->mi_col_end, tile_row);
  av1_init_above_context(&cm->above_contexts, av1_num_planes(cm), tile_row,
                         &td->mb.e_mbd);

  if (cpi->oxcf.intra_mode_cfg.enable_cfl_intra)
    cfl_init(&td->mb.e_mbd.cfl, cm->seq_params);

  if (td->mb.txfm_search_info.mb_rd_record != NULL) {
    av1_crc32c_calculator_init(
        &td->mb.txfm_search_info.mb_rd_record->crc_calculator);
  }

  for (int mi_row = tile_info->mi_row_start; mi_row < tile_info->mi_row_end;
       mi_row += cm->seq_params->mib_size) {
    av1_encode_sb_row(cpi, td, tile_row, tile_col, mi_row);
  }
  this_tile->abs_sum_level = td->abs_sum_level;
}

/*!\brief Break one frame into tiles and encode the tiles
 *
 * \ingroup partition_search
 *
 * \param[in]    cpi    Top-level encoder structure
 */
static AOM_INLINE void encode_tiles(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int tile_cols = cm->tiles.cols;
  const int tile_rows = cm->tiles.rows;
  int tile_col, tile_row;

  MACROBLOCK *const mb = &cpi->td.mb;
  assert(IMPLIES(cpi->tile_data == NULL,
                 cpi->allocated_tiles < tile_cols * tile_rows));
  if (cpi->allocated_tiles < tile_cols * tile_rows) av1_alloc_tile_data(cpi);

  av1_init_tile_data(cpi);
  av1_alloc_mb_data(cm, mb, cpi->sf.rt_sf.use_nonrd_pick_mode,
                    cpi->sf.rd_sf.use_mb_rd_hash);

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileDataEnc *const this_tile =
          &cpi->tile_data[tile_row * cm->tiles.cols + tile_col];
      cpi->td.intrabc_used = 0;
      cpi->td.deltaq_used = 0;
      cpi->td.abs_sum_level = 0;
      cpi->td.mb.e_mbd.tile_ctx = &this_tile->tctx;
      cpi->td.mb.tile_pb_ctx = &this_tile->tctx;
      // Reset cyclic refresh counters.
      av1_init_cyclic_refresh_counters(&cpi->td.mb);

      av1_encode_tile(cpi, &cpi->td, tile_row, tile_col);
      // Accumulate cyclic refresh params.
      if (cpi->oxcf.q_cfg.aq_mode == CYCLIC_REFRESH_AQ &&
          !frame_is_intra_only(&cpi->common))
        av1_accumulate_cyclic_refresh_counters(cpi->cyclic_refresh,
                                               &cpi->td.mb);
      cpi->intrabc_used |= cpi->td.intrabc_used;
      cpi->deltaq_used |= cpi->td.deltaq_used;
    }
  }

  av1_dealloc_mb_data(cm, mb);
}

// Set the relative distance of a reference frame w.r.t. current frame
static AOM_INLINE void set_rel_frame_dist(
    const AV1_COMMON *const cm, RefFrameDistanceInfo *const ref_frame_dist_info,
    const int ref_frame_flags) {
  MV_REFERENCE_FRAME ref_frame;
  int min_past_dist = INT32_MAX, min_future_dist = INT32_MAX;
  ref_frame_dist_info->nearest_past_ref = NONE_FRAME;
  ref_frame_dist_info->nearest_future_ref = NONE_FRAME;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    ref_frame_dist_info->ref_relative_dist[ref_frame - LAST_FRAME] = 0;
    if (ref_frame_flags & av1_ref_frame_flag_list[ref_frame]) {
      int dist = av1_encoder_get_relative_dist(
          cm->cur_frame->ref_display_order_hint[ref_frame - LAST_FRAME],
          cm->current_frame.display_order_hint);
      ref_frame_dist_info->ref_relative_dist[ref_frame - LAST_FRAME] = dist;
      // Get the nearest ref_frame in the past
      if (abs(dist) < min_past_dist && dist < 0) {
        ref_frame_dist_info->nearest_past_ref = ref_frame;
        min_past_dist = abs(dist);
      }
      // Get the nearest ref_frame in the future
      if (dist < min_future_dist && dist > 0) {
        ref_frame_dist_info->nearest_future_ref = ref_frame;
        min_future_dist = dist;
      }
    }
  }
}

static INLINE int refs_are_one_sided(const AV1_COMMON *cm) {
  assert(!frame_is_intra_only(cm));

  int one_sided_refs = 1;
  const int cur_display_order_hint = cm->current_frame.display_order_hint;
  for (int ref = LAST_FRAME; ref <= ALTREF_FRAME; ++ref) {
    const RefCntBuffer *const buf = get_ref_frame_buf(cm, ref);
    if (buf == NULL) continue;
    if (av1_encoder_get_relative_dist(buf->display_order_hint,
                                      cur_display_order_hint) > 0) {
      one_sided_refs = 0;  // bwd reference
      break;
    }
  }
  return one_sided_refs;
}

static INLINE void get_skip_mode_ref_offsets(const AV1_COMMON *cm,
                                             int ref_order_hint[2]) {
  const SkipModeInfo *const skip_mode_info = &cm->current_frame.skip_mode_info;
  ref_order_hint[0] = ref_order_hint[1] = 0;
  if (!skip_mode_info->skip_mode_allowed) return;

  const RefCntBuffer *const buf_0 =
      get_ref_frame_buf(cm, LAST_FRAME + skip_mode_info->ref_frame_idx_0);
  const RefCntBuffer *const buf_1 =
      get_ref_frame_buf(cm, LAST_FRAME + skip_mode_info->ref_frame_idx_1);
  assert(buf_0 != NULL && buf_1 != NULL);

  ref_order_hint[0] = buf_0->order_hint;
  ref_order_hint[1] = buf_1->order_hint;
}

static int check_skip_mode_enabled(AV1_COMP *const cpi) {
  AV1_COMMON *const cm = &cpi->common;

  av1_setup_skip_mode_allowed(cm);
  if (!cm->current_frame.skip_mode_info.skip_mode_allowed) return 0;

  // Turn off skip mode if the temporal distances of the reference pair to the
  // current frame are different by more than 1 frame.
  const int cur_offset = (int)cm->current_frame.order_hint;
  int ref_offset[2];
  get_skip_mode_ref_offsets(cm, ref_offset);
  const int cur_to_ref0 = get_relative_dist(&cm->seq_params->order_hint_info,
                                            cur_offset, ref_offset[0]);
  const int cur_to_ref1 = abs(get_relative_dist(
      &cm->seq_params->order_hint_info, cur_offset, ref_offset[1]));
  if (abs(cur_to_ref0 - cur_to_ref1) > 1) return 0;

  // High Latency: Turn off skip mode if all refs are fwd.
  if (cpi->all_one_sided_refs && cpi->oxcf.gf_cfg.lag_in_frames > 0) return 0;

  static const int flag_list[REF_FRAMES] = { 0,
                                             AOM_LAST_FLAG,
                                             AOM_LAST2_FLAG,
                                             AOM_LAST3_FLAG,
                                             AOM_GOLD_FLAG,
                                             AOM_BWD_FLAG,
                                             AOM_ALT2_FLAG,
                                             AOM_ALT_FLAG };
  const int ref_frame[2] = {
    cm->current_frame.skip_mode_info.ref_frame_idx_0 + LAST_FRAME,
    cm->current_frame.skip_mode_info.ref_frame_idx_1 + LAST_FRAME
  };
  if (!(cpi->ref_frame_flags & flag_list[ref_frame[0]]) ||
      !(cpi->ref_frame_flags & flag_list[ref_frame[1]]))
    return 0;

  return 1;
}

static AOM_INLINE void set_default_interp_skip_flags(
    const AV1_COMMON *cm, InterpSearchFlags *interp_search_flags) {
  const int num_planes = av1_num_planes(cm);
  interp_search_flags->default_interp_skip_flags =
      (num_planes == 1) ? INTERP_SKIP_LUMA_EVAL_CHROMA
                        : INTERP_SKIP_LUMA_SKIP_CHROMA;
}

static AOM_INLINE void setup_prune_ref_frame_mask(AV1_COMP *cpi) {
  if ((!cpi->oxcf.ref_frm_cfg.enable_onesided_comp ||
       cpi->sf.inter_sf.disable_onesided_comp) &&
      cpi->all_one_sided_refs) {
    // Disable all compound references
    cpi->prune_ref_frame_mask = (1 << MODE_CTX_REF_FRAMES) - (1 << REF_FRAMES);
  } else if (!cpi->sf.rt_sf.use_nonrd_pick_mode &&
             cpi->sf.inter_sf.selective_ref_frame >= 2) {
    AV1_COMMON *const cm = &cpi->common;
    const int cur_frame_display_order_hint =
        cm->current_frame.display_order_hint;
    unsigned int *ref_display_order_hint =
        cm->cur_frame->ref_display_order_hint;
    const int arf2_dist = av1_encoder_get_relative_dist(
        ref_display_order_hint[ALTREF2_FRAME - LAST_FRAME],
        cur_frame_display_order_hint);
    const int bwd_dist = av1_encoder_get_relative_dist(
        ref_display_order_hint[BWDREF_FRAME - LAST_FRAME],
        cur_frame_display_order_hint);

    for (int ref_idx = REF_FRAMES; ref_idx < MODE_CTX_REF_FRAMES; ++ref_idx) {
      MV_REFERENCE_FRAME rf[2];
      av1_set_ref_frame(rf, ref_idx);
      if (!(cpi->ref_frame_flags & av1_ref_frame_flag_list[rf[0]]) ||
          !(cpi->ref_frame_flags & av1_ref_frame_flag_list[rf[1]])) {
        continue;
      }

      if (!cpi->all_one_sided_refs) {
        int ref_dist[2];
        for (int i = 0; i < 2; ++i) {
          ref_dist[i] = av1_encoder_get_relative_dist(
              ref_display_order_hint[rf[i] - LAST_FRAME],
              cur_frame_display_order_hint);
        }

        // One-sided compound is used only when all reference frames are
        // one-sided.
        if ((ref_dist[0] > 0) == (ref_dist[1] > 0)) {
          cpi->prune_ref_frame_mask |= 1 << ref_idx;
        }
      }

      if (cpi->sf.inter_sf.selective_ref_frame >= 4 &&
          (rf[0] == ALTREF2_FRAME || rf[1] == ALTREF2_FRAME) &&
          (cpi->ref_frame_flags & av1_ref_frame_flag_list[BWDREF_FRAME])) {
        // Check if both ALTREF2_FRAME and BWDREF_FRAME are future references.
        if (arf2_dist > 0 && bwd_dist > 0 && bwd_dist <= arf2_dist) {
          // Drop ALTREF2_FRAME as a reference if BWDREF_FRAME is a closer
          // reference to the current frame than ALTREF2_FRAME
          cpi->prune_ref_frame_mask |= 1 << ref_idx;
        }
      }
    }
  }
}

/*!\brief Encoder setup(only for the current frame), encoding, and recontruction
 * for a single frame
 *
 * \ingroup high_level_algo
 */
static AOM_INLINE void encode_frame_internal(AV1_COMP *cpi) {
  ThreadData *const td = &cpi->td;
  MACROBLOCK *const x = &td->mb;
  AV1_COMMON *const cm = &cpi->common;
  CommonModeInfoParams *const mi_params = &cm->mi_params;
  FeatureFlags *const features = &cm->features;
  MACROBLOCKD *const xd = &x->e_mbd;
  RD_COUNTS *const rdc = &cpi->td.rd_counts;
#if CONFIG_FRAME_PARALLEL_ENCODE && CONFIG_FPMT_TEST
  FrameProbInfo *const temp_frame_probs = &cpi->ppi->temp_frame_probs;
  FrameProbInfo *const temp_frame_probs_simulation =
      &cpi->ppi->temp_frame_probs_simulation;
#endif
  FrameProbInfo *const frame_probs = &cpi->ppi->frame_probs;
  IntraBCHashInfo *const intrabc_hash_info = &x->intrabc_hash_info;
  MultiThreadInfo *const mt_info = &cpi->mt_info;
  AV1EncRowMultiThreadInfo *const enc_row_mt = &mt_info->enc_row_mt;
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;
  const DELTAQ_MODE deltaq_mode = oxcf->q_cfg.deltaq_mode;
  int i;

  if (!cpi->sf.rt_sf.use_nonrd_pick_mode) {
    mi_params->setup_mi(mi_params);
  }

  set_mi_offsets(mi_params, xd, 0, 0);

  av1_zero(*td->counts);
  av1_zero(rdc->tx_type_used);
  av1_zero(rdc->obmc_used);
  av1_zero(rdc->warped_used);

  // Reset the flag.
  cpi->intrabc_used = 0;
  // Need to disable intrabc when superres is selected
  if (av1_superres_scaled(cm)) {
    features->allow_intrabc = 0;
  }

  features->allow_intrabc &= (oxcf->kf_cfg.enable_intrabc);

  if (features->allow_warped_motion &&
      cpi->sf.inter_sf.prune_warped_prob_thresh > 0) {
    const FRAME_UPDATE_TYPE update_type =
        get_frame_update_type(&cpi->ppi->gf_group, cpi->gf_frame_index);
    int warped_probability =
#if CONFIG_FRAME_PARALLEL_ENCODE && CONFIG_FPMT_TEST
        cpi->ppi->fpmt_unit_test_cfg == PARALLEL_SIMULATION_ENCODE
            ? temp_frame_probs->warped_probs[update_type]
            :
#endif  // CONFIG_FRAME_PARALLEL_ENCODE && CONFIG_FPMT_TEST
            frame_probs->warped_probs[update_type];
    if (warped_probability < cpi->sf.inter_sf.prune_warped_prob_thresh)
      features->allow_warped_motion = 0;
  }

  int hash_table_created = 0;
  if (!is_stat_generation_stage(cpi) && av1_use_hash_me(cpi) &&
      !cpi->sf.rt_sf.use_nonrd_pick_mode) {
    // TODO(any): move this outside of the recoding loop to avoid recalculating
    // the hash table.
    // add to hash table
    const int pic_width = cpi->source->y_crop_width;
    const int pic_height = cpi->source->y_crop_height;
    uint32_t *block_hash_values[2][2];
    int8_t *is_block_same[2][3];
    int k, j;

    for (k = 0; k < 2; k++) {
      for (j = 0; j < 2; j++) {
        CHECK_MEM_ERROR(cm, block_hash_values[k][j],
                        aom_malloc(sizeof(uint32_t) * pic_width * pic_height));
      }

      for (j = 0; j < 3; j++) {
        CHECK_MEM_ERROR(cm, is_block_same[k][j],
                        aom_malloc(sizeof(int8_t) * pic_width * pic_height));
      }
    }

    av1_hash_table_init(intrabc_hash_info);
    av1_hash_table_create(&intrabc_hash_info->intrabc_hash_table);
    hash_table_created = 1;
    av1_generate_block_2x2_hash_value(intrabc_hash_info, cpi->source,
                                      block_hash_values[0], is_block_same[0]);
    // Hash data generated for screen contents is used for intraBC ME
    const int min_alloc_size = block_size_wide[mi_params->mi_alloc_bsize];
    const int max_sb_size =
        (1 << (cm->seq_params->mib_size_log2 + MI_SIZE_LOG2));
    int src_idx = 0;
    for (int size = 4; size <= max_sb_size; size *= 2, src_idx = !src_idx) {
      const int dst_idx = !src_idx;
      av1_generate_block_hash_value(
          intrabc_hash_info, cpi->source, size, block_hash_values[src_idx],
          block_hash_values[dst_idx], is_block_same[src_idx],
          is_block_same[dst_idx]);
      if (size >= min_alloc_size) {
        av1_add_to_hash_map_by_row_with_precal_data(
            &intrabc_hash_info->intrabc_hash_table, block_hash_values[dst_idx],
            is_block_same[dst_idx][2], pic_width, pic_height, size);
      }
    }

    for (k = 0; k < 2; k++) {
      for (j = 0; j < 2; j++) {
        aom_free(block_hash_values[k][j]);
      }

      for (j = 0; j < 3; j++) {
        aom_free(is_block_same[k][j]);
      }
    }
  }

  const CommonQuantParams *quant_params = &cm->quant_params;
  for (i = 0; i < MAX_SEGMENTS; ++i) {
    const int qindex =
        cm->seg.enabled ? av1_get_qindex(&cm->seg, i, quant_params->base_qindex)
                        : quant_params->base_qindex;
    xd->lossless[i] =
        qindex == 0 && quant_params->y_dc_delta_q == 0 &&
        quant_params->u_dc_delta_q == 0 && quant_params->u_ac_delta_q == 0 &&
        quant_params->v_dc_delta_q == 0 && quant_params->v_ac_delta_q == 0;
    if (xd->lossless[i]) cpi->enc_seg.has_lossless_segment = 1;
    xd->qindex[i] = qindex;
    if (xd->lossless[i]) {
      cpi->optimize_seg_arr[i] = NO_TRELLIS_OPT;
    } else {
      cpi->optimize_seg_arr[i] = cpi->sf.rd_sf.optimize_coefficients;
    }
  }
  features->coded_lossless = is_coded_lossless(cm, xd);
  features->all_lossless = features->coded_lossless && !av1_superres_scaled(cm);

  // Fix delta q resolution for the moment
  cm->delta_q_info.delta_q_res = 0;
  if (cpi->oxcf.q_cfg.aq_mode != CYCLIC_REFRESH_AQ) {
    if (deltaq_mode == DELTA_Q_OBJECTIVE)
      cm->delta_q_info.delta_q_res = DEFAULT_DELTA_Q_RES_OBJECTIVE;
    else if (deltaq_mode == DELTA_Q_PERCEPTUAL)
      cm->delta_q_info.delta_q_res = DEFAULT_DELTA_Q_RES_PERCEPTUAL;
    else if (deltaq_mode == DELTA_Q_PERCEPTUAL_AI)
      cm->delta_q_info.delta_q_res = DEFAULT_DELTA_Q_RES_PERCEPTUAL;
    else if (deltaq_mode == DELTA_Q_USER_RATING_BASED)
      cm->delta_q_info.delta_q_res = DEFAULT_DELTA_Q_RES_PERCEPTUAL;
    else if (deltaq_mode == DELTA_Q_HDR)
      cm->delta_q_info.delta_q_res = DEFAULT_DELTA_Q_RES_PERCEPTUAL;
    // Set delta_q_present_flag before it is used for the first time
    cm->delta_q_info.delta_lf_res = DEFAULT_DELTA_LF_RES;
    cm->delta_q_info.delta_q_present_flag = deltaq_mode != NO_DELTA_Q;

    // Turn off cm->delta_q_info.delta_q_present_flag if objective delta_q
    // is used for ineligible frames. That effectively will turn off row_mt
    // usage. Note objective delta_q and tpl eligible frames are only altref
    // frames currently.
    const GF_GROUP *gf_group = &cpi->ppi->gf_group;
    if (cm->delta_q_info.delta_q_present_flag) {
      if (deltaq_mode == DELTA_Q_OBJECTIVE &&
          !is_frame_tpl_eligible(gf_group, cpi->gf_frame_index))
        cm->delta_q_info.delta_q_present_flag = 0;
    }

    // Reset delta_q_used flag
    cpi->deltaq_used = 0;

    cm->delta_q_info.delta_lf_present_flag =
        cm->delta_q_info.delta_q_present_flag &&
        oxcf->tool_cfg.enable_deltalf_mode;
    cm->delta_q_info.delta_lf_multi = DEFAULT_DELTA_LF_MULTI;

    // update delta_q_present_flag and delta_lf_present_flag based on
    // base_qindex
    cm->delta_q_info.delta_q_present_flag &= quant_params->base_qindex > 0;
    cm->delta_q_info.delta_lf_present_flag &= quant_params->base_qindex > 0;
  } else {
    cpi->cyclic_refresh->actual_num_seg1_blocks = 0;
    cpi->cyclic_refresh->actual_num_seg2_blocks = 0;
    cpi->cyclic_refresh->cnt_zeromv = 0;
  }

  av1_frame_init_quantizer(cpi);

  init_encode_frame_mb_context(cpi);
  set_default_interp_skip_flags(cm, &cpi->interp_search_flags);
  if (cm->prev_frame && cm->prev_frame->seg.enabled)
    cm->last_frame_seg_map = cm->prev_frame->seg_map;
  else
    cm->last_frame_seg_map = NULL;
  if (features->allow_intrabc || features->coded_lossless) {
    av1_set_default_ref_deltas(cm->lf.ref_deltas);
    av1_set_default_mode_deltas(cm->lf.mode_deltas);
  } else if (cm->prev_frame) {
    memcpy(cm->lf.ref_deltas, cm->prev_frame->ref_deltas, REF_FRAMES);
    memcpy(cm->lf.mode_deltas, cm->prev_frame->mode_deltas, MAX_MODE_LF_DELTAS);
  }
  memcpy(cm->cur_frame->ref_deltas, cm->lf.ref_deltas, REF_FRAMES);
  memcpy(cm->cur_frame->mode_deltas, cm->lf.mode_deltas, MAX_MODE_LF_DELTAS);

  cpi->all_one_sided_refs =
      frame_is_intra_only(cm) ? 0 : refs_are_one_sided(cm);

  cpi->prune_ref_frame_mask = 0;
  // Figure out which ref frames can be skipped at frame level.
  setup_prune_ref_frame_mask(cpi);

  x->txfm_search_info.txb_split_count = 0;
#if CONFIG_SPEED_STATS
  x->txfm_search_info.tx_search_count = 0;
#endif  // CONFIG_SPEED_STATS

#if !CONFIG_REALTIME_ONLY
#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_compute_global_motion_time);
#endif
  av1_compute_global_motion_facade(cpi);
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_compute_global_motion_time);
#endif
#endif  // !CONFIG_REALTIME_ONLY

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_setup_motion_field_time);
#endif
  av1_calculate_ref_frame_side(cm);
  if (features->allow_ref_frame_mvs) av1_setup_motion_field(cm);
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_setup_motion_field_time);
#endif

  cm->current_frame.skip_mode_info.skip_mode_flag =
      check_skip_mode_enabled(cpi);

  // Initialization of skip mode cost depends on the value of
  // 'skip_mode_flag'. This initialization happens in the function
  // av1_fill_mode_rates(), which is in turn called in
  // av1_initialize_rd_consts(). Thus, av1_initialize_rd_consts()
  // has to be called after 'skip_mode_flag' is initialized.
  av1_initialize_rd_consts(cpi);
  av1_set_sad_per_bit(cpi, &x->sadperbit, quant_params->base_qindex);

  enc_row_mt->sync_read_ptr = av1_row_mt_sync_read_dummy;
  enc_row_mt->sync_write_ptr = av1_row_mt_sync_write_dummy;
  mt_info->row_mt_enabled = 0;

  if (oxcf->row_mt && (mt_info->num_workers > 1)) {
    mt_info->row_mt_enabled = 1;
    enc_row_mt->sync_read_ptr = av1_row_mt_sync_read;
    enc_row_mt->sync_write_ptr = av1_row_mt_sync_write;
    av1_encode_tiles_row_mt(cpi);
  } else {
    if (AOMMIN(mt_info->num_workers, cm->tiles.cols * cm->tiles.rows) > 1)
      av1_encode_tiles_mt(cpi);
    else
      encode_tiles(cpi);
  }

  // If intrabc is allowed but never selected, reset the allow_intrabc flag.
  if (features->allow_intrabc && !cpi->intrabc_used) {
    features->allow_intrabc = 0;
  }
  if (features->allow_intrabc) {
    cm->delta_q_info.delta_lf_present_flag = 0;
  }

  if (cm->delta_q_info.delta_q_present_flag && cpi->deltaq_used == 0) {
    cm->delta_q_info.delta_q_present_flag = 0;
  }

  // Set the transform size appropriately before bitstream creation
  const MODE_EVAL_TYPE eval_type =
      cpi->sf.winner_mode_sf.enable_winner_mode_for_tx_size_srch
          ? WINNER_MODE_EVAL
          : DEFAULT_EVAL;
  const TX_SIZE_SEARCH_METHOD tx_search_type =
      cpi->winner_mode_params.tx_size_search_methods[eval_type];
  assert(oxcf->txfm_cfg.enable_tx64 || tx_search_type != USE_LARGESTALL);
  features->tx_mode = select_tx_mode(cm, tx_search_type);

#if CONFIG_FRAME_PARALLEL_ENCODE
  // Retain the frame level probability update conditions for parallel frames.
  // These conditions will be consumed during postencode stage to update the
  // probability.
  if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] > 0) {
    cpi->do_update_frame_probs_txtype[cpi->num_frame_recode] =
        cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats;
    cpi->do_update_frame_probs_obmc[cpi->num_frame_recode] =
        (cpi->sf.inter_sf.prune_obmc_prob_thresh > 0 &&
         cpi->sf.inter_sf.prune_obmc_prob_thresh < INT_MAX);
    cpi->do_update_frame_probs_warp[cpi->num_frame_recode] =
        (features->allow_warped_motion &&
         cpi->sf.inter_sf.prune_warped_prob_thresh > 0);
    cpi->do_update_frame_probs_interpfilter[cpi->num_frame_recode] =
        (cm->current_frame.frame_type != KEY_FRAME &&
         cpi->sf.interp_sf.adaptive_interp_filter_search == 2 &&
         features->interp_filter == SWITCHABLE);
  }
#endif

  if (cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats ||
      ((cpi->sf.tx_sf.tx_type_search.fast_inter_tx_type_prob_thresh !=
        INT_MAX) &&
       (cpi->sf.tx_sf.tx_type_search.fast_inter_tx_type_prob_thresh != 0))) {
    const FRAME_UPDATE_TYPE update_type =
        get_frame_update_type(&cpi->ppi->gf_group, cpi->gf_frame_index);
    for (i = 0; i < TX_SIZES_ALL; i++) {
      int sum = 0;
      int j;
      int left = MAX_TX_TYPE_PROB;

      for (j = 0; j < TX_TYPES; j++)
        sum += cpi->td.rd_counts.tx_type_used[i][j];

      for (j = TX_TYPES - 1; j >= 0; j--) {
        int update_txtype_frameprobs = 1;
        const int new_prob =
            sum ? MAX_TX_TYPE_PROB * cpi->td.rd_counts.tx_type_used[i][j] / sum
                : (j ? 0 : MAX_TX_TYPE_PROB);
#if CONFIG_FRAME_PARALLEL_ENCODE
#if CONFIG_FPMT_TEST
        if (cpi->ppi->fpmt_unit_test_cfg == PARALLEL_SIMULATION_ENCODE) {
          if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] ==
              0) {
            int prob =
                (temp_frame_probs_simulation->tx_type_probs[update_type][i][j] +
                 new_prob) >>
                1;
            left -= prob;
            if (j == 0) prob += left;
            temp_frame_probs_simulation->tx_type_probs[update_type][i][j] =
                prob;
            // Copy temp_frame_probs_simulation to temp_frame_probs
            for (int update_type_idx = 0; update_type_idx < FRAME_UPDATE_TYPES;
                 update_type_idx++) {
              temp_frame_probs->tx_type_probs[update_type_idx][i][j] =
                  temp_frame_probs_simulation
                      ->tx_type_probs[update_type_idx][i][j];
            }
          }
          update_txtype_frameprobs = 0;
        }
#endif  // CONFIG_FPMT_TEST
        // Track the frame probabilities of parallel encode frames to update
        // during postencode stage.
        if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] > 0) {
          update_txtype_frameprobs = 0;
          cpi->frame_new_probs[cpi->num_frame_recode]
              .tx_type_probs[update_type][i][j] = new_prob;
        }
#endif  // CONFIG_FRAME_PARALLEL_ENCODE
        if (update_txtype_frameprobs) {
          int prob =
              (frame_probs->tx_type_probs[update_type][i][j] + new_prob) >> 1;
          left -= prob;
          if (j == 0) prob += left;
          frame_probs->tx_type_probs[update_type][i][j] = prob;
        }
      }
    }
  }

  if (cpi->sf.inter_sf.prune_obmc_prob_thresh > 0 &&
      cpi->sf.inter_sf.prune_obmc_prob_thresh < INT_MAX) {
    const FRAME_UPDATE_TYPE update_type =
        get_frame_update_type(&cpi->ppi->gf_group, cpi->gf_frame_index);

    for (i = 0; i < BLOCK_SIZES_ALL; i++) {
      int sum = 0;
      int update_obmc_frameprobs = 1;
      for (int j = 0; j < 2; j++) sum += cpi->td.rd_counts.obmc_used[i][j];

      const int new_prob =
          sum ? 128 * cpi->td.rd_counts.obmc_used[i][1] / sum : 0;
#if CONFIG_FRAME_PARALLEL_ENCODE
#if CONFIG_FPMT_TEST
      if (cpi->ppi->fpmt_unit_test_cfg == PARALLEL_SIMULATION_ENCODE) {
        if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] == 0) {
          temp_frame_probs_simulation->obmc_probs[update_type][i] =
              (temp_frame_probs_simulation->obmc_probs[update_type][i] +
               new_prob) >>
              1;
          // Copy temp_frame_probs_simulation to temp_frame_probs
          for (int update_type_idx = 0; update_type_idx < FRAME_UPDATE_TYPES;
               update_type_idx++) {
            temp_frame_probs->obmc_probs[update_type_idx][i] =
                temp_frame_probs_simulation->obmc_probs[update_type_idx][i];
          }
        }
        update_obmc_frameprobs = 0;
      }
#endif  // CONFIG_FPMT_TEST
      // Track the frame probabilities of parallel encode frames to update
      // during postencode stage.
      if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] > 0) {
        update_obmc_frameprobs = 0;
        cpi->frame_new_probs[cpi->num_frame_recode].obmc_probs[update_type][i] =
            new_prob;
      }
#endif  // CONFIG_FRAME_PARALLEL_ENCODE
      if (update_obmc_frameprobs) {
        frame_probs->obmc_probs[update_type][i] =
            (frame_probs->obmc_probs[update_type][i] + new_prob) >> 1;
      }
    }
  }

  if (features->allow_warped_motion &&
      cpi->sf.inter_sf.prune_warped_prob_thresh > 0) {
    const FRAME_UPDATE_TYPE update_type =
        get_frame_update_type(&cpi->ppi->gf_group, cpi->gf_frame_index);
    int update_warp_frameprobs = 1;
    int sum = 0;
    for (i = 0; i < 2; i++) sum += cpi->td.rd_counts.warped_used[i];
    const int new_prob = sum ? 128 * cpi->td.rd_counts.warped_used[1] / sum : 0;
#if CONFIG_FRAME_PARALLEL_ENCODE
#if CONFIG_FPMT_TEST
    if (cpi->ppi->fpmt_unit_test_cfg == PARALLEL_SIMULATION_ENCODE) {
      if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] == 0) {
        temp_frame_probs_simulation->warped_probs[update_type] =
            (temp_frame_probs_simulation->warped_probs[update_type] +
             new_prob) >>
            1;
        // Copy temp_frame_probs_simulation to temp_frame_probs
        for (int update_type_idx = 0; update_type_idx < FRAME_UPDATE_TYPES;
             update_type_idx++) {
          temp_frame_probs->warped_probs[update_type_idx] =
              temp_frame_probs_simulation->warped_probs[update_type_idx];
        }
      }
      update_warp_frameprobs = 0;
    }
#endif  // CONFIG_FPMT_TEST
    // Track the frame probabilities of parallel encode frames to update
    // during postencode stage.
    if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] > 0) {
      update_warp_frameprobs = 0;
      cpi->frame_new_probs[cpi->num_frame_recode].warped_probs[update_type] =
          new_prob;
    }
#endif  // CONFIG_FRAME_PARALLEL_ENCODE
    if (update_warp_frameprobs) {
      frame_probs->warped_probs[update_type] =
          (frame_probs->warped_probs[update_type] + new_prob) >> 1;
    }
  }

  if (cm->current_frame.frame_type != KEY_FRAME &&
      cpi->sf.interp_sf.adaptive_interp_filter_search == 2 &&
      features->interp_filter == SWITCHABLE) {
    const FRAME_UPDATE_TYPE update_type =
        get_frame_update_type(&cpi->ppi->gf_group, cpi->gf_frame_index);

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; i++) {
      int sum = 0;
      int j;
      int left = 1536;

      for (j = 0; j < SWITCHABLE_FILTERS; j++) {
        sum += cpi->td.counts->switchable_interp[i][j];
      }

      for (j = SWITCHABLE_FILTERS - 1; j >= 0; j--) {
        int update_interpfilter_frameprobs = 1;
        const int new_prob =
            sum ? 1536 * cpi->td.counts->switchable_interp[i][j] / sum
                : (j ? 0 : 1536);
#if CONFIG_FRAME_PARALLEL_ENCODE
#if CONFIG_FPMT_TEST
        if (cpi->ppi->fpmt_unit_test_cfg == PARALLEL_SIMULATION_ENCODE) {
          if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] ==
              0) {
            int prob = (temp_frame_probs_simulation
                            ->switchable_interp_probs[update_type][i][j] +
                        new_prob) >>
                       1;
            left -= prob;
            if (j == 0) prob += left;
            temp_frame_probs_simulation
                ->switchable_interp_probs[update_type][i][j] = prob;
            // Copy temp_frame_probs_simulation to temp_frame_probs
            for (int update_type_idx = 0; update_type_idx < FRAME_UPDATE_TYPES;
                 update_type_idx++) {
              temp_frame_probs->switchable_interp_probs[update_type_idx][i][j] =
                  temp_frame_probs_simulation
                      ->switchable_interp_probs[update_type_idx][i][j];
            }
          }
          update_interpfilter_frameprobs = 0;
        }
#endif  // CONFIG_FPMT_TEST
        // Track the frame probabilities of parallel encode frames to update
        // during postencode stage.
        if (cpi->ppi->gf_group.frame_parallel_level[cpi->gf_frame_index] > 0) {
          update_interpfilter_frameprobs = 0;
          cpi->frame_new_probs[cpi->num_frame_recode]
              .switchable_interp_probs[update_type][i][j] = new_prob;
        }
#endif  // CONFIG_FRAME_PARALLEL_ENCODE
        if (update_interpfilter_frameprobs) {
          int prob = (frame_probs->switchable_interp_probs[update_type][i][j] +
                      new_prob) >>
                     1;
          left -= prob;
          if (j == 0) prob += left;
          frame_probs->switchable_interp_probs[update_type][i][j] = prob;
        }
      }
    }
  }
  if (hash_table_created) {
    av1_hash_table_destroy(&intrabc_hash_info->intrabc_hash_table);
  }
}

/*!\brief Setup reference frame buffers and encode a frame
 *
 * \ingroup high_level_algo
 * \callgraph
 * \callergraph
 *
 * \param[in]    cpi    Top-level encoder structure
 */
void av1_encode_frame(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  CurrentFrame *const current_frame = &cm->current_frame;
  FeatureFlags *const features = &cm->features;
  const int num_planes = av1_num_planes(cm);
  RD_COUNTS *const rdc = &cpi->td.rd_counts;
  // Indicates whether or not to use a default reduced set for ext-tx
  // rather than the potential full set of 16 transforms
  features->reduced_tx_set_used = cpi->oxcf.txfm_cfg.reduced_tx_type_set;

  // Make sure segment_id is no larger than last_active_segid.
  if (cm->seg.enabled && cm->seg.update_map) {
    const int mi_rows = cm->mi_params.mi_rows;
    const int mi_cols = cm->mi_params.mi_cols;
    const int last_active_segid = cm->seg.last_active_segid;
    uint8_t *map = cpi->enc_seg.map;
    for (int mi_row = 0; mi_row < mi_rows; ++mi_row) {
      for (int mi_col = 0; mi_col < mi_cols; ++mi_col) {
        map[mi_col] = AOMMIN(map[mi_col], last_active_segid);
      }
      map += mi_cols;
    }
  }

  av1_setup_frame_buf_refs(cm);
  enforce_max_ref_frames(cpi, &cpi->ref_frame_flags,
                         cm->cur_frame->ref_display_order_hint,
                         cm->current_frame.display_order_hint);
  set_rel_frame_dist(&cpi->common, &cpi->ref_frame_dist_info,
                     cpi->ref_frame_flags);
  av1_setup_frame_sign_bias(cm);

#if CONFIG_MISMATCH_DEBUG
  mismatch_reset_frame(num_planes);
#else
  (void)num_planes;
#endif

  rdc->newmv_or_intra_blocks = 0;

  if (cpi->sf.hl_sf.frame_parameter_update ||
      cpi->sf.rt_sf.use_comp_ref_nonrd) {
    if (frame_is_intra_only(cm))
      current_frame->reference_mode = SINGLE_REFERENCE;
    else
      current_frame->reference_mode = REFERENCE_MODE_SELECT;

    features->interp_filter = SWITCHABLE;
    if (cm->tiles.large_scale) features->interp_filter = EIGHTTAP_REGULAR;

    features->switchable_motion_mode = 1;

    rdc->compound_ref_used_flag = 0;
    rdc->skip_mode_used_flag = 0;

    encode_frame_internal(cpi);

    if (current_frame->reference_mode == REFERENCE_MODE_SELECT) {
      // Use a flag that includes 4x4 blocks
      if (rdc->compound_ref_used_flag == 0) {
        current_frame->reference_mode = SINGLE_REFERENCE;
#if CONFIG_ENTROPY_STATS
        av1_zero(cpi->td.counts->comp_inter);
#endif  // CONFIG_ENTROPY_STATS
      }
    }
    // Re-check on the skip mode status as reference mode may have been
    // changed.
    SkipModeInfo *const skip_mode_info = &current_frame->skip_mode_info;
    if (frame_is_intra_only(cm) ||
        current_frame->reference_mode == SINGLE_REFERENCE) {
      skip_mode_info->skip_mode_allowed = 0;
      skip_mode_info->skip_mode_flag = 0;
    }
    if (skip_mode_info->skip_mode_flag && rdc->skip_mode_used_flag == 0)
      skip_mode_info->skip_mode_flag = 0;

    if (!cm->tiles.large_scale) {
      if (features->tx_mode == TX_MODE_SELECT &&
          cpi->td.mb.txfm_search_info.txb_split_count == 0)
        features->tx_mode = TX_MODE_LARGEST;
    }
  } else {
    // This is needed if real-time speed setting is changed on the fly
    // from one using compound prediction to one using single reference.
    if (current_frame->reference_mode == REFERENCE_MODE_SELECT)
      current_frame->reference_mode = SINGLE_REFERENCE;
    encode_frame_internal(cpi);
  }
}
