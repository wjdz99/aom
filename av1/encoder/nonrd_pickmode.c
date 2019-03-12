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

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/blend.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/aom_timer.h"
#include "aom_ports/mem.h"
#include "aom_ports/system_state.h"

#include "av1/common/cfl.h"
#include "av1/common/common.h"
#include "av1/common/common_data.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/idct.h"
#include "av1/common/mvref_common.h"
#include "av1/common/obmc.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/pred_common.h"
#include "av1/common/quant_common.h"
#include "av1/common/reconinter.h"
#include "av1/common/reconintra.h"
#include "av1/common/scan.h"
#include "av1/common/seg_common.h"
#include "av1/common/txb_common.h"
#include "av1/common/warped_motion.h"

#include "av1/encoder/aq_variance.h"
#include "av1/encoder/av1_quantize.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/encodemb.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/hybrid_fwd_txfm.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/ml.h"
#include "av1/encoder/palette.h"
#include "av1/encoder/pustats.h"
#include "av1/encoder/random.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/tokenize.h"

extern int g_pick_inter_mode_cnt;
typedef struct {
  uint8_t *data;
  int stride;
  int in_use;
} PRED_BUFFER;

typedef struct {
  PRED_BUFFER *best_pred;
  PREDICTION_MODE best_mode;
  TX_SIZE best_tx_size;
  TX_SIZE best_intra_tx_size;
  MV_REFERENCE_FRAME best_ref_frame;
  MV_REFERENCE_FRAME best_second_ref_frame;
  uint8_t best_mode_skip_txfm;
  InterpFilters best_pred_filter;
} BEST_PICKMODE;

typedef struct {
  MV_REFERENCE_FRAME ref_frame;
  PREDICTION_MODE pred_mode;
} REF_MODE;

#define RT_INTER_MODES 9
static const REF_MODE ref_mode_set[RT_INTER_MODES] = {
  { LAST_FRAME, NEARESTMV },   { LAST_FRAME, NEARMV },
  { LAST_FRAME, NEWMV },       { GOLDEN_FRAME, NEARESTMV },
  { GOLDEN_FRAME, NEARMV },    { GOLDEN_FRAME, NEWMV },
  { ALTREF_FRAME, NEARESTMV }, { ALTREF_FRAME, NEARMV },
  { ALTREF_FRAME, NEWMV }
};

static INLINE void init_best_pickmode(BEST_PICKMODE *bp) {
  bp->best_mode = NEARESTMV;
  bp->best_ref_frame = LAST_FRAME;
  bp->best_tx_size = TX_SIZES;
  bp->best_intra_tx_size = TX_SIZES;
  bp->best_pred_filter = EIGHTTAP_REGULAR;
  bp->best_mode_skip_txfm = 0;
  bp->best_second_ref_frame = NONE_FRAME;
  bp->best_pred = NULL;
}

static int combined_motion_search(AV1_COMP *cpi, MACROBLOCK *x,
                                  BLOCK_SIZE bsize, int mi_row, int mi_col,
                                  int_mv *tmp_mv, int *rate_mv,
                                  int64_t best_rd_sofar, int use_base_mv) {
  MACROBLOCKD *xd = &x->e_mbd;
  const AV1_COMMON *cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MB_MODE_INFO *mi = xd->mi[0];
  struct buf_2d backup_yv12[MAX_MB_PLANE] = { { 0, 0, 0, 0, 0 } };
  int step_param = cpi->sf.mv.reduce_first_step_size;
  const int sadpb = x->sadperbit16;
  MV mvp_full;
  const int ref = mi->ref_frame[0];
  const MV ref_mv = av1_get_ref_mv(x, 0).as_mv;
  MV center_mv;
  int dis;
  int rate_mode;
  const MvLimits tmp_mv_limits = x->mv_limits;
  int rv = 0;
  int cost_list[5];
  int search_subpel = 1;
  const YV12_BUFFER_CONFIG *scaled_ref_frame =
      av1_get_scaled_ref_frame(cpi, ref);

  step_param = AOMMIN(step_param, cpi->mv_step_param);
  if (scaled_ref_frame) {
    int i;
    // Swap out the reference frame for a version that's been scaled to
    // match the resolution of the current frame, allowing the existing
    // motion search code to be used without additional modifications.
    for (i = 0; i < MAX_MB_PLANE; i++) backup_yv12[i] = xd->plane[i].pre[0];
    av1_setup_pre_planes(xd, 0, scaled_ref_frame, mi_row, mi_col, NULL,
                         num_planes);
  }
  av1_set_mv_search_range(&x->mv_limits, &ref_mv);

  // Limit motion vector for large lightning change.
  if (cpi->oxcf.speed > 5 /*&& x->lowvar_highsumdiff*/) {
    x->mv_limits.col_min = AOMMAX(x->mv_limits.col_min, -10);
    x->mv_limits.row_min = AOMMAX(x->mv_limits.row_min, -10);
    x->mv_limits.col_max = AOMMIN(x->mv_limits.col_max, 10);
    x->mv_limits.row_max = AOMMIN(x->mv_limits.row_max, 10);
  }

  mvp_full = x->pred_mv[ref];

  mvp_full.col >>= 3;
  mvp_full.row >>= 3;

  if (!use_base_mv)
    center_mv = ref_mv;
  else
    center_mv = tmp_mv->as_mv;

  av1_full_pixel_search(
      cpi, x, bsize, &mvp_full, step_param, cpi->sf.mv.search_method, 0, sadpb,
      cond_cost_list(cpi, cost_list), &center_mv, INT_MAX, 0,
      (MI_SIZE * mi_col), (MI_SIZE * mi_row), 0, &cpi->ss_cfg[SS_CFG_SRC]);

  x->mv_limits = tmp_mv_limits;
  *tmp_mv = x->best_mv;
  // calculate the bit cost on motion vector
  mvp_full.row = tmp_mv->as_mv.row * 8;
  mvp_full.col = tmp_mv->as_mv.col * 8;

  *rate_mv = av1_mv_bit_cost(&mvp_full, &ref_mv, x->nmv_vec_cost,
                             x->mv_cost_stack, MV_COST_WEIGHT);

  rate_mode = x->newmv_mode_cost[x->mbmi_ext->mode_context[ref]][0];
  rv = !(RDCOST(x->rdmult, (*rate_mv + rate_mode), 0) > best_rd_sofar);

  if (rv && search_subpel) {
    SUBPEL_FORCE_STOP subpel_force_stop = cpi->sf.mv.subpel_force_stop;
    if (use_base_mv) subpel_force_stop = HALF_PEL;
    cpi->find_fractional_mv_step(
        x, cm, mi_row, mi_col, &ref_mv, cpi->common.allow_high_precision_mv,
        x->errorperbit, &cpi->fn_ptr[bsize], subpel_force_stop,
        cpi->sf.mv.subpel_iters_per_step, cond_cost_list(cpi, cost_list),
        x->nmv_vec_cost, x->mv_cost_stack, &dis, &x->pred_sse[ref], NULL, NULL,
        0, 0, 0, 0, 0, 1);
    *tmp_mv = x->best_mv;
    *rate_mv = av1_mv_bit_cost(&tmp_mv->as_mv, &ref_mv, x->nmv_vec_cost,
                               x->mv_cost_stack, MV_COST_WEIGHT);
  }

  if (scaled_ref_frame) {
    int i;
    for (i = 0; i < MAX_MB_PLANE; i++) xd->plane[i].pre[0] = backup_yv12[i];
  }
  return rv;
}

static int search_new_mv(AV1_COMP *cpi, MACROBLOCK *x,
                         int_mv frame_mv[][REF_FRAMES],
                         MV_REFERENCE_FRAME ref_frame, int gf_temporal_ref,
                         BLOCK_SIZE bsize, int mi_row, int mi_col,
                         int best_pred_sad, int *rate_mv,
                         unsigned int best_sse_sofar, RD_STATS *best_rdc) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mi = xd->mi[0];
  AV1_COMMON *cm = &cpi->common;
  (void)best_sse_sofar;
  if (ref_frame > LAST_FRAME && gf_temporal_ref &&
      cpi->oxcf.rc_mode == AOM_CBR) {
    int tmp_sad;
    int dis;
    int cost_list[5] = { INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX };

    if (bsize < BLOCK_16X16) return -1;

    tmp_sad = av1_int_pro_motion_estimation(
        cpi, x, bsize, mi_row, mi_col,
        &x->mbmi_ext->ref_mv_stack[ref_frame][0].this_mv.as_mv);

    if (tmp_sad > x->pred_mv_sad[LAST_FRAME]) return -1;
    if (tmp_sad + (num_pels_log2_lookup[bsize] << 4) > best_pred_sad) return -1;

    frame_mv[NEWMV][ref_frame].as_int = mi->mv[0].as_int;
    MV ref_mv = av1_get_ref_mv(x, 0).as_mv;

    *rate_mv =
        av1_mv_bit_cost(&frame_mv[NEWMV][ref_frame].as_mv, &ref_mv,
                        x->nmv_vec_cost, x->mv_cost_stack, MV_COST_WEIGHT);
    frame_mv[NEWMV][ref_frame].as_mv.row >>= 3;
    frame_mv[NEWMV][ref_frame].as_mv.col >>= 3;

    cpi->find_fractional_mv_step(
        x, cm, mi_row, mi_col, &ref_mv, cm->allow_high_precision_mv,
        x->errorperbit, &cpi->fn_ptr[bsize], cpi->sf.mv.subpel_force_stop,
        cpi->sf.mv.subpel_iters_per_step, cond_cost_list(cpi, cost_list),
        x->nmv_vec_cost, x->mv_cost_stack, &dis, &x->pred_sse[ref_frame], NULL,
        NULL, 0, 0, 0, 0, 0, 1);
  } else if (!combined_motion_search(cpi, x, bsize, mi_row, mi_col,
                                     &frame_mv[NEWMV][ref_frame], rate_mv,
                                     best_rdc->rdcost, 0)) {
    return -1;
  }

  return 0;
}

static INLINE void find_predictors(
    AV1_COMP *cpi, MACROBLOCK *x, MV_REFERENCE_FRAME ref_frame,
    int_mv frame_mv[MB_MODE_COUNT][REF_FRAMES], int const_motion[REF_FRAMES],
    int *ref_frame_skip_mask, const int flag_list[4], TileDataEnc *tile_data,
    int mi_row, int mi_col, struct buf_2d yv12_mb[4][MAX_MB_PLANE],
    BLOCK_SIZE bsize, int force_skip_low_temp_var, int comp_pred_allowed) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  MB_MODE_INFO_EXT *const mbmi_ext = x->mbmi_ext;
  const YV12_BUFFER_CONFIG *yv12 = get_ref_frame_yv12_buf(cm, ref_frame);
  const int num_planes = av1_num_planes(cm);
  (void)tile_data;
  (void)const_motion;
  (void)comp_pred_allowed;

  // TODO(jingning) placeholder for inter-frame non-RD mode decision.
  x->pred_mv_sad[ref_frame] = INT_MAX;
  frame_mv[NEWMV][ref_frame].as_int = INVALID_MV;
  // this needs various further optimizations. to be continued..
  if ((cpi->ref_frame_flags & flag_list[ref_frame]) && (yv12 != NULL)) {
    const struct scale_factors *const sf =
        get_ref_scale_factors_const(cm, ref_frame);
    av1_setup_pred_block(xd, yv12_mb[ref_frame], yv12, mi_row, mi_col, sf, sf,
                         num_planes);
    av1_find_mv_refs(cm, xd, mbmi, ref_frame, mbmi_ext->ref_mv_count,
                     mbmi_ext->ref_mv_stack, NULL, mbmi_ext->global_mvs, mi_row,
                     mi_col, mbmi_ext->mode_context);
    av1_find_best_ref_mvs_from_stack(0, mbmi_ext, ref_frame,
                                     &frame_mv[NEARESTMV][ref_frame],
                                     &frame_mv[NEARMV][ref_frame], 0);
    // Early exit for golden frame if force_skip_low_temp_var is set.
    if (!av1_is_scaled(sf) && bsize >= BLOCK_8X8 &&
        !(force_skip_low_temp_var && ref_frame == GOLDEN_FRAME)) {
      av1_mv_pred(cpi, x, yv12_mb[ref_frame][0].buf, yv12->y_stride, ref_frame,
                  bsize);
    }
  } else {
    *ref_frame_skip_mask |= (1 << ref_frame);
  }
}

static void estimate_ref_frame_costs(
    const AV1_COMMON *cm, const MACROBLOCKD *xd, const MACROBLOCK *x,
    int segment_id, unsigned int *ref_costs_single,
    unsigned int (*ref_costs_comp)[REF_FRAMES]) {
  int seg_ref_active =
      segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME);
  if (seg_ref_active) {
    memset(ref_costs_single, 0, REF_FRAMES * sizeof(*ref_costs_single));
    int ref_frame;
    for (ref_frame = 0; ref_frame < REF_FRAMES; ++ref_frame)
      memset(ref_costs_comp[ref_frame], 0,
             REF_FRAMES * sizeof((*ref_costs_comp)[0]));
  } else {
    int intra_inter_ctx = av1_get_intra_inter_context(xd);
    ref_costs_single[INTRA_FRAME] = x->intra_inter_cost[intra_inter_ctx][0];
    unsigned int base_cost = x->intra_inter_cost[intra_inter_ctx][1];

    for (int i = LAST_FRAME; i <= ALTREF_FRAME; ++i)
      ref_costs_single[i] = base_cost;

    const int ctx_p1 = av1_get_pred_context_single_ref_p1(xd);
    const int ctx_p2 = av1_get_pred_context_single_ref_p2(xd);
    const int ctx_p3 = av1_get_pred_context_single_ref_p3(xd);
    const int ctx_p4 = av1_get_pred_context_single_ref_p4(xd);
    const int ctx_p5 = av1_get_pred_context_single_ref_p5(xd);
    const int ctx_p6 = av1_get_pred_context_single_ref_p6(xd);

    // Determine cost of a single ref frame, where frame types are represented
    // by a tree:
    // Level 0: add cost whether this ref is a forward or backward ref
    ref_costs_single[LAST_FRAME] += x->single_ref_cost[ctx_p1][0][0];
    ref_costs_single[LAST2_FRAME] += x->single_ref_cost[ctx_p1][0][0];
    ref_costs_single[LAST3_FRAME] += x->single_ref_cost[ctx_p1][0][0];
    ref_costs_single[GOLDEN_FRAME] += x->single_ref_cost[ctx_p1][0][0];
    ref_costs_single[BWDREF_FRAME] += x->single_ref_cost[ctx_p1][0][1];
    ref_costs_single[ALTREF2_FRAME] += x->single_ref_cost[ctx_p1][0][1];
    ref_costs_single[ALTREF_FRAME] += x->single_ref_cost[ctx_p1][0][1];

    // Level 1: if this ref is forward ref,
    // add cost whether it is last/last2 or last3/golden
    ref_costs_single[LAST_FRAME] += x->single_ref_cost[ctx_p3][2][0];
    ref_costs_single[LAST2_FRAME] += x->single_ref_cost[ctx_p3][2][0];
    ref_costs_single[LAST3_FRAME] += x->single_ref_cost[ctx_p3][2][1];
    ref_costs_single[GOLDEN_FRAME] += x->single_ref_cost[ctx_p3][2][1];

    // Level 1: if this ref is backward ref
    // then add cost whether this ref is altref or backward ref
    ref_costs_single[BWDREF_FRAME] += x->single_ref_cost[ctx_p2][1][0];
    ref_costs_single[ALTREF2_FRAME] += x->single_ref_cost[ctx_p2][1][0];
    ref_costs_single[ALTREF_FRAME] += x->single_ref_cost[ctx_p2][1][1];

    // Level 2: further add cost whether this ref is last or last2
    ref_costs_single[LAST_FRAME] += x->single_ref_cost[ctx_p4][3][0];
    ref_costs_single[LAST2_FRAME] += x->single_ref_cost[ctx_p4][3][1];

    // Level 2: last3 or golden
    ref_costs_single[LAST3_FRAME] += x->single_ref_cost[ctx_p5][4][0];
    ref_costs_single[GOLDEN_FRAME] += x->single_ref_cost[ctx_p5][4][1];

    // Level 2: bwdref or altref2
    ref_costs_single[BWDREF_FRAME] += x->single_ref_cost[ctx_p6][5][0];
    ref_costs_single[ALTREF2_FRAME] += x->single_ref_cost[ctx_p6][5][1];

    if (cm->current_frame.reference_mode != SINGLE_REFERENCE) {
      // Similar to single ref, determine cost of compound ref frames.
      // cost_compound_refs = cost_first_ref + cost_second_ref
      const int bwdref_comp_ctx_p = av1_get_pred_context_comp_bwdref_p(xd);
      const int bwdref_comp_ctx_p1 = av1_get_pred_context_comp_bwdref_p1(xd);
      const int ref_comp_ctx_p = av1_get_pred_context_comp_ref_p(xd);
      const int ref_comp_ctx_p1 = av1_get_pred_context_comp_ref_p1(xd);
      const int ref_comp_ctx_p2 = av1_get_pred_context_comp_ref_p2(xd);

      const int comp_ref_type_ctx = av1_get_comp_reference_type_context(xd);
      unsigned int ref_bicomp_costs[REF_FRAMES] = { 0 };

      ref_bicomp_costs[LAST_FRAME] = ref_bicomp_costs[LAST2_FRAME] =
          ref_bicomp_costs[LAST3_FRAME] = ref_bicomp_costs[GOLDEN_FRAME] =
              base_cost + x->comp_ref_type_cost[comp_ref_type_ctx][1];
      ref_bicomp_costs[BWDREF_FRAME] = ref_bicomp_costs[ALTREF2_FRAME] = 0;
      ref_bicomp_costs[ALTREF_FRAME] = 0;

      // cost of first ref frame
      ref_bicomp_costs[LAST_FRAME] += x->comp_ref_cost[ref_comp_ctx_p][0][0];
      ref_bicomp_costs[LAST2_FRAME] += x->comp_ref_cost[ref_comp_ctx_p][0][0];
      ref_bicomp_costs[LAST3_FRAME] += x->comp_ref_cost[ref_comp_ctx_p][0][1];
      ref_bicomp_costs[GOLDEN_FRAME] += x->comp_ref_cost[ref_comp_ctx_p][0][1];

      ref_bicomp_costs[LAST_FRAME] += x->comp_ref_cost[ref_comp_ctx_p1][1][0];
      ref_bicomp_costs[LAST2_FRAME] += x->comp_ref_cost[ref_comp_ctx_p1][1][1];

      ref_bicomp_costs[LAST3_FRAME] += x->comp_ref_cost[ref_comp_ctx_p2][2][0];
      ref_bicomp_costs[GOLDEN_FRAME] += x->comp_ref_cost[ref_comp_ctx_p2][2][1];

      // cost of second ref frame
      ref_bicomp_costs[BWDREF_FRAME] +=
          x->comp_bwdref_cost[bwdref_comp_ctx_p][0][0];
      ref_bicomp_costs[ALTREF2_FRAME] +=
          x->comp_bwdref_cost[bwdref_comp_ctx_p][0][0];
      ref_bicomp_costs[ALTREF_FRAME] +=
          x->comp_bwdref_cost[bwdref_comp_ctx_p][0][1];

      ref_bicomp_costs[BWDREF_FRAME] +=
          x->comp_bwdref_cost[bwdref_comp_ctx_p1][1][0];
      ref_bicomp_costs[ALTREF2_FRAME] +=
          x->comp_bwdref_cost[bwdref_comp_ctx_p1][1][1];

      // cost: if one ref frame is forward ref, the other ref is backward ref
      int ref0, ref1;
      for (ref0 = LAST_FRAME; ref0 <= GOLDEN_FRAME; ++ref0) {
        for (ref1 = BWDREF_FRAME; ref1 <= ALTREF_FRAME; ++ref1) {
          ref_costs_comp[ref0][ref1] =
              ref_bicomp_costs[ref0] + ref_bicomp_costs[ref1];
        }
      }

      // cost: if both ref frames are the same side.
      const int uni_comp_ref_ctx_p = av1_get_pred_context_uni_comp_ref_p(xd);
      const int uni_comp_ref_ctx_p1 = av1_get_pred_context_uni_comp_ref_p1(xd);
      const int uni_comp_ref_ctx_p2 = av1_get_pred_context_uni_comp_ref_p2(xd);
      ref_costs_comp[LAST_FRAME][LAST2_FRAME] =
          base_cost + x->comp_ref_type_cost[comp_ref_type_ctx][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p][0][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p1][1][0];
      ref_costs_comp[LAST_FRAME][LAST3_FRAME] =
          base_cost + x->comp_ref_type_cost[comp_ref_type_ctx][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p][0][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p1][1][1] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p2][2][0];
      ref_costs_comp[LAST_FRAME][GOLDEN_FRAME] =
          base_cost + x->comp_ref_type_cost[comp_ref_type_ctx][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p][0][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p1][1][1] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p2][2][1];
      ref_costs_comp[BWDREF_FRAME][ALTREF_FRAME] =
          base_cost + x->comp_ref_type_cost[comp_ref_type_ctx][0] +
          x->uni_comp_ref_cost[uni_comp_ref_ctx_p][0][1];
    } else {
      int ref0, ref1;
      for (ref0 = LAST_FRAME; ref0 <= GOLDEN_FRAME; ++ref0) {
        for (ref1 = BWDREF_FRAME; ref1 <= ALTREF_FRAME; ++ref1)
          ref_costs_comp[ref0][ref1] = 512;
      }
      ref_costs_comp[LAST_FRAME][LAST2_FRAME] = 512;
      ref_costs_comp[LAST_FRAME][LAST3_FRAME] = 512;
      ref_costs_comp[LAST_FRAME][GOLDEN_FRAME] = 512;
      ref_costs_comp[BWDREF_FRAME][ALTREF_FRAME] = 512;
    }
  }
}

enum {
  //  INTER_ALL = (1 << NEARESTMV) | (1 << NEARMV) | (1 << NEWMV),
  INTER_NEAREST = (1 << NEARESTMV),
  INTER_NEAREST_NEW = (1 << NEARESTMV) | (1 << NEWMV),
  INTER_NEAREST_NEAR = (1 << NEARESTMV) | (1 << NEARMV),
  INTER_NEAR_NEW = (1 << NEARMV) | (1 << NEWMV),
};

static void model_rd_from_sse(const AV1_COMP *const cpi,
                              const MACROBLOCK *const x, BLOCK_SIZE plane_bsize,
                              int plane, int64_t sse, int num_samples,
                              int *rate, int64_t *dist) {
  (void)num_samples;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const int dequant_shift = (is_cur_buf_hbd(xd)) ? xd->bd - 5 : 3;

  // Fast approximate the modelling function.
  if (cpi->sf.simple_model_rd_from_var) {
    const int64_t square_error = sse;
    int quantizer = pd->dequant_Q3[1] >> dequant_shift;
    if (quantizer < 120)
      *rate = (int)AOMMIN(
          (square_error * (280 - quantizer)) >> (16 - AV1_PROB_COST_SHIFT),
          INT_MAX);
    else
      *rate = 0;
    assert(*rate >= 0);
    *dist = (square_error * quantizer) >> 8;
  } else {
    av1_model_rd_from_var_lapndz(sse, num_pels_log2_lookup[plane_bsize],
                                 pd->dequant_Q3[1] >> dequant_shift, rate,
                                 dist);
  }
  *dist <<= 4;
}

static void model_rd_for_sb(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                            MACROBLOCK *x, MACROBLOCKD *xd, int plane_from,
                            int plane_to, int mi_row, int mi_col,
                            int *out_rate_sum, int64_t *out_dist_sum,
                            int *skip_txfm_sb, int64_t *skip_sse_sb,
                            int *plane_rate, int64_t *plane_sse,
                            int64_t *plane_dist) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  int plane;
  (void)mi_row;
  (void)mi_col;
  const int ref = xd->mi[0]->ref_frame[0];

  int64_t rate_sum = 0;
  int64_t dist_sum = 0;
  int64_t total_sse = 0;

  for (plane = plane_from; plane <= plane_to; ++plane) {
    struct macroblock_plane *const p = &x->plane[plane];
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    const int bw = block_size_wide[plane_bsize];
    const int bh = block_size_high[plane_bsize];
    int64_t sse;
    int rate;
    int64_t dist;

    if (x->skip_chroma_rd && plane) continue;

    if (is_cur_buf_hbd(xd)) {
      sse = aom_highbd_sse(p->src.buf, p->src.stride, pd->dst.buf,
                           pd->dst.stride, bw, bh);
    } else {
      sse = aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw,
                    bh);
    }
    sse = ROUND_POWER_OF_TWO(sse, (xd->bd - 8) * 2);

    model_rd_from_sse(cpi, x, plane_bsize, plane, sse, bw * bh, &rate, &dist);

    if (plane == 0) x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

    total_sse += sse;
    rate_sum += rate;
    dist_sum += dist;
    if (plane_rate) plane_rate[plane] = rate;
    if (plane_sse) plane_sse[plane] = sse;
    if (plane_dist) plane_dist[plane] = dist;
    assert(rate_sum >= 0);
  }

  if (skip_txfm_sb) *skip_txfm_sb = total_sse == 0;
  if (skip_sse_sb) *skip_sse_sb = total_sse << 4;
  rate_sum = AOMMIN(rate_sum, INT_MAX);
  *out_rate_sum = (int)rate_sum;
  *out_dist_sum = dist_sum;
}
// The factor to scale from cost in bits to cost in vp9_prob_cost units.
#define AV1_PROB_COST_SHIFT 9

static void block_yrd(AV1_COMP *cpi, MACROBLOCK *x, int mi_row, int mi_col,
                      RD_STATS *this_rdc, int *skippable, int64_t *sse,
                      BLOCK_SIZE bsize, TX_SIZE tx_size, int rd_computed) {
  MACROBLOCKD *xd = &x->e_mbd;
  AV1_COMMON *const cm = &cpi->common;
  const struct macroblockd_plane *pd = &xd->plane[0];
  struct macroblock_plane *const p = &x->plane[0];
  const int num_4x4_w = mi_size_wide[bsize];
  const int num_4x4_h = mi_size_high[bsize];
  const int step = 1 << (tx_size << 1);
  const int block_step = (1 << tx_size);
  int block = 0, r, c;
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >> 5);
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 : xd->mb_to_bottom_edge >> 5);
  int eob_cost = 0;
  const int bw = 4 * num_4x4_w;
  const int bh = 4 * num_4x4_h;
  const int planes = av1_num_planes(cm);

  assert(tx_size > 0 && tx_size <= 4);
  (void)mi_row;
  (void)mi_col;
  (void)rd_computed;
  /*
    if (cpi->common.current_frame.frame_type != KEY_FRAME &&  bsize <
    BLOCK_32X32) { (void)tx_size; if (!rd_computed) model_rd_for_sb(cpi, bsize,
    x, xd, AOM_PLANE_Y, AOM_PLANE_Y, mi_row, mi_col, &this_rdc->rate,
                      &this_rdc->dist, &this_rdc->skip, NULL, NULL, NULL, NULL);
      *sse = INT_MAX;
      *skippable = 0;
      return;
    }
  */
  (void)cpi;

  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
  *skippable = 1;
  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (r = 0; r < max_blocks_high; r += block_step) {
    for (c = 0; c < num_4x4_w; c += block_step) {
      if (c < max_blocks_wide) {
        const SCAN_ORDER *const scan_order = &av1_default_scan_orders[tx_size];
        tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
        tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
        tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
        uint16_t *const eob = &p->eobs[block];
        const int diff_stride = bw;
        const int16_t *src_diff;
        src_diff = &p->src_diff[(r * diff_stride + c) << 2];

        switch (tx_size) {
          case TX_64X64:
            assert(0);
            // aom_hadamard_64x64(src_diff, diff_stride, coeff);
            // av1_quantize_fp(coeff, 64*64, p->zbin_QTX, p->round_fp_QTX,
            //             p->quant_fp_QTX, p->quant_shift_QTX, qcoeff, dqcoeff,
            //             p->dequant_QTX, eob, scan_order->scan,
            //             scan_order->iscan);
            break;
          case TX_32X32:
            aom_hadamard_32x32(src_diff, diff_stride, coeff);
            av1_quantize_fp(coeff, 32 * 32, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob, scan_order->scan,
                            scan_order->iscan);
            break;
          case TX_16X16:
            aom_hadamard_16x16(src_diff, diff_stride, coeff);
            av1_quantize_fp(coeff, 16 * 16, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob, scan_order->scan,
                            scan_order->iscan);
            break;
          case TX_8X8:
            aom_hadamard_8x8(src_diff, diff_stride, coeff);
            av1_quantize_fp(coeff, 8 * 8, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob, scan_order->scan,
                            scan_order->iscan);
            break;
          default: assert(0); break;
        }
        *skippable &= (*eob == 0);
        eob_cost += 1;
        const int blk_idx =
            r * (block_size_wide[bsize] >> tx_size_wide_log2[0]) + c;
        for (int i = 0; i < planes; i++)
          set_blk_skip(x, i, blk_idx,
                       *eob == 0);  // this is wrong. need extra care about UV
      }
      block += step;
    }
  }

  this_rdc->rate = 0;
  if (*sse < INT64_MAX) {
    *sse = (*sse << 6) >> 2;
    if (*skippable) {
      this_rdc->dist = *sse;
      return;
    }
  }

  block = 0;
  this_rdc->dist = 0;
  for (r = 0; r < max_blocks_high; r += block_step) {
    for (c = 0; c < num_4x4_w; c += block_step) {
      if (c < max_blocks_wide) {
        int64_t dummy;
        tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
        tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
        tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
        uint16_t *const eob = &p->eobs[block];

        if (*eob == 1)
          this_rdc->rate += (int)abs(qcoeff[0]);
        else if (*eob > 1)
          this_rdc->rate += aom_satd(qcoeff, step << 4);

        this_rdc->dist +=
            av1_block_error(coeff, dqcoeff, step << 4, &dummy) >> 2;
      }
      block += step;
    }
  }

  // If skippable is set, rate gets clobbered later.
  this_rdc->rate <<= (2 + AV1_PROB_COST_SHIFT);
  this_rdc->rate += (eob_cost << AV1_PROB_COST_SHIFT);
}

static const THR_MODES mode_idx[REF_FRAMES][4] = {
  { THR_DC, THR_V_PRED, THR_H_PRED, THR_SMOOTH },
  { THR_NEARESTMV, THR_NEARMV, THR_GLOBALMV, THR_NEWMV },
  { THR_NEARESTG, THR_NEARG, THR_GLOBALMV, THR_NEWG },
  { THR_NEARESTA, THR_NEARA, THR_GLOBALMV, THR_NEWA },
};

typedef struct {
  PREDICTION_MODE mode;
  MV_REFERENCE_FRAME ref_frame[2];
} MODE_DEFINITION;

static const MODE_DEFINITION av1_mode_order[MAX_MODES] = {
  { NEARESTMV, { LAST_FRAME, NONE_FRAME } },
  { NEARESTMV, { LAST2_FRAME, NONE_FRAME } },
  { NEARESTMV, { LAST3_FRAME, NONE_FRAME } },
  { NEARESTMV, { BWDREF_FRAME, NONE_FRAME } },
  { NEARESTMV, { ALTREF2_FRAME, NONE_FRAME } },
  { NEARESTMV, { ALTREF_FRAME, NONE_FRAME } },
  { NEARESTMV, { GOLDEN_FRAME, NONE_FRAME } },

  { NEWMV, { LAST_FRAME, NONE_FRAME } },
  { NEWMV, { LAST2_FRAME, NONE_FRAME } },
  { NEWMV, { LAST3_FRAME, NONE_FRAME } },
  { NEWMV, { BWDREF_FRAME, NONE_FRAME } },
  { NEWMV, { ALTREF2_FRAME, NONE_FRAME } },
  { NEWMV, { ALTREF_FRAME, NONE_FRAME } },
  { NEWMV, { GOLDEN_FRAME, NONE_FRAME } },

  { NEARMV, { LAST_FRAME, NONE_FRAME } },
  { NEARMV, { LAST2_FRAME, NONE_FRAME } },
  { NEARMV, { LAST3_FRAME, NONE_FRAME } },
  { NEARMV, { BWDREF_FRAME, NONE_FRAME } },
  { NEARMV, { ALTREF2_FRAME, NONE_FRAME } },
  { NEARMV, { ALTREF_FRAME, NONE_FRAME } },
  { NEARMV, { GOLDEN_FRAME, NONE_FRAME } },

  { GLOBALMV, { LAST_FRAME, NONE_FRAME } },
  { GLOBALMV, { LAST2_FRAME, NONE_FRAME } },
  { GLOBALMV, { LAST3_FRAME, NONE_FRAME } },
  { GLOBALMV, { BWDREF_FRAME, NONE_FRAME } },
  { GLOBALMV, { ALTREF2_FRAME, NONE_FRAME } },
  { GLOBALMV, { GOLDEN_FRAME, NONE_FRAME } },
  { GLOBALMV, { ALTREF_FRAME, NONE_FRAME } },

  // TODO(zoeliu): May need to reconsider the order on the modes to check

  { NEAREST_NEARESTMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEAREST_NEARESTMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEAREST_NEARESTMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEAREST_NEARESTMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEARESTMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEARESTMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEARESTMV, { GOLDEN_FRAME, ALTREF2_FRAME } },

  { NEAREST_NEARESTMV, { LAST_FRAME, LAST2_FRAME } },
  { NEAREST_NEARESTMV, { LAST_FRAME, LAST3_FRAME } },
  { NEAREST_NEARESTMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEAREST_NEARESTMV, { BWDREF_FRAME, ALTREF_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, ALTREF_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, ALTREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, ALTREF_FRAME } },

  { NEAR_NEARMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEW_NEARESTMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEAREST_NEWMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEW_NEARMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEAR_NEWMV, { LAST2_FRAME, ALTREF_FRAME } },
  { NEW_NEWMV, { LAST2_FRAME, ALTREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST2_FRAME, ALTREF_FRAME } },

  { NEAR_NEARMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEW_NEARESTMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEAREST_NEWMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEW_NEARMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEAR_NEWMV, { LAST3_FRAME, ALTREF_FRAME } },
  { NEW_NEWMV, { LAST3_FRAME, ALTREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST3_FRAME, ALTREF_FRAME } },

  { NEAR_NEARMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEW_NEARESTMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEAREST_NEWMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEW_NEARMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEAR_NEWMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { NEW_NEWMV, { GOLDEN_FRAME, ALTREF_FRAME } },
  { GLOBAL_GLOBALMV, { GOLDEN_FRAME, ALTREF_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, BWDREF_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, BWDREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, BWDREF_FRAME } },

  { NEAR_NEARMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEW_NEARESTMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEAREST_NEWMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEW_NEARMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEAR_NEWMV, { LAST2_FRAME, BWDREF_FRAME } },
  { NEW_NEWMV, { LAST2_FRAME, BWDREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST2_FRAME, BWDREF_FRAME } },

  { NEAR_NEARMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEW_NEARESTMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEAREST_NEWMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEW_NEARMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEAR_NEWMV, { LAST3_FRAME, BWDREF_FRAME } },
  { NEW_NEWMV, { LAST3_FRAME, BWDREF_FRAME } },
  { GLOBAL_GLOBALMV, { LAST3_FRAME, BWDREF_FRAME } },

  { NEAR_NEARMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEW_NEARESTMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEAREST_NEWMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEW_NEARMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEAR_NEWMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { NEW_NEWMV, { GOLDEN_FRAME, BWDREF_FRAME } },
  { GLOBAL_GLOBALMV, { GOLDEN_FRAME, BWDREF_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, ALTREF2_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, ALTREF2_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, ALTREF2_FRAME } },

  { NEAR_NEARMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEW_NEARESTMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEWMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEW_NEARMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEAR_NEWMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { NEW_NEWMV, { LAST2_FRAME, ALTREF2_FRAME } },
  { GLOBAL_GLOBALMV, { LAST2_FRAME, ALTREF2_FRAME } },

  { NEAR_NEARMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEW_NEARESTMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEWMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEW_NEARMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEAR_NEWMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { NEW_NEWMV, { LAST3_FRAME, ALTREF2_FRAME } },
  { GLOBAL_GLOBALMV, { LAST3_FRAME, ALTREF2_FRAME } },

  { NEAR_NEARMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { NEW_NEARESTMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { NEAREST_NEWMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { NEW_NEARMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { NEAR_NEWMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { NEW_NEWMV, { GOLDEN_FRAME, ALTREF2_FRAME } },
  { GLOBAL_GLOBALMV, { GOLDEN_FRAME, ALTREF2_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, LAST2_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, LAST2_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, LAST2_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, LAST2_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, LAST2_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, LAST2_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, LAST2_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, LAST3_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, LAST3_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, LAST3_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, LAST3_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, LAST3_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, LAST3_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, LAST3_FRAME } },

  { NEAR_NEARMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEW_NEARESTMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEAREST_NEWMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEW_NEARMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEAR_NEWMV, { LAST_FRAME, GOLDEN_FRAME } },
  { NEW_NEWMV, { LAST_FRAME, GOLDEN_FRAME } },
  { GLOBAL_GLOBALMV, { LAST_FRAME, GOLDEN_FRAME } },

  { NEAR_NEARMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { NEW_NEARESTMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { NEAREST_NEWMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { NEW_NEARMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { NEAR_NEWMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { NEW_NEWMV, { BWDREF_FRAME, ALTREF_FRAME } },
  { GLOBAL_GLOBALMV, { BWDREF_FRAME, ALTREF_FRAME } },

  // intra modes
  { DC_PRED, { INTRA_FRAME, NONE_FRAME } },
  { PAETH_PRED, { INTRA_FRAME, NONE_FRAME } },
  { SMOOTH_PRED, { INTRA_FRAME, NONE_FRAME } },
  { SMOOTH_V_PRED, { INTRA_FRAME, NONE_FRAME } },
  { SMOOTH_H_PRED, { INTRA_FRAME, NONE_FRAME } },
  { H_PRED, { INTRA_FRAME, NONE_FRAME } },
  { V_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D135_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D203_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D157_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D67_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D113_PRED, { INTRA_FRAME, NONE_FRAME } },
  { D45_PRED, { INTRA_FRAME, NONE_FRAME } },
};

static INLINE void init_mbmi(MB_MODE_INFO *mbmi, int mode_index,
                             const AV1_COMMON *cm) {
  PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  PREDICTION_MODE this_mode = av1_mode_order[mode_index].mode;
  mbmi->ref_mv_idx = 0;
  mbmi->mode = this_mode;
  mbmi->uv_mode = UV_DC_PRED;
  mbmi->ref_frame[0] = av1_mode_order[mode_index].ref_frame[0];
  mbmi->ref_frame[1] = av1_mode_order[mode_index].ref_frame[1];
  pmi->palette_size[0] = 0;
  pmi->palette_size[1] = 0;
  mbmi->filter_intra_mode_info.use_filter_intra = 0;
  mbmi->mv[0].as_int = mbmi->mv[1].as_int = 0;
  mbmi->motion_mode = SIMPLE_TRANSLATION;
  mbmi->interintra_mode = (INTERINTRA_MODE)(II_DC_PRED - 1);
  set_default_interp_filters(mbmi, cm->interp_filter);
}

static void store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                                 int mode_index, int skippable) {
  MACROBLOCKD *const xd = &x->e_mbd;

  // Take a snapshot of the coding context so it can be
  // restored if we decide to encode this way
  ctx->skip = x->skip;
  memcpy(ctx->blk_skip, x->blk_skip, sizeof(x->blk_skip[0]) * ctx->num_4x4_blk);
  ctx->skippable = skippable;
  ctx->best_mode_index = mode_index;
  ctx->mic = *xd->mi[0];
  ctx->mbmi_ext = *x->mbmi_ext;
}

static INLINE int bsize_to_num_blk(BLOCK_SIZE bsize) {
  int num_blk = 1 << (num_pels_log2_lookup[bsize] - 2 * tx_size_wide_log2[0]);
  return num_blk;
}

// Used to set proper context for early termination with skip = 1.
static void set_skip_flag(MACROBLOCK *x, RD_STATS *rd_stats, int bsize,
                          int64_t dist) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int n4 = bsize_to_num_blk(bsize);
  const TX_SIZE tx_size = max_txsize_rect_lookup[bsize];
  memset(mbmi->txk_type, DCT_DCT, sizeof(mbmi->txk_type[0]) * TXK_TYPE_BUF_LEN);
  memset(mbmi->inter_tx_size, tx_size, sizeof(mbmi->inter_tx_size));
  mbmi->tx_size = tx_size;
  for (int i = 0; i < n4; ++i) {
    set_blk_skip(x, 0, i, 1);
    set_blk_skip(x, 1, i, 1);
    set_blk_skip(x, 2, i, 1);
  }
  rd_stats->skip = 1;
  if (is_cur_buf_hbd(xd)) dist = ROUND_POWER_OF_TWO(dist, (xd->bd - 8) * 2);
  rd_stats->dist = rd_stats->sse = (dist << 4);
  // Though decision is to make the block as skip based on luma stats,
  // it is possible that block becomes non skip after chroma rd. In addition
  // intermediate non skip costs calculated by caller function will be
  // incorrect, if rate is set as  zero (i.e., if zero_blk_rate is not
  // accounted). Hence intermediate rate is populated to code the luma tx blks
  // as skip, the caller function based on final rd decision (i.e., skip vs
  // non-skip) sets the final rate accordingly. Here the rate populated
  // corresponds to coding all the tx blocks with zero_blk_rate (based on max tx
  // size possible) in the current block. Eg: For 128*128 block, rate would be
  // 4 * zero_blk_rate where zero_blk_rate corresponds to coding of one 64x64 tx
  // block as 'all zeros'
  ENTROPY_CONTEXT ctxa[MAX_MIB_SIZE];
  ENTROPY_CONTEXT ctxl[MAX_MIB_SIZE];
  av1_get_entropy_contexts(bsize, &xd->plane[0], ctxa, ctxl);
  ENTROPY_CONTEXT *ta = ctxa;
  ENTROPY_CONTEXT *tl = ctxl;
  const TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);
  TXB_CTX txb_ctx;
  get_txb_ctx(bsize, tx_size, 0, ta, tl, &txb_ctx);
  const int zero_blk_rate = x->coeff_costs[txs_ctx][PLANE_TYPE_Y]
                                .txb_skip_cost[txb_ctx.txb_skip_ctx][1];
  rd_stats->rate = zero_blk_rate *
                   (block_size_wide[bsize] >> tx_size_wide_log2[tx_size]) *
                   (block_size_high[bsize] >> tx_size_high_log2[tx_size]);
}

void av1_nonrd_pick_inter_mode_sb(AV1_COMP *cpi, TileDataEnc *tile_data,
                                  MACROBLOCK *x, int mi_row, int mi_col,
                                  RD_STATS *rd_cost, BLOCK_SIZE bsize,
                                  PICK_MODE_CONTEXT *ctx,
                                  int64_t best_rd_so_far) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mi = xd->mi[0];

  BEST_PICKMODE best_pickmode;
  int inter_mode_mask[BLOCK_SIZES];

  MV_REFERENCE_FRAME ref_frame;
  MV_REFERENCE_FRAME usable_ref_frame, second_ref_frame;
  int_mv frame_mv[MB_MODE_COUNT][REF_FRAMES];
  uint8_t mode_checked[MB_MODE_COUNT][REF_FRAMES];
  struct buf_2d yv12_mb[4][MAX_MB_PLANE];
  static const int flag_list[5] = { 0, AOM_LAST_FLAG, AOM_LAST2_FLAG,
                                    AOM_LAST3_FLAG, AOM_GOLD_FLAG };
  RD_STATS this_rdc, best_rdc;
  // var_y and sse_y are saved to be used in skipping checking
  unsigned int sse_y = UINT_MAX;
  const int *const rd_threshes = cpi->rd.threshes[mi->segment_id][bsize];
  const int *const rd_thresh_freq_fact = tile_data->thresh_freq_fact[bsize];

  InterpFilter filter_ref;
  int const_motion[REF_FRAMES] = { 0 };
  int ref_frame_skip_mask = 0;
  int idx;
  int best_pred_sad = INT_MAX;
  int best_early_term = 0;
  unsigned int ref_costs_single[REF_FRAMES],
      ref_costs_comp[REF_FRAMES][REF_FRAMES];
  int use_golden_nonzeromv = 1;
  int force_skip_low_temp_var = 0;
  int skip_ref_find_pred[5] = { 0 };
  unsigned int sse_zeromv_normalized = UINT_MAX;
  unsigned int best_sse_sofar = UINT_MAX;
  int gf_temporal_ref = 0;
  int force_test_gf_zeromv = 0;
  const struct segmentation *const seg = &cm->seg;
  int comp_modes = 0;
  int num_inter_modes = RT_INTER_MODES;
  unsigned int thresh_skip_golden = (bsize >= BLOCK_32X32) ? 2500 : 500;
  unsigned char segment_id = mi->segment_id;

  (void)best_rd_so_far;

  init_best_pickmode(&best_pickmode);

  for (int i = 0; i < BLOCK_SIZES; ++i) inter_mode_mask[i] = INTER_ALL;
  inter_mode_mask[BLOCK_32X32] = INTER_NEAREST_NEW;
  inter_mode_mask[BLOCK_32X64] = INTER_NEAREST_NEW;
  inter_mode_mask[BLOCK_64X32] = INTER_NEAREST_NEW;
  inter_mode_mask[BLOCK_64X64] = INTER_NEAREST_NEW;
  inter_mode_mask[BLOCK_128X128] = INTER_NEAREST_NEAR;

  x->source_variance = UINT_MAX;

  struct scale_factors *const sf_last = get_ref_scale_factors(cm, LAST_FRAME);
  struct scale_factors *const sf_golden =
      get_ref_scale_factors(cm, GOLDEN_FRAME);
  gf_temporal_ref = 1;
  // For temporal long term prediction, check that the golden reference
  // is same scale as last reference, otherwise disable.
  if ((sf_last->x_scale_fp != sf_golden->x_scale_fp) ||
      (sf_last->y_scale_fp != sf_golden->y_scale_fp)) {
    gf_temporal_ref = 0;
  }

  estimate_ref_frame_costs(cm, xd, x, segment_id, ref_costs_single,
                           ref_costs_comp);

  memset(&mode_checked[0][0], 0, MB_MODE_COUNT * REF_FRAMES);

  x->skip = 0;

  // Instead of using av1_get_pred_context_switchable_interp(xd) to assign
  // filter_ref, we use a less strict condition on assigning filter_ref.
  // This is to reduce the probabily of entering the flow of not assigning
  // filter_ref and then skip filter search.
  filter_ref = cm->interp_filter;

  // initialize mode decisions
  av1_invalid_rd_stats(&best_rdc);
  av1_invalid_rd_stats(rd_cost);
  mi->sb_type = bsize;
  mi->ref_frame[0] = NONE_FRAME;
  mi->ref_frame[1] = NONE_FRAME;

  mi->tx_size = AOMMIN(
      AOMMIN(max_txsize_lookup[bsize], tx_mode_to_biggest_tx_size[cm->tx_mode]),
      TX_32X32);

  x->source_variance =
      av1_get_sby_perpixel_variance(cpi, &x->plane[0].src, bsize);

  if (cpi->rc.frames_since_golden == 0 && gf_temporal_ref) {
    usable_ref_frame = LAST_FRAME;
  } else {
    usable_ref_frame = GOLDEN_FRAME;
  }

  // For svc mode, on spatial_layer_id > 0: if the reference has different scale
  // constrain the inter mode to only test zero motion.

  force_skip_low_temp_var = 1;
  //        get_force_skip_low_temp_var(&x->variance_low[0], mi_row, mi_col,
  //        bsize);
  // If force_skip_low_temp_var is set, and for short circuit mode = 1 and 3,
  // skip golden reference.
  if (force_skip_low_temp_var) {
    usable_ref_frame = LAST_FRAME;
  }

  if (!((cpi->ref_frame_flags & flag_list[GOLDEN_FRAME]) &&
        !force_skip_low_temp_var))
    use_golden_nonzeromv = 0;

  // If the segment reference frame feature is enabled and it's set to GOLDEN
  // reference, then make sure we don't skip checking GOLDEN, this is to
  // prevent possibility of not picking any mode.
  if (segfeature_active(seg, mi->segment_id, SEG_LVL_REF_FRAME) &&
      get_segdata(seg, mi->segment_id, SEG_LVL_REF_FRAME) == GOLDEN_FRAME) {
    usable_ref_frame = GOLDEN_FRAME;
    skip_ref_find_pred[GOLDEN_FRAME] = 0;
  }

  for (ref_frame = LAST_FRAME; ref_frame <= usable_ref_frame; ++ref_frame) {
    // Skip find_predictor if the reference frame is not in the
    // ref_frame_flags (i.e., not used as a reference for this frame).
    skip_ref_find_pred[ref_frame] =
        !(cpi->ref_frame_flags & flag_list[ref_frame]);
    if (!skip_ref_find_pred[ref_frame]) {
      find_predictors(cpi, x, ref_frame, frame_mv, const_motion,
                      &ref_frame_skip_mask, flag_list, tile_data, mi_row,
                      mi_col, yv12_mb, bsize, force_skip_low_temp_var,
                      comp_modes > 0);
    }
  }

  for (int i = 0; i < bsize_to_num_blk(bsize); ++i) {
    const int planes = av1_num_planes(cm);
    for (int j = 0; j < planes; j++) set_blk_skip(x, j, i, 0);
  }
  for (idx = 0; idx < num_inter_modes; ++idx) {
    int rate_mv = 0;
    int mode_rd_thresh;
    int mode_index;
    int i;
    int64_t this_sse;
    int is_skippable;
    int this_early_term = 0;
    int rd_computed = 0;
    int inter_mv_mode = 0;
    int skip_this_mv = 0;
    int comp_pred = 0;
    int force_mv_inter_layer = 0;
    PREDICTION_MODE this_mode;
    second_ref_frame = NONE_FRAME;

    this_mode = ref_mode_set[idx].pred_mode;
    ref_frame = ref_mode_set[idx].ref_frame;
    init_mbmi(mi, idx, cm);

    if (ref_frame > usable_ref_frame) continue;
    if (skip_ref_find_pred[ref_frame]) continue;

    // If the segment reference frame feature is enabled then do nothing if the
    // current ref frame is not allowed.
    if (segfeature_active(seg, mi->segment_id, SEG_LVL_REF_FRAME) &&
        get_segdata(seg, mi->segment_id, SEG_LVL_REF_FRAME) != (int)ref_frame)
      continue;

    // For CBR mode: skip the golden reference search if sse of zeromv_last is
    // below threshold.
    if (ref_frame == GOLDEN_FRAME && cpi->oxcf.rc_mode == AOM_CBR &&
        sse_zeromv_normalized < thresh_skip_golden)
      continue;

    if (!(cpi->ref_frame_flags & flag_list[ref_frame])) continue;

    //    if (x->source_variance == 0 && frame_mv[this_mode][ref_frame].as_int
    //    != 0) {
    //      continue;
    //    }

    if (!(inter_mode_mask[bsize] & (1 << this_mode))) continue;

    if (const_motion[ref_frame] && this_mode == NEARMV) continue;

    // Skip non-zeromv mode search for golden frame if force_skip_low_temp_var
    // is set. If nearestmv for golden frame is 0, zeromv mode will be skipped
    // later.
    if (!force_mv_inter_layer && force_skip_low_temp_var &&
        ref_frame == GOLDEN_FRAME &&
        frame_mv[this_mode][ref_frame].as_int != 0) {
      continue;
    }

    /*
        if (x->content_state_sb != kVeryHighSad &&
            (cpi->sf.short_circuit_low_temp_var >= 2 ||
             (cpi->sf.short_circuit_low_temp_var == 1 && bsize == BLOCK_64X64))
       && force_skip_low_temp_var && ref_frame == LAST_FRAME && this_mode ==
       NEWMV) { continue;
        }

        // Disable this drop out case if the ref frame segment level feature is
        // enabled for this segment. This is to prevent the possibility that we
       end
        // up unable to pick any mode.
        if (!segfeature_active(seg, mi->segment_id, SEG_LVL_REF_FRAME)) {
          if (sf->reference_masking &&
              !(frame_mv[this_mode][ref_frame].as_int == 0 &&
                ref_frame == LAST_FRAME)) {
            if (usable_ref_frame < ALTREF_FRAME) {
              if (!force_skip_low_temp_var && usable_ref_frame > LAST_FRAME) {
                i = (ref_frame == LAST_FRAME) ? GOLDEN_FRAME : LAST_FRAME;
                if ((cpi->ref_frame_flags & flag_list[i]))
                  if (x->pred_mv_sad[ref_frame] > (x->pred_mv_sad[i] << 1))
                    ref_frame_skip_mask |= (1 << ref_frame);
              }
            } else if (!cpi->rc.is_src_frame_alt_ref &&
                       !(frame_mv[this_mode][ref_frame].as_int == 0 &&
                         ref_frame == ALTREF_FRAME)) {
              int ref1 = (ref_frame == GOLDEN_FRAME) ? LAST_FRAME :
       GOLDEN_FRAME; int ref2 = (ref_frame == ALTREF_FRAME) ? LAST_FRAME :
       ALTREF_FRAME; if (((cpi->ref_frame_flags & flag_list[ref1]) &&
                   (x->pred_mv_sad[ref_frame] > (x->pred_mv_sad[ref1] << 1))) ||
                  ((cpi->ref_frame_flags & flag_list[ref2]) &&
                   (x->pred_mv_sad[ref_frame] > (x->pred_mv_sad[ref2] << 1))))
                ref_frame_skip_mask |= (1 << ref_frame);
            }
          }
          if (ref_frame_skip_mask & (1 << ref_frame)) continue;
        }
    */

    // Select prediction reference frames.
    for (i = 0; i < MAX_MB_PLANE; i++) {
      xd->plane[i].pre[0] = yv12_mb[ref_frame][i];
      if (comp_pred) xd->plane[i].pre[1] = yv12_mb[second_ref_frame][i];
    }

    mi->ref_frame[0] = ref_frame;
    mi->ref_frame[1] = second_ref_frame;
    set_ref_ptrs(cm, xd, ref_frame, second_ref_frame);

    mode_index = mode_idx[ref_frame][INTER_OFFSET(this_mode)];
    mode_rd_thresh = best_pickmode.best_mode_skip_txfm
                         ? rd_threshes[mode_index] << 1
                         : rd_threshes[mode_index];

    if (rd_less_than_thresh(best_rdc.rdcost, mode_rd_thresh,
                            rd_thresh_freq_fact[mode_index]))
      if (frame_mv[this_mode][ref_frame].as_int != 0) continue;

    if (this_mode == NEWMV && !force_mv_inter_layer) {
      if (search_new_mv(cpi, x, frame_mv, ref_frame, gf_temporal_ref, bsize,
                        mi_row, mi_col, best_pred_sad, &rate_mv, best_sse_sofar,
                        &best_rdc))
        continue;
    }

    // TODO(jianj): Skipping the testing of (duplicate) non-zero motion vector
    // causes some regression, leave it for duplicate zero-mv for now, until
    // regression issue is resolved.
    for (inter_mv_mode = NEARESTMV; inter_mv_mode <= NEWMV; inter_mv_mode++) {
      if (inter_mv_mode == this_mode || comp_pred) continue;
      if (mode_checked[inter_mv_mode][ref_frame] &&
          frame_mv[this_mode][ref_frame].as_int ==
              frame_mv[inter_mv_mode][ref_frame].as_int &&
          frame_mv[inter_mv_mode][ref_frame].as_int == 0) {
        skip_this_mv = 1;
        break;
      }
    }

    if (skip_this_mv) continue;

    // If use_golden_nonzeromv is false, NEWMV mode is skipped for golden, no
    // need to compute best_pred_sad which is only used to skip golden NEWMV.
    if (use_golden_nonzeromv && this_mode == NEWMV && ref_frame == LAST_FRAME &&
        frame_mv[NEWMV][LAST_FRAME].as_int != INVALID_MV) {
      const int pre_stride = xd->plane[0].pre[0].stride;
      const uint8_t *const pre_buf =
          xd->plane[0].pre[0].buf +
          (frame_mv[NEWMV][LAST_FRAME].as_mv.row >> 3) * pre_stride +
          (frame_mv[NEWMV][LAST_FRAME].as_mv.col >> 3);
      best_pred_sad = cpi->fn_ptr[bsize].sdf(
          x->plane[0].src.buf, x->plane[0].src.stride, pre_buf, pre_stride);
      x->pred_mv_sad[LAST_FRAME] = best_pred_sad;
    }

    if (this_mode != NEARESTMV && !comp_pred &&
        frame_mv[this_mode][ref_frame].as_int ==
            frame_mv[NEARESTMV][ref_frame].as_int)
      continue;

    mi->mode = this_mode;
    mi->mv[0].as_int = frame_mv[this_mode][ref_frame].as_int;
    mi->mv[1].as_int = 0;

    /*
        if ((this_mode == NEWMV || filter_ref == SWITCHABLE) &&
            pred_filter_search &&
            (ref_frame == LAST_FRAME ||
             (ref_frame == GOLDEN_FRAME && !force_mv_inter_layer &&
              (cpi->use_svc || cpi->oxcf.rc_mode == VPX_VBR))) &&
            (((mi->mv[0].as_mv.row | mi->mv[0].as_mv.col) & 0x07) != 0)) {
          rd_computed = 1;
          search_filter_ref(cpi, x, &this_rdc, mi_row, mi_col, tmp, bsize,
                            reuse_inter_pred, &this_mode_pred, &var_y, &sse_y);
        } else
    */
    {
      mi->interp_filters =
          (filter_ref == SWITCHABLE) ? EIGHTTAP_REGULAR : filter_ref;
      av1_enc_build_inter_predictor(cm, xd, mi_row, mi_col, NULL, bsize,
                                    AOM_PLANE_Y, AOM_PLANE_Y);

      // For large partition blocks, extra testing is done.
      rd_computed = 1;
      model_rd_for_sb(cpi, bsize, x, xd, AOM_PLANE_Y, AOM_PLANE_Y, mi_row,
                      mi_col, &this_rdc.rate, &this_rdc.dist, &this_rdc.skip,
                      NULL, NULL, NULL, NULL);
      // Save normalized sse (between current and last frame) for (0, 0) motion.
      if (ref_frame == LAST_FRAME &&
          frame_mv[this_mode][ref_frame].as_int == 0) {
        sse_zeromv_normalized =
            sse_y >> (block_size_wide[bsize] + block_size_high[bsize]);
      }
      if (sse_y < best_sse_sofar) best_sse_sofar = sse_y;
    }

    const int skip_mode_ctx = av1_get_skip_mode_context(xd);
    const int skip_mode_cost = x->skip_mode_cost[skip_mode_ctx][1];
    if (!this_early_term) {
      this_sse = (int64_t)sse_y;

      block_yrd(cpi, x, mi_row, mi_col, &this_rdc, &is_skippable, &this_sse,
                bsize, mi->tx_size, rd_computed);

      x->skip = is_skippable;
      if (is_skippable) {
        set_skip_flag(x, &this_rdc, bsize, this_sse);
      } else {
        if (RDCOST(x->rdmult, this_rdc.rate, this_rdc.dist) >=
            RDCOST(x->rdmult, skip_mode_cost, this_sse)) {
          set_skip_flag(x, &this_rdc, bsize, this_sse);
        }
      }

    } else {
      this_rdc.rate += skip_mode_cost;
    }
    /*
        if (!this_early_term &&
            (x->color_sensitivity[0] || x->color_sensitivity[1])) {
          RD_COST rdc_uv;
          const BLOCK_SIZE uv_bsize = get_plane_block_size(bsize,
       &xd->plane[1]); if (x->color_sensitivity[0] && !flag_preduv_computed[0])
       { av1_build_inter_predictors_sbp(xd, mi_row, mi_col, bsize, 1);
            flag_preduv_computed[0] = 1;
          }
          if (x->color_sensitivity[1] && !flag_preduv_computed[1]) {
            av1_build_inter_predictors_sbp(xd, mi_row, mi_col, bsize, 2);
            flag_preduv_computed[1] = 1;
          }
          model_rd_for_sb_uv(cpi, uv_bsize, x, xd, &rdc_uv, &var_y, &sse_y, 1,
       2); this_rdc.rate += rdc_uv.rate; this_rdc.dist += rdc_uv.dist;
        }
    */
    this_rdc.rate += rate_mv;
    //    this_rdc.rate +=
    //    cpi->inter_mode_cost[x->mbmi_ext->mode_context[ref_frame]]
    //                                         [INTER_OFFSET(this_mode)];
    // TODO(marpan): Add costing for compound mode.
    this_rdc.rate += ref_costs_single[ref_frame];
    this_rdc.rdcost = RDCOST(x->rdmult, this_rdc.rate, this_rdc.dist);
    /*
        // Bias against NEWMV that is very different from its neighbors, and
       bias
        // to small motion-lastref for noisy input.
        if (cpi->oxcf.rc_mode == VPX_CBR && cpi->oxcf.speed >= 5 &&
            cpi->oxcf.content != VP9E_CONTENT_SCREEN) {
          av1_NEWMV_diff_bias(&cpi->noise_estimate, xd, this_mode, &this_rdc,
       bsize, frame_mv[this_mode][ref_frame].as_mv.row,
                              frame_mv[this_mode][ref_frame].as_mv.col,
                              ref_frame == LAST_FRAME, x->lowvar_highsumdiff,
                              x->sb_is_skin);
        }

    */

    mode_checked[this_mode][ref_frame] = 1;

    if (this_rdc.rdcost < best_rdc.rdcost || x->skip) {
      best_rdc = this_rdc;
      best_early_term = this_early_term;
      best_pickmode.best_mode = this_mode;
      best_pickmode.best_pred_filter = mi->interp_filters;
      best_pickmode.best_tx_size = mi->tx_size;
      best_pickmode.best_ref_frame = ref_frame;
      best_pickmode.best_mode_skip_txfm = x->skip;
      best_pickmode.best_second_ref_frame = second_ref_frame;
    }
    if (x->skip && !force_test_gf_zeromv) break;

    // If early termination flag is 1 and at least 2 modes are checked,
    // the mode search is terminated.
    if (best_early_term && idx > 0) {
      x->skip = 1;
      break;
    }
  }

  mi->mode = best_pickmode.best_mode;
  mi->interp_filters = best_pickmode.best_pred_filter;
  mi->tx_size = best_pickmode.best_tx_size;
  mi->ref_frame[0] = best_pickmode.best_ref_frame;
  mi->mv[0].as_int =
      frame_mv[best_pickmode.best_mode][best_pickmode.best_ref_frame].as_int;
  x->skip = best_pickmode.best_mode_skip_txfm;
  mi->ref_frame[1] = best_pickmode.best_second_ref_frame;

  // If the segment reference frame feature is enabled and set then
  // skip the intra prediction.
  /*
  if (segfeature_active(seg, mi->segment_id, SEG_LVL_REF_FRAME) &&
      get_segdata(seg, mi->segment_id, SEG_LVL_REF_FRAME) > 0)
    perform_intra_pred = 0;
    // Perform intra prediction search, if the best SAD is above a certain
    // threshold.
    if (best_rdc.rdcost == INT64_MAX ||
        (cpi->oxcf.content == VP9E_CONTENT_SCREEN && x->source_variance == 0 &&
         !x->zero_temp_sad_source) ||
        (scene_change_detected && perform_intra_pred) ||
        ((!force_skip_low_temp_var || bsize < BLOCK_32X32 ||
          x->content_state_sb == kVeryHighSad) &&
         perform_intra_pred && !x->skip && best_rdc.rdcost > inter_mode_thresh
  && bsize <= cpi->sf.max_intra_bsize && !x->skip_low_source_sad &&
         !x->lowvar_highsumdiff)) {
      struct estimate_block_intra_args args = { cpi, x, DC_PRED, 1, 0 };
      int i;
      PRED_BUFFER *const best_pred = best_pickmode.best_pred;
      TX_SIZE intra_tx_size =
          VPXMIN(max_txsize_lookup[bsize],
                 tx_mode_to_biggest_tx_size[cpi->common.tx_mode]);
      if (cpi->oxcf.content != VP9E_CONTENT_SCREEN && intra_tx_size > TX_16X16)
        intra_tx_size = TX_16X16;

      if (reuse_inter_pred && best_pred != NULL) {
        if (best_pred->data == orig_dst.buf) {
          this_mode_pred = &tmp[get_pred_buffer(tmp, 3)];
  #if CONFIG_VP9_HIGHBITDEPTH
          if (cm->use_highbitdepth)
            vpx_highbd_convolve_copy(
                CONVERT_TO_SHORTPTR(best_pred->data), best_pred->stride,
                CONVERT_TO_SHORTPTR(this_mode_pred->data),
  this_mode_pred->stride, NULL, 0, 0, 0, 0, bw, bh, xd->bd); else
            vpx_convolve_copy(best_pred->data, best_pred->stride,
                              this_mode_pred->data, this_mode_pred->stride,
  NULL, 0, 0, 0, 0, bw, bh); #else vpx_convolve_copy(best_pred->data,
  best_pred->stride, this_mode_pred->data, this_mode_pred->stride, NULL, 0, 0,
  0, 0, bw, bh); #endif  // CONFIG_VP9_HIGHBITDEPTH best_pickmode.best_pred =
  this_mode_pred;
        }
      }
      pd->dst = orig_dst;

      for (i = 0; i < 4; ++i) {
        const PREDICTION_MODE this_mode = intra_mode_list[i];
        THR_MODES mode_index = mode_idx[INTRA_FRAME][mode_offset(this_mode)];
        int mode_rd_thresh = rd_threshes[mode_index];
        // For spatially flat blocks, under short_circuit_flat_blocks flag:
        // only check DC mode for stationary blocks, otherwise also check
        // H and V mode.
        if (sf->short_circuit_flat_blocks && x->source_variance == 0 &&
            ((x->zero_temp_sad_source && this_mode != DC_PRED) || i > 2)) {
          continue;
        }

        if (!((1 << this_mode) & cpi->sf.intra_y_mode_bsize_mask[bsize]))
          continue;

        if ((cpi->sf.adaptive_rd_thresh_row_mt &&
             rd_less_than_thresh_row_mt(best_rdc.rdcost, mode_rd_thresh,
                                        &rd_thresh_freq_fact[mode_index])) ||
            (!cpi->sf.adaptive_rd_thresh_row_mt &&
             rd_less_than_thresh(best_rdc.rdcost, mode_rd_thresh,
                                 &rd_thresh_freq_fact[mode_index])))
          continue;

        mi->mode = this_mode;
        mi->ref_frame[0] = INTRA_FRAME;
        this_rdc.dist = this_rdc.rate = 0;
        args.mode = this_mode;
        args.skippable = 1;
        args.rdc = &this_rdc;
        mi->tx_size = intra_tx_size;
        av1_foreach_transformed_block_in_plane(xd, bsize, 0,
  estimate_block_intra, &args);
        // Check skip cost here since skippable is not set for for uv, this
        // mirrors the behavior used by inter
        if (args.skippable) {
          x->skip_txfm[0] = SKIP_TXFM_AC_DC;
          this_rdc.rate = av1_cost_bit(av1_get_skip_prob(&cpi->common, xd), 1);
        } else {
          x->skip_txfm[0] = SKIP_TXFM_NONE;
          this_rdc.rate += av1_cost_bit(av1_get_skip_prob(&cpi->common, xd), 0);
        }
        // Inter and intra RD will mismatch in scale for non-screen content.
        if (cpi->oxcf.content == VP9E_CONTENT_SCREEN) {
          if (x->color_sensitivity[0])
            av1_foreach_transformed_block_in_plane(xd, bsize, 1,
                                                   estimate_block_intra, &args);
          if (x->color_sensitivity[1])
            av1_foreach_transformed_block_in_plane(xd, bsize, 2,
                                                   estimate_block_intra, &args);
        }
        this_rdc.rate += cpi->mbmode_cost[this_mode];
        this_rdc.rate += ref_frame_cost[INTRA_FRAME];
        this_rdc.rate += intra_cost_penalty;
        this_rdc.rdcost =
            RDCOST(x->rdmult, x->rddiv, this_rdc.rate, this_rdc.dist);

        if (this_rdc.rdcost < best_rdc.rdcost) {
          best_rdc = this_rdc;
          best_pickmode.best_mode = this_mode;
          best_pickmode.best_intra_tx_size = mi->tx_size;
          best_pickmode.best_ref_frame = INTRA_FRAME;
          best_pickmode.best_second_ref_frame = NONE;
          mi->uv_mode = this_mode;
          mi->mv[0].as_int = INVALID_MV;
          mi->mv[1].as_int = INVALID_MV;
          best_pickmode.best_mode_skip_txfm = x->skip_txfm[0];
        }
      }

      // Reset mb_mode_info to the best inter mode.
      if (best_pickmode.best_ref_frame != INTRA_FRAME) {
        mi->tx_size = best_pickmode.best_tx_size;
      } else {
        mi->tx_size = best_pickmode.best_intra_tx_size;
      }
    }
    pd->dst = orig_dst;
    mi->mode = best_pickmode.best_mode;
    mi->ref_frame[0] = best_pickmode.best_ref_frame;
    mi->ref_frame[1] = best_pickmode.best_second_ref_frame;
    x->skip[0] = best_pickmode.best_mode_skip_txfm;

    if (!is_inter_block(mi)) {
      mi->interp_filters = SWITCHABLE_FILTERS;
    }


    if (best_pickmode.best_ref_frame == ALTREF_FRAME ||
        best_pickmode.best_second_ref_frame == ALTREF_FRAME)
      x->arf_frame_usage++;
    else if (best_pickmode.best_ref_frame != INTRA_FRAME)
      x->lastgolden_frame_usage++;

    if (cpi->sf.adaptive_rd_thresh) {
      THR_MODES best_mode_idx =
          mode_idx[best_pickmode.best_ref_frame][mode_offset(mi->mode)];

      if (best_pickmode.best_ref_frame == INTRA_FRAME) {
        // Only consider the modes that are included in the intra_mode_list.
        int intra_modes = sizeof(intra_mode_list) / sizeof(PREDICTION_MODE);
        int i;

        // TODO(yunqingwang): Check intra mode mask and only update freq_fact
        // for those valid modes.
        for (i = 0; i < intra_modes; i++) {
          update_thresh_freq_fact(cpi, tile_data, x->source_variance, bsize,
                                    INTRA_FRAME, best_mode_idx,
                                    intra_mode_list[i]);
        }
      } else {
        for (ref_frame = LAST_FRAME; ref_frame <= GOLDEN_FRAME; ++ref_frame) {
          PREDICTION_MODE this_mode;
          if (best_pickmode.best_ref_frame != ref_frame) continue;
          for (this_mode = NEARESTMV; this_mode <= NEWMV; ++this_mode) {
            update_thresh_freq_fact(cpi, tile_data, x->source_variance, bsize,
                                      ref_frame, best_mode_idx, this_mode);
          }
        }
      }
    }
    */

  store_coding_context(x, ctx, mi->mode, x->skip);
  *rd_cost = best_rdc;
}
