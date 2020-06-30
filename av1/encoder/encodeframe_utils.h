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

#ifndef AOM_AV1_ENCODER_ENCODEFRAME_UTILS_H_
#define AOM_AV1_ENCODER_ENCODEFRAME_UTILS_H_

#include "av1/common/reconinter.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/partition_strategy.h"
#include "av1/encoder/rdopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!\cond */

enum { PICK_MODE_RD = 0, PICK_MODE_NONRD };

enum {
  SB_SINGLE_PASS,  // Single pass encoding: all ctxs get updated normally
  SB_DRY_PASS,     // First pass of multi-pass: does not update the ctxs
  SB_WET_PASS      // Second pass of multi-pass: finalize and update the ctx
} UENUM1BYTE(SB_MULTI_PASS_MODE);

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

// This struct is used to store the statistics used by sb-level multi-pass
// encoding. Currently, this is only used to make a copy of the state before we
// perform the first pass
typedef struct SB_FIRST_PASS_STATS {
  RD_SEARCH_MACROBLOCK_CONTEXT x_ctx;
  RD_COUNTS rd_count;

  int split_count;
  FRAME_COUNTS fc;
  InterModeRdModel inter_mode_rd_models[BLOCK_SIZES_ALL];
  int thresh_freq_fact[BLOCK_SIZES_ALL][MAX_MODES];
  int current_qindex;

#if CONFIG_INTERNAL_STATS
  unsigned int mode_chosen_counts[MAX_MODES];
#endif  // CONFIG_INTERNAL_STATS
} SB_FIRST_PASS_STATS;
/*!\endcond */

/*!\cond */
// This structure contains block size related
// variables for use in rd_pick_partition().
typedef struct {
  // Half of block width to determine block edge.
  int mi_step;

  // Block row and column indices.
  int mi_row;
  int mi_col;

  // Block edge row and column indices.
  int mi_row_edge;
  int mi_col_edge;

  // Block width of current partition block.
  int width;

  // Block width of minimum partition size allowed.
  int min_partition_size_1d;

  // Flag to indicate if partition is 8x8 or higher size.
  int bsize_at_least_8x8;

  // Indicates edge blocks in frame.
  int has_rows;
  int has_cols;

  // Block size of current partition.
  BLOCK_SIZE bsize;

  // Size of current sub-partition.
  BLOCK_SIZE subsize;

  // Size of split partition.
  BLOCK_SIZE split_bsize2;
} PartitionBlkParams;

// Structure holding state variables for partition search.
typedef struct {
  // Intra partitioning related info.
  PartitionSearchInfo *intra_part_info;

  // Parameters related to partition block size.
  PartitionBlkParams part_blk_params;

  // Win flags for HORZ and VERT partition evaluations.
  RD_RECT_PART_WIN_INFO split_part_rect_win[4];

  // RD cost for the current block of given partition type.
  RD_STATS this_rdc;

  // RD cost summed across all blocks of partition type.
  RD_STATS sum_rdc;

  // Array holding partition type cost.
  int tmp_partition_cost[PARTITION_TYPES];

  // Pointer to partition cost buffer
  int *partition_cost;

  // RD costs for different partition types.
  int64_t none_rd;
  int64_t split_rd[4];
  // RD costs for rectangular partitions.
  // rect_part_rd[0][i] is the RD cost of ith partition index of PARTITION_HORZ.
  // rect_part_rd[1][i] is the RD cost of ith partition index of PARTITION_VERT.
  int64_t rect_part_rd[NUM_RECT_PARTS][2];

  // Flags indicating if the corresponding partition was winner or not.
  // Used to bypass similar blocks during AB partition evaluation.
  int is_split_ctx_is_ready[2];
  int is_rect_ctx_is_ready[NUM_RECT_PARTS];

  // Flags to prune/skip particular partition size evaluation.
  int terminate_partition_search;
  int partition_none_allowed;
  int partition_rect_allowed[NUM_RECT_PARTS];
  int do_rectangular_split;
  int do_square_split;
  int prune_rect_part[NUM_RECT_PARTS];

  // Chroma subsampling in x and y directions.
  int ss_x;
  int ss_y;

  // Partition plane context index.
  int pl_ctx_idx;

  // This flag will be set if best partition is found from the search.
  bool found_best_partition;
} PartitionSearchState;
/*!\endcond */

static AOM_INLINE void update_global_motion_used(PREDICTION_MODE mode,
                                                 BLOCK_SIZE bsize,
                                                 const MB_MODE_INFO *mbmi,
                                                 RD_COUNTS *rdc) {
  if (mode == GLOBALMV || mode == GLOBAL_GLOBALMV) {
    const int num_4x4s = mi_size_wide[bsize] * mi_size_high[bsize];
    int ref;
    for (ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
      rdc->global_motion_used[mbmi->ref_frame[ref]] += num_4x4s;
    }
  }
}

static AOM_INLINE void update_filter_type_cdf(const MACROBLOCKD *xd,
                                              const MB_MODE_INFO *mbmi) {
  int dir;
  for (dir = 0; dir < 2; ++dir) {
    const int ctx = av1_get_pred_context_switchable_interp(xd, dir);
    InterpFilter filter = av1_extract_interp_filter(mbmi->interp_filters, dir);
    update_cdf(xd->tile_ctx->switchable_interp_cdf[ctx], filter,
               SWITCHABLE_FILTERS);
  }
}

static AOM_INLINE int set_segment_rdmult(const AV1_COMP *const cpi,
                                         MACROBLOCK *const x,
                                         int8_t segment_id) {
  const AV1_COMMON *const cm = &cpi->common;
  av1_init_plane_quantizers(cpi, x, segment_id);
  aom_clear_system_state();
  const int segment_qindex =
      av1_get_qindex(&cm->seg, segment_id, cm->quant_params.base_qindex);
  return av1_compute_rd_mult(cpi,
                             segment_qindex + cm->quant_params.y_dc_delta_q);
}

static AOM_INLINE int do_slipt_check(BLOCK_SIZE bsize) {
  return (bsize == BLOCK_16X16 || bsize == BLOCK_32X32);
}

#if !CONFIG_REALTIME_ONLY
static AOM_INLINE const FIRSTPASS_STATS *read_one_frame_stats(const TWO_PASS *p,
                                                              int frm) {
  assert(frm >= 0);
  if (frm < 0 ||
      p->stats_buf_ctx->stats_in_start + frm > p->stats_buf_ctx->stats_in_end) {
    return NULL;
  }

  return &p->stats_buf_ctx->stats_in_start[frm];
}

static int get_rdmult_delta(AV1_COMP *cpi, BLOCK_SIZE bsize, int analysis_type,
                            int mi_row, int mi_col, int orig_rdmult) {
  AV1_COMMON *const cm = &cpi->common;
  const GF_GROUP *const gf_group = &cpi->gf_group;
  assert(IMPLIES(cpi->gf_group.size > 0,
                 cpi->gf_group.index < cpi->gf_group.size));
  const int tpl_idx = cpi->gf_group.index;
  TplParams *const tpl_data = &cpi->tpl_data;
  TplDepFrame *tpl_frame = &tpl_data->tpl_frame[tpl_idx];
  TplDepStats *tpl_stats = tpl_frame->tpl_stats_ptr;
  const uint8_t block_mis_log2 = tpl_data->tpl_stats_block_mis_log2;
  int tpl_stride = tpl_frame->stride;
  int64_t intra_cost = 0;
  int64_t mc_dep_cost = 0;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  if (tpl_frame->is_valid == 0) return orig_rdmult;

  if (!is_frame_tpl_eligible(gf_group)) return orig_rdmult;

  if (cpi->gf_group.index >= MAX_TPL_FRAME_IDX) return orig_rdmult;

  int64_t mc_count = 0, mc_saved = 0;
  int mi_count = 0;
  const int mi_col_sr =
      coded_to_superres_mi(mi_col, cm->superres_scale_denominator);
  const int mi_col_end_sr =
      coded_to_superres_mi(mi_col + mi_wide, cm->superres_scale_denominator);
  const int mi_cols_sr = av1_pixels_to_mi(cm->superres_upscaled_width);
  const int step = 1 << block_mis_log2;
  for (int row = mi_row; row < mi_row + mi_high; row += step) {
    for (int col = mi_col_sr; col < mi_col_end_sr; col += step) {
      if (row >= cm->mi_params.mi_rows || col >= mi_cols_sr) continue;
      TplDepStats *this_stats =
          &tpl_stats[av1_tpl_ptr_pos(row, col, tpl_stride, block_mis_log2)];
      int64_t mc_dep_delta =
          RDCOST(tpl_frame->base_rdmult, this_stats->mc_dep_rate,
                 this_stats->mc_dep_dist);
      intra_cost += this_stats->recrf_dist << RDDIV_BITS;
      mc_dep_cost += (this_stats->recrf_dist << RDDIV_BITS) + mc_dep_delta;
      mc_count += this_stats->mc_count;
      mc_saved += this_stats->mc_saved;
      mi_count++;
    }
  }

  aom_clear_system_state();

  double beta = 1.0;
  if (analysis_type == 0) {
    if (mc_dep_cost > 0 && intra_cost > 0) {
      const double r0 = cpi->rd.r0;
      const double rk = (double)intra_cost / mc_dep_cost;
      beta = (r0 / rk);
    }
  } else if (analysis_type == 1) {
    const double mc_count_base = (mi_count * cpi->rd.mc_count_base);
    beta = (mc_count + 1.0) / (mc_count_base + 1.0);
    beta = pow(beta, 0.5);
  } else if (analysis_type == 2) {
    const double mc_saved_base = (mi_count * cpi->rd.mc_saved_base);
    beta = (mc_saved + 1.0) / (mc_saved_base + 1.0);
    beta = pow(beta, 0.5);
  }

  int rdmult = av1_get_adaptive_rdmult(cpi, beta);

  aom_clear_system_state();

  rdmult = AOMMIN(rdmult, orig_rdmult * 3 / 2);
  rdmult = AOMMAX(rdmult, orig_rdmult * 1 / 2);

  rdmult = AOMMAX(1, rdmult);

  return rdmult;
}

static AOM_INLINE void adjust_rdmult_tpl_model(AV1_COMP *cpi, MACROBLOCK *x,
                                               int mi_row, int mi_col) {
  const BLOCK_SIZE sb_size = cpi->common.seq_params.sb_size;
  const int orig_rdmult = cpi->rd.RDMULT;

  assert(IMPLIES(cpi->gf_group.size > 0,
                 cpi->gf_group.index < cpi->gf_group.size));
  const int gf_group_index = cpi->gf_group.index;
  if (cpi->oxcf.enable_tpl_model && cpi->oxcf.q_cfg.aq_mode == NO_AQ &&
      cpi->oxcf.q_cfg.deltaq_mode == NO_DELTA_Q && gf_group_index > 0 &&
      cpi->gf_group.update_type[gf_group_index] == ARF_UPDATE) {
    const int dr =
        get_rdmult_delta(cpi, sb_size, 0, mi_row, mi_col, orig_rdmult);
    x->rdmult = dr;
  }
}

static BLOCK_SIZE dim_to_size(int dim) {
  switch (dim) {
    case 4: return BLOCK_4X4;
    case 8: return BLOCK_8X8;
    case 16: return BLOCK_16X16;
    case 32: return BLOCK_32X32;
    case 64: return BLOCK_64X64;
    case 128: return BLOCK_128X128;
    default: assert(0); return 0;
  }
}

static AOM_INLINE void set_max_min_partition_size(SuperBlockEnc *sb_enc,
                                                  AV1_COMP *cpi, MACROBLOCK *x,
                                                  const SPEED_FEATURES *sf,
                                                  BLOCK_SIZE sb_size,
                                                  int mi_row, int mi_col) {
  const AV1_COMMON *cm = &cpi->common;

  sb_enc->max_partition_size =
      AOMMIN(sf->part_sf.default_max_partition_size,
             dim_to_size(cpi->oxcf.part_cfg.max_partition_size));
  sb_enc->min_partition_size =
      AOMMAX(sf->part_sf.default_min_partition_size,
             dim_to_size(cpi->oxcf.part_cfg.min_partition_size));
  sb_enc->max_partition_size =
      AOMMIN(sb_enc->max_partition_size, cm->seq_params.sb_size);
  sb_enc->min_partition_size =
      AOMMIN(sb_enc->min_partition_size, cm->seq_params.sb_size);

  if (use_auto_max_partition(cpi, sb_size, mi_row, mi_col)) {
    float features[FEATURE_SIZE_MAX_MIN_PART_PRED] = { 0.0f };

    av1_get_max_min_partition_features(cpi, x, mi_row, mi_col, features);
    sb_enc->max_partition_size =
        AOMMAX(AOMMIN(av1_predict_max_partition(cpi, x, features),
                      sb_enc->max_partition_size),
               sb_enc->min_partition_size);
  }
}

int av1_active_h_edge(const AV1_COMP *cpi, int mi_row, int mi_step);

int av1_active_v_edge(const AV1_COMP *cpi, int mi_col, int mi_step);

void av1_get_tpl_stats_sb(AV1_COMP *cpi, BLOCK_SIZE bsize, int mi_row,
                          int mi_col, SuperBlockEnc *sb_enc);

int av1_get_q_for_deltaq_objective(AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                   int mi_row, int mi_col);
#endif  // !CONFIG_REALTIME_ONLY

void av1_set_ssim_rdmult(const AV1_COMP *const cpi, MvCosts *const mv_costs,
                         const BLOCK_SIZE bsize, const int mi_row,
                         const int mi_col, int *const rdmult);

int av1_get_hier_tpl_rdmult(const AV1_COMP *const cpi, MACROBLOCK *const x,
                            const BLOCK_SIZE bsize, const int mi_row,
                            const int mi_col, int orig_rdmult);

void av1_update_state(const AV1_COMP *const cpi, ThreadData *td,
                      const PICK_MODE_CONTEXT *const ctx, int mi_row,
                      int mi_col, BLOCK_SIZE bsize, RUN_TYPE dry_run);

void av1_update_inter_mode_stats(FRAME_CONTEXT *fc, FRAME_COUNTS *counts,
                                 PREDICTION_MODE mode, int16_t mode_context);

void av1_sum_intra_stats(const AV1_COMMON *const cm, FRAME_COUNTS *counts,
                         MACROBLOCKD *xd, const MB_MODE_INFO *const mbmi,
                         const MB_MODE_INFO *above_mi,
                         const MB_MODE_INFO *left_mi, const int intraonly);

void av1_restore_context(MACROBLOCK *x, const RD_SEARCH_MACROBLOCK_CONTEXT *ctx,
                         int mi_row, int mi_col, BLOCK_SIZE bsize,
                         const int num_planes);

void av1_save_context(const MACROBLOCK *x, RD_SEARCH_MACROBLOCK_CONTEXT *ctx,
                      int mi_row, int mi_col, BLOCK_SIZE bsize,
                      const int num_planes);

void av1_set_fixed_partitioning(AV1_COMP *cpi, const TileInfo *const tile,
                                MB_MODE_INFO **mib, int mi_row, int mi_col,
                                BLOCK_SIZE bsize);

int av1_is_leaf_split_partition(AV1_COMMON *cm, int mi_row, int mi_col,
                                BLOCK_SIZE bsize);

void av1_reset_simple_motion_tree_partition(SIMPLE_MOTION_DATA_TREE *sms_tree,
                                            BLOCK_SIZE bsize);

void av1_update_picked_ref_frames_mask(MACROBLOCK *const x, int ref_type,
                                       BLOCK_SIZE bsize, int mib_size,
                                       int mi_row, int mi_col);

void av1_avg_cdf_symbols(FRAME_CONTEXT *ctx_left, FRAME_CONTEXT *ctx_tr,
                         int wt_left, int wt_tr);

void av1_source_content_sb(AV1_COMP *cpi, MACROBLOCK *x, int offset);

void av1_reset_mbmi(CommonModeInfoParams *const mi_params, BLOCK_SIZE sb_size,
                    int mi_row, int mi_col);

void av1_backup_sb_state(SB_FIRST_PASS_STATS *sb_fp_stats, const AV1_COMP *cpi,
                         ThreadData *td, const TileDataEnc *tile_data,
                         int mi_row, int mi_col);

void av1_restore_sb_state(const SB_FIRST_PASS_STATS *sb_fp_stats, AV1_COMP *cpi,
                          ThreadData *td, TileDataEnc *tile_data, int mi_row,
                          int mi_col);

void av1_set_cost_upd_freq(AV1_COMP *cpi, ThreadData *td,
                           const TileInfo *const tile_info, const int mi_row,
                           const int mi_col);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_ENCODEFRAME_UTILS_H_
