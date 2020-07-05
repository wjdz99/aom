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
#include "aom_ports/system_state.h"

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
#include "av1/encoder/ml.h"
#include "av1/encoder/motion_search_facade.h"
#include "av1/encoder/partition_strategy.h"
#if !CONFIG_REALTIME_ONLY
#include "av1/encoder/partition_model_weights.h"
#endif
#include "av1/encoder/partitioning.h"
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
      cpi->fn_ptr[bs].vf(ref->buf, ref->stride, AV1_VAR_OFFS, 0, &sse);
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
  var =
      cpi->fn_ptr[bs].vf(ref->buf, ref->stride,
                         CONVERT_TO_BYTEPTR(high_var_offs[off_index]), 0, &sse);
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

#define AVG_CDF_WEIGHT_LEFT 3
#define AVG_CDF_WEIGHT_TOP_RIGHT 1

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
  const BLOCK_SIZE sb_size = cm->seq_params.sb_size;
  const int mib_size = cm->seq_params.mib_size;
  const int mib_size_log2 = cm->seq_params.mib_size_log2;
  const int sb_row = (mi_row - tile_info->mi_row_start) >> mib_size_log2;
  const int use_nonrd_mode = cpi->sf.rt_sf.use_nonrd_pick_mode;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, encode_sb_time);
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
    (*(enc_row_mt->sync_read_ptr))(row_mt_sync, sb_row, sb_col_in_tile);

    if (tile_data->allow_update_cdf && row_mt_enabled &&
        (tile_info->mi_row_start != mi_row)) {
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
    x->color_sensitivity[0] = 0;
    x->color_sensitivity[1] = 0;
    x->content_state_sb = 0;

    xd->cur_frame_force_integer_mv = cm->features.cur_frame_force_integer_mv;
    x->source_variance = UINT_MAX;
    x->simple_motion_pred_sse = UINT_MAX;
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

    // encode the superblock
    if (use_nonrd_mode) {
      av1_encode_nonrd_sb(cpi, td, tile_data, tp, mi_row, mi_col, seg_skip);
    } else {
      av1_encode_rd_sb(cpi, td, tile_data, tp, mi_row, mi_col, seg_skip);
    }

    // Update the top-right context in row_mt coding
    if (tile_data->allow_update_cdf && row_mt_enabled &&
        (tile_info->mi_row_end > (mi_row + mib_size))) {
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
  end_timing(cpi, encode_sb_time);
#endif
}

static AOM_INLINE void init_encode_frame_mb_context(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCK *const x = &cpi->td.mb;
  MACROBLOCKD *const xd = &x->e_mbd;

  // Copy data over into macro block data structures.
  av1_setup_src_planes(x, cpi->source, 0, 0, num_planes,
                       cm->seq_params.sb_size);

  av1_setup_block_planes(xd, cm->seq_params.subsampling_x,
                         cm->seq_params.subsampling_y, num_planes);
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

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileDataEnc *const tile_data =
          &cpi->tile_data[tile_row * tile_cols + tile_col];
      TileInfo *const tile_info = &tile_data->tile_info;
      av1_tile_init(tile_info, cm, tile_row, tile_col);
      tile_data->firstpass_top_mv = kZeroMv;

      if (pre_tok != NULL && tplist != NULL) {
        token_info->tile_tok[tile_row][tile_col] = pre_tok + tile_tok;
        pre_tok = token_info->tile_tok[tile_row][tile_col];
        tile_tok = allocated_tokens(*tile_info,
                                    cm->seq_params.mib_size_log2 + MI_SIZE_LOG2,
                                    num_planes);
        token_info->tplist[tile_row][tile_col] = tplist + tplist_count;
        tplist = token_info->tplist[tile_row][tile_col];
        tplist_count = av1_get_sb_rows_in_tile(cm, tile_data->tile_info);
      }
      tile_data->allow_update_cdf = !cm->tiles.large_scale;
      tile_data->allow_update_cdf =
          tile_data->allow_update_cdf && !cm->features.disable_cdf_update;
      tile_data->tctx = *cm->fc;
    }
  }
}

/*!\brief Encode a superblock row
 *
 * \ingroup partition_search
 */
void av1_encode_sb_row(AV1_COMP *cpi, ThreadData *td, int tile_row,
                       int tile_col, int mi_row) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  const int tile_cols = cm->tiles.cols;
  TileDataEnc *this_tile = &cpi->tile_data[tile_row * tile_cols + tile_col];
  const TileInfo *const tile_info = &this_tile->tile_info;
  TokenExtra *tok = NULL;
  TokenList *const tplist = cpi->token_info.tplist[tile_row][tile_col];
  const int sb_row_in_tile =
      (mi_row - tile_info->mi_row_start) >> cm->seq_params.mib_size_log2;
  const int tile_mb_cols =
      (tile_info->mi_col_end - tile_info->mi_col_start + 2) >> 2;
  const int num_mb_rows_in_sb =
      ((1 << (cm->seq_params.mib_size_log2 + MI_SIZE_LOG2)) + 8) >> 4;

  get_start_tok(cpi, tile_row, tile_col, mi_row, &tok,
                cm->seq_params.mib_size_log2 + MI_SIZE_LOG2, num_planes);
  tplist[sb_row_in_tile].start = tok;

  encode_sb_row(cpi, td, this_tile, mi_row, &tok);

  tplist[sb_row_in_tile].count =
      (unsigned int)(tok - tplist[sb_row_in_tile].start);

  assert((unsigned int)(tok - tplist[sb_row_in_tile].start) <=
         get_token_alloc(num_mb_rows_in_sb, tile_mb_cols,
                         cm->seq_params.mib_size_log2 + MI_SIZE_LOG2,
                         num_planes));

  (void)tile_mb_cols;
  (void)num_mb_rows_in_sb;
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
    cfl_init(&td->mb.e_mbd.cfl, &cm->seq_params);

  av1_crc32c_calculator_init(
      &td->mb.txfm_search_info.mb_rd_record.crc_calculator);

  for (int mi_row = tile_info->mi_row_start; mi_row < tile_info->mi_row_end;
       mi_row += cm->seq_params.mib_size) {
    av1_encode_sb_row(cpi, td, tile_row, tile_col, mi_row);
  }
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

  assert(IMPLIES(cpi->tile_data == NULL,
                 cpi->allocated_tiles < tile_cols * tile_rows));
  if (cpi->allocated_tiles < tile_cols * tile_rows) av1_alloc_tile_data(cpi);

  av1_init_tile_data(cpi);

  for (tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (tile_col = 0; tile_col < tile_cols; ++tile_col) {
      TileDataEnc *const this_tile =
          &cpi->tile_data[tile_row * cm->tiles.cols + tile_col];
      cpi->td.intrabc_used = 0;
      cpi->td.deltaq_used = 0;
      cpi->td.mb.e_mbd.tile_ctx = &this_tile->tctx;
      cpi->td.mb.tile_pb_ctx = &this_tile->tctx;
      av1_encode_tile(cpi, &cpi->td, tile_row, tile_col);
      cpi->intrabc_used |= cpi->td.intrabc_used;
      cpi->deltaq_used |= cpi->td.deltaq_used;
    }
  }
}

// Set the relative distance of a reference frame w.r.t. current frame
static AOM_INLINE void set_rel_frame_dist(
    const AV1_COMMON *const cm, RefFrameDistanceInfo *const ref_frame_dist_info,
    const int ref_frame_flags) {
  const OrderHintInfo *const order_hint_info = &cm->seq_params.order_hint_info;
  MV_REFERENCE_FRAME ref_frame;
  int min_past_dist = INT32_MAX, min_future_dist = INT32_MAX;
  ref_frame_dist_info->nearest_past_ref = NONE_FRAME;
  ref_frame_dist_info->nearest_future_ref = NONE_FRAME;
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    ref_frame_dist_info->ref_relative_dist[ref_frame - LAST_FRAME] = 0;
    if (ref_frame_flags & av1_ref_frame_flag_list[ref_frame]) {
      int dist = av1_encoder_get_relative_dist(
          order_hint_info,
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
  for (int ref = LAST_FRAME; ref <= ALTREF_FRAME; ++ref) {
    const RefCntBuffer *const buf = get_ref_frame_buf(cm, ref);
    if (buf == NULL) continue;

    const int ref_display_order_hint = buf->display_order_hint;
    if (av1_encoder_get_relative_dist(
            &cm->seq_params.order_hint_info, ref_display_order_hint,
            (int)cm->current_frame.display_order_hint) > 0) {
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
  const int cur_to_ref0 = get_relative_dist(&cm->seq_params.order_hint_info,
                                            cur_offset, ref_offset[0]);
  const int cur_to_ref1 = abs(get_relative_dist(&cm->seq_params.order_hint_info,
                                                cur_offset, ref_offset[1]));
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
  if (!cpi->sf.rt_sf.use_nonrd_pick_mode &&
      cpi->sf.inter_sf.selective_ref_frame >= 2) {
    AV1_COMMON *const cm = &cpi->common;
    const OrderHintInfo *const order_hint_info =
        &cm->seq_params.order_hint_info;
    const int cur_frame_display_order_hint =
        cm->current_frame.display_order_hint;
    unsigned int *ref_display_order_hint =
        cm->cur_frame->ref_display_order_hint;
    const int arf2_dist = av1_encoder_get_relative_dist(
        order_hint_info, ref_display_order_hint[ALTREF2_FRAME - LAST_FRAME],
        cur_frame_display_order_hint);
    const int bwd_dist = av1_encoder_get_relative_dist(
        order_hint_info, ref_display_order_hint[BWDREF_FRAME - LAST_FRAME],
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
              order_hint_info, ref_display_order_hint[rf[i] - LAST_FRAME],
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
  FrameProbInfo *const frame_probs = &cpi->frame_probs;
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
  av1_zero(rdc->comp_pred_diff);
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
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);
    if (frame_probs->warped_probs[update_type] <
        cpi->sf.inter_sf.prune_warped_prob_thresh)
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
        (1 << (cm->seq_params.mib_size_log2 + MI_SIZE_LOG2));
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
    // Set delta_q_present_flag before it is used for the first time
    cm->delta_q_info.delta_lf_res = DEFAULT_DELTA_LF_RES;
    cm->delta_q_info.delta_q_present_flag = deltaq_mode != NO_DELTA_Q;

    // Turn off cm->delta_q_info.delta_q_present_flag if objective delta_q
    // is used for ineligible frames. That effectively will turn off row_mt
    // usage. Note objective delta_q and tpl eligible frames are only altref
    // frames currently.
    const GF_GROUP *gf_group = &cpi->gf_group;
    if (cm->delta_q_info.delta_q_present_flag) {
      if (deltaq_mode == DELTA_Q_OBJECTIVE && !is_frame_tpl_eligible(gf_group))
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
  }

  av1_frame_init_quantizer(cpi);
  av1_initialize_rd_consts(cpi);
  av1_set_sad_per_bit(cpi, &x->mv_costs, quant_params->base_qindex);

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

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_compute_global_motion_time);
#endif
  av1_compute_global_motion_facade(cpi);
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_compute_global_motion_time);
#endif

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_setup_motion_field_time);
#endif
  if (features->allow_ref_frame_mvs) av1_setup_motion_field(cm);
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_setup_motion_field_time);
#endif

  cm->current_frame.skip_mode_info.skip_mode_flag =
      check_skip_mode_enabled(cpi);

  enc_row_mt->sync_read_ptr = av1_row_mt_sync_read_dummy;
  enc_row_mt->sync_write_ptr = av1_row_mt_sync_write_dummy;
  mt_info->row_mt_enabled = 0;

  if (oxcf->row_mt && (oxcf->max_threads > 1)) {
    mt_info->row_mt_enabled = 1;
    enc_row_mt->sync_read_ptr = av1_row_mt_sync_read;
    enc_row_mt->sync_write_ptr = av1_row_mt_sync_write;
    av1_encode_tiles_row_mt(cpi);
  } else {
    if (AOMMIN(oxcf->max_threads, cm->tiles.cols * cm->tiles.rows) > 1)
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

  if (cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats) {
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);

    for (i = 0; i < TX_SIZES_ALL; i++) {
      int sum = 0;
      int j;
      int left = 1024;

      for (j = 0; j < TX_TYPES; j++)
        sum += cpi->td.rd_counts.tx_type_used[i][j];

      for (j = TX_TYPES - 1; j >= 0; j--) {
        const int new_prob =
            sum ? 1024 * cpi->td.rd_counts.tx_type_used[i][j] / sum
                : (j ? 0 : 1024);
        int prob =
            (frame_probs->tx_type_probs[update_type][i][j] + new_prob) >> 1;
        left -= prob;
        if (j == 0) prob += left;
        frame_probs->tx_type_probs[update_type][i][j] = prob;
      }
    }
  }

  if (!cpi->sf.inter_sf.disable_obmc &&
      cpi->sf.inter_sf.prune_obmc_prob_thresh > 0) {
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);

    for (i = 0; i < BLOCK_SIZES_ALL; i++) {
      int sum = 0;
      for (int j = 0; j < 2; j++) sum += cpi->td.rd_counts.obmc_used[i][j];

      const int new_prob =
          sum ? 128 * cpi->td.rd_counts.obmc_used[i][1] / sum : 0;
      frame_probs->obmc_probs[update_type][i] =
          (frame_probs->obmc_probs[update_type][i] + new_prob) >> 1;
    }
  }

  if (features->allow_warped_motion &&
      cpi->sf.inter_sf.prune_warped_prob_thresh > 0) {
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);
    int sum = 0;
    for (i = 0; i < 2; i++) sum += cpi->td.rd_counts.warped_used[i];
    const int new_prob = sum ? 128 * cpi->td.rd_counts.warped_used[1] / sum : 0;
    frame_probs->warped_probs[update_type] =
        (frame_probs->warped_probs[update_type] + new_prob) >> 1;
  }

  if (cm->current_frame.frame_type != KEY_FRAME &&
      cpi->sf.interp_sf.adaptive_interp_filter_search == 2 &&
      features->interp_filter == SWITCHABLE) {
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);

    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; i++) {
      int sum = 0;
      int j;
      int left = 1536;

      for (j = 0; j < SWITCHABLE_FILTERS; j++) {
        sum += cpi->td.counts->switchable_interp[i][j];
      }

      for (j = SWITCHABLE_FILTERS - 1; j >= 0; j--) {
        const int new_prob =
            sum ? 1536 * cpi->td.counts->switchable_interp[i][j] / sum
                : (j ? 0 : 1536);
        int prob = (frame_probs->switchable_interp_probs[update_type][i][j] +
                    new_prob) >>
                   1;
        left -= prob;
        if (j == 0) prob += left;
        frame_probs->switchable_interp_probs[update_type][i][j] = prob;
      }
    }
  }

  if ((!is_stat_generation_stage(cpi) && av1_use_hash_me(cpi) &&
       !cpi->sf.rt_sf.use_nonrd_pick_mode) ||
      hash_table_created) {
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
  enforce_max_ref_frames(cpi, &cpi->ref_frame_flags);
  set_rel_frame_dist(&cpi->common, &cpi->ref_frame_dist_info,
                     cpi->ref_frame_flags);
  av1_setup_frame_sign_bias(cm);

#if CONFIG_MISMATCH_DEBUG
  mismatch_reset_frame(num_planes);
#else
  (void)num_planes;
#endif

  if (cpi->sf.hl_sf.frame_parameter_update) {
    RD_COUNTS *const rdc = &cpi->td.rd_counts;

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
    encode_frame_internal(cpi);
  }
}

static AOM_INLINE void update_txfm_count(MACROBLOCK *x, MACROBLOCKD *xd,
                                         FRAME_COUNTS *counts, TX_SIZE tx_size,
                                         int depth, int blk_row, int blk_col,
                                         uint8_t allow_update_cdf) {
  MB_MODE_INFO *mbmi = xd->mi[0];
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int max_blocks_high = max_block_high(xd, bsize, 0);
  const int max_blocks_wide = max_block_wide(xd, bsize, 0);
  int ctx = txfm_partition_context(xd->above_txfm_context + blk_col,
                                   xd->left_txfm_context + blk_row,
                                   mbmi->sb_type, tx_size);
  const int txb_size_index = av1_get_txb_size_index(bsize, blk_row, blk_col);
  const TX_SIZE plane_tx_size = mbmi->inter_tx_size[txb_size_index];

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;
  assert(tx_size > TX_4X4);

  if (depth == MAX_VARTX_DEPTH) {
    // Don't add to counts in this case
    mbmi->tx_size = tx_size;
    txfm_partition_update(xd->above_txfm_context + blk_col,
                          xd->left_txfm_context + blk_row, tx_size, tx_size);
    return;
  }

  if (tx_size == plane_tx_size) {
#if CONFIG_ENTROPY_STATS
    ++counts->txfm_partition[ctx][0];
#endif
    if (allow_update_cdf)
      update_cdf(xd->tile_ctx->txfm_partition_cdf[ctx], 0, 2);
    mbmi->tx_size = tx_size;
    txfm_partition_update(xd->above_txfm_context + blk_col,
                          xd->left_txfm_context + blk_row, tx_size, tx_size);
  } else {
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsw = tx_size_wide_unit[sub_txs];
    const int bsh = tx_size_high_unit[sub_txs];

#if CONFIG_ENTROPY_STATS
    ++counts->txfm_partition[ctx][1];
#endif
    if (allow_update_cdf)
      update_cdf(xd->tile_ctx->txfm_partition_cdf[ctx], 1, 2);
    ++x->txfm_search_info.txb_split_count;

    if (sub_txs == TX_4X4) {
      mbmi->inter_tx_size[txb_size_index] = TX_4X4;
      mbmi->tx_size = TX_4X4;
      txfm_partition_update(xd->above_txfm_context + blk_col,
                            xd->left_txfm_context + blk_row, TX_4X4, tx_size);
      return;
    }

    for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh) {
      for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
        int offsetr = row;
        int offsetc = col;

        update_txfm_count(x, xd, counts, sub_txs, depth + 1, blk_row + offsetr,
                          blk_col + offsetc, allow_update_cdf);
      }
    }
  }
}

static AOM_INLINE void tx_partition_count_update(const AV1_COMMON *const cm,
                                                 MACROBLOCK *x,
                                                 BLOCK_SIZE plane_bsize,
                                                 FRAME_COUNTS *td_counts,
                                                 uint8_t allow_update_cdf) {
  MACROBLOCKD *xd = &x->e_mbd;
  const int mi_width = mi_size_wide[plane_bsize];
  const int mi_height = mi_size_high[plane_bsize];
  const TX_SIZE max_tx_size = get_vartx_max_txsize(xd, plane_bsize, 0);
  const int bh = tx_size_high_unit[max_tx_size];
  const int bw = tx_size_wide_unit[max_tx_size];

  xd->above_txfm_context =
      cm->above_contexts.txfm[xd->tile.tile_row] + xd->mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (xd->mi_row & MAX_MIB_MASK);

  for (int idy = 0; idy < mi_height; idy += bh) {
    for (int idx = 0; idx < mi_width; idx += bw) {
      update_txfm_count(x, xd, td_counts, max_tx_size, 0, idy, idx,
                        allow_update_cdf);
    }
  }
}

static AOM_INLINE void set_txfm_context(MACROBLOCKD *xd, TX_SIZE tx_size,
                                        int blk_row, int blk_col) {
  MB_MODE_INFO *mbmi = xd->mi[0];
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int max_blocks_high = max_block_high(xd, bsize, 0);
  const int max_blocks_wide = max_block_wide(xd, bsize, 0);
  const int txb_size_index = av1_get_txb_size_index(bsize, blk_row, blk_col);
  const TX_SIZE plane_tx_size = mbmi->inter_tx_size[txb_size_index];

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  if (tx_size == plane_tx_size) {
    mbmi->tx_size = tx_size;
    txfm_partition_update(xd->above_txfm_context + blk_col,
                          xd->left_txfm_context + blk_row, tx_size, tx_size);

  } else {
    if (tx_size == TX_8X8) {
      mbmi->inter_tx_size[txb_size_index] = TX_4X4;
      mbmi->tx_size = TX_4X4;
      txfm_partition_update(xd->above_txfm_context + blk_col,
                            xd->left_txfm_context + blk_row, TX_4X4, tx_size);
      return;
    }
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsw = tx_size_wide_unit[sub_txs];
    const int bsh = tx_size_high_unit[sub_txs];
    for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh) {
      for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
        const int offsetr = blk_row + row;
        const int offsetc = blk_col + col;
        if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide) continue;
        set_txfm_context(xd, sub_txs, offsetr, offsetc);
      }
    }
  }
}

static AOM_INLINE void tx_partition_set_contexts(const AV1_COMMON *const cm,
                                                 MACROBLOCKD *xd,
                                                 BLOCK_SIZE plane_bsize) {
  const int mi_width = mi_size_wide[plane_bsize];
  const int mi_height = mi_size_high[plane_bsize];
  const TX_SIZE max_tx_size = get_vartx_max_txsize(xd, plane_bsize, 0);
  const int bh = tx_size_high_unit[max_tx_size];
  const int bw = tx_size_wide_unit[max_tx_size];

  xd->above_txfm_context =
      cm->above_contexts.txfm[xd->tile.tile_row] + xd->mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (xd->mi_row & MAX_MIB_MASK);

  for (int idy = 0; idy < mi_height; idy += bh) {
    for (int idx = 0; idx < mi_width; idx += bw) {
      set_txfm_context(xd, max_tx_size, idy, idx);
    }
  }
}

static void update_zeromv_cnt(const AV1_COMP *const cpi,
                              const MB_MODE_INFO *const mi, int mi_row,
                              int mi_col, BLOCK_SIZE bsize) {
  const AV1_COMMON *const cm = &cpi->common;
  MV mv = mi->mv[0].as_mv;
  const int bw = mi_size_wide[bsize] >> 1;
  const int bh = mi_size_high[bsize] >> 1;
  const int xmis = AOMMIN((cm->mi_params.mi_cols - mi_col) >> 1, bw);
  const int ymis = AOMMIN((cm->mi_params.mi_rows - mi_row) >> 1, bh);
  const int block_index =
      (mi_row >> 1) * (cm->mi_params.mi_cols >> 1) + (mi_col >> 1);
  for (int y = 0; y < ymis; y++)
    for (int x = 0; x < xmis; x++) {
      // consec_zero_mv is in the scale of 8x8 blocks
      const int map_offset = block_index + y * (cm->mi_params.mi_cols >> 1) + x;
      if (mi->ref_frame[0] == LAST_FRAME && is_inter_block(mi) &&
          mi->segment_id <= CR_SEGMENT_ID_BOOST2) {
        if (abs(mv.row) < 10 && abs(mv.col) < 10) {
          if (cpi->consec_zero_mv[map_offset] < 255)
            cpi->consec_zero_mv[map_offset]++;
        } else {
          cpi->consec_zero_mv[map_offset] = 0;
        }
      }
    }
}

void av1_encode_superblock(const AV1_COMP *const cpi, TileDataEnc *tile_data,
                           ThreadData *td, TokenExtra **t, RUN_TYPE dry_run,
                           BLOCK_SIZE bsize, int *rate) {
  const AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO **mi_4x4 = xd->mi;
  MB_MODE_INFO *mbmi = mi_4x4[0];
  const int seg_skip =
      segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP);
  const int mis = cm->mi_params.mi_stride;
  const int mi_width = mi_size_wide[bsize];
  const int mi_height = mi_size_high[bsize];
  const int is_inter = is_inter_block(mbmi);

  // Initialize tx_mode and tx_size_search_method
  TxfmSearchParams *txfm_params = &x->txfm_search_params;
  set_tx_size_search_method(
      cm, &cpi->winner_mode_params, txfm_params,
      cpi->sf.winner_mode_sf.enable_winner_mode_for_tx_size_srch, 1);

  const int mi_row = xd->mi_row;
  const int mi_col = xd->mi_col;
  if (!is_inter) {
    xd->cfl.store_y = store_cfl_required(cm, xd);
    mbmi->skip_txfm = 1;
    for (int plane = 0; plane < num_planes; ++plane) {
      av1_encode_intra_block_plane(cpi, x, bsize, plane, dry_run,
                                   cpi->optimize_seg_arr[mbmi->segment_id]);
    }

    // If there is at least one lossless segment, force the skip for intra
    // block to be 0, in order to avoid the segment_id to be changed by in
    // write_segment_id().
    if (!cpi->common.seg.segid_preskip && cpi->common.seg.update_map &&
        cpi->enc_seg.has_lossless_segment)
      mbmi->skip_txfm = 0;

    xd->cfl.store_y = 0;
    if (av1_allow_palette(cm->features.allow_screen_content_tools, bsize)) {
      for (int plane = 0; plane < AOMMIN(2, num_planes); ++plane) {
        if (mbmi->palette_mode_info.palette_size[plane] > 0) {
          if (!dry_run) {
            av1_tokenize_color_map(x, plane, t, bsize, mbmi->tx_size,
                                   PALETTE_MAP, tile_data->allow_update_cdf,
                                   td->counts);
          } else if (dry_run == DRY_RUN_COSTCOEFFS) {
            rate +=
                av1_cost_color_map(x, plane, bsize, mbmi->tx_size, PALETTE_MAP);
          }
        }
      }
    }

    av1_update_txb_context(cpi, td, dry_run, bsize,
                           tile_data->allow_update_cdf);
  } else {
    int ref;
    const int is_compound = has_second_ref(mbmi);

    set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);
    for (ref = 0; ref < 1 + is_compound; ++ref) {
      const YV12_BUFFER_CONFIG *cfg =
          get_ref_frame_yv12_buf(cm, mbmi->ref_frame[ref]);
      assert(IMPLIES(!is_intrabc_block(mbmi), cfg));
      av1_setup_pre_planes(xd, ref, cfg, mi_row, mi_col,
                           xd->block_ref_scale_factors[ref], num_planes);
    }
    int start_plane = (cpi->sf.rt_sf.reuse_inter_pred_nonrd) ? 1 : 0;
    av1_enc_build_inter_predictor(cm, xd, mi_row, mi_col, NULL, bsize,
                                  start_plane, av1_num_planes(cm) - 1);
    if (mbmi->motion_mode == OBMC_CAUSAL) {
      assert(cpi->oxcf.motion_mode_cfg.enable_obmc);
      av1_build_obmc_inter_predictors_sb(cm, xd);
    }

#if CONFIG_MISMATCH_DEBUG
    if (dry_run == OUTPUT_ENABLED) {
      for (int plane = 0; plane < num_planes; ++plane) {
        const struct macroblockd_plane *pd = &xd->plane[plane];
        int pixel_c, pixel_r;
        mi_to_pixel_loc(&pixel_c, &pixel_r, mi_col, mi_row, 0, 0,
                        pd->subsampling_x, pd->subsampling_y);
        if (!is_chroma_reference(mi_row, mi_col, bsize, pd->subsampling_x,
                                 pd->subsampling_y))
          continue;
        mismatch_record_block_pre(pd->dst.buf, pd->dst.stride,
                                  cm->current_frame.order_hint, plane, pixel_c,
                                  pixel_r, pd->width, pd->height,
                                  xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH);
      }
    }
#else
    (void)num_planes;
#endif

    av1_encode_sb(cpi, x, bsize, dry_run);
    av1_tokenize_sb_vartx(cpi, td, dry_run, bsize, rate,
                          tile_data->allow_update_cdf);
  }

  if (!dry_run) {
    if (av1_allow_intrabc(cm) && is_intrabc_block(mbmi)) td->intrabc_used = 1;
    if (txfm_params->tx_mode_search_type == TX_MODE_SELECT &&
        !xd->lossless[mbmi->segment_id] && mbmi->sb_type > BLOCK_4X4 &&
        !(is_inter && (mbmi->skip_txfm || seg_skip))) {
      if (is_inter) {
        tx_partition_count_update(cm, x, bsize, td->counts,
                                  tile_data->allow_update_cdf);
      } else {
        if (mbmi->tx_size != max_txsize_rect_lookup[bsize])
          ++x->txfm_search_info.txb_split_count;
        if (block_signals_txsize(bsize)) {
          const int tx_size_ctx = get_tx_size_context(xd);
          const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);
          const int depth = tx_size_to_depth(mbmi->tx_size, bsize);
          const int max_depths = bsize_to_max_depth(bsize);

          if (tile_data->allow_update_cdf)
            update_cdf(xd->tile_ctx->tx_size_cdf[tx_size_cat][tx_size_ctx],
                       depth, max_depths + 1);
#if CONFIG_ENTROPY_STATS
          ++td->counts->intra_tx_size[tx_size_cat][tx_size_ctx][depth];
#endif
        }
      }
      assert(IMPLIES(is_rect_tx(mbmi->tx_size), is_rect_tx_allowed(xd, mbmi)));
    } else {
      int i, j;
      TX_SIZE intra_tx_size;
      // The new intra coding scheme requires no change of transform size
      if (is_inter) {
        if (xd->lossless[mbmi->segment_id]) {
          intra_tx_size = TX_4X4;
        } else {
          intra_tx_size =
              tx_size_from_tx_mode(bsize, txfm_params->tx_mode_search_type);
        }
      } else {
        intra_tx_size = mbmi->tx_size;
      }

      for (j = 0; j < mi_height; j++)
        for (i = 0; i < mi_width; i++)
          if (mi_col + i < cm->mi_params.mi_cols &&
              mi_row + j < cm->mi_params.mi_rows)
            mi_4x4[mis * j + i]->tx_size = intra_tx_size;

      if (intra_tx_size != max_txsize_rect_lookup[bsize])
        ++x->txfm_search_info.txb_split_count;
    }
  }

  if (txfm_params->tx_mode_search_type == TX_MODE_SELECT &&
      block_signals_txsize(mbmi->sb_type) && is_inter &&
      !(mbmi->skip_txfm || seg_skip) && !xd->lossless[mbmi->segment_id]) {
    if (dry_run) tx_partition_set_contexts(cm, xd, bsize);
  } else {
    TX_SIZE tx_size = mbmi->tx_size;
    // The new intra coding scheme requires no change of transform size
    if (is_inter) {
      if (xd->lossless[mbmi->segment_id]) {
        tx_size = TX_4X4;
      } else {
        tx_size = tx_size_from_tx_mode(bsize, txfm_params->tx_mode_search_type);
      }
    } else {
      tx_size = (bsize > BLOCK_4X4) ? tx_size : TX_4X4;
    }
    mbmi->tx_size = tx_size;
    set_txfm_ctxs(tx_size, xd->width, xd->height,
                  (mbmi->skip_txfm || seg_skip) && is_inter_block(mbmi), xd);
  }

  if (is_inter_block(mbmi) && !xd->is_chroma_ref && is_cfl_allowed(xd)) {
    cfl_store_block(xd, mbmi->sb_type, mbmi->tx_size);
  }
  if (!dry_run) {
    if (cpi->oxcf.pass == 0 && cpi->svc.temporal_layer_id == 0 &&
        cpi->sf.rt_sf.use_temporal_noise_estimate &&
        (!cpi->use_svc ||
         (cpi->use_svc &&
          !cpi->svc.layer_context[cpi->svc.temporal_layer_id].is_key_frame &&
          cpi->svc.spatial_layer_id == cpi->svc.number_spatial_layers - 1)))
      update_zeromv_cnt(cpi, mbmi, mi_row, mi_col, bsize);
  }
}
