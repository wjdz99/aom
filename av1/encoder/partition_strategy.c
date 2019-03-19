/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom_ports/system_state.h"

#include "av1/common/enums.h"
#include "av1/common/reconinter.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/partition_model_weights.h"
#include "av1/encoder/partition_strategy.h"
#include "av1/encoder/reconinter_enc.h"

void av1_simple_motion_search(AV1_COMP *const cpi, MACROBLOCK *x, int mi_row,
                              int mi_col, BLOCK_SIZE bsize, int ref,
                              MV ref_mv_full, int num_planes,
                              int use_subpixel) {
  assert(num_planes == 1 &&
         "Currently simple_motion_search only supports luma plane");
  assert(!frame_is_intra_only(&cpi->common) &&
         "Simple motion search only enabled for non-key frames");
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;

  set_offsets_for_motion_search(cpi, x, mi_row, mi_col, bsize);

  MB_MODE_INFO *mbmi = xd->mi[0];
  mbmi->sb_type = bsize;
  mbmi->ref_frame[0] = ref;
  mbmi->ref_frame[1] = NONE_FRAME;
  mbmi->motion_mode = SIMPLE_TRANSLATION;

  const YV12_BUFFER_CONFIG *yv12 = get_ref_frame_yv12_buf(cm, ref);
  const YV12_BUFFER_CONFIG *scaled_ref_frame =
      av1_get_scaled_ref_frame(cpi, ref);
  struct buf_2d backup_yv12;
  // ref_mv is used to code the motion vector. ref_mv_full is the initial point.
  // ref_mv is in units of 1/8 pel whereas ref_mv_full is in units of pel.
  MV ref_mv = { 0, 0 };
  const int step_param = cpi->mv_step_param;
  const MvLimits tmp_mv_limits = x->mv_limits;
  const SEARCH_METHODS search_methods = NSTEP;
  const int do_mesh_search = 0;
  const int sadpb = x->sadperbit16;
  int cost_list[5];
  const int ref_idx = 0;
  int var;

  if (scaled_ref_frame) {
    backup_yv12 = xd->plane[AOM_PLANE_Y].pre[ref_idx];
    av1_setup_pre_planes(xd, ref_idx, scaled_ref_frame, mi_row, mi_col, NULL,
                         num_planes);
  } else {
    av1_setup_pre_planes(xd, ref_idx, yv12, mi_row, mi_col,
                         get_ref_scale_factors(cm, ref), num_planes);
  }

  // This overwrites the mv_limits so we will need to restore it later.
  av1_set_mv_search_range(&x->mv_limits, &ref_mv);
  var = av1_full_pixel_search(
      cpi, x, bsize, &ref_mv_full, step_param, search_methods, do_mesh_search,
      sadpb, cond_cost_list(cpi, cost_list), &ref_mv, INT_MAX, 1,
      mi_col * MI_SIZE, mi_row * MI_SIZE, 0, &cpi->ss_cfg[SS_CFG_SRC]);
  // Restore
  x->mv_limits = tmp_mv_limits;

  const int use_subpel_search =
      var < INT_MAX && !cpi->common.cur_frame_force_integer_mv && use_subpixel;
  if (use_subpel_search) {
    int not_used = 0;
    if (cpi->sf.use_accurate_subpel_search) {
      const int pw = block_size_wide[bsize];
      const int ph = block_size_high[bsize];
      cpi->find_fractional_mv_step(
          x, cm, mi_row, mi_col, &ref_mv, cm->allow_high_precision_mv,
          x->errorperbit, &cpi->fn_ptr[bsize], cpi->sf.mv.subpel_force_stop,
          cpi->sf.mv.subpel_iters_per_step, cond_cost_list(cpi, cost_list),
          x->nmv_vec_cost, x->mv_cost_stack, &not_used, &x->pred_sse[ref], NULL,
          NULL, 0, 0, pw, ph, cpi->sf.use_accurate_subpel_search, 1);
    } else {
      cpi->find_fractional_mv_step(
          x, cm, mi_row, mi_col, &ref_mv, cm->allow_high_precision_mv,
          x->errorperbit, &cpi->fn_ptr[bsize], cpi->sf.mv.subpel_force_stop,
          cpi->sf.mv.subpel_iters_per_step, cond_cost_list(cpi, cost_list),
          x->nmv_vec_cost, x->mv_cost_stack, &not_used, &x->pred_sse[ref], NULL,
          NULL, 0, 0, 0, 0, 0, 1);
    }
  } else {
    // Manually convert from units of pixel to 1/8-pixels if we are not doing
    // subpel search
    x->best_mv.as_mv.row *= 8;
    x->best_mv.as_mv.col *= 8;
  }

  mbmi->mv[0].as_mv = x->best_mv.as_mv;

  // Get a copy of the prediction output
  set_ref_ptrs(cm, xd, mbmi->ref_frame[0], mbmi->ref_frame[1]);
  av1_enc_build_inter_predictor(cm, xd, mi_row, mi_col, NULL, bsize,
                                AOM_PLANE_Y, AOM_PLANE_Y);

  aom_clear_system_state();

  if (scaled_ref_frame) {
    xd->plane[AOM_PLANE_Y].pre[ref_idx] = backup_yv12;
  }
}

void av1_simple_motion_sse_var(AV1_COMP *cpi, MACROBLOCK *x, int mi_row,
                               int mi_col, BLOCK_SIZE bsize,
                               const MV ref_mv_full, int use_subpixel,
                               unsigned int *sse, unsigned int *var) {
  MACROBLOCKD *xd = &x->e_mbd;
  const MV_REFERENCE_FRAME ref =
      cpi->rc.is_src_frame_alt_ref ? ALTREF_FRAME : LAST_FRAME;

  av1_simple_motion_search(cpi, x, mi_row, mi_col, bsize, ref, ref_mv_full, 1,
                           use_subpixel);

  const uint8_t *src = x->plane[0].src.buf;
  const int src_stride = x->plane[0].src.stride;
  const uint8_t *dst = xd->plane[0].dst.buf;
  const int dst_stride = xd->plane[0].dst.stride;

  *var = cpi->fn_ptr[bsize].vf(src, src_stride, dst, dst_stride, sse);
}

// Performs a simple_motion_search with a single reference frame and extract
// the variance of residues. Here features is assumed to be a length 6 array.
// After this function is called, we will store the following in to features:
// features[0] = log(1 + dc_q**2/256)
// features[1] = log(1 + variance_of_residue)
// for i in [2, 3, 4, 5]:
//  features[i] = log(1 + variance_of_residue_in_block[i]/variance_of_residue)
static void get_res_var_features(AV1_COMP *const cpi, MACROBLOCK *x, int mi_row,
                                 int mi_col, BLOCK_SIZE bsize,
                                 float *features) {
  // TODO(chiyotsai@google.com): The data this model trained on did not also use
  // SIMPLE_TRANSLATION to build the inter_predictor. Retraining and tuning the
  // model with the correct data should give better performance.
  assert(mi_size_wide[bsize] == mi_size_high[bsize]);

  MACROBLOCKD *xd = &x->e_mbd;

  // Perform a single motion search in Y_PLANE to make a prediction
  const int use_subpixel = 0;

  // Start getting the features
  int f_idx = 0;

  // Q_INDEX
  const int dc_q = av1_dc_quant_QTX(x->qindex, 0, xd->bd) >> (xd->bd - 8);
  aom_clear_system_state();
  features[f_idx++] = logf(1.0f + (float)(dc_q * dc_q) / 256.0f);

  // VARIANCE
  unsigned int sse = 0;
  unsigned int var = 0;
  const MV ref_mv_full = { .row = 0, .col = 0 };
  av1_simple_motion_sse_var(cpi, x, mi_row, mi_col, bsize, ref_mv_full,
                            use_subpixel, &sse, &var);
  aom_clear_system_state();
  features[f_idx++] = logf(1.0f + (float)var);

  // Regional
  const uint8_t *src = x->plane[0].src.buf;
  const int src_stride = x->plane[0].src.stride;
  const uint8_t *dst = xd->plane[0].dst.buf;
  const int dst_stride = xd->plane[0].dst.stride;
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, PARTITION_SPLIT);
  int r_idx = 0;
  for (r_idx = 0; r_idx < 4; r_idx++) {
    const int x_idx = (r_idx & 1) * bw / 2;
    const int y_idx = (r_idx >> 1) * bh / 2;
    const int src_offset = y_idx * src_stride + x_idx;
    const int dst_offset = y_idx * dst_stride + x_idx;
    const unsigned int sub_var = cpi->fn_ptr[subsize].vf(
        src + src_offset, src_stride, dst + dst_offset, dst_stride, &sse);
    aom_clear_system_state();
    const float var_ratio = (1.0f + (float)sub_var) / (4.0f + (float)var);
    features[f_idx++] = var_ratio;
  }
}

void av1_simple_motion_search_based_split(
    AV1_COMP *const cpi, MACROBLOCK *x, int mi_row, int mi_col,
    BLOCK_SIZE bsize, int *partition_none_allowed, int *partition_horz_allowed,
    int *partition_vert_allowed, int *do_rectangular_split) {
  const NN_CONFIG *nn_config = NULL;
  float split_only_thresh = 0.0f;
  if (bsize == BLOCK_128X128) {
    nn_config = &av1_simple_motion_search_based_split_nn_config_128;
    split_only_thresh = av1_simple_motion_search_based_split_thresh_128;
  } else if (bsize == BLOCK_64X64) {
    nn_config = &av1_simple_motion_search_based_split_nn_config_64;
    split_only_thresh = av1_simple_motion_search_based_split_thresh_64;
  } else if (bsize == BLOCK_32X32) {
    nn_config = &av1_simple_motion_search_based_split_nn_config_32;
    split_only_thresh = av1_simple_motion_search_based_split_thresh_32;
  } else if (bsize == BLOCK_16X16) {
    nn_config = &av1_simple_motion_search_based_split_nn_config_16;
    split_only_thresh = av1_simple_motion_search_based_split_thresh_16;
  } else if (bsize == BLOCK_8X8) {
    // Disable BLOCK_8X8 for now
#if !CONFIG_DISABLE_FULL_PIXEL_SPLIT_8X8
    nn_config = &av1_simple_motion_search_based_split_nn_config_8;
    split_only_thresh = av1_simple_motion_search_based_split_thresh_8;
#endif
  } else {
    assert(0 && "Unexpected block size in simple_motion_based_split");
  }
  if (nn_config) {
    float features[6] = { 0 };
    float score = 0;
    get_res_var_features(cpi, x, mi_row, mi_col, bsize, features);
    av1_nn_predict(features, nn_config, &score);

    if (score > split_only_thresh) {
      *partition_none_allowed = 0;
      *partition_horz_allowed = 0;
      *partition_vert_allowed = 0;
      *do_rectangular_split = 0;
    }
  }
}
