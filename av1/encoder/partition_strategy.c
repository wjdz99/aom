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

#include "av1/encoder/encodeframe.h"
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

// Given a list of ref frames in refs, performs simple_motion_search on each of
// the refs and returns the ref with the smallest sse. Returns -1 if none of the
// ref in the list is available. Also stores the best sse and var in best_sse,
// best_var, respectively. If save_mv_code is -1, don't update mv_ref_fulls in
// pc_tree. If save_mv_code is between 0 and 3, update mv_ref_fulls under
// pc_tree->split[i]. If save_mv_code is 4, update mv_ref_fulls under pc_tree.
static int simple_motion_search_get_best_ref(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, const int *const refs, int num_refs,
    int use_subpixel, int save_mv_code, unsigned int *best_sse,
    unsigned int *best_var) {
  // TODO(chiyotsai@google.com): The calculation of variance currently uses
  // bsize, so we might take area outside of the image into account. We need to
  // modify the SIMD functions to fix this later.
  const AV1_COMMON *const cm = &cpi->common;
  int best_ref = -1;

  if (mi_col >= cm->mi_cols || mi_row >= cm->mi_rows) {
    // If the whole block is outside of the image, set the var and sse to 0.
    *best_var = 0;
    *best_sse = 0;

    return best_ref;
  }

  // Otherwise do loop through the reference frames and find the one with the
  // minimum SSE
  const MACROBLOCKD *xd = &x->e_mbd;
  const MV *mv_ref_fulls = pc_tree->mv_ref_fulls;

  const int num_planes = 1;

  *best_sse = INT_MAX;

  for (int ref_idx = 0; ref_idx < num_refs; ref_idx++) {
    const int ref = refs[ref_idx];

    if (cpi->ref_frame_flags & av1_ref_frame_flag_list[ref]) {
      unsigned int curr_sse = 0, curr_var = 0;
      av1_simple_motion_search(cpi, x, mi_row, mi_col, bsize, ref,
                               mv_ref_fulls[ref], num_planes, use_subpixel);
      curr_var = cpi->fn_ptr[bsize].vf(
          x->plane[0].src.buf, x->plane[0].src.stride, xd->plane[0].dst.buf,
          xd->plane[0].dst.stride, &curr_sse);
      if (curr_sse < *best_sse) {
        *best_sse = curr_sse;
        *best_var = curr_var;
        best_ref = ref;
      }

      const int new_mv_row = x->best_mv.as_mv.row / 8;
      const int new_mv_col = x->best_mv.as_mv.col / 8;
      if (save_mv_code == 4) {
        pc_tree->mv_ref_fulls[ref].row = new_mv_row;
        pc_tree->mv_ref_fulls[ref].col = new_mv_col;
      } else if (save_mv_code >= 0 && save_mv_code < 4) {
        // Propagate the new motion vectors to a lower level
        pc_tree->split[save_mv_code]->mv_ref_fulls[ref].row = new_mv_row;
        pc_tree->split[save_mv_code]->mv_ref_fulls[ref].col = new_mv_col;
      } else {
        assert(save_mv_code == -1 &&
               "Unknown code in simple_motion_search_get_best_ref.");
      }
    }
  }

  return best_ref;
}

// Performs fullpixel simple_motion_search with LAST_FRAME and ALTREF_FRAME on
// each subblock and extract the variance and sse of residues. Then store the
// var and sse from each partition subblock to features. The DC qindex is also
// stored in features.
// Here features is assumed to be a length 19 array.
// After this function is called, we will store the following to features:
// features[0:17] = var and sse from subblocks
// features[18] = DC q_index
#define NUM_FEATURES 25
static void simple_motion_search_prune_part_features(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, float *features) {
  // TODO(chiyotsai@google.com): Cache the result of the motion search from the
  // larger bsize.
  const int w_mi = mi_size_wide[bsize];
  const int h_mi = mi_size_high[bsize];
  int f_idx = 0;
  assert(mi_size_wide[bsize] == mi_size_high[bsize]);
  assert(cpi->ref_frame_flags & av1_ref_frame_flag_list[LAST_FRAME] ||
         cpi->ref_frame_flags & av1_ref_frame_flag_list[ALTREF_FRAME]);

  // Setting up motion search
  const int ref_list[] = { LAST_FRAME, ALTREF_FRAME };
  const int num_refs = 2;
  const int use_subpixel = 1;

  unsigned int int_features[NUM_FEATURES - 1];

  // Doing whole block first to update the mv
  simple_motion_search_get_best_ref(
      cpi, x, pc_tree, mi_row, mi_col, bsize, ref_list, num_refs, use_subpixel,
      4, &int_features[f_idx], &int_features[f_idx + 1]);
  f_idx += 2;

  // Split subblocks
  BLOCK_SIZE subsize = get_partition_subsize(bsize, PARTITION_SPLIT);
  int r_idx = 0;
  for (r_idx = 0; r_idx < 4; r_idx++) {
    const int sub_mi_col = mi_col + (r_idx & 1) * w_mi / 2;
    const int sub_mi_row = mi_row + (r_idx >> 1) * h_mi / 2;

    simple_motion_search_get_best_ref(
        cpi, x, pc_tree, sub_mi_row, sub_mi_col, subsize, ref_list, num_refs,
        use_subpixel, r_idx, &int_features[f_idx], &int_features[f_idx + 1]);
    f_idx += 2;
  }

  // Horz subblocks
  subsize = get_partition_subsize(bsize, PARTITION_HORZ);
  for (r_idx = 0; r_idx < 2; r_idx++) {
    const int sub_mi_col = mi_col + 0;
    const int sub_mi_row = mi_row + r_idx * h_mi / 2;

    simple_motion_search_get_best_ref(
        cpi, x, pc_tree, sub_mi_row, sub_mi_col, subsize, ref_list, num_refs,
        use_subpixel, -1, &int_features[f_idx], &int_features[f_idx + 1]);

    f_idx += 2;
  }

  // Vert subblock
  subsize = get_partition_subsize(bsize, PARTITION_VERT);
  for (r_idx = 0; r_idx < 2; r_idx++) {
    const int sub_mi_col = mi_col + r_idx * w_mi / 2;
    const int sub_mi_row = mi_row + 0;

    simple_motion_search_get_best_ref(
        cpi, x, pc_tree, sub_mi_row, sub_mi_col, subsize, ref_list, num_refs,
        use_subpixel, -1, &int_features[f_idx], &int_features[f_idx + 1]);

    f_idx += 2;
  }

  aom_clear_system_state();
  for (int idx = 0; idx < f_idx; idx++) {
    features[idx] = logf(1.0f + (float)int_features[idx]);
  }

  const MACROBLOCKD *xd = &x->e_mbd;
  set_offsets_for_motion_search(cpi, x, mi_row, mi_col, bsize);

  // Q_INDEX
  const int dc_q = av1_dc_quant_QTX(x->qindex, 0, xd->bd) >> (xd->bd - 8);
  features[f_idx++] = logf(1.0f + (float)(dc_q * dc_q) / 256.0f);

  // Neighbor stuff
  const int has_above = !!xd->above_mbmi;
  const int has_left = !!xd->left_mbmi;
  const BLOCK_SIZE above_bsize = has_above ? xd->above_mbmi->sb_type : bsize;
  const BLOCK_SIZE left_bsize = has_left ? xd->left_mbmi->sb_type : bsize;
  features[f_idx++] = (float)has_above;
  features[f_idx++] = (float)mi_size_wide_log2[above_bsize];
  features[f_idx++] = (float)mi_size_high_log2[above_bsize];
  features[f_idx++] = (float)has_left;
  features[f_idx++] = (float)mi_size_wide_log2[left_bsize];
  features[f_idx++] = (float)mi_size_high_log2[left_bsize];

  assert(f_idx == NUM_FEATURES);
}

#define MAX_NUM_CLASSES 10
void av1_simple_motion_search_prune_part(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, int *partition_none_allowed,
    int *partition_horz_allowed, int *partition_vert_allowed,
    int *do_square_split, int *do_rectangular_split, int *prune_horz,
    int *prune_vert, float *features, int *valid) {
  const AV1_COMMON *const cm = &cpi->common;
  // Get model parameters
  const NN_CONFIG *nn_config = NULL;
  const float *prune_thresh = NULL, *only_thresh = NULL;
  const float *ml_mean = NULL, *ml_std = NULL;
  float normalized_features[NUM_FEATURES] = { 0.0f };

  if (bsize == BLOCK_128X128) {
    nn_config = &av1_simple_motion_search_prune_part_nn_config_128;
    ml_mean = av1_simple_motion_search_prune_part_mean_128;
    ml_std = av1_simple_motion_search_prune_part_std_128;
    prune_thresh = av1_simple_motion_search_prune_part_prune_thresh_128;
    only_thresh = av1_simple_motion_search_prune_part_only_thresh_128;
  } else if (bsize == BLOCK_64X64) {
    nn_config = &av1_simple_motion_search_prune_part_nn_config_64;
    ml_mean = av1_simple_motion_search_prune_part_mean_64;
    ml_std = av1_simple_motion_search_prune_part_std_64;
    prune_thresh = av1_simple_motion_search_prune_part_prune_thresh_64;
    only_thresh = av1_simple_motion_search_prune_part_only_thresh_64;
  } else if (bsize == BLOCK_32X32) {
    nn_config = &av1_simple_motion_search_prune_part_nn_config_32;
    ml_mean = av1_simple_motion_search_prune_part_mean_32;
    ml_std = av1_simple_motion_search_prune_part_std_32;
    prune_thresh = av1_simple_motion_search_prune_part_prune_thresh_32;
    only_thresh = av1_simple_motion_search_prune_part_only_thresh_32;
  } else if (bsize == BLOCK_16X16) {
    nn_config = &av1_simple_motion_search_prune_part_nn_config_16;
    ml_mean = av1_simple_motion_search_prune_part_mean_16;
    ml_std = av1_simple_motion_search_prune_part_std_16;
    prune_thresh = av1_simple_motion_search_prune_part_prune_thresh_16;
    only_thresh = av1_simple_motion_search_prune_part_only_thresh_16;
  } else if (bsize == BLOCK_8X8) {
    nn_config = &av1_simple_motion_search_prune_part_nn_config_8;
    ml_mean = av1_simple_motion_search_prune_part_mean_8;
    ml_std = av1_simple_motion_search_prune_part_std_8;
    prune_thresh = av1_simple_motion_search_prune_part_prune_thresh_8;
    only_thresh = av1_simple_motion_search_prune_part_only_thresh_8;
  } else {
    assert(0 && "Unexpected block size in simple_motion_prune_part");
  }

  // If there is no valid threshold, return immediately.
  if (!nn_config || (prune_thresh[PARTITION_HORZ] == 0.0f &&
                     prune_thresh[PARTITION_VERT] == 0.0f)) {
    return;
  }
  if (bsize < BLOCK_8X8) {
    return;
  }

  // Get features
  simple_motion_search_prune_part_features(cpi, x, pc_tree, mi_row, mi_col,
                                           bsize, features);
  *valid = 1;
  for (int f_idx = 0; f_idx < NUM_FEATURES; f_idx++) {
    normalized_features[f_idx] =
        (features[f_idx] - ml_mean[f_idx]) / ml_std[f_idx];
  }

  // Get probabilities
  float scores[MAX_NUM_CLASSES] = { 0.0f }, probs[MAX_NUM_CLASSES] = { 0.0f };
  const int num_classes =
      (bsize == BLOCK_128X128 || bsize == BLOCK_8X8) ? 4 : 10;

  av1_nn_predict(normalized_features, nn_config, scores);
  aom_clear_system_state();

  av1_nn_softmax(scores, probs, num_classes);

  // Determine if we should prune rectangular partitions.
  if (cpi->sf.simple_motion_search_prune_rect && !frame_is_intra_only(cm) &&
      (*partition_horz_allowed || *partition_vert_allowed) &&
      bsize >= BLOCK_8X8 && !av1_superres_scaled(cm)) {
    *prune_horz = probs[PARTITION_HORZ] <= prune_thresh[PARTITION_HORZ];
    *prune_vert = probs[PARTITION_VERT] <= prune_thresh[PARTITION_VERT];
  }

  // Silence compiler warnings
  (void)only_thresh;
  (void)partition_none_allowed;
  (void)do_square_split;
  (void)do_rectangular_split;
}
#undef MAX_NUM_CLASSES
#undef NUM_FEATURES

// Early terminates PARTITION_NONE using simple_motion_search features and the
// rate, distortion, and rdcost of PARTITION_NONE. This is only called when:
//  - The frame is a show frame
//  - The frame is not intra only
//  - The current bsize is > BLOCK_8X8
//  - blk_row + blk_height/2 < total_rows and blk_col + blk_width/2 < total_cols
#define NUM_FEATURES 28
void av1_simple_motion_search_early_term_none(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, const RD_STATS *none_rdc,
    int *early_terminate, float *simple_motion_features,
    int *simple_motion_features_are_valid) {
  // TODO(chiyotsai@google.com): There are other features we can extract from
  // PARTITION_NONE. Play with this later.
  int f_idx = 0;
  if (!*simple_motion_features_are_valid) {
    simple_motion_search_prune_part_features(cpi, x, pc_tree, mi_row, mi_col,
                                             bsize, simple_motion_features);
    *simple_motion_features_are_valid = 1;
  }
  f_idx = 25;

  simple_motion_features[f_idx++] = logf(1.0f + (float)none_rdc->rate);
  simple_motion_features[f_idx++] = logf(1.0f + (float)none_rdc->dist);
  simple_motion_features[f_idx++] = logf(1.0f + (float)none_rdc->rdcost);

  assert(f_idx == NUM_FEATURES);

  const float *ml_mean = NULL;
  const float *ml_std = NULL;
  const float *ml_model = NULL;

  if (bsize == BLOCK_128X128) {
    ml_mean = av1_simple_motion_search_term_none_mean_128;
    ml_std = av1_simple_motion_search_term_none_std_128;
    ml_model = av1_simple_motion_search_term_none_model_128;
  } else if (bsize == BLOCK_64X64) {
    ml_mean = av1_simple_motion_search_term_none_mean_64;
    ml_std = av1_simple_motion_search_term_none_std_64;
    ml_model = av1_simple_motion_search_term_none_model_64;
  } else if (bsize == BLOCK_32X32) {
    ml_mean = av1_simple_motion_search_term_none_mean_32;
    ml_std = av1_simple_motion_search_term_none_std_32;
    ml_model = av1_simple_motion_search_term_none_model_32;
  } else if (bsize == BLOCK_16X16) {
    ml_mean = av1_simple_motion_search_term_none_mean_16;
    ml_std = av1_simple_motion_search_term_none_std_16;
    ml_model = av1_simple_motion_search_term_none_model_16;
  } else {
    assert(0 && "Unexpected block size in simple_motion_term_none");
  }

  if (ml_model) {
    float score = 0.0f;
    for (f_idx = 0; f_idx < NUM_FEATURES; f_idx++) {
      score += ml_model[f_idx] *
               (simple_motion_features[f_idx] - ml_mean[f_idx]) / ml_std[f_idx];
    }
    score += ml_model[NUM_FEATURES];

    if (score >= 0.0f) {
      *early_terminate = 1;
    }
  }
}
#undef NUM_FEATURES

static void firstpass_simple_motion_search_features(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, float *features) {
  assert(mi_size_wide[bsize] == mi_size_high[bsize]);
  assert(cpi->ref_frame_flags & av1_ref_frame_flag_list[LAST_FRAME] ||
         cpi->ref_frame_flags & av1_ref_frame_flag_list[ALTREF_FRAME]);

  // Setting up motion search
  const int ref_list[] = { LAST_FRAME, ALTREF_FRAME };
  const int num_refs = 2;
  const int use_subpixel = 0;

  unsigned int int_features[10] = { 0 };

  int f_idx = 0;
  // Doing whole block first to update the mv
  simple_motion_search_get_best_ref(
      cpi, x, pc_tree, mi_row, mi_col, bsize, ref_list, num_refs, use_subpixel,
      4, &int_features[f_idx], &int_features[f_idx + 1]);
  f_idx += 2;

  // Split subblocks
  const BLOCK_SIZE subsize = get_partition_subsize(bsize, PARTITION_SPLIT);
  const int w_mi = mi_size_wide[bsize];
  const int h_mi = mi_size_high[bsize];
  for (int r_idx = 0; r_idx < 4; r_idx++) {
    const int sub_mi_col = mi_col + (r_idx & 1) * w_mi / 2;
    const int sub_mi_row = mi_row + (r_idx >> 1) * h_mi / 2;

    simple_motion_search_get_best_ref(
        cpi, x, pc_tree, sub_mi_row, sub_mi_col, subsize, ref_list, num_refs,
        use_subpixel, r_idx, &int_features[f_idx], &int_features[f_idx + 1]);
    f_idx += 2;
  }

  aom_clear_system_state();
  for (int idx = 0; idx < f_idx; idx++) {
    features[idx] = logf(1.0f + (float)int_features[idx]);
  }

  const MACROBLOCKD *xd = &x->e_mbd;
  set_offsets_for_motion_search(cpi, x, mi_row, mi_col, bsize);

  // Q_INDEX
  const int dc_q = av1_dc_quant_QTX(x->qindex, 0, xd->bd) >> (xd->bd - 8);
  features[f_idx++] = logf(1.0f + (float)(dc_q * dc_q) / 256.0f);

  // Neighbor stuff
  const int has_above = !!xd->above_mbmi;
  const int has_left = !!xd->left_mbmi;
  const BLOCK_SIZE above_bsize = has_above ? xd->above_mbmi->sb_type : bsize;
  const BLOCK_SIZE left_bsize = has_left ? xd->left_mbmi->sb_type : bsize;
  features[f_idx++] = (float)has_above;
  features[f_idx++] = (float)mi_size_wide_log2[above_bsize];
  features[f_idx++] = (float)mi_size_high_log2[above_bsize];
  features[f_idx++] = (float)has_left;
  features[f_idx++] = (float)mi_size_wide_log2[left_bsize];
  features[f_idx++] = (float)mi_size_high_log2[left_bsize];
}

#define NUM_FEATURES 20
void av1_firstpass_simple_motion_search_early_term(AV1_COMP *const cpi,
                                                   MACROBLOCK *x,
                                                   PC_TREE *pc_tree, int mi_row,
                                                   int mi_col, BLOCK_SIZE bsize,
                                                   const RD_STATS *none_rdc,
                                                   int *do_square_split) {
  const NN_CONFIG *nn_config = NULL;
  float thresh = 0.0f;
  const float *ml_mean = NULL, *ml_std = NULL;
  if (bsize == BLOCK_32X32) {
    nn_config = &av1_fp_simple_motion_search_term_none_nn_config_32;
    ml_mean = av1_fp_simple_motion_search_term_none_mean_32;
    ml_std = av1_fp_simple_motion_search_term_none_std_32;
    thresh = av1_fp_simple_motion_search_term_none_thresh_32;
  } else if (bsize == BLOCK_16X16) {
    nn_config = &av1_fp_simple_motion_search_term_none_nn_config_16;
    ml_mean = av1_fp_simple_motion_search_term_none_mean_16;
    ml_std = av1_fp_simple_motion_search_term_none_std_16;
    thresh = av1_fp_simple_motion_search_term_none_thresh_16;
  } else if (bsize == BLOCK_8X8) {
    nn_config = &av1_fp_simple_motion_search_term_none_nn_config_8;
    ml_mean = av1_fp_simple_motion_search_term_none_mean_8;
    ml_std = av1_fp_simple_motion_search_term_none_std_8;
    thresh = av1_fp_simple_motion_search_term_none_thresh_8;
  } else {
    assert(0 &&
           "Unexpected bsize in firstpass_simple_motion_search_early_term");
    return;
  }

  float ml_features[NUM_FEATURES] = { 0.0f };

  firstpass_simple_motion_search_features(cpi, x, pc_tree, mi_row, mi_col,
                                          bsize, ml_features);
  int f_idx = 17;

  ml_features[f_idx++] = logf(1.0f + (float)none_rdc->rate);
  ml_features[f_idx++] = logf(1.0f + (float)none_rdc->dist);
  ml_features[f_idx++] = logf(1.0f + (float)none_rdc->rdcost);

  for (f_idx = 0; f_idx < 20; f_idx++) {
    ml_features[f_idx] = (ml_features[f_idx] - ml_mean[f_idx]) / ml_std[f_idx];
  }

  // Get probabilities
  float score = 0.0f;

  av1_nn_predict(ml_features, nn_config, &score);
  aom_clear_system_state();

  // Determine if we should prune square partitions.
  if (score < thresh) {
    *do_square_split = 0;
  }
}
#undef NUM_FEATURES
