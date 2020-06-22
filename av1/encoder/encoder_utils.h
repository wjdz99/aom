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

#ifndef AOM_AV1_ENCODER_ENCODER_UTILS_H_
#define AOM_AV1_ENCODER_ENCODER_UTILS_H_

#include "config/aom_scale_rtcd.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define AM_SEGMENT_ID_INACTIVE 7
#define AM_SEGMENT_ID_ACTIVE 0

extern const int default_tx_type_probs[FRAME_UPDATE_TYPES][TX_SIZES_ALL]
                                      [TX_TYPES];

extern const int default_obmc_probs[FRAME_UPDATE_TYPES][BLOCK_SIZES_ALL];

extern const int default_warped_probs[FRAME_UPDATE_TYPES];

// TODO(yunqing): the default probs can be trained later from better
// performance.
extern const int default_switchable_interp_probs[FRAME_UPDATE_TYPES]
                                                [SWITCHABLE_FILTER_CONTEXTS]
                                                [SWITCHABLE_FILTERS];

// Mark all inactive blocks as active. Other segmentation features may be set
// so memset cannot be used, instead only inactive blocks should be reset.
static AOM_INLINE void suppress_active_map(AV1_COMP *cpi) {
  unsigned char *const seg_map = cpi->enc_seg.map;
  int i;
  if (cpi->active_map.enabled || cpi->active_map.update)
    for (i = 0;
         i < cpi->common.mi_params.mi_rows * cpi->common.mi_params.mi_cols; ++i)
      if (seg_map[i] == AM_SEGMENT_ID_INACTIVE)
        seg_map[i] = AM_SEGMENT_ID_ACTIVE;
}

static AOM_INLINE void set_mb_mi(CommonModeInfoParams *mi_params, int width,
                                 int height) {
  // Ensure that the decoded width and height are both multiples of
  // 8 luma pixels (note: this may only be a multiple of 4 chroma pixels if
  // subsampling is used).
  // This simplifies the implementation of various experiments,
  // eg. cdef, which operates on units of 8x8 luma pixels.
  const int aligned_width = ALIGN_POWER_OF_TWO(width, 3);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, 3);

  mi_params->mi_cols = aligned_width >> MI_SIZE_LOG2;
  mi_params->mi_rows = aligned_height >> MI_SIZE_LOG2;
  mi_params->mi_stride = calc_mi_size(mi_params->mi_cols);

  mi_params->mb_cols = (mi_params->mi_cols + 2) >> 2;
  mi_params->mb_rows = (mi_params->mi_rows + 2) >> 2;
  mi_params->MBs = mi_params->mb_rows * mi_params->mb_cols;

  const int mi_alloc_size_1d = mi_size_wide[mi_params->mi_alloc_bsize];
  mi_params->mi_alloc_stride =
      (mi_params->mi_stride + mi_alloc_size_1d - 1) / mi_alloc_size_1d;

  assert(mi_size_wide[mi_params->mi_alloc_bsize] ==
         mi_size_high[mi_params->mi_alloc_bsize]);

#if CONFIG_LPF_MASK
  av1_alloc_loop_filter_mask(mi_params);
#endif
}

static AOM_INLINE void enc_free_mi(CommonModeInfoParams *mi_params) {
  aom_free(mi_params->mi_alloc);
  mi_params->mi_alloc = NULL;
  aom_free(mi_params->mi_grid_base);
  mi_params->mi_grid_base = NULL;
  mi_params->mi_alloc_size = 0;
  aom_free(mi_params->tx_type_map);
  mi_params->tx_type_map = NULL;
}

static AOM_INLINE void enc_set_mb_mi(CommonModeInfoParams *mi_params, int width,
                                     int height) {
  const int is_4k_or_larger = AOMMIN(width, height) >= 2160;
  mi_params->mi_alloc_bsize = is_4k_or_larger ? BLOCK_8X8 : BLOCK_4X4;

  set_mb_mi(mi_params, width, height);
}

static AOM_INLINE void stat_stage_set_mb_mi(CommonModeInfoParams *mi_params,
                                            int width, int height) {
  mi_params->mi_alloc_bsize = BLOCK_16X16;

  set_mb_mi(mi_params, width, height);
}

static AOM_INLINE void enc_setup_mi(CommonModeInfoParams *mi_params) {
  const int mi_grid_size =
      mi_params->mi_stride * calc_mi_size(mi_params->mi_rows);
  memset(mi_params->mi_alloc, 0,
         mi_params->mi_alloc_size * sizeof(*mi_params->mi_alloc));
  memset(mi_params->mi_grid_base, 0,
         mi_grid_size * sizeof(*mi_params->mi_grid_base));
  memset(mi_params->tx_type_map, 0,
         mi_grid_size * sizeof(*mi_params->tx_type_map));
}

static AOM_INLINE void init_buffer_indices(
    ForceIntegerMVInfo *const force_intpel_info, int *const remapped_ref_idx) {
  int fb_idx;
  for (fb_idx = 0; fb_idx < REF_FRAMES; ++fb_idx)
    remapped_ref_idx[fb_idx] = fb_idx;
  force_intpel_info->rate_index = 0;
  force_intpel_info->rate_size = 0;
}

static AOM_INLINE void copy_frame_prob_info(AV1_COMP *cpi) {
  FrameProbInfo *const frame_probs = &cpi->frame_probs;
  if (cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats) {
    av1_copy(frame_probs->tx_type_probs, default_tx_type_probs);
  }
  if (!cpi->sf.inter_sf.disable_obmc &&
      cpi->sf.inter_sf.prune_obmc_prob_thresh > 0) {
    av1_copy(frame_probs->obmc_probs, default_obmc_probs);
  }
  if (cpi->sf.inter_sf.prune_warped_prob_thresh > 0) {
    av1_copy(frame_probs->warped_probs, default_warped_probs);
  }
  if (cpi->sf.interp_sf.adaptive_interp_filter_search == 2) {
    av1_copy(frame_probs->switchable_interp_probs,
             default_switchable_interp_probs);
  }
}

static AOM_INLINE void restore_cur_buf(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;
  aom_yv12_copy_frame(&cc->copy_buffer, &cm->cur_frame->buf,
                      av1_num_planes(cm));
}

// Coding context that only needs to be restored when recode loop includes
// filtering (deblocking, CDEF, superres post-encode upscale and/or loop
// restoraton).
static AOM_INLINE void restore_extra_coding_context(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;
  cm->lf = cc->lf;
  cm->cdef_info = cc->cdef_info;
  cpi->rc = cc->rc;
}

static AOM_INLINE void release_copy_buffer(CODING_CONTEXT *cc) {
  aom_free_frame_buffer(&cc->copy_buffer);
}

static AOM_INLINE int equal_dimensions_and_border(const YV12_BUFFER_CONFIG *a,
                                                  const YV12_BUFFER_CONFIG *b) {
  return a->y_height == b->y_height && a->y_width == b->y_width &&
         a->uv_height == b->uv_height && a->uv_width == b->uv_width &&
         a->y_stride == b->y_stride && a->uv_stride == b->uv_stride &&
         a->border == b->border &&
         (a->flags & YV12_FLAG_HIGHBITDEPTH) ==
             (b->flags & YV12_FLAG_HIGHBITDEPTH);
}

static AOM_INLINE int update_entropy(bool *ext_refresh_frame_context,
                                     bool *ext_refresh_frame_context_pending,
                                     bool update) {
  *ext_refresh_frame_context = update;
  *ext_refresh_frame_context_pending = 1;
  return 0;
}

#if !CONFIG_REALTIME_ONLY
static AOM_INLINE int combine_prior_with_tpl_boost(double min_factor,
                                                   double max_factor,
                                                   int prior_boost,
                                                   int tpl_boost,
                                                   int frames_to_key) {
  double factor = sqrt((double)frames_to_key);
  double range = max_factor - min_factor;
  factor = AOMMIN(factor, max_factor);
  factor = AOMMAX(factor, min_factor);
  factor -= min_factor;
  int boost =
      (int)((factor * prior_boost + (range - factor) * tpl_boost) / range);
  return boost;
}
#endif

#if CONFIG_AV1_HIGHBITDEPTH
void highbd_set_var_fns(AV1_COMP *const cpi);
#endif

#if !CONFIG_REALTIME_ONLY
void process_tpl_stats_frame(AV1_COMP *cpi);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_ENCODER_UTILS_H_
