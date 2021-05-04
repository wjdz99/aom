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

#ifndef AOM_AV1_ENCODER_RECONINTER_ENC_H_
#define AOM_AV1_ENCODER_RECONINTER_ENC_H_

#include "aom/aom_integer.h"
#include "av1/common/av1_common_int.h"
#include "av1/common/blockd.h"
#include "av1/common/convolve.h"
#include "av1/common/filter.h"
#include "av1/common/reconinter.h"
#include "av1/common/warped_motion.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD   (de3a75 Set mv_step_param based on previous frame in speed >= 2)
void aom_highbd_comp_mask_upsampled_pred(
    MACROBLOCKD *xd, const struct AV1Common *const cm, int mi_row, int mi_col,
    const MV *const mv, uint8_t *comp_pred8, const uint8_t *pred8, int width,
    int height, int subpel_x_q3, int subpel_y_q3, const uint8_t *ref8,
    int ref_stride, const uint8_t *mask, int mask_stride, int invert_mask,
    int bd, int subpel_search);
=======
static INLINE void enc_calc_subpel_params(const MV *const src_mv,
                                   InterPredParams *const inter_pred_params,
                                   MACROBLOCKD *xd, int mi_x, int mi_y, int ref,
                                   uint8_t **mc_buf, uint8_t **pre,
                                   SubpelParams *subpel_params,
                                   int *src_stride) {
  // These are part of the function signature to use this function through a
  // function pointer. See typedef of 'CalcSubpelParamsFunc'.
  (void)xd;
  (void)mi_x;
  (void)mi_y;
  (void)ref;
  (void)mc_buf;

  const struct scale_factors *sf = inter_pred_params->scale_factors;

  struct buf_2d *pre_buf = &inter_pred_params->ref_frame_buf;
  int ssx = inter_pred_params->subsampling_x;
  int ssy = inter_pred_params->subsampling_y;
  int orig_pos_y = inter_pred_params->pix_row << SUBPEL_BITS;
  orig_pos_y += src_mv->row * (1 << (1 - ssy));
  int orig_pos_x = inter_pred_params->pix_col << SUBPEL_BITS;
  orig_pos_x += src_mv->col * (1 << (1 - ssx));
  int pos_y = sf->scale_value_y(orig_pos_y, sf);
  int pos_x = sf->scale_value_x(orig_pos_x, sf);
  pos_x += SCALE_EXTRA_OFF;
  pos_y += SCALE_EXTRA_OFF;

  const int top = -AOM_LEFT_TOP_MARGIN_SCALED(ssy);
  const int left = -AOM_LEFT_TOP_MARGIN_SCALED(ssx);
  const int bottom = (pre_buf->height + AOM_INTERP_EXTEND) << SCALE_SUBPEL_BITS;
  const int right = (pre_buf->width + AOM_INTERP_EXTEND) << SCALE_SUBPEL_BITS;
  pos_y = clamp(pos_y, top, bottom);
  pos_x = clamp(pos_x, left, right);

  subpel_params->subpel_x = pos_x & SCALE_SUBPEL_MASK;
  subpel_params->subpel_y = pos_y & SCALE_SUBPEL_MASK;
  subpel_params->xs = sf->x_step_q4;
  subpel_params->ys = sf->y_step_q4;
  *pre = pre_buf->buf0 + (pos_y >> SCALE_SUBPEL_BITS) * pre_buf->stride +
         (pos_x >> SCALE_SUBPEL_BITS);
  *src_stride = pre_buf->stride;
}


>>>>>>> CHANGE (4a2a54 wip: rot in tf)

// Build single or compound reference inter predictors for all planes.
// Can build inter-intra predictors, masked predictors etc as well.
void av1_enc_build_inter_predictor(const AV1_COMMON *cm, MACROBLOCKD *xd,
                                   int mi_row, int mi_col,
                                   const BUFFER_SET *ctx, BLOCK_SIZE bsize,
                                   int plane_from, int plane_to);

void av1_enc_build_inter_predictor_y(MACROBLOCKD *xd, int mi_row, int mi_col);

// Build one inter predictor. It is called for building predictor for single
// reference case, or just the 1st or 2nd reference in compound reference case.
// Can build both regular and masked predictors.
void av1_enc_build_one_inter_predictor(uint8_t *dst, int dst_stride,
                                       const MV *src_mv,
                                       InterPredParams *inter_pred_params);

void av1_build_prediction_by_above_preds(const AV1_COMMON *cm, MACROBLOCKD *xd,
                                         uint8_t *tmp_buf[MAX_MB_PLANE],
                                         int tmp_width[MAX_MB_PLANE],
                                         int tmp_height[MAX_MB_PLANE],
                                         int tmp_stride[MAX_MB_PLANE]);

void av1_build_prediction_by_left_preds(const AV1_COMMON *cm, MACROBLOCKD *xd,
                                        uint8_t *tmp_buf[MAX_MB_PLANE],
                                        int tmp_width[MAX_MB_PLANE],
                                        int tmp_height[MAX_MB_PLANE],
                                        int tmp_stride[MAX_MB_PLANE]);

void av1_build_obmc_inter_predictors_sb(const AV1_COMMON *cm, MACROBLOCKD *xd);

void av1_build_inter_predictors_for_planes_single_buf(
    MACROBLOCKD *xd, BLOCK_SIZE bsize, int plane_from, int plane_to, int ref,
    uint8_t *ext_dst[3], int ext_dst_stride[3]);

void av1_build_wedge_inter_predictor_from_buf(MACROBLOCKD *xd, BLOCK_SIZE bsize,
                                              int plane_from, int plane_to,
                                              uint8_t *ext_dst0[3],
                                              int ext_dst_stride0[3],
                                              uint8_t *ext_dst1[3],
                                              int ext_dst_stride1[3]);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_RECONINTER_ENC_H_
