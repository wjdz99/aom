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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "config/aom_dsp_rtcd.h"

#include "aom_dsp/binary_codes_writer.h"
#include "aom_ports/system_state.h"

#include "av1/common/convolve.h"
#include "av1/common/resize.h"
#include "av1/common/warped_motion.h"

#include "av1/encoder/corner_detect.h"
#include "av1/encoder/corner_match.h"
#include "av1/encoder/global_motion.h"
#include "av1/encoder/ransac.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/segmentation.h"

#define MIN_INLIER_PROB 0.1

#define MIN_TRANS_THRESH (1 * GM_TRANS_DECODE_FACTOR)

// Border over which to compute the global motion
#define ERRORADV_BORDER 0

// Number of pyramid levels in disflow computation
#define N_LEVELS 2
// Size of square patches in the disflow dense grid
#define PATCH_SIZE 8
// Center point of square patch
#define PATCH_CENTER ((PATCH_SIZE + 1) >> 1)
// Step size between patches, lower value means greater patch overlap
#define PATCH_STEP 1
// Minimum size of border padding for disflow
#define MIN_PAD 7
// Warp error convergence threshold for disflow
#define DISFLOW_ERROR_TR 0.01
// Max number of iterations if warp convergence is not found
#define DISFLOW_MAX_ITR 10
// Highest motion model to search
#define GLOBAL_TRANS_TYPES_ENC 3

// Struct for an image pyramid
typedef struct {
  int n_levels;
  int pad_size;
  int has_gradient;
  int widths[N_LEVELS];
  int heights[N_LEVELS];
  int strides[N_LEVELS];
  int level_loc[N_LEVELS];
  unsigned char *level_buffer;
  double *level_dx_buffer;
  double *level_dy_buffer;
} ImagePyramid;

int av1_is_enough_erroradvantage(double best_erroradvantage, int params_cost,
                                 int erroradv_type) {
  assert(erroradv_type < GM_ERRORADV_TR_TYPES);
  return best_erroradvantage < erroradv_tr[erroradv_type] &&
         best_erroradvantage * params_cost < erroradv_prod_tr[erroradv_type];
}

static void convert_to_params(const double *params, int32_t *model) {
  int i;
  int alpha_present = 0;
  model[0] = (int32_t)floor(params[0] * (1 << GM_TRANS_PREC_BITS) + 0.5);
  model[1] = (int32_t)floor(params[1] * (1 << GM_TRANS_PREC_BITS) + 0.5);
  model[0] = (int32_t)clamp(model[0], GM_TRANS_MIN, GM_TRANS_MAX) *
             GM_TRANS_DECODE_FACTOR;
  model[1] = (int32_t)clamp(model[1], GM_TRANS_MIN, GM_TRANS_MAX) *
             GM_TRANS_DECODE_FACTOR;

  for (i = 2; i < 6; ++i) {
    const int diag_value = ((i == 2 || i == 5) ? (1 << GM_ALPHA_PREC_BITS) : 0);
    model[i] = (int32_t)floor(params[i] * (1 << GM_ALPHA_PREC_BITS) + 0.5);
    model[i] =
        (int32_t)clamp(model[i] - diag_value, GM_ALPHA_MIN, GM_ALPHA_MAX);
    alpha_present |= (model[i] != 0);
    model[i] = (model[i] + diag_value) * GM_ALPHA_DECODE_FACTOR;
  }
  for (; i < 8; ++i) {
    model[i] = (int32_t)floor(params[i] * (1 << GM_ROW3HOMO_PREC_BITS) + 0.5);
    model[i] = (int32_t)clamp(model[i], GM_ROW3HOMO_MIN, GM_ROW3HOMO_MAX) *
               GM_ROW3HOMO_DECODE_FACTOR;
    alpha_present |= (model[i] != 0);
  }

  if (!alpha_present) {
    if (abs(model[0]) < MIN_TRANS_THRESH && abs(model[1]) < MIN_TRANS_THRESH) {
      model[0] = 0;
      model[1] = 0;
    }
  }
}

void av1_convert_model_to_params(const double *params,
                                 WarpedMotionParams *model) {
  convert_to_params(params, model->wmmat);
  model->wmtype = get_wmtype(model);
  model->invalid = 0;
}

// Adds some offset to a global motion parameter and handles
// all of the necessary precision shifts, clamping, and
// zero-centering.
static int32_t add_param_offset(int param_index, int32_t param_value,
                                int32_t offset) {
  const int scale_vals[3] = { GM_TRANS_PREC_DIFF, GM_ALPHA_PREC_DIFF,
                              GM_ROW3HOMO_PREC_DIFF };
  const int clamp_vals[3] = { GM_TRANS_MAX, GM_ALPHA_MAX, GM_ROW3HOMO_MAX };
  // type of param: 0 - translation, 1 - affine, 2 - homography
  const int param_type = (param_index < 2 ? 0 : (param_index < 6 ? 1 : 2));
  const int is_one_centered = (param_index == 2 || param_index == 5);

  // Make parameter zero-centered and offset the shift that was done to make
  // it compatible with the warped model
  param_value = (param_value - (is_one_centered << WARPEDMODEL_PREC_BITS)) >>
                scale_vals[param_type];
  // Add desired offset to the rescaled/zero-centered parameter
  param_value += offset;
  // Clamp the parameter so it does not overflow the number of bits allotted
  // to it in the bitstream
  param_value = (int32_t)clamp(param_value, -clamp_vals[param_type],
                               clamp_vals[param_type]);
  // Rescale the parameter to WARPEDMODEL_PRECISION_BITS so it is compatible
  // with the warped motion library
  param_value *= (1 << scale_vals[param_type]);

  // Undo the zero-centering step if necessary
  return param_value + (is_one_centered << WARPEDMODEL_PREC_BITS);
}

static void force_wmtype(WarpedMotionParams *wm, TransformationType wmtype) {
  switch (wmtype) {
    case IDENTITY:
      wm->wmmat[0] = 0;
      wm->wmmat[1] = 0;
      AOM_FALLTHROUGH_INTENDED;
    case TRANSLATION:
      wm->wmmat[2] = 1 << WARPEDMODEL_PREC_BITS;
      wm->wmmat[3] = 0;
      AOM_FALLTHROUGH_INTENDED;
    case ROTZOOM:
      wm->wmmat[4] = -wm->wmmat[3];
      wm->wmmat[5] = wm->wmmat[2];
      AOM_FALLTHROUGH_INTENDED;
    case AFFINE: wm->wmmat[6] = wm->wmmat[7] = 0; break;
    default: assert(0);
  }
  wm->wmtype = wmtype;
}

#if CONFIG_AV1_HIGHBITDEPTH
static int64_t highbd_warp_error(
    WarpedMotionParams *wm, const uint16_t *const ref, int width, int height,
    int stride, const uint16_t *const dst, int p_col, int p_row, int p_width,
    int p_height, int p_stride, int subsampling_x, int subsampling_y, int bd,
    int64_t best_error, uint8_t *segment_map, int segment_map_stride) {
  int64_t gm_sumerr = 0;
  const int error_bsize_w = AOMMIN(p_width, WARP_ERROR_BLOCK);
  const int error_bsize_h = AOMMIN(p_height, WARP_ERROR_BLOCK);
  uint16_t tmp[WARP_ERROR_BLOCK * WARP_ERROR_BLOCK];

  ConvolveParams conv_params = get_conv_params(0, 0, bd);
  conv_params.use_dist_wtd_comp_avg = 0;
  for (int i = p_row; i < p_row + p_height; i += WARP_ERROR_BLOCK) {
    for (int j = p_col; j < p_col + p_width; j += WARP_ERROR_BLOCK) {
      int seg_x = j >> WARP_ERROR_BLOCK_LOG;
      int seg_y = i >> WARP_ERROR_BLOCK_LOG;
      // Only compute the error if this block contains inliers from the motion
      // model
      if (!segment_map[seg_y * segment_map_stride + seg_x]) continue;
      // avoid warping extra 8x8 blocks in the padded region of the frame
      // when p_width and p_height are not multiples of WARP_ERROR_BLOCK
      const int warp_w = AOMMIN(error_bsize_w, p_col + p_width - j);
      const int warp_h = AOMMIN(error_bsize_h, p_row + p_height - i);
      highbd_warp_plane(wm, ref, width, height, stride, tmp, j, i, warp_w,
                        warp_h, WARP_ERROR_BLOCK, subsampling_x, subsampling_y,
                        bd, &conv_params);
      gm_sumerr += av1_calc_highbd_frame_error(tmp, WARP_ERROR_BLOCK,
                                               dst + j + i * p_stride, warp_w,
                                               warp_h, p_stride, bd);
      if (gm_sumerr > best_error) return INT64_MAX;
    }
  }
  return gm_sumerr;
}
#endif

static int64_t warp_error(WarpedMotionParams *wm, const uint8_t *const ref,
                          int width, int height, int stride,
                          const uint8_t *const dst, int p_col, int p_row,
                          int p_width, int p_height, int p_stride,
                          int subsampling_x, int subsampling_y,
                          int64_t best_error, uint8_t *segment_map,
                          int segment_map_stride) {
  int64_t gm_sumerr = 0;
  int warp_w, warp_h;
  const int error_bsize_w = AOMMIN(p_width, WARP_ERROR_BLOCK);
  const int error_bsize_h = AOMMIN(p_height, WARP_ERROR_BLOCK);
  DECLARE_ALIGNED(16, uint8_t, tmp[WARP_ERROR_BLOCK * WARP_ERROR_BLOCK]);
  ConvolveParams conv_params = get_conv_params(0, 0, 8);
  conv_params.use_dist_wtd_comp_avg = 0;

  for (int i = p_row; i < p_row + p_height; i += WARP_ERROR_BLOCK) {
    for (int j = p_col; j < p_col + p_width; j += WARP_ERROR_BLOCK) {
      int seg_x = j >> WARP_ERROR_BLOCK_LOG;
      int seg_y = i >> WARP_ERROR_BLOCK_LOG;
      // Only compute the error if this block contains inliers from the motion
      // model
      if (!segment_map[seg_y * segment_map_stride + seg_x]) continue;
      // avoid warping extra 8x8 blocks in the padded region of the frame
      // when p_width and p_height are not multiples of WARP_ERROR_BLOCK
      warp_w = AOMMIN(error_bsize_w, p_col + p_width - j);
      warp_h = AOMMIN(error_bsize_h, p_row + p_height - i);
      warp_plane(wm, ref, width, height, stride, tmp, j, i, warp_w, warp_h,
                 WARP_ERROR_BLOCK, subsampling_x, subsampling_y, &conv_params);

      gm_sumerr +=
          av1_calc_frame_error(tmp, WARP_ERROR_BLOCK, dst + j + i * p_stride,
                               warp_w, warp_h, p_stride);
      if (gm_sumerr > best_error) return INT64_MAX;
    }
  }
  return gm_sumerr;
}

int64_t av1_warp_error(WarpedMotionParams *wm, int use_hbd, int bd,
                       const uint8_t *ref, int width, int height, int stride,
                       uint8_t *dst, int p_col, int p_row, int p_width,
                       int p_height, int p_stride, int subsampling_x,
                       int subsampling_y, int64_t best_error,
                       uint8_t *segment_map, int segment_map_stride) {
  if (wm->wmtype <= AFFINE)
    if (!av1_get_shear_params(wm)) return INT64_MAX;
#if CONFIG_AV1_HIGHBITDEPTH
  if (use_hbd)
    return highbd_warp_error(wm, CONVERT_TO_SHORTPTR(ref), width, height,
                             stride, CONVERT_TO_SHORTPTR(dst), p_col, p_row,
                             p_width, p_height, p_stride, subsampling_x,
                             subsampling_y, bd, best_error, segment_map,
                             segment_map_stride);
#endif
  (void)use_hbd;
  (void)bd;
  return warp_error(wm, ref, width, height, stride, dst, p_col, p_row, p_width,
                    p_height, p_stride, subsampling_x, subsampling_y,
                    best_error, segment_map, segment_map_stride);
}

// Factors used to calculate the thresholds for av1_warp_error
static double thresh_factors[GM_REFINEMENT_COUNT] = { 1.25, 1.20, 1.15, 1.10,
                                                      1.05 };

static INLINE int64_t calc_approx_erroradv_threshold(
    double scaling_factor, int64_t erroradv_threshold) {
  return erroradv_threshold <
                 (int64_t)(((double)INT64_MAX / scaling_factor) + 0.5)
             ? (int64_t)(scaling_factor * erroradv_threshold + 0.5)
             : INT64_MAX;
}

int64_t av1_refine_integerized_param(
    WarpedMotionParams *wm, TransformationType wmtype, int use_hbd, int bd,
    uint8_t *ref, int r_width, int r_height, int r_stride, uint8_t *dst,
    int d_width, int d_height, int d_stride, int n_refinements,
    int64_t best_frame_error, uint8_t *segment_map, int segment_map_stride,
    int64_t erroradv_threshold) {
  static const int max_trans_model_params[TRANS_TYPES] = { 0, 2, 4, 6 };
  const int border = ERRORADV_BORDER;
  int i = 0, p;
  int n_params = max_trans_model_params[wmtype];
  int32_t *param_mat = wm->wmmat;
  int64_t step_error, best_error;
  int32_t step;
  int32_t *param;
  int32_t curr_param;
  int32_t best_param;

  force_wmtype(wm, wmtype);
  best_error =
      av1_warp_error(wm, use_hbd, bd, ref, r_width, r_height, r_stride,
                     dst + border * d_stride + border, border, border,
                     d_width - 2 * border, d_height - 2 * border, d_stride, 0,
                     0, best_frame_error, segment_map, segment_map_stride);
  best_error = AOMMIN(best_error, best_frame_error);
  step = 1 << (n_refinements - 1);
  for (i = 0; i < n_refinements; i++, step >>= 1) {
    int64_t error_adv_thresh =
        calc_approx_erroradv_threshold(thresh_factors[i], erroradv_threshold);
    for (p = 0; p < n_params; ++p) {
      int step_dir = 0;
      // Skip searches for parameters that are forced to be 0
      param = param_mat + p;
      curr_param = *param;
      best_param = curr_param;
      // look to the left
      *param = add_param_offset(p, curr_param, -step);
      step_error =
          av1_warp_error(wm, use_hbd, bd, ref, r_width, r_height, r_stride,
                         dst + border * d_stride + border, border, border,
                         d_width - 2 * border, d_height - 2 * border, d_stride,
                         0, 0, AOMMIN(best_error, error_adv_thresh),
                         segment_map, segment_map_stride);
      if (step_error < best_error) {
        best_error = step_error;
        best_param = *param;
        step_dir = -1;
      }

      // look to the right
      *param = add_param_offset(p, curr_param, step);
      step_error =
          av1_warp_error(wm, use_hbd, bd, ref, r_width, r_height, r_stride,
                         dst + border * d_stride + border, border, border,
                         d_width - 2 * border, d_height - 2 * border, d_stride,
                         0, 0, AOMMIN(best_error, error_adv_thresh),
                         segment_map, segment_map_stride);
      if (step_error < best_error) {
        best_error = step_error;
        best_param = *param;
        step_dir = 1;
      }
      *param = best_param;

      // look to the direction chosen above repeatedly until error increases
      // for the biggest step size
      while (step_dir) {
        *param = add_param_offset(p, best_param, step * step_dir);
        step_error =
            av1_warp_error(wm, use_hbd, bd, ref, r_width, r_height, r_stride,
                           dst + border * d_stride + border, border, border,
                           d_width - 2 * border, d_height - 2 * border,
                           d_stride, 0, 0, AOMMIN(best_error, error_adv_thresh),
                           segment_map, segment_map_stride);
        if (step_error < best_error) {
          best_error = step_error;
          best_param = *param;
        } else {
          *param = best_param;
          step_dir = 0;
        }
      }
    }
  }
  force_wmtype(wm, wmtype);
  wm->wmtype = get_wmtype(wm);
  return best_error;
}

unsigned char *av1_downconvert_frame(YV12_BUFFER_CONFIG *frm, int bit_depth) {
  int i, j;
  uint16_t *orig_buf = CONVERT_TO_SHORTPTR(frm->y_buffer);
  uint8_t *buf_8bit = frm->y_buffer_8bit;
  assert(buf_8bit);
  if (!frm->buf_8bit_valid) {
    for (i = 0; i < frm->y_height; ++i) {
      for (j = 0; j < frm->y_width; ++j) {
        buf_8bit[i * frm->y_stride + j] =
            orig_buf[i * frm->y_stride + j] >> (bit_depth - 8);
      }
    }
    frm->buf_8bit_valid = 1;
  }
  return buf_8bit;
}

static void get_inliers_from_indices(MotionModel *params,
                                     int *correspondences) {
  int *inliers_tmp = (int *)aom_malloc(2 * MAX_CORNERS * sizeof(*inliers_tmp));
  memset(inliers_tmp, 0, 2 * MAX_CORNERS * sizeof(*inliers_tmp));

  for (int i = 0; i < params->num_inliers; i++) {
    int index = params->inliers[i];
    inliers_tmp[2 * i] = correspondences[4 * index];
    inliers_tmp[2 * i + 1] = correspondences[4 * index + 1];
  }
  memcpy(params->inliers, inliers_tmp, sizeof(*inliers_tmp) * 2 * MAX_CORNERS);
  aom_free(inliers_tmp);
}

#define FEAT_COUNT_TR 3
#define SEG_COUNT_TR 0.40
void av1_compute_feature_segmentation_map(uint8_t *segment_map, int width,
                                          int height, int *inliers,
                                          int num_inliers) {
  int seg_count = 0;
  memset(segment_map, 0, sizeof(*segment_map) * width * height);

  for (int i = 0; i < num_inliers; i++) {
    int x = inliers[i * 2];
    int y = inliers[i * 2 + 1];
    int seg_x = x >> WARP_ERROR_BLOCK_LOG;
    int seg_y = y >> WARP_ERROR_BLOCK_LOG;
    segment_map[seg_y * width + seg_x] += 1;
  }

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      uint8_t feat_count = segment_map[i * width + j];
      segment_map[i * width + j] = (feat_count >= FEAT_COUNT_TR);
      seg_count += (segment_map[i * width + j]);
    }
  }

  // If this motion does not make up a large enough portion of the frame,
  // use the unsegmented version of the error metric
  if (seg_count < (width * height * SEG_COUNT_TR))
    memset(segment_map, 1, width * height * sizeof(*segment_map));
}

static int compute_global_motion_feature_based(
    TransformationType type, unsigned char *src_buffer, int src_width,
    int src_height, int src_stride, int *src_corners, int num_src_corners,
    YV12_BUFFER_CONFIG *ref, int bit_depth, int *num_inliers_by_motion,
    MotionModel *params_by_motion, int num_motions) {
  int i;
  int num_ref_corners;
  int num_correspondences;
  int *correspondences;
  int ref_corners[2 * MAX_CORNERS];
  unsigned char *ref_buffer = ref->y_buffer;
  RansacFunc ransac = av1_get_ransac_type(type);

  if (ref->flags & YV12_FLAG_HIGHBITDEPTH) {
    ref_buffer = av1_downconvert_frame(ref, bit_depth);
  }

  num_ref_corners =
      av1_fast_corner_detect(ref_buffer, ref->y_width, ref->y_height,
                             ref->y_stride, ref_corners, MAX_CORNERS);

  // find correspondences between the two images
  correspondences =
      (int *)malloc(num_src_corners * 4 * sizeof(*correspondences));
  num_correspondences = av1_determine_correspondence(
      src_buffer, (int *)src_corners, num_src_corners, ref_buffer,
      (int *)ref_corners, num_ref_corners, src_width, src_height, src_stride,
      ref->y_stride, correspondences);

  ransac(correspondences, num_correspondences, num_inliers_by_motion,
         params_by_motion, num_motions);

  // Set num_inliers = 0 for motions with too few inliers so they are ignored.
  for (i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] < MIN_INLIER_PROB * num_correspondences ||
        num_correspondences == 0) {
      num_inliers_by_motion[i] = 0;
    } else {
      get_inliers_from_indices(&params_by_motion[i], correspondences);
    }
  }

  free(correspondences);

  // Return true if any one of the motions has inliers.
  for (i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] > 0) return 1;
  }
  return 0;
}

// Don't use points around the frame border since they are less reliable
static INLINE int valid_point(int x, int y, int width, int height) {
  return (x > (PATCH_SIZE + PATCH_CENTER)) &&
         (x < (width - PATCH_SIZE - PATCH_CENTER)) &&
         (y > (PATCH_SIZE + PATCH_CENTER)) &&
         (y < (height - PATCH_SIZE - PATCH_CENTER));
}

static int determine_disflow_correspondence(int *frm_corners,
                                            int num_frm_corners, double *flow_u,
                                            double *flow_v, int width,
                                            int height, int stride,
                                            double *correspondences) {
  int num_correspondences = 0;
  int x, y;
  for (int i = 0; i < num_frm_corners; ++i) {
    x = frm_corners[2 * i];
    y = frm_corners[2 * i + 1];
    if (valid_point(x, y, width, height)) {
      correspondences[4 * num_correspondences] = x;
      correspondences[4 * num_correspondences + 1] = y;
      correspondences[4 * num_correspondences + 2] = x + flow_u[y * stride + x];
      correspondences[4 * num_correspondences + 3] = y + flow_v[y * stride + x];
      num_correspondences++;
    }
  }
  return num_correspondences;
}

static double getCubicValue(double p[4], double x) {
  return p[1] + 0.5 * x *
                    (p[2] - p[0] +
                     x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                          x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

static void get_subcolumn(unsigned char *ref, double col[4], int stride, int x,
                          int y_start) {
  int i;
  for (i = 0; i < 4; ++i) {
    col[i] = ref[(i + y_start) * stride + x];
  }
}

static double bicubic(unsigned char *ref, double x, double y, int stride) {
  double arr[4];
  int k;
  int i = (int)x;
  int j = (int)y;
  for (k = 0; k < 4; ++k) {
    double arr_temp[4];
    get_subcolumn(ref, arr_temp, stride, i + k - 1, j - 1);
    arr[k] = getCubicValue(arr_temp, y - j);
  }
  return getCubicValue(arr, x - i);
}

// Interpolate a warped block using bicubic interpolation when possible
static unsigned char interpolate(unsigned char *ref, double x, double y,
                                 int width, int height, int stride) {
  if (x < 0 && y < 0)
    return ref[0];
  else if (x < 0 && y > height - 1)
    return ref[(height - 1) * stride];
  else if (x > width - 1 && y < 0)
    return ref[width - 1];
  else if (x > width - 1 && y > height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (x < 0) {
    int v;
    int i = (int)y;
    double a = y - i;
    if (y > 1 && y < height - 2) {
      double arr[4];
      get_subcolumn(ref, arr, stride, 0, i - 1);
      return clamp((int)(getCubicValue(arr, a) + 0.5), 0, 255);
    }
    v = (int)(ref[i * stride] * (1 - a) + ref[(i + 1) * stride] * a + 0.5);
    return clamp(v, 0, 255);
  } else if (y < 0) {
    int v;
    int j = (int)x;
    double b = x - j;
    if (x > 1 && x < width - 2) {
      double arr[4] = { ref[j - 1], ref[j], ref[j + 1], ref[j + 2] };
      return clamp((int)(getCubicValue(arr, b) + 0.5), 0, 255);
    }
    v = (int)(ref[j] * (1 - b) + ref[j + 1] * b + 0.5);
    return clamp(v, 0, 255);
  } else if (x > width - 1) {
    int v;
    int i = (int)y;
    double a = y - i;
    if (y > 1 && y < height - 2) {
      double arr[4];
      get_subcolumn(ref, arr, stride, width - 1, i - 1);
      return clamp((int)(getCubicValue(arr, a) + 0.5), 0, 255);
    }
    v = (int)(ref[i * stride + width - 1] * (1 - a) +
              ref[(i + 1) * stride + width - 1] * a + 0.5);
    return clamp(v, 0, 255);
  } else if (y > height - 1) {
    int v;
    int j = (int)x;
    double b = x - j;
    if (x > 1 && x < width - 2) {
      int row = (height - 1) * stride;
      double arr[4] = { ref[row + j - 1], ref[row + j], ref[row + j + 1],
                        ref[row + j + 2] };
      return clamp((int)(getCubicValue(arr, b) + 0.5), 0, 255);
    }
    v = (int)(ref[(height - 1) * stride + j] * (1 - b) +
              ref[(height - 1) * stride + j + 1] * b + 0.5);
    return clamp(v, 0, 255);
  } else if (x > 1 && y > 1 && x < width - 2 && y < height - 2) {
    return clamp((int)(bicubic(ref, x, y, stride) + 0.5), 0, 255);
  } else {
    int i = (int)y;
    int j = (int)x;
    double a = y - i;
    double b = x - j;
    int v = (int)(ref[i * stride + j] * (1 - a) * (1 - b) +
                  ref[i * stride + j + 1] * (1 - a) * b +
                  ref[(i + 1) * stride + j] * a * (1 - b) +
                  ref[(i + 1) * stride + j + 1] * a * b);
    return clamp(v, 0, 255);
  }
}

// Warps a block using flow vector [u, v] and computes the mse
static double compute_warp_and_error(unsigned char *ref, unsigned char *frm,
                                     int width, int height, int stride, int x,
                                     int y, double u, double v, int16_t *dt) {
  int i, j;
  unsigned char warped;
  double x_w, y_w;
  double mse = 0;
  int16_t err = 0;
  for (i = y; i < y + PATCH_SIZE; ++i)
    for (j = x; j < x + PATCH_SIZE; ++j) {
      x_w = (double)j + u;
      y_w = (double)i + v;
      warped = interpolate(ref, x_w, y_w, width, height, stride);
      err = warped - frm[j + i * stride];
      mse += err * err;
      dt[(i - y) * PATCH_SIZE + (j - x)] = err;
    }

  mse /= (PATCH_SIZE * PATCH_SIZE);
  return mse;
}

// Computes the components of the system of equations used to solve for
// a flow vector. This includes:
// 1.) The hessian matrix for optical flow. This matrix is in the
// form of:
//
//       M = |sum(dx * dx)  sum(dx * dy)|
//           |sum(dx * dy)  sum(dy * dy)|
//
// 2.)   b = |sum(dx * dt)|
//           |sum(dy * dt)|
// Where the sums are computed over a square window of PATCH_SIZE.
static INLINE void compute_flow_system(const double *dx, int dx_stride,
                                       const double *dy, int dy_stride,
                                       const int16_t *dt, int dt_stride,
                                       double *M, double *b) {
  for (int i = 0; i < PATCH_SIZE; i++) {
    for (int j = 0; j < PATCH_SIZE; j++) {
      M[0] += dx[i * dx_stride + j] * dx[i * dx_stride + j];
      M[1] += dx[i * dx_stride + j] * dy[i * dy_stride + j];
      M[3] += dy[i * dy_stride + j] * dy[i * dy_stride + j];

      b[0] += dx[i * dx_stride + j] * dt[i * dt_stride + j];
      b[1] += dy[i * dy_stride + j] * dt[i * dt_stride + j];
    }
  }

  M[2] = M[1];
}

// Solves a general Mx = b where M is a 2x2 matrix and b is a 2x1 matrix
static INLINE void solve_2x2_system(const double *M, const double *b,
                                    double *output_vec) {
  double M_0 = M[0];
  double M_3 = M[3];
  double det = (M_0 * M_3) - (M[1] * M[2]);
  if (det < 1e-5) {
    // Handle singular matrix
    // TODO(sarahparker) compare results using pseudo inverse instead
    M_0 += 1e-10;
    M_3 += 1e-10;
    det = (M_0 * M_3) - (M[1] * M[2]);
  }
  const double det_inv = 1 / det;
  const double mult_b0 = det_inv * b[0];
  const double mult_b1 = det_inv * b[1];
  output_vec[0] = M_3 * mult_b0 - M[1] * mult_b1;
  output_vec[1] = -M[2] * mult_b0 + M_0 * mult_b1;
}

/*
static INLINE void image_difference(const uint8_t *src, int src_stride,
                                    const uint8_t *ref, int ref_stride,
                                    int16_t *dst, int dst_stride, int height,
                                    int width) {
  const int block_unit = 8;
  // Take difference in 8x8 blocks to make use of optimized diff function
  for (int i = 0; i < height; i += block_unit) {
    for (int j = 0; j < width; j += block_unit) {
      aom_subtract_block(block_unit, block_unit, dst + i * dst_stride + j,
                         dst_stride, src + i * src_stride + j, src_stride,
                         ref + i * ref_stride + j, ref_stride);
    }
  }
}
*/

// Compute an image gradient using a sobel filter.
// If dir == 1, compute the x gradient. If dir == 0, compute y. This function
// assumes the images have been padded so that they can be processed in units
// of 8.
static INLINE void sobel_xy_image_gradient(const uint8_t *src, int src_stride,
                                           double *dst, int dst_stride,
                                           int height, int width, int dir) {
  double norm = 1.0;
  // TODO(sarahparker) experiment with doing this over larger block sizes
  const int block_unit = 8;
  // Filter in 8x8 blocks to eventually make use of optimized convolve function
  for (int i = 0; i < height; i += block_unit) {
    for (int j = 0; j < width; j += block_unit) {
      av1_convolve_2d_sobel_y_c(src + i * src_stride + j, src_stride,
                                dst + i * dst_stride + j, dst_stride,
                                block_unit, block_unit, dir, norm);
    }
  }
}

static ImagePyramid *alloc_pyramid(int width, int height, int pad_size,
                                   int compute_gradient) {
  ImagePyramid *pyr = aom_malloc(sizeof(*pyr));
  pyr->has_gradient = compute_gradient;
  // 2 * width * height is the upper bound for a buffer that fits
  // all pyramid levels + padding for each level
  const int buffer_size = sizeof(*pyr->level_buffer) * 2 * width * height +
                          (width + 2 * pad_size) * 2 * pad_size * N_LEVELS;
  pyr->level_buffer = aom_malloc(buffer_size);
  memset(pyr->level_buffer, 0, buffer_size);

  if (compute_gradient) {
    const int gradient_size =
        sizeof(*pyr->level_dx_buffer) * 2 * width * height +
        (width + 2 * pad_size) * 2 * pad_size * N_LEVELS;
    pyr->level_dx_buffer = aom_malloc(gradient_size);
    pyr->level_dy_buffer = aom_malloc(gradient_size);
    memset(pyr->level_dx_buffer, 0, gradient_size);
    memset(pyr->level_dy_buffer, 0, gradient_size);
  }
  return pyr;
}

static void free_pyramid(ImagePyramid *pyr) {
  aom_free(pyr->level_buffer);
  if (pyr->has_gradient) {
    aom_free(pyr->level_dx_buffer);
    aom_free(pyr->level_dy_buffer);
  }
  aom_free(pyr);
}

static INLINE void update_level_dims(ImagePyramid *frm_pyr, int level) {
  frm_pyr->widths[level] = frm_pyr->widths[level - 1] >> 1;
  frm_pyr->heights[level] = frm_pyr->heights[level - 1] >> 1;
  frm_pyr->strides[level] = frm_pyr->widths[level] + 2 * frm_pyr->pad_size;
  // Point the beginning of the next level buffer to the correct location inside
  // the padded border
  frm_pyr->level_loc[level] =
      frm_pyr->level_loc[level - 1] +
      frm_pyr->strides[level - 1] *
          (2 * frm_pyr->pad_size + frm_pyr->heights[level - 1]);
}

// Compute coarse to fine pyramids for a frame
static void compute_flow_pyramids(unsigned char *frm, const int frm_width,
                                  const int frm_height, const int frm_stride,
                                  int n_levels, int pad_size, int compute_grad,
                                  ImagePyramid *frm_pyr) {
  int cur_width, cur_height, cur_stride, cur_loc;
  assert((frm_width >> n_levels) > 0);
  assert((frm_height >> n_levels) > 0);

  // Initialize first level
  frm_pyr->n_levels = n_levels;
  frm_pyr->pad_size = pad_size;
  frm_pyr->widths[0] = frm_width;
  frm_pyr->heights[0] = frm_height;
  frm_pyr->strides[0] = frm_width + 2 * frm_pyr->pad_size;
  // Point the beginning of the level buffer to the location inside
  // the padded border
  frm_pyr->level_loc[0] =
      frm_pyr->strides[0] * frm_pyr->pad_size + frm_pyr->pad_size;
  // This essentially copies the original buffer into the pyramid buffer
  // without the original padding
  av1_resize_plane(frm, frm_height, frm_width, frm_stride,
                   frm_pyr->level_buffer + frm_pyr->level_loc[0],
                   frm_pyr->heights[0], frm_pyr->widths[0],
                   frm_pyr->strides[0]);

  if (compute_grad) {
    cur_width = frm_pyr->widths[0];
    cur_height = frm_pyr->heights[0];
    cur_stride = frm_pyr->strides[0];
    cur_loc = frm_pyr->level_loc[0];
    assert(frm_pyr->has_gradient && frm_pyr->level_dx_buffer != NULL &&
           frm_pyr->level_dy_buffer != NULL);
    // Computation x gradient
    sobel_xy_image_gradient(frm_pyr->level_buffer + cur_loc, cur_stride,
                            frm_pyr->level_dx_buffer + cur_loc, cur_stride,
                            cur_height, cur_width, 1);

    // Computation y gradient
    sobel_xy_image_gradient(frm_pyr->level_buffer + cur_loc, cur_stride,
                            frm_pyr->level_dy_buffer + cur_loc, cur_stride,
                            cur_height, cur_width, 0);
  }

  // Start at the finest level and resize down to the coarsest level
  for (int level = 1; level < n_levels; ++level) {
    update_level_dims(frm_pyr, level);
    cur_width = frm_pyr->widths[level];
    cur_height = frm_pyr->heights[level];
    cur_stride = frm_pyr->strides[level];
    cur_loc = frm_pyr->level_loc[level];

    av1_resize_plane(frm_pyr->level_buffer + frm_pyr->level_loc[level - 1],
                     frm_pyr->heights[level - 1], frm_pyr->widths[level - 1],
                     frm_pyr->strides[level - 1],
                     frm_pyr->level_buffer + cur_loc, cur_height, cur_width,
                     cur_stride);

    if (compute_grad) {
      assert(frm_pyr->has_gradient && frm_pyr->level_dx_buffer != NULL &&
             frm_pyr->level_dy_buffer != NULL);
      // Computation x gradient
      sobel_xy_image_gradient(frm_pyr->level_buffer + cur_loc, cur_stride,
                              frm_pyr->level_dx_buffer + cur_loc, cur_stride,
                              cur_height, cur_width, 1);

      // Computation y gradient
      sobel_xy_image_gradient(frm_pyr->level_buffer + cur_loc, cur_stride,
                              frm_pyr->level_dy_buffer + cur_loc, cur_stride,
                              cur_height, cur_width, 0);
    }
  }
}

static INLINE void compute_flow_at_point(unsigned char *frm, unsigned char *ref,
                                         double *dx, double *dy, int x, int y,
                                         int width, int height, int stride,
                                         double *u, double *v) {
  double M[4] = { 0 };
  double b[2] = { 0 };
  double tmp_output_vec[2] = { 0 };
  double error = 0;
  int16_t dt[PATCH_SIZE * PATCH_SIZE];
  double o_u = *u;
  double o_v = *v;

  for (int itr = 0; itr < DISFLOW_MAX_ITR; itr++) {
    error = compute_warp_and_error(ref, frm, width, height, stride, x, y, *u,
                                   *v, dt);
    if (error <= DISFLOW_ERROR_TR) break;
    compute_flow_system(dx, stride, dy, stride, dt, PATCH_SIZE, M, b);
    solve_2x2_system(M, b, tmp_output_vec);
    *u += tmp_output_vec[0];
    *v += tmp_output_vec[1];
  }
  if (fabs(*u - o_u) > PATCH_SIZE || fabs(*v - o_u) > PATCH_SIZE) {
    *u = o_u;
    *v = o_v;
  }
}

// make sure flow_u and flow_v start at 0
static void compute_flow_field(ImagePyramid *frm_pyr, ImagePyramid *ref_pyr,
                               double *flow_u, double *flow_v) {
  int cur_width, cur_height, cur_stride, cur_loc, patch_loc, patch_center;
  double *u_upscale =
      aom_malloc(frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_u));
  double *v_upscale =
      aom_malloc(frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_v));

  assert(frm_pyr->n_levels == ref_pyr->n_levels);

  // Compute flow field from coarsest to finest level of the pyramid
  for (int level = frm_pyr->n_levels - 1; level >= 0; --level) {
    cur_width = frm_pyr->widths[level];
    cur_height = frm_pyr->heights[level];
    cur_stride = frm_pyr->strides[level];
    cur_loc = frm_pyr->level_loc[level];

    for (int i = PATCH_SIZE; i < cur_height - PATCH_SIZE; i += PATCH_STEP) {
      for (int j = PATCH_SIZE; j < cur_width - PATCH_SIZE; j += PATCH_STEP) {
        patch_loc = i * cur_stride + j;
        patch_center = patch_loc + PATCH_CENTER * cur_stride + PATCH_CENTER;
        compute_flow_at_point(frm_pyr->level_buffer + cur_loc,
                              ref_pyr->level_buffer + cur_loc,
                              frm_pyr->level_dx_buffer + cur_loc + patch_loc,
                              frm_pyr->level_dy_buffer + cur_loc + patch_loc, j,
                              i, cur_width, cur_height, cur_stride,
                              flow_u + patch_center, flow_v + patch_center);
      }
    }
    // TODO(sarahparker) Replace this with upscale function in resize.c
    if (level > 0) {
      int h_upscale = frm_pyr->heights[level - 1];
      int w_upscale = frm_pyr->widths[level - 1];
      int s_upscale = frm_pyr->strides[level - 1];
      for (int i = 0; i < h_upscale; ++i) {
        for (int j = 0; j < w_upscale; ++j) {
          u_upscale[j + i * s_upscale] =
              flow_u[(int)(j >> 1) + (int)(i >> 1) * cur_stride];
          v_upscale[j + i * s_upscale] =
              flow_v[(int)(j >> 1) + (int)(i >> 1) * cur_stride];
        }
      }
      memcpy(flow_u, u_upscale,
             frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_u));
      memcpy(flow_v, v_upscale,
             frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_v));
    }
  }
  aom_free(u_upscale);
  aom_free(v_upscale);
}

static int compute_global_motion_disflow_based(
    TransformationType type, unsigned char *frm_buffer, int frm_width,
    int frm_height, int frm_stride, int *frm_corners, int num_frm_corners,
    YV12_BUFFER_CONFIG *ref, int bit_depth, int *num_inliers_by_motion,
    MotionModel *params_by_motion, int num_motions) {
  unsigned char *ref_buffer = ref->y_buffer;
  const int ref_width = ref->y_width;
  const int ref_height = ref->y_height;
  const int pad_size = AOMMAX(PATCH_SIZE, MIN_PAD);
  int num_correspondences;
  double *correspondences;
  RansacFuncDouble ransac = av1_get_ransac_double_prec_type(type);
  assert(frm_width == ref_width);
  assert(frm_height == ref_height);

  // Ensure the number of pyramid levels will work with the frame resolution
  const int msb =
      frm_width < frm_height ? get_msb(frm_width) : get_msb(frm_height);
  const int n_levels = AOMMIN(msb, N_LEVELS);

  if (ref->flags & YV12_FLAG_HIGHBITDEPTH) {
    ref_buffer = av1_downconvert_frame(ref, bit_depth);
  }

  // TODO(sarahparker) We will want to do the source pyramid computation
  // outside of this function so it doesn't get recomputed for every
  // reference. We also don't need to compute every pyramid level for the
  // reference in advance, since lower levels can be overwritten once their
  // flow field is computed and upscaled. I'll add these optimizations
  // once the full implementation is working.
  // Allocate frm image pyramids
  int compute_gradient = 1;
  ImagePyramid *frm_pyr =
      alloc_pyramid(frm_width, frm_height, pad_size, compute_gradient);
  compute_flow_pyramids(frm_buffer, frm_width, frm_height, frm_stride, n_levels,
                        pad_size, compute_gradient, frm_pyr);
  // Allocate ref image pyramids
  compute_gradient = 0;
  ImagePyramid *ref_pyr =
      alloc_pyramid(ref_width, ref_height, pad_size, compute_gradient);
  compute_flow_pyramids(ref_buffer, ref_width, ref_height, ref->y_stride,
                        n_levels, pad_size, compute_gradient, ref_pyr);

  double *flow_u =
      aom_malloc(frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_u));
  double *flow_v =
      aom_malloc(frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_v));

  memset(flow_u, 0,
         frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_u));
  memset(flow_v, 0,
         frm_pyr->strides[0] * frm_pyr->heights[0] * sizeof(*flow_v));

  compute_flow_field(frm_pyr, ref_pyr, flow_u, flow_v);

  // find correspondences between the two images using the flow field
  correspondences = aom_malloc(num_frm_corners * 4 * sizeof(*correspondences));
  num_correspondences = determine_disflow_correspondence(
      frm_corners, num_frm_corners, flow_u, flow_v, frm_width, frm_height,
      frm_pyr->strides[0], correspondences);
  ransac(correspondences, num_correspondences, num_inliers_by_motion,
         params_by_motion, num_motions);

  free_pyramid(frm_pyr);
  free_pyramid(ref_pyr);
  aom_free(correspondences);
  aom_free(flow_u);
  aom_free(flow_v);
  // Set num_inliers = 0 for motions with too few inliers so they are ignored.
  for (int i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] < MIN_INLIER_PROB * num_correspondences) {
      num_inliers_by_motion[i] = 0;
    }
  }

  // Return true if any one of the motions has inliers.
  for (int i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] > 0) return 1;
  }
  return 0;
}

/*
Computes "num_motions" candidate global motion parameters between two frames.
The array "params_by_motion" should be length 8 * "num_motions". The ordering
of each set of parameters is best described  by the homography:

[x'     (m2 m3 m0   [x
z .  y'  =   m4 m5 m1 *  y
1]      m6 m7 1)    1]

where m{i} represents the ith value in any given set of parameters.

"num_inliers" should be length "num_motions", and will be populated with the
number of inlier feature points for each motion. Params for which the
num_inliers entry is 0 should be ignored by the caller.
*/
static AOM_INLINE int compute_global_motion(
    TransformationType type, unsigned char *src_buffer, int src_width,
    int src_height, int src_stride, int *src_corners, int num_src_corners,
    YV12_BUFFER_CONFIG *ref, int bit_depth,
    GlobalMotionEstimationType gm_estimation_type, int *num_inliers_by_motion,
    MotionModel *params_by_motion, int num_motions) {
  switch (gm_estimation_type) {
    case GLOBAL_MOTION_FEATURE_BASED:
      return compute_global_motion_feature_based(
          type, src_buffer, src_width, src_height, src_stride, src_corners,
          num_src_corners, ref, bit_depth, num_inliers_by_motion,
          params_by_motion, num_motions);
    case GLOBAL_MOTION_DISFLOW_BASED:
      return compute_global_motion_disflow_based(
          type, src_buffer, src_width, src_height, src_stride, src_corners,
          num_src_corners, ref, bit_depth, num_inliers_by_motion,
          params_by_motion, num_motions);
    default: assert(0 && "Unknown global motion estimation type");
  }
  return 0;
}

static int gm_get_params_cost(const WarpedMotionParams *gm,
                              const WarpedMotionParams *ref_gm, int allow_hp) {
  int params_cost = 0;
  int trans_bits, trans_prec_diff;
  switch (gm->wmtype) {
    case AFFINE:
    case ROTZOOM:
      params_cost += aom_count_signed_primitive_refsubexpfin(
          GM_ALPHA_MAX + 1, SUBEXPFIN_K,
          (ref_gm->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS),
          (gm->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
      params_cost += aom_count_signed_primitive_refsubexpfin(
          GM_ALPHA_MAX + 1, SUBEXPFIN_K,
          (ref_gm->wmmat[3] >> GM_ALPHA_PREC_DIFF),
          (gm->wmmat[3] >> GM_ALPHA_PREC_DIFF));
      if (gm->wmtype >= AFFINE) {
        params_cost += aom_count_signed_primitive_refsubexpfin(
            GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            (ref_gm->wmmat[4] >> GM_ALPHA_PREC_DIFF),
            (gm->wmmat[4] >> GM_ALPHA_PREC_DIFF));
        params_cost += aom_count_signed_primitive_refsubexpfin(
            GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            (ref_gm->wmmat[5] >> GM_ALPHA_PREC_DIFF) -
                (1 << GM_ALPHA_PREC_BITS),
            (gm->wmmat[5] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
      }
      AOM_FALLTHROUGH_INTENDED;
    case TRANSLATION:
      trans_bits = (gm->wmtype == TRANSLATION)
                       ? GM_ABS_TRANS_ONLY_BITS - !allow_hp
                       : GM_ABS_TRANS_BITS;
      trans_prec_diff = (gm->wmtype == TRANSLATION)
                            ? GM_TRANS_ONLY_PREC_DIFF + !allow_hp
                            : GM_TRANS_PREC_DIFF;
      params_cost += aom_count_signed_primitive_refsubexpfin(
          (1 << trans_bits) + 1, SUBEXPFIN_K,
          (ref_gm->wmmat[0] >> trans_prec_diff),
          (gm->wmmat[0] >> trans_prec_diff));
      params_cost += aom_count_signed_primitive_refsubexpfin(
          (1 << trans_bits) + 1, SUBEXPFIN_K,
          (ref_gm->wmmat[1] >> trans_prec_diff),
          (gm->wmmat[1] >> trans_prec_diff));
      AOM_FALLTHROUGH_INTENDED;
    case IDENTITY: break;
    default: assert(0);
  }
  return (params_cost << AV1_PROB_COST_SHIFT);
}

// TODO(Remya): Can include erroradv_prod_tr[] for threshold calculation
static INLINE int64_t calc_erroradv_threshold(AV1_COMP *cpi,
                                              int64_t ref_frame_error) {
  if (!cpi->sf.gm_sf.disable_adaptive_warp_error_thresh)
    return (int64_t)(
        ref_frame_error * erroradv_tr[cpi->sf.gm_sf.gm_erroradv_type] + 0.5);
  else
    return INT64_MAX;
}

static void compute_global_motion_for_ref_frame(
    AV1_COMP *cpi, YV12_BUFFER_CONFIG *ref_buf[REF_FRAMES], int frame,
    int num_src_corners, int *src_corners, unsigned char *src_buffer,
    MotionModel *params_by_motion, uint8_t *segment_map,
    const int segment_map_w, const int segment_map_h,
    const WarpedMotionParams *ref_params) {
  ThreadData *const td = &cpi->td;
  MACROBLOCK *const x = &td->mb;
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  int i;
  int src_width = cpi->source->y_width;
  int src_height = cpi->source->y_height;
  int src_stride = cpi->source->y_stride;
  // clang-format off
	static const double kIdentityParams[MAX_PARAMDIM - 1] = {
		0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
	};
  // clang-format on
  WarpedMotionParams tmp_wm_params;
  const double *params_this_motion;
  int inliers_by_motion[RANSAC_NUM_MOTIONS];
  assert(ref_buf[frame] != NULL);
  TransformationType model;

  aom_clear_system_state();

  // TODO(sarahparker, debargha): Explore do_adaptive_gm_estimation = 1
  const int do_adaptive_gm_estimation = 0;

  const int ref_frame_dist = get_relative_dist(
      &cm->seq_params.order_hint_info, cm->current_frame.order_hint,
      cm->cur_frame->ref_order_hints[frame - LAST_FRAME]);
  const GlobalMotionEstimationType gm_estimation_type =
      cm->seq_params.order_hint_info.enable_order_hint &&
              abs(ref_frame_dist) <= 2 && do_adaptive_gm_estimation
          ? GLOBAL_MOTION_DISFLOW_BASED
          : GLOBAL_MOTION_FEATURE_BASED;
  for (model = ROTZOOM; model < GLOBAL_TRANS_TYPES_ENC; ++model) {
    int64_t best_warp_error = INT64_MAX;
    // Initially set all params to identity.
    for (i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
      memcpy(params_by_motion[i].params, kIdentityParams,
             (MAX_PARAMDIM - 1) * sizeof(*(params_by_motion[i].params)));
      params_by_motion[i].num_inliers = 0;
    }

    compute_global_motion(model, src_buffer, src_width, src_height, src_stride,
                          src_corners, num_src_corners, ref_buf[frame],
                          cpi->common.seq_params.bit_depth, gm_estimation_type,
                          inliers_by_motion, params_by_motion,
                          RANSAC_NUM_MOTIONS);
    int64_t ref_frame_error = 0;
    for (i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
      if (inliers_by_motion[i] == 0) continue;

      params_this_motion = params_by_motion[i].params;
      av1_convert_model_to_params(params_this_motion, &tmp_wm_params);

      if (tmp_wm_params.wmtype != IDENTITY) {
        av1_compute_feature_segmentation_map(
            segment_map, segment_map_w, segment_map_h,
            params_by_motion[i].inliers, params_by_motion[i].num_inliers);

        ref_frame_error = av1_segmented_frame_error(
            is_cur_buf_hbd(xd), xd->bd, ref_buf[frame]->y_buffer,
            ref_buf[frame]->y_stride, cpi->source->y_buffer, src_width,
            src_height, src_stride, segment_map, segment_map_w);

        int64_t erroradv_threshold =
            calc_erroradv_threshold(cpi, ref_frame_error);

        const int64_t warp_error = av1_refine_integerized_param(
            &tmp_wm_params, tmp_wm_params.wmtype, is_cur_buf_hbd(xd), xd->bd,
            ref_buf[frame]->y_buffer, ref_buf[frame]->y_width,
            ref_buf[frame]->y_height, ref_buf[frame]->y_stride,
            cpi->source->y_buffer, src_width, src_height, src_stride,
            GM_REFINEMENT_COUNT, best_warp_error, segment_map, segment_map_w,
            erroradv_threshold);

        if (warp_error < best_warp_error) {
          best_warp_error = warp_error;
          // Save the wm_params modified by
          // av1_refine_integerized_param() rather than motion index to
          // avoid rerunning refine() below.
          memcpy(&(cm->global_motion[frame]), &tmp_wm_params,
                 sizeof(WarpedMotionParams));
        }
      }
    }
    if (cm->global_motion[frame].wmtype <= AFFINE)
      if (!av1_get_shear_params(&cm->global_motion[frame]))
        cm->global_motion[frame] = default_warp_params;

    if (cm->global_motion[frame].wmtype == TRANSLATION) {
      cm->global_motion[frame].wmmat[0] =
          convert_to_trans_prec(cm->features.allow_high_precision_mv,
                                cm->global_motion[frame].wmmat[0]) *
          GM_TRANS_ONLY_DECODE_FACTOR;
      cm->global_motion[frame].wmmat[1] =
          convert_to_trans_prec(cm->features.allow_high_precision_mv,
                                cm->global_motion[frame].wmmat[1]) *
          GM_TRANS_ONLY_DECODE_FACTOR;
    }

    if (cm->global_motion[frame].wmtype == IDENTITY) continue;

    if (ref_frame_error == 0) continue;

    // If the best error advantage found doesn't meet the threshold for
    // this motion type, revert to IDENTITY.
    if (!av1_is_enough_erroradvantage(
            (double)best_warp_error / ref_frame_error,
            gm_get_params_cost(&cm->global_motion[frame], ref_params,
                               cm->features.allow_high_precision_mv),
            cpi->sf.gm_sf.gm_erroradv_type)) {
      cm->global_motion[frame] = default_warp_params;
    }

    if (cm->global_motion[frame].wmtype != IDENTITY) break;
  }

  aom_clear_system_state();
}

static INLINE void compute_gm_for_valid_ref_frames(
    AV1_COMP *cpi, YV12_BUFFER_CONFIG *ref_buf[REF_FRAMES], int frame,
    int num_src_corners, int *src_corners, unsigned char *src_buffer,
    MotionModel *params_by_motion, uint8_t *segment_map,
    const int segment_map_w, const int segment_map_h) {
  AV1_COMMON *const cm = &cpi->common;
  GlobalMotionInfo *const gm_info = &cpi->gm_info;
  const WarpedMotionParams *ref_params =
      cm->prev_frame ? &cm->prev_frame->global_motion[frame]
                     : &default_warp_params;

  compute_global_motion_for_ref_frame(
      cpi, ref_buf, frame, num_src_corners, src_corners, src_buffer,
      params_by_motion, segment_map, segment_map_w, segment_map_h, ref_params);

  gm_info->params_cost[frame] =
      gm_get_params_cost(&cm->global_motion[frame], ref_params,
                         cm->features.allow_high_precision_mv) +
      gm_info->type_cost[cm->global_motion[frame].wmtype] -
      gm_info->type_cost[IDENTITY];
}

static INLINE void compute_global_motion_for_references(
    AV1_COMP *cpi, YV12_BUFFER_CONFIG *ref_buf[REF_FRAMES],
    FrameDistPair reference_frame[REF_FRAMES - 1], int num_ref_frames,
    int num_src_corners, int *src_corners, unsigned char *src_buffer,
    MotionModel *params_by_motion, uint8_t *segment_map,
    const int segment_map_w, const int segment_map_h) {
  // Computation of frame corners for the source frame will be done already.
  assert(num_src_corners != -1);
  AV1_COMMON *const cm = &cpi->common;
  // Compute global motion w.r.t. reference frames starting from the nearest ref
  // frame in a given direction.
  for (int frame = 0; frame < num_ref_frames; frame++) {
    int ref_frame = reference_frame[frame].frame;
    compute_gm_for_valid_ref_frames(cpi, ref_buf, ref_frame, num_src_corners,
                                    src_corners, src_buffer, params_by_motion,
                                    segment_map, segment_map_w, segment_map_h);
    // If global motion w.r.t. current ref frame is
    // INVALID/TRANSLATION/IDENTITY, skip the evaluation of global motion w.r.t
    // the remaining ref frames in that direction. The below exit is disabled
    // when ref frame distance w.r.t. current frame is zero. E.g.:
    // source_alt_ref_frame w.r.t. ARF frames.
    if (cpi->sf.gm_sf.prune_ref_frame_for_gm_search &&
        reference_frame[frame].distance != 0 &&
        cm->global_motion[ref_frame].wmtype != ROTZOOM)
      break;
  }
}

static int compare_distance(const void *a, const void *b) {
  const int diff =
      ((FrameDistPair *)a)->distance - ((FrameDistPair *)b)->distance;
  if (diff > 0)
    return 1;
  else if (diff < 0)
    return -1;
  return 0;
}

// Function to decide if we can skip the global motion parameter computation
// for a particular ref frame
static INLINE int skip_gm_frame(AV1_COMMON *const cm, int ref_frame) {
  if ((ref_frame == LAST3_FRAME || ref_frame == LAST2_FRAME) &&
      cm->global_motion[GOLDEN_FRAME].wmtype != IDENTITY) {
    return get_relative_dist(
               &cm->seq_params.order_hint_info,
               cm->cur_frame->ref_order_hints[ref_frame - LAST_FRAME],
               cm->cur_frame->ref_order_hints[GOLDEN_FRAME - LAST_FRAME]) <= 0;
  }
  return 0;
}

static int do_gm_search_logic(SPEED_FEATURES *const sf, int frame) {
  (void)frame;
  switch (sf->gm_sf.gm_search_type) {
    case GM_FULL_SEARCH: return 1;
    case GM_REDUCED_REF_SEARCH_SKIP_L2_L3:
      return !(frame == LAST2_FRAME || frame == LAST3_FRAME);
    case GM_REDUCED_REF_SEARCH_SKIP_L2_L3_ARF2:
      return !(frame == LAST2_FRAME || frame == LAST3_FRAME ||
               (frame == ALTREF2_FRAME));
    case GM_DISABLE_SEARCH: return 0;
    default: assert(0);
  }
  return 1;
}

static INLINE void update_valid_ref_frames_for_gm(
    AV1_COMP *cpi, YV12_BUFFER_CONFIG *ref_buf[REF_FRAMES],
    FrameDistPair reference_frames[MAX_DIRECTIONS][REF_FRAMES - 1],
    int *num_ref_frames) {
  AV1_COMMON *const cm = &cpi->common;
  const OrderHintInfo *const order_hint_info = &cm->seq_params.order_hint_info;
  int *num_past_ref_frames = &num_ref_frames[0];
  int *num_future_ref_frames = &num_ref_frames[1];
  for (int frame = ALTREF_FRAME; frame >= LAST_FRAME; --frame) {
    const MV_REFERENCE_FRAME ref_frame[2] = { frame, NONE_FRAME };
    RefCntBuffer *buf = get_ref_frame_buf(cm, frame);
    const int ref_disabled =
        !(cpi->ref_frame_flags & av1_ref_frame_flag_list[frame]);
    ref_buf[frame] = NULL;
    cm->global_motion[frame] = default_warp_params;
    // Skip global motion estimation for invalid ref frames
    if (buf == NULL ||
        (ref_disabled && cpi->sf.hl_sf.recode_loop != DISALLOW_RECODE)) {
      cpi->gm_info.params_cost[frame] = 0;
      continue;
    } else {
      ref_buf[frame] = &buf->buf;
    }

    if (ref_buf[frame]->y_crop_width == cpi->source->y_crop_width &&
        ref_buf[frame]->y_crop_height == cpi->source->y_crop_height &&
        do_gm_search_logic(&cpi->sf, frame) &&
        !prune_ref_by_selective_ref_frame(
            cpi, NULL, ref_frame, cm->cur_frame->ref_display_order_hint) &&
        !(cpi->sf.gm_sf.selective_ref_gm && skip_gm_frame(cm, frame))) {
      assert(ref_buf[frame] != NULL);
      int relative_frame_dist = av1_encoder_get_relative_dist(
          order_hint_info, buf->display_order_hint,
          cm->cur_frame->display_order_hint);
      // Populate past and future ref frames.
      // reference_frames[0][] indicates past direction and
      // reference_frames[1][] indicates future direction.
      if (relative_frame_dist <= 0) {
        reference_frames[0][*num_past_ref_frames].distance =
            abs(relative_frame_dist);
        reference_frames[0][*num_past_ref_frames].frame = frame;
        (*num_past_ref_frames)++;
      } else {
        reference_frames[1][*num_future_ref_frames].distance =
            abs(relative_frame_dist);
        reference_frames[1][*num_future_ref_frames].frame = frame;
        (*num_future_ref_frames)++;
      }
    }
  }
}

// Allocates and initializes memory for segment_map and MotionModel.
static AOM_INLINE void alloc_global_motion_data(MotionModel *params_by_motion,
                                                uint8_t **segment_map,
                                                const int segment_map_w,
                                                const int segment_map_h) {
  for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) {
    av1_zero(params_by_motion[m]);
    params_by_motion[m].inliers =
        aom_malloc(sizeof(*(params_by_motion[m].inliers)) * 2 * MAX_CORNERS);
  }

  *segment_map = (uint8_t *)aom_malloc(sizeof(*segment_map) * segment_map_w *
                                       segment_map_h);
  av1_zero_array(*segment_map, segment_map_w * segment_map_h);
}

// Deallocates segment_map and inliers.
static AOM_INLINE void dealloc_global_motion_data(MotionModel *params_by_motion,
                                                  uint8_t *segment_map) {
  aom_free(segment_map);

  for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) {
    aom_free(params_by_motion[m].inliers);
  }
}

void av1_compute_global_motion(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  GlobalMotionInfo *const gm_info = &cpi->gm_info;

  av1_zero(cpi->td.rd_counts.global_motion_used);
  av1_zero(gm_info->params_cost);
  if (cpi->common.current_frame.frame_type == INTER_FRAME && cpi->source &&
      cpi->oxcf.enable_global_motion && !gm_info->search_done) {
    YV12_BUFFER_CONFIG *ref_buf[REF_FRAMES];
    MotionModel params_by_motion[RANSAC_NUM_MOTIONS];
    uint8_t *segment_map = NULL;

    int num_src_corners = -1;
    int src_corners[2 * MAX_CORNERS];
    unsigned char *src_buffer = cpi->source->y_buffer;
    if (cpi->source->flags & YV12_FLAG_HIGHBITDEPTH) {
      // The source buffer is 16-bit, so we need to convert to 8 bits for the
      // following code. We cache the result until the source frame is released.
      src_buffer =
          av1_downconvert_frame(cpi->source, cpi->common.seq_params.bit_depth);
    }
    const int segment_map_w =
        (cpi->source->y_width + WARP_ERROR_BLOCK) >> WARP_ERROR_BLOCK_LOG;
    const int segment_map_h =
        (cpi->source->y_height + WARP_ERROR_BLOCK) >> WARP_ERROR_BLOCK_LOG;

    alloc_global_motion_data(params_by_motion, &segment_map, segment_map_w,
                             segment_map_h);

    FrameDistPair reference_frames[MAX_DIRECTIONS][REF_FRAMES - 1];
    memset(reference_frames, -1,
           sizeof(reference_frames[0][0]) * MAX_DIRECTIONS * (REF_FRAMES - 1));

    int num_ref_frames[MAX_DIRECTIONS] = { 0 };

    // Populate ref_buf for valid ref frames in global motion.
    update_valid_ref_frames_for_gm(cpi, ref_buf, reference_frames,
                                   num_ref_frames);

    // Sort the past and future ref frames in the ascending order of their
    // distance from the current frame. reference_frames[0] => past direction
    // and reference_frames[1] => future direction.
    qsort(reference_frames[0], num_ref_frames[0],
          sizeof(reference_frames[0][0]), compare_distance);
    qsort(reference_frames[1], num_ref_frames[1],
          sizeof(reference_frames[1][0]), compare_distance);

    // If atleast one valid reference frame exists in past/future directions,
    // compute interest points of source frame using FAST features.
    if (num_ref_frames[0] > 0 || num_ref_frames[1] > 0) {
      num_src_corners = av1_fast_corner_detect(
          src_buffer, cpi->source->y_width, cpi->source->y_height,
          cpi->source->y_stride, src_corners, MAX_CORNERS);
    }

    // Compute global motion w.r.t. past reference frames and future reference
    // frames.
    for (int dir = 0; dir < MAX_DIRECTIONS; dir++) {
      if (num_ref_frames[dir] > 0)
        compute_global_motion_for_references(
            cpi, ref_buf, reference_frames[dir], num_ref_frames[dir],
            num_src_corners, src_corners, src_buffer, params_by_motion,
            segment_map, segment_map_w, segment_map_h);
    }

    dealloc_global_motion_data(params_by_motion, segment_map);

    gm_info->search_done = 1;
  }
  memcpy(cm->cur_frame->global_motion, cm->global_motion,
         REF_FRAMES * sizeof(WarpedMotionParams));
}
