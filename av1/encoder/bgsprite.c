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

#include <float.h>
#include <math.h>

#include "av1/encoder/bgsprite.h"

#include "av1/common/mv.h"
#include "av1/common/warped_motion.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/global_motion.h"
#include "av1/encoder/mathutils.h"

// Helper function to convert the parameter array to a 3x3 matrix form
// target is modified
static void params_to_matrix(double *const params, double *target) {
  target[0] = params[2];
  target[1] = params[3];
  target[2] = params[0];
  target[3] = params[4];
  target[4] = params[5];
  target[5] = params[1];
  target[6] = params[6];
  target[7] = params[7];
  target[8] = 1;
}


// Helper function to convert a 3x3 matrix to a parameter array form
// target is modified
static void matrix_to_params(double *const matrix, double *target) {
  target[2] = matrix[0];
  target[3] = matrix[1];
  target[0] = matrix[2];
  target[4] = matrix[3];
  target[5] = matrix[4];
  target[1] = matrix[5];
  target[6] = matrix[6];
  target[7] = matrix[7];
}

// Helper function to do matrix multiplication on params
// target is modified,
// m1 and m2 are unmodified (unless m1 or m2 is also equal to target)
static void multiply_params(double *const m1, double *const m2,
                            double *target) {
  double m1_matrix[9];
  double m2_matrix[9];
  double result[9];

  params_to_matrix(m1, m1_matrix);
  params_to_matrix(m2, m2_matrix);
  multiply_mat(m2_matrix, m1_matrix, result, 3, 3, 3);
  matrix_to_params(result, target);
}

// Helper function to find limits of a single transformed image
// x_min, x_max, y_min, and y_max are modified
// width and height are the size of the input video
static void find_limits(int width, int height, double *const transform,
                        double *x_min, double *x_max, double *y_min,
                        double *y_max) {
  double transform_matrix[9];
  double xy_matrix[3] = {0, 0, 1};
  double uv_matrix[3] = {0};

  params_to_matrix(transform, transform_matrix);
  xy_matrix[0] = 0;
  xy_matrix[1] = 0;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, 3, 3, 1);
  *x_max = uv_matrix[0];
  *x_min = uv_matrix[0];
  *y_max = uv_matrix[1];
  *y_min = uv_matrix[1];

  xy_matrix[0] = width;
  xy_matrix[1] = 0;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, 3, 3, 1);
  if (uv_matrix[0] > *x_max) {
    *x_max = uv_matrix[0];
  }
  if (uv_matrix[0] < *x_min) {
    *x_min = uv_matrix[0];
  }
  if (uv_matrix[1] > *y_max) {
    *y_max = uv_matrix[1];
  }
  if (uv_matrix[1] < *y_min) {
    *y_min = uv_matrix[1];
  }

  xy_matrix[0] = width;
  xy_matrix[1] = height;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, 3, 3, 1);
  if (uv_matrix[0] > *x_max) {
    *x_max = uv_matrix[0];
  }
  if (uv_matrix[0] < *x_min) {
    *x_min = uv_matrix[0];
  }
  if (uv_matrix[1] > *y_max) {
    *y_max = uv_matrix[1];
  }
  if (uv_matrix[1] < *y_min) {
    *y_min = uv_matrix[1];
  }

  xy_matrix[0] = 0;
  xy_matrix[1] = height;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, 3, 3, 1);
  if (uv_matrix[0] > *x_max) {
    *x_max = uv_matrix[0];
  }
  if (uv_matrix[0] < *x_min) {
    *x_min = uv_matrix[0];
  }
  if (uv_matrix[1] > *y_max) {
    *y_max = uv_matrix[1];
  }
  if (uv_matrix[1] < *y_min) {
    *y_min = uv_matrix[1];
  }
}

// Helper function to invert a 3x3 matrix that is in the parameter form.
// target is changed, returns inverse of params in the parameter form.
static void invert_params(double *const params, double *target) {
  double temp[9] = {0};
  params_to_matrix(params, temp);

  // Find determinant of matrix (expansion by minors)
  double det = temp[0] * ((temp[4] * temp[8]) - (temp[5] * temp[7])) -
               temp[1] * ((temp[3] * temp[8]) - (temp[5] * temp[6])) +
               temp[2] * ((temp[3] * temp[7]) - (temp[4] * temp[6]));

  // inverse is transpose of cofactor * 1/det
  double inverse[9] = {0};
  inverse[0] = (temp[4] * temp[8] - temp[7] * temp[5]) / det;
  inverse[1] = (temp[2] * temp[7] - temp[1] * temp[8]) / det;
  inverse[2] = (temp[1] * temp[5] - temp[2] * temp[4]) / det;
  inverse[3] = (temp[5] * temp[6] - temp[3] * temp[8]) / det;
  inverse[4] = (temp[0] * temp[8] - temp[2] * temp[6]) / det;
  inverse[5] = (temp[3] * temp[2] - temp[0] * temp[5]) / det;
  inverse[6] = (temp[3] * temp[7] - temp[6] * temp[4]) / det;
  inverse[7] = (temp[6] * temp[1] - temp[0] * temp[7]) / det;
  inverse[8] = (temp[0] * temp[4] - temp[3] * temp[1]) / det;

  matrix_to_params(inverse, target);
}

int av1_background_sprite(AV1_COMP *cpi, int distance){
  int i;
  int frame;
  YV12_BUFFER_CONFIG *frames[MAX_LAG_BUFFERS] = {NULL};
  int inliers_by_motion[RANSAC_NUM_MOTIONS];
  static const double kIdentityParams[MAX_PARAMDIM - 1] = {0.0, 0.0, 1.0, 0.0,
                                                           0.0, 1.0, 0.0, 0.0};

  // Get frames to be included in background sprite
  frames[0] = cpi->source;
  for (frame = 0; frame < distance; ++frame) {
    struct lookahead_entry *buf = av1_lookahead_peek(cpi->lookahead, frame);
    frames[frame + 1] = &buf->img;
  }

  // Allocate empty arrays for parameters between frames
  double **params = malloc((distance + 1) * sizeof(double *));
  for (i = 0; i < distance + 1; ++i) {
    params[i] = malloc((MAX_PARAMDIM - 1) * sizeof(double));
    memcpy(params[i], kIdentityParams, (MAX_PARAMDIM - 1) * sizeof(double));
  }

  // Use global motion to find affine transformations between frames
  // params[i] will have the transfrom from frame[i] to frame[i-1]
  // params[0] will have the identity matrix because it has no previous frame
  TransformationType model = AFFINE;
  for (frame = 0; frame < distance; ++frame) {
    compute_global_motion_feature_based(model, frames[frame + 1], frames[frame],
#if CONFIG_HIGHBITDEPTH
                                        cpi->common.bit_depth,
#endif  // CONFIG_HIGHBITDEPTH
                                        inliers_by_motion, params[frame + 1],
                                        RANSAC_NUM_MOTIONS);
  }

  // Compound the transformation parameters
  for (i = 1; i < distance + 1; ++i) {
    multiply_params(params[i - 1], params[i], params[i]);
  }

  // Compute frame limits for final stitched images
  double * x_max = malloc((distance + 1) * sizeof(double));
  double * x_min = malloc((distance + 1) * sizeof(double));
  double * y_max = malloc((distance + 1) * sizeof(double));
  double * y_min = malloc((distance + 1) * sizeof(double));
  double pano_x_max = 0;
  double pano_x_min = 0;
  double pano_y_max = 0;
  double pano_y_min = 0;
  for (i = 0; i < distance + 1; ++i) {
    find_limits(cpi->initial_width, cpi->initial_height, params[i], &x_min[i],
                &x_max[i], &y_min[i], &y_max[i]);
    if (x_max[i] > pano_x_max) {
      pano_x_max = x_max[i];
    }
    if (x_min[i] < pano_x_min) {
      pano_x_min = x_min[i];
    }
    if (y_max[i] > pano_y_max) {
      pano_y_max = y_max[i];
    }
    if (y_min[i] > pano_y_min) {
      pano_y_min = y_min[i];
    }
  }

  // Estimate center image based on frame limits
  double pano_center_x = (pano_x_max + pano_x_min) / 2;
  double pano_center_y = (pano_y_max + pano_y_min) / 2;
  double nearest_distance = DBL_MAX;
  int center_idx = -1;
  for (i = 0; i < distance + 1; ++i) {
    double image_center_x = (x_max[i] + x_min[i]) / 2;
    double image_center_y = (y_max[i] + y_min[i]) / 2;
    double distance_from_center = pow(pano_center_x - image_center_x, 2) +
                                  pow(pano_center_y + image_center_y, 2);
    if (distance_from_center < nearest_distance) {
      center_idx = i;
      nearest_distance = distance_from_center;
    }
  }

  // Recompute transformations to adjust to center image
  // Invert center image's transform
  double inverse[8] = {0};
  invert_params(params[center_idx], inverse);

  // Multiply the inverse to all transformation parameters
  for (i = 0; i < distance + 1; ++i) {
    multiply_params(inverse, params[i], params[i]);
  }

  // Recompute frame limits for new adjusted center
  pano_x_max = 0;
  pano_x_min = 0;
  pano_y_max = 0;
  pano_y_min = 0;
  for (i = 0; i < distance + 1; ++i) {
    find_limits(cpi->initial_width, cpi->initial_height, params[i], &x_min[i],
                &x_max[i], &y_min[i], &y_max[i]);
    if (x_max[i] > pano_x_max) {
      pano_x_max = x_max[i];
    }
    if (x_min[i] < pano_x_min) {
      pano_x_min = x_min[i];
    }
    if (y_max[i] > pano_y_max) {
      pano_y_max = y_max[i];
    }
    if (y_min[i] > pano_y_min) {
      pano_y_min = y_min[i];
    }
  }

  // Stitch Images
  //TODO: Allocate panorama and stitch images together

  // Free memory
  for (i = 0; i < distance + 1; ++i) {
    free(params[i]);
    params[i] = NULL;
  }
  free(params);
  free(x_max);
  free(x_min);
  free(y_max);
  free(y_min);
  params = NULL;
  x_max = NULL;
  x_min = NULL;
  y_max = NULL;
  y_min = NULL;

  return distance;
}
