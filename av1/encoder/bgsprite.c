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
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "av1/encoder/bgsprite.h"

#include "aom_mem/aom_mem.h"
#include "av1/common/mv.h"
#include "av1/common/warped_motion.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/global_motion.h"
#include "av1/encoder/mathutils.h"

#define TRANSFORM_MAT_DIM 3

typedef struct {
  // TODO(toddnguyen): Support high bit depth
  uint8_t y;
  uint8_t u;
  uint8_t v;
} yuv_pixel;

// Maps to convert from matrix form to param vector form.
static const int params_to_matrix_map[] = { 2, 3, 0, 4, 5, 1, 6, 7 };
static const int matrix_to_params_map[] = { 2, 5, 0, 1, 3, 4, 6, 7 };

// Convert the parameter array to a 3x3 matrix form.
static void params_to_matrix(double *const params, double *target) {
  for (int i = 0; i < MAX_PARAMDIM - 1; i++) {
    assert(params_to_matrix_map[i] < MAX_PARAMDIM - 1);
    target[i] = params[params_to_matrix_map[i]];
  }
  target[8] = 1;
}

// Convert a 3x3 matrix to a parameter array form.
static void matrix_to_params(double *const matrix, double *target) {
  for (int i = 0; i < MAX_PARAMDIM - 1; i++) {
    assert(matrix_to_params_map[i] < MAX_PARAMDIM - 1);
    target[i] = matrix[matrix_to_params_map[i]];
  }
}

// Do matrix multiplication on params.
static void multiply_params(double *const m1, double *const m2,
                            double *target) {
  double m1_matrix[MAX_PARAMDIM];
  double m2_matrix[MAX_PARAMDIM];
  double result[MAX_PARAMDIM];

  params_to_matrix(m1, m1_matrix);
  params_to_matrix(m2, m2_matrix);
  multiply_mat(m2_matrix, m1_matrix, result, TRANSFORM_MAT_DIM,
               TRANSFORM_MAT_DIM, TRANSFORM_MAT_DIM);
  matrix_to_params(result, target);
}

// Finds x and y limits of a single transformed image.
// Width and height are the size of the input video.
static void find_frame_limit(int width, int height, double *const transform,
                             double *x_min, double *x_max, double *y_min,
                             double *y_max) {
  double transform_matrix[MAX_PARAMDIM];
  double xy_matrix[3] = { 0, 0, 1 };
  double uv_matrix[3] = { 0 };

  params_to_matrix(transform, transform_matrix);
  xy_matrix[0] = 0;
  xy_matrix[1] = 0;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, TRANSFORM_MAT_DIM,
               TRANSFORM_MAT_DIM, 1);
  *x_max = uv_matrix[0];
  *x_min = uv_matrix[0];
  *y_max = uv_matrix[1];
  *y_min = uv_matrix[1];

  xy_matrix[0] = width;
  xy_matrix[1] = 0;
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, TRANSFORM_MAT_DIM,
               TRANSFORM_MAT_DIM, 1);
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
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, TRANSFORM_MAT_DIM,
               TRANSFORM_MAT_DIM, 1);
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
  multiply_mat(transform_matrix, xy_matrix, uv_matrix, TRANSFORM_MAT_DIM,
               TRANSFORM_MAT_DIM, 1);
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

// Finds x and y limits for arrays. Also finds the overall max and minimums
static void find_limits(int width, int height, double **const params,
                        int num_frames, double *x_min, double *x_max,
                        double *y_min, double *y_max, double *pano_x_min,
                        double *pano_x_max, double *pano_y_min,
                        double *pano_y_max) {
  *pano_x_max = DBL_MIN;
  *pano_x_min = DBL_MAX;
  *pano_y_max = DBL_MIN;
  *pano_y_min = DBL_MAX;
  for (int i = 0; i < num_frames; ++i) {
    find_frame_limit(width, height, params[i], &x_min[i], &x_max[i], &y_min[i],
                     &y_max[i]);
    if (x_max[i] > *pano_x_max) {
      *pano_x_max = x_max[i];
    }
    if (x_min[i] < *pano_x_min) {
      *pano_x_min = x_min[i];
    }
    if (y_max[i] > *pano_y_max) {
      *pano_y_max = y_max[i];
    }
    if (y_min[i] < *pano_y_min) {
      *pano_y_min = y_min[i];
    }
  }
}

// Inverts a 3x3 matrix that is in the parameter form.
static void invert_params(double *const params, double *target) {
  double temp[MAX_PARAMDIM] = { 0 };
  params_to_matrix(params, temp);

  // Find determinant of matrix (expansion by minors).
  double det = temp[0] * ((temp[4] * temp[8]) - (temp[5] * temp[7])) -
               temp[1] * ((temp[3] * temp[8]) - (temp[5] * temp[6])) +
               temp[2] * ((temp[3] * temp[7]) - (temp[4] * temp[6]));
  assert(det != 0);

  // inverse is transpose of cofactor * 1/det.
  double inverse[MAX_PARAMDIM] = { 0 };
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

// Swaps two yuv_pixels.
static void swap(yuv_pixel *a, yuv_pixel *b) {
  yuv_pixel temp;
  temp.y = b->y;
  temp.u = b->u;
  temp.v = b->v;

  b->y = a->y;
  b->u = a->u;
  b->v = a->v;

  a->y = temp.y;
  a->u = temp.u;
  a->v = temp.v;
}

// Partitions array to find pivot index in qselect.
static int partition(yuv_pixel arr[], int left, int right, int pivot_idx) {
  yuv_pixel pivot = arr[pivot_idx];

  // Move pivot to the end.
  swap(&arr[pivot_idx], &arr[right]);

  int p_idx = left;
  for (int i = left; i < right; ++i) {
    if (arr[i].y <= pivot.y) {
      swap(&arr[i], &arr[p_idx]);
      p_idx++;
    }
  }

  swap(&arr[p_idx], &arr[right]);

  return p_idx;
}

// Returns the kth element in array, partially sorted in place (quickselect).
static yuv_pixel qselect(yuv_pixel arr[], int left, int right, int k) {
  if (left == right) {
    return arr[left];
  }

  int pivot_idx = left + rand() % (right - left + 1);
  pivot_idx = partition(arr, left, right, pivot_idx);

  if (k == pivot_idx) {
    return arr[k];
  } else if (k < pivot_idx) {
    return qselect(arr, left, pivot_idx - 1, k);
  } else {
    return qselect(arr, pivot_idx + 1, right, k);
  }
}

// Stitches images together to create ARF
static void stitch_images(YV12_BUFFER_CONFIG **const frames, int num_frames,
                          double **const params, YV12_BUFFER_CONFIG *panorama,
                          double *const x_min, double *const x_max,
                          double *const y_min, double *const y_max,
                          double pano_x_min, double pano_x_max,
                          double pano_y_min, double pano_y_max) {
  int width = ceil(pano_x_max) - floor(pano_x_min) + 1;
  int height = ceil(pano_y_max) - floor(pano_y_min) + 1;

  // Create temp_pano[y][x][num_frames] stack of pixel values
  yuv_pixel ***temp_pano = aom_malloc(height * sizeof(*temp_pano));
  for (int i = 0; i < height; ++i) {
    temp_pano[i] = aom_malloc(width * sizeof(**temp_pano));
    for (int j = 0; j < width; ++j) {
      temp_pano[i][j] = aom_malloc(num_frames * sizeof(***temp_pano));
    }
  }
  // Create count[y][x] to count how many values in stack for median filtering
  int **count = aom_malloc(height * sizeof(*count));
  for (int i = 0; i < height; ++i) {
    count[i] = aom_calloc(width, sizeof(**count));  // counts initialized to 0
  }

  // Re-sample images onto panorama (pre-median filtering).
  int x_offset = -floor(pano_x_min);
  int y_offset = -floor(pano_y_min);
  int frame_width = frames[0]->y_width;
  int frame_height = frames[0]->y_height;
  for (int i = 0; i < num_frames; ++i) {
    // Find transforms from panorama coordinate system back to single image
    // coordinate system for sampling.
    int transformed_width = ceil(x_max[i]) - floor(x_min[i]) + 1;
    int transformed_height = ceil(y_max[i]) - floor(y_min[i]) + 1;

    double transform_matrix[MAX_PARAMDIM];
    double transform_params[MAX_PARAMDIM - 1];
    invert_params(params[i], transform_params);
    params_to_matrix(transform_params, transform_matrix);

    for (int y = 0; y < transformed_height; ++y) {
      for (int x = 0; x < transformed_width; ++x) {
        // Do transform.
        double xy_matrix[3] = { x + floor(x_min[i]), y + floor(y_min[i]), 1 };
        double uv_matrix[3] = { 0 };
        multiply_mat(transform_matrix, xy_matrix, uv_matrix, TRANSFORM_MAT_DIM,
                     TRANSFORM_MAT_DIM, 1);
        int image_x = round(uv_matrix[0]);
        int image_y = round(uv_matrix[1]);

        // Check if valid point in original image.
        if (image_x >= 0 && image_x < frame_width && image_y >= 0 &&
            image_y < frame_height) {
          // Place in panorama stack.
          int pano_x = x + floor(x_min[i]) + x_offset;
          int pano_y = y + floor(y_min[i]) + y_offset;

          int ychannel_idx = image_y * frames[i]->y_stride + image_x;
          int uvchannel_idx =
              (image_y >> frames[i]->subsampling_y) * frames[i]->uv_stride +
              (image_x >> frames[i]->subsampling_x);

          temp_pano[pano_y][pano_x][count[pano_y][pano_x]].y =
              frames[i]->y_buffer[ychannel_idx];
          temp_pano[pano_y][pano_x][count[pano_y][pano_x]].u =
              frames[i]->u_buffer[uvchannel_idx];
          temp_pano[pano_y][pano_x][count[pano_y][pano_x]].v =
              frames[i]->v_buffer[uvchannel_idx];
          // Update count.
          count[pano_y][pano_x]++;
        }
      }
    }
  }

  // Apply median filtering using quickselect.
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      if (count[y][x] == 0) {
        // Just make the pixel black.
        // TODO(toddnguyen): Color the pixel with nearest neighbor
      } else {
        // Find
        int median_idx = floor(count[y][x] / 2);
        yuv_pixel median =
            qselect(temp_pano[y][x], 0, count[y][x] - 1, median_idx);

        // Make the median value the 0th index for UV subsampling later
        swap(&temp_pano[y][x][median_idx], &temp_pano[y][x][0]);
        assert(median.y == temp_pano[y][x][0].y &&
               median.u == temp_pano[y][x][0].u &&
               median.v == temp_pano[y][x][0].v);

        // TODO(toddnguyen): Get rid of this (supress compiler warnings)
        (void)median;
      }
    }
  }

  // UV subsampling

  // TODO(toddnguyen): Get rid of this (supress compiler warnings)
  (void)panorama;

  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      aom_free(temp_pano[i][j]);
    }
    aom_free(temp_pano[i]);
    aom_free(count[i]);
  }
  aom_free(count);
  aom_free(temp_pano);
}

int av1_background_sprite(AV1_COMP *cpi, int distance) {
  int i;
  int frame;
  YV12_BUFFER_CONFIG *frames[MAX_LAG_BUFFERS] = { NULL };
  int inliers_by_motion[RANSAC_NUM_MOTIONS];
  static const double identity_params[MAX_PARAMDIM - 1] = {
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
  };

  // Get frames to be included in background sprite.
  frames[0] = cpi->source;
  for (frame = 0; frame < distance; ++frame) {
    struct lookahead_entry *buf = av1_lookahead_peek(cpi->lookahead, frame);
    frames[frame + 1] = &buf->img;
  }

  // Allocate empty arrays for parameters between frames.
  double **params = aom_malloc((distance + 1) * sizeof(*params));
  for (i = 0; i < distance + 1; ++i) {
    params[i] = aom_malloc(sizeof(identity_params));
    memcpy(params[i], identity_params, sizeof(identity_params));
  }

  // Use global motion to find affine transformations between frames.
  // params[i] will have the transform from frame[i] to frame[i-1].
  // params[0] will have the identity matrix because it has no previous frame.
  TransformationType model = AFFINE;
  for (frame = 0; frame < distance; ++frame) {
    int global_motion_ret = compute_global_motion_feature_based(
        model, frames[frame + 1], frames[frame],
#if CONFIG_HIGHBITDEPTH
        cpi->common.bit_depth,
#endif  // CONFIG_HIGHBITDEPTH
        inliers_by_motion, params[frame + 1], RANSAC_NUM_MOTIONS);

    // Quit if global motion had an error.
    if (global_motion_ret == 0) {
      for (i = 0; i < distance + 1; ++i) {
        aom_free(params[i]);
      }
      aom_free(params);
      return EXIT_FAILURE;
    }
  }

  // Compound the transformation parameters.
  for (i = 1; i < distance + 1; ++i) {
    multiply_params(params[i - 1], params[i], params[i]);
  }

  // Compute frame limits for final stitched images.
  double pano_x_max = DBL_MIN;
  double pano_x_min = DBL_MAX;
  double pano_y_max = DBL_MIN;
  double pano_y_min = DBL_MAX;
  double *x_max = aom_malloc((distance + 1) * sizeof(*x_max));
  double *x_min = aom_malloc((distance + 1) * sizeof(*x_min));
  double *y_max = aom_malloc((distance + 1) * sizeof(*y_max));
  double *y_min = aom_malloc((distance + 1) * sizeof(*y_min));

  find_limits(cpi->initial_width, cpi->initial_height, params, distance + 1,
              x_min, x_max, y_min, y_max, &pano_x_min, &pano_x_max, &pano_y_min,
              &pano_y_max);

  // Estimate center image based on frame limits.
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
  assert(center_idx != -1);

  // Recompute transformations to adjust to center image.
  // Invert center image's transform.
  double inverse[MAX_PARAMDIM - 1] = { 0 };
  invert_params(params[center_idx], inverse);

  // Multiply the inverse to all transformation parameters.
  for (i = 0; i < distance + 1; ++i) {
    multiply_params(inverse, params[i], params[i]);
  }

  // Recompute frame limits for new adjusted center.
  find_limits(cpi->initial_width, cpi->initial_height, params, distance + 1,
              x_min, x_max, y_min, y_max, &pano_x_min, &pano_x_max, &pano_y_min,
              &pano_y_max);

  // Stitch Images.
  stitch_images(frames, distance + 1, params, &cpi->alt_ref_buffer, x_min,
                x_max, y_min, y_max, pano_x_min, pano_x_max, pano_y_min,
                pano_y_max);

  // Free memory.
  for (i = 0; i < distance + 1; ++i) {
    aom_free(params[i]);
  }
  aom_free(params);
  aom_free(x_max);
  aom_free(x_min);
  aom_free(y_max);
  aom_free(y_min);

  return EXIT_SUCCESS;
}
