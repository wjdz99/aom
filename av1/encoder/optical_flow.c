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
#if DCONFIG_OPTICAL_FLOW_API

#include <math.h>
#include <limits.h>

#include "config/aom_config.h"
#include "av1/common/av1_common_int.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/mathutils.h"
#include "av1/encoder/optical_flow.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/temporal_filter.h"
#include "aom_mem/aom_mem.h"

// Helper function to determine whether a frame is encoded with high bit-depth.
static INLINE int is_frame_high_bitdepth(const YV12_BUFFER_CONFIG *frame) {
  return (frame->flags & YV12_FLAG_HIGHBITDEPTH) ? 1 : 0;
}

typedef struct LOCALMV {
  double row;
  double col;
} LOCALMV;

// coefficients for bilinear interpolation on unit square
int pixel_interp(const double x, const double y, const int b00, const int b01,
                 const int b10, const int b11) {
  int xint = (int)x;
  int yint = (int)y;
  double xdec = x - xint;
  double ydec = y - yint;
  double a = (1 - xdec) * (1 - ydec);
  double b = xdec * (1 - ydec);
  double c = (1 - xdec) * ydec;
  double d = xdec * ydec;
  // if x, y are already integers, this results to b00
  int interp = a * b00 + b * b01 + c * b10 + d * b11;
  return interp;
}
// bilinear interpolation to find subpixel values
int get_subpixels(const YV12_BUFFER_CONFIG *frame, int *pred, const int w,
                  const int h, LOCALMV mv, const double x_coord,
                  const double y_coord) {
  double left = x_coord + mv.row;
  double top = y_coord + mv.col;
  int fromedge = 2;
  int height = frame->y_crop_height;
  int width = frame->y_crop_width;
  if (left < 1) left = 1;
  if (top < 1) top = 1;
  // could use elements past boundary where stride > width
  if (top > height - fromedge) top = height - fromedge;
  if (left > width - fromedge) left = width - fromedge;
  uint8_t *buf = frame->y_buffer;
  int s = frame->y_stride;
  int prev = -1;

  int xint;
  int yint;
  int idx = 0;
  for (int y = prev; y < prev + h; y++) {
    for (int x = prev; x < prev + w; x++) {
      double xx = left + x;
      double yy = top + y;
      xint = (int)xx;
      yint = (int)yy;
      int interp = pixel_interp(
          xx, yy, buf[yint * s + xint], buf[yint * s + (xint + 1)],
          buf[(yint + 1) * s + xint], buf[(yint + 1) * s + (xint + 1)]);
      pred[idx++] = interp;
    }
  }
  return 0;
}
// Scharr filter to compute spatial gradient
void spatial_gradient(const YV12_BUFFER_CONFIG *frame, const int x_coord,
                      const int y_coord, const int direction,
                      double *derivative) {
  double filter[9];
  // Scharr filters
  double gx[9] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
  double gy[9] = { -3, -10, -3, 0, 0, 0, 3, 10, 3 };
  if (direction == 0) {  // x direction
    memcpy(filter, gx, sizeof(filter));
  } else {  // y direction
    memcpy(filter, gy, sizeof(filter));
  }
  int idx = 0;
  double d = 0;
  for (int yy = -1; yy <= 1; yy++) {
    for (int xx = -1; xx <= 1; xx++) {
      d += filter[idx] *
           frame->y_buffer[(y_coord + yy) * frame->y_stride + (x_coord + xx)];
      idx++;
    }
  }
  *derivative = d;
}
// Determine the spatial gradient at subpixel locations
// For example, when reducing images for pyramidal LK,
// corners found in original image may be at subpixel locations.
void gradient_interp(int *fullpel_deriv, const double x_coord,
                     const double y_coord, const int w, const int h,
                     double *derivative) {
  int xint = (int)x_coord;
  int yint = (int)y_coord;
  int interp;
  if (xint + 1 > w - 1 || yint + 1 > h - 1) {
    interp = fullpel_deriv[yint * w + xint];
  } else {
    interp = pixel_interp(x_coord, y_coord, fullpel_deriv[yint * w + xint],
                          fullpel_deriv[yint * w + (xint + 1)],
                          fullpel_deriv[(yint + 1) * w + xint],
                          fullpel_deriv[(yint + 1) * w + (xint + 1)]);
  }

  *derivative = interp;
}

void temporal_gradient(const YV12_BUFFER_CONFIG *frame,
                       const YV12_BUFFER_CONFIG *frame2, const double x_coord,
                       const double y_coord, const int bit_depth,
                       double *derivative, LOCALMV *mv) {
  int idx = 0;
  const int w = 5;
  const int h = 5;
  uint8_t pred1[25];
  uint8_t pred2[25];

  const int y = y_coord;
  const int x = x_coord;
  double ydec = y_coord - y;
  double xdec = x_coord - x;
  const int is_intrabc = 0;  // Is intra-copied?
  const int is_high_bitdepth = is_frame_high_bitdepth(frame2);
  const int subsampling_x = 0, subsampling_y = 0;  // for y-buffer
  const int_interpfilters interp_filters =
      av1_broadcast_interp_filter(MULTITAP_SHARP);
  const int plane = 0;  // y-plane
  const struct buf_2d ref_buf2 = { NULL, frame2->y_buffer, frame2->y_crop_width,
                                   frame2->y_crop_height, frame2->y_stride };
  struct scale_factors scale;
  av1_setup_scale_factors_for_frame(&scale, frame->y_crop_width,
                                    frame->y_crop_height, frame->y_crop_width,
                                    frame->y_crop_height);
  InterPredParams inter_pred_params;
  av1_init_inter_params(&inter_pred_params, w, h, y, x, subsampling_x,
                        subsampling_y, bit_depth, is_high_bitdepth, is_intrabc,
                        &scale, &ref_buf2, interp_filters);
  inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
  MV newmv = { .row = round((mv->row + xdec) * 8),
               .col = round((mv->col + ydec) * 8) };
  av1_enc_build_one_inter_predictor(pred2, w, &newmv, &inter_pred_params);
  pred2[3] = pred2[0];
  const struct buf_2d ref_buf1 = { NULL, frame->y_buffer, frame->y_crop_width,
                                   frame->y_crop_height, frame->y_stride };
  av1_init_inter_params(&inter_pred_params, w, h, y, x, subsampling_x,
                        subsampling_y, bit_depth, is_high_bitdepth, is_intrabc,
                        &scale, &ref_buf1, interp_filters);
  inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
  MV zeroMV = { .row = round(xdec * 8), .col = round(ydec * 8) };
  av1_enc_build_one_inter_predictor(pred1, w, &zeroMV, &inter_pred_params);
  pred1[3] = pred1[0];
  double d1 = 0, d2 = 0;
  d1 = pred1[3];
  d2 = pred2[3];
  *derivative = d2 - d1;
}
// Numerical differentiate over window_size x window_size surrounding (x,y)
// location. Alters ix, iy, it to contain numerical partial derivatives
// TODO(lpartin) if LK is computed over all pixels instead of on corners,
// it would be more efficient to compute gradients on entire image (instead of
// window)
void gradients_over_window(const YV12_BUFFER_CONFIG *frame,
                           const YV12_BUFFER_CONFIG *ref_frame,
                           const double x_coord, const double y_coord,
                           const int window_size, const int bit_depth,
                           double *ix, double *iy, double *it, LOCALMV *mv) {
  const double left = x_coord - window_size / 2;
  const double top = y_coord - window_size / 2;
  // gradient operators need pixel before and after (start at 1)
  const double x_start = AOMMAX(1, left);
  const double y_start = AOMMAX(1, top);
  const int frame_height = frame->y_crop_height;
  const int frame_width = frame->y_crop_width;
  const int ystride = frame->y_stride;
  const int border = frame->border;
  double deriv_x;
  double deriv_y;
  double deriv_t;

  const double x_end = AOMMIN(x_coord + window_size / 2, frame_width - 2);
  const double y_end = AOMMIN(y_coord + window_size / 2, frame_height - 2);
  int xs = AOMMAX(1, x_start - 1);
  int ys = AOMMAX(1, y_start - 1);
  int xe = AOMMIN(x_end + 2, frame_width - 2);
  int ye = AOMMIN(y_end + 2, frame_height - 2);
  // assuming that derivatives are integers (which is the case
  // for both sobel and scharr filters on greyscale int values)
  int *fullpel_dx = malloc((ye - ys) * (xe - xs) * sizeof(int));
  int *fullpel_dy = malloc((ye - ys) * (xe - xs) * sizeof(int));
  // TODO: This could be more efficient in the case that x_coord
  // and y_coord are integers.. but it may look more messy.
  for (int j = ys; j < ye; j++) {
    for (int i = xs; i < xe; i++) {
      spatial_gradient(frame, i, j, 0, &deriv_x);
      spatial_gradient(frame, i, j, 1, &deriv_y);
      int idx = (j - ys) * (xe - xs) + (i - xs);
      fullpel_dx[idx] = (int)deriv_x;
      fullpel_dy[idx] = (int)deriv_y;
    }
  }
  // compute numerical differentiation for every pixel in window
  for (double j = y_start; j < y_end; j++) {
    for (double i = x_start; i < x_end; i++) {
      temporal_gradient(frame, ref_frame, i, j, bit_depth, &deriv_t, mv);
      gradient_interp(fullpel_dx, i - xs, j - ys, xe - xs, ye - ys, &deriv_x);
      gradient_interp(fullpel_dy, i - xs, j - ys, xe - xs, ye - ys, &deriv_y);
      int idx = (int)(j - top) * window_size + (int)(i - left);
      // TODO: move scaling into spatial gradient function
      // keeping here for now, since fullpel array type is int
      ix[idx] = deriv_x;  // / 32.0; //scaling for scharr filter
      iy[idx] = deriv_y;  // / 32.0;
      it[idx] = deriv_t;
    }
  }
  // TODO: to avoid setting deriv arrays to zero for every iteration,
  // could instead pass these two values back through function call
  // int first_idx = (int)(y_start - top) * window_size + (int)(x_start - left);
  // int width = window_size - ((int)(x_start - left) + (int)(left + window_size
  // - x_end));

  free(fullpel_dx);
  free(fullpel_dy);
}

// To compute eigenvalues of 2x2 matrix: Solve for lambda where
// Determinant of matrix - lambda*identity == 0
void eigenvalues_2x2(const double *matrix, double *eig) {
  double a = 1;
  double b = -1 * matrix[0] - matrix[3];
  double c = -1 * matrix[1] * matrix[2] + matrix[0] * matrix[3];
  // quadratic formula
  double discriminant = b * b - 4 * a * c;
  eig[0] = (-b - sqrt(discriminant)) / (2.0 * a);
  eig[1] = (-b + sqrt(discriminant)) / (2.0 * a);
  // double check that eigenvalues are ordered by magnitude
  if (fabs(eig[0]) > fabs(eig[1])) {
    double tmp = eig[0];
    eig[0] = eig[1];
    eig[1] = tmp;
  }
}
// forward declaration for corner score
void gradients_over_window(const YV12_BUFFER_CONFIG *frame,
                           const YV12_BUFFER_CONFIG *ref_frame,
                           const double x_coord, const double y_coord,
                           const int window_size, const int bit_depth,
                           double *ix, double *iy, double *it, LOCALMV *mv);
// Shi-Tomasi corner detection criteria
double corner_score(const YV12_BUFFER_CONFIG *frame_to_filter,
                    const YV12_BUFFER_CONFIG *ref_frame, const int x,
                    const int y, double *i_x, double *i_y, double *i_t,
                    const int n, const int bit_depth) {
  double eig[2];
  LOCALMV mv = { .row = 0, .col = 0 };
  // TODO: technically, ref_frame and i_t are not used by corner score
  // so these could be replaced by dummy variables,
  // or change this to spatial gradient function over window only
  gradients_over_window(frame_to_filter, ref_frame, x, y, n, bit_depth, i_x,
                        i_y, i_t, &mv);
  double Mres1[1] = { 0 }, Mres2[1] = { 0 }, Mres3[1] = { 0 };
  multiply_mat(i_x, i_x, Mres1, 1, n * n, 1);
  multiply_mat(i_x, i_y, Mres2, 1, n * n, 1);
  multiply_mat(i_y, i_y, Mres3, 1, n * n, 1);
  double M[4] = { Mres1[0], Mres2[0], Mres2[0], Mres3[0] };
  eigenvalues_2x2(M, eig);
  return fabs(eig[0]);
}
// Finds corners in frame_to_filter
// For less strict requirements (i.e. more corners), decrease threshold
int detect_corners(const YV12_BUFFER_CONFIG *frame_to_filter,
                   const YV12_BUFFER_CONFIG *ref_frame, const int maxcorners,
                   int *ref_corners, const int bit_depth) {
  // TODO: currently if maxcorners is decreased, then it only means
  // corners will be omited from bottom-right of image. if maxcorners
  // is actually used, then this algorithm would need to re-iterate
  // and choose threshold based on that
  int frame_height = frame_to_filter->y_crop_height;
  int frame_width = frame_to_filter->y_crop_width;
  int countcorners = 0;
  int fromedge = 10;
  double max_score;
  double threshold = 0.15;
  double score;
  int n = 7;
  double i_x[49];
  double i_y[49];
  double i_t[49];
  // rough estimate of max corner score in image
  for (int x = fromedge; x < frame_width - fromedge; x += 1) {
    for (int y = fromedge; y < frame_height - fromedge; y += frame_height / 5) {
      for (int i = 0; i < 49; i++) {
        i_x[i] = 0;
        i_y[i] = 0;
        i_t[i] = 0;
      }
      score = corner_score(frame_to_filter, ref_frame, x, y, i_x, i_y, i_t, n,
                           bit_depth);
      if (x == fromedge && y == fromedge) {
        max_score = score;
      }
      if (score > max_score) {
        max_score = score;
      }
    }
  }
  // score all the points and choose corners over threshold
  for (int x = fromedge; x < frame_width - fromedge; x += 1) {
    for (int y = fromedge;
         (y < frame_height - fromedge) && countcorners < maxcorners; y += 1) {
      for (int i = 0; i < 49; i++) {
        i_x[i] = 0;
        i_y[i] = 0;
        i_t[i] = 0;
      }
      score = corner_score(frame_to_filter, ref_frame, x, y, i_x, i_y, i_t, n,
                           bit_depth);
      if (score > threshold * max_score) {
        ref_corners[countcorners * 2] = x;
        ref_corners[countcorners * 2 + 1] = y;
        countcorners++;
      }
    }
  }
  return countcorners;
}
void gaussian(const double sigma, const int n, const int normalize,
              double *weights) {
  double total_weight = 0;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      double distance = sqrt(pow(n / 2 - i, 2) + pow(n / 2 - j, 2));
      double weight = exp(-0.5 * pow(distance / sigma, 2));
      weights[j * n + i] = weight;
      total_weight += weight;
    }
  }
  if (normalize == 1) {
    for (int j = 0; j < n; j++) {
      weights[j] = weights[j] / total_weight;
    }
  }
}

double convolve(double *filter, int *img, int size) {
  double result = 0;
  for (int i = 0; i < size; i++) {
    result += filter[i] * img[i];
  }
  return result;
}
// Applies a Gaussian low-pass smoothing filter to produce
// a corresponding lower resolution image with halved dimensions
void reduce(uint8_t *img, int height, int width, int stride,
            uint8_t *reduced_img) {
  int new_width = width / 2;
  int new_height = height / 2;
  int window_size = 5;
  double gaussian_filter[25] = {
    1. / 256, 1.0 / 64, 3. / 128, 1. / 64,  1. / 256, 1. / 64, 1. / 16,
    3. / 32,  1. / 16,  1. / 64,  3. / 128, 3. / 32,  9. / 64, 3. / 32,
    3. / 128, 1. / 64,  1. / 16,  3. / 32,  1. / 16,  1. / 64, 1. / 256,
    1. / 64,  3. / 128, 1. / 64,  1. / 256
  };
  // filter is 5x5 so need prev and forward 2 pixels
  int img_section[25];
  for (int y = 0; y < height - 1; y += 2) {
    for (int x = 0; x < width - 1; x += 2) {
      int i = 0;
      for (int yy = y - window_size / 2; yy <= y + window_size / 2; yy++) {
        for (int xx = x - window_size / 2; xx <= x + window_size / 2; xx++) {
          int yvalue = yy;
          int xvalue = xx;
          // copied pixels outside the boundary
          if (yvalue < 0) yvalue = 0;
          if (xvalue < 0) xvalue = 0;
          if (yvalue >= height) yvalue = height - 1;
          if (xvalue >= width) xvalue = width - 1;
          img_section[i++] = img[yvalue * stride + xvalue];
        }
      }
      reduced_img[(y / 2) * new_width + (x / 2)] =
          (uint8_t)convolve(gaussian_filter, img_section, pow(window_size, 2));
    }
  }
}
int cmpfunc(const void *a, const void *b) { return (*(int *)a - *(int *)b); }
void filter_mvs(const mv_filter_type mv_filter, LOCALMV *localmvs, MV *mvs) {
  // for smoothing filter
  double gaussian_filter[25] = {
    1. / 256, 1.0 / 64, 3. / 128, 1. / 64,  1. / 256, 1. / 64, 1. / 16,
    3. / 32,  1. / 16,  1. / 64,  3. / 128, 3. / 32,  9. / 64, 3. / 32,
    3. / 128, 1. / 64,  1. / 16,  3. / 32,  1. / 16,  1. / 64, 1. / 256,
    1. / 64,  3. / 128, 1. / 64,  1. / 256
  };
  // for median filter
  int mvrows[25];
  int mvcols[25];
  if (mv_filter != None) {
    for (int y = 0; y < frame_height; y++) {
      for (int x = 0; x < frame_width; x++) {
        if (fabs(localmvs[y * frame_width + x].row) > 0 ||
            fabs(localmvs[y * frame_width + x].col) > 0) {
          int i = 0;
          double filtered_row = 0;
          double filtered_col = 0;
          // TODO: might be worth applying filter only to surrounding
          // mvs within the pixel's block, since the next block's mvs
          // may change drastically
          for (int yy = y - 5 / 2; yy <= y + 5 / 2; yy++) {
            for (int xx = x - 5 / 2; xx <= x + 5 / 2; xx++) {
              int yvalue = yy + y;
              int xvalue = xx + x;
              // copied pixels outside the boundary
              if (yvalue < 0) yvalue = 0;
              if (xvalue < 0) xvalue = 0;
              if (yvalue >= frame_height) yvalue = frame_height - 1;
              if (xvalue >= frame_width) xvalue = frame_width - 1;
              int index = yvalue * frame_width + xvalue;
              if (mv_filter == Smooth) {
                filtered_row += mvs[index].row * gaussian_filter[i];
                filtered_col += mvs[index].col * gaussian_filter[i];
              } else if (mv_filter == Median) {
                mvrows[i] = mvs[index].row;
                mvcols[i] = mvs[index].col;
              }
              i++;
            }
          }
          MV mv;
          if (mv_filter == Smooth) {
            mv.row = filtered_row;
            mv.col = filtered_col;
          } else if (mv_filter == Median) {
            qsort(mvrows, 25, sizeof(int), cmpfunc);
            qsort(mvcols, 25, sizeof(int), cmpfunc);
            mv.row = mvrows[25 / 2];
            mv.col = mvcols[25 / 2];
          }
          LOCALMV localmv = { .row = mv.row / 8.0, .col = mv.row / 8.0 };
          localmvs[y * frame_width + x] = localmv;
          // if mvs array is immediately updated here, then the result may
          // propagate to other pixels.
        }
      }
    }
    for (int i = 0; i < frame_height * frame_width; i++) {
      if (fabs(localmvs[i].row) > 0 || fabs(localmvs[i].col) > 0) {
        MV mv = { .row = round(8 * localmvs[i].row),
                  .col = round(8 * localmvs[i].col) };
        mvs[i] = mv;
      }
    }
  }
}
// Computes optical flow by applying algorithm at
// multiple pyramid levels of images (lower-resolution, smoothed images)
// This accounts for larger motions.
// Inputs:
//   from_frame Frame buffer.
//   to_frame: Frame buffer. MVs point from_frame -> to_frame.
//   from_frame_idx: Index of from_frame.
//   to_frame_idx: Index of to_frame. Return all zero MVs when idx are equal.
//   bit_depth:
//   levels: total pyramid levels. Must be in [1,5].
//   window_size: as used by specific algorithm.
//   mv_filter: None, Smooth, or Median.
//   method: LucasKanade,
//   mvs: pointer to MVs. Contains initialization, and modified
//   based on optical flow.
void optical_flow(const YV12_BUFFER_CONFIG *from_frame,
                  const YV12_BUFFER_CONFIG *to_frame, const int from_frame_idx,
                  const int to_frame_idx, const int bit_depth, int levels,
                  int window_size, const mv_filter_type mv_filter,
                  const optflow_method method, MV *mvs) {
  if (levels < 1 || levels > 5) {
    printf("Pyramid levels out of bounds. Choose a value within [%d, %d].\n", 1,
           5);
    printf("Resetting to default value %d\n", OPFL_PYRAMID_LEVELS);
    levels = OPFL_PYRAMID_LEVELS;
  }
  if (window_size < 1 || window_size % 2 != 1) {
    printf("Window size is not positive or not even.\n");
    printf("Resetting to default value %d\n", OPFL_WINDOW_SIZE);
    window_size = OPFL_WINDOW_SIZE;
  }
  const int frame_height = from_frame->y_crop_height;
  const int frame_width = from_frame->y_crop_width;
  if ((frame_height / pow(2.0, levels - 1) < 50 ||
       frame_height / pow(2.0, levels - 1) < 50) &&
      levels > 1)
    levels = levels - 1;
  LOCALMV *localmvs = malloc(frame_height * frame_width * sizeof(LOCALMV));
  if (from_frame_idx == to_frame_idx) {
    // immediately return all zero mvs when frame indices are equal
    for (int yy = 0; yy < frame_height; yy++) {
      for (int xx = 0; xx < frame_width; xx++) {
        MV mv = { .row = 0, .col = 0 };
        mvs[yy * frame_width + xx] = mv;
      }
    }
    return;
  }

  // Initialize double mvs based on input parameter mvs array
  for (int i = 0; i < frame_width * frame_height; i++) {
    MV mv = mvs[i];
    LOCALMV localmv = { .row = mv.row / 8.0, .col = mv.col / 8.0 };
    localmvs[i] = localmv;
  }
  // Apply optical flow algorithm

  // Update original mvs array
  for (int j = 0; j < frame_height; j++) {
    for (int i = 0; i < frame_width; i++) {
      int idx = j * frame_width + i;
      int new_x = localmvs[idx].row + i;
      int new_y = localmvs[idx].col + j;
      if ((fabs(localmvs[idx].row) >= 0.125 ||
           fabs(localmvs[idx].col) >= 0.125)) {
        // if mv points outside of frame (lost feature), keep old mv.
        if (new_x < frame_width && new_x >= 0 && new_y < frame_height &&
            new_y >= 0) {
          MV mv = { .row = round(8 * localmvs[idx].row),
                    .col = round(8 * localmvs[idx].col) };
          mvs[idx] = mv;
        }
      }
    }
  }

  filter_mvs(mv_filter, localmvs, mvs);

  for (int i = 1; i < levels; i++) {
    free(images1[i]);
    free(images2[i]);
  }
  free(localmvs);
}
#endif
