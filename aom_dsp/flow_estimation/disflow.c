/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 3-Clause Clear License
 * and the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear
 * License was not distributed with this source code in the LICENSE file, you
 * can obtain it at aomedia.org/license/software-license/bsd-3-c-c/.  If the
 * Alliance for Open Media Patent License 1.0 was not distributed with this
 * source code in the PATENTS file, you can obtain it at
 * aomedia.org/license/patent-license/.
 */

// Dense Inverse Search flow algorithm
// Paper: https://arxiv.org/abs/1603.03590

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/flow_estimation/disflow.h"
#include "aom_dsp/flow_estimation/corner_detect.h"
#include "aom_dsp/flow_estimation/pyramid.h"
#include "aom_dsp/flow_estimation/ransac.h"
#include "aom_mem/aom_mem.h"

#include "config/aom_dsp_rtcd.h"

// TODO(rachelbarker):
// Implement specialized functions for upscaling flow fields,
// replacing av1_upscale_plane_double_prec().
// Then we can avoid needing to include code from av1/
#include "av1/common/resize.h"

#include <assert.h>

// Amount to downsample the flow field by.
// eg. DOWNSAMPLE_SHIFT = 2 (DOWNSAMPLE_FACTOR == 4) means we calculate
// one flow point for each 4x4 pixel region of the frame
// Must be a power of 2
#define DOWNSAMPLE_SHIFT 3
#define DOWNSAMPLE_FACTOR (1 << DOWNSAMPLE_SHIFT)
// Number of outermost flow field entries (on each edge) which can't be
// computed, because the patch they correspond to extends outside of the
// frame
// The border is (DISFLOW_PATCH_SIZE >> 1) pixels, which is
// (DISFLOW_PATCH_SIZE >> 1) >> DOWNSAMPLE_SHIFT many flow field entries
#define FLOW_BORDER ((DISFLOW_PATCH_SIZE >> 1) >> DOWNSAMPLE_SHIFT)
// When downsampling the flow field, each flow field entry covers a square
// region of pixels in the image pyramid. This value is equal to the position
// of the center of that region, as an offset from the top/left edge.
//
// Note: Using ((DOWNSAMPLE_FACTOR - 1) / 2) is equivalent to the more
// natural expression ((DOWNSAMPLE_FACTOR / 2) - 1),
// unless DOWNSAMPLE_FACTOR == 1 (ie, no downsampling), in which case
// this gives the correct offset of 0 instead of -1.
#define UPSAMPLE_CENTER_OFFSET ((DOWNSAMPLE_FACTOR - 1) / 2)

// Internal precision of cubic interpolation filters
// The limiting factor here is that:
// * Before integerizing, the maximum value of any kernel tap is 1.0
// * After integerizing, each tap must fit into an int16_t.
// Thus the largest multiplier we can get away with is 2^14 = 16384,
// as 2^15 = 32768 is too large to fit in an int16_t.
#define CUBIC_PREC_BITS 14

static INLINE void getCubicKernel_dbl(double x, double *kernel) {
  assert(0 <= x && x < 1);
  double x2 = x * x;
  double x3 = x2 * x;
  kernel[0] = -0.5 * x + x2 - 0.5 * x3;
  kernel[1] = 1.0 - 2.5 * x2 + 1.5 * x3;
  kernel[2] = 0.5 * x + 2.0 * x2 - 1.5 * x3;
  kernel[3] = -0.5 * x2 + 0.5 * x3;
}

static INLINE void getCubicKernel_int(double x, int *kernel) {
  assert(0 <= x && x < 1);
  double x2 = x * x;
  double x3 = x2 * x;

  double k0 = -0.5 * x + x2 - 0.5 * x3;
  double k1 = 1.0 - 2.5 * x2 + 1.5 * x3;
  double k2 = 0.5 * x + 2.0 * x2 - 1.5 * x3;
  double k3 = -0.5 * x2 + 0.5 * x3;

  kernel[0] = (int)rint(k0 * (1 << CUBIC_PREC_BITS));
  kernel[1] = (int)rint(k1 * (1 << CUBIC_PREC_BITS));
  kernel[2] = (int)rint(k2 * (1 << CUBIC_PREC_BITS));
  kernel[3] = (int)rint(k3 * (1 << CUBIC_PREC_BITS));
}

static INLINE double getCubicValue_dbl(const double *p, const double *kernel) {
  return kernel[0] * p[0] + kernel[1] * p[1] + kernel[2] * p[2] +
         kernel[3] * p[3];
}

static INLINE int getCubicValue_int(const int *p, const int *kernel) {
  return kernel[0] * p[0] + kernel[1] * p[1] + kernel[2] * p[2] +
         kernel[3] * p[3];
}

static INLINE double bicubic_interp_one(const double *arr, int stride,
                                        double *h_kernel, double *v_kernel) {
  double tmp[1 * 4];

  // Horizontal convolution
  for (int i = -1; i < 3; ++i) {
    tmp[i + 1] = getCubicValue_dbl(&arr[i * stride - 1], h_kernel);
  }

  // Vertical convolution
  return getCubicValue_dbl(tmp, v_kernel);
}

static int determine_disflow_correspondence(int *frm_corners,
                                            int num_frm_corners,
                                            const FlowField *flow,
                                            Correspondence *correspondences) {
  int width = flow->width;
  int height = flow->height;
  int stride = flow->stride;

  int num_correspondences = 0;
  for (int i = 0; i < num_frm_corners; ++i) {
    int x0 = frm_corners[2 * i];
    int y0 = frm_corners[2 * i + 1];

    // Offset points, to compensate for the fact that (say) a flow field entry
    // at horizontal index i, is nominally associated with the pixel at
    // horizontal coordinate (i << DOWNSAMPLE_FACTOR) + UPSAMPLE_CENTER_OFFSET
    // This offset must be applied before we split the coordinate into integer
    // and fractional parts, in order for the interpolation to be correct.
    int x = x0 - UPSAMPLE_CENTER_OFFSET;
    int y = y0 - UPSAMPLE_CENTER_OFFSET;

    // Split the pixel coordinates into integer flow field coordinates and
    // an offset for interpolation
    int flow_x = x >> DOWNSAMPLE_SHIFT;
    double flow_sub_x =
        (x & (DOWNSAMPLE_FACTOR - 1)) / (double)DOWNSAMPLE_FACTOR;
    int flow_y = y >> DOWNSAMPLE_SHIFT;
    double flow_sub_y =
        (y & (DOWNSAMPLE_FACTOR - 1)) / (double)DOWNSAMPLE_FACTOR;

    // Make sure that bicubic interpolation won't read outside of the flow field
    if (flow_x < 1 || (flow_x + 2) >= width) continue;
    if (flow_y < 1 || (flow_y + 2) >= height) continue;

    double h_kernel[4];
    double v_kernel[4];
    getCubicKernel_dbl(flow_sub_x, h_kernel);
    getCubicKernel_dbl(flow_sub_y, v_kernel);

    double flow_u = bicubic_interp_one(&flow->u[flow_y * stride + flow_x],
                                       stride, h_kernel, v_kernel);
    double flow_v = bicubic_interp_one(&flow->v[flow_y * stride + flow_x],
                                       stride, h_kernel, v_kernel);

    // Use original points (without offsets) when filling in correspondence
    // array
    correspondences[num_correspondences].x = x0;
    correspondences[num_correspondences].y = y0;
    correspondences[num_correspondences].rx = x0 + flow_u;
    correspondences[num_correspondences].ry = y0 + flow_v;
    num_correspondences++;
  }
  return num_correspondences;
}

// Compare two regions of width x height pixels, one rooted at position
// (x, y) in frm and the other at (x + u, y + v) in ref.
// This function returns the sum of squared pixel differences between
// the two regions.
static INLINE int compute_flow_error(const uint8_t *ref, const uint8_t *frm,
                                     int width, int height, int stride, int x,
                                     int y, double u, double v, int16_t *dt) {
  // Split offset into integer and fractional parts, and compute cubic
  // interpolation kernels
  int u_int = (int)floor(u);
  int v_int = (int)floor(v);
  double u_frac = u - floor(u);
  double v_frac = v - floor(v);

  int h_kernel[4];
  int v_kernel[4];
  getCubicKernel_int(u_frac, h_kernel);
  getCubicKernel_int(v_frac, v_kernel);

  // Storage for intermediate values between the two convolution directions
  int tmp_[DISFLOW_PATCH_SIZE * (DISFLOW_PATCH_SIZE + 3)];
  int *tmp = tmp_ + DISFLOW_PATCH_SIZE;  // Offset by one row

  // Clamp coordinates so that all pixels we fetch will remain within the
  // allocated border region, but allow them to go far enough out that
  // the border pixels' values do not change.
  // Since we are calculating an 8x8 block, the bottom-right pixel
  // in the block has coordinates (x0 + 7, y0 + 7). Then, the cubic
  // interpolation has 4 taps, meaning that the output of pixel
  // (x_w, y_w) depends on the pixels in the range
  // ([x_w - 1, x_w + 2], [y_w - 1, y_w + 2]).
  //
  // Thus the most extreme coordinates which will be fetched are
  // (x0 - 1, y0 - 1) and (x0 + 9, y0 + 9).
  int x0 = clamp(x + u_int, -9, width);
  int y0 = clamp(y + v_int, -9, height);

  // Horizontal convolution
  for (int i = -1; i < DISFLOW_PATCH_SIZE + 2; ++i) {
    int y_w = y0 + i;
    for (int j = 0; j < DISFLOW_PATCH_SIZE; ++j) {
      int x_w = x0 + j;
      int arr[4];

      arr[0] = (int)ref[y_w * stride + (x_w - 1)];
      arr[1] = (int)ref[y_w * stride + (x_w + 0)];
      arr[2] = (int)ref[y_w * stride + (x_w + 1)];
      arr[3] = (int)ref[y_w * stride + (x_w + 2)];

      // Apply kernel and round, keeping 6 extra bits of precision.
      //
      // 6 is the maximum allowable number of extra bits which will avoid
      // the intermediate values overflowing an int16_t. The most extreme
      // intermediate value occurs when:
      // * The input pixels are [0, 255, 255, 0]
      // * u_frac = 0.5
      // In this case, the un-scaled output is 255 * 1.125 = 286.875.
      // As an integer with 6 fractional bits, that is 18360, which fits
      // in an int16_t. But with 7 fractional bits it would be 36720,
      // which is too large.
      tmp[i * DISFLOW_PATCH_SIZE + j] = ROUND_POWER_OF_TWO(
          getCubicValue_int(arr, h_kernel), CUBIC_PREC_BITS - 6);
    }
  }

  // Vertical convolution
  int sse = 0;
  for (int i = 0; i < DISFLOW_PATCH_SIZE; ++i) {
    for (int j = 0; j < DISFLOW_PATCH_SIZE; ++j) {
      int *p = &tmp[i * DISFLOW_PATCH_SIZE + j];
      int arr[4] = { p[-DISFLOW_PATCH_SIZE], p[0], p[DISFLOW_PATCH_SIZE],
                     p[2 * DISFLOW_PATCH_SIZE] };
      int result = getCubicValue_int(arr, v_kernel);

      // Apply kernel and round.
      // This time, we have to round off the 6 extra bits which were kept
      // earlier, but we also want to keep DISFLOW_DERIV_SCALE_LOG2 extra bits
      // of precision to match the scale of the dx and dy arrays.
      const int round_bits = CUBIC_PREC_BITS + 6 - DISFLOW_DERIV_SCALE_LOG2;
      int warped = ROUND_POWER_OF_TWO(result, round_bits);
      int src_px = frm[(x + j) + (y + i) * stride] << 3;
      int err = warped - src_px;
      sse += err * err;
      dt[i * DISFLOW_PATCH_SIZE + j] = err;
    }
  }

  return sse;
}

static INLINE void sobel_filter(const uint8_t *src, int src_stride,
                                int16_t *dst, int dst_stride, int dir) {
  int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];

  // Sobel filter kernel
  // This must have an overall scale factor equal to DISFLOW_DERIV_SCALE,
  // in order to produce correctly scaled outputs.
  // To work out the scale factor, we multiply two factors:
  //
  // * For the derivative filter (sobel_a), comparing our filter
  //    image[x - 1] - image[x + 1]
  //   to the standard form
  //    d/dx image[x] = image[x+1] - image[x]
  //   tells us that we're actually calculating -2 * d/dx image[2]
  //
  // * For the smoothing filter (sobel_b), all coefficients are positive
  //   so the scale factor is just the sum of the coefficients
  //
  // Thus we need to make sure that DISFLOW_DERIV_SCALE = 2 * sum(sobel_b)
  // (and take care of the - sign from sobel_a elsewhere)
  static const int16_t sobel_a[3] = { 1, 0, -1 };
  static const int16_t sobel_b[3] = { 1, 2, 1 };
  const int taps = 3;

  int im_h = DISFLOW_PATCH_SIZE + taps - 1;
  int im_stride = DISFLOW_PATCH_SIZE;
  const int fo_vert = 1;
  const int fo_horiz = 1;

  // horizontal filter
  const uint8_t *src_horiz = src - fo_vert * src_stride;
  const int16_t *h_kernel = dir ? sobel_a : sobel_b;

  for (int y = 0; y < im_h; ++y) {
    for (int x = 0; x < DISFLOW_PATCH_SIZE; ++x) {
      int sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += h_kernel[k] * src_horiz[y * src_stride + x - fo_horiz + k];
      }
      im_block[y * im_stride + x] = sum;
    }
  }

  // vertical filter
  int16_t *src_vert = im_block + fo_vert * im_stride;
  const int16_t *v_kernel = dir ? sobel_b : sobel_a;

  for (int y = 0; y < DISFLOW_PATCH_SIZE; ++y) {
    for (int x = 0; x < DISFLOW_PATCH_SIZE; ++x) {
      int sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += v_kernel[k] * src_vert[(y - fo_vert + k) * im_stride + x];
      }
      dst[y * dst_stride + x] = sum;
    }
  }
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
// Where the sums are computed over a square window of DISFLOW_PATCH_SIZE.
static INLINE void compute_hessian(const int16_t *dx, int dx_stride,
                                   const int16_t *dy, int dy_stride,
                                   double *M) {
  int tmp[4] = { 0 };

  for (int i = 0; i < DISFLOW_PATCH_SIZE; i++) {
    for (int j = 0; j < DISFLOW_PATCH_SIZE; j++) {
      tmp[0] += dx[i * dx_stride + j] * dx[i * dx_stride + j];
      tmp[1] += dx[i * dx_stride + j] * dy[i * dy_stride + j];
      // Don't compute tmp[2], as it should be equal to tmp[1]
      tmp[3] += dy[i * dy_stride + j] * dy[i * dy_stride + j];
    }
  }

  tmp[2] = tmp[1];

  M[0] = (double)tmp[0];
  M[1] = (double)tmp[1];
  M[2] = (double)tmp[2];
  M[3] = (double)tmp[3];
}

static INLINE void compute_flow_vector(const int16_t *dx, int dx_stride,
                                       const int16_t *dy, int dy_stride,
                                       const int16_t *dt, int dt_stride,
                                       int *b) {
  memset(b, 0, 2 * sizeof(*b));

  for (int i = 0; i < DISFLOW_PATCH_SIZE; i++) {
    for (int j = 0; j < DISFLOW_PATCH_SIZE; j++) {
      b[0] += dx[i * dx_stride + j] * dt[i * dt_stride + j];
      b[1] += dy[i * dy_stride + j] * dt[i * dt_stride + j];
    }
  }
}

static INLINE void invert_2x2(const double *M, double *M_inv) {
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

  // TODO(rachelbarker): Is using regularized values
  // or original values better here?
  M_inv[0] = M_3 * det_inv;
  M_inv[1] = -M[1] * det_inv;
  M_inv[2] = -M[2] * det_inv;
  M_inv[3] = M_0 * det_inv;
}

static INLINE void compute_flow_at_point(const uint8_t *frm, const uint8_t *ref,
                                         int x, int y, int width, int height,
                                         int stride, double *u, double *v) {
  double M[4];
  double M_inv[4];
  int b[2];
  int16_t dt[DISFLOW_PATCH_SIZE * DISFLOW_PATCH_SIZE];
  int16_t dx[DISFLOW_PATCH_SIZE * DISFLOW_PATCH_SIZE];
  int16_t dy[DISFLOW_PATCH_SIZE * DISFLOW_PATCH_SIZE];
  const double o_u = *u;
  const double o_v = *v;

  // Compute gradients within this patch
  const uint8_t *frm_patch = &frm[y * stride + x];
  sobel_filter(frm_patch, stride, dx, DISFLOW_PATCH_SIZE, 1);
  sobel_filter(frm_patch, stride, dy, DISFLOW_PATCH_SIZE, 0);

  compute_hessian(dx, DISFLOW_PATCH_SIZE, dy, DISFLOW_PATCH_SIZE, M);
  invert_2x2(M, M_inv);

  for (int itr = 0; itr < DISFLOW_MAX_ITR; itr++) {
    int sse =
        compute_flow_error(ref, frm, width, height, stride, x, y, *u, *v, dt);
    if (sse <= DISFLOW_SSE_TR) break;

    compute_flow_vector(dx, DISFLOW_PATCH_SIZE, dy, DISFLOW_PATCH_SIZE, dt,
                        DISFLOW_PATCH_SIZE, b);

    // Solve flow equations to find a better estimate for the flow vector
    // at this point
    double delta_u = M_inv[0] * b[0] + M_inv[1] * b[1];
    double delta_v = M_inv[2] * b[0] + M_inv[3] * b[1];

    *u += delta_u * DISFLOW_STEP_SIZE;
    *v += delta_v * DISFLOW_STEP_SIZE;
  }
  if (fabs(*u - o_u) > DISFLOW_PATCH_SIZE ||
      fabs(*v - o_v) > DISFLOW_PATCH_SIZE) {
    *u = o_u;
    *v = o_v;
  }
}

static void fill_flow_field_borders(double *flow, int width, int height,
                                    int stride) {
  // Calculate the bounds of the rectangle which was filled in by
  // compute_flow_field() before calling this function.
  // These indices are inclusive on both ends.
  const int left_index = FLOW_BORDER;
  const int right_index = (width - FLOW_BORDER - 1);
  const int top_index = FLOW_BORDER;
  const int bottom_index = (height - FLOW_BORDER - 1);

  // Left area
  for (int i = top_index; i <= bottom_index; i += 1) {
    double *row = flow + i * stride;
    double left = row[left_index];
    for (int j = 0; j < left_index; j++) {
      row[j] = left;
    }
  }

  // Right area
  for (int i = top_index; i <= bottom_index; i += 1) {
    double *row = flow + i * stride;
    double right = row[right_index];
    for (int j = right_index + 1; j < width; j++) {
      row[j] = right;
    }
  }

  // Top area
  double *top_row = flow + top_index * stride;
  for (int i = 0; i < top_index; i++) {
    double *row = flow + i * stride;
    memcpy(row, top_row, width * sizeof(double));
  }

  // Bottom area
  double *bottom_row = flow + bottom_index * stride;
  for (int i = bottom_index + 1; i < height; i++) {
    double *row = flow + i * stride;
    memcpy(row, bottom_row, width * sizeof(double));
  }
}

// make sure flow_u and flow_v start at 0
static void compute_flow_field(const ImagePyramid *frm_pyr,
                               const ImagePyramid *ref_pyr, FlowField *flow) {
  assert(frm_pyr->n_levels == ref_pyr->n_levels);

  double *flow_u = flow->u;
  double *flow_v = flow->v;

  size_t flow_size = flow->stride * (size_t)flow->height;
  double *u_upscale = aom_malloc(flow_size * sizeof(*u_upscale));
  double *v_upscale = aom_malloc(flow_size * sizeof(*v_upscale));

  // Compute flow field from coarsest to finest level of the pyramid
  for (int level = frm_pyr->n_levels - 1; level >= 0; --level) {
    const PyramidLayer *cur_layer = &frm_pyr->layers[level];
    int cur_width = cur_layer->width;
    int cur_height = cur_layer->height;
    int cur_stride = cur_layer->stride;

    const uint8_t *src_buffer = cur_layer->buffer;
    const uint8_t *ref_buffer = ref_pyr->layers[level].buffer;

    int cur_flow_width = cur_width >> DOWNSAMPLE_SHIFT;
    int cur_flow_height = cur_height >> DOWNSAMPLE_SHIFT;
    int cur_flow_stride = flow->stride;

    for (int i = FLOW_BORDER; i < cur_flow_height - FLOW_BORDER; i += 1) {
      for (int j = FLOW_BORDER; j < cur_flow_width - FLOW_BORDER; j += 1) {
        int flow_field_idx = i * cur_flow_stride + j;  // In flow field entries

        // Calculate the position of a patch of size DISFLOW_PATCH_SIZE pixels,
        // which is centered on the region covered by this flow field entry
        int patch_center_x =
            (j << DOWNSAMPLE_SHIFT) + UPSAMPLE_CENTER_OFFSET;  // In pixels
        int patch_center_y =
            (i << DOWNSAMPLE_SHIFT) + UPSAMPLE_CENTER_OFFSET;  // In pixels
        int patch_tl_x = patch_center_x - DISFLOW_PATCH_CENTER;
        int patch_tl_y = patch_center_y - DISFLOW_PATCH_CENTER;
        assert(patch_tl_x >= 0);
        assert(patch_tl_y >= 0);

        compute_flow_at_point(src_buffer, ref_buffer, patch_tl_x, patch_tl_y,
                              cur_width, cur_height, cur_stride,
                              &flow_u[flow_field_idx], &flow_v[flow_field_idx]);
      }
    }

    // Fill in the areas which we haven't explicitly computed, with copies
    // of the outermost values which we did compute
    fill_flow_field_borders(flow_u, cur_flow_width, cur_flow_height,
                            cur_flow_stride);
    fill_flow_field_borders(flow_v, cur_flow_width, cur_flow_height,
                            cur_flow_stride);

    if (level > 0) {
      int upscale_flow_width = cur_flow_width << 1;
      int upscale_flow_height = cur_flow_height << 1;
      int upscale_stride = flow->stride;

      av1_upscale_plane_double_prec(
          flow_u, cur_flow_height, cur_flow_width, cur_flow_stride, u_upscale,
          upscale_flow_height, upscale_flow_width, upscale_stride);
      av1_upscale_plane_double_prec(
          flow_v, cur_flow_height, cur_flow_width, cur_flow_stride, v_upscale,
          upscale_flow_height, upscale_flow_width, upscale_stride);

      // Multiply all flow vectors by 2.
      // When we move down a pyramid level, the image resolution doubles.
      // Thus we need to double all vectors in order for them to represent
      // the same translation at the next level down
      for (int i = 0; i < upscale_flow_height; i++) {
        for (int j = 0; j < upscale_flow_width; j++) {
          int index = i * upscale_stride + j;
          flow_u[index] = u_upscale[index] * 2.0;
          flow_v[index] = v_upscale[index] * 2.0;
        }
      }

      // If we didn't fill in the rightmost column or bottommost row during
      // upsampling (in order to keep the ratio to exactly 2), fill them
      // in here by copying the next closest column/row
      const PyramidLayer *next_layer = &frm_pyr->layers[level - 1];
      int next_flow_width = next_layer->width >> DOWNSAMPLE_SHIFT;
      int next_flow_height = next_layer->height >> DOWNSAMPLE_SHIFT;

      // Rightmost column
      if (next_flow_width > upscale_flow_width) {
        assert(next_flow_width == upscale_flow_width + 1);
        for (int i = 0; i < upscale_flow_height; i++) {
          int index = i * upscale_stride + upscale_flow_width;
          flow_u[index] = flow_u[index - 1];
          flow_v[index] = flow_v[index - 1];
        }
      }

      // Bottommost row
      if (next_flow_height > upscale_flow_height) {
        assert(next_flow_height == upscale_flow_height + 1);
        for (int j = 0; j < next_flow_width; j++) {
          int index = upscale_flow_height * upscale_stride + j;
          flow_u[index] = flow_u[index - upscale_stride];
          flow_v[index] = flow_v[index - upscale_stride];
        }
      }
    }
  }
  aom_free(u_upscale);
  aom_free(v_upscale);
}

FlowField *aom_alloc_flow_field(int frame_width, int frame_height) {
  FlowField *flow = (FlowField *)aom_malloc(sizeof(FlowField));
  if (flow == NULL) return NULL;

  // Calculate the size of the bottom (largest) layer of the flow pyramid
  flow->width = frame_width >> DOWNSAMPLE_SHIFT;
  flow->height = frame_height >> DOWNSAMPLE_SHIFT;
  flow->stride = flow->width;

  size_t flow_size = flow->stride * (size_t)flow->height;
  flow->u = aom_calloc(flow_size, sizeof(double));
  flow->v = aom_calloc(flow_size, sizeof(double));

  if (flow->u == NULL || flow->v == NULL) {
    aom_free(flow->u);
    aom_free(flow->v);
    aom_free(flow);
    return NULL;
  }

  return flow;
}

void aom_free_flow_field(FlowField *flow) {
  aom_free(flow->u);
  aom_free(flow->v);
  aom_free(flow);
}

FlowField *aom_compute_flow_field(YV12_BUFFER_CONFIG *frm,
                                  YV12_BUFFER_CONFIG *ref, int bit_depth) {
  const int frm_width = frm->y_width;
  const int frm_height = frm->y_height;
  assert(frm->y_width == ref->y_width);
  assert(frm->y_height == ref->y_height);

  const ImagePyramid *frm_pyr =
      aom_compute_pyramid(frm, bit_depth, MAX_PYRAMID_LEVELS);
  const ImagePyramid *ref_pyr =
      aom_compute_pyramid(ref, bit_depth, MAX_PYRAMID_LEVELS);

  FlowField *flow = aom_alloc_flow_field(frm_width, frm_height);

  compute_flow_field(frm_pyr, ref_pyr, flow);

  return flow;
}

bool aom_fit_global_model_to_flow_field(const FlowField *flow,
                                        TransformationType type,
                                        YV12_BUFFER_CONFIG *frm, int bit_depth,
                                        MotionModel *params_by_motion,
                                        int num_motions) {
  int num_correspondences;

  aom_find_corners_in_frame(frm, bit_depth);

  // find correspondences between the two images using the flow field
  Correspondence *correspondences =
      aom_malloc(frm->num_corners * sizeof(*correspondences));
  num_correspondences = determine_disflow_correspondence(
      frm->corners, frm->num_corners, flow, correspondences);

  ransac(correspondences, num_correspondences, type, params_by_motion,
         num_motions);

  aom_free(correspondences);

  // Set num_inliers = 0 for motions with too few inliers so they are ignored.
  for (int i = 0; i < num_motions; ++i) {
    if (params_by_motion[i].num_inliers <
        MIN_INLIER_PROB * num_correspondences) {
      params_by_motion[i].num_inliers = 0;
    }
  }

  // Return true if any one of the motions has inliers.
  for (int i = 0; i < num_motions; ++i) {
    if (params_by_motion[i].num_inliers > 0) return true;
  }
  return false;
}

bool aom_fit_local_model_to_flow_field(const FlowField *flow,
                                       const PixelRect *rect,
                                       TransformationType type, double *mat) {
  // Transform input rectangle to flow-field space
  // Generally `rect` will be the rectangle of a single coding block,
  // so the edges will be aligned to multiples of DOWNSAMPLE_FACTOR already.
  PixelRect downsampled_rect = { .left = rect->left >> DOWNSAMPLE_SHIFT,
                                 .right = rect->right >> DOWNSAMPLE_SHIFT,
                                 .top = rect->top >> DOWNSAMPLE_SHIFT,
                                 .bottom = rect->bottom >> DOWNSAMPLE_SHIFT };

  // Generate one point for each flow field entry covered by the rectangle
  int width = rect_height(&downsampled_rect);
  int height = rect_width(&downsampled_rect);

  int num_points = width * height;

  double *pts1 = aom_malloc(num_points * 2 * sizeof(double));
  double *pts2 = aom_malloc(num_points * 2 * sizeof(double));

  int flow_stride = flow->stride;

  int index = 0;
  for (int i = rect->top; i < rect->bottom; i++) {
    for (int j = rect->left; j < rect->right; j++) {
      int flow_pos = i * flow_stride + j;
      // Associate each flow field entry with the center-most pixel that
      // it covers
      int patch_center_x = (j << DOWNSAMPLE_SHIFT) + UPSAMPLE_CENTER_OFFSET;
      int patch_center_y = (i << DOWNSAMPLE_SHIFT) + UPSAMPLE_CENTER_OFFSET;
      pts1[2 * index + 0] = (double)patch_center_x;
      pts1[2 * index + 1] = (double)patch_center_y;
      pts2[2 * index + 0] = (double)patch_center_x + flow->u[flow_pos];
      pts2[2 * index + 1] = (double)patch_center_y + flow->v[flow_pos];
      index++;
    }
  }

  // Check that we filled the expected number of points
  assert(index == num_points);

  bool result = aom_fit_motion_model(type, num_points, pts1, pts2, mat);

  aom_free(pts1);
  aom_free(pts2);
  return result;
}
