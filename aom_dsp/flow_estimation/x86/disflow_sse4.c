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

#include <math.h>
#include <smmintrin.h>

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/flow_estimation/disflow.h"

#include "config/aom_dsp_rtcd.h"

// Internal cross-check against C code
// If you set this to 1 and compile in debug mode, then the outputs of the two
// convolution stages will be checked against the plain C version of the code,
// and an assertion will be fired if the results differ.
#define CHECK_RESULTS 0

// Internal precision of cubic interpolation filters
// The limiting factor here is that:
// * Before integerizing, the maximum value of any kernel tap is 1.0
// * After integerizing, each tap must fit into an int16_t.
// Thus the largest multiplier we can get away with is 2^14 = 16384,
// as 2^15 = 32768 is too large to fit in an int16_t.
#define CUBIC_PREC_BITS 14

// Note: Max sum(+ve coefficients) = 1.125 * scale
static INLINE void getCubicKernel_int(double x, int16_t *kernel) {
  assert(0 <= x && x < 1);
  double x2 = x * x;
  double x3 = x2 * x;

  double k0 = -0.5 * x + x2 - 0.5 * x3;
  double k1 = 1.0 - 2.5 * x2 + 1.5 * x3;
  double k2 = 0.5 * x + 2.0 * x2 - 1.5 * x3;
  double k3 = -0.5 * x2 + 0.5 * x3;

  kernel[0] = (int16_t)rint(k0 * (1 << CUBIC_PREC_BITS));
  kernel[1] = (int16_t)rint(k1 * (1 << CUBIC_PREC_BITS));
  kernel[2] = (int16_t)rint(k2 * (1 << CUBIC_PREC_BITS));
  kernel[3] = (int16_t)rint(k3 * (1 << CUBIC_PREC_BITS));
}

#if CHECK_RESULTS
static INLINE int getCubicValue_int(int *p, int16_t *kernel) {
  return kernel[0] * p[0] + kernel[1] * p[1] + kernel[2] * p[2] +
         kernel[3] * p[3];
}
#endif  // CHECK_RESULTS

// Compare two regions of width x height pixels, one rooted at position
// (x, y) in frm and the other at (x + u, y + v) in ref.
// This function returns the sum of squared pixel differences between
// the two regions.
//
// TODO(rachelbarker): Explore the speed vs. accuracy tradeoff space
// for this function. Eg:
//
// * Keep multiplication by 256, but round after first stage
//   This would allow us to use 16 x 16 -> 32 bit multiplies in both stages
//
// * Use linear interpolation instead of bicubic interpolation
//   This would halve the number of multiplies needed, as well as simplifying
//   the non-vectorized initial part of the function
static INLINE int compute_flow_error(unsigned char *ref, unsigned char *frm,
                                     int width, int height, int stride, int x,
                                     int y, double u, double v, int16_t *dt) {
  // This function is written to do 8x8 convolutions only
  assert(DISFLOW_PATCH_SIZE == 8);

  // Split offset into integer and fractional parts, and compute cubic
  // interpolation kernels
  int u_int = (int)floor(u);
  int v_int = (int)floor(v);
  double u_frac = u - floor(u);
  double v_frac = v - floor(v);

  int16_t h_kernel[4];
  int16_t v_kernel[4];
  getCubicKernel_int(u_frac, h_kernel);
  getCubicKernel_int(v_frac, v_kernel);

  // Storage for intermediate values between the two convolution directions
  int16_t tmp_[DISFLOW_PATCH_SIZE * (DISFLOW_PATCH_SIZE + 3)];
  int16_t *tmp = tmp_ + DISFLOW_PATCH_SIZE;  // Offset by one row

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

  // Prepare the kernel vectors
  // We split the kernel into two vectors with kernel indices:
  // 0, 1, 0, 1, 0, 1, 0, 1, and
  // 2, 3, 2, 3, 2, 3, 2, 3
  __m128i h_kernel_01 =
      _mm_set1_epi32((h_kernel[0] & 0xFFFF) | (h_kernel[1] << 16));
  __m128i h_kernel_23 =
      _mm_set1_epi32((h_kernel[2] & 0xFFFF) | (h_kernel[3] << 16));

  __m128i round_const_h = _mm_set1_epi32(1 << (CUBIC_PREC_BITS - 6 - 1));

  for (int i = -1; i < DISFLOW_PATCH_SIZE + 2; ++i) {
    int y_w = y0 + i;
    const unsigned char *ref_row = &ref[y_w * stride + (x0 - 1)];
    int16_t *tmp_row = &tmp[i * DISFLOW_PATCH_SIZE];

    // Load this row of pixels.
    // For an 8x8 patch, we need to load the 8 image pixels + 3 extras,
    // for a total of 11 pixels. Here we load 16 pixels, but only use
    // the first 11.
    __m128i row = _mm_loadu_si128((__m128i *)ref_row);

    // Expand pixels to int16s
    __m128i px_0to7_i16 = _mm_cvtepu8_epi16(row);
    __m128i px_4to10_i16 = _mm_cvtepu8_epi16(_mm_srli_si128(row, 4));

    // Relevant multiply instruction
    // This multiplies pointwise, then sums in pairs.
    //_mm_madd_epi16();

    // Compute first four outputs
    // input pixels 0, 1, 1, 2, 2, 3, 3, 4
    // * kernel     0, 1, 0, 1, 0, 1, 0, 1
    __m128i px0 =
        _mm_unpacklo_epi16(px_0to7_i16, _mm_srli_si128(px_0to7_i16, 2));
    // input pixels 2, 3, 3, 4, 4, 5, 5, 6
    // * kernel     2, 3, 2, 3, 2, 3, 2, 3
    __m128i px1 = _mm_unpacklo_epi16(_mm_srli_si128(px_0to7_i16, 4),
                                     _mm_srli_si128(px_0to7_i16, 6));
    // Convolve with kernel and sum 2x2 boxes to form first 4 outputs
    __m128i sum0 = _mm_add_epi32(_mm_madd_epi16(px0, h_kernel_01),
                                 _mm_madd_epi16(px1, h_kernel_23));

    __m128i out0 =
        _mm_srai_epi32(_mm_add_epi32(sum0, round_const_h), CUBIC_PREC_BITS - 6);

    // Compute second four outputs
    __m128i px2 =
        _mm_unpacklo_epi16(px_4to10_i16, _mm_srli_si128(px_4to10_i16, 2));
    __m128i px3 = _mm_unpacklo_epi16(_mm_srli_si128(px_4to10_i16, 4),
                                     _mm_srli_si128(px_4to10_i16, 6));
    __m128i sum1 = _mm_add_epi32(_mm_madd_epi16(px2, h_kernel_01),
                                 _mm_madd_epi16(px3, h_kernel_23));

    // Round by just enough bits that the result is
    // guaranteed to fit into an i16. Then the next stage can use 16 x 16 -> 32
    // bit multiplies, which should be a fair bit faster than 32 x 32 -> 32
    // as it does now
    // This means shifting down so we have 6 extra bits, for a maximum value
    // of +18360, which can occur if u_frac == 0.5 and the input pixels are
    // {0, 255, 255, 0}.
    __m128i out1 =
        _mm_srai_epi32(_mm_add_epi32(sum1, round_const_h), CUBIC_PREC_BITS - 6);

    _mm_storeu_si128((__m128i *)tmp_row, _mm_packs_epi32(out0, out1));

#if CHECK_RESULTS
    // Cross-check
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
      int c_value = ROUND_POWER_OF_TWO(getCubicValue_int(arr, h_kernel), 4);
      (void)c_value;  // Suppress warnings
      assert(tmp_row[j] == c_value);
    }
#endif  // CHECK_RESULTS
  }

  // Vertical convolution
  // Accumulator for summed squared errors.
  // This is split into four partial sums which are 32 bits wide each,
  // which must be summed at the end to find the true SSE.
  __m128i partial_sse = _mm_setzero_si128();
#if CHECK_RESULTS
  int c_sse = 0;
#endif  // CHECK_RESULTS

  __m128i round_const_v = _mm_set1_epi32(1 << (CUBIC_PREC_BITS + 6 - 1));
  __m128i clamp_const = _mm_set1_epi16(255);

  __m128i v_kernel_01 =
      _mm_set1_epi32((v_kernel[0] & 0xFFFF) | (v_kernel[1] << 16));
  __m128i v_kernel_23 =
      _mm_set1_epi32((v_kernel[2] & 0xFFFF) | (v_kernel[3] << 16));

  for (int i = 0; i < DISFLOW_PATCH_SIZE; ++i) {
    int16_t *tmp_row = &tmp[i * DISFLOW_PATCH_SIZE];

    // Load 4 rows of 8 x 16-bit values
    __m128i px0 = _mm_loadu_si128((__m128i *)(tmp_row - DISFLOW_PATCH_SIZE));
    __m128i px1 = _mm_loadu_si128((__m128i *)tmp_row);
    __m128i px2 = _mm_loadu_si128((__m128i *)(tmp_row + DISFLOW_PATCH_SIZE));
    __m128i px3 =
        _mm_loadu_si128((__m128i *)(tmp_row + 2 * DISFLOW_PATCH_SIZE));

    // We want to calculate px0 * v_kernel[0] + px1 * v_kernel[1] + ... ,
    // but each multiply expands its output to 32 bits. So we need to be
    // a little clever about how we do this
    __m128i sum0 = _mm_add_epi32(
        _mm_madd_epi16(_mm_unpacklo_epi16(px0, px1), v_kernel_01),
        _mm_madd_epi16(_mm_unpacklo_epi16(px2, px3), v_kernel_23));
    __m128i sum1 = _mm_add_epi32(
        _mm_madd_epi16(_mm_unpackhi_epi16(px0, px1), v_kernel_01),
        _mm_madd_epi16(_mm_unpackhi_epi16(px2, px3), v_kernel_23));

    __m128i sum0_rounded =
        _mm_srai_epi32(_mm_add_epi32(sum0, round_const_v), CUBIC_PREC_BITS + 6);
    __m128i sum1_rounded =
        _mm_srai_epi32(_mm_add_epi32(sum1, round_const_v), CUBIC_PREC_BITS + 6);

    // Pack the outputs into one vector of 8 16-bit values,
    // and clamp to the range [0, 255].
    // The easiest way to do this, given the intrinsics available,
    // is to reduce from 32 to 16 bits with clamping to [0, 65535],
    // then using min() to clamp to [0, 255].
    __m128i out_i16 = _mm_packus_epi32(sum0_rounded, sum1_rounded);
    __m128i warped = _mm_min_epi16(out_i16, clamp_const);

    // Calculate delta from the target patch
    __m128i frm_pixels = _mm_loadu_si64((__m128i *)&frm[(y + i) * stride + x]);
    __m128i err = _mm_sub_epi16(warped, _mm_cvtepu8_epi16(frm_pixels));
    _mm_storeu_si128((__m128i *)&dt[i * DISFLOW_PATCH_SIZE], err);

    // Calculate squared errors, sum in pairs, and add to accumulator.
    // Each value in the accumulator must be 32 bits wide, as the individual
    // error values are in the range [-255, 255] and thus a single squared
    // error value just about fits in 16 bits. So a sum of many couple overflow
    // 16 bits, but not 32 bits.
    partial_sse = _mm_add_epi32(partial_sse, _mm_madd_epi16(err, err));

#if CHECK_RESULTS
    for (int j = 0; j < DISFLOW_PATCH_SIZE; ++j) {
      int16_t *p = &tmp[i * DISFLOW_PATCH_SIZE + j];
      int arr[4] = { p[-DISFLOW_PATCH_SIZE], p[0], p[DISFLOW_PATCH_SIZE],
                     p[2 * DISFLOW_PATCH_SIZE] };
      int result = getCubicValue_int(arr, v_kernel);

      // Apply kernel and round.
      // This time, we want the result to be at the same precision as the input
      // pixels, so we have to round off the kernel precision plus the 6 extra
      // bits used in the intermediate array.
      int c_warped = clamp(ROUND_POWER_OF_TWO(result, 16), 0, 255);
      int c_err = c_warped - frm[(x + j) + (y + i) * stride];
      c_sse += c_err * c_err;
      assert(dt[i * DISFLOW_PATCH_SIZE + j] == c_err);
    }
#endif  // CHECK_RESULTS
  }

  // Finish accumulating errors
  int sse =
      _mm_extract_epi32(partial_sse, 0) + _mm_extract_epi32(partial_sse, 1) +
      _mm_extract_epi32(partial_sse, 2) + _mm_extract_epi32(partial_sse, 3);

#if CHECK_RESULTS
  (void)sse;  // Suppress warnings
  assert(sse == c_sse);
#endif  // CHECK_RESULTS

  return sse;
}

static INLINE void compute_flow_vector(const int16_t *dx, int dx_stride,
                                       const int16_t *dy, int dy_stride,
                                       const int16_t *dt, int dt_stride,
                                       int *b) {
  __m128i b0_acc = _mm_setzero_si128();
  __m128i b1_acc = _mm_setzero_si128();

  for (int i = 0; i < DISFLOW_PATCH_SIZE; i++) {
    // Need to load 8 values of dx, 8 of dy, 8 of dt, which conveniently
    // works out to one register each. Then just calculate dx * dt, dy * dt,
    // and (implicitly) sum horizontally in pairs.
    // This gives four 32-bit partial sums for each of b[0] and b[1],
    // which can be accumulated and summed at the end.
    __m128i dx_row = _mm_loadu_si128((__m128i *)&dx[i * dx_stride]);
    __m128i dy_row = _mm_loadu_si128((__m128i *)&dy[i * dy_stride]);
    __m128i dt_row = _mm_loadu_si128((__m128i *)&dt[i * dt_stride]);

    b0_acc = _mm_add_epi32(b0_acc, _mm_madd_epi16(dx_row, dt_row));
    b1_acc = _mm_add_epi32(b1_acc, _mm_madd_epi16(dy_row, dt_row));
  }

  // We need to set b[0] = sum(b0_acc), b[1] = sum(b1_acc).
  // We might as well use a `hadd` instruction to do 4 of the additions
  // needed here. Then that just leaves two more additions, which can be
  // done in scalar code
  __m128i partial_sum = _mm_hadd_epi32(b0_acc, b1_acc);
  b[0] = _mm_extract_epi32(partial_sum, 0) + _mm_extract_epi32(partial_sum, 1);
  b[1] = _mm_extract_epi32(partial_sum, 2) + _mm_extract_epi32(partial_sum, 3);

#if CHECK_RESULTS
  int c_result[2];
  aom_compute_flow_vector_c(dx, dx_stride, dy, dy_stride, dt, dt_stride,
                            c_result);
  assert(b[0] == c_result[0]);
  assert(b[1] == c_result[1]);
#endif  // CHECK_RESULTS
}

static INLINE void sobel_filter(const uint8_t *src, int src_stride,
                                int16_t *dst, int dst_stride, int dir) {
  int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];
  DECLARE_ALIGNED(256, static const int16_t, sobel_a[4]) = { 1, 0, -1, 0 };
  DECLARE_ALIGNED(256, static const int16_t, sobel_b[4]) = { 1, 2, 1, 0 };
  const int taps = 4;
  int im_h = DISFLOW_PATCH_SIZE + taps - 1;
  int im_stride = DISFLOW_PATCH_SIZE;
  const int fo_vert = 1;
  const int fo_horiz = 1;

  // horizontal filter
  const uint8_t *src_horiz = src - fo_vert * src_stride - fo_horiz;
  const int16_t *h_kernel = dir ? sobel_a : sobel_b;

  __m128i h_kernel_01 =
      _mm_set1_epi32((h_kernel[0] & 0xFFFF) | (h_kernel[1] << 16));
  __m128i h_kernel_23 =
      _mm_set1_epi32((h_kernel[2] & 0xFFFF) | (h_kernel[3] << 16));

  for (int y = 0; y < im_h; ++y) {
    const uint8_t *src_row = src_horiz + y * src_stride;
    int16_t *im_row = im_block + y * im_stride;

    // Load pixels and expand to 16 bits
    __m128i row = _mm_loadu_si128((__m128i *)src_row);
    __m128i px_0to7_i16 = _mm_cvtepu8_epi16(row);
    __m128i px_4to10_i16 = _mm_cvtepu8_epi16(_mm_srli_si128(row, 4));

    // Compute first four outputs
    // input pixels 0, 1, 1, 2, 2, 3, 3, 4
    // * kernel     0, 1, 0, 1, 0, 1, 0, 1
    __m128i px0 =
        _mm_unpacklo_epi16(px_0to7_i16, _mm_srli_si128(px_0to7_i16, 2));
    // input pixels 2, 3, 3, 4, 4, 5, 5, 6
    // * kernel     2, 3, 2, 3, 2, 3, 2, 3
    __m128i px1 = _mm_unpacklo_epi16(_mm_srli_si128(px_0to7_i16, 4),
                                     _mm_srli_si128(px_0to7_i16, 6));
    // Convolve with kernel and sum 2x2 boxes to form first 4 outputs
    __m128i sum0 = _mm_add_epi32(_mm_madd_epi16(px0, h_kernel_01),
                                 _mm_madd_epi16(px1, h_kernel_23));

    // Compute second four outputs
    __m128i px2 =
        _mm_unpacklo_epi16(px_4to10_i16, _mm_srli_si128(px_4to10_i16, 2));
    __m128i px3 = _mm_unpacklo_epi16(_mm_srli_si128(px_4to10_i16, 4),
                                     _mm_srli_si128(px_4to10_i16, 6));
    __m128i sum1 = _mm_add_epi32(_mm_madd_epi16(px2, h_kernel_01),
                                 _mm_madd_epi16(px3, h_kernel_23));

    // Store to intermediate array
    _mm_storeu_si128((__m128i *)im_row, _mm_packs_epi32(sum0, sum1));

#if CHECK_RESULTS
    // Cross-check
    for (int x = 0; x < DISFLOW_PATCH_SIZE; ++x) {
      int sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += h_kernel[k] * src_row[x + k];
      }
      (void)sum;
      assert(im_row[x] == sum);
    }
#endif  // CHECK_RESULTS
  }

  // vertical filter
  const int16_t *v_kernel = dir ? sobel_b : sobel_a;

  __m128i v_kernel_01 =
      _mm_set1_epi32((v_kernel[0] & 0xFFFF) | (v_kernel[1] << 16));
  __m128i v_kernel_23 =
      _mm_set1_epi32((v_kernel[2] & 0xFFFF) | (v_kernel[3] << 16));

  for (int y = 0; y < DISFLOW_PATCH_SIZE; ++y) {
    // TODO: Rework this to look more like the other function
    const int16_t *im_row = im_block + (y + 1) * im_stride;
    int16_t *dst_row = dst + y * dst_stride;

    // Load 4 rows of 8 x 16-bit values
    __m128i px0 = _mm_loadu_si128((__m128i *)(im_row - im_stride));
    __m128i px1 = _mm_loadu_si128((__m128i *)im_row);
    __m128i px2 = _mm_loadu_si128((__m128i *)(im_row + im_stride));
    __m128i px3 = _mm_loadu_si128((__m128i *)(im_row + 2 * im_stride));

    // We want to calculate px0 * v_kernel[0] + px1 * v_kernel[1] + ... ,
    // but each multiply expands its output to 32 bits. So we need to be
    // a little clever about how we do this
    __m128i sum0 = _mm_add_epi32(
        _mm_madd_epi16(_mm_unpacklo_epi16(px0, px1), v_kernel_01),
        _mm_madd_epi16(_mm_unpacklo_epi16(px2, px3), v_kernel_23));
    __m128i sum1 = _mm_add_epi32(
        _mm_madd_epi16(_mm_unpackhi_epi16(px0, px1), v_kernel_01),
        _mm_madd_epi16(_mm_unpackhi_epi16(px2, px3), v_kernel_23));

    _mm_storeu_si128((__m128i *)dst_row, _mm_packs_epi32(sum0, sum1));

#if CHECK_RESULTS
    for (int x = 0; x < DISFLOW_PATCH_SIZE; ++x) {
      int sum = 0;
      for (int k = 0; k < taps; ++k) {
        sum += v_kernel[k] * im_row[(k - 1) * im_stride + x];
      }
      (void)sum;
      assert(dst_row[x] == sum);
    }
#endif  // CHECK_RESULTS
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
  __m128i acc[4] = { 0 };

  for (int i = 0; i < DISFLOW_PATCH_SIZE; i++) {
    __m128i dx_row = _mm_loadu_si128((__m128i *)&dx[i * dx_stride]);
    __m128i dy_row = _mm_loadu_si128((__m128i *)&dy[i * dy_stride]);

    acc[0] = _mm_add_epi32(acc[0], _mm_madd_epi16(dx_row, dx_row));
    acc[1] = _mm_add_epi32(acc[1], _mm_madd_epi16(dx_row, dy_row));
    // Don't compute acc[2], as it should be equal to acc[1]
    acc[3] = _mm_add_epi32(acc[3], _mm_madd_epi16(dy_row, dy_row));
  }

  // Condense sums
  __m128i partial_sum_0 = _mm_hadd_epi32(acc[0], acc[1]);
  __m128i partial_sum_1 = _mm_hadd_epi32(acc[1], acc[3]);
  __m128i result = _mm_hadd_epi32(partial_sum_0, partial_sum_1);

#if CHECK_RESULTS
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

  assert(tmp[0] == _mm_extract_epi32(result, 0));
  assert(tmp[1] == _mm_extract_epi32(result, 1));
  assert(tmp[2] == _mm_extract_epi32(result, 2));
  assert(tmp[3] == _mm_extract_epi32(result, 3));
#endif  // CHECK_RESULTS

  // Convert results to doubles and store
  _mm_storeu_pd(M, _mm_cvtepi32_pd(result));
  _mm_storeu_pd(M + 2, _mm_cvtepi32_pd(_mm_srli_si128(result, 8)));
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

void aom_compute_flow_at_point_sse4_1(unsigned char *frm, unsigned char *ref,
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
  unsigned char *frm_patch = &frm[y * stride + x];
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
    // Multiply flow vector by the inverse of the Hessian to find
    // what size and direction of step we should take
    *u += M_inv[0] * b[0] + M_inv[1] * b[1];
    *v += M_inv[2] * b[0] + M_inv[3] * b[1];
  }
  if (fabs(*u - o_u) > DISFLOW_PATCH_SIZE ||
      fabs(*v - o_v) > DISFLOW_PATCH_SIZE) {
    *u = o_u;
    *v = o_v;
  }
}
