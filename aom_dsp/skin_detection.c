/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom_dsp/skin_detection.h"

#define MODEL_MODE 1

// Fixed-point skin color model parameters.
static const int skin_mean[5][2] = { { 7463, 9614 },
                                     { 6400, 10240 },
                                     { 7040, 10240 },
                                     { 8320, 9280 },
                                     { 6800, 9614 } };
static const int skin_inv_cov[4] = { 4107, 1663, 1663, 2157 };  // q16
static const int skin_threshold[6] = { 1570636, 1400000, 800000,
                                       800000,  800000,  800000 };  // q18
// Thresholds on luminance.
static const int y_low = 40;
static const int y_high = 220;

// Evaluates the Mahalanobis distance measure for the input CbCr values.
static int aom_evaluate_skin_color_difference(const int cb, const int cr,
                                              const int idx) {
  const int cb_q6 = cb << 6;
  const int cr_q6 = cr << 6;
  const int cb_diff_q12 =
      (cb_q6 - skin_mean[idx][0]) * (cb_q6 - skin_mean[idx][0]);
  const int cbcr_diff_q12 =
      (cb_q6 - skin_mean[idx][0]) * (cr_q6 - skin_mean[idx][1]);
  const int cr_diff_q12 =
      (cr_q6 - skin_mean[idx][1]) * (cr_q6 - skin_mean[idx][1]);
  const int cb_diff_q2 = (cb_diff_q12 + (1 << 9)) >> 10;
  const int cbcr_diff_q2 = (cbcr_diff_q12 + (1 << 9)) >> 10;
  const int cr_diff_q2 = (cr_diff_q12 + (1 << 9)) >> 10;
  const int skin_diff =
      skin_inv_cov[0] * cb_diff_q2 + skin_inv_cov[1] * cbcr_diff_q2 +
      skin_inv_cov[2] * cbcr_diff_q2 + skin_inv_cov[3] * cr_diff_q2;
  return skin_diff;
}

// Checks if the input yCbCr values corresponds to skin color.
int aom_skin_pixel(const int y, const int cb, const int cr, int motion) {
  if (y < y_low || y > y_high) {
    return 0;
  } else if (MODEL_MODE == 0) {
    return (aom_evaluate_skin_color_difference(cb, cr, 0) < skin_threshold[0]);
  } else {
    int i = 0;
    // Exit on grey.
    if (cb == 128 && cr == 128) return 0;
    // Exit on very strong cb.
    if (cb > 150 && cr < 110) return 0;
    for (; i < 5; ++i) {
      int skin_color_diff = aom_evaluate_skin_color_difference(cb, cr, i);
      if (skin_color_diff < skin_threshold[i + 1]) {
        if (y < 60 && skin_color_diff > 3 * (skin_threshold[i + 1] >> 2)) {
          return 0;
        } else if (motion == 0 &&
                   skin_color_diff > (skin_threshold[i + 1] >> 1)) {
          return 0;
        } else {
          return 1;
        }
      }
      // Exit if difference is much large than the threshold.
      if (skin_color_diff > (skin_threshold[i + 1] << 3)) {
        return 0;
      }
    }
    return 0;
  }
}
