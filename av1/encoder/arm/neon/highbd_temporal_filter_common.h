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
#ifndef AOM_AV1_ENCODER_ARM_NEON_HIGHBD_TEMPORAL_FILTER_COMMON_H
#define AOM_AV1_ENCODER_ARM_NEON_HIGHBD_TEMPORAL_FILTER_COMMON_H

#include "config/aom_config.h"
#include "config/av1_rtcd.h"

#include "aom_dsp/mathutils.h"

static AOM_FORCE_INLINE void do_filtering(
    const uint16_t *frame, int stride, uint32_t acc_5x5_neon[32][32],
    const uint32_t *luma_sse_sum, const double inv_num_ref_pixels,
    const uint32_t block_height, const uint32_t block_width,
    const int *subblock_mses, const double weight_factor,
    const double inv_factor, const double decay_factor, const double *d_factor,
    unsigned int *accumulator, uint16_t *count, int tf_wgt_calc_lvl, int bd) {
  if (tf_wgt_calc_lvl == 0) {
    for (unsigned int i = 0, k = 0; i < block_height; i++) {
      for (unsigned int j = 0; j < block_width; j++, k++) {
        const int pixel_value = frame[i * stride + j];
        // Scale down the difference for high bit depth input.
        const uint32_t diff_sse =
            (acc_5x5_neon[i][j] + luma_sse_sum[i * BW + j]) >> ((bd - 8) * 2);

        const double window_error = diff_sse * inv_num_ref_pixels;
        const int subblock_idx =
            (i >= block_height / 2) * 2 + (j >= block_width / 2);
        const double block_error = (double)subblock_mses[subblock_idx];
        const double combined_error =
            weight_factor * window_error + block_error * inv_factor;
        // Compute filter weight.
        double scaled_error =
            combined_error * d_factor[subblock_idx] * decay_factor;
        scaled_error = AOMMIN(scaled_error, 7);
        const int weight = (int)(exp(-scaled_error) * TF_WEIGHT_SCALE);
        accumulator[k] += weight * pixel_value;
        count[k] += weight;
      }
    }
  } else {
    for (unsigned int i = 0, k = 0; i < block_height; i++) {
      for (unsigned int j = 0; j < block_width; j++, k++) {
        const int pixel_value = frame[i * stride + j];
        // Scale down the difference for high bit depth input.
        const uint32_t diff_sse =
            (acc_5x5_neon[i][j] + luma_sse_sum[i * BW + j]) >> ((bd - 8) * 2);

        const double window_error = diff_sse * inv_num_ref_pixels;
        const int subblock_idx =
            (i >= block_height / 2) * 2 + (j >= block_width / 2);
        const double block_error = (double)subblock_mses[subblock_idx];
        const double combined_error =
            weight_factor * window_error + block_error * inv_factor;
        // Compute filter weight.
        double scaled_error =
            combined_error * d_factor[subblock_idx] * decay_factor;
        scaled_error = AOMMIN(scaled_error, 7);
        const float fweight =
            approx_exp((float)-scaled_error) * TF_WEIGHT_SCALE;
        const int weight = iroundpf(fweight);
        accumulator[k] += weight * pixel_value;
        count[k] += weight;
      }
    }
  }
}

#endif  // AOM_AV1_ENCODER_ARM_NEON_HIGHBD_TEMPORAL_FILTER_COMMON_H
