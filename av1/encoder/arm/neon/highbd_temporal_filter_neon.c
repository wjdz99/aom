/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>

#include "config/aom_config.h"
#include "config/av1_rtcd.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/temporal_filter.h"
#include "aom_dsp/mathutils.h"
#include "aom_dsp/arm/mem_neon.h"
#include "aom_dsp/arm/sum_neon.h"

// For the squared error buffer, add padding for 4 samples.
#define SSE_STRIDE (BW + 4)

static INLINE void get_squared_error(
    const uint16_t *frame1, const uint32_t stride1, const uint16_t *frame2,
    const uint32_t stride2, const uint32_t block_width,
    const uint32_t block_height, uint32_t *frame_sse,
    const unsigned int dst_stride) {
  uint32_t *dst = frame_sse;

  uint32_t i = 0;
  do {
    uint32_t j = 0;
    do {
      uint16x8_t s = vld1q_u16(frame1 + i * stride1 + j);
      uint16x8_t r = vld1q_u16(frame2 + i * stride2 + j);

      uint16x8_t abs_diff = vabdq_u16(s, r);
      uint32x4_t sse_lo =
          vmull_u16(vget_low_u16(abs_diff), vget_low_u16(abs_diff));
      uint32x4_t sse_hi =
          vmull_u16(vget_high_u16(abs_diff), vget_high_u16(abs_diff));

      vst1q_u32(dst + j, sse_lo);
      vst1q_u32(dst + j + 4, sse_hi);

      j += 8;
    } while (j < block_width);

    dst += dst_stride;
    i++;
  } while (i < block_height);
}

static INLINE void load_and_pad(const uint32_t *src, uint32x4_t *dstvec,
                                const uint32_t col,
                                const uint32_t block_width) {
  uint32x4_t vtmp1 = vld1q_u32(src);
  uint32x4_t vtmp2 = vld1q_u32(src + 4);

  if (col == 0) {
    // For the first column, replicate the first element twice to the left
    const uint32_t lane2 = vgetq_lane_u32(vtmp1, 2);
    vtmp1 = vsetq_lane_u32(lane2, vtmp1, 0);
    vtmp1 = vsetq_lane_u32(lane2, vtmp1, 1);
  } else if (col >= block_width - 4) {
    // For the last column, replicate the last element twice to the right
    const uint32_t lane5 = vgetq_lane_u32(vtmp2, 1);
    vtmp2 = vsetq_lane_u32(lane5, vtmp2, 2);
    vtmp2 = vsetq_lane_u32(lane5, vtmp2, 3);
  }
  dstvec[0] = vtmp1;
  dstvec[1] = vtmp2;
}

static void highbd_apply_temporal_filter(
    const uint16_t *frame, const unsigned int stride,
    const uint32_t block_width, const uint32_t block_height,
    const int *subblock_mses, unsigned int *accumulator, uint16_t *count,
    const uint32_t *frame_sse, const uint32_t *luma_sse_sum,
    const double inv_num_ref_pixels, const double decay_factor,
    const double inv_factor, const double weight_factor, const double *d_factor,
    int tf_wgt_calc_lvl) {
  assert(((block_width == 16) || (block_width == 32)) &&
         ((block_height == 16) || (block_height == 32)));

  uint32_t acc_5x5_neon[BH][BW];

  // Traverse 4 columns at a time - first and last two columns need padding.
  for (uint32_t col = 0; col < block_width; col += 4) {
    uint32x4_t vsrc[5][2];
    const uint32_t *src = frame_sse + col;

    // Load and pad (for first and last two columns) 3 rows from the top.
    for (int i = 2; i < 5; i++) {
      load_and_pad(src, vsrc[i], col, block_width);
      src += SSE_STRIDE;
    }

    // Pad the top 2 rows.
    vsrc[0][0] = vsrc[2][0];
    vsrc[0][1] = vsrc[2][1];
    vsrc[1][0] = vsrc[2][0];
    vsrc[1][1] = vsrc[2][1];

    const uint32x4_t vzero = vdupq_n_u32(0);
    const uint32x4_t _1111 = vdupq_n_u32(1);
    const uint32x4_t _0111 = vextq_u32(vzero, _1111, 3);
    const uint32x4_t _0011 = vextq_u32(vzero, _1111, 2);
    const uint32x4_t _0001 = vextq_u32(vzero, _1111, 1);
    const uint32x4_t _1100 = vextq_u32(_1111, vzero, 2);
    const uint32x4_t _1000 = vextq_u32(_1111, vzero, 3);

    uint32x4_t vsum = vsrc[0][0];
    vsum = vmlaq_u32(vsum, vsrc[0][1], _1000);
    vsum = vmlaq_u32(vsum, vsrc[1][0], _0111);
    vsum = vmlaq_u32(vsum, vsrc[1][1], _1100);
    vsum = vmlaq_u32(vsum, vsrc[2][0], _0011);
    vsum = vmlaq_u32(vsum, vsrc[2][1], _1000);
    vsum = vmlaq_u32(vsum, vsrc[3][0], _0001);
    vsum = vaddq_u32(vsum, vsrc[3][1]);
    uint32_t sum = horizontal_add_u32x4(vsum);

    for (unsigned int row = 0; row < block_height; row++) {
      for (int i = 0; i < 4; i++) {
        acc_5x5_neon[row][col + i] = sum;
      }
    }
  }

  // Perform filtering.
  if (tf_wgt_calc_lvl == 0) {
    for (unsigned int i = 0, k = 0; i < block_height; i++) {
      for (unsigned int j = 0; j < block_width; j++, k++) {
        const int pixel_value = frame[i * stride + j];
        const uint32_t diff_sse = acc_5x5_neon[i][j] + luma_sse_sum[i * BW + j];

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
        const uint32_t diff_sse = acc_5x5_neon[i][j] + luma_sse_sum[i * BW + j];

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

void av1_highbd_apply_temporal_filter_neon(
    const YV12_BUFFER_CONFIG *frame_to_filter, const MACROBLOCKD *mbd,
    const BLOCK_SIZE block_size, const int mb_row, const int mb_col,
    const int num_planes, const double *noise_levels, const MV *subblock_mvs,
    const int *subblock_mses, const int q_factor, const int filter_strength,
    int tf_wgt_calc_lvl, const uint8_t *pred8, uint32_t *accum,
    uint16_t *count) {
  const int is_high_bitdepth = frame_to_filter->flags & YV12_FLAG_HIGHBITDEPTH;
  assert(TF_WINDOW_LENGTH == 5 && "Only support window length 5 with Neon!");
  assert(num_planes >= 1 && num_planes <= MAX_MB_PLANE);
  (void)is_high_bitdepth;
  assert(is_high_bitdepth);

  // Block information.
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  // Frame information.
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  const int min_frame_size = AOMMIN(frame_height, frame_width);
  // Variables to simplify combined error calculation.
  const double inv_factor = 1.0 / ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) *
                                   TF_SEARCH_ERROR_NORM_WEIGHT);
  const double weight_factor =
      (double)TF_WINDOW_BLOCK_BALANCE_WEIGHT * inv_factor;
  // Adjust filtering based on q.
  // Larger q -> stronger filtering -> larger weight.
  // Smaller q -> weaker filtering -> smaller weight.
  double q_decay = pow((double)q_factor / TF_Q_DECAY_THRESHOLD, 2);
  q_decay = CLIP(q_decay, 1e-5, 1);
  if (q_factor >= TF_QINDEX_CUTOFF) {
    // Max q_factor is 255, therefore the upper bound of q_decay is 8.
    // We do not need a clip here.
    q_decay = 0.5 * pow((double)q_factor / 64, 2);
  }
  // Smaller strength -> smaller filtering weight.
  double s_decay = pow((double)filter_strength / TF_STRENGTH_THRESHOLD, 2);
  s_decay = CLIP(s_decay, 1e-5, 1);
  double d_factor[4] = { 0 };
  uint32_t frame_sse[SSE_STRIDE * BH] = { 0 };
  uint32_t luma_sse_sum[BW * BH] = { 0 };
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);

  for (int subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
    // Larger motion vector -> smaller filtering weight.
    const MV mv = subblock_mvs[subblock_idx];
    const double distance = sqrt(pow(mv.row, 2) + pow(mv.col, 2));
    double distance_threshold = min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD;
    distance_threshold = AOMMAX(distance_threshold, 1);
    d_factor[subblock_idx] = distance / distance_threshold;
    d_factor[subblock_idx] = AOMMAX(d_factor[subblock_idx], 1);
  }

  // Handle planes in sequence.
  int plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const uint32_t plane_h = mb_height >> mbd->plane[plane].subsampling_y;
    const uint32_t plane_w = mb_width >> mbd->plane[plane].subsampling_x;
    const uint32_t frame_stride =
        frame_to_filter->strides[plane == AOM_PLANE_Y ? 0 : 1];
    const int frame_offset = mb_row * plane_h * frame_stride + mb_col * plane_w;

    const uint16_t *ref =
        CONVERT_TO_SHORTPTR(frame_to_filter->buffers[plane]) + frame_offset;
    const int ss_x_shift =
        mbd->plane[plane].subsampling_x - mbd->plane[AOM_PLANE_Y].subsampling_x;
    const int ss_y_shift =
        mbd->plane[plane].subsampling_y - mbd->plane[AOM_PLANE_Y].subsampling_y;
    const int num_ref_pixels = TF_WINDOW_LENGTH * TF_WINDOW_LENGTH +
                               ((plane) ? (1 << (ss_x_shift + ss_y_shift)) : 0);
    const double inv_num_ref_pixels = 1.0 / num_ref_pixels;
    // Larger noise -> larger filtering weight.
    const double n_decay = 0.5 + log(2 * noise_levels[plane] + 5.0);
    // Decay factors for non-local mean approach.
    const double decay_factor = 1 / (n_decay * q_decay * s_decay);

    // Filter U-plane and V-plane using Y-plane. This is because motion
    // search is only done on Y-plane, so the information from Y-plane
    // will be more accurate. The luma sse sum is reused in both chroma
    // planes.
    if (plane == AOM_PLANE_U) {
      for (unsigned int i = 0; i < plane_h; i++) {
        for (unsigned int j = 0; j < plane_w; j++) {
          for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
            for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
              const int yy = (i << ss_y_shift) + ii;  // Y-coord on Y-plane.
              const int xx = (j << ss_x_shift) + jj;  // X-coord on Y-plane.
              luma_sse_sum[i * BW + j] += frame_sse[yy * SSE_STRIDE + xx + 2];
            }
          }
        }
      }
    }

    get_squared_error(ref, frame_stride, pred + plane_offset, plane_w, plane_w,
                      plane_h, frame_sse, SSE_STRIDE);

    highbd_apply_temporal_filter(pred + plane_offset, plane_w, plane_w, plane_h,
                                 subblock_mses, accum + plane_offset,
                                 count + plane_offset, frame_sse, luma_sse_sum,
                                 inv_num_ref_pixels, decay_factor, inv_factor,
                                 weight_factor, d_factor, tf_wgt_calc_lvl);

    plane_offset += plane_h * plane_w;
  }
}
