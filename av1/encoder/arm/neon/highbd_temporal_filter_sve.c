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

#include <arm_neon.h>
#include <arm_neon_sve_bridge.h>
#include <assert.h>

#include "config/aom_config.h"
#include "config/av1_rtcd.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/temporal_filter.h"
#include "av1/encoder/arm/neon/highbd_temporal_filter_common.h"
#include "aom_dsp/arm/dot_sve.h"
#include "aom_dsp/arm/mem_neon.h"
#include "aom_dsp/arm/sum_neon.h"

// For the squared error buffer, add padding for 4 samples.
#define SSE_STRIDE (BW + 4)

static INLINE void get_abs_diff(const uint16_t *frame1, const uint32_t stride1,
                                const uint16_t *frame2, const uint32_t stride2,
                                const uint32_t block_width,
                                const uint32_t block_height,
                                uint16_t *frame_abs_diff,
                                const unsigned int frame_abs_diff_stride) {
  uint16_t *dst = frame_abs_diff;

  uint32_t i = 0;
  do {
    uint32_t j = 0;
    do {
      uint16x8_t s = vld1q_u16(frame1 + i * stride1 + j);
      uint16x8_t r = vld1q_u16(frame2 + i * stride2 + j);
      uint16x8_t abs_diff = vabdq_u16(s, r);
      vst1q_u16(dst + j + 2, abs_diff);
      j += 8;
    } while (j < block_width);

    dst += frame_abs_diff_stride;
  } while (++i < block_height);
}

// Table used to pad the first and last columns and apply the sliding window.
DECLARE_ALIGNED(16, static const uint16_t, kLoadPad[4][16]) = {
  { 2, 2, 2, 3, 4, 255, 255, 255, 0, 1, 2, 3, 4, 255, 255, 255 },
  { 255, 2, 2, 3, 4, 5, 255, 255, 255, 1, 2, 3, 4, 5, 255, 255 },
  { 255, 255, 2, 3, 4, 5, 6, 255, 255, 255, 2, 3, 4, 5, 5, 255 },
  { 255, 255, 255, 3, 4, 5, 6, 7, 255, 255, 255, 3, 4, 5, 5, 5 }
};

static INLINE uint16x8_t aom_tbl_u16(uint16x8_t src, uint16x8_t table) {
  return svget_neonq_u16(svtbl_u16(svset_neonq_u16(svundef_u16(), src),
                                   svset_neonq_u16(svundef_u16(), table)));
}

static void apply_temporal_filter(
    const uint16_t *frame, const unsigned int stride,
    const uint32_t block_width, const uint32_t block_height,
    const int *subblock_mses, unsigned int *accumulator, uint16_t *count,
    const uint16_t *frame_abs_diff, const uint32_t *luma_sse_sum,
    const double inv_num_ref_pixels, const double decay_factor,
    const double inv_factor, const double weight_factor, const double *d_factor,
    int tf_wgt_calc_lvl, int bd) {
  assert((block_width == 16 || block_width == 32) &&
         (block_height == 16 || block_height == 32));

  uint32_t acc_5x5_neon[BH][BW];
  const uint16x8x2_t pad_tbl0 = vld1q_u16_x2(kLoadPad[0]);
  const uint16x8x2_t pad_tbl1 = vld1q_u16_x2(kLoadPad[1]);
  const uint16x8x2_t pad_tbl2 = vld1q_u16_x2(kLoadPad[2]);
  const uint16x8x2_t pad_tbl3 = vld1q_u16_x2(kLoadPad[3]);

  // For columns that don't need to be padded we can just use predicated loads
  // to apply the sliding window directly. Construct one predicate per line of
  // the sliding window (similar to kLoadPad declared above) using pre-defined
  // patterns. SV_VLn means a predicate pattern with the first n elements true.
  svbool_t sliding_window_mask0 = svptrue_pat_b16(SV_VL5);
  svbool_t p6 = svptrue_pat_b16(SV_VL6);
  svbool_t p1 = svptrue_pat_b16(SV_VL1);
  svbool_t sliding_window_mask1 = sveor_b_z(p6, p6, p1);
  svbool_t p7 = svptrue_pat_b16(SV_VL7);
  svbool_t p2 = svptrue_pat_b16(SV_VL2);
  svbool_t sliding_window_mask2 = sveor_b_z(p7, p7, p2);
  svbool_t p8 = svptrue_pat_b16(SV_VL8);
  svbool_t p3 = svptrue_pat_b16(SV_VL3);
  svbool_t sliding_window_mask3 = sveor_b_z(p8, p8, p3);

  // Traverse 4 columns at a time - first and last two columns need padding.
  for (uint32_t col = 0; col < block_width; col += 4) {
    uint16x8_t vsrc[5][4];
    const uint16_t *src = frame_abs_diff + col;

    // Load, pad (for first and last two columns) and mask 3 rows from the top.
    for (int i = 2; i < 5; i++) {
      if (col == 0) {
        uint16x8_t s = vld1q_u16(src);
        vsrc[i][0] = aom_tbl_u16(s, pad_tbl0.val[0]);
        vsrc[i][1] = aom_tbl_u16(s, pad_tbl1.val[0]);
        vsrc[i][2] = aom_tbl_u16(s, pad_tbl2.val[0]);
        vsrc[i][3] = aom_tbl_u16(s, pad_tbl3.val[0]);
      } else if (col >= block_width - 4) {
        uint16x8_t s = vld1q_u16(src);
        vsrc[i][0] = aom_tbl_u16(s, pad_tbl0.val[1]);
        vsrc[i][1] = aom_tbl_u16(s, pad_tbl1.val[1]);
        vsrc[i][2] = aom_tbl_u16(s, pad_tbl2.val[1]);
        vsrc[i][3] = aom_tbl_u16(s, pad_tbl3.val[1]);
      } else {
        vsrc[i][0] = svget_neonq_u16(svld1_u16(sliding_window_mask0, src));
        vsrc[i][1] = svget_neonq_u16(svld1_u16(sliding_window_mask1, src));
        vsrc[i][2] = svget_neonq_u16(svld1_u16(sliding_window_mask2, src));
        vsrc[i][3] = svget_neonq_u16(svld1_u16(sliding_window_mask3, src));
      }
      src += SSE_STRIDE;
    }

    // Pad the top 2 rows.
    vsrc[0][0] = vsrc[2][0];
    vsrc[0][1] = vsrc[2][1];
    vsrc[0][2] = vsrc[2][2];
    vsrc[0][3] = vsrc[2][3];
    vsrc[1][0] = vsrc[2][0];
    vsrc[1][1] = vsrc[2][1];
    vsrc[1][2] = vsrc[2][2];
    vsrc[1][3] = vsrc[2][3];

    for (uint32_t row = 0; row < block_height; ++row) {
      uint64x2_t sum[4] = { vdupq_n_u64(0), vdupq_n_u64(0), vdupq_n_u64(0),
                            vdupq_n_u64(0) };

      sum[0] = aom_udotq_u16(sum[0], vsrc[0][0], vsrc[0][0]);
      sum[0] = aom_udotq_u16(sum[0], vsrc[1][0], vsrc[1][0]);
      sum[0] = aom_udotq_u16(sum[0], vsrc[2][0], vsrc[2][0]);
      sum[0] = aom_udotq_u16(sum[0], vsrc[3][0], vsrc[3][0]);
      sum[0] = aom_udotq_u16(sum[0], vsrc[4][0], vsrc[4][0]);

      sum[1] = aom_udotq_u16(sum[1], vsrc[0][1], vsrc[0][1]);
      sum[1] = aom_udotq_u16(sum[1], vsrc[1][1], vsrc[1][1]);
      sum[1] = aom_udotq_u16(sum[1], vsrc[2][1], vsrc[2][1]);
      sum[1] = aom_udotq_u16(sum[1], vsrc[3][1], vsrc[3][1]);
      sum[1] = aom_udotq_u16(sum[1], vsrc[4][1], vsrc[4][1]);

      sum[2] = aom_udotq_u16(sum[2], vsrc[0][2], vsrc[0][2]);
      sum[2] = aom_udotq_u16(sum[2], vsrc[1][2], vsrc[1][2]);
      sum[2] = aom_udotq_u16(sum[2], vsrc[2][2], vsrc[2][2]);
      sum[2] = aom_udotq_u16(sum[2], vsrc[3][2], vsrc[3][2]);
      sum[2] = aom_udotq_u16(sum[2], vsrc[4][2], vsrc[4][2]);

      sum[3] = aom_udotq_u16(sum[3], vsrc[0][3], vsrc[0][3]);
      sum[3] = aom_udotq_u16(sum[3], vsrc[1][3], vsrc[1][3]);
      sum[3] = aom_udotq_u16(sum[3], vsrc[2][3], vsrc[2][3]);
      sum[3] = aom_udotq_u16(sum[3], vsrc[3][3], vsrc[3][3]);
      sum[3] = aom_udotq_u16(sum[3], vsrc[4][3], vsrc[4][3]);

      uint32x4_t sum_0123 = vcombine_u32(vmovn_u64(vpaddq_u64(sum[0], sum[1])),
                                         vmovn_u64(vpaddq_u64(sum[2], sum[3])));

      // Narrow back to 32 bit to store
      vst1q_u32(&acc_5x5_neon[row][col], sum_0123);

      // Push all rows in the sliding window up one.
      for (int i = 0; i < 4; i++) {
        vsrc[i][0] = vsrc[i + 1][0];
        vsrc[i][1] = vsrc[i + 1][1];
        vsrc[i][2] = vsrc[i + 1][2];
        vsrc[i][3] = vsrc[i + 1][3];
      }

      if (row <= block_height - 4) {
        if (col == 0) {
          uint16x8_t s = vld1q_u16(src);
          vsrc[4][0] = aom_tbl_u16(s, pad_tbl0.val[0]);
          vsrc[4][1] = aom_tbl_u16(s, pad_tbl1.val[0]);
          vsrc[4][2] = aom_tbl_u16(s, pad_tbl2.val[0]);
          vsrc[4][3] = aom_tbl_u16(s, pad_tbl3.val[0]);
        } else if (col >= block_width - 4) {
          uint16x8_t s = vld1q_u16(src);
          vsrc[4][0] = aom_tbl_u16(s, pad_tbl0.val[1]);
          vsrc[4][1] = aom_tbl_u16(s, pad_tbl1.val[1]);
          vsrc[4][2] = aom_tbl_u16(s, pad_tbl2.val[1]);
          vsrc[4][3] = aom_tbl_u16(s, pad_tbl3.val[1]);
        } else {
          vsrc[4][0] = svget_neonq_u16(svld1_u16(sliding_window_mask0, src));
          vsrc[4][1] = svget_neonq_u16(svld1_u16(sliding_window_mask1, src));
          vsrc[4][2] = svget_neonq_u16(svld1_u16(sliding_window_mask2, src));
          vsrc[4][3] = svget_neonq_u16(svld1_u16(sliding_window_mask3, src));
        }
        src += SSE_STRIDE;
      } else {
        // Pad the bottom 2 rows.
        vsrc[4][0] = vsrc[3][0];
        vsrc[4][1] = vsrc[3][1];
        vsrc[4][2] = vsrc[3][2];
        vsrc[4][3] = vsrc[3][3];
      }
    }
  }

  // Perform filtering.
  apply_temporal_filtering(
      frame, stride, acc_5x5_neon, luma_sse_sum, inv_num_ref_pixels,
      block_height, block_width, subblock_mses, weight_factor, inv_factor,
      decay_factor, d_factor, accumulator, count, tf_wgt_calc_lvl, bd);
}

void av1_highbd_apply_temporal_filter_sve(
    const YV12_BUFFER_CONFIG *frame_to_filter, const MACROBLOCKD *mbd,
    const BLOCK_SIZE block_size, const int mb_row, const int mb_col,
    const int num_planes, const double *noise_levels, const MV *subblock_mvs,
    const int *subblock_mses, const int q_factor, const int filter_strength,
    int tf_wgt_calc_lvl, const uint8_t *pred, uint32_t *accum,
    uint16_t *count) {
  assert(block_size == BLOCK_32X32 && "Only support 32x32 block with Neon!");
  assert(TF_WINDOW_LENGTH == 5 && "Only support window length 5 with Neon!");
  assert(num_planes >= 1 && num_planes <= MAX_MB_PLANE);

  const uint16_t *pred16 = CONVERT_TO_SHORTPTR(pred);
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
  uint16_t frame_abs_diff[SSE_STRIDE * BH] = { 0 };
  uint32_t luma_sse_sum[BW * BH] = { 0 };

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
              luma_sse_sum[i * BW + j] +=
                  (frame_abs_diff[yy * SSE_STRIDE + xx + 2] *
                   frame_abs_diff[yy * SSE_STRIDE + xx + 2]);
            }
          }
        }
      }
    }

    get_abs_diff(ref, frame_stride, pred16 + plane_offset, plane_w, plane_w,
                 plane_h, frame_abs_diff, SSE_STRIDE);

    apply_temporal_filter(pred16 + plane_offset, plane_w, plane_w, plane_h,
                          subblock_mses, accum + plane_offset,
                          count + plane_offset, frame_abs_diff, luma_sse_sum,
                          inv_num_ref_pixels, decay_factor, inv_factor,
                          weight_factor, d_factor, tf_wgt_calc_lvl, mbd->bd);

    plane_offset += plane_h * plane_w;
  }
}
