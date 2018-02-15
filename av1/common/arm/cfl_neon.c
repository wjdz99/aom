/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <arm_neon.h>

#include "./av1_rtcd.h"

#include "av1/common/cfl.h"

static INLINE void subtract_average_neon(int16_t *pred_buf, int width,
                                         int height, int round_offset,
                                         int num_pel_log2) {
  (void)round_offset;
  const int16_t *const end = pred_buf + height * CFL_BUF_LINE;
  const int step = (1 << ((width < 16) + (width == 4))) * CFL_BUF_LINE;

  const int16_t *sum_buf = pred_buf;
  int32x4_t sum_32x4 = { 0, 0, 0, 0 };
  do {
    int16x8_t row;
    if (width == 4) {
      const int16x8_t c0 =
          vcombine_s16(vld1_s16(sum_buf), vld1_s16(sum_buf + CFL_BUF_LINE));
      const int16x8_t c1 = vcombine_s16(vld1_s16(sum_buf + 2 * CFL_BUF_LINE),
                                        vld1_s16(sum_buf + 3 * CFL_BUF_LINE));
      row = vaddq_s16(c0, c1);
    } else if (width == 8) {
      row = vaddq_s16(vld1q_s16(sum_buf), vld1q_s16(sum_buf + CFL_BUF_LINE));
    } else {
      row = vaddq_s16(vld1q_s16(sum_buf), vld1q_s16(sum_buf + 8));
      if (width == 32) {
        sum_32x4 = vpadalq_s16(sum_32x4, vaddq_s16(vld1q_s16(sum_buf + 16),
                                                   vld1q_s16(sum_buf + 24)));
      }
    }
    sum_32x4 = vpadalq_s16(sum_32x4, row);
  } while ((sum_buf += step) < end);

  int32x4_t flip =
      vcombine_s32(vget_high_s32(sum_32x4), vget_low_s32(sum_32x4));

  sum_32x4 = vaddq_s32(sum_32x4, flip);
  sum_32x4 = vaddq_s32(sum_32x4, vrev64q_s32(sum_32x4));

  const int32x4_t shift = vdupq_n_s32(-num_pel_log2);
  int32x4_t avg = vqrshlq_s32(sum_32x4, shift);
  int16x4_t avg_16x4 = vqmovn_s32(avg);
  const int16x8_t avg_16x8 = vcombine_s16(avg_16x4, avg_16x4);

  // TODO(ltrudeau) SWITCH FROM INT TO UINT IN SUM ONLY
  do {
    if (width == 4) {
      vst1_s16(pred_buf, vsub_s16(vld1_s16(pred_buf), avg_16x4));
    } else {
      vst1q_s16(pred_buf, vsubq_s16(vld1q_s16(pred_buf), avg_16x8));
      if (width > 8) {
        vst1q_s16(pred_buf + 8, vsubq_s16(vld1q_s16(pred_buf + 8), avg_16x8));
      }
      if (width == 32) {
        vst1q_s16(pred_buf + 16, vsubq_s16(vld1q_s16(pred_buf + 16), avg_16x8));
        vst1q_s16(pred_buf + 24, vsubq_s16(vld1q_s16(pred_buf + 24), avg_16x8));
      }
    }
  } while ((pred_buf += CFL_BUF_LINE) < end);
}

CFL_SUB_AVG_FN(neon)
