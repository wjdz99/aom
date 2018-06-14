/*
 *
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>
#include <assert.h>

#include "aom/aom_integer.h"
#include "aom_dsp/blend.h"
#include "aom_ports/mem.h"
#include "av1/common/arm/mem_neon.h"
#include "av1/common/blockd.h"

void av1_build_compound_diffwtd_mask_d16_neon(
    uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0,
    int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w,
    ConvolveParams *conv_params, int bd) {
  int round =
      2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1 + (bd - 8);
  uint16x8_t diff_q, tmp0, tmp1;
  uint8x8_t diff_d;
  const CONV_BUF_TYPE *src0_1, *src1_1;
  const int16x8_t dup_round = vdupq_n_s16((int16_t)(-round));
  const uint16x8_t dup_38 = vdupq_n_u16(38);
  const uint16x8_t dup_64 = vdupq_n_u16(64);
  if (!(w & 0x07)) {
    for (int i = 0; i < h; ++i) {
      src0_1 = src0;
      src1_1 = src1;
      for (int j = 0; j < w; j += 8) {
        __builtin_prefetch(src0_1);
        __builtin_prefetch(src1_1);
        diff_q = vabdq_u16(vld1q_u16(src0_1), vld1q_u16(src1_1));
        diff_q = vrshlq_u16(diff_q, dup_round);
        diff_q = vshrq_n_u16(diff_q, DIFF_FACTOR_LOG2);
        diff_q = vminq_u16(vaddq_u16(diff_q, dup_38), dup_64);
        diff_d = vmovn_u16(diff_q);
        if (DIFFWTD_38 == mask_type) {
          vst1_u8(mask, diff_d);
        } else if (DIFFWTD_38_INV == mask_type) {
          diff_d = vsub_u8(vmovn_u16(dup_64), diff_d);
          vst1_u8(mask, diff_d);
        } else {
          assert(0);
        }
        src0_1 += 8;
        src1_1 += 8;
        mask += 8;
      }
      src0 += src0_stride;
      src1 += src1_stride;
    }
  } else if (!(w & 0x03)) {
    for (int i = 0; i < h; i += 2) {
      src0_1 = src0;
      src1_1 = src1;
      __builtin_prefetch(src0_1 + 0 * src0_stride);
      __builtin_prefetch(src0_1 + 1 * src0_stride);
      __builtin_prefetch(src1_1 + 0 * src1_stride);
      __builtin_prefetch(src1_1 + 1 * src1_stride);
      tmp0 = vcombine_u16(vld1_u16(src0_1 + (0 * src0_stride)),
                          vld1_u16(src0_1 + (1 * src0_stride)));
      tmp1 = vcombine_u16(vld1_u16(src1_1 + (0 * src1_stride)),
                          vld1_u16(src1_1 + (1 * src1_stride)));
      diff_q = vabdq_u16(tmp0, tmp1);
      diff_q = vrshlq_u16(diff_q, dup_round);
      diff_q = vshrq_n_u16(diff_q, DIFF_FACTOR_LOG2);
      diff_q = vminq_u16(vaddq_u16(diff_q, dup_38), dup_64);
      diff_d = vmovn_u16(diff_q);
      if (DIFFWTD_38 == mask_type) {
        vst1_u8(mask, diff_d);
      } else if (DIFFWTD_38_INV == mask_type) {
        diff_d = vsub_u8(vmovn_u16(dup_64), diff_d);
        vst1_u8(mask, diff_d);
      } else {
        assert(0);
      }
      src0 += src0_stride * 2;
      src1 += src1_stride * 2;
      mask += w * 2;
    }
  }
}
