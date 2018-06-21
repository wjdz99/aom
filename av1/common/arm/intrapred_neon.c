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

#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"
#include "av1/common/arm/mem_neon.h"
#include "config/aom_dsp_rtcd.h"

static INLINE void highbd_dc_predictor_neon(uint16_t *dst, ptrdiff_t stride,
                                            int bw, int bh,
                                            const uint16_t *above,
                                            const uint16_t *left, int bd) {
  assert(bw == bh);
  assert(bw >= 4);
  assert(IS_POWER_OF_TWO(bw));
  int expected_dc, sum = 0;
  const int count = bw + bh;
  (void)bd;
  uint32x4_t sum_q = vdupq_n_u32(0);
  uint32x2_t sum_d;
  uint16x8_t tmp_q;
  uint16_t *dst_1;
  if (!(bw & 0x07)) {
    for (int i = 0; i < bw; i += 8) {
      tmp_q = vld1q_u16(above);
      sum_q = vpadalq_u16(sum_q, tmp_q);
      above += 8;
    }
    for (int i = 0; i < bh; i += 8) {
      tmp_q = vld1q_u16(left);
      sum_q = vpadalq_u16(sum_q, tmp_q);
      left += 8;
    }
    sum_d = vadd_u32(vget_low_u32(sum_q), vget_high_u32(sum_q));
    sum = vget_lane_s32(vreinterpret_s32_u64(vpaddl_u32(sum_d)), 0);
    expected_dc = (sum + (count >> 1)) / count;
    const uint16x8_t dc = vdupq_n_u16((uint16_t)expected_dc);
    for (int r = 0; r < bh; r++) {
      dst_1 = dst;
      for (int i = 0; i < bw; i += 8) {
        vst1q_u16(dst_1, dc);
        dst_1 += 8;
      }
      dst += stride;
    }
  } else {
    sum_q = vaddl_u16(vld1_u16(above), vld1_u16(left));
    sum_d = vadd_u32(vget_low_u32(sum_q), vget_high_u32(sum_q));
    sum = vget_lane_s32(vreinterpret_s32_u64(vpaddl_u32(sum_d)), 0);
    expected_dc = (sum + (count >> 1)) / count;
    const uint16x4_t dc = vdup_n_u16((uint16_t)expected_dc);
    for (int r = 0; r < bh; r++) {
      vst1_u16(dst, dc);
      dst += stride;
    }
  }
}

#define intra_pred_highbd_sized(type, width, height)                        \
  void aom_highbd_##type##_predictor_##width##x##height##_neon(             \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,               \
      const uint16_t *left, int bd) {                                       \
    highbd_##type##_predictor_neon(dst, stride, width, height, above, left, \
                                   bd);                                     \
  }

#define intra_pred_square(type)          \
  intra_pred_highbd_sized(type, 4, 4);   \
  intra_pred_highbd_sized(type, 8, 8);   \
  intra_pred_highbd_sized(type, 16, 16); \
  intra_pred_highbd_sized(type, 32, 32); \
  intra_pred_highbd_sized(type, 64, 64);

intra_pred_square(dc)

#undef intra_pred_square
