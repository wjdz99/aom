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
#include <assert.h>

#include "aom_dsp/arm/sum_neon.h"

#define MAX_UPSAMPLE_SZ 16

void av1_upsample_intra_edge_high_neon(uint16_t *p, int sz, int bd) {
  if (!sz) return;

  assert(sz <= MAX_UPSAMPLE_SZ);

  uint16_t edge[MAX_UPSAMPLE_SZ + 3];
  const uint16_t *src = edge;

  // Copy p[-1..(sz-1)] and pad out both ends.
  edge[0] = p[-1];
  edge[1] = p[-1];
  memcpy(edge + 2, p, sz * 2);
  edge[sz + 2] = p[sz - 1];
  p[-2] = p[-1];

  uint16x8_t pixel_val_max =
      bd == 8 ? vdupq_n_u16(255)
              : (bd == 10 ? vdupq_n_u16(1023) : vdupq_n_u16(4095));

  uint16_t *dst = p - 1;

  if (bd == 12) {

    do {
      uint16x8_t p0 = vld1q_u16(src);
      uint16x8_t p1 = vld1q_u16(src + 1);
      uint16x8_t p2 = vld1q_u16(src + 2);
      uint16x8_t p3 = vld1q_u16(src + 3);

      uint16x8_t t0 = vaddq_u16(p1, p2);
      uint16x8_t t1 = vaddq_u16(p0, p3);
      uint32x4_t acc0 = vmull_n_u16(vget_low_u16(t0), 9);
      acc0 = vqsubq_u32(acc0, vmovl_u16(vget_low_u16(t1)));
      uint32x4_t acc1 = vmull_n_u16(vget_high_u16(t0), 9);
      acc1 = vqsubq_u32(acc1, vmovl_u16(vget_high_u16(t1)));

      uint16x8x2_t res;
      res.val[0] = vcombine_u16(vrshrn_n_u32(acc0, 4), vrshrn_n_u32(acc1, 4));
      // Clip pixel values at their bit depth's maximum.
      res.val[0] = vminq_u16(res.val[0], pixel_val_max);
      res.val[1] = p2;

      vst2q_u16(dst, res);

      src += 8;
      dst += 16;
      sz -= 8;
    } while (sz > 0);
  } else {  // Bit depth is 8 or 10.

    do {
      uint16x8_t p0 = vld1q_u16(src);
      uint16x8_t p1 = vld1q_u16(src + 1);
      uint16x8_t p2 = vld1q_u16(src + 2);
      uint16x8_t p3 = vld1q_u16(src + 3);

      uint16x8_t t0 = vaddq_u16(p0, p3);
      uint16x8_t t1 = vaddq_u16(p1, p2);
      t1 = vmulq_n_u16(t1, 9);
      t1 = vqsubq_u16(t1, t0);

      uint16x8x2_t res;
      res.val[0] = vrshrq_n_u16(t1, 4);
      // Clip pixel values at their bit depth's maximum.
      res.val[0] = vminq_u16(res.val[0], pixel_val_max);
      res.val[1] = p2;

      vst2q_u16(dst, res);

      src += 8;
      dst += 16;
      sz -= 8;
    } while (sz > 0);
  }
}
