/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <stdio.h>

#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "aom_dsp/aom_dsp_common.h"

#include "./aom_dsp_rtcd.h"

#define MASK_BITS 6

void aom_blend_mask6_c(uint8_t *dst, uint32_t dst_stride, uint8_t *src0,
                       uint32_t src0_stride, uint8_t *src1,
                       uint32_t src1_stride, const uint8_t *mask,
                       uint32_t mask_stride, int h, int w, int subh, int subw) {
  int i, j;

  assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
  assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

  assert(h >= 4);
  assert(w >= 4);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  if (subw == 0 && subh == 0) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 = mask[i * mask_stride + j];
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else if (subw == 1 && subh == 1) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 = ROUND_POWER_OF_TWO(
            mask[(2 * i) * mask_stride + (2 * j)] +
                mask[(2 * i + 1) * mask_stride + (2 * j)] +
                mask[(2 * i) * mask_stride + (2 * j + 1)] +
                mask[(2 * i + 1) * mask_stride + (2 * j + 1)],
            2);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else if (subw == 1 && subh == 0) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 =
            ROUND_POWER_OF_TWO(mask[i * mask_stride + (2 * j)] +
                                   mask[i * mask_stride + (2 * j + 1)],
                               1);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 =
            ROUND_POWER_OF_TWO(mask[(2 * i) * mask_stride + j] +
                                   mask[(2 * i + 1) * mask_stride + j],
                               1);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
void aom_highbd_blend_mask6_c(uint8_t *dst_8, uint32_t dst_stride,
                              uint8_t *src0_8, uint32_t src0_stride,
                              uint8_t *src1_8, uint32_t src1_stride,
                              const uint8_t *mask, uint32_t mask_stride, int h,
                              int w, int subh, int subw, int bd) {
  int i, j;
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);
  uint16_t *src0 = CONVERT_TO_SHORTPTR(src0_8);
  uint16_t *src1 = CONVERT_TO_SHORTPTR(src1_8);

  assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
  assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

  assert(h >= 4);
  assert(w >= 4);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  assert(bd == 8 || bd == 10 || bd == 12);

  if (subw == 0 && subh == 0) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 = mask[i * mask_stride + j];
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else if (subw == 1 && subh == 1) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 = ROUND_POWER_OF_TWO(
            mask[(2 * i) * mask_stride + (2 * j)] +
                mask[(2 * i + 1) * mask_stride + (2 * j)] +
                mask[(2 * i) * mask_stride + (2 * j + 1)] +
                mask[(2 * i + 1) * mask_stride + (2 * j + 1)],
            2);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else if (subw == 1 && subh == 0) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 =
            ROUND_POWER_OF_TWO(mask[i * mask_stride + (2 * j)] +
                                   mask[i * mask_stride + (2 * j + 1)],
                               1);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  } else {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        const int m0 =
            ROUND_POWER_OF_TWO(mask[(2 * i) * mask_stride + j] +
                                   mask[(2 * i + 1) * mask_stride + j],
                               1);
        const int m1 = ((1 << MASK_BITS) - m0);
        dst[i * dst_stride + j] = ROUND_POWER_OF_TWO(
            src0[i * src0_stride + j] * m0 + src1[i * src1_stride + j] * m1,
            MASK_BITS);
      }
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH
