/*
 * Copyright (c) 2014 The WebM project authors. All Rights Reserved.
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
#include <string.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/aom_filter.h"
#include "aom_dsp/arm/mem_neon.h"
#include "aom_dsp/arm/transpose_neon.h"
#include "aom_ports/mem.h"

static INLINE int32x4_t highbd_convolve8_4_s32(
    const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
    const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
    const int16x4_t s6, const int16x4_t s7, const int16x8_t y_filter,
    const int32x4_t offset) {
  const int16x4_t y_filter_lo = vget_low_s16(y_filter);
  const int16x4_t y_filter_hi = vget_high_s16(y_filter);

  int32x4_t sum = vmlal_lane_s16(offset, s0, y_filter_lo, 0);
  sum = vmlal_lane_s16(sum, s1, y_filter_lo, 1);
  sum = vmlal_lane_s16(sum, s2, y_filter_lo, 2);
  sum = vmlal_lane_s16(sum, s3, y_filter_lo, 3);
  sum = vmlal_lane_s16(sum, s4, y_filter_hi, 0);
  sum = vmlal_lane_s16(sum, s5, y_filter_hi, 1);
  sum = vmlal_lane_s16(sum, s6, y_filter_hi, 2);
  sum = vmlal_lane_s16(sum, s7, y_filter_hi, 3);

  return sum;
}

static INLINE uint16x4_t highbd_convolve8_4_s32_s16(
    const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
    const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
    const int16x4_t s6, const int16x4_t s7, const int16x8_t y_filter,
    const int32x4_t offset) {
  int32x4_t sum =
      highbd_convolve8_4_s32(s0, s1, s2, s3, s4, s5, s6, s7, y_filter, offset);

  return vqrshrun_n_s32(sum, FILTER_BITS);
}

static INLINE int32x4_t highbd_convolve8_horiz4_s32(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x8_t x_filter_0_7, const int32x4_t offset) {
  const int16x4_t s0_lo = vget_low_s16(s0);
  const int16x4_t s1_lo = vget_low_s16(s1);
  const int16x4_t s2_lo = vget_low_s16(s2);
  const int16x4_t s3_lo = vget_low_s16(s3);
  const int16x4_t s4_lo = vget_high_s16(s0);
  const int16x4_t s5_lo = vget_high_s16(s1);
  const int16x4_t s6_lo = vget_high_s16(s2);
  const int16x4_t s7_lo = vget_high_s16(s3);

  return highbd_convolve8_4_s32(s0_lo, s1_lo, s2_lo, s3_lo, s4_lo, s5_lo, s6_lo,
                                s7_lo, x_filter_0_7, offset);
}

static INLINE uint16x4_t highbd_convolve8_horiz4_s32_s16(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x8_t x_filter_0_7, const int32x4_t offset) {
  int32x4_t sum =
      highbd_convolve8_horiz4_s32(s0, s1, s2, s3, x_filter_0_7, offset);

  return vqrshrun_n_s32(sum, FILTER_BITS);
}

static void highbd_convolve_horiz_neon(const uint16_t *src_ptr,
                                       ptrdiff_t src_stride, uint16_t *dst_ptr,
                                       ptrdiff_t dst_stride,
                                       const int16_t *x_filter_ptr,
                                       int x_step_q4, int w, int h, int bd) {
  const int32x4_t idx = { 0, 1, 2, 3 };
  const uint16x4_t max = vdup_n_u16((1 << bd) - 1);

  int height = h;
  int16x8_t s0, s1, s2, s3;
  uint16x4_t d0;

  do {
    int width = w;
    int x_q4 = 0;
    uint16_t *d = dst_ptr;
    const uint16_t *s = src_ptr;

    do {
      // Load 4 src vectors at a time, they might be the same, but we have to
      // calculate the indices anyway. Doing it in SIMD and then storing the
      // indices is faster than having to calculate the expression
      // &src_ptr[((x_q4 + i*x_step_q4) >> SUBPEL_BITS)] 4 times
      // Ideally this should be a gather using the indices, but NEON does not
      // have that, so have to emulate
      const int32x4_t xq4_idx = vmlaq_n_s32(vdupq_n_s32(x_q4), idx, x_step_q4);
      // We have to multiply x2 to get the actual pointer as sizeof(uint16_t)
      // = 2
      const int32x4_t src_idx =
          vshlq_n_s32(vshrq_n_s32(xq4_idx, SUBPEL_BITS), 1);

#if AOM_ARCH_AARCH64
      uint64x2_t tmp4[2];
      tmp4[0] = vreinterpretq_u64_s64(
          vaddw_s32(vdupq_n_s64((const int64_t)s), vget_low_s32(src_idx)));
      tmp4[1] = vreinterpretq_u64_s64(
          vaddw_s32(vdupq_n_s64((const int64_t)s), vget_high_s32(src_idx)));
      int16_t *src4_ptr[4];
      uint64_t *tmp_ptr = (uint64_t *)&src4_ptr;
      vst1q_u64(tmp_ptr, tmp4[0]);
      vst1q_u64(tmp_ptr + 2, tmp4[1]);
#else
      uint32x4_t tmp4;
      tmp4 = vreinterpretq_u32_s32(
          vaddq_s32(vdupq_n_s32((const int32_t)s), src_idx));
      int16_t *src4_ptr[4];
      uint32_t *tmp_ptr = (uint32_t *)&src4_ptr;
      vst1q_u32(tmp_ptr, tmp4);
#endif  // AOM_ARCH_AARCH64
      // Load source
      s0 = vld1q_s16(src4_ptr[0]);
      s1 = vld1q_s16(src4_ptr[1]);
      s2 = vld1q_s16(src4_ptr[2]);
      s3 = vld1q_s16(src4_ptr[3]);

      // Actually load the filters
      const int16x8_t x_filter = vld1q_s16(x_filter_ptr);

      d0 = highbd_convolve8_horiz4_s32_s16(s0, s1, s2, s3, x_filter,
                                           vdupq_n_s32(0));

      d0 = vmin_u16(d0, max);
      vst1_u16(d, d0);

      x_q4 += 4 * x_step_q4;
      d += 4;
      width -= 4;
    } while (width > 0);

    src_ptr += src_stride;
    dst_ptr += dst_stride;
    height--;
  } while (height > 0);
}

static void highbd_convolve_vert_neon(const uint16_t *src_ptr,
                                      ptrdiff_t src_stride, uint16_t *dst_ptr,
                                      ptrdiff_t dst_stride,
                                      const int16_t *y_filter_ptr,
                                      int y_step_q4, int w, int h, int bd) {
  const int32x4_t idx = { 0, 1, 2, 3 };
  const uint16x4_t max = vdup_n_u16((1 << bd) - 1);

  int width = w;
  int16x4_t s0, s1, s2, s3, s4, s5, s6, s7;
  uint16x4_t d0;

  do {
    int height = h;
    int y_q4 = 0;
    uint16_t *d = dst_ptr;
    const uint16_t *s = src_ptr;

    do {
      // Load 4 src vectors at a time, they might be the same, but we have to
      // calculate the indices anyway. Doing it in SIMD and then storing the
      // indices is faster than having to calculate the expression
      // &src_ptr[((x_q4 + i*x_step_q4) >> SUBPEL_BITS)] 4 times
      // Ideally this should be a gather using the indices, but NEON does not
      // have that, so have to emulate
      const int32x4_t yq4_idx = vmlaq_n_s32(vdupq_n_s32(y_q4), idx, y_step_q4);
      // We have to multiply x2 to get the actual pointer as sizeof(uint16_t)
      // = 2
      const int32x4_t src_idx =
          vshlq_n_s32(vshrq_n_s32(yq4_idx, SUBPEL_BITS), 1);
#if AOM_ARCH_AARCH64
      uint64x2_t tmp4[2];
      tmp4[0] = vreinterpretq_u64_s64(
          vaddw_s32(vdupq_n_s64((const int64_t)s), vget_low_s32(src_idx)));
      tmp4[1] = vreinterpretq_u64_s64(
          vaddw_s32(vdupq_n_s64((const int64_t)s), vget_high_s32(src_idx)));
      const int16_t *src4_ptr[4];
      uint64_t *tmp_ptr = (uint64_t *)&src4_ptr;
      vst1q_u64(tmp_ptr, tmp4[0]);
      vst1q_u64(tmp_ptr + 2, tmp4[1]);
#else
      uint32x4_t tmp4;
      tmp4 = vreinterpretq_u32_s32(
          vaddq_s32(vdupq_n_s32((const int32_t)s), src_idx));
      int16_t *src4_ptr[4];
      uint32_t *tmp_ptr = (uint32_t *)&src4_ptr;
      vst1q_u32(tmp_ptr, tmp4);
#endif  // AOM_ARCH_AARCH64

      // Load source
      load_s16_4x8(src4_ptr[0], src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6,
                   &s7);

      // Actually load the filters
      const int16x8_t y_filter = vld1q_s16(y_filter_ptr);

      // Run the convolution
      d0 = highbd_convolve8_4_s32_s16(s0, s1, s2, s3, s4, s5, s6, s7, y_filter,
                                      vdupq_n_s32(0));
      d0 = vmin_u16(d0, max);
      vst1_u16(d, d0);

      s += src_stride;
      d += dst_stride;
      height--;
    } while (height > 0);

    y_q4 += 4 * y_step_q4;
    src_ptr += 4;
    dst_ptr += 4;
    width -= 4;
  } while (width > 0);
}

void aom_highbd_convolve8_horiz_neon(const uint8_t *src8, ptrdiff_t src_stride,
                                     uint8_t *dst8, ptrdiff_t dst_stride,
                                     const int16_t *filter_x, int x_step_q4,
                                     const int16_t *filter_y, int y_step_q4,
                                     int w, int h, int bd) {
  (void)filter_y;
  (void)y_step_q4;

  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);

  src -= SUBPEL_TAPS / 2 - 1;
  highbd_convolve_horiz_neon(src, src_stride, dst, dst_stride, filter_x,
                             x_step_q4, w, h, bd);
}

void aom_highbd_convolve8_vert_neon(const uint8_t *src8, ptrdiff_t src_stride,
                                    uint8_t *dst8, ptrdiff_t dst_stride,
                                    const int16_t *filter_x, int x_step_q4,
                                    const int16_t *filter_y, int y_step_q4,
                                    int w, int h, int bd) {
  (void)filter_x;
  (void)x_step_q4;

  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);

  src -= (SUBPEL_TAPS / 2 - 1) * src_stride;
  highbd_convolve_vert_neon(src, src_stride, dst, dst_stride, filter_y,
                            y_step_q4, w, h, bd);
}

void aom_highbd_convolve_copy_neon(const uint16_t *src, ptrdiff_t src_stride,
                                   uint16_t *dst, ptrdiff_t dst_stride, int w,
                                   int h) {
  for (int y = 0; y < h; ++y) {
    int x;
    for (x = 0; x + 8 <= w; x += 8) {
      vst1q_u16(dst + x, vld1q_u16(src + x));
    }
    for (; x + 4 <= w; x += 4) {
      vst1_u16(dst + x, vld1_u16(src + x));
    }
    for (; x + 2 <= w; x += 2) {
      vst1_lane_u32((uint32_t *)(dst + x),
                    vreinterpret_u32_u16(vld1_u16(src + x)), 0);
    }
    src += src_stride;
    dst += dst_stride;
  }
}
