/*
 *  Copyright (c) 2020, Alliance for Open Media. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <arm_neon.h>

#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/aom_filter.h"
#include "aom_ports/mem.h"
#include "av1/common/arm/convolve_neon.h"
#include "av1/common/arm/mem_neon.h"
#include "av1/common/arm/transpose_neon.h"

static const InterpKernel *get_filter_base(const int16_t *filter) {
  // NOTE: This assumes that the filter table is 256-byte aligned.
  return (const InterpKernel *)(((intptr_t)filter) & ~((intptr_t)0xFF));
}

static int get_filter_offset(const int16_t *f, const InterpKernel *base) {
  return (int)((const InterpKernel *)(intptr_t)f - base);
}

void aom_convolve8_horiz_neon(const uint8_t *src, ptrdiff_t src_stride,
                              uint8_t *dst, ptrdiff_t dst_stride,
                              const int16_t *filter_x, int x_step_q4,
                              const int16_t *filter_y, int y_step_q4, int w,
                              int h) {
  const InterpKernel *const filters_x = get_filter_base(filter_x);
  const int x0_q4 = get_filter_offset(filter_x, filters_x);

  uint8x8_t t0, t1, t2, t3;

  assert(!((intptr_t)dst & 3));
  assert(!(dst_stride & 3));
  assert(x_step_q4 == 16);

  (void)y_step_q4;
  (void)filter_y;

  src -= 3;

  if (h == 4) {
    uint8x8_t d01, d23;
    int16x4_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, d0, d1, d2, d3;
    int16x8_t tt0, tt1, tt2, tt3;

    __builtin_prefetch(src + 0 * src_stride);
    __builtin_prefetch(src + 1 * src_stride);
    __builtin_prefetch(src + 2 * src_stride);
    __builtin_prefetch(src + 3 * src_stride);
    load_u8_8x4(src, src_stride, &t0, &t1, &t2, &t3);
    transpose_u8_8x4(&t0, &t1, &t2, &t3);
    tt0 = vreinterpretq_s16_u16(vmovl_u8(t0));
    tt1 = vreinterpretq_s16_u16(vmovl_u8(t1));
    tt2 = vreinterpretq_s16_u16(vmovl_u8(t2));
    tt3 = vreinterpretq_s16_u16(vmovl_u8(t3));
    s0 = vget_low_s16(tt0);
    s1 = vget_low_s16(tt1);
    s2 = vget_low_s16(tt2);
    s3 = vget_low_s16(tt3);
    s4 = vget_high_s16(tt0);
    s5 = vget_high_s16(tt1);
    s6 = vget_high_s16(tt2);
    __builtin_prefetch(dst + 0 * dst_stride);
    __builtin_prefetch(dst + 1 * dst_stride);
    __builtin_prefetch(dst + 2 * dst_stride);
    __builtin_prefetch(dst + 3 * dst_stride);
    src += 7;

    do {
      load_u8_8x4(src, src_stride, &t0, &t1, &t2, &t3);
      transpose_u8_8x4(&t0, &t1, &t2, &t3);
      tt0 = vreinterpretq_s16_u16(vmovl_u8(t0));
      tt1 = vreinterpretq_s16_u16(vmovl_u8(t1));
      tt2 = vreinterpretq_s16_u16(vmovl_u8(t2));
      tt3 = vreinterpretq_s16_u16(vmovl_u8(t3));
      s7 = vget_low_s16(tt0);
      s8 = vget_low_s16(tt1);
      s9 = vget_low_s16(tt2);
      s10 = vget_low_s16(tt3);

      d0 = convolve8_4x4(s0, s1, s2, s3, s4, s5, s6, s7, filter_x);
      d1 = convolve8_4x4(s1, s2, s3, s4, s5, s6, s7, s8, filter_x);
      d2 = convolve8_4x4(s2, s3, s4, s5, s6, s7, s8, s9, filter_x);
      d3 = convolve8_4x4(s3, s4, s5, s6, s7, s8, s9, s10, filter_x);

      d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), 7);
      d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), 7);
      transpose_u8_4x4(&d01, &d23);

      vst1_lane_u32((uint32_t *)(dst + 0 * dst_stride),
                    vreinterpret_u32_u8(d01), 0);
      vst1_lane_u32((uint32_t *)(dst + 1 * dst_stride),
                    vreinterpret_u32_u8(d23), 0);
      vst1_lane_u32((uint32_t *)(dst + 2 * dst_stride),
                    vreinterpret_u32_u8(d01), 1);
      vst1_lane_u32((uint32_t *)(dst + 3 * dst_stride),
                    vreinterpret_u32_u8(d23), 1);

      s0 = s4;
      s1 = s5;
      s2 = s6;
      s3 = s7;
      s4 = s8;
      s5 = s9;
      s6 = s10;
      src += 4;
      dst += 4;
      w -= 4;
    } while (w > 0);
  } else {
    int width;
    const uint8_t *s;
    uint8x8_t t4, t5, t6, t7;
    int16x8_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10;

    if (w == 4) {
      do {
        load_u8_8x8(src, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
        transpose_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
        s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
        s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
        s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
        s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
        s4 = vreinterpretq_s16_u16(vmovl_u8(t4));
        s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
        s6 = vreinterpretq_s16_u16(vmovl_u8(t6));

        load_u8_8x8(src + 7, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6,
                    &t7);
        src += 8 * src_stride;
        __builtin_prefetch(dst + 0 * dst_stride);
        __builtin_prefetch(dst + 1 * dst_stride);
        __builtin_prefetch(dst + 2 * dst_stride);
        __builtin_prefetch(dst + 3 * dst_stride);
        __builtin_prefetch(dst + 4 * dst_stride);
        __builtin_prefetch(dst + 5 * dst_stride);
        __builtin_prefetch(dst + 6 * dst_stride);
        __builtin_prefetch(dst + 7 * dst_stride);
        transpose_u8_4x8(&t0, &t1, &t2, &t3, t4, t5, t6, t7);
        s7 = vreinterpretq_s16_u16(vmovl_u8(t0));
        s8 = vreinterpretq_s16_u16(vmovl_u8(t1));
        s9 = vreinterpretq_s16_u16(vmovl_u8(t2));
        s10 = vreinterpretq_s16_u16(vmovl_u8(t3));

        __builtin_prefetch(src + 0 * src_stride);
        __builtin_prefetch(src + 1 * src_stride);
        __builtin_prefetch(src + 2 * src_stride);
        __builtin_prefetch(src + 3 * src_stride);
        __builtin_prefetch(src + 4 * src_stride);
        __builtin_prefetch(src + 5 * src_stride);
        __builtin_prefetch(src + 6 * src_stride);
        __builtin_prefetch(src + 7 * src_stride);
        t0 = convolve8_8x8(s0, s1, s2, s3, s4, s5, s6, s7, filter_x);
        t1 = convolve8_8x8(s1, s2, s3, s4, s5, s6, s7, s8, filter_x);
        t2 = convolve8_8x8(s2, s3, s4, s5, s6, s7, s8, s9, filter_x);
        t3 = convolve8_8x8(s3, s4, s5, s6, s7, s8, s9, s10, filter_x);

        transpose_u8_8x4(&t0, &t1, &t2, &t3);
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t0), 0);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t1), 0);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t2), 0);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t3), 0);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t0), 1);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t1), 1);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t2), 1);
        dst += dst_stride;
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(t3), 1);
        dst += dst_stride;
        h -= 8;
      } while (h > 0);
    } else {
      uint8_t *d;
      int16x8_t s11, s12, s13, s14;

      do {
        __builtin_prefetch(src + 0 * src_stride);
        __builtin_prefetch(src + 1 * src_stride);
        __builtin_prefetch(src + 2 * src_stride);
        __builtin_prefetch(src + 3 * src_stride);
        __builtin_prefetch(src + 4 * src_stride);
        __builtin_prefetch(src + 5 * src_stride);
        __builtin_prefetch(src + 6 * src_stride);
        __builtin_prefetch(src + 7 * src_stride);
        load_u8_8x8(src, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
        transpose_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
        s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
        s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
        s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
        s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
        s4 = vreinterpretq_s16_u16(vmovl_u8(t4));
        s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
        s6 = vreinterpretq_s16_u16(vmovl_u8(t6));

        width = w;
        s = src + 7;
        d = dst;
        __builtin_prefetch(dst + 0 * dst_stride);
        __builtin_prefetch(dst + 1 * dst_stride);
        __builtin_prefetch(dst + 2 * dst_stride);
        __builtin_prefetch(dst + 3 * dst_stride);
        __builtin_prefetch(dst + 4 * dst_stride);
        __builtin_prefetch(dst + 5 * dst_stride);
        __builtin_prefetch(dst + 6 * dst_stride);
        __builtin_prefetch(dst + 7 * dst_stride);

        do {
          load_u8_8x8(s, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
          transpose_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
          s7 = vreinterpretq_s16_u16(vmovl_u8(t0));
          s8 = vreinterpretq_s16_u16(vmovl_u8(t1));
          s9 = vreinterpretq_s16_u16(vmovl_u8(t2));
          s10 = vreinterpretq_s16_u16(vmovl_u8(t3));
          s11 = vreinterpretq_s16_u16(vmovl_u8(t4));
          s12 = vreinterpretq_s16_u16(vmovl_u8(t5));
          s13 = vreinterpretq_s16_u16(vmovl_u8(t6));
          s14 = vreinterpretq_s16_u16(vmovl_u8(t7));

          t0 = convolve8_8x8(s0, s1, s2, s3, s4, s5, s6, s7, filter_x);
          t1 = convolve8_8x8(s1, s2, s3, s4, s5, s6, s7, s8, filter_x);
          t2 = convolve8_8x8(s2, s3, s4, s5, s6, s7, s8, s9, filter_x);
          t3 = convolve8_8x8(s3, s4, s5, s6, s7, s8, s9, s10, filter_x);
          t4 = convolve8_8x8(s4, s5, s6, s7, s8, s9, s10, s11, filter_x);
          t5 = convolve8_8x8(s5, s6, s7, s8, s9, s10, s11, s12, filter_x);
          t6 = convolve8_8x8(s6, s7, s8, s9, s10, s11, s12, s13, filter_x);
          t7 = convolve8_8x8(s7, s8, s9, s10, s11, s12, s13, s14, filter_x);

          transpose_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
          store_u8_8x8(d, dst_stride, t0, t1, t2, t3, t4, t5, t6, t7);

          s0 = s8;
          s1 = s9;
          s2 = s10;
          s3 = s11;
          s4 = s12;
          s5 = s13;
          s6 = s14;
          s += 8;
          d += 8;
          width -= 8;
        } while (width > 0);
        src += 8 * src_stride;
        dst += 8 * dst_stride;
        h -= 8;
      } while (h > 0);
    }
  }
}

void aom_convolve8_vert_neon(const uint8_t *src, ptrdiff_t src_stride,
                             uint8_t *dst, ptrdiff_t dst_stride,
                             const int16_t *filter_x, int x_step_q4,
                             const int16_t *filter_y, int y_step_q4, int w,
                             int h) {
  const InterpKernel *const filters_y = get_filter_base(filter_y);
  const int x0_q4 = get_filter_offset(filter_y, filters_y);
  const int16x8_t filters = vld1q_s16(filter_y);

  assert(!((intptr_t)dst & 3));
  assert(!(dst_stride & 3));
  assert(y_step_q4 == 16);

  (void)filter_x;
  (void)x_step_q4;

  src -= 3 * src_stride;

  if (w == 4) {
    uint8x8_t d01, d23;
    int16x4_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, d0, d1, d2, d3;

    s0 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s1 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s2 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s3 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s4 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s5 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;
    s6 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
    src += src_stride;

    do {
      s7 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
      src += src_stride;
      s8 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
      src += src_stride;
      s9 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
      src += src_stride;
      s10 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(vld1_u8(src))));
      src += src_stride;

      __builtin_prefetch(dst + 0 * dst_stride);
      __builtin_prefetch(dst + 1 * dst_stride);
      __builtin_prefetch(dst + 2 * dst_stride);
      __builtin_prefetch(dst + 3 * dst_stride);
      __builtin_prefetch(src + 0 * src_stride);
      __builtin_prefetch(src + 1 * src_stride);
      __builtin_prefetch(src + 2 * src_stride);
      __builtin_prefetch(src + 3 * src_stride);
      d0 = convolve8_4x4(s0, s1, s2, s3, s4, s5, s6, s7, filter_y);
      d1 = convolve8_4x4(s1, s2, s3, s4, s5, s6, s7, s8, filter_y);
      d2 = convolve8_4x4(s2, s3, s4, s5, s6, s7, s8, s9, filter_y);
      d3 = convolve8_4x4(s3, s4, s5, s6, s7, s8, s9, s10, filter_y);

      d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), 7);
      d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), 7);
      vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d01), 0);
      dst += dst_stride;
      vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d01), 1);
      dst += dst_stride;
      vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d23), 0);
      dst += dst_stride;
      vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(d23), 1);
      dst += dst_stride;

      s0 = s4;
      s1 = s5;
      s2 = s6;
      s3 = s7;
      s4 = s8;
      s5 = s9;
      s6 = s10;
      h -= 4;
    } while (h > 0);
  } else {
    int height;
    const uint8_t *s;
    uint8_t *d;
    uint8x8_t t0, t1, t2, t3;
    int16x8_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10;

    do {
      __builtin_prefetch(src + 0 * src_stride);
      __builtin_prefetch(src + 1 * src_stride);
      __builtin_prefetch(src + 2 * src_stride);
      __builtin_prefetch(src + 3 * src_stride);
      __builtin_prefetch(src + 4 * src_stride);
      __builtin_prefetch(src + 5 * src_stride);
      __builtin_prefetch(src + 6 * src_stride);
      s = src;
      s0 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s1 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s2 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s3 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s4 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s5 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      s6 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
      s += src_stride;
      d = dst;
      height = h;

      do {
        s7 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
        s += src_stride;
        s8 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
        s += src_stride;
        s9 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
        s += src_stride;
        s10 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(s)));
        s += src_stride;

        __builtin_prefetch(d + 0 * dst_stride);
        __builtin_prefetch(d + 1 * dst_stride);
        __builtin_prefetch(d + 2 * dst_stride);
        __builtin_prefetch(d + 3 * dst_stride);
        __builtin_prefetch(s + 0 * src_stride);
        __builtin_prefetch(s + 1 * src_stride);
        __builtin_prefetch(s + 2 * src_stride);
        __builtin_prefetch(s + 3 * src_stride);
        t0 = convolve8_8x8(s0, s1, s2, s3, s4, s5, s6, s7, filter_y);
        t1 = convolve8_8x8(s1, s2, s3, s4, s5, s6, s7, s8, filter_y);
        t2 = convolve8_8x8(s2, s3, s4, s5, s6, s7, s8, s9, filter_y);
        t3 = convolve8_8x8(s3, s4, s5, s6, s7, s8, s9, s10, filter_y);

        vst1_u8(d, t0);
        d += dst_stride;
        vst1_u8(d, t1);
        d += dst_stride;
        vst1_u8(d, t2);
        d += dst_stride;
        vst1_u8(d, t3);
        d += dst_stride;

        s0 = s4;
        s1 = s5;
        s2 = s6;
        s3 = s7;
        s4 = s8;
        s5 = s9;
        s6 = s10;
        height -= 4;
      } while (height > 0);
      src += 8;
      dst += 8;
      w -= 8;
    } while (w > 0);
  }
}

void aom_convolve8_neon(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                        ptrdiff_t dst_stride, const InterpKernel *filter,
                        int x0_q4, int x_step_q4, int y0_q4, int y_step_q4,
                        int w, int h) {
  /* Given our constraints: w <= 64, h <= 64, taps == 8 we can reduce the
   * maximum buffer size to 64 * 64 + 7 (+ 1 to make it divisible by 4).
   */
  uint8_t temp[64 * 72];

  // Account for the vertical phase needing 3 lines prior and 4 lines post
  // (+ 1 to make it divisible by 4).
  const int intermediate_height = h + 8;

  assert(y_step_q4 == 16);
  assert(x_step_q4 == 16);

  /* Filter starting 3 lines back. The neon implementation will ignore the given
   * height and filter a multiple of 4 lines. Since this goes in to the temp
   * buffer which has lots of extra room and is subsequently discarded this is
   * safe if somewhat less than ideal.   */
  aom_convolve8_horiz_neon(src - src_stride * 3, src_stride, temp, w,
                           filter[x0_q4], x_step_q4, NULL, y_step_q4, w,
                           intermediate_height);

  /* Step into the temp buffer 3 lines to get the actual frame data */
  aom_convolve8_vert_neon(temp + w * 3, w, dst, dst_stride, NULL, x_step_q4,
                          filter[y0_q4], y_step_q4, w, h);
}

void aom_convolve_copy_neon(const uint8_t *src, ptrdiff_t src_stride,
                            uint8_t *dst, ptrdiff_t dst_stride, int w, int h) {
  const uint8_t *src1;
  uint8_t *dst1;
  int y;

  if (!(w & 0x0F)) {
    for (y = 0; y < h; ++y) {
      src1 = src;
      dst1 = dst;
      for (int x = 0; x < (w >> 4); ++x) {
        vst1q_u8(dst1, vld1q_u8(src1));
        src1 += 16;
        dst1 += 16;
      }
      src += src_stride;
      dst += dst_stride;
    }
  } else if (!(w & 0x07)) {
    for (y = 0; y < h; ++y) {
      vst1_u8(dst, vld1_u8(src));
      src += src_stride;
      dst += dst_stride;
    }
  } else if (!(w & 0x03)) {
    for (y = 0; y < h; ++y) {
      vst1_lane_u32((uint32_t *)(dst), vreinterpret_u32_u8(vld1_u8(src)), 0);
      src += src_stride;
      dst += dst_stride;
    }
  } else if (!(w & 0x01)) {
    for (y = 0; y < h; ++y) {
      vst1_lane_u16((uint16_t *)(dst), vreinterpret_u16_u8(vld1_u8(src)), 0);
      src += src_stride;
      dst += dst_stride;
    }
  }
}
