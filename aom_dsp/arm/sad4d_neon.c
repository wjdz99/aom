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

#include <arm_neon.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"

static INLINE unsigned int horizontal_long_add_16x8(const uint16x8_t vec_lo,
                                                    const uint16x8_t vec_hi) {
  const uint32x4_t vec_l_lo =
      vaddl_u16(vget_low_u16(vec_lo), vget_high_u16(vec_lo));
  const uint32x4_t vec_l_hi =
      vaddl_u16(vget_low_u16(vec_hi), vget_high_u16(vec_hi));
  const uint32x4_t a = vaddq_u32(vec_l_lo, vec_l_hi);
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

// Calculate the absolute difference of 64 bytes from vec_src_00, vec_src_16,
// vec_src_32, vec_src_48 and ref. Accumulate partial sums in vec_sum_ref_lo
// and vec_sum_ref_hi.
static void sad_neon_64(const uint8x16_t vec_src_00,
                        const uint8x16_t vec_src_16,
                        const uint8x16_t vec_src_32,
                        const uint8x16_t vec_src_48, const uint8_t *ref,
                        uint16x8_t *vec_sum_ref_lo,
                        uint16x8_t *vec_sum_ref_hi) {
  const uint8x16_t vec_ref_00 = vld1q_u8(ref);
  const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);
  const uint8x16_t vec_ref_32 = vld1q_u8(ref + 32);
  const uint8x16_t vec_ref_48 = vld1q_u8(ref + 48);

  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_00),
                             vget_low_u8(vec_ref_00));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_00),
                             vget_high_u8(vec_ref_00));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_16),
                             vget_low_u8(vec_ref_16));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_16),
                             vget_high_u8(vec_ref_16));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_32),
                             vget_low_u8(vec_ref_32));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_32),
                             vget_high_u8(vec_ref_32));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_48),
                             vget_low_u8(vec_ref_48));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_48),
                             vget_high_u8(vec_ref_48));
}

// Calculate the absolute difference of 32 bytes from vec_src_00, vec_src_16,
// and ref. Accumulate partial sums in vec_sum_ref_lo and vec_sum_ref_hi.
static void sad_neon_32(const uint8x16_t vec_src_00,
                        const uint8x16_t vec_src_16, const uint8_t *ref,
                        uint16x8_t *vec_sum_ref_lo,
                        uint16x8_t *vec_sum_ref_hi) {
  const uint8x16_t vec_ref_00 = vld1q_u8(ref);
  const uint8x16_t vec_ref_16 = vld1q_u8(ref + 16);

  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_00),
                             vget_low_u8(vec_ref_00));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_00),
                             vget_high_u8(vec_ref_00));
  *vec_sum_ref_lo = vabal_u8(*vec_sum_ref_lo, vget_low_u8(vec_src_16),
                             vget_low_u8(vec_ref_16));
  *vec_sum_ref_hi = vabal_u8(*vec_sum_ref_hi, vget_high_u8(vec_src_16),
                             vget_high_u8(vec_ref_16));
}

void aom_sad64x64x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 64; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);
    const uint8x16_t vec_src_32 = vld1q_u8(src + 32);
    const uint8x16_t vec_src_48 = vld1q_u8(src + 48);

    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref0,
                &vec_sum_ref0_lo, &vec_sum_ref0_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref1,
                &vec_sum_ref1_lo, &vec_sum_ref1_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref2,
                &vec_sum_ref2_lo, &vec_sum_ref2_hi);
    sad_neon_64(vec_src_00, vec_src_16, vec_src_32, vec_src_48, ref3,
                &vec_sum_ref3_lo, &vec_sum_ref3_hi);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}

void aom_sad32x32x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 32; ++i) {
    const uint8x16_t vec_src_00 = vld1q_u8(src);
    const uint8x16_t vec_src_16 = vld1q_u8(src + 16);

    sad_neon_32(vec_src_00, vec_src_16, ref0, &vec_sum_ref0_lo,
                &vec_sum_ref0_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref1, &vec_sum_ref1_lo,
                &vec_sum_ref1_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref2, &vec_sum_ref2_lo,
                &vec_sum_ref2_hi);
    sad_neon_32(vec_src_00, vec_src_16, ref3, &vec_sum_ref3_lo,
                &vec_sum_ref3_hi);

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}

void aom_sad16x16x4d_neon(const uint8_t *src, int src_stride,
                          const uint8_t *const ref[4], int ref_stride,
                          uint32_t *res) {
  int i;
  uint16x8_t vec_sum_ref0_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref0_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref1_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref2_hi = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_lo = vdupq_n_u16(0);
  uint16x8_t vec_sum_ref3_hi = vdupq_n_u16(0);
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];

  for (i = 0; i < 16; ++i) {
    const uint8x16_t vec_src = vld1q_u8(src);
    const uint8x16_t vec_ref0 = vld1q_u8(ref0);
    const uint8x16_t vec_ref1 = vld1q_u8(ref1);
    const uint8x16_t vec_ref2 = vld1q_u8(ref2);
    const uint8x16_t vec_ref3 = vld1q_u8(ref3);

    vec_sum_ref0_lo =
        vabal_u8(vec_sum_ref0_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref0));
    vec_sum_ref0_hi = vabal_u8(vec_sum_ref0_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref0));
    vec_sum_ref1_lo =
        vabal_u8(vec_sum_ref1_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref1));
    vec_sum_ref1_hi = vabal_u8(vec_sum_ref1_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref1));
    vec_sum_ref2_lo =
        vabal_u8(vec_sum_ref2_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref2));
    vec_sum_ref2_hi = vabal_u8(vec_sum_ref2_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref2));
    vec_sum_ref3_lo =
        vabal_u8(vec_sum_ref3_lo, vget_low_u8(vec_src), vget_low_u8(vec_ref3));
    vec_sum_ref3_hi = vabal_u8(vec_sum_ref3_hi, vget_high_u8(vec_src),
                               vget_high_u8(vec_ref3));

    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }

  res[0] = horizontal_long_add_16x8(vec_sum_ref0_lo, vec_sum_ref0_hi);
  res[1] = horizontal_long_add_16x8(vec_sum_ref1_lo, vec_sum_ref1_hi);
  res[2] = horizontal_long_add_16x8(vec_sum_ref2_lo, vec_sum_ref2_hi);
  res[3] = horizontal_long_add_16x8(vec_sum_ref3_lo, vec_sum_ref3_hi);
}

static INLINE unsigned int horizontal_add_16x8(const uint16x8_t vec_16x8) {
  const uint32x4_t a = vpaddlq_u16(vec_16x8);
  const uint64x2_t b = vpaddlq_u32(a);
  const uint32x2_t c = vadd_u32(vreinterpret_u32_u64(vget_low_u64(b)),
                                vreinterpret_u32_u64(vget_high_u64(b)));
  return vget_lane_u32(c, 0);
}

static void sad_row16_neon(uint16x8_t *vec_src, const uint8x16_t *q0,
                           const uint8_t *ref_ptr) {
  uint8x16_t q1 = vld1q_u8(ref_ptr);
  uint8x16_t q2 = vabdq_u8(*q0, q1);
  *vec_src = vpadalq_u8(*vec_src, q2);
}

void aom_sadMxNx4d_neon(int width, int height, const uint8_t *src,
                        int src_stride, const uint8_t *ref[4], int ref_stride,
                        uint32_t res[4]) {
  uint8x16_t q4;

  res[0] = 0;
  res[1] = 0;
  res[2] = 0;
  res[3] = 0;

  switch (width) {
    case 32: {
      for (int i = 0; i < height; i++) {
        uint16x8_t q0 = vdupq_n_u16(0);
        uint16x8_t q1 = vdupq_n_u16(0);
        uint16x8_t q2 = vdupq_n_u16(0);
        uint16x8_t q3 = vdupq_n_u16(0);

        q4 = vld1q_u8(src);

        sad_row16_neon(&q0, &q4, ref[0]);
        sad_row16_neon(&q1, &q4, ref[1]);
        sad_row16_neon(&q2, &q4, ref[2]);
        sad_row16_neon(&q3, &q4, ref[3]);

        q4 = vld1q_u8(src + 16);

        sad_row16_neon(&q0, &q4, ref[0] + 16);
        sad_row16_neon(&q1, &q4, ref[1] + 16);
        sad_row16_neon(&q2, &q4, ref[2] + 16);
        sad_row16_neon(&q3, &q4, ref[3] + 16);

        src += src_stride;
        ref[0] += ref_stride;
        ref[1] += ref_stride;
        ref[2] += ref_stride;
        ref[3] += ref_stride;

        res[0] += horizontal_add_16x8(q0);
        res[1] += horizontal_add_16x8(q1);
        res[2] += horizontal_add_16x8(q2);
        res[3] += horizontal_add_16x8(q3);
      }
      break;
    }
    case 64: {
      for (int i = 0; i < height; i++) {
        uint16x8_t q0 = vdupq_n_u16(0);
        uint16x8_t q1 = vdupq_n_u16(0);
        uint16x8_t q2 = vdupq_n_u16(0);
        uint16x8_t q3 = vdupq_n_u16(0);

        q4 = vld1q_u8(src);

        sad_row16_neon(&q0, &q4, ref[0]);
        sad_row16_neon(&q1, &q4, ref[1]);
        sad_row16_neon(&q2, &q4, ref[2]);
        sad_row16_neon(&q3, &q4, ref[3]);

        q4 = vld1q_u8(src + 16);

        sad_row16_neon(&q0, &q4, ref[0] + 16);
        sad_row16_neon(&q1, &q4, ref[1] + 16);
        sad_row16_neon(&q2, &q4, ref[2] + 16);
        sad_row16_neon(&q3, &q4, ref[3] + 16);

        q4 = vld1q_u8(src + 32);

        sad_row16_neon(&q0, &q4, ref[0] + 32);
        sad_row16_neon(&q1, &q4, ref[1] + 32);
        sad_row16_neon(&q2, &q4, ref[2] + 32);
        sad_row16_neon(&q3, &q4, ref[3] + 32);

        q4 = vld1q_u8(src + 48);

        sad_row16_neon(&q0, &q4, ref[0] + 48);
        sad_row16_neon(&q1, &q4, ref[1] + 48);
        sad_row16_neon(&q2, &q4, ref[2] + 48);
        sad_row16_neon(&q3, &q4, ref[3] + 48);

        src += src_stride;
        ref[0] += ref_stride;
        ref[1] += ref_stride;
        ref[2] += ref_stride;
        ref[3] += ref_stride;

        res[0] += horizontal_add_16x8(q0);
        res[1] += horizontal_add_16x8(q1);
        res[2] += horizontal_add_16x8(q2);
        res[3] += horizontal_add_16x8(q3);
      }
      break;
    }
    case 128: {
      for (int i = 0; i < height; i++) {
        uint16x8_t q0 = vdupq_n_u16(0);
        uint16x8_t q1 = vdupq_n_u16(0);
        uint16x8_t q2 = vdupq_n_u16(0);
        uint16x8_t q3 = vdupq_n_u16(0);

        q4 = vld1q_u8(src);

        sad_row16_neon(&q0, &q4, ref[0]);
        sad_row16_neon(&q1, &q4, ref[1]);
        sad_row16_neon(&q2, &q4, ref[2]);
        sad_row16_neon(&q3, &q4, ref[3]);

        q4 = vld1q_u8(src + 16);

        sad_row16_neon(&q0, &q4, ref[0] + 16);
        sad_row16_neon(&q1, &q4, ref[1] + 16);
        sad_row16_neon(&q2, &q4, ref[2] + 16);
        sad_row16_neon(&q3, &q4, ref[3] + 16);

        q4 = vld1q_u8(src + 32);

        sad_row16_neon(&q0, &q4, ref[0] + 32);
        sad_row16_neon(&q1, &q4, ref[1] + 32);
        sad_row16_neon(&q2, &q4, ref[2] + 32);
        sad_row16_neon(&q3, &q4, ref[3] + 32);

        q4 = vld1q_u8(src + 48);

        sad_row16_neon(&q0, &q4, ref[0] + 48);
        sad_row16_neon(&q1, &q4, ref[1] + 48);
        sad_row16_neon(&q2, &q4, ref[2] + 48);
        sad_row16_neon(&q3, &q4, ref[3] + 48);

        q4 = vld1q_u8(src + 64);

        sad_row16_neon(&q0, &q4, ref[0] + 64);
        sad_row16_neon(&q1, &q4, ref[1] + 64);
        sad_row16_neon(&q2, &q4, ref[2] + 64);
        sad_row16_neon(&q3, &q4, ref[3] + 64);

        q4 = vld1q_u8(src + 80);

        sad_row16_neon(&q0, &q4, ref[0] + 80);
        sad_row16_neon(&q1, &q4, ref[1] + 80);
        sad_row16_neon(&q2, &q4, ref[2] + 80);
        sad_row16_neon(&q3, &q4, ref[3] + 80);

        q4 = vld1q_u8(src + 96);

        sad_row16_neon(&q0, &q4, ref[0] + 96);
        sad_row16_neon(&q1, &q4, ref[1] + 96);
        sad_row16_neon(&q2, &q4, ref[2] + 96);
        sad_row16_neon(&q3, &q4, ref[3] + 96);

        q4 = vld1q_u8(src + 112);

        sad_row16_neon(&q0, &q4, ref[0] + 112);
        sad_row16_neon(&q1, &q4, ref[1] + 112);
        sad_row16_neon(&q2, &q4, ref[2] + 112);
        sad_row16_neon(&q3, &q4, ref[3] + 112);

        src += src_stride;
        ref[0] += ref_stride;
        ref[1] += ref_stride;
        ref[2] += ref_stride;
        ref[3] += ref_stride;

        res[0] += horizontal_add_16x8(q0);
        res[1] += horizontal_add_16x8(q1);
        res[2] += horizontal_add_16x8(q2);
        res[3] += horizontal_add_16x8(q3);
      }
    }
  }
}

#define sad_skip_MxN_neon(m, n)                                                \
  void aom_sad_skip_##m##x##n##x4d_neon(const uint8_t *src, int src_stride,    \
                                        const uint8_t *ref[4], int ref_stride, \
                                        uint32_t res[4]) {                     \
    aom_sadMxNx4d_neon(m, ((n) >> 1), src, 2 * src_stride, ref,                \
                       2 * ref_stride, res);                                   \
    res[0] <<= 1;                                                              \
    res[1] <<= 1;                                                              \
    res[2] <<= 1;                                                              \
    res[3] <<= 1;                                                              \
  }

sad_skip_MxN_neon(32, 8);
sad_skip_MxN_neon(32, 16);
sad_skip_MxN_neon(32, 32);
sad_skip_MxN_neon(32, 64);

sad_skip_MxN_neon(64, 16);
sad_skip_MxN_neon(64, 32);
sad_skip_MxN_neon(64, 64);
sad_skip_MxN_neon(64, 128);

sad_skip_MxN_neon(128, 64);
sad_skip_MxN_neon(128, 128);
#undef sad_skip_MxN_neon
