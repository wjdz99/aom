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

#ifndef AV1_COMMON_ARM_NEON_COMMON_NEON_H_
#define AV1_COMMON_ARM_NEON_COMMON_NEON_H_

#include <arm_neon.h>

#include "./aom_config.h"

static INLINE void av1_transpose_8x8_neon(int16x8_t *q8s16, int16x8_t *q9s16,
                                          int16x8_t *q10s16, int16x8_t *q11s16,
                                          int16x8_t *q12s16, int16x8_t *q13s16,
                                          int16x8_t *q14s16,
                                          int16x8_t *q15s16) {
  int16x4_t d16s16, d17s16, d18s16, d19s16, d20s16, d21s16, d22s16, d23s16;
  int16x4_t d24s16, d25s16, d26s16, d27s16, d28s16, d29s16, d30s16, d31s16;
  int32x4x2_t q0x2s32, q1x2s32, q2x2s32, q3x2s32;
  int16x8x2_t q0x2s16, q1x2s16, q2x2s16, q3x2s16;

  d16s16 = vget_low_s16(*q8s16);
  d17s16 = vget_high_s16(*q8s16);
  d18s16 = vget_low_s16(*q9s16);
  d19s16 = vget_high_s16(*q9s16);
  d20s16 = vget_low_s16(*q10s16);
  d21s16 = vget_high_s16(*q10s16);
  d22s16 = vget_low_s16(*q11s16);
  d23s16 = vget_high_s16(*q11s16);
  d24s16 = vget_low_s16(*q12s16);
  d25s16 = vget_high_s16(*q12s16);
  d26s16 = vget_low_s16(*q13s16);
  d27s16 = vget_high_s16(*q13s16);
  d28s16 = vget_low_s16(*q14s16);
  d29s16 = vget_high_s16(*q14s16);
  d30s16 = vget_low_s16(*q15s16);
  d31s16 = vget_high_s16(*q15s16);

  *q8s16 = vcombine_s16(d16s16, d24s16);   // vswp d17, d24
  *q9s16 = vcombine_s16(d18s16, d26s16);   // vswp d19, d26
  *q10s16 = vcombine_s16(d20s16, d28s16);  // vswp d21, d28
  *q11s16 = vcombine_s16(d22s16, d30s16);  // vswp d23, d30
  *q12s16 = vcombine_s16(d17s16, d25s16);
  *q13s16 = vcombine_s16(d19s16, d27s16);
  *q14s16 = vcombine_s16(d21s16, d29s16);
  *q15s16 = vcombine_s16(d23s16, d31s16);

  q0x2s32 =
      vtrnq_s32(vreinterpretq_s32_s16(*q8s16), vreinterpretq_s32_s16(*q10s16));
  q1x2s32 =
      vtrnq_s32(vreinterpretq_s32_s16(*q9s16), vreinterpretq_s32_s16(*q11s16));
  q2x2s32 =
      vtrnq_s32(vreinterpretq_s32_s16(*q12s16), vreinterpretq_s32_s16(*q14s16));
  q3x2s32 =
      vtrnq_s32(vreinterpretq_s32_s16(*q13s16), vreinterpretq_s32_s16(*q15s16));

  q0x2s16 = vtrnq_s16(vreinterpretq_s16_s32(q0x2s32.val[0]),   // q8
                      vreinterpretq_s16_s32(q1x2s32.val[0]));  // q9
  q1x2s16 = vtrnq_s16(vreinterpretq_s16_s32(q0x2s32.val[1]),   // q10
                      vreinterpretq_s16_s32(q1x2s32.val[1]));  // q11
  q2x2s16 = vtrnq_s16(vreinterpretq_s16_s32(q2x2s32.val[0]),   // q12
                      vreinterpretq_s16_s32(q3x2s32.val[0]));  // q13
  q3x2s16 = vtrnq_s16(vreinterpretq_s16_s32(q2x2s32.val[1]),   // q14
                      vreinterpretq_s16_s32(q3x2s32.val[1]));  // q15

  *q8s16 = q0x2s16.val[0];
  *q9s16 = q0x2s16.val[1];
  *q10s16 = q1x2s16.val[0];
  *q11s16 = q1x2s16.val[1];
  *q12s16 = q2x2s16.val[0];
  *q13s16 = q2x2s16.val[1];
  *q14s16 = q3x2s16.val[0];
  *q15s16 = q3x2s16.val[1];
  return;
}

static INLINE void av1_idct_8x8_1d_neon(int16x8_t *q8s16, int16x8_t *q9s16,
                                        int16x8_t *q10s16, int16x8_t *q11s16,
                                        int16x8_t *q12s16, int16x8_t *q13s16,
                                        int16x8_t *q14s16, int16x8_t *q15s16) {
  int16x4_t d0s16, d1s16, d2s16, d3s16;
  int16x4_t d8s16, d9s16, d10s16, d11s16, d12s16, d13s16, d14s16, d15s16;
  int16x4_t d16s16, d17s16, d18s16, d19s16, d20s16, d21s16, d22s16, d23s16;
  int16x4_t d24s16, d25s16, d26s16, d27s16, d28s16, d29s16, d30s16, d31s16;
  int16x8_t q0s16, q1s16, q2s16, q3s16, q4s16, q5s16, q6s16, q7s16;
  int32x4_t q2s32, q3s32, q5s32, q6s32, q8s32, q9s32;
  int32x4_t q10s32, q11s32, q12s32, q13s32, q15s32;

  d0s16 = vdup_n_s16((int16_t)cospi_28_64);
  d1s16 = vdup_n_s16((int16_t)cospi_4_64);
  d2s16 = vdup_n_s16((int16_t)cospi_12_64);
  d3s16 = vdup_n_s16((int16_t)cospi_20_64);

  d16s16 = vget_low_s16(*q8s16);
  d17s16 = vget_high_s16(*q8s16);
  d18s16 = vget_low_s16(*q9s16);
  d19s16 = vget_high_s16(*q9s16);
  d20s16 = vget_low_s16(*q10s16);
  d21s16 = vget_high_s16(*q10s16);
  d22s16 = vget_low_s16(*q11s16);
  d23s16 = vget_high_s16(*q11s16);
  d24s16 = vget_low_s16(*q12s16);
  d25s16 = vget_high_s16(*q12s16);
  d26s16 = vget_low_s16(*q13s16);
  d27s16 = vget_high_s16(*q13s16);
  d28s16 = vget_low_s16(*q14s16);
  d29s16 = vget_high_s16(*q14s16);
  d30s16 = vget_low_s16(*q15s16);
  d31s16 = vget_high_s16(*q15s16);

  q2s32 = vmull_s16(d18s16, d0s16);
  q3s32 = vmull_s16(d19s16, d0s16);
  q5s32 = vmull_s16(d26s16, d2s16);
  q6s32 = vmull_s16(d27s16, d2s16);

  q2s32 = vmlsl_s16(q2s32, d30s16, d1s16);
  q3s32 = vmlsl_s16(q3s32, d31s16, d1s16);
  q5s32 = vmlsl_s16(q5s32, d22s16, d3s16);
  q6s32 = vmlsl_s16(q6s32, d23s16, d3s16);

  d8s16 = vqrshrn_n_s32(q2s32, 14);
  d9s16 = vqrshrn_n_s32(q3s32, 14);
  d10s16 = vqrshrn_n_s32(q5s32, 14);
  d11s16 = vqrshrn_n_s32(q6s32, 14);
  q4s16 = vcombine_s16(d8s16, d9s16);
  q5s16 = vcombine_s16(d10s16, d11s16);

  q2s32 = vmull_s16(d18s16, d1s16);
  q3s32 = vmull_s16(d19s16, d1s16);
  q9s32 = vmull_s16(d26s16, d3s16);
  q13s32 = vmull_s16(d27s16, d3s16);

  q2s32 = vmlal_s16(q2s32, d30s16, d0s16);
  q3s32 = vmlal_s16(q3s32, d31s16, d0s16);
  q9s32 = vmlal_s16(q9s32, d22s16, d2s16);
  q13s32 = vmlal_s16(q13s32, d23s16, d2s16);

  d14s16 = vqrshrn_n_s32(q2s32, 14);
  d15s16 = vqrshrn_n_s32(q3s32, 14);
  d12s16 = vqrshrn_n_s32(q9s32, 14);
  d13s16 = vqrshrn_n_s32(q13s32, 14);
  q6s16 = vcombine_s16(d12s16, d13s16);
  q7s16 = vcombine_s16(d14s16, d15s16);

  d0s16 = vdup_n_s16((int16_t)cospi_16_64);

  q2s32 = vmull_s16(d16s16, d0s16);
  q3s32 = vmull_s16(d17s16, d0s16);
  q13s32 = vmull_s16(d16s16, d0s16);
  q15s32 = vmull_s16(d17s16, d0s16);

  q2s32 = vmlal_s16(q2s32, d24s16, d0s16);
  q3s32 = vmlal_s16(q3s32, d25s16, d0s16);
  q13s32 = vmlsl_s16(q13s32, d24s16, d0s16);
  q15s32 = vmlsl_s16(q15s32, d25s16, d0s16);

  d0s16 = vdup_n_s16((int16_t)cospi_24_64);
  d1s16 = vdup_n_s16((int16_t)cospi_8_64);

  d18s16 = vqrshrn_n_s32(q2s32, 14);
  d19s16 = vqrshrn_n_s32(q3s32, 14);
  d22s16 = vqrshrn_n_s32(q13s32, 14);
  d23s16 = vqrshrn_n_s32(q15s32, 14);
  *q9s16 = vcombine_s16(d18s16, d19s16);
  *q11s16 = vcombine_s16(d22s16, d23s16);

  q2s32 = vmull_s16(d20s16, d0s16);
  q3s32 = vmull_s16(d21s16, d0s16);
  q8s32 = vmull_s16(d20s16, d1s16);
  q12s32 = vmull_s16(d21s16, d1s16);

  q2s32 = vmlsl_s16(q2s32, d28s16, d1s16);
  q3s32 = vmlsl_s16(q3s32, d29s16, d1s16);
  q8s32 = vmlal_s16(q8s32, d28s16, d0s16);
  q12s32 = vmlal_s16(q12s32, d29s16, d0s16);

  d26s16 = vqrshrn_n_s32(q2s32, 14);
  d27s16 = vqrshrn_n_s32(q3s32, 14);
  d30s16 = vqrshrn_n_s32(q8s32, 14);
  d31s16 = vqrshrn_n_s32(q12s32, 14);
  *q13s16 = vcombine_s16(d26s16, d27s16);
  *q15s16 = vcombine_s16(d30s16, d31s16);

  q0s16 = vaddq_s16(*q9s16, *q15s16);
  q1s16 = vaddq_s16(*q11s16, *q13s16);
  q2s16 = vsubq_s16(*q11s16, *q13s16);
  q3s16 = vsubq_s16(*q9s16, *q15s16);

  *q13s16 = vsubq_s16(q4s16, q5s16);
  q4s16 = vaddq_s16(q4s16, q5s16);
  *q14s16 = vsubq_s16(q7s16, q6s16);
  q7s16 = vaddq_s16(q7s16, q6s16);
  d26s16 = vget_low_s16(*q13s16);
  d27s16 = vget_high_s16(*q13s16);
  d28s16 = vget_low_s16(*q14s16);
  d29s16 = vget_high_s16(*q14s16);

  d16s16 = vdup_n_s16((int16_t)cospi_16_64);

  q9s32 = vmull_s16(d28s16, d16s16);
  q10s32 = vmull_s16(d29s16, d16s16);
  q11s32 = vmull_s16(d28s16, d16s16);
  q12s32 = vmull_s16(d29s16, d16s16);

  q9s32 = vmlsl_s16(q9s32, d26s16, d16s16);
  q10s32 = vmlsl_s16(q10s32, d27s16, d16s16);
  q11s32 = vmlal_s16(q11s32, d26s16, d16s16);
  q12s32 = vmlal_s16(q12s32, d27s16, d16s16);

  d10s16 = vqrshrn_n_s32(q9s32, 14);
  d11s16 = vqrshrn_n_s32(q10s32, 14);
  d12s16 = vqrshrn_n_s32(q11s32, 14);
  d13s16 = vqrshrn_n_s32(q12s32, 14);
  q5s16 = vcombine_s16(d10s16, d11s16);
  q6s16 = vcombine_s16(d12s16, d13s16);

  *q8s16 = vaddq_s16(q0s16, q7s16);
  *q9s16 = vaddq_s16(q1s16, q6s16);
  *q10s16 = vaddq_s16(q2s16, q5s16);
  *q11s16 = vaddq_s16(q3s16, q4s16);
  *q12s16 = vsubq_s16(q3s16, q4s16);
  *q13s16 = vsubq_s16(q2s16, q5s16);
  *q14s16 = vsubq_s16(q1s16, q6s16);
  *q15s16 = vsubq_s16(q0s16, q7s16);
  return;
}

#endif  // AV1_COMMON_ARM_NEON_COMMON_NEON_H_
