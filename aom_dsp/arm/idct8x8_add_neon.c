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

#include "./aom_config.h"
#include "./av1_rtcd.h"
#include "aom_dsp/txfm_common.h"
#include "av1/common/arm/neon/common_neon.h"

void aom_idct8x8_64_add_neon(int16_t *input, uint8_t *dest, int dest_stride) {
  av1_iht8x8_64_add_neon(input, dest, dest_stride, DCT_DCT);
}

void aom_idct8x8_12_add_neon(int16_t *input, uint8_t *dest, int dest_stride) {
  uint8_t *d1, *d2;
  uint8x8_t d0u8, d1u8, d2u8, d3u8;
  int16x4_t d10s16, d11s16, d12s16, d13s16, d16s16;
  int16x4_t d26s16, d27s16, d28s16, d29s16;
  uint64x1_t d0u64, d1u64, d2u64, d3u64;
  int16x8_t q0s16, q1s16, q2s16, q3s16, q4s16, q5s16, q6s16, q7s16;
  int16x8_t q8s16, q9s16, q10s16, q11s16, q12s16, q13s16, q14s16, q15s16;
  uint16x8_t q8u16, q9u16, q10u16, q11u16;
  int32x4_t q9s32, q10s32, q11s32, q12s32;

  q8s16 = vld1q_s16(input);
  q9s16 = vld1q_s16(input + 8);
  q10s16 = vld1q_s16(input + 16);
  q11s16 = vld1q_s16(input + 24);
  q12s16 = vld1q_s16(input + 32);
  q13s16 = vld1q_s16(input + 40);
  q14s16 = vld1q_s16(input + 48);
  q15s16 = vld1q_s16(input + 56);

  av1_transpose_8x8_neon(&q8s16, &q9s16, &q10s16, &q11s16, &q12s16, &q13s16,
                         &q14s16, &q15s16);

  // First transform rows
  // stage 1
  q0s16 = vdupq_n_s16((int16_t)cospi_28_64 * 2);
  q1s16 = vdupq_n_s16((int16_t)cospi_4_64 * 2);

  q4s16 = vqrdmulhq_s16(q9s16, q0s16);

  q0s16 = vdupq_n_s16(-(int16_t)cospi_20_64 * 2);

  q7s16 = vqrdmulhq_s16(q9s16, q1s16);

  q1s16 = vdupq_n_s16((int16_t)cospi_12_64 * 2);

  q5s16 = vqrdmulhq_s16(q11s16, q0s16);

  q0s16 = vdupq_n_s16((int16_t)cospi_16_64 * 2);

  q6s16 = vqrdmulhq_s16(q11s16, q1s16);

  // stage 2 & stage 3 - even half
  q1s16 = vdupq_n_s16((int16_t)cospi_24_64 * 2);

  q9s16 = vqrdmulhq_s16(q8s16, q0s16);

  q0s16 = vdupq_n_s16((int16_t)cospi_8_64 * 2);

  q13s16 = vqrdmulhq_s16(q10s16, q1s16);

  q15s16 = vqrdmulhq_s16(q10s16, q0s16);

  // stage 3 -odd half
  q0s16 = vaddq_s16(q9s16, q15s16);
  q1s16 = vaddq_s16(q9s16, q13s16);
  q2s16 = vsubq_s16(q9s16, q13s16);
  q3s16 = vsubq_s16(q9s16, q15s16);

  // stage 2 - odd half
  q13s16 = vsubq_s16(q4s16, q5s16);
  q4s16 = vaddq_s16(q4s16, q5s16);
  q14s16 = vsubq_s16(q7s16, q6s16);
  q7s16 = vaddq_s16(q7s16, q6s16);
  d26s16 = vget_low_s16(q13s16);
  d27s16 = vget_high_s16(q13s16);
  d28s16 = vget_low_s16(q14s16);
  d29s16 = vget_high_s16(q14s16);

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

  // stage 4
  q8s16 = vaddq_s16(q0s16, q7s16);
  q9s16 = vaddq_s16(q1s16, q6s16);
  q10s16 = vaddq_s16(q2s16, q5s16);
  q11s16 = vaddq_s16(q3s16, q4s16);
  q12s16 = vsubq_s16(q3s16, q4s16);
  q13s16 = vsubq_s16(q2s16, q5s16);
  q14s16 = vsubq_s16(q1s16, q6s16);
  q15s16 = vsubq_s16(q0s16, q7s16);

  av1_transpose_8x8_neon(&q8s16, &q9s16, &q10s16, &q11s16, &q12s16, &q13s16,
                         &q14s16, &q15s16);

  av1_idct_8x8_1d_neon(&q8s16, &q9s16, &q10s16, &q11s16, &q12s16, &q13s16,
                       &q14s16, &q15s16);

  q8s16 = vrshrq_n_s16(q8s16, 5);
  q9s16 = vrshrq_n_s16(q9s16, 5);
  q10s16 = vrshrq_n_s16(q10s16, 5);
  q11s16 = vrshrq_n_s16(q11s16, 5);
  q12s16 = vrshrq_n_s16(q12s16, 5);
  q13s16 = vrshrq_n_s16(q13s16, 5);
  q14s16 = vrshrq_n_s16(q14s16, 5);
  q15s16 = vrshrq_n_s16(q15s16, 5);

  d1 = d2 = dest;

  d0u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d1u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d2u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d3u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;

  q8u16 = vaddw_u8(vreinterpretq_u16_s16(q8s16), vreinterpret_u8_u64(d0u64));
  q9u16 = vaddw_u8(vreinterpretq_u16_s16(q9s16), vreinterpret_u8_u64(d1u64));
  q10u16 = vaddw_u8(vreinterpretq_u16_s16(q10s16), vreinterpret_u8_u64(d2u64));
  q11u16 = vaddw_u8(vreinterpretq_u16_s16(q11s16), vreinterpret_u8_u64(d3u64));

  d0u8 = vqmovun_s16(vreinterpretq_s16_u16(q8u16));
  d1u8 = vqmovun_s16(vreinterpretq_s16_u16(q9u16));
  d2u8 = vqmovun_s16(vreinterpretq_s16_u16(q10u16));
  d3u8 = vqmovun_s16(vreinterpretq_s16_u16(q11u16));

  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d0u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d1u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d2u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d3u8));
  d2 += dest_stride;

  q8s16 = q12s16;
  q9s16 = q13s16;
  q10s16 = q14s16;
  q11s16 = q15s16;

  d0u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d1u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d2u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;
  d3u64 = vld1_u64((uint64_t *)d1);
  d1 += dest_stride;

  q8u16 = vaddw_u8(vreinterpretq_u16_s16(q8s16), vreinterpret_u8_u64(d0u64));
  q9u16 = vaddw_u8(vreinterpretq_u16_s16(q9s16), vreinterpret_u8_u64(d1u64));
  q10u16 = vaddw_u8(vreinterpretq_u16_s16(q10s16), vreinterpret_u8_u64(d2u64));
  q11u16 = vaddw_u8(vreinterpretq_u16_s16(q11s16), vreinterpret_u8_u64(d3u64));

  d0u8 = vqmovun_s16(vreinterpretq_s16_u16(q8u16));
  d1u8 = vqmovun_s16(vreinterpretq_s16_u16(q9u16));
  d2u8 = vqmovun_s16(vreinterpretq_s16_u16(q10u16));
  d3u8 = vqmovun_s16(vreinterpretq_s16_u16(q11u16));

  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d0u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d1u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d2u8));
  d2 += dest_stride;
  vst1_u64((uint64_t *)d2, vreinterpret_u64_u8(d3u8));
  d2 += dest_stride;
  return;
}
