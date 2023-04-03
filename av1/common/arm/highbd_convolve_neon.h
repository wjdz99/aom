/*
 *  Copyright (c) 2023, Alliance for Open Media. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AOM_AV1_COMMON_ARM_HIGHBD_CONVOLVE_NEON_H_
#define AOM_AV1_COMMON_ARM_HIGHBD_CONVOLVE_NEON_H_

#include <arm_neon.h>

static INLINE int32x4_t highbd_convolve8_4_s32(
    const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
    const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
    const int16x4_t s6, const int16x4_t s7, const int16x8_t y_filter,
    const int32x4_t offset_const) {
  const int16x4_t y_filter_lo = vget_low_s16(y_filter);
  const int16x4_t y_filter_hi = vget_high_s16(y_filter);

  int32x4_t sum = offset_const;
  sum = vmlal_lane_s16(sum, s0, y_filter_lo, 0);
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
    const int32x4_t offset_const) {
  int32x4_t sum = offset_const;
  sum = highbd_convolve8_4_s32(s0, s1, s2, s3, s4, s5, s6, s7, y_filter, sum);

  return vqrshrun_n_s32(sum, COMPOUND_ROUND1_BITS);
}

static INLINE void highbd_convolve8_8_s32(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
    const int16x8_t s6, const int16x8_t s7, const int16x8_t y_filter,
    int32x4_t *sum0, int32x4_t *sum1) {
  int16x4_t s0_lo = vget_low_s16(s0);
  int16x4_t s1_lo = vget_low_s16(s1);
  int16x4_t s2_lo = vget_low_s16(s2);
  int16x4_t s3_lo = vget_low_s16(s3);
  int16x4_t s4_lo = vget_low_s16(s4);
  int16x4_t s5_lo = vget_low_s16(s5);
  int16x4_t s6_lo = vget_low_s16(s6);
  int16x4_t s7_lo = vget_low_s16(s7);
  int16x4_t s0_hi = vget_high_s16(s0);
  int16x4_t s1_hi = vget_high_s16(s1);
  int16x4_t s2_hi = vget_high_s16(s2);
  int16x4_t s3_hi = vget_high_s16(s3);
  int16x4_t s4_hi = vget_high_s16(s4);
  int16x4_t s5_hi = vget_high_s16(s5);
  int16x4_t s6_hi = vget_high_s16(s6);
  int16x4_t s7_hi = vget_high_s16(s7);

  *sum0 = highbd_convolve8_4_s32(s0_lo, s1_lo, s2_lo, s3_lo, s4_lo, s5_lo,
                                 s6_lo, s7_lo, y_filter, *sum0);
  *sum1 = highbd_convolve8_4_s32(s0_hi, s1_hi, s2_hi, s3_hi, s4_hi, s5_hi,
                                 s6_hi, s7_hi, y_filter, *sum1);
}

static INLINE uint16x8_t highbd_convolve8_8_s32_s16(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
    const int16x8_t s6, const int16x8_t s7, const int16x8_t y_filter,
    const int32x4_t offset_const) {
  int32x4_t sum0 = offset_const;
  int32x4_t sum1 = offset_const;
  highbd_convolve8_8_s32(s0, s1, s2, s3, s4, s5, s6, s7, y_filter, &sum0,
                         &sum1);

  return vcombine_u16(vqrshrun_n_s32(sum0, COMPOUND_ROUND1_BITS),
                      vqrshrun_n_s32(sum1, COMPOUND_ROUND1_BITS));
}

#endif  // AOM_AV1_COMMON_ARM_HIGHBD_CONVOLVE_NEON_H_
