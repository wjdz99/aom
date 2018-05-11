/*
 *  Copyright (c) 2018, Alliance for Open Media. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_CONVOLVE_NEON_H_
#define AV1_CONVOLVE_NEON_H_

#include <arm_neon.h>

static INLINE void load_s16_8x8(const uint8_t *s, const ptrdiff_t p,
                                uint8x8_t *const s0, uint8x8_t *const s1,
                                uint8x8_t *const s2, uint8x8_t *const s3) {
  *s0 = vld1_u8(s);
  s += p;
  *s1 = vld1_u8(s);
  s += p;
  *s2 = vld1_u8(s);
  s += p;
  *s3 = vld1_u8(s);
}

static INLINE void transpose_u16_4x8(uint16x4_t *a0, uint16x4_t *a1,
                                     uint16x4_t *a2, uint16x4_t *a3,
                                     uint16x4_t *a4, uint16x4_t *a5,
                                     uint16x4_t *a6, uint16x4_t *a7,
                                     uint16x8_t *o0, uint16x8_t *o1,
                                     uint16x8_t *o2, uint16x8_t *o3) {
  // Swap 16 bit elements. Goes from:
  // a0: 00 01 02 03
  // a1: 10 11 12 13
  // a2: 20 21 22 23
  // a3: 30 31 32 33
  // a4: 40 41 42 43
  // a5: 50 51 52 53
  // a6: 60 61 62 63
  // a7: 70 71 72 73
  // to:
  // b0.val[0]: 00 10 02 12
  // b0.val[1]: 01 11 03 13
  // b1.val[0]: 20 30 22 32
  // b1.val[1]: 21 31 23 33
  // b2.val[0]: 40 50 42 52
  // b2.val[1]: 41 51 43 53
  // b3.val[0]: 60 70 62 72
  // b3.val[1]: 61 71 63 73

  uint16x4x2_t b0 = vtrn_u16(*a0, *a1);
  uint16x4x2_t b1 = vtrn_u16(*a2, *a3);
  uint16x4x2_t b2 = vtrn_u16(*a4, *a5);
  uint16x4x2_t b3 = vtrn_u16(*a6, *a7);

  // Swap 32 bit elements resulting in:
  // c0.val[0]: 00 10 20 30
  // c0.val[1]: 02 12 22 32
  // c1.val[0]: 01 11 21 31
  // c1.val[1]: 03 13 23 33
  // c2.val[0]: 40 50 60 70
  // c2.val[1]: 42 52 62 72
  // c3.val[0]: 41 51 61 71
  // c3.val[1]: 43 53 63 73

  uint32x2x2_t c0 = vtrn_u32(vreinterpret_u32_u16(b0.val[0]),
                             vreinterpret_u32_u16(b1.val[0]));
  uint32x2x2_t c1 = vtrn_u32(vreinterpret_u32_u16(b0.val[1]),
                             vreinterpret_u32_u16(b1.val[1]));
  uint32x2x2_t c2 = vtrn_u32(vreinterpret_u32_u16(b2.val[0]),
                             vreinterpret_u32_u16(b3.val[0]));
  uint32x2x2_t c3 = vtrn_u32(vreinterpret_u32_u16(b2.val[1]),
                             vreinterpret_u32_u16(b3.val[1]));

  // Swap 64 bit elements resulting in:
  // o0: 00 10 20 30 40 50 60 70
  // o1: 01 11 21 31 41 51 61 71
  // o2: 02 12 22 32 42 52 62 72
  // o3: 03 13 23 33 43 53 63 73

  *o0 = vcombine_u16(vreinterpret_u16_u32(c0.val[0]),
                     vreinterpret_u16_u32(c2.val[0]));
  *o1 = vcombine_u16(vreinterpret_u16_u32(c1.val[0]),
                     vreinterpret_u16_u32(c3.val[0]));
  *o2 = vcombine_u16(vreinterpret_u16_u32(c0.val[1]),
                     vreinterpret_u16_u32(c2.val[1]));
  *o3 = vcombine_u16(vreinterpret_u16_u32(c1.val[1]),
                     vreinterpret_u16_u32(c3.val[1]));
}

static INLINE void transpose_u16_8x8(uint16x8_t *a0, uint16x8_t *a1,
                                     uint16x8_t *a2, uint16x8_t *a3,
                                     uint16x8_t *a4, uint16x8_t *a5,
                                     uint16x8_t *a6, uint16x8_t *a7) {
  // Swap 16 bit elements. Goes from:
  // a0: 00 01 02 03 04 05 06 07
  // a1: 10 11 12 13 14 15 16 17
  // a2: 20 21 22 23 24 25 26 27
  // a3: 30 31 32 33 34 35 36 37
  // a4: 40 41 42 43 44 45 46 47
  // a5: 50 51 52 53 54 55 56 57
  // a6: 60 61 62 63 64 65 66 67
  // a7: 70 71 72 73 74 75 76 77
  // to:
  // b0.val[0]: 00 10 02 12 04 14 06 16
  // b0.val[1]: 01 11 03 13 05 15 07 17
  // b1.val[0]: 20 30 22 32 24 34 26 36
  // b1.val[1]: 21 31 23 33 25 35 27 37
  // b2.val[0]: 40 50 42 52 44 54 46 56
  // b2.val[1]: 41 51 43 53 45 55 47 57
  // b3.val[0]: 60 70 62 72 64 74 66 76
  // b3.val[1]: 61 71 63 73 65 75 67 77

  const uint16x8x2_t b0 = vtrnq_u16(*a0, *a1);
  const uint16x8x2_t b1 = vtrnq_u16(*a2, *a3);
  const uint16x8x2_t b2 = vtrnq_u16(*a4, *a5);
  const uint16x8x2_t b3 = vtrnq_u16(*a6, *a7);

  // Swap 32 bit elements resulting in:
  // c0.val[0]: 00 10 20 30 04 14 24 34
  // c0.val[1]: 02 12 22 32 06 16 26 36
  // c1.val[0]: 01 11 21 31 05 15 25 35
  // c1.val[1]: 03 13 23 33 07 17 27 37
  // c2.val[0]: 40 50 60 70 44 54 64 74
  // c2.val[1]: 42 52 62 72 46 56 66 76
  // c3.val[0]: 41 51 61 71 45 55 65 75
  // c3.val[1]: 43 53 63 73 47 57 67 77

  const uint32x4x2_t c0 = vtrnq_u32(vreinterpretq_u32_u16(b0.val[0]),
                                    vreinterpretq_u32_u16(b1.val[0]));
  const uint32x4x2_t c1 = vtrnq_u32(vreinterpretq_u32_u16(b0.val[1]),
                                    vreinterpretq_u32_u16(b1.val[1]));
  const uint32x4x2_t c2 = vtrnq_u32(vreinterpretq_u32_u16(b2.val[0]),
                                    vreinterpretq_u32_u16(b3.val[0]));
  const uint32x4x2_t c3 = vtrnq_u32(vreinterpretq_u32_u16(b2.val[1]),
                                    vreinterpretq_u32_u16(b3.val[1]));

  *a0 = vcombine_u16(vget_low_u16(vreinterpretq_u16_u32(c0.val[0])),
                     vget_low_u16(vreinterpretq_u16_u32(c2.val[0])));
  *a4 = vcombine_u16(vget_high_u16(vreinterpretq_u16_u32(c0.val[0])),
                     vget_high_u16(vreinterpretq_u16_u32(c2.val[0])));

  *a2 = vcombine_u16(vget_low_u16(vreinterpretq_u16_u32(c0.val[1])),
                     vget_low_u16(vreinterpretq_u16_u32(c2.val[1])));
  *a6 = vcombine_u16(vget_high_u16(vreinterpretq_u16_u32(c0.val[1])),
                     vget_high_u16(vreinterpretq_u16_u32(c2.val[1])));

  *a1 = vcombine_u16(vget_low_u16(vreinterpretq_u16_u32(c1.val[0])),
                     vget_low_u16(vreinterpretq_u16_u32(c3.val[0])));
  *a5 = vcombine_u16(vget_high_u16(vreinterpretq_u16_u32(c1.val[0])),
                     vget_high_u16(vreinterpretq_u16_u32(c3.val[0])));

  *a3 = vcombine_u16(vget_low_u16(vreinterpretq_u16_u32(c1.val[1])),
                     vget_low_u16(vreinterpretq_u16_u32(c3.val[1])));
  *a7 = vcombine_u16(vget_high_u16(vreinterpretq_u16_u32(c1.val[1])),
                     vget_high_u16(vreinterpretq_u16_u32(c3.val[1])));
}

static INLINE uint8x8_t wiener_convolve_s16_8_8(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
    const int16x8_t s6, const int16x8_t s7, const int16x8_t filters,
    const int bd, const int round1_bits) {
  const int16x4_t filters_lo = vget_low_s16(filters);
  int16x8_t ss0, ss1, ss2;
  int32x4_t sum0, sum1, zero;
  uint16x4_t tmp0, tmp1;
  int16x4_t filter3;
  uint16x8_t tmp;
  uint8x8_t res;
  int32x4_t round_vec;
  int32_t round_const = (1 << (bd + round1_bits - 1));
  const int x = -round1_bits;
  int32x4_t round_bits = vdupq_n_s32(x);
  (void)s7;

  zero = vdupq_n_s32(0);
  round_vec = vdupq_n_s32(round_const);
  filter3 =
      vadd_s16(vdup_lane_s16(filters_lo, 3), vdup_n_s16(1 << FILTER_BITS));

  ss0 = vaddq_s16(s0, s6);
  ss1 = vaddq_s16(s1, s5);
  ss2 = vaddq_s16(s2, s4);

  sum0 = vmull_lane_s16(vget_low_s16(ss0), filters_lo, 0);
  sum0 = vmlal_lane_s16(sum0, vget_low_s16(ss1), filters_lo, 1);
  sum0 = vmlal_lane_s16(sum0, vget_low_s16(ss2), filters_lo, 2);
  sum0 = vmlal_lane_s16(sum0, vget_low_s16(s3), filter3, 3);

  sum1 = vmull_lane_s16(vget_high_s16(ss0), filters_lo, 0);
  sum1 = vmlal_lane_s16(sum1, vget_high_s16(ss1), filters_lo, 1);
  sum1 = vmlal_lane_s16(sum1, vget_high_s16(ss2), filters_lo, 2);
  sum1 = vmlal_lane_s16(sum1, vget_high_s16(s3), filter3, 3);

  sum0 = vsubq_s32(sum0, round_vec);
  sum1 = vsubq_s32(sum1, round_vec);

  /*right shift & rounding & saturating*/
  // sum0 = vqrshlq_s32(sum0, round_bits);
  // sum1 = vqrshlq_s32(sum1, round_bits);

  sum0 = vrshlq_s32(sum0, round_bits);
  sum1 = vrshlq_s32(sum1, round_bits);

  sum0 = vmaxq_s32(sum0, zero);
  sum1 = vmaxq_s32(sum1, zero);

  /*Converting from int32x4_t to uint8x8_t*/
  tmp0 = vqmovn_u32(vreinterpretq_u32_s32(sum0));
  tmp1 = vqmovn_u32(vreinterpretq_u32_s32(sum1));
  tmp = vcombine_u16(tmp0, tmp1);
  res = vqmovn_u16(tmp);

  return res;
}

static INLINE uint16x8_t wiener_convolve_u8_8_8(
    const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
    const int16x8_t s3, const int16x4_t filters_lo, const int bd,
    const int round0_bits) {
  int16x8_t sum;
  uint16x8_t res;
  int32x4_t sum_0, sum_1;
  int32x4_t s3_0, s3_1;
  int32x4_t round_vec_0, round_vec_1;
  int32_t round_const_0 = (1 << (bd + FILTER_BITS - 1));
  int32_t round_const_1 = (1 << ((bd) + 1 + FILTER_BITS - round0_bits));

  /*for the purpose of right shift by { conv_params->round_0 }*/
  const int x = -round0_bits;
  int32x4_t round_bits = vdupq_n_s32(x);

  /*adding 128 to the coeef_3 as it is being used as rounding factor*/
  const int16x4_t filter3 =
      vadd_s16(vdup_lane_s16(filters_lo, 3), vdup_n_s16(1 << FILTER_BITS));

  round_vec_0 = vdupq_n_s32(round_const_0);
  round_vec_1 = vdupq_n_s32(round_const_1);

  sum = vmulq_lane_s16(s0, filters_lo, 0);
  sum = vmlaq_lane_s16(sum, s1, filters_lo, 1);
  sum = vmlaq_lane_s16(sum, s2, filters_lo, 2);

  /* converting sum from 16x8 to 2 32x4 registers*/
  sum_0 = vmovl_s16(vget_low_s16(sum));
  sum_1 = vmovl_s16(vget_high_s16(sum));

  /*  s[3]*128 -- and filter coff max can be 128
          them max value possible = 128*128*255 exceeding 16 bit*/

  s3_0 = vmull_s16(vget_low_s16(s3), filter3);
  s3_1 = vmull_s16(vget_high_s16(s3), filter3);
  sum_0 = vaddq_s32(sum_0, s3_0);
  sum_1 = vaddq_s32(sum_1, s3_1);

  /*Adding the rounding value*/
  sum_0 = vaddq_s32(sum_0, round_vec_0);
  sum_1 = vaddq_s32(sum_1, round_vec_0);

  /*right shift & rounding & saturating*/
  sum_0 = vqrshlq_s32(sum_0, round_bits);
  sum_1 = vqrshlq_s32(sum_1, round_bits);

  /*Clipping to max value*/
  sum_0 = vminq_s32(sum_0, round_vec_1);
  sum_1 = vminq_s32(sum_1, round_vec_1);

  res = vcombine_u16(vqmovun_s32(sum_0), vqmovun_s32(sum_1));
  return res;
}

static INLINE uint16x4_t wiener_convolve_u8_8_4(
    const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
    const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
    const int16x4_t s6, const int16x4_t s7, const int16x4_t filters_lo,
    const int bd, const int round0_bits) {
  int32x4_t sum_0, zero, round_vec_0, round_vec_1, s3_0;
  int16x4_t sum, temp0, temp1, temp2;
  uint16x4_t res;

  int32_t round_const_0 = (1 << (bd + FILTER_BITS - 1));
  int32_t round_const_1 = (1 << ((bd) + 1 + FILTER_BITS - round0_bits));
  const int x = -round0_bits;
  int32x4_t round_bits = vdupq_n_s32(x);
  const int16x4_t filter3 =
      vadd_s16(vdup_lane_s16(filters_lo, 3), vdup_n_s16(1 << FILTER_BITS));
  (void)s7;

  temp0 = vadd_s16(s0, s6);
  temp1 = vadd_s16(s1, s5);
  temp2 = vadd_s16(s2, s4);

  round_vec_0 = vdupq_n_s32(round_const_0);
  round_vec_1 = vdupq_n_s32(round_const_1);
  zero = vdupq_n_s32(0);

  sum = vmul_lane_s16(temp0, filters_lo, 0);
  sum = vmla_lane_s16(sum, temp1, filters_lo, 1);
  sum = vmla_lane_s16(sum, temp2, filters_lo, 2);
  sum_0 = vmovl_s16(sum);

  /*  s[3]*128 -- and filter coff max can be 128
      128*128*255 */
  s3_0 = vmull_s16(s3, filter3);
  sum_0 = vaddq_s32(sum_0, s3_0);

  /*Adding the rounding value*/
  sum_0 = vaddq_s32(sum_0, round_vec_0);

  /*right & rounding*/
  sum_0 = vrshlq_s32(sum_0, round_bits);

  sum_0 = vmaxq_s32(sum_0, zero);
  sum_0 = vminq_s32(sum_0, round_vec_1);
  res = vqmovun_s32(sum_0);
  return res;
}

static INLINE void store_u16_8x8(uint16_t *s, ptrdiff_t dst_stride,
                                 const uint16x8_t s0, const uint16x8_t s1,
                                 const uint16x8_t s2, const uint16x8_t s3,
                                 const uint16x8_t s4, const uint16x8_t s5,
                                 const uint16x8_t s6, const uint16x8_t s7) {
  vst1q_u16(s, s0);
  s += dst_stride;
  vst1q_u16(s, s1);
  s += dst_stride;
  vst1q_u16(s, s2);
  s += dst_stride;
  vst1q_u16(s, s3);
  s += dst_stride;
  vst1q_u16(s, s4);
  s += dst_stride;
  vst1q_u16(s, s5);
  s += dst_stride;
  vst1q_u16(s, s6);
  s += dst_stride;
  vst1q_u16(s, s7);
}

static INLINE void store_u16_8x4(uint16_t *s, ptrdiff_t dst_stride,
                                 const uint16x8_t s0, const uint16x8_t s1,
                                 const uint16x8_t s2, const uint16x8_t s3) {
  vst1q_u16(s, s0);
  s += dst_stride;
  vst1q_u16(s, s1);
  s += dst_stride;
  vst1q_u16(s, s2);
  s += dst_stride;
  vst1q_u16(s, s3);
}

#endif  // AV1_CONVOLVE_NEON_H_