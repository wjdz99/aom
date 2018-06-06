/*
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

#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/blend.h"
#include "aom_ports/mem.h"
#include "av1/common/arm/mem_neon.h"
#include "config/aom_dsp_rtcd.h"

static INLINE void blend_w8(uint8_t *dst, uint32_t dst_stride,
                            const CONV_BUF_TYPE *src0, uint32_t src0_stride,
                            const CONV_BUF_TYPE *src1, uint32_t src1_stride,
                            int16x8_t mask0, int16x8_t mask1, int16x8_t mask2,
                            int16x8_t mask3, const int16x4_t v_maxval,
                            const uint16x8_t vec_round_offset,
                            const int16x8_t vec_round_bits) {
  int16x8_t src0_0, src0_1, src0_2, src0_3;
  int16x8_t src1_0, src1_1, src1_2, src1_3;
  int16x4_t max_minus_mask;

  load_s16_8x4((int16_t *)src0, (int32_t)src0_stride, &src0_0, &src0_1, &src0_2,
               &src0_3);
  load_s16_8x4((int16_t *)src1, (int32_t)src1_stride, &src1_0, &src1_1, &src1_2,
               &src1_3);

  int32x4_t m0_lowxsrc0_0_low =
      vmull_s16(vget_low_s16(mask0), vget_low_s16(src0_0));

  max_minus_mask = vsub_s16(v_maxval, vget_low_s16(mask0));

  m0_lowxsrc0_0_low =
      vmlal_s16(m0_lowxsrc0_0_low, max_minus_mask, vget_low_s16(src1_0));

  int32x4_t m1_lowxsrc0_1_low =
      vmull_s16(vget_low_s16(mask1), vget_low_s16(src0_1));

  max_minus_mask = vsub_s16(v_maxval, vget_low_s16(mask1));

  m1_lowxsrc0_1_low =
      vmlal_s16(m1_lowxsrc0_1_low, max_minus_mask, vget_low_s16(src1_1));

  int32x4_t m2_lowxsrc0_2_low =
      vmull_s16(vget_low_s16(mask2), vget_low_s16(src0_2));

  max_minus_mask = vsub_s16(v_maxval, vget_low_s16(mask2));

  m2_lowxsrc0_2_low =
      vmlal_s16(m2_lowxsrc0_2_low, max_minus_mask, vget_low_s16(src1_2));

  int32x4_t m3_lowxsrc0_3_low =
      vmull_s16(vget_low_s16(mask3), vget_low_s16(src0_3));

  max_minus_mask = vsub_s16(v_maxval, vget_low_s16(mask3));

  m3_lowxsrc0_3_low =
      vmlal_s16(m3_lowxsrc0_3_low, max_minus_mask, vget_low_s16(src1_3));

  int32x4_t m0_highxsrc0_0_high =
      vmull_s16(vget_high_s16(mask0), vget_high_s16(src0_0));

  max_minus_mask = vsub_s16(v_maxval, vget_high_s16(mask0));

  m0_highxsrc0_0_high =
      vmlal_s16(m0_highxsrc0_0_high, max_minus_mask, vget_high_s16(src1_0));

  int32x4_t m1_highxsrc0_1_high =
      vmull_s16(vget_high_s16(mask1), vget_high_s16(src0_1));

  max_minus_mask = vsub_s16(v_maxval, vget_high_s16(mask1));

  m1_highxsrc0_1_high =
      vmlal_s16(m1_highxsrc0_1_high, max_minus_mask, vget_high_s16(src1_1));

  int32x4_t m2_highxsrc0_2_high =
      vmull_s16(vget_high_s16(mask2), vget_high_s16(src0_2));

  max_minus_mask = vsub_s16(v_maxval, vget_high_s16(mask2));

  m2_highxsrc0_2_high =
      vmlal_s16(m2_highxsrc0_2_high, max_minus_mask, vget_high_s16(src1_2));

  int32x4_t m3_highxsrc0_3_high =
      vmull_s16(vget_high_s16(mask3), vget_high_s16(src0_3));

  max_minus_mask = vsub_s16(v_maxval, vget_high_s16(mask3));

  m3_highxsrc0_3_high =
      vmlal_s16(m3_highxsrc0_3_high, max_minus_mask, vget_high_s16(src1_3));

  m0_lowxsrc0_0_low = vshrq_n_s32(m0_lowxsrc0_0_low, AOM_BLEND_A64_ROUND_BITS);
  m1_lowxsrc0_1_low = vshrq_n_s32(m1_lowxsrc0_1_low, AOM_BLEND_A64_ROUND_BITS);
  m2_lowxsrc0_2_low = vshrq_n_s32(m2_lowxsrc0_2_low, AOM_BLEND_A64_ROUND_BITS);
  m3_lowxsrc0_3_low = vshrq_n_s32(m3_lowxsrc0_3_low, AOM_BLEND_A64_ROUND_BITS);

  m0_highxsrc0_0_high =
      vshrq_n_s32(m0_highxsrc0_0_high, AOM_BLEND_A64_ROUND_BITS);
  m1_highxsrc0_1_high =
      vshrq_n_s32(m1_highxsrc0_1_high, AOM_BLEND_A64_ROUND_BITS);
  m2_highxsrc0_2_high =
      vshrq_n_s32(m2_highxsrc0_2_high, AOM_BLEND_A64_ROUND_BITS);
  m3_highxsrc0_3_high =
      vshrq_n_s32(m3_highxsrc0_3_high, AOM_BLEND_A64_ROUND_BITS);

  int16x8_t im_res_0 = vcombine_s16(vmovn_s32(m0_lowxsrc0_0_low),
                                    vmovn_s32(m0_highxsrc0_0_high));
  int16x8_t im_res_1 = vcombine_s16(vmovn_s32(m1_lowxsrc0_1_low),
                                    vmovn_s32(m1_highxsrc0_1_high));
  int16x8_t im_res_2 = vcombine_s16(vmovn_s32(m2_lowxsrc0_2_low),
                                    vmovn_s32(m2_highxsrc0_2_high));
  int16x8_t im_res_3 = vcombine_s16(vmovn_s32(m3_lowxsrc0_3_low),
                                    vmovn_s32(m3_highxsrc0_3_high));

  uint16x8_t im_res1_0 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_0), vec_round_offset);
  uint16x8_t im_res1_1 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_1), vec_round_offset);
  uint16x8_t im_res1_2 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_2), vec_round_offset);
  uint16x8_t im_res1_3 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_3), vec_round_offset);

  im_res_0 = vshlq_s16(vreinterpretq_s16_u16(im_res1_0), vec_round_bits);
  im_res_1 = vshlq_s16(vreinterpretq_s16_u16(im_res1_1), vec_round_bits);
  im_res_2 = vshlq_s16(vreinterpretq_s16_u16(im_res1_2), vec_round_bits);
  im_res_3 = vshlq_s16(vreinterpretq_s16_u16(im_res1_3), vec_round_bits);

  vst1_u8((dst + 0 * dst_stride), vqmovun_s16(im_res_0));
  vst1_u8((dst + 1 * dst_stride), vqmovun_s16(im_res_1));
  vst1_u8((dst + 2 * dst_stride), vqmovun_s16(im_res_2));
  vst1_u8((dst + 3 * dst_stride), vqmovun_s16(im_res_3));
}

static INLINE void blend_w4(uint8_t *dst, uint32_t dst_stride,
                            const CONV_BUF_TYPE *src0, uint32_t src0_stride,
                            const CONV_BUF_TYPE *src1, uint32_t src1_stride,
                            int16x4_t mask0, int16x4_t mask1, int16x4_t mask2,
                            int16x4_t mask3, const int16x4_t v_maxval,
                            const uint16x8_t vec_round_offset,
                            const int16x8_t vec_round_bits) {
  const int16x4_t vec_zer0 = vdup_n_s16(0);
  int16x8_t src0_0, src0_1;
  int16x8_t src1_0, src1_1;
  uint64x2_t tu0, tu1, tu2, tu3;
  int16x4_t max_minus_mask;

  load_unaligned_u16_4x4(src0, src0_stride, &tu0, &tu1);
  load_unaligned_u16_4x4(src1, src1_stride, &tu2, &tu3);

  src0_0 = vreinterpretq_s16_u64(tu0);
  src0_1 = vreinterpretq_s16_u64(tu1);

  src1_0 = vreinterpretq_s16_u64(tu2);
  src1_1 = vreinterpretq_s16_u64(tu3);

  int32x4_t m0_lowxsrc0_0_low = vmull_s16(mask0, vget_low_s16(src0_0));
  max_minus_mask = vsub_s16(v_maxval, mask0);
  m0_lowxsrc0_0_low =
      vmlal_s16(m0_lowxsrc0_0_low, max_minus_mask, vget_low_s16(src1_0));

  int32x4_t m1_lowxsrc0_1_low = vmull_s16(mask1, vget_high_s16(src0_0));
  max_minus_mask = vsub_s16(v_maxval, mask1);
  m1_lowxsrc0_1_low =
      vmlal_s16(m1_lowxsrc0_1_low, max_minus_mask, vget_high_s16(src1_0));

  int32x4_t m2_lowxsrc0_2_low = vmull_s16(mask2, vget_low_s16(src0_1));
  max_minus_mask = vsub_s16(v_maxval, mask2);
  m2_lowxsrc0_2_low =
      vmlal_s16(m2_lowxsrc0_2_low, max_minus_mask, vget_low_s16(src1_1));

  int32x4_t m3_lowxsrc0_3_low = vmull_s16(mask3, vget_high_s16(src0_1));
  max_minus_mask = vsub_s16(v_maxval, mask3);
  m3_lowxsrc0_3_low =
      vmlal_s16(m3_lowxsrc0_3_low, max_minus_mask, vget_high_s16(src1_1));

  m0_lowxsrc0_0_low = vshrq_n_s32(m0_lowxsrc0_0_low, AOM_BLEND_A64_ROUND_BITS);
  m1_lowxsrc0_1_low = vshrq_n_s32(m1_lowxsrc0_1_low, AOM_BLEND_A64_ROUND_BITS);
  m2_lowxsrc0_2_low = vshrq_n_s32(m2_lowxsrc0_2_low, AOM_BLEND_A64_ROUND_BITS);
  m3_lowxsrc0_3_low = vshrq_n_s32(m3_lowxsrc0_3_low, AOM_BLEND_A64_ROUND_BITS);

  int16x8_t im_res_0 = vcombine_s16(vmovn_s32(m0_lowxsrc0_0_low), vec_zer0);
  int16x8_t im_res_1 = vcombine_s16(vmovn_s32(m1_lowxsrc0_1_low), vec_zer0);
  int16x8_t im_res_2 = vcombine_s16(vmovn_s32(m2_lowxsrc0_2_low), vec_zer0);
  int16x8_t im_res_3 = vcombine_s16(vmovn_s32(m3_lowxsrc0_3_low), vec_zer0);

  uint16x8_t im_res1_0 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_0), vec_round_offset);
  uint16x8_t im_res1_1 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_1), vec_round_offset);
  uint16x8_t im_res1_2 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_2), vec_round_offset);
  uint16x8_t im_res1_3 =
      vqsubq_u16(vreinterpretq_u16_s16(im_res_3), vec_round_offset);

  im_res_0 = vshlq_s16(vreinterpretq_s16_u16(im_res1_0), vec_round_bits);
  im_res_1 = vshlq_s16(vreinterpretq_s16_u16(im_res1_1), vec_round_bits);
  im_res_2 = vshlq_s16(vreinterpretq_s16_u16(im_res1_2), vec_round_bits);
  im_res_3 = vshlq_s16(vreinterpretq_s16_u16(im_res1_3), vec_round_bits);

  uint8x8_t res_0 = vqmovun_s16(im_res_0);
  uint8x8_t res_1 = vqmovun_s16(im_res_1);
  uint8x8_t res_2 = vqmovun_s16(im_res_2);
  uint8x8_t res_3 = vqmovun_s16(im_res_3);

  vst1_lane_u32((uint32_t *)(dst + 0 * dst_stride), vreinterpret_u32_u8(res_0),
                0);
  vst1_lane_u32((uint32_t *)(dst + 1 * dst_stride), vreinterpret_u32_u8(res_1),
                0);
  vst1_lane_u32((uint32_t *)(dst + 2 * dst_stride), vreinterpret_u32_u8(res_2),
                0);
  vst1_lane_u32((uint32_t *)(dst + 3 * dst_stride), vreinterpret_u32_u8(res_3),
                0);
}

void aom_lowbd_blend_a64_d16_mask_neon(
    uint8_t *dst, uint32_t dst_stride, const CONV_BUF_TYPE *src0,
    uint32_t src0_stride, const CONV_BUF_TYPE *src1, uint32_t src1_stride,
    const uint8_t *mask, uint32_t mask_stride, int h, int w, int subh, int subw,
    ConvolveParams *conv_params) {
  int i = 0, j = 0;
  const int bd = 8;
  int w_ = w;
  const uint8_t *mask_tmp = mask;
  const CONV_BUF_TYPE *src0_tmp = src0;
  const CONV_BUF_TYPE *src1_tmp = src1;
  uint8_t *dst_tmp = dst;

  const int offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
  const int round_offset = (1 << (offset_bits - conv_params->round_1)) +
                           (1 << (offset_bits - conv_params->round_1 - 1));
  const int round_bits =
      2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;

  assert(IMPLIES((void *)src0 == dst, src0_stride == dst_stride));
  assert(IMPLIES((void *)src1 == dst, src1_stride == dst_stride));

  assert(h >= 4);
  assert(w >= 4);
  assert(IS_POWER_OF_TWO(h));
  assert(IS_POWER_OF_TWO(w));

  uint8x8_t s0, s1, s2, s3;
  uint32x2_t tu0, tu1, tu2, tu3;
  uint8x16_t t0, t1, t2, t3, t4, t5, t6, t7;
  int16x8_t mask0, mask1, mask2, mask3;
  int16x8_t mask4, mask5, mask6, mask7;
  int32x4_t m0_32, m1_32, m2_32, m3_32;
  int32x4_t m4_32, m5_32, m6_32, m7_32;
  uint8x8_t mask0_l, mask1_l, mask2_l, mask3_l;
  uint8x8_t mask4_l, mask5_l, mask6_l, mask7_l;
  int16x4_t mask0_low, mask1_low, mask2_low, mask3_low;
  const uint16x4_t vec_zer0 = vdup_n_u16(0);
  const uint16_t offset = round_offset - (1 << (round_bits - 1));
  const int16x4_t v_maxval = vdup_n_s16(AOM_BLEND_A64_MAX_ALPHA);
  const int16x8_t vec_round_bits = vdupq_n_s16(-round_bits);
  const uint16x8_t vec_offset = vdupq_n_u16(offset);

  if (subw == 0 && subh == 0) {
    if (w_ > 7) {
      do {
        w_ = w;

        do {
          load_u8_8x4(mask_tmp, mask_stride, &s0, &s1, &s2, &s3);

          mask0 = vmovl_s8(vreinterpret_s8_u8(s0));
          mask1 = vmovl_s8(vreinterpret_s8_u8(s1));
          mask2 = vmovl_s8(vreinterpret_s8_u8(s2));
          mask3 = vmovl_s8(vreinterpret_s8_u8(s3));

          blend_w8(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0, mask1, mask2, mask3, v_maxval,
                   vec_offset, vec_round_bits);

          w_ -= 8;
          j += 8;
          mask_tmp += 8;
          dst_tmp += 8;
          src0_tmp += 8;
          src1_tmp += 8;
        } while (w_ > 7);
        j = 0;
        i += 4;
        mask_tmp = mask + (i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }

    if (w_ > 3) {
      do {
        w_ = w;

        do {
          load_unaligned_u8_4x4(mask_tmp, mask_stride, &tu0, &tu1);

          mask0 = vreinterpretq_s16_u16(vmovl_u8(vreinterpret_u8_u32(tu0)));
          mask1 = vreinterpretq_s16_u16(vmovl_u8(vreinterpret_u8_u32(tu1)));

          mask0_low = vget_low_s16(mask0);
          mask1_low = vget_high_s16(mask0);
          mask2_low = vget_low_s16(mask1);
          mask3_low = vget_high_s16(mask1);

          blend_w4(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0_low, mask1_low, mask2_low, mask3_low,
                   v_maxval, vec_offset, vec_round_bits);
          w_ -= 4;
          j += 4;
          mask_tmp += 4;
          dst_tmp += 4;
          src0_tmp += 4;
          src1_tmp += 4;
        } while (w_ > 3);
        j = 0;
        i += 4;
        mask_tmp = mask + (i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
  } else if (subw == 1 && subh == 1) {
    if (w_ > 7) {
      do {
        w_ = w;

        do {
          load_u8_16x8(mask_tmp, mask_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6,
                       &t7);

          mask0 =
              vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t0), vget_low_u8(t1)));
          mask1 =
              vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t2), vget_low_u8(t3)));
          mask2 =
              vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t4), vget_low_u8(t5)));
          mask3 =
              vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t6), vget_low_u8(t7)));

          mask4 = vreinterpretq_s16_u16(
              vaddl_u8(vget_high_u8(t0), vget_high_u8(t1)));
          mask5 = vreinterpretq_s16_u16(
              vaddl_u8(vget_high_u8(t2), vget_high_u8(t3)));
          mask6 = vreinterpretq_s16_u16(
              vaddl_u8(vget_high_u8(t4), vget_high_u8(t5)));
          mask7 = vreinterpretq_s16_u16(
              vaddl_u8(vget_high_u8(t6), vget_high_u8(t7)));

          m0_32 = vpaddlq_s16(mask0);
          m1_32 = vpaddlq_s16(mask1);
          m2_32 = vpaddlq_s16(mask2);
          m3_32 = vpaddlq_s16(mask3);

          m4_32 = vpaddlq_s16(mask4);
          m5_32 = vpaddlq_s16(mask5);
          m6_32 = vpaddlq_s16(mask6);
          m7_32 = vpaddlq_s16(mask7);

          mask0 =
              vcombine_s16(vqrshrn_n_s32(m0_32, 2), vqrshrn_n_s32(m4_32, 2));
          mask1 =
              vcombine_s16(vqrshrn_n_s32(m1_32, 2), vqrshrn_n_s32(m5_32, 2));
          mask2 =
              vcombine_s16(vqrshrn_n_s32(m2_32, 2), vqrshrn_n_s32(m6_32, 2));
          mask3 =
              vcombine_s16(vqrshrn_n_s32(m3_32, 2), vqrshrn_n_s32(m7_32, 2));

          blend_w8(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0, mask1, mask2, mask3, v_maxval,
                   vec_offset, vec_round_bits);

          w_ -= 8;
          j += 8;
          mask_tmp += 16;
          dst_tmp += 8;
          src0_tmp += 8;
          src1_tmp += 8;
        } while (w_ > 7);
        j = 0;
        i += 4;
        mask_tmp = mask + (2 * i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
    if (w_ > 3) {
      do {
        w_ = w;

        do {
          load_u8_8x8(mask_tmp, mask_stride, &mask0_l, &mask1_l, &mask2_l,
                      &mask3_l, &mask4_l, &mask5_l, &mask6_l, &mask7_l);

          mask0 = vreinterpretq_s16_u16(vaddl_u8(mask0_l, mask1_l));
          mask1 = vreinterpretq_s16_u16(vaddl_u8(mask2_l, mask3_l));
          mask2 = vreinterpretq_s16_u16(vaddl_u8(mask4_l, mask5_l));
          mask3 = vreinterpretq_s16_u16(vaddl_u8(mask6_l, mask7_l));

          m0_32 = vpaddlq_s16(mask0);
          m1_32 = vpaddlq_s16(mask1);
          m2_32 = vpaddlq_s16(mask2);
          m3_32 = vpaddlq_s16(mask3);

          mask0_low = vqrshrn_n_s32(m0_32, 2);
          mask1_low = vqrshrn_n_s32(m1_32, 2);
          mask2_low = vqrshrn_n_s32(m2_32, 2);
          mask3_low = vqrshrn_n_s32(m3_32, 2);

          blend_w4(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0_low, mask1_low, mask2_low, mask3_low,
                   v_maxval, vec_offset, vec_round_bits);

          w_ -= 4;
          j += 4;
          mask_tmp += 8;
          dst_tmp += 4;
          src0_tmp += 4;
          src1_tmp += 4;
        } while (w_ > 3);
        j = 0;
        i += 4;
        mask_tmp = mask + (2 * i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
  } else if (subw == 1 && subh == 0) {
    if (w_ > 7) {
      do {
        w_ = w;
        do {
          load_u8_16x4(mask_tmp, mask_stride, &t0, &t1, &t2, &t3);

          mask0 = vreinterpretq_s16_u16(vcombine_u16(
              vpaddl_u8(vget_low_u8(t0)), vpaddl_u8(vget_high_u8(t0))));
          mask1 = vreinterpretq_s16_u16(vcombine_u16(
              vpaddl_u8(vget_low_u8(t1)), vpaddl_u8(vget_high_u8(t1))));
          mask2 = vreinterpretq_s16_u16(vcombine_u16(
              vpaddl_u8(vget_low_u8(t2)), vpaddl_u8(vget_high_u8(t2))));
          mask3 = vreinterpretq_s16_u16(vcombine_u16(
              vpaddl_u8(vget_low_u8(t3)), vpaddl_u8(vget_high_u8(t3))));

          mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
          mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));
          mask2 = vmovl_s8(vqrshrn_n_s16(mask2, 1));
          mask3 = vmovl_s8(vqrshrn_n_s16(mask3, 1));

          blend_w8(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0, mask1, mask2, mask3, v_maxval,
                   vec_offset, vec_round_bits);
          w_ -= 8;
          j += 8;
          mask_tmp += 16;
          dst_tmp += 8;
          src0_tmp += 8;
          src1_tmp += 8;
        } while (w_ > 7);
        j = 0;
        i += 4;
        mask_tmp = mask + (i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }

    if (w_ > 3) {
      do {
        w_ = w;
        do {
          load_u8_8x4(mask_tmp, mask_stride, &mask0_l, &mask1_l, &mask2_l,
                      &mask3_l);

          mask0 =
              vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask0_l), vec_zer0));
          mask1 =
              vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask1_l), vec_zer0));
          mask2 =
              vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask2_l), vec_zer0));
          mask3 =
              vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask3_l), vec_zer0));

          mask0_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask0, 1)));
          mask1_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask1, 1)));
          mask2_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask2, 1)));
          mask3_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask3, 1)));

          blend_w4(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0_low, mask1_low, mask2_low, mask3_low,
                   v_maxval, vec_offset, vec_round_bits);

          w_ -= 4;
          j += 4;
          mask_tmp += 8;
          dst_tmp += 4;
          src0_tmp += 4;
          src1_tmp += 4;
        } while (w_ > 3);
        j = 0;
        i += 4;
        mask_tmp = mask + (i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
  } else {
    if (w_ > 7) {
      do {
        w_ = w;
        do {
          load_u8_8x8(mask_tmp, mask_stride, &mask0_l, &mask1_l, &mask2_l,
                      &mask3_l, &mask4_l, &mask5_l, &mask6_l, &mask7_l);

          mask0 = vreinterpretq_s16_u16(vaddl_u8(mask0_l, mask1_l));
          mask1 = vreinterpretq_s16_u16(vaddl_u8(mask2_l, mask3_l));
          mask2 = vreinterpretq_s16_u16(vaddl_u8(mask4_l, mask5_l));
          mask3 = vreinterpretq_s16_u16(vaddl_u8(mask6_l, mask7_l));

          mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
          mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));
          mask2 = vmovl_s8(vqrshrn_n_s16(mask2, 1));
          mask3 = vmovl_s8(vqrshrn_n_s16(mask3, 1));

          blend_w8(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0, mask1, mask2, mask3, v_maxval,
                   vec_offset, vec_round_bits);

          w_ -= 8;
          j += 8;
          mask_tmp += 8;
          dst_tmp += 8;
          src0_tmp += 8;
          src1_tmp += 8;
        } while (w_ > 7);
        j = 0;
        i += 4;
        mask_tmp = mask + (2 * i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
    if (w_ > 3) {
      do {
        w_ = w;
        do {
          load_unaligned_u8_4x4(mask_tmp, 2 * mask_stride, &tu0, &tu1);
          load_unaligned_u8_4x4(mask_tmp + mask_stride, 2 * mask_stride, &tu2,
                                &tu3);

          s0 = vreinterpret_u8_u32(tu0);
          s1 = vreinterpret_u8_u32(tu1);
          s2 = vreinterpret_u8_u32(tu2);
          s3 = vreinterpret_u8_u32(tu3);

          mask0 = vreinterpretq_s16_u16(vaddl_u8(s0, s2));
          mask1 = vreinterpretq_s16_u16(vaddl_u8(s1, s3));

          mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
          mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));

          mask0_low = vget_low_s16(mask0);
          mask1_low = vget_high_s16(mask0);
          mask2_low = vget_low_s16(mask1);
          mask3_low = vget_high_s16(mask1);

          blend_w4(dst_tmp, dst_stride, src0_tmp, src0_stride, src1_tmp,
                   src1_stride, mask0_low, mask1_low, mask2_low, mask3_low,
                   v_maxval, vec_offset, vec_round_bits);

          w_ -= 4;
          j += 4;
          mask_tmp += 4;
          dst_tmp += 4;
          src0_tmp += 4;
          src1_tmp += 4;
        } while (w_ > 3);
        j = 0;
        i += 4;
        mask_tmp = mask + (2 * i * mask_stride);
        dst_tmp = dst + (i * dst_stride);
        src0_tmp = src0 + (i * src0_stride);
        src1_tmp = src1 + (i * src1_stride);
      } while (i < h);
    }
  }
}
