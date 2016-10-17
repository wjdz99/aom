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

#include "./aom_dsp_rtcd.h"
#include "aom_ports/mem.h"

// delta = 4/16 * clamp(a - o, -s, s) + 1/16 * clamp(b - o, -s, s) +
//         3/16 * clamp(c - o, -s, s) + 3/16 * clamp(d - o, -s, s) +
//         1/16 * clamp(e - o, -s, s) + 4/16 * clamp(f - o, -s, s)
SIMD_INLINE v128 calc_delta(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                            v128 f, v128 sp, v128 sm) {
  // The difference will be 9 bit, offset by 128 so we can use saturated
  // sub to avoid going to 16 bit temporarily before "strength" clipping.
  const v128 c128 = v128_dup_8(128);
  const v128 x = v128_add_8(c128, o);
  const v128 c8 = v128_dup_8(8);
  const v128 tmp = v128_add_8(
      v128_max_s8(v128_min_s8(v128_ssub_s8(v128_add_8(c128, c), x), sp), sm),
      v128_max_s8(v128_min_s8(v128_ssub_s8(v128_add_8(c128, d), x), sp), sm));
  const v128 delta = v128_add_8(
      v128_add_8(
          v128_shl_8(
              v128_add_8(
                  v128_max_s8(
                      v128_min_s8(v128_ssub_s8(v128_add_8(c128, a), x), sp),
                      sm),
                  v128_max_s8(
                      v128_min_s8(v128_ssub_s8(v128_add_8(c128, f), x), sp),
                      sm)),
              2),
          v128_add_8(
              v128_max_s8(v128_min_s8(v128_ssub_s8(v128_add_8(c128, b), x), sp),
                          sm),
              v128_max_s8(v128_min_s8(v128_ssub_s8(v128_add_8(c128, e), x), sp),
                          sm))),
      v128_add_8(v128_add_8(tmp, tmp), tmp));
  return v128_add_8(
      o,
      v128_shr_s8(
          v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
          4));
}

// v256 version of calc_delta()
SIMD_INLINE v256 calc_delta256(v256 o, v256 a, v256 b, v256 c, v256 d, v256 e,
                               v256 f, v256 sp, v256 sm) {
  // The difference will be 9 bit, offset by 128 so we can use saturated
  // sub to avoid going to 16 bit temporarily before "strength" clipping.
  const v256 c128 = v256_dup_8(128);
  const v256 x = v256_add_8(c128, o);
  const v256 c8 = v256_dup_8(8);
  const v256 tmp = v256_add_8(
      v256_max_s8(v256_min_s8(v256_ssub_s8(v256_add_8(c128, c), x), sp), sm),
      v256_max_s8(v256_min_s8(v256_ssub_s8(v256_add_8(c128, d), x), sp), sm));
  const v256 delta = v256_add_8(
      v256_add_8(
          v256_shl_8(
              v256_add_8(
                  v256_max_s8(
                      v256_min_s8(v256_ssub_s8(v256_add_8(c128, a), x), sp),
                      sm),
                  v256_max_s8(
                      v256_min_s8(v256_ssub_s8(v256_add_8(c128, f), x), sp),
                      sm)),
              2),
          v256_add_8(
              v256_max_s8(v256_min_s8(v256_ssub_s8(v256_add_8(c128, b), x), sp),
                          sm),
              v256_max_s8(v256_min_s8(v256_ssub_s8(v256_add_8(c128, e), x), sp),
                          sm))),
      v256_add_8(v256_add_8(tmp, tmp), tmp));
  return v256_add_8(
      o,
      v256_shr_s8(
          v256_add_8(c8, v256_add_8(delta, v256_cmplt_s8(delta, v256_zero()))),
          4));
}

// Process blocks of width 4, four lines at a time, 8 bit.
static void clpf_block4x4(const uint8_t *src, uint8_t *dst, int sstride,
                          int dstride, int x0, int y0, int sizey, int width,
                          int height, unsigned int strength) {
  const v128 sp = v128_dup_8(strength);
  const v128 sm = v128_dup_8(-(int)strength);
  const int right = width - 4 - x0;
  const int bottom = height - 4 - y0;
  DECLARE_ALIGNED(16, static const uint64_t,
                  b_shuff[]) = { 0x0504040401000000LL, 0x0d0c0c0c09080808LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  c_shuff[]) = { 0x0605040402010000LL, 0x0e0d0c0c0a090808LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  d_shuff[]) = { 0x0707060503030201LL, 0x0f0f0e0d0b0b0a09LL };
  DECLARE_ALIGNED(16, static const uint64_t,
                  e_shuff[]) = { 0x0707070603030302LL, 0x0f0f0f0e0b0b0b0aLL };
  int y;

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  for (y = 0; y < sizey; y += 4) {
    const uint32_t l0 = u32_load_aligned(src - (y != -y0) * sstride);
    const uint32_t l1 = u32_load_aligned(src);
    const uint32_t l2 = u32_load_aligned(src + sstride);
    const uint32_t l3 = u32_load_aligned(src + 2 * sstride);
    const uint32_t l4 = u32_load_aligned(src + 3 * sstride);
    const uint32_t l5 = u32_load_aligned(src + ((y != bottom) + 3) * sstride);
    v128 o = v128_from_32(l1, l2, l3, l4);
    const v128 a = v128_from_32(l0, l1, l2, l3);
    const v128 f = v128_from_32(l2, l3, l4, l5);
    v128 b, c, d, e;

    if (x0) {
      b = v128_from_32(u32_load_unaligned(src - 2),
                       u32_load_unaligned(src + sstride - 2),
                       u32_load_unaligned(src + 2 * sstride - 2),
                       u32_load_unaligned(src + 3 * sstride - 2));
      c = v128_from_32(u32_load_unaligned(src - 1),
                       u32_load_unaligned(src + sstride - 1),
                       u32_load_unaligned(src + 2 * sstride - 1),
                       u32_load_unaligned(src + 3 * sstride - 1));
    } else {  // Left clipping
      b = v128_shuffle_8(o, v128_load_aligned(b_shuff));
      c = v128_shuffle_8(o, v128_load_aligned(c_shuff));
    }
    if (right) {
      d = v128_from_32(u32_load_unaligned(src + 1),
                       u32_load_unaligned(src + sstride + 1),
                       u32_load_unaligned(src + 2 * sstride + 1),
                       u32_load_unaligned(src + 3 * sstride + 1));
      e = v128_from_32(u32_load_unaligned(src + 2 * !!right),
                       u32_load_unaligned(src + sstride + 2),
                       u32_load_unaligned(src + 2 * sstride + 2),
                       u32_load_unaligned(src + 3 * sstride + 2));
    } else {  // Right clipping
      d = v128_shuffle_8(o, v128_load_aligned(d_shuff));
      e = v128_shuffle_8(o, v128_load_aligned(e_shuff));
    }

    o = calc_delta(o, a, b, c, d, e, f, sp, sm);
    u32_store_aligned(dst, v128_low_u32(v128_shr_n_byte(o, 12)));
    u32_store_aligned(dst + dstride, v128_low_u32(v128_shr_n_byte(o, 8)));
    u32_store_aligned(dst + 2 * dstride, v128_low_u32(v128_shr_n_byte(o, 4)));
    u32_store_aligned(dst + 3 * dstride, v128_low_u32(o));

    dst += 4 * dstride;
    src += 4 * sstride;
  }
}

// Process blocks of width 8, four lines at a time, 8 bit.
static void clpf_block8x4(const uint8_t *src, uint8_t *dst, int sstride,
                          int dstride, int x0, int y0, int sizey, int width,
                          int height, unsigned int strength) {
  const v256 sp = v256_dup_8(strength);
  const v256 sm = v256_dup_8(-(int)strength);
  const int right = width - 8 - x0;
  const int bottom = height - 4 - y0;
  DECLARE_ALIGNED(32, static const uint64_t,
                  b_shuff[]) = { 0x0504030201000000LL, 0x0d0c0b0a09080808LL,
                                 0x0504030201000000LL, 0x0d0c0b0a09080808LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  c_shuff[]) = { 0x0605040302010000LL, 0x0e0d0c0b0a090808LL,
                                 0x0605040302010000LL, 0x0e0d0c0b0a090808LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  d_shuff[]) = { 0x0707060504030201LL, 0x0f0f0e0d0c0b0a09LL,
                                 0x0707060504030201LL, 0x0f0f0e0d0c0b0a09LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  e_shuff[]) = { 0x0707070605040302LL, 0x0f0f0f0e0d0c0b0aLL,
                                 0x0707070605040302LL, 0x0f0f0f0e0d0c0b0aLL };
  int y;

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  for (y = 0; y < sizey; y += 4) {
    const v64 l0 = v64_load_aligned(src - (y != -y0) * sstride);
    const v64 l1 = v64_load_aligned(src);
    const v64 l2 = v64_load_aligned(src + sstride);
    const v64 l3 = v64_load_aligned(src + 2 * sstride);
    const v64 l4 = v64_load_aligned(src + 3 * sstride);
    const v64 l5 = v64_load_aligned(src + ((y != bottom) + 3) * sstride);
    v256 o = v256_from_v64(l1, l2, l3, l4);
    const v256 a = v256_from_v64(l0, l1, l2, l3);
    v256 b = v256_from_v64(v64_load_unaligned(src - 2 * !!x0),
                           v64_load_unaligned(src + sstride - 2 * !!x0),
                           v64_load_unaligned(src + 2 * sstride - 2 * !!x0),
                           v64_load_unaligned(src + 3 * sstride - 2 * !!x0));
    v256 c = v256_from_v64(v64_load_unaligned(src - !!x0),
                           v64_load_unaligned(src + sstride - !!x0),
                           v64_load_unaligned(src + 2 * sstride - !!x0),
                           v64_load_unaligned(src + 3 * sstride - !!x0));
    v256 d = v256_from_v64(v64_load_unaligned(src + !!right),
                           v64_load_unaligned(src + sstride + !!right),
                           v64_load_unaligned(src + 2 * sstride + !!right),
                           v64_load_unaligned(src + 3 * sstride + !!right));
    v256 e = v256_from_v64(v64_load_unaligned(src + 2 * !!right),
                           v64_load_unaligned(src + sstride + 2 * !!right),
                           v64_load_unaligned(src + 2 * sstride + 2 * !!right),
                           v64_load_unaligned(src + 3 * sstride + 2 * !!right));
    const v256 f = v256_from_v64(l2, l3, l4, l5);

    if (!x0) {  // Left clipping
      b = v256_pshuffle_8(b, v256_load_aligned(b_shuff));
      c = v256_pshuffle_8(c, v256_load_aligned(c_shuff));
    }
    if (!right) {  // Right clipping
      d = v256_pshuffle_8(d, v256_load_aligned(d_shuff));
      e = v256_pshuffle_8(e, v256_load_aligned(e_shuff));
    }

    o = calc_delta256(o, a, b, c, d, e, f, sp, sm);
    v64_store_aligned(dst, v256_low_v64(v256_shr_n_byte(o, 24)));
    v64_store_aligned(dst + dstride, v256_low_v64(v256_shr_n_byte(o, 16)));
    v64_store_aligned(dst + 2 * dstride, v256_low_v64(v256_shr_n_byte(o, 8)));
    v64_store_aligned(dst + 3 * dstride, v256_low_v64(o));

    dst += 4 * dstride;
    src += 4 * sstride;
  }
}

void SIMD_FUNC(aom_clpf_block)(const uint8_t *src, uint8_t *dst, int sstride,
                               int dstride, int x0, int y0, int sizex,
                               int sizey, int width, int height,
                               unsigned int strength) {
  if ((sizex != 4 && sizex != 8) || ((sizey & 3) && sizex == 4)) {
    // Fallback to C for odd sizes:
    // * block widths not 4 or 8
    // * block heights not a multiple of 4 if the block width is 4
    aom_clpf_block_c(src, dst, sstride, dstride, x0, y0, sizex, sizey, width,
                     height, strength);
  } else {
    (sizex == 4 ? clpf_block4x4 : clpf_block8x4)(
        src, dst, sstride, dstride, x0, y0, sizey, width, height, strength);
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
// delta = 4/16 * clamp(a - o, -s, s) + 1/16 * clamp(b - o, -s, s) +
//         3/16 * clamp(c - o, -s, s) + 3/16 * clamp(d - o, -s, s) +
//         1/16 * clamp(e - o, -s, s) + 4/16 * clamp(f - o, -s, s)
SIMD_INLINE v128 calc_delta_hbd(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                                v128 f, v128 sp, v128 sm) {
  const v128 c8 = v128_dup_16(8);
  const v128 tmp =
      v128_add_16(v128_max_s16(v128_min_s16(v128_sub_16(c, o), sp), sm),
                  v128_max_s16(v128_min_s16(v128_sub_16(d, o), sp), sm));
  const v128 delta = v128_add_16(
      v128_add_16(
          v128_shl_16(
              v128_add_16(
                  v128_max_s16(v128_min_s16(v128_sub_16(a, o), sp), sm),
                  v128_max_s16(v128_min_s16(v128_sub_16(f, o), sp), sm)),
              2),
          v128_add_16(v128_max_s16(v128_min_s16(v128_sub_16(b, o), sp), sm),
                      v128_max_s16(v128_min_s16(v128_sub_16(e, o), sp), sm))),
      v128_add_16(v128_add_16(tmp, tmp), tmp));
  return v128_add_16(
      o, v128_shr_s16(
             v128_add_16(
                 c8, v128_add_16(delta, v128_cmplt_s16(delta, v128_zero()))),
             4));
}

// 256 bit version
SIMD_INLINE v256 calc_delta256_hbd(v256 o, v256 a, v256 b, v256 c, v256 d,
                                   v256 e, v256 f, v256 sp, v256 sm) {
  const v256 c8 = v256_dup_16(8);
  const v256 tmp =
      v256_add_16(v256_max_s16(v256_min_s16(v256_sub_16(c, o), sp), sm),
                  v256_max_s16(v256_min_s16(v256_sub_16(d, o), sp), sm));
  const v256 delta = v256_add_16(
      v256_add_16(
          v256_shl_16(
              v256_add_16(
                  v256_max_s16(v256_min_s16(v256_sub_16(a, o), sp), sm),
                  v256_max_s16(v256_min_s16(v256_sub_16(f, o), sp), sm)),
              2),
          v256_add_16(v256_max_s16(v256_min_s16(v256_sub_16(b, o), sp), sm),
                      v256_max_s16(v256_min_s16(v256_sub_16(e, o), sp), sm))),
      v256_add_16(v256_add_16(tmp, tmp), tmp));
  return v256_add_16(
      o, v256_shr_s16(
             v256_add_16(
                 c8, v256_add_16(delta, v256_cmplt_s16(delta, v256_zero()))),
             4));
}

// Process blocks of width 4, four lines at a time.
static void clpf_block4x4_hbd(const uint16_t *src, uint16_t *dst, int sstride,
                              int dstride, int x0, int y0, int sizey, int width,
                              int height, unsigned int strength) {
  const v256 sp = v256_dup_16(strength);
  const v256 sm = v256_dup_16(-(int)strength);
  const int right = width - 4 - x0;
  const int bottom = height - 4 - y0;
  DECLARE_ALIGNED(32, static const uint64_t,
                  b_shuff[]) = { 0x0302010001000100LL, 0x0b0a090809080908LL,
                                 0x0302010001000100LL, 0x0b0a090809080908LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  c_shuff[]) = { 0x0504030201000100LL, 0x0d0c0b0a09080908LL,
                                 0x0504030201000100LL, 0x0d0c0b0a09080908LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  d_shuff[]) = { 0x0706070605040302LL, 0x0f0e0f0e0d0c0b0aLL,
                                 0x0706070605040302LL, 0x0f0e0f0e0d0c0b0aLL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  e_shuff[]) = { 0x0706070607060504LL, 0x0f0e0f0e0f0e0d0cLL,
                                 0x0706070607060504LL, 0x0f0e0f0e0f0e0d0cLL };
  int y;

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  for (y = 0; y < sizey; y += 4) {
    const v64 l0 = v64_load_aligned(src - (y != -y0) * sstride);
    const v64 l1 = v64_load_aligned(src);
    const v64 l2 = v64_load_aligned(src + sstride);
    const v64 l3 = v64_load_aligned(src + 2 * sstride);
    const v64 l4 = v64_load_aligned(src + 3 * sstride);
    const v64 l5 = v64_load_aligned(src + ((y != bottom) + 3) * sstride);
    v256 o = v256_from_v64(l1, l2, l3, l4);
    const v256 a = v256_from_v64(l0, l1, l2, l3);
    const v256 f = v256_from_v64(l2, l3, l4, l5);
    v256 b, c, d, e;

    if (x0) {
      b = v256_from_v64(v64_load_unaligned(src - 2),
                        v64_load_unaligned(src + sstride - 2),
                        v64_load_unaligned(src + 2 * sstride - 2),
                        v64_load_unaligned(src + 3 * sstride - 2));
      c = v256_from_v64(v64_load_unaligned(src - 1),
                        v64_load_unaligned(src + sstride - 1),
                        v64_load_unaligned(src + 2 * sstride - 1),
                        v64_load_unaligned(src + 3 * sstride - 1));
    } else {  // Left clipping
      b = v256_pshuffle_8(o, v256_load_aligned(b_shuff));
      c = v256_pshuffle_8(o, v256_load_aligned(c_shuff));
    }
    if (right) {
      d = v256_from_v64(v64_load_unaligned(src + 1),
                        v64_load_unaligned(src + sstride + 1),
                        v64_load_unaligned(src + 2 * sstride + 1),
                        v64_load_unaligned(src + 3 * sstride + 1));
      e = v256_from_v64(v64_load_unaligned(src + 2 * !!right),
                        v64_load_unaligned(src + sstride + 2),
                        v64_load_unaligned(src + 2 * sstride + 2),
                        v64_load_unaligned(src + 3 * sstride + 2));
    } else {  // Right clipping
      d = v256_pshuffle_8(o, v256_load_aligned(d_shuff));
      e = v256_pshuffle_8(o, v256_load_aligned(e_shuff));
    }

    o = calc_delta256_hbd(o, a, b, c, d, e, f, sp, sm);
    v64_store_aligned(dst, v256_low_v64(v256_shr_n_byte(o, 24)));
    v64_store_aligned(dst + dstride, v256_low_v64(v256_shr_n_byte(o, 16)));
    v64_store_aligned(dst + 2 * dstride, v256_low_v64(v256_shr_n_byte(o, 8)));
    v64_store_aligned(dst + 3 * dstride, v256_low_v64(o));

    dst += 4 * dstride;
    src += 4 * sstride;
  }
}

// Process blocks of width 8, two lines at a time.
static void clpf_block8x2_hbd(const uint16_t *src, uint16_t *dst, int sstride,
                              int dstride, int x0, int y0, int sizey, int width,
                              int height, unsigned int strength) {
  const int right = width - 8 - x0;
  const int bottom = height - 2 - y0;
  const v256 sp = v256_dup_16(strength);
  const v256 sm = v256_dup_16(-(int)strength);
  DECLARE_ALIGNED(32, static const uint64_t,
                  b_shuff[]) = { 0x0302010001000100LL, 0x0b0a090807060504LL,
                                 0x0302010001000100LL, 0x0b0a090807060504LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  c_shuff[]) = { 0x0504030201000100LL, 0x0d0c0b0a09080706LL,
                                 0x0504030201000100LL, 0x0d0c0b0a09080706LL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  d_shuff[]) = { 0x0908070605040302LL, 0x0f0e0f0e0d0c0b0aLL,
                                 0x0908070605040302LL, 0x0f0e0f0e0d0c0b0aLL };
  DECLARE_ALIGNED(32, static const uint64_t,
                  e_shuff[]) = { 0x0b0a090807060504LL, 0x0f0e0f0e0f0e0d0cLL,
                                 0x0b0a090807060504LL, 0x0f0e0f0e0f0e0d0cLL };
  int y;

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  for (y = 0; y < sizey; y += 2) {
    const v128 l1 = v128_load_aligned(src);
    const v128 l2 = v128_load_aligned(src + sstride);
    v256 o = v256_from_v128(l1, l2);
    const v256 a =
        v256_from_v128(v128_load_aligned(src - (y != -y0) * sstride), l1);
    const v256 f = v256_from_v128(
        l2, v128_load_aligned(src + ((y != bottom) + 1) * sstride));
    v256 b, c, d, e;

    if (x0) {
      b = v256_from_v128(v128_load_unaligned(src - 2),
                         v128_load_unaligned(src - 2 + sstride));
      c = v256_from_v128(v128_load_unaligned(src - 1),
                         v128_load_unaligned(src - 1 + sstride));
    } else {  // Left clipping
      b = v256_pshuffle_8(o, v256_load_aligned(b_shuff));
      c = v256_pshuffle_8(o, v256_load_aligned(c_shuff));
    }
    if (right) {
      d = v256_from_v128(v128_load_unaligned(src + 1),
                         v128_load_unaligned(src + 1 + sstride));
      e = v256_from_v128(v128_load_unaligned(src + 2),
                         v128_load_unaligned(src + 2 + sstride));
    } else {  // Right clipping
      d = v256_pshuffle_8(o, v256_load_aligned(d_shuff));
      e = v256_pshuffle_8(o, v256_load_aligned(e_shuff));
    }

    o = calc_delta256_hbd(o, a, b, c, d, e, f, sp, sm);
    v128_store_aligned(dst, v256_high_v128(o));
    v128_store_aligned(dst + dstride, v256_low_v128(o));
    src += sstride * 2;
    dst += dstride * 2;
  }
}

void SIMD_FUNC(aom_clpf_block_hbd)(const uint16_t *src, uint16_t *dst,
                                   int sstride, int dstride, int x0, int y0,
                                   int sizex, int sizey, int width, int height,
                                   unsigned int strength) {
  if ((sizex != 4 && sizex != 8) || ((sizey & 1) && sizex == 4)) {
    // Fallback to C for odd sizes:
    // * block width not 4 or 8
    // * block heights not a multiple of 2 if the block width is 4
    aom_clpf_block_hbd_c(src, dst, sstride, dstride, x0, y0, sizex, sizey,
                         width, height, strength);
  } else {
    (sizex == 4 ? clpf_block4x4_hbd : clpf_block8x2_hbd)(
        src, dst, sstride, dstride, x0, y0, sizey, width, height, strength);
  }
}
#endif
