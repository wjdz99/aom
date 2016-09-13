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

// delta = 4/16 * clamp(a - o, -s, s) + 1/16 * clamp(b - o, -s, s) +
//         3/16 * clamp(c - o, -s, s) + 3/16 * clamp(d - o, -s, s) +
//         1/16 * clamp(e - o, -s, s) + 4/16 * clamp(f - o, -s, s)
SIMD_INLINE void calc_delta(v128 o, v128 x, v128 a, v128 b, v128 c, v128 d,
                            v128 e, v128 f, uint8_t *dst, v128 sp, v128 sm,
                            int dstride) {
  const v128 c8 = v128_dup_8(8);
  const v128 tmp =
      v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                 v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
  const v128 delta = v128_add_8(
      v128_add_8(
          v128_shl_8(
              v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                         v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
              2),
          v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                     v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
      v128_add_8(v128_add_8(tmp, tmp), tmp));
  o = v128_add_8(
      o,
      v128_shr_s8(
          v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
          4));
  v64_store_aligned(dst, v128_high_v64(o));
  v64_store_aligned(dst + dstride, v128_low_v64(o));
}

// Process 4 lines at a time, 8 bit.
static void clpf_block4(const uint8_t *src, uint8_t *dst, int sstride,
                        int dstride, int x0, int y0, int sizey, int width,
                        int height, unsigned int strength) {
  const int right = width - 4 - x0;
  const int bottom = height - 4 - y0;
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
    const v128 sp = v128_dup_8(strength);
    const v128 sm = v128_dup_8(-(int)strength);
    // The difference will be 9 bit, offset by 128 so we can use saturated
    // sub to avoid going to 16 bit temporarily before "strength" clipping.
    const v128 c128 = v128_dup_8(128);
    const v128 o = v128_from_32(l1, l2, l3, l4);
    const v128 x = v128_add_8(c128, o);
    const v128 a = v128_add_8(c128, v128_from_32(l0, l1, l2, l3));
    v128 b = v128_add_8(
        c128, v128_from_32(u32_load_unaligned(src - 2 * !!x0),
                           u32_load_unaligned(src + sstride - 2 * !!x0),
                           u32_load_unaligned(src + 2 * sstride - 2 * !!x0),
                           u32_load_unaligned(src + 3 * sstride - 2 * !!x0)));
    v128 c = v128_add_8(
        c128, v128_from_32(u32_load_unaligned(src - !!x0),
                           u32_load_unaligned(src + sstride - !!x0),
                           u32_load_unaligned(src + 2 * sstride - !!x0),
                           u32_load_unaligned(src + 3 * sstride - !!x0)));
    v128 d = v128_add_8(
        c128, v128_from_32(u32_load_unaligned(src + !!right),
                           u32_load_unaligned(src + sstride + !!right),
                           u32_load_unaligned(src + 2 * sstride + !!right),
                           u32_load_unaligned(src + 3 * sstride + !!right)));
    v128 e = v128_add_8(
        c128,
        v128_from_32(u32_load_unaligned(src + 2 * !!right),
                     u32_load_unaligned(src + sstride + 2 * !!right),
                     u32_load_unaligned(src + 2 * sstride + 2 * !!right),
                     u32_load_unaligned(src + 3 * sstride + 2 * !!right)));
    const v128 f = v128_add_8(c128, v128_from_32(l2, l3, l4, l5));
    v128 tmp, delta, r;

    if (!x0) {
      b = v128_shuffle_8(b, v128_from_v64(v64_from_64(0x0d0c0c0c09080808LL),
                                          v64_from_64(0x0504040401000000LL)));
      c = v128_shuffle_8(c, v128_from_v64(v64_from_64(0x0e0d0c0c0a090808LL),
                                          v64_from_64(0x0605040402010000LL)));
    }
    if (!right) {
      d = v128_shuffle_8(d, v128_from_v64(v64_from_64(0x0f0f0e0d0b0b0a09LL),
                                          v64_from_64(0x0707060503030201LL)));
      e = v128_shuffle_8(e, v128_from_v64(v64_from_64(0x0f0f0f0e0b0b0b0aLL),
                                          v64_from_64(0x0707070603030302LL)));
    }

    // delta = 4/16 * clamp(a - x, -s, s) + 1/16 * clamp(b - x, -s, s) +
    //         3/16 * clamp(c - x, -s, s) + 3/16 * clamp(d - x, -s, s) +
    //         1/16 * clamp(e - x, -s, s) + 4/16 * clamp(f - x, -s, s)
    // x is o with an offset of 128 so 8 bit saturated sub can be used.
    tmp = v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                     v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
    delta = v128_add_8(
        v128_add_8(
            v128_shl_8(
                v128_add_8(
                    v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                    v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
                2),
            v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                       v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
        v128_add_8(v128_add_8(tmp, tmp), tmp));
    delta = v128_shr_s8(
        v128_add_8(v128_dup_8(8),
                   v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
        4);
    r = v128_add_8(o, delta);
    u32_store_aligned(dst, v128_low_u32(v128_shr_n_byte(r, 12)));
    u32_store_aligned(dst + dstride, v128_low_u32(v128_shr_n_byte(r, 8)));
    u32_store_aligned(dst + 2 * dstride, v128_low_u32(v128_shr_n_byte(r, 4)));
    u32_store_aligned(dst + 3 * dstride, v128_low_u32(r));

    dst += 4 * dstride;
    src += 4 * sstride;
  }
}

// Process 2 lines at a time, 8 bit.
static void clpf_block(const uint8_t *src, uint8_t *dst, int sstride,
                       int dstride, int x0, int y0, int sizey, int width,
                       int height, unsigned int strength) {
  int bottom = height - 2 - y0;
  const v128 sp = v128_dup_8(strength);
  const v128 sm = v128_dup_8(-(int)strength);
  // The difference will be 9 bit, offset by 128 so we can use saturated
  // sub to avoid going to 16 bit temporarily before "strength" clipping.
  const v128 c128 = v128_dup_8(128);
  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  if (!x0) {  // Clip left
    const v128 b_shuff = v128_from_v64(v64_from_64(0x0d0c0b0a09080808LL),
                                       v64_from_64(0x0504030201000000LL));
    const v128 c_shuff = v128_from_v64(v64_from_64(0x0e0d0c0b0a090808LL),
                                       v64_from_64(0x0605040302010000LL));
    int y;

    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_shuffle_8(x, b_shuff);
      const v128 c = v128_shuffle_8(x, c_shuff);
      const v128 d = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 1),
                              v64_load_unaligned(src + 1 + sstride)));
      const v128 e = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 2),
                              v64_load_unaligned(src + 2 + sstride)));
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    const v128 d_shuff = v128_from_v64(v64_from_64(0x0f0f0e0d0c0b0a09LL),
                                       v64_from_64(0x0707060504030201LL));
    const v128 e_shuff = v128_from_v64(v64_from_64(0x0f0f0f0e0d0c0b0aLL),
                                       v64_from_64(0x0707070605040302LL));
    int y;

    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 2),
                              v64_load_unaligned(src - 2 + sstride)));
      const v128 c = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 1),
                              v64_load_unaligned(src - 1 + sstride)));
      const v128 d = v128_shuffle_8(x, d_shuff);
      const v128 e = v128_shuffle_8(x, e_shuff);
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else {  // No left/right clipping
    int y;
    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 x = v128_add_8(c128, o);
      const v128 a = v128_add_8(
          c128,
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1));
      const v128 b = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 2),
                              v64_load_unaligned(src - 2 + sstride)));
      const v128 c = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src - 1),
                              v64_load_unaligned(src - 1 + sstride)));
      const v128 d = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 1),
                              v64_load_unaligned(src + 1 + sstride)));
      const v128 e = v128_add_8(
          c128, v128_from_v64(v64_load_unaligned(src + 2),
                              v64_load_unaligned(src + 2 + sstride)));
      const v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride)));
      calc_delta(o, x, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  }
}

// SIMD implementation for aom_clpf_block_c()
void SIMD_FUNC(aom_clpf_block)(const uint8_t *src, uint8_t *dst, int sstride,
                               int dstride, int x0, int y0, int sizex,
                               int sizey, int width, int height,
                               unsigned int strength) {
  if ((sizex != 4 && sizex != 8) || width < 16 || y0 + 4 > height ||
      x0 + 4 > width) {
    // Fallback to C for odd sizes
    aom_clpf_block_c(src, dst, sstride, dstride, x0, y0, sizex, sizey, width,
                     height, strength);
  } else {
    (sizex == 4 ? clpf_block4 : clpf_block)(src, dst, sstride, dstride, x0, y0,
                                            sizey, width, height, strength);
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
// delta = 4/16 * clamp(a - o, -s, s) + 1/16 * clamp(b - o, -s, s) +
//         3/16 * clamp(c - o, -s, s) + 3/16 * clamp(d - o, -s, s) +
//         1/16 * clamp(e - o, -s, s) + 4/16 * clamp(f - o, -s, s)
static void calc_delta_hbd4(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                            v128 f, uint16_t *dst, v128 sp, v128 sm,
                            int dstride) {
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
  o = v128_add_16(
      o, v128_shr_s16(
             v128_add_16(
                 c8, v128_add_16(delta, v128_cmplt_s16(delta, v128_zero()))),
             4));
  v64_store_aligned(dst, v128_high_v64(o));
  v64_store_aligned(dst + dstride, v128_low_v64(o));
}

// As above but the result is stored in just one line.
static void calc_delta_hbd(v128 o, v128 a, v128 b, v128 c, v128 d, v128 e,
                           v128 f, uint16_t *dst, v128 sp, v128 sm) {
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
  v128_store_aligned(
      dst,
      v128_add_16(
          o, v128_shr_s16(
                 v128_add_16(c8, v128_add_16(delta, v128_cmplt_s16(
                                                        delta, v128_zero()))),
                 4)));
}

// Modelled like clpf_block, processing two lines at a time, 16 bit.
SIMD_INLINE void clpf_block_hbd4(const uint16_t *src, uint16_t *dst,
                                 int sstride, int dstride, int x0, int y0,
                                 int sizey, int width, int height,
                                 unsigned int strength) {
  int bottom = height - 2 - y0;
  const v128 sp = v128_dup_16(strength);
  const v128 sm = v128_dup_16(-(int)strength);
  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  if (!x0) {  // Clip left
    const v128 b_shuff = v128_from_v64(v64_from_64(0x0b0a090809080908LL),
                                       v64_from_64(0x0302010001000100LL));
    const v128 c_shuff = v128_from_v64(v64_from_64(0x0d0c0b0a09080908LL),
                                       v64_from_64(0x0504030201000100LL));
    int y;

    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 a =
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1);
      const v128 b = v128_shuffle_8(o, b_shuff);
      const v128 c = v128_shuffle_8(o, c_shuff);
      const v128 d = v128_from_v64(v64_load_unaligned(src + 1),
                                   v64_load_unaligned(src + 1 + sstride));
      const v128 e = v128_from_v64(v64_load_unaligned(src + 2),
                                   v64_load_unaligned(src + 2 + sstride));
      const v128 f = v128_from_v64(
          l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride));
      calc_delta_hbd4(o, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else if (!(width - x0 - 4)) {  // Clip right
    const v128 d_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0d0c0b0aLL),
                                       v64_from_64(0x0706070605040302LL));
    const v128 e_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0f0e0d0cLL),
                                       v64_from_64(0x0706070607060504LL));
    int y;
    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 a =
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1);
      const v128 b = v128_from_v64(v64_load_unaligned(src - 2),
                                   v64_load_unaligned(src - 2 + sstride));
      const v128 c = v128_from_v64(v64_load_unaligned(src - 1),
                                   v64_load_unaligned(src - 1 + sstride));
      const v128 d = v128_shuffle_8(o, d_shuff);
      const v128 e = v128_shuffle_8(o, e_shuff);
      const v128 f = v128_from_v64(
          l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride));
      calc_delta_hbd4(o, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  } else {  // No left/right clipping
    int y;
    for (y = 0; y < sizey; y += 2) {
      const v64 l1 = v64_load_aligned(src);
      const v64 l2 = v64_load_aligned(src + sstride);
      v128 o = v128_from_v64(l1, l2);
      const v128 a =
          v128_from_v64(v64_load_aligned(src - (y != -y0) * sstride), l1);
      const v128 b = v128_from_v64(v64_load_unaligned(src - 2),
                                   v64_load_unaligned(src - 2 + sstride));
      const v128 c = v128_from_v64(v64_load_unaligned(src - 1),
                                   v64_load_unaligned(src - 1 + sstride));
      const v128 d = v128_from_v64(v64_load_unaligned(src + 1),
                                   v64_load_unaligned(src + 1 + sstride));
      const v128 e = v128_from_v64(v64_load_unaligned(src + 2),
                                   v64_load_unaligned(src + 2 + sstride));
      const v128 f = v128_from_v64(
          l2, v64_load_aligned(src + ((y != bottom) + 1) * sstride));
      calc_delta_hbd4(o, a, b, c, d, e, f, dst, sp, sm, dstride);
      src += sstride * 2;
      dst += dstride * 2;
    }
  }
}

// The most simple case.  Start here if you need to understand the functions.
SIMD_INLINE void clpf_block_hbd(const uint16_t *src, uint16_t *dst, int sstride,
                                int dstride, int x0, int y0, int sizey,
                                int width, int height, unsigned int strength) {
  int y;
  int bottom = height - 2 - y0;
  const v128 sp = v128_dup_16(strength);
  const v128 sm = v128_dup_16(-(int)strength);

  dst += x0 + y0 * dstride;
  src += x0 + y0 * sstride;

  // Read 8 set of pixels at a time.  Clipping along upper and lower
  // edges is handled by reading the upper or lower line twice.
  // Clipping along the left and right edges is handled by shuffle
  // instructions doing shift and pad.
  if (!x0) {  // Clip left
    const v128 b_shuff = v128_from_v64(v64_from_64(0x0b0a090807060504LL),
                                       v64_from_64(0x0302010001000100LL));
    const v128 c_shuff = v128_from_v64(v64_from_64(0x0d0c0b0a09080706LL),
                                       v64_from_64(0x0504030201000100LL));
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_shuffle_8(o, b_shuff);
      const v128 c = v128_shuffle_8(o, c_shuff);
      const v128 d = v128_load_unaligned(src + 1);
      const v128 e = v128_load_unaligned(src + 2);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    const v128 d_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0d0c0b0aLL),
                                       v64_from_64(0x0908070605040302LL));
    const v128 e_shuff = v128_from_v64(v64_from_64(0x0f0e0f0e0f0e0d0cLL),
                                       v64_from_64(0x0b0a090807060504LL));
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_load_unaligned(src - 2);
      const v128 c = v128_load_unaligned(src - 1);
      const v128 d = v128_shuffle_8(o, d_shuff);
      const v128 e = v128_shuffle_8(o, e_shuff);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  } else {  // No left/right clipping
    for (y = 0; y < sizey; y++) {
      const v128 o = v128_load_aligned(src);
      const v128 a = v128_load_aligned(src - (y != -y0) * sstride);
      const v128 b = v128_load_unaligned(src - 2);
      const v128 c = v128_load_unaligned(src - 1);
      const v128 d = v128_load_unaligned(src + 1);
      const v128 e = v128_load_unaligned(src + 2);
      const v128 f = v128_load_aligned(src + (y - 1 != bottom) * sstride);
      calc_delta_hbd(o, a, b, c, d, e, f, dst, sp, sm);
      src += sstride;
      dst += dstride;
    }
  }
}

// SIMD implementation for aom_clpf_block_hbd_c()
void SIMD_FUNC(aom_clpf_block_hbd)(const uint16_t *src, uint16_t *dst,
                                   int sstride, int dstride, int x0, int y0,
                                   int sizex, int sizey, int width, int height,
                                   unsigned int strength) {
  if ((sizex != 4 && sizex != 8) || width < 16 || y0 + 4 > height ||
      x0 + 4 > width) {
    // Fallback to C for odd sizes
    aom_clpf_block_hbd_c(src, dst, sstride, dstride, x0, y0, sizex, sizey,
                         width, height, strength);
  } else {
    (sizex == 4 ? clpf_block_hbd4 : clpf_block_hbd)(
        src, dst, sstride, dstride, x0, y0, sizey, width, height, strength);
  }
}
#endif
