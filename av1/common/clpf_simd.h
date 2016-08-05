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

static void clpf_block(const uint8_t *src, uint8_t *dst, int stride, int x0,
                       int y0, int width, int height, int sizey,
                       unsigned int strength) {
  dst += x0 + y0 * stride;
  src += x0 + y0 * stride;
  {
    int bottom = height - 2 - y0;
    v128 sp = v128_dup_8(strength);
    v128 sm = v128_dup_8(-(int)strength);
    v128 c8 = v128_dup_8(8);
    v128 c128 = v128_dup_8(128);

    if (!x0) {  // Clip left
      v128 s1 = v128_from_v64(v64_from_64(0x0e0d0c0b0a090808LL),
                              v64_from_64(0x0605040302010000LL));
      v128 s2 = v128_from_v64(v64_from_64(0x0d0c0b0a09080808LL),
                              v64_from_64(0x0504030201000000LL));
      int y;

      for (y = 0; y < sizey; y += 2) {
        v64 l1 = v64_load_aligned(src);
        v64 l2 = v64_load_aligned(src + stride);
        v128 o = v128_from_v64(l1, l2);
        v128 x = v128_add_8(c128, o);
        v128 a = v128_add_8(
            c128,
            v128_from_v64(v64_load_aligned(src - (y != -y0) * stride), l1));
        v128 b = v128_add_8(c128, v128_shuffle_8(o, s2));
        v128 c = v128_add_8(c128, v128_shuffle_8(o, s1));
        v128 d = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src + 1),
                                v64_load_unaligned(src + 1 + stride)));
        v128 e = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src + 2),
                                v64_load_unaligned(src + 2 + stride)));
        v128 f = v128_add_8(
            c128, v128_from_v64(l2, v64_load_aligned(
                                        src + ((y != bottom) + 1) * stride)));

        v128 tmp =
            v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                       v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
        v128 delta = v128_add_8(
            v128_add_8(
                v128_shl_8(
                    v128_add_8(
                        v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                        v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
                    2),
                v128_add_8(
                    v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                    v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
            v128_add_8(v128_add_8(tmp, tmp), tmp));
        o = v128_add_8(
            o, v128_shr_s8(
                   v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(
                                                        delta, v128_zero()))),
                   4));
        v64_store_aligned(dst, v128_high_v64(o));
        v64_store_aligned(dst + stride, v128_low_v64(o));
        src += stride * 2;
        dst += stride * 2;
      }
    } else if (!(width - x0 - 8)) {  // Clip right
      v128 s1 = v128_from_v64(v64_from_64(0x0f0f0e0d0c0b0a09LL),
                              v64_from_64(0x0707060504030201LL));
      v128 s2 = v128_from_v64(v64_from_64(0x0f0f0f0e0d0c0b0aLL),
                              v64_from_64(0x0707070605040302LL));
      int y;

      for (y = 0; y < sizey; y += 2) {
        v64 l1 = v64_load_aligned(src);
        v64 l2 = v64_load_aligned(src + stride);
        v128 o = v128_from_v64(l1, l2);
        v128 x = v128_add_8(c128, o);
        v128 a = v128_add_8(
            c128,
            v128_from_v64(v64_load_aligned(src - (y != -y0) * stride), l1));
        v128 b = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src - 2),
                                v64_load_unaligned(src - 2 + stride)));
        v128 c = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src - 1),
                                v64_load_unaligned(src - 1 + stride)));
        v128 d = v128_add_8(c128, v128_shuffle_8(o, s1));
        v128 e = v128_add_8(c128, v128_shuffle_8(o, s2));
        v128 f = v128_add_8(
            c128, v128_from_v64(l2, v64_load_aligned(
                                        src + ((y != bottom) + 1) * stride)));

        v128 tmp =
            v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                       v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
        v128 delta = v128_add_8(
            v128_add_8(
                v128_shl_8(
                    v128_add_8(
                        v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                        v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
                    2),
                v128_add_8(
                    v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                    v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
            v128_add_8(v128_add_8(tmp, tmp), tmp));
        o = v128_add_8(
            o, v128_shr_s8(
                   v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(
                                                        delta, v128_zero()))),
                   4));
        v64_store_aligned(dst, v128_high_v64(o));
        v64_store_aligned(dst + stride, v128_low_v64(o));
        src += stride * 2;
        dst += stride * 2;
      }
    } else {  // No left/right clipping
      int y;
      for (y = 0; y < sizey; y += 2) {
        v64 l1 = v64_load_aligned(src);
        v64 l2 = v64_load_aligned(src + stride);
        v128 o = v128_from_v64(l1, l2);
        v128 x = v128_add_8(c128, o);
        v128 a = v128_add_8(
            c128,
            v128_from_v64(v64_load_aligned(src - (y != -y0) * stride), l1));
        v128 b = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src - 2),
                                v64_load_unaligned(src - 2 + stride)));
        v128 c = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src - 1),
                                v64_load_unaligned(src - 1 + stride)));
        v128 d = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src + 1),
                                v64_load_unaligned(src + 1 + stride)));
        v128 e = v128_add_8(
            c128, v128_from_v64(v64_load_unaligned(src + 2),
                                v64_load_unaligned(src + 2 + stride)));
        v128 f = v128_add_8(
            c128, v128_from_v64(l2, v64_load_aligned(
                                        src + ((y != bottom) + 1) * stride)));

        v128 tmp =
            v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(c, x), sp), sm),
                       v128_max_s8(v128_min_s8(v128_ssub_s8(d, x), sp), sm));
        v128 delta = v128_add_8(
            v128_add_8(
                v128_shl_8(
                    v128_add_8(
                        v128_max_s8(v128_min_s8(v128_ssub_s8(a, x), sp), sm),
                        v128_max_s8(v128_min_s8(v128_ssub_s8(f, x), sp), sm)),
                    2),
                v128_add_8(
                    v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                    v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
            v128_add_8(v128_add_8(tmp, tmp), tmp));
        o = v128_add_8(
            o, v128_shr_s8(
                   v128_add_8(c8, v128_add_8(delta, v128_cmplt_s8(
                                                        delta, v128_zero()))),
                   4));
        v64_store_aligned(dst, v128_high_v64(o));
        v64_store_aligned(dst + stride, v128_low_v64(o));
        src += stride * 2;
        dst += stride * 2;
      }
    }
  }
}

void SIMD_FUNC(aom_clpf_block)(const uint8_t *src, uint8_t *dst, int stride,
                               int x0, int y0, int sizex, int sizey, int width,
                               int height, unsigned int strength) {
  // TODO(stemidts):
  // A sizex different from 8 will only be needed if CLPF is extended to chroma.
  // This will only be used if 4:2:0 and width not a multiple of 16 and along
  // the right edge only, so we can fall back to the plain C implementation in
  // this case.  If not extended to chroma, this test will be redundant.
  if (sizex != 8) {
    aom_clpf_block_c(src, dst, stride, x0, y0, sizex, sizey, width, height,
                     strength);
  } else {
    clpf_block(src, dst, stride, x0, y0, width, height, sizey, strength);
  }
}
