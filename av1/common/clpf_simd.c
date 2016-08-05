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

static void clpf_block4(const uint8_t *src, uint8_t *dst, int stride, int x0,
                        int y0, int width, int height, unsigned int strength) {
  dst += x0 + y0 * stride;
  src += x0 + y0 * stride;

  {
    int right = width - x0 - 4;
    int bottom = height - y0 - 4;
    uint32_t l0 = *(uint32_t *)(src - !!y0 * stride);
    uint32_t l1 = *(uint32_t *)(src);
    uint32_t l2 = *(uint32_t *)(src + stride);
    uint32_t l3 = *(uint32_t *)(src + 2 * stride);
    uint32_t l4 = *(uint32_t *)(src + 3 * stride);
    uint32_t l5 = *(uint32_t *)(src + (3 + !!bottom) * stride);
    v128 c128 = v128_dup_8(128);
    v128 o = v128_from_32(l1, l2, l3, l4);
    v128 x = v128_add_8(c128, o);
    v128 a = v128_add_8(c128, v128_from_32(l0, l1, l2, l3));
    v128 b = v128_add_8(
        c128, v128_from_32(*(uint32_t *)(src - 2 * !!x0),
                           *(uint32_t *)(src + stride - 2 * !!x0),
                           *(uint32_t *)(src + 2 * stride - 2 * !!x0),
                           *(uint32_t *)(src + 3 * stride - 2 * !!x0)));
    v128 c =
        v128_add_8(c128, v128_from_32(*(uint32_t *)(src - !!x0),
                                      *(uint32_t *)(src + stride - !!x0),
                                      *(uint32_t *)(src + 2 * stride - !!x0),
                                      *(uint32_t *)(src + 3 * stride - !!x0)));
    v128 d = v128_add_8(
        c128, v128_from_32(*(uint32_t *)(src + !!right),
                           *(uint32_t *)(src + stride + !!right),
                           *(uint32_t *)(src + 2 * stride + !!right),
                           *(uint32_t *)(src + 3 * stride + !!right)));
    v128 e = v128_add_8(
        c128, v128_from_32(*(uint32_t *)(src + 2 * !!right),
                           *(uint32_t *)(src + stride + 2 * !!right),
                           *(uint32_t *)(src + 2 * stride + 2 * !!right),
                           *(uint32_t *)(src + 3 * stride + 2 * !!right)));
    v128 f = v128_add_8(c128, v128_from_32(l2, l3, l4, l5));
    if (!x0) {
      b = v128_shuffle_8(c, v128_from_v64(v64_from_64(0x0d0c0c0c09080808LL),
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
    {
      v128 r;
      v128 sp = v128_dup_8(strength);
      v128 sm = v128_dup_8(-(int)strength);

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
              v128_add_8(v128_max_s8(v128_min_s8(v128_ssub_s8(b, x), sp), sm),
                         v128_max_s8(v128_min_s8(v128_ssub_s8(e, x), sp), sm))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      delta = v128_shr_s8(
          v128_add_8(v128_dup_8(8),
                     v128_add_8(delta, v128_cmplt_s8(delta, v128_zero()))),
          4);
      r = v128_add_8(o, delta);
      *(uint32_t *)dst = v128_low_u32(v128_shr_n_byte(r, 12));
      *(uint32_t *)(dst + stride) = v128_low_u32(v128_shr_n_byte(r, 8));
      *(uint32_t *)(dst + 2 * stride) = v128_low_u32(v128_shr_n_byte(r, 4));
      *(uint32_t *)(dst + 3 * stride) = v128_low_u32(r);
    }
  }
}

static void clpf_block8(const uint8_t *src, uint8_t *dst, int stride, int x0,
                        int y0, int width, int height, unsigned int strength) {
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

      for (int y = 0; y < 8; y += 2) {
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

      for (int y = 0; y < 8; y += 2) {
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
      for (int y = 0; y < 8; y += 2) {
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
  if (sizex == 4 && sizey == 4)  // chroma 420
    clpf_block4(src, dst, stride, x0, y0, width, height, strength);
  else if (sizex == 4) {  // chroma 422
    clpf_block4(src, dst, stride, x0, y0, width, height, strength);
    clpf_block4(src + 4 * stride, dst + 4 * stride, stride, x0, y0, width,
                height, strength);
  } else  // luma 444 or luma 420
    clpf_block8(src, dst, stride, x0, y0, width, height, strength);
}
