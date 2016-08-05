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

#include "aom_dsp/aom_simd.h"

void SIMD_FUNC(aom_detect_clpf)(const uint8_t *rec, const uint8_t *org, int x0,
                                int y0, int width, int height, int so,
                                int stride, int *sum0, int *sum1,
                                unsigned int strength) {
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  v128 c128 = v128_dup_8(128);
  v128 sp = v128_dup_8(strength);
  v128 sm = v128_dup_8(-(int)strength);
  int bottom = height - 2 - y0;

  rec += x0 + y0 * stride;
  org += x0 + y0 * so;

  if (!x0) {  // Clip left
    v128 s1 = v128_from_v64(v64_from_64(0x0e0d0c0b0a090808LL),
                            v64_from_64(0x0605040302010000LL));
    v128 s2 = v128_from_v64(v64_from_64(0x0d0c0b0a09080808LL),
                            v64_from_64(0x0504030201000000LL));

    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b = v128_add_8(c128, v128_shuffle_8(q, s2));
      v128 c = v128_add_8(c128, v128_shuffle_8(q, s1));
      v128 d =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 1),
                                         v64_load_unaligned(rec + 1 + stride)));
      v128 e =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 2),
                                         v64_load_unaligned(rec + 2 + stride)));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));

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
      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(ssd1, o, v128_add_8(q, delta));
      rec += stride * 2;
      org += so * 2;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    v128 s1 = v128_from_v64(v64_from_64(0x0f0f0e0d0c0b0a09LL),
                            v64_from_64(0x0707060504030201LL));
    v128 s2 = v128_from_v64(v64_from_64(0x0f0f0f0e0d0c0b0aLL),
                            v64_from_64(0x0707070605040302LL));

    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 2),
                                         v64_load_unaligned(rec - 2 + stride)));
      v128 c =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 1),
                                         v64_load_unaligned(rec - 1 + stride)));
      v128 d = v128_add_8(c128, v128_shuffle_8(q, s1));
      v128 e = v128_add_8(c128, v128_shuffle_8(q, s2));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));

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
      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(ssd1, o, v128_add_8(q, delta));
      rec += stride * 2;
      org += so * 2;
    }
  } else {  // No left/right clipping
    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 2),
                                         v64_load_unaligned(rec - 2 + stride)));
      v128 c =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 1),
                                         v64_load_unaligned(rec - 1 + stride)));
      v128 d =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 1),
                                         v64_load_unaligned(rec + 1 + stride)));
      v128 e =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 2),
                                         v64_load_unaligned(rec + 2 + stride)));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));

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

      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(ssd1, o, v128_add_8(q, delta));
      rec += stride * 2;
      org += so * 2;
    }
  }
  *sum0 += v128_ssd_u8_sum(ssd0);
  *sum1 += v128_ssd_u8_sum(ssd1);
}

// Test multiple filter strengths at once.  Use a simpler filter (4 tap, every
// second line).
void SIMD_FUNC(aom_detect_multi_clpf)(const uint8_t *rec, const uint8_t *org,
                                      int x0, int y0, int width, int height,
                                      int so, int stride, int *sum) {
  v128 c128 = v128_dup_8(128);
  v128 cp1 = v128_dup_8(1);
  v128 cm1 = v128_dup_8(-1);
  v128 cp2 = v128_dup_8(2);
  v128 cm2 = v128_dup_8(-2);
  v128 cp4 = v128_dup_8(4);
  v128 cm4 = v128_dup_8(-4);
  v128 c8 = v128_dup_8(8);
  int bottom = height - 2 - y0;
  ssd128_internal ssd0 = v128_ssd_u8_init();
  ssd128_internal ssd1 = v128_ssd_u8_init();
  ssd128_internal ssd2 = v128_ssd_u8_init();
  ssd128_internal ssd3 = v128_ssd_u8_init();

  rec += x0 + y0 * stride;
  org += x0 + y0 * so;

  if (!x0) {  // Clip left
    v128 s1 = v128_from_v64(v64_from_64(0x0e0d0c0b0a090808LL),
                            v64_from_64(0x0605040302010000LL));
    v128 s2 = v128_from_v64(v64_from_64(0x0d0c0b0a09080808LL),
                            v64_from_64(0x0504030201000000LL));

    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b = v128_add_8(c128, v128_shuffle_8(q, s2));
      v128 c = v128_add_8(c128, v128_shuffle_8(q, s1));
      v128 d =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 1),
                                         v64_load_unaligned(rec + 1 + stride)));
      v128 e =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 2),
                                         v64_load_unaligned(rec + 2 + stride)));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));
      v128 tmp, delta1, delta2, delta3;

      a = v128_ssub_s8(a, x);
      b = v128_ssub_s8(b, x);
      c = v128_ssub_s8(c, x);
      d = v128_ssub_s8(d, x);
      e = v128_ssub_s8(e, x);
      f = v128_ssub_s8(f, x);
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp1), cm1),
                       v128_max_s8(v128_min_s8(d, cp1), cm1));
      delta1 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp1), cm1),
                                    v128_max_s8(v128_min_s8(f, cp1), cm1)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp1), cm1),
                         v128_max_s8(v128_min_s8(e, cp1), cm1))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp2), cm2),
                       v128_max_s8(v128_min_s8(d, cp2), cm2));
      delta2 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp2), cm2),
                                    v128_max_s8(v128_min_s8(f, cp2), cm2)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp2), cm2),
                         v128_max_s8(v128_min_s8(e, cp2), cm2))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp4), cm4),
                       v128_max_s8(v128_min_s8(d, cp4), cm4));
      delta3 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp4), cm4),
                                    v128_max_s8(v128_min_s8(f, cp4), cm4)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp4), cm4),
                         v128_max_s8(v128_min_s8(e, cp4), cm4))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));

      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(
          ssd1, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta1, v128_cmplt_s8(
                                                        delta1, v128_zero()))),
                  4)));
      ssd2 = v128_ssd_u8(
          ssd2, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta2, v128_cmplt_s8(
                                                        delta2, v128_zero()))),
                  4)));
      ssd3 = v128_ssd_u8(
          ssd3, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta3, v128_cmplt_s8(
                                                        delta3, v128_zero()))),
                  4)));
      rec += 2 * stride;
      org += 2 * so;
    }
  } else if (!(width - x0 - 8)) {  // Clip right
    v128 s1 = v128_from_v64(v64_from_64(0x0f0f0e0d0c0b0a09LL),
                            v64_from_64(0x0707060504030201LL));
    v128 s2 = v128_from_v64(v64_from_64(0x0f0f0f0e0d0c0b0aLL),
                            v64_from_64(0x0707070605040302LL));

    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 2),
                                         v64_load_unaligned(rec - 2 + stride)));
      v128 c =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 1),
                                         v64_load_unaligned(rec - 1 + stride)));
      v128 d = v128_add_8(c128, v128_shuffle_8(q, s1));
      v128 e = v128_add_8(c128, v128_shuffle_8(q, s2));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));
      v128 tmp, delta1, delta2, delta3;

      a = v128_ssub_s8(a, x);
      b = v128_ssub_s8(b, x);
      c = v128_ssub_s8(c, x);
      d = v128_ssub_s8(d, x);
      e = v128_ssub_s8(e, x);
      f = v128_ssub_s8(f, x);
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp1), cm1),
                       v128_max_s8(v128_min_s8(d, cp1), cm1));
      delta1 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp1), cm1),
                                    v128_max_s8(v128_min_s8(f, cp1), cm1)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp1), cm1),
                         v128_max_s8(v128_min_s8(e, cp1), cm1))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp2), cm2),
                       v128_max_s8(v128_min_s8(d, cp2), cm2));
      delta2 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp2), cm2),
                                    v128_max_s8(v128_min_s8(f, cp2), cm2)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp2), cm2),
                         v128_max_s8(v128_min_s8(e, cp2), cm2))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp4), cm4),
                       v128_max_s8(v128_min_s8(d, cp4), cm4));
      delta3 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp4), cm4),
                                    v128_max_s8(v128_min_s8(f, cp4), cm4)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp4), cm4),
                         v128_max_s8(v128_min_s8(e, cp4), cm4))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));

      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(
          ssd1, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta1, v128_cmplt_s8(
                                                        delta1, v128_zero()))),
                  4)));
      ssd2 = v128_ssd_u8(
          ssd2, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta2, v128_cmplt_s8(
                                                        delta2, v128_zero()))),
                  4)));
      ssd3 = v128_ssd_u8(
          ssd3, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta3, v128_cmplt_s8(
                                                        delta3, v128_zero()))),
                  4)));
      rec += 2 * stride;
      org += 2 * so;
    }
  } else {  // No left/right clipping
    for (int y = 0; y < 8; y += 2) {
      v64 k1 = v64_load_aligned(org);
      v64 k2 = v64_load_aligned(org + so);
      v64 l1 = v64_load_aligned(rec);
      v64 l2 = v64_load_aligned(rec + stride);
      v128 o = v128_from_v64(k1, k2);
      v128 q = v128_from_v64(l1, l2);
      v128 x = v128_add_8(c128, q);
      v128 a = v128_add_8(
          c128, v128_from_v64(v64_load_aligned(rec - (y != -y0) * stride), l1));
      v128 b =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 2),
                                         v64_load_unaligned(rec - 2 + stride)));
      v128 c =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec - 1),
                                         v64_load_unaligned(rec - 1 + stride)));
      v128 d =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 1),
                                         v64_load_unaligned(rec + 1 + stride)));
      v128 e =
          v128_add_8(c128, v128_from_v64(v64_load_unaligned(rec + 2),
                                         v64_load_unaligned(rec + 2 + stride)));
      v128 f = v128_add_8(
          c128, v128_from_v64(
                    l2, v64_load_aligned(rec + ((y != bottom) + 1) * stride)));
      v128 tmp, delta1, delta2, delta3;

      a = v128_ssub_s8(a, x);
      b = v128_ssub_s8(b, x);
      c = v128_ssub_s8(c, x);
      d = v128_ssub_s8(d, x);
      e = v128_ssub_s8(e, x);
      f = v128_ssub_s8(f, x);
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp1), cm1),
                       v128_max_s8(v128_min_s8(d, cp1), cm1));
      delta1 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp1), cm1),
                                    v128_max_s8(v128_min_s8(f, cp1), cm1)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp1), cm1),
                         v128_max_s8(v128_min_s8(e, cp1), cm1))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp2), cm2),
                       v128_max_s8(v128_min_s8(d, cp2), cm2));
      delta2 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp2), cm2),
                                    v128_max_s8(v128_min_s8(f, cp2), cm2)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp2), cm2),
                         v128_max_s8(v128_min_s8(e, cp2), cm2))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));
      tmp = v128_add_8(v128_max_s8(v128_min_s8(c, cp4), cm4),
                       v128_max_s8(v128_min_s8(d, cp4), cm4));
      delta3 = v128_add_8(
          v128_add_8(
              v128_shl_8(v128_add_8(v128_max_s8(v128_min_s8(a, cp4), cm4),
                                    v128_max_s8(v128_min_s8(f, cp4), cm4)),
                         2),
              v128_add_8(v128_max_s8(v128_min_s8(b, cp4), cm4),
                         v128_max_s8(v128_min_s8(e, cp4), cm4))),
          v128_add_8(v128_add_8(tmp, tmp), tmp));

      ssd0 = v128_ssd_u8(ssd0, o, q);
      ssd1 = v128_ssd_u8(
          ssd1, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta1, v128_cmplt_s8(
                                                        delta1, v128_zero()))),
                  4)));
      ssd2 = v128_ssd_u8(
          ssd2, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta2, v128_cmplt_s8(
                                                        delta2, v128_zero()))),
                  4)));
      ssd3 = v128_ssd_u8(
          ssd3, o,
          v128_add_8(
              q,
              v128_shr_s8(
                  v128_add_8(c8, v128_add_8(delta3, v128_cmplt_s8(
                                                        delta3, v128_zero()))),
                  4)));
      rec += 2 * stride;
      org += 2 * so;
    }
  }
  sum[0] += v128_ssd_u8_sum(ssd0);
  sum[1] += v128_ssd_u8_sum(ssd1);
  sum[2] += v128_ssd_u8_sum(ssd2);
  sum[3] += v128_ssd_u8_sum(ssd3);
}
