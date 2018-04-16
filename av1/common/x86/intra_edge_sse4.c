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

#include <assert.h>
#include <smmintrin.h>

#include "./aom_config.h"
#include "./av1_rtcd.h"

void av1_filter_intra_edge_sse4_1(uint8_t *p, int sz, int strength) {
  if (!strength) return;

  DECLARE_ALIGNED(16, static const int8_t, kern[3][16]) = {
    { 4, 8, 4, 0, 4, 8, 4, 0, 4, 8, 4, 0, 4, 8, 4, 0 },  // strength 1: 4,8,4
    { 5, 6, 5, 0, 5, 6, 5, 0, 5, 6, 5, 0, 5, 6, 5, 0 },  // strength 2: 5,6,5
    { 2, 4, 4, 4, 2, 0, 0, 0, 2, 4, 4, 4, 2, 0, 0, 0 }  // strength 3: 2,4,4,4,2
  };

  DECLARE_ALIGNED(16, static const int8_t, v_const[5][16]) = {
    { 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6 },
    { 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  };

  // Extend the first and last samples to simplify the loop for the 5-tap case
  p[-1] = p[0];
  __m128i last = _mm_set1_epi8(p[sz - 1]);
  _mm_storeu_si128((__m128i *)&p[sz], last);

  // Adjust input pointer for filter support area
  uint8_t *in = (strength == 3) ? p - 1 : p;

  // Avoid modifying first sample
  uint8_t *out = p + 1;
  int len = sz - 1;

  const int use_3tap_filter = (strength < 3);

  if (use_3tap_filter) {
    __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
    __m128i shuf0 = _mm_lddqu_si128((__m128i const *)v_const[0]);
    __m128i shuf1 = _mm_lddqu_si128((__m128i const *)v_const[1]);
    __m128i iden = _mm_lddqu_si128((__m128i *)v_const[3]);
    __m128i in0 = _mm_lddqu_si128((__m128i *)in);
    while (len > 0) {
      int n_out = (len < 8) ? len : 8;
      __m128i d0 = _mm_shuffle_epi8(in0, shuf0);
      __m128i d1 = _mm_shuffle_epi8(in0, shuf1);
      d0 = _mm_maddubs_epi16(d0, coef0);
      d1 = _mm_maddubs_epi16(d1, coef0);
      d0 = _mm_hadd_epi16(d0, d1);
      __m128i eight = _mm_set1_epi16(8);
      d0 = _mm_add_epi16(d0, eight);
      d0 = _mm_srai_epi16(d0, 4);
      d0 = _mm_packus_epi16(d0, d0);
      __m128i out0 = _mm_lddqu_si128((__m128i *)out);
      __m128i n0 = _mm_set1_epi8(n_out);
      __m128i mask = _mm_cmpgt_epi8(n0, iden);
      out0 = _mm_blendv_epi8(out0, d0, mask);
      _mm_storel_epi64((__m128i *)out, out0);
      __m128i in1 = _mm_lddqu_si128((__m128i *)(in + 16));
      in0 = _mm_alignr_epi8(in1, in0, 8);
      in += 8;
      out += 8;
      len -= n_out;
    }
  } else {  // 5-tap filter
    __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
    __m128i two = _mm_set1_epi8(2);
    __m128i shuf_a = _mm_lddqu_si128((__m128i const *)v_const[2]);
    __m128i shuf_b = _mm_add_epi8(shuf_a, two);
    __m128i shuf_c = _mm_add_epi8(shuf_b, two);
    __m128i shuf_d = _mm_add_epi8(shuf_c, two);
    __m128i iden = _mm_lddqu_si128((__m128i *)v_const[3]);
    __m128i in0 = _mm_lddqu_si128((__m128i *)in);
    while (len > 0) {
      int n_out = (len < 8) ? len : 8;
      __m128i d0 = _mm_shuffle_epi8(in0, shuf_a);
      __m128i d1 = _mm_shuffle_epi8(in0, shuf_b);
      __m128i d2 = _mm_shuffle_epi8(in0, shuf_c);
      __m128i d3 = _mm_shuffle_epi8(in0, shuf_d);
      d0 = _mm_maddubs_epi16(d0, coef0);
      d1 = _mm_maddubs_epi16(d1, coef0);
      d2 = _mm_maddubs_epi16(d2, coef0);
      d3 = _mm_maddubs_epi16(d3, coef0);
      d0 = _mm_hadd_epi16(d0, d1);
      d2 = _mm_hadd_epi16(d2, d3);
      d0 = _mm_hadd_epi16(d0, d2);
      __m128i eight = _mm_set1_epi16(8);
      d0 = _mm_add_epi16(d0, eight);
      d0 = _mm_srai_epi16(d0, 4);
      d0 = _mm_packus_epi16(d0, d0);
      __m128i out0 = _mm_lddqu_si128((__m128i *)out);
      __m128i n0 = _mm_set1_epi8(n_out);
      __m128i mask = _mm_cmpgt_epi8(n0, iden);
      out0 = _mm_blendv_epi8(out0, d0, mask);
      _mm_storel_epi64((__m128i *)out, out0);
      __m128i in1 = _mm_lddqu_si128((__m128i *)(in + 16));
      in0 = _mm_alignr_epi8(in1, in0, 8);
      in += 8;
      out += 8;
      len -= n_out;
    }
  }
}

void av1_filter_intra_edge_high_sse4_1(uint16_t *p, int sz, int strength) {
  if (!strength) return;

  DECLARE_ALIGNED(16, static const int16_t, kern[3][8]) = {
    { 4, 8, 4, 8, 4, 8, 4, 8 },  // strength 1: 4,8,4
    { 5, 6, 5, 6, 5, 6, 5, 6 },  // strength 2: 5,6,5
    { 2, 4, 2, 4, 2, 4, 2, 4 }   // strength 3: 2,4,4,4,2
  };

  DECLARE_ALIGNED(16, static const int16_t,
                  v_const[1][8]) = { { 0, 1, 2, 3, 4, 5, 6, 7 } };

  // Extend the first and last samples to simplify the loop for the 5-tap case
  p[-1] = p[0];
  __m128i last = _mm_set1_epi16(p[sz - 1]);
  _mm_storeu_si128((__m128i *)&p[sz], last);

  // Adjust input pointer for filter support area
  uint16_t *in = (strength == 3) ? p - 1 : p;

  // Avoid modifying first sample
  uint16_t *out = p + 1;
  int len = sz - 1;

  const int use_3tap_filter = (strength < 3);

  if (use_3tap_filter) {
    __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
    __m128i iden = _mm_lddqu_si128((__m128i *)v_const[0]);
    __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
    __m128i in8 = _mm_lddqu_si128((__m128i *)&in[8]);
    while (len > 0) {
      int n_out = (len < 8) ? len : 8;
      __m128i in1 = _mm_alignr_epi8(in8, in0, 2);
      __m128i in2 = _mm_alignr_epi8(in8, in0, 4);
      __m128i in02 = _mm_add_epi16(in0, in2);
      __m128i d0 = _mm_unpacklo_epi16(in02, in1);
      __m128i d1 = _mm_unpackhi_epi16(in02, in1);
      d0 = _mm_mullo_epi16(d0, coef0);
      d1 = _mm_mullo_epi16(d1, coef0);
      d0 = _mm_hadd_epi16(d0, d1);
      __m128i eight = _mm_set1_epi16(8);
      d0 = _mm_add_epi16(d0, eight);
      d0 = _mm_srli_epi16(d0, 4);
      __m128i out0 = _mm_lddqu_si128((__m128i *)out);
      __m128i n0 = _mm_set1_epi16(n_out);
      __m128i mask = _mm_cmpgt_epi16(n0, iden);
      out0 = _mm_blendv_epi8(out0, d0, mask);
      _mm_storeu_si128((__m128i *)out, out0);
      in += 8;
      in0 = in8;
      in8 = _mm_lddqu_si128((__m128i *)&in[8]);
      out += 8;
      len -= n_out;
    }
  } else {  // 5-tap filter
    __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
    __m128i iden = _mm_lddqu_si128((__m128i *)v_const[0]);
    __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
    __m128i in8 = _mm_lddqu_si128((__m128i *)&in[8]);
    while (len > 0) {
      int n_out = (len < 8) ? len : 8;
      __m128i in1 = _mm_alignr_epi8(in8, in0, 2);
      __m128i in2 = _mm_alignr_epi8(in8, in0, 4);
      __m128i in3 = _mm_alignr_epi8(in8, in0, 6);
      __m128i in4 = _mm_alignr_epi8(in8, in0, 8);
      __m128i in04 = _mm_add_epi16(in0, in4);
      __m128i in123 = _mm_add_epi16(in1, in2);
      in123 = _mm_add_epi16(in123, in3);
      __m128i d0 = _mm_unpacklo_epi16(in04, in123);
      __m128i d1 = _mm_unpackhi_epi16(in04, in123);
      d0 = _mm_mullo_epi16(d0, coef0);
      d1 = _mm_mullo_epi16(d1, coef0);
      d0 = _mm_hadd_epi16(d0, d1);
      __m128i eight = _mm_set1_epi16(8);
      d0 = _mm_add_epi16(d0, eight);
      d0 = _mm_srli_epi16(d0, 4);
      __m128i out0 = _mm_lddqu_si128((__m128i *)out);
      __m128i n0 = _mm_set1_epi16(n_out);
      __m128i mask = _mm_cmpgt_epi16(n0, iden);
      out0 = _mm_blendv_epi8(out0, d0, mask);
      _mm_storeu_si128((__m128i *)out, out0);
      in += 8;
      in0 = in8;
      in8 = _mm_lddqu_si128((__m128i *)&in[8]);
      out += 8;
      len -= n_out;
    }
  }
}

void av1_upsample_intra_edge_sse4_1(uint8_t *p, int sz) {
  // interpolate half-sample positions
  assert(sz <= 24);

  DECLARE_ALIGNED(16, static const int8_t, kernel[1][16]) = {
    { -1, 9, 9, -1, -1, 9, 9, -1, -1, 9, 9, -1, -1, 9, 9, -1 }
  };

  DECLARE_ALIGNED(16, static const int8_t, v_const[2][16]) = {
    { 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6 },
    { 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10 }
  };

  // Extend first/last samples (upper-left p[-1], last p[sz-1])
  // to support 4-tap filter
  p[-2] = p[-1];
  p[sz] = p[sz - 1];

  uint8_t *in = &p[-2];
  uint8_t *out = &p[-2];

  int n = sz + 1;  // Input length including upper-left sample

  __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
  __m128i in16 = _mm_lddqu_si128((__m128i *)&in[16]);

  __m128i coef0 = _mm_lddqu_si128((__m128i *)kernel[0]);
  __m128i shuf0 = _mm_lddqu_si128((__m128i *)v_const[0]);
  __m128i shuf1 = _mm_lddqu_si128((__m128i *)v_const[1]);

  while (n > 0) {
    __m128i in8 = _mm_alignr_epi8(in16, in0, 8);
    __m128i d0 = _mm_shuffle_epi8(in0, shuf0);
    __m128i d1 = _mm_shuffle_epi8(in0, shuf1);
    __m128i d2 = _mm_shuffle_epi8(in8, shuf0);
    __m128i d3 = _mm_shuffle_epi8(in8, shuf1);
    d0 = _mm_maddubs_epi16(d0, coef0);
    d1 = _mm_maddubs_epi16(d1, coef0);
    d2 = _mm_maddubs_epi16(d2, coef0);
    d3 = _mm_maddubs_epi16(d3, coef0);
    d0 = _mm_hadd_epi16(d0, d1);
    d2 = _mm_hadd_epi16(d2, d3);
    __m128i eight = _mm_set1_epi16(8);
    d0 = _mm_add_epi16(d0, eight);
    d2 = _mm_add_epi16(d2, eight);
    d0 = _mm_srai_epi16(d0, 4);
    d2 = _mm_srai_epi16(d2, 4);
    d0 = _mm_packus_epi16(d0, d2);
    __m128i in1 = _mm_alignr_epi8(in16, in0, 1);
    __m128i out0 = _mm_unpacklo_epi8(in1, d0);
    __m128i out1 = _mm_unpackhi_epi8(in1, d0);
    _mm_storeu_si128((__m128i *)&out[0], out0);
    _mm_storeu_si128((__m128i *)&out[16], out1);
    in0 = in16;
    in16 = _mm_setzero_si128();
    out += 32;
    n -= 16;
  }
}

void av1_upsample_intra_edge_high_sse4_1(uint16_t *p, int sz, int bd) {
  // interpolate half-sample positions
  assert(sz <= 24);

  DECLARE_ALIGNED(16, static const int16_t,
                  kernel[1][8]) = { { -1, 9, -1, 9, -1, 9, -1, 9 } };

  // Extend first/last samples (upper-left p[-1], last p[sz-1])
  // to support 4-tap filter
  p[-2] = p[-1];
  p[sz] = p[sz - 1];

  uint16_t *in = &p[-2];
  uint16_t *out = in;
  int n = sz + 1;

  __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
  __m128i in8 = _mm_lddqu_si128((__m128i *)&in[8]);
  __m128i in16 = _mm_lddqu_si128((__m128i *)&in[16]);
  __m128i in24 = _mm_lddqu_si128((__m128i *)&in[24]);

  while (n > 0) {
    __m128i in1 = _mm_alignr_epi8(in8, in0, 2);
    __m128i in2 = _mm_alignr_epi8(in8, in0, 4);
    __m128i in3 = _mm_alignr_epi8(in8, in0, 6);
    __m128i sum0 = _mm_add_epi16(in0, in3);
    __m128i sum1 = _mm_add_epi16(in1, in2);
    __m128i d0 = _mm_unpacklo_epi16(sum0, sum1);
    __m128i d1 = _mm_unpackhi_epi16(sum0, sum1);
    __m128i coef0 = _mm_lddqu_si128((__m128i *)kernel[0]);
    d0 = _mm_madd_epi16(d0, coef0);
    d1 = _mm_madd_epi16(d1, coef0);
    __m128i eight = _mm_set1_epi32(8);
    d0 = _mm_add_epi32(d0, eight);
    d1 = _mm_add_epi32(d1, eight);
    d0 = _mm_srai_epi32(d0, 4);
    d1 = _mm_srai_epi32(d1, 4);
    d0 = _mm_packus_epi32(d0, d1);
    __m128i max0 = _mm_set1_epi16((1 << bd) - 1);
    d0 = _mm_min_epi16(d0, max0);
    __m128i out0 = _mm_unpacklo_epi16(in1, d0);
    __m128i out1 = _mm_unpackhi_epi16(in1, d0);
    _mm_storeu_si128((__m128i *)&out[0], out0);
    _mm_storeu_si128((__m128i *)&out[8], out1);
    in0 = in8;
    in8 = in16;
    in16 = in24;
    in24 = _mm_setzero_si128();
    out += 16;
    n -= 8;
  }
}

#include "aom_dsp/x86/synonyms.h"

// Get the shift (up-scaled by 256) in X w.r.t a unit change in Y.
// If angle > 0 && angle < 90, dx = -((int)(256 / t));
// If angle > 90 && angle < 180, dx = (int)(256 / t);
// If angle > 180 && angle < 270, dx = 1;
static INLINE int get_dx(int angle) {
  if (angle > 0 && angle < 90) {
    return dr_intra_derivative[angle];
  } else if (angle > 90 && angle < 180) {
    return dr_intra_derivative[180 - angle];
  } else {
    // In this case, we are not really going to use dx. We may return any value.
    return 1;
  }
}

// Get the shift (up-scaled by 256) in Y w.r.t a unit change in X.
// If angle > 0 && angle < 90, dy = 1;
// If angle > 90 && angle < 180, dy = (int)(256 * t);
// If angle > 180 && angle < 270, dy = -((int)(256 * t));
static INLINE int get_dy(int angle) {
  if (angle > 90 && angle < 180) {
    return dr_intra_derivative[angle - 90];
  } else if (angle > 180 && angle < 270) {
    return dr_intra_derivative[270 - angle];
  } else {
    // In this case, we are not really going to use dy. We may return any value.
    return 1;
  }
}

void test_z1(int ang, int dx, int dy, int upsample_above, int upsample_left) {
  (void) upsample_left;
  (void) dy;
  int bh = 16;
  int bw = 64;

  const int max_base_x = ((bw + bh) - 1) << upsample_above;
  const int frac_bits = 6 - upsample_above;
  const int base_inc = 1 << upsample_above;


  if (dx == 0) {
    printf(">>>>>ang: %d dx: %d \n", ang, dx);
    return;
  }

  int x = dx;

  for (int r = 0; r < bh; ++r, x += dx) {
    int base = x >> frac_bits;
    printf("%3d: %8d ", r, base);
    for (int c = 0; c < bw; ++c, base += base_inc) {
      if (base < max_base_x) {
        printf(" %8d", base);
      } else {
//        printf(" (%6d)", max_base_x);
        printf(" (%6d)", base);
      }
    }
    printf("\n");
  }
  printf("ang: %d dx: %d dy: %d\n", ang, dx, dy);
}


// Directional prediction, zone 1: 0 < angle < 90
void av1_dr_prediction_z1_sse4_1(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                                 const uint8_t *above, const uint8_t *left,
                                 int upsample_above, int dx, int dy) {
  (void)left;
  (void)dy;
  assert(dy == 1);
  assert(dx > 0);

#if 0
  for (int ang = 1; ang < 90; ++ang) {
    test_z1(ang, get_dx(ang), get_dy(ang), 0, 0);
    printf("~~~~~~~~~~~~~~~~~ bw: %d bh: %d\n", bw, bh);
  }
  exit(0);
#endif

  const int max_base_x = ((bw + bh) - 1) << upsample_above;
  const int frac_bits = 6 - upsample_above;
  const __m128i dup16 =
      _mm_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100);
  const __m128i mbx = _mm_cvtsi32_si128(max_base_x);
  const __m128i v_max_base_x = _mm_shuffle_epi8(mbx, dup16);
  const __m128i ambx = _mm_cvtsi32_si128(above[max_base_x]);
  const __m128i above_max_base_x = _mm_shuffle_epi8(ambx, dup16);
  const __m128i offsets =
      _mm_set_epi32(0x00070006, 0x00050004, 0x00030002, 0x00010000);

  int x = dx;

  if (bw > 4) {
    if(!upsample_above) {
      const int base_inc = 8;
      const __m128i v_base_inc = _mm_set1_epi16(8);

      for (int h = 0; h < bh; ++h, dst += stride, x += dx) {
        int base = x >> frac_bits;
        const int shift = (x & 0x3F) >> 1;
if (1)
        if (base >= max_base_x) {
          for (int i = h; i < bh; ++i) {
            memset(dst, above[max_base_x], bw * sizeof(dst[0]));
            dst += stride;
          }
          return;
        }

        const __m128i base0 = _mm_cvtsi32_si128(base);
        const __m128i base1 = _mm_shuffle_epi8(base0, dup16);
        const __m128i v_32mshift_shift =
            _mm_shuffle_epi8(_mm_cvtsi32_si128((32 - shift) | (shift << 8)), dup16);
        __m128i v_base = _mm_add_epi16(base1, offsets);
        for (int w = 0; w < bw; w += 8, base += base_inc) {
          const __m128i ab0 = xx_loadu_128(&above[base]);
          const __m128i ab1 = _mm_srli_si128(ab0, 1);  // above[base + 1]
          const __m128i ab = _mm_unpacklo_epi8(ab0, ab1);
          const __m128i e =_mm_maddubs_epi16(ab, v_32mshift_shift);
          const __m128i f = xx_roundn_epu16(e, 5);
          const __m128i _m = _mm_cmplt_epi16(v_base, v_max_base_x);
          const __m128i final0 = _mm_blendv_epi8(above_max_base_x, f, _m);
          const __m128i final1 = _mm_packus_epi16(final0, final0);
          const __m128i final = final1; // _mm_shuffle_epi8(final1, gat);
          v_base = _mm_add_epi16(v_base, v_base_inc);
          xx_storel_64(&dst[w], final);
        }
      }
    } else {
      const __m128i gat_e = _mm_set_epi32(0, 0, 0x0e0c0a08, 0x06040200);
      const __m128i gat_o = _mm_set_epi32(0, 0, 0x0f0d0b09, 0x07050301);
      const int base_inc = (1 << upsample_above) * 8;
      const __m128i v_base_inc = _mm_set1_epi16(16);
      const __m128i offsets_us = _mm_slli_epi16(offsets, 1);

      for (int h = 0; h < bh; ++h, dst += stride, x += dx) {
        int base = x >> frac_bits;
//        const int shift = ((x << upsample_above) & 0x3F) >> 1;
        const int shift = x & 0x1F;

        if (base >= max_base_x) {
          for (int i = h; i < bh; ++i) {
            memset(dst, above[max_base_x], bw * sizeof(dst[0]));
            dst += stride;
          }
          return;
        }

        const __m128i base0 = _mm_cvtsi32_si128(base);
        const __m128i base1 = _mm_shuffle_epi8(base0, dup16);
        const __m128i v_32mshift_shift =
            _mm_shuffle_epi8(_mm_cvtsi32_si128((32 - shift) | (shift << 8)), dup16);
        __m128i v_base = _mm_add_epi16(base1, offsets_us);

        for (int w = 0; w < bw; w += 8, base += base_inc) {
          const __m128i ab0 = xx_loadu_128(&above[base]);
          const __m128i ab0_e = _mm_shuffle_epi8(ab0, gat_e);
          const __m128i ab0_o = _mm_shuffle_epi8(ab0, gat_o);
          const __m128i ab = _mm_unpacklo_epi8(ab0_e, ab0_o);
          const __m128i e =_mm_maddubs_epi16(ab, v_32mshift_shift);
          const __m128i f = xx_roundn_epu16(e, 5);
          const __m128i _m = _mm_cmplt_epi16(v_base, v_max_base_x);
          const __m128i final0 = _mm_blendv_epi8(above_max_base_x, f, _m);
          const __m128i final1 = _mm_packus_epi16(final0, final0);
          v_base = _mm_add_epi16(v_base, v_base_inc);
          xx_storel_64(&dst[w], final1);
        }
      }
    }
  } else {
    const __m128i gat = upsample_above
                            ? _mm_set_epi32(0, 0, 0x0e0c0a08, 0x06040200)
                            : _mm_set_epi32(0, 0, 0x07060504, 0x03020100);
    for (int h = 0; h < bh; ++h, dst += stride, x += dx) {
      int base = x >> frac_bits;
      const int shift = ((x << upsample_above) & 0x3F) >> 1;

      if (base >= max_base_x) {
        for (int i = h; i < bh; ++i) {
          memset(dst, above[max_base_x], 4 * sizeof(dst[0]));
          dst += stride;
        }
        return;
      }

      const __m128i base0 = _mm_cvtsi32_si128(base);
      const __m128i base1 = _mm_shuffle_epi8(base0, dup16);
      const __m128i v_shift = _mm_shuffle_epi8(_mm_cvtsi32_si128(shift), dup16);
      const __m128i v_32mshift =
          _mm_shuffle_epi8(_mm_cvtsi32_si128(32 - shift), dup16);
      __m128i v_base = _mm_add_epi16(base1, offsets);
      {
        const __m128i ab0 = xx_loadl_64(&above[base]);
        const __m128i a = _mm_cvtepu8_epi16(ab0);
        const __m128i b = _mm_srli_si128(a, 2);  // above[base + 1]
        const __m128i c = _mm_mullo_epi16(a, v_32mshift);
        const __m128i d = _mm_mullo_epi16(b, v_shift);
        const __m128i e = _mm_add_epi16(c, d);
        const __m128i f = xx_roundn_epu16(e, 5);
        const __m128i _m = _mm_cmplt_epi16(v_base, v_max_base_x);
        const __m128i final0 = _mm_blendv_epi8(above_max_base_x, f, _m);
        const __m128i final1 = _mm_packus_epi16(final0, final0);
        const __m128i final = _mm_shuffle_epi8(final1, gat);
        xx_storel_32(&dst[0], final);
      }
    }

  }
}

static void get_base_y(__m128i *base0, __m128i *base1, __m128i *v_y,
                       const __m128i *_m, const uint8_t *left,
                       int frac_bits_y) {
  // NOTE: the calculated base index can be invalid.  ie exceed in the negative
  // direction.  Could cause out of bounds reads.
  // Possible solution: keep base in simd and use mask to blend with zero

#if 1
  const __m128i zero = _mm_setzero_si128();
  const __m128i v_base2 = _mm_sra_epi16(*v_y, _mm_cvtsi32_si128(frac_bits_y));
  // if (base1 >= min_base_x), replace calculated y (invalid) with zero to
  // prevent out-ranging.
  const __m128i v_base2_0 = _mm_blendv_epi8(zero, v_base2, *_m);
  // sign extended left base indexes
  const __m128i v_base2_s32 = _mm_cvtepi16_epi32(v_base2_0);

  int base2 = _mm_extract_epi32(v_base2_s32, 0);

  uint16_t *left16 = (uint16_t *)&left[base2];
  __m128i a = zero;
  a = _mm_insert_epi16 (a, *left16, 0);
  base2 = _mm_extract_epi32(v_base2_s32, 1);
  left16 = (uint16_t *)&left[base2];
  a = _mm_insert_epi16 (a, *left16, 1);
  base2 = _mm_extract_epi32(v_base2_s32, 2);
  left16 = (uint16_t *)&left[base2];
  a = _mm_insert_epi16 (a, *left16, 2);
  base2 = _mm_extract_epi32(v_base2_s32, 3);
  left16 = (uint16_t *)&left[base2];
  a = _mm_insert_epi16 (a, *left16, 3);

  *base1 = _mm_srli_epi16(a, 8);  // base + 1 expanded to 16 bits
  a = _mm_slli_epi16(a, 8);       // clear out base + 1
  *base0 = _mm_srli_epi16(a, 8);  // base + 0 expanded to 16 bits
#else
  // NOTE: broken... :-(

//  const int base_end = y >> frac_bits_y;
  y -= dy * 4;  // y start
  const int base_start = y >> frac_bits_y;

  __m128i left_a = xx_loadu_128(&left[base_start]);
  const __m128i dup16 =
      _mm_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100);

  const __m128i v_y = _mm_shuffle_epi8(_mm_cvtsi32_si128(y), dup16);
  const __m128i v_dy = _mm_shuffle_epi8(_mm_cvtsi32_si128(dy), dup16);
  const __m128i v_dy_off =
      _mm_set_epi32(0x0, 0x0, 0x00000001, 0x00020003);
  const __m128i v_dy_reverse = _mm_mullo_epi16(v_dy, v_dy_off);
  const __m128i v_y_dy = _mm_add_epi16(v_y, v_dy_reverse);
  const __m128i v_y_dy_shift = _mm_sra_epi16(v_y_dy, _mm_cvtsi32_si128(frac_bits_y));

  // calculate left index used in final shuffle
  const __m128i v_left_index0 = _mm_sub_epi16(v_y_dy_shift, _mm_shuffle_epi8(_mm_cvtsi32_si128(base_start), dup16));
  const __m128i v_left_index1 = _mm_slli_epi16(v_left_index0, 8);
  const __m128i v_left_index = _mm_add_epi16(_mm_srli_epi16(v_left_index1, 8), _mm_set1_epi16(0x8000));
  const __m128i v_left_indexp1 = _mm_add_epi16(v_left_index, _mm_set1_epi16(0x8001));

  const __m128i base_a = _mm_shuffle_epi8(left_a, v_left_index);
  const __m128i base_b = _mm_shuffle_epi8(left_a, v_left_indexp1);

  *base0 = base_a;
  *base1 = base_b;

  // - read 8/16 bytes left[base2] <- using y starting point
  // - save base starting index. use this to reset index
  // - simd y >> frac_bits_y <- should be reversed order
  // - simd sub (sat) y starting point
  // - simd pshufb left,new index
  // - simd expand to 16bits using shifts
#endif
}

void test_z2(int ang, int dx, int dy, int upsample_above, int upsample_left) {
  const int min_base_x = -(1 << upsample_above);
  const int frac_bits_x = 6 - upsample_above;
  const int frac_bits_y = 6 - upsample_left;
  const int base_inc_x = 1 << upsample_above;
  int x = -dx;

  int bh = 16;
  int bw = 16;

  if (dx == 0 && dy == 0) {
    printf(">>>>>ang: %d dx: %d dy: %d\n", ang, dx, dy);
    return;
  }

  for (int r = 0; r < bh; ++r, x -= dx) {
    int base1 = x >> frac_bits_x;
    int y = (r << 6) - dy;
    printf("%3d: %8d ", r, y);
    for (int w = 0; w < bw; w += 1, base1 += base_inc_x, y -= dy) {
      if (base1 >= min_base_x) {
        int base2 = y >> frac_bits_y;
        printf(" (%6d)", base2);
      } else {
        int base2 = y >> frac_bits_y;
        printf(" %8d", base2);
      }
    }
    printf("\n");
  }
  printf("ang: %d dx: %d dy: %d\n", ang, dx, dy);
}

// Directional prediction, zone 2: 90 < angle < 180
void av1_dr_prediction_z2_sse4_1(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                                 const uint8_t *above, const uint8_t *left,
                                 int upsample_above, int upsample_left, int dx,
                                 int dy) {
  assert(dx > 0);
  assert(dy > 0);

  // FIXME:
  if (0) {
//  if (upsample_above || upsample_left) {
//  if ( upsample_above) {
//printf("upsample_above || upsample_left\n");
    av1_dr_prediction_z2_c(dst, stride, bw, bh,
                                     above, left,
                                     upsample_above, upsample_left, dx,
                                     dy);
    return;
  }

#if 0
  for (int ang = 91; ang < 180; ++ang) {
    test_z2(ang, get_dx(ang), get_dy(ang), 0, 1);
    printf("~~~~~~~~~~~~~~~~~ bw: %d bh: %d\n", bw, bh);
  }
  exit(0);
#endif

  const int min_base_x = -(1 << upsample_above);
  const int frac_bits_x = 6 - upsample_above;
  const int frac_bits_y = 6 - upsample_left;
  const int base_inc_x = 1 << upsample_above;

  const __m128i gat = upsample_above
                          ? _mm_set_epi32(0, 0, 0x80068004, 0x80028000)
                          : _mm_set_epi32(0, 0, 0x80038002, 0x80018000);
  const __m128i v_base_inc =
      upsample_above ? _mm_set1_epi16(8) : _mm_set1_epi16(4);
  const __m128i dup16 =
      _mm_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100);
  const __m128i v_min_base_x =
      _mm_shuffle_epi8(_mm_cvtsi32_si128(min_base_x), dup16);
  const __m128i v_dy = _mm_shuffle_epi8(_mm_cvtsi32_si128(dy), dup16);
  const __m128i v_dy_scale =
      _mm_set_epi32(0x00070006, 0x00050004, 0x00030002, 0x00010000);

  // dy * 0, dy * 1, dy * 2, dy * 3, ...
  const __m128i v_dy_offsets = _mm_mullo_epi16(v_dy, v_dy_scale);
  const __m128i v_dy4 = _mm_shuffle_epi8(_mm_cvtsi32_si128(dy*4), dup16);

  int x = -dx;

  const __m128i offsets = upsample_above
                          ? _mm_set_epi32(0, 0, 0x00060004, 0x00020000)
                          : _mm_set_epi32(0, 0, 0x00030002, 0x00010000);

  for (int r = 0; r < bh; ++r, x -= dx, dst += stride) {
    int base1 = x >> frac_bits_x;
    int y = (r << 6) - dy;

    __m128i v_base = _mm_shuffle_epi8(_mm_cvtsi32_si128(base1), dup16);
    v_base = _mm_add_epi16(v_base, offsets);

    const __m128i v_y_base = _mm_shuffle_epi8(_mm_cvtsi32_si128(y), dup16);
    const int shift_x = ((x << upsample_above) & 0x3F) >> 1;
    const __m128i v_shift_x =
        _mm_shuffle_epi8(_mm_cvtsi32_si128(shift_x), dup16);

    __m128i v_y = _mm_sub_epi16(v_y_base, v_dy_offsets);

    for (int w = 0; w < bw; w += 4, base1 += base_inc_x * 4) {
      // mask to blend 16 bit values
      const __m128i _m = _mm_cmplt_epi16(v_base, v_min_base_x);

      v_base = _mm_add_epi16(v_base, v_base_inc);
      const __m128i v_base0_x0 = xx_loadl_64(&above[base1]);

//      const __m128i v_base0_x0 = xx_loadl_64(&above[base1]);
      const __m128i v_base1_x0 = _mm_srli_si128(v_base0_x0, 1);
      const __m128i v_base0_x = _mm_shuffle_epi8(v_base0_x0, gat);
      const __m128i v_base1_x = _mm_shuffle_epi8(v_base1_x0, gat);

      // FIXME: left + base index can outrange....
      __m128i v_base0_y;
      __m128i v_base1_y;
       get_base_y(&v_base0_y, &v_base1_y, &v_y, &_m, left, frac_bits_y);

       // blend bases (above, left) here
       const __m128i v_base0 = _mm_blendv_epi8(v_base0_x, v_base0_y, _m);
       const __m128i v_base1 = _mm_blendv_epi8(v_base1_x, v_base1_y, _m);

      const __m128i v_shift_y0 =
          _mm_sll_epi16(v_y, _mm_cvtsi32_si128(upsample_left));
      const __m128i v_shift_y1 = _mm_slli_epi16(v_shift_y0, 10);  // & 0x3f
      const __m128i v_shift_y = _mm_srli_epi16(v_shift_y1, 11);

      // blend v_shift_y, v_shift_x here
      const __m128i v_shift = _mm_blendv_epi8(v_shift_x, v_shift_y, _m);
      const __m128i v_32mshift = _mm_sub_epi16(_mm_set1_epi16(32), v_shift);

      const __m128i a = v_base0;
      const __m128i b = v_base1;
      const __m128i c = _mm_mullo_epi16(a, v_32mshift);
      const __m128i d = _mm_mullo_epi16(b, v_shift);
      const __m128i e = _mm_add_epi16(c, d);
      const __m128i f = xx_roundn_epu16(e, 5);
      const __m128i g = _mm_packus_epi16(f, f);

      v_y = _mm_sub_epi16(v_y, v_dy4);

      xx_storel_32(&dst[w], g);
    }
  }
}

void test_z3a(int ang, int dy, int upsample_left, int bw, int bh) {
  const int frac_bits = 6 - upsample_left;
  const int base_inc = 1 << upsample_left;
  int y = dy;

  if (dy == 0) {
    printf(">>>>>ang: %d dy: %d\n", ang, dy);
    return;
  }

  const int max_base_y = (bw + bh - 1) << upsample_left;
  int y_start = y;

  for (int c = 0; c < bw; ++c, y += dy) {
    int base = y >> frac_bits;
    printf("%3d: %8d ", c, y_start);

    for (int r = 0; r < bh; ++r, base += base_inc) {
      if (base < max_base_y) {
        printf(" %8d", base);
      } else {
//        base = max_base_y;
        printf(" (%6d)", base);
      }
    }
    printf("\n");
  }

  y_start += bh - 0;
  int base_max_diff = (y_start + dy * bw)  >> frac_bits;
  if (base_max_diff >= max_base_y)
    base_max_diff = max_base_y;
  base_max_diff -= (y_start >> frac_bits);

  printf("ang: %d dy: %d base_max_diff %d\n", ang, dy, base_max_diff);
}

// Directional prediction, zone 3: 180 < angle < 270
void av1_dr_prediction_z3_sse4_1(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                                 const uint8_t *above, const uint8_t *left,
                                 int upsample_left, int dx, int dy) {
  int r, c;

  (void)above;
  (void)dx;

  assert(dx == 1);
  assert(dy > 0);

  const int max_base_y = (bw + bh - 1) << upsample_left;
  const int frac_bits = 6 - upsample_left;
//  const int base_inc = 1 << upsample_left;

  int y = dy;

//  uint8_t left_buf[64 + 64];
//  uint8_t *left = &left_buf[0];
//  memcpy(left_buf, left_ - 0, 1 * (64 * 2 + 16) * sizeof(*left_));

#if 0
  for (int ang = 181; ang < 270; ++ang) {
    test_z3a(ang, get_dy(ang), 0, 64, 64);
    printf("~~~~~~~~~~~~~~~~~ bw: %d bh: %d\n", bw, bh);
  }
  exit(0);
#endif

  const __m128i v_base_inc =
      upsample_left ? _mm_set1_epi16(8) : _mm_set1_epi16(4);
  const __m128i dup16 =
      _mm_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100);
  const __m128i v_max_base_y =
      _mm_shuffle_epi8(_mm_cvtsi32_si128(max_base_y), dup16);
  const __m128i lmby = _mm_cvtsi32_si128(left[max_base_y]);
  const __m128i left_max_base_x = _mm_shuffle_epi8(lmby, dup16);
  const __m128i gat = upsample_left
                          ? _mm_set_epi32(0, 0, 0x0e0c0a08, 0x06040200)
                          : _mm_set_epi32(0, 0, 0x07060504, 0x03020100);
  const __m128i offsets =
      _mm_set_epi32(0x00070006, 0x00050004, 0x00030002, 0x00010000);

  const __m128i zero = _mm_setzero_si128();

  for (c = 0; c < bw; ++c, y += dy) {
    int base = y >> frac_bits;
    const int shift = ((y << upsample_left) & 0x3F) >> 1;
#if 0
    const __m128i v_shift = _mm_shuffle_epi8(_mm_cvtsi32_si128(shift), dup16);
    const __m128i v_32mshift = _mm_sub_epi16(_mm_set1_epi16(32), v_shift);
#else
    const __m128i v_32mshift_shift =
        _mm_shuffle_epi8(_mm_cvtsi32_si128((32 - shift) | (shift << 8)), dup16);
#endif

    __m128i v_base = _mm_shuffle_epi8(_mm_cvtsi32_si128(base), dup16);
    v_base = _mm_add_epi16(v_base, offsets);

    if (1 && base >= max_base_y) {
      for (int i = 0; i < bh; ++i) {
        memset(dst + c, left[max_base_y], (bw - c) * sizeof(dst[0]));
        dst += stride;
      }
      return;
    }

    for (r = 0; r < bh; r += 4) {
      const __m128i _m = _mm_cmplt_epi16(v_base, v_max_base_y);

      // If base exceeds max_base, force to 0 to prevent OOB reads.
      const __m128i v_new_base = v_base; //_mm_blendv_epi8(zero, v_base, _m);
      const int new_base = _mm_extract_epi16(v_new_base, 0);
      v_base = _mm_add_epi16(v_base, v_base_inc);

      const __m128i ab0 = xx_loadl_64(&left[new_base]);
#if 0
      const __m128i a = _mm_cvtepu8_epi16(ab0);
      const __m128i b = _mm_srli_si128(a, 2);  // left[base + 1]
      const __m128i cc = _mm_mullo_epi16(a, v_32mshift);
      const __m128i d = _mm_mullo_epi16(b, v_shift);
      const __m128i e = _mm_add_epi16(cc, d);
#else
      const __m128i ab1 = _mm_srli_si128(ab0, 1);  // left[base + 1]
      const __m128i ab = _mm_unpacklo_epi8(ab0, ab1);
      const __m128i e =_mm_maddubs_epi16(ab, v_32mshift_shift);
#endif
      const __m128i f = xx_roundn_epu16(e, 5);
      const __m128i final0 = _mm_blendv_epi8(left_max_base_x, f, _m);
      const __m128i final1 = _mm_packus_epi16(final0, final0);
      const __m128i final = _mm_shuffle_epi8(final1, gat);

      dst[(r + 0) * stride + c] = (uint8_t)_mm_extract_epi8(final, 0);
      dst[(r + 1) * stride + c] = (uint8_t)_mm_extract_epi8(final, 1);
      dst[(r + 2) * stride + c] = (uint8_t)_mm_extract_epi8(final, 2);
      dst[(r + 3) * stride + c] = (uint8_t)_mm_extract_epi8(final, 3);
    }
  }
}



/*
--- see recontintra.c

upsample_above =
    use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
upsample_left =
    use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);

static int use_intra_edge_upsample(int bs0, int bs1, int delta, int type) {
  const int d = abs(delta);
  const int blk_wh = bs0 + bs1;
  if (d <= 0 || d >= 40) return 0;
  return type ? (blk_wh <= 8) : (blk_wh <= 16);
}
*/


/*
const int16_t dr_intra_derivative[90] = {
  // More evenly spread out angles and limited to 10-bit
  // Values that are 0 will never be used
  //                    Approx angle
  0,    0, 0,        //
  1023, 0, 0,        // 3, ...
  547,  0, 0,        // 6, ...
  372,  0, 0, 0, 0,  // 9, ...
  273,  0, 0,        // 14, ...
  215,  0, 0,        // 17, ...
  178,  0, 0,        // 20, ...
  151,  0, 0,        // 23, ... (113 & 203 are base angles)
  132,  0, 0,        // 26, ...
  116,  0, 0,        // 29, ...
  102,  0, 0, 0,     // 32, ...
  90,   0, 0,        // 36, ...
  80,   0, 0,        // 39, ...
  71,   0, 0,        // 42, ...
  64,   0, 0,        // 45, ... (45 & 135 are base angles)
  57,   0, 0,        // 48, ...
  51,   0, 0,        // 51, ...
  45,   0, 0, 0,     // 54, ...
  40,   0, 0,        // 58, ...
  35,   0, 0,        // 61, ...
  31,   0, 0,        // 64, ...
  27,   0, 0,        // 67, ... (67 & 157 are base angles)
  23,   0, 0,        // 70, ...
  19,   0, 0,        // 73, ...
  15,   0, 0, 0, 0,  // 76, ...
  11,   0, 0,        // 81, ...
  7,    0, 0,        // 84, ...
  3,    0, 0,        // 87, ...
};
*/
