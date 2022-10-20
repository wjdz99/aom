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

#include <emmintrin.h>  // SSE2

#include "config/aom_dsp_rtcd.h"

#include "aom_dsp/txfm_common.h"
#include "aom_dsp/x86/fwd_txfm_sse2.h"
#include "aom_dsp/x86/txfm_common_sse2.h"
#include "aom_ports/mem.h"

#define ADD_EPI16 _mm_add_epi16
#define SUB_EPI16 _mm_sub_epi16

static void FDCT4x4_2D_HELPER(const int16_t *input, int stride, __m128i *in0,
                              __m128i *in1) {
  // Constants
  // These are the coefficients used for the multiplies.
  // In the comments, pN means cos(N pi /64) and mN is -cos(N pi /64),
  // where cospi_N_64 = cos(N pi /64)
  const __m128i k__cospi_A =
      octa_set_epi16(cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64,
                     cospi_16_64, -cospi_16_64, cospi_16_64, -cospi_16_64);
  const __m128i k__cospi_B =
      octa_set_epi16(cospi_16_64, -cospi_16_64, cospi_16_64, -cospi_16_64,
                     cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64);
  const __m128i k__cospi_C =
      octa_set_epi16(cospi_8_64, cospi_24_64, cospi_8_64, cospi_24_64,
                     cospi_24_64, -cospi_8_64, cospi_24_64, -cospi_8_64);
  const __m128i k__cospi_D =
      octa_set_epi16(cospi_24_64, -cospi_8_64, cospi_24_64, -cospi_8_64,
                     cospi_8_64, cospi_24_64, cospi_8_64, cospi_24_64);
  const __m128i k__cospi_E =
      octa_set_epi16(cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64,
                     cospi_16_64, cospi_16_64, cospi_16_64, cospi_16_64);
  const __m128i k__cospi_F =
      octa_set_epi16(cospi_16_64, -cospi_16_64, cospi_16_64, -cospi_16_64,
                     cospi_16_64, -cospi_16_64, cospi_16_64, -cospi_16_64);
  const __m128i k__cospi_G =
      octa_set_epi16(cospi_8_64, cospi_24_64, cospi_8_64, cospi_24_64,
                     -cospi_8_64, -cospi_24_64, -cospi_8_64, -cospi_24_64);
  const __m128i k__cospi_H =
      octa_set_epi16(cospi_24_64, -cospi_8_64, cospi_24_64, -cospi_8_64,
                     -cospi_24_64, cospi_8_64, -cospi_24_64, cospi_8_64);

  const __m128i k__DCT_CONST_ROUNDING = _mm_set1_epi32(DCT_CONST_ROUNDING);
  // This second rounding constant saves doing some extra adds at the end
  const __m128i k__DCT_CONST_ROUNDING2 =
      _mm_set1_epi32(DCT_CONST_ROUNDING + (DCT_CONST_ROUNDING << 1));
  const int DCT_CONST_BITS2 = DCT_CONST_BITS + 2;
  const __m128i k__nonzero_bias_a = _mm_setr_epi16(0, 1, 1, 1, 1, 1, 1, 1);
  const __m128i k__nonzero_bias_b = _mm_setr_epi16(1, 0, 0, 0, 0, 0, 0, 0);

  // Load inputs.
  *in0 = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
  *in1 = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
  *in1 = _mm_unpacklo_epi64(
      *in1, _mm_loadl_epi64((const __m128i *)(input + 2 * stride)));
  *in0 = _mm_unpacklo_epi64(
      *in0, _mm_loadl_epi64((const __m128i *)(input + 3 * stride)));
  // in0 = [i0 i1 i2 i3 iC iD iE iF]
  // in1 = [i4 i5 i6 i7 i8 i9 iA iB]
  // multiply by 16 to give some extra precision
  *in0 = _mm_slli_epi16(*in0, 4);
  *in1 = _mm_slli_epi16(*in1, 4);
  // if (i == 0 && input[0]) input[0] += 1;
  // add 1 to the upper left pixel if it is non-zero, which helps reduce
  // the round-trip error
  {
    // The mask will only contain whether the first value is zero, all
    // other comparison will fail as something shifted by 4 (above << 4)
    // can never be equal to one. To increment in the non-zero case, we
    // add the mask and one for the first element:
    //   - if zero, mask = -1, v = v - 1 + 1 = v
    //   - if non-zero, mask = 0, v = v + 0 + 1 = v + 1
    __m128i mask = _mm_cmpeq_epi16(*in0, k__nonzero_bias_a);
    *in0 = _mm_add_epi16(*in0, mask);
    *in0 = _mm_add_epi16(*in0, k__nonzero_bias_b);
  }
  // There are 4 total stages, alternating between an add/subtract stage
  // followed by an multiply-and-add stage.
  {
    // Stage 1: Add/subtract

    // in0 = [i0 i1 i2 i3 iC iD iE iF]
    // in1 = [i4 i5 i6 i7 i8 i9 iA iB]
    const __m128i r0 = _mm_unpacklo_epi16(*in0, *in1);
    const __m128i r1 = _mm_unpackhi_epi16(*in0, *in1);
    // r0 = [i0 i4 i1 i5 i2 i6 i3 i7]
    // r1 = [iC i8 iD i9 iE iA iF iB]
    const __m128i r2 = _mm_shuffle_epi32(r0, 0xB4);
    const __m128i r3 = _mm_shuffle_epi32(r1, 0xB4);
    // r2 = [i0 i4 i1 i5 i3 i7 i2 i6]
    // r3 = [iC i8 iD i9 iF iB iE iA]

    const __m128i t0 = _mm_add_epi16(r2, r3);
    const __m128i t1 = _mm_sub_epi16(r2, r3);
    // t0 = [a0 a4 a1 a5 a3 a7 a2 a6]
    // t1 = [aC a8 aD a9 aF aB aE aA]

    // Stage 2: multiply by constants (which gets us into 32 bits).
    // The constants needed here are:
    // k__cospi_A = [p16 p16 p16 p16 p16 m16 p16 m16]
    // k__cospi_B = [p16 m16 p16 m16 p16 p16 p16 p16]
    // k__cospi_C = [p08 p24 p08 p24 p24 m08 p24 m08]
    // k__cospi_D = [p24 m08 p24 m08 p08 p24 p08 p24]
    const __m128i u0 = _mm_madd_epi16(t0, k__cospi_A);
    const __m128i u2 = _mm_madd_epi16(t0, k__cospi_B);
    const __m128i u1 = _mm_madd_epi16(t1, k__cospi_C);
    const __m128i u3 = _mm_madd_epi16(t1, k__cospi_D);
    // Then add and right-shift to get back to 16-bit range
    const __m128i v0 = _mm_add_epi32(u0, k__DCT_CONST_ROUNDING);
    const __m128i v1 = _mm_add_epi32(u1, k__DCT_CONST_ROUNDING);
    const __m128i v2 = _mm_add_epi32(u2, k__DCT_CONST_ROUNDING);
    const __m128i v3 = _mm_add_epi32(u3, k__DCT_CONST_ROUNDING);
    const __m128i w0 = _mm_srai_epi32(v0, DCT_CONST_BITS);
    const __m128i w1 = _mm_srai_epi32(v1, DCT_CONST_BITS);
    const __m128i w2 = _mm_srai_epi32(v2, DCT_CONST_BITS);
    const __m128i w3 = _mm_srai_epi32(v3, DCT_CONST_BITS);
    // w0 = [b0 b1 b7 b6]
    // w1 = [b8 b9 bF bE]
    // w2 = [b4 b5 b3 b2]
    // w3 = [bC bD bB bA]
    const __m128i x0 = _mm_packs_epi32(w0, w1);
    const __m128i x1 = _mm_packs_epi32(w2, w3);

    // x0 = [b0 b1 b7 b6 b8 b9 bF bE]
    // x1 = [b4 b5 b3 b2 bC bD bB bA]
    *in0 = _mm_shuffle_epi32(x0, 0xD8);
    *in1 = _mm_shuffle_epi32(x1, 0x8D);
    // in0 = [b0 b1 b8 b9 b7 b6 bF bE]
    // in1 = [b3 b2 bB bA b4 b5 bC bD]
  }
  {
    // vertical DCTs finished. Now we do the horizontal DCTs.
    // Stage 3: Add/subtract

    const __m128i t0 = ADD_EPI16(*in0, *in1);
    const __m128i t1 = SUB_EPI16(*in0, *in1);

    // Stage 4: multiply by constants (which gets us into 32 bits).
    {
      // The constants needed here are:
      // k__cospi_E = [p16 p16 p16 p16 p16 p16 p16 p16]
      // k__cospi_F = [p16 m16 p16 m16 p16 m16 p16 m16]
      // k__cospi_G = [p08 p24 p08 p24 m08 m24 m08 m24]
      // k__cospi_H = [p24 m08 p24 m08 m24 p08 m24 p08]
      const __m128i u0 = _mm_madd_epi16(t0, k__cospi_E);
      const __m128i u1 = _mm_madd_epi16(t0, k__cospi_F);
      const __m128i u2 = _mm_madd_epi16(t1, k__cospi_G);
      const __m128i u3 = _mm_madd_epi16(t1, k__cospi_H);
      // Then add and right-shift to get back to 16-bit range
      // but this combines the final right-shift as well to save operations
      // This unusual rounding operations is to maintain bit-accurate
      // compatibility with the c version of this function which has two
      // rounding steps in a row.
      const __m128i v0 = _mm_add_epi32(u0, k__DCT_CONST_ROUNDING2);
      const __m128i v1 = _mm_add_epi32(u1, k__DCT_CONST_ROUNDING2);
      const __m128i v2 = _mm_add_epi32(u2, k__DCT_CONST_ROUNDING2);
      const __m128i v3 = _mm_add_epi32(u3, k__DCT_CONST_ROUNDING2);
      const __m128i w0 = _mm_srai_epi32(v0, DCT_CONST_BITS2);
      const __m128i w1 = _mm_srai_epi32(v1, DCT_CONST_BITS2);
      const __m128i w2 = _mm_srai_epi32(v2, DCT_CONST_BITS2);
      const __m128i w3 = _mm_srai_epi32(v3, DCT_CONST_BITS2);
      *in0 = _mm_packs_epi32(w0, w2);
      *in1 = _mm_packs_epi32(w1, w3);
    }
  }
}

void FDCT4x4_2D(const int16_t *input, tran_low_t *output, int stride) {
  // This 2D transform implements 4 vertical 1D transforms followed
  // by 4 horizontal 1D transforms.  The multiplies and adds are as given
  // by Chen, Smith and Fralick ('77).  The commands for moving the data
  // around have been minimized by hand.
  // For the purposes of the comments, the 16 inputs are referred to at i0
  // through iF (in raster order), intermediate variables are a0, b0, c0
  // through f, and correspond to the in-place computations mapped to input
  // locations.  The outputs, o0 through oF are labeled according to the
  // output locations.
  __m128i in0, in1;
  FDCT4x4_2D_HELPER(input, stride, &in0, &in1);

  // Post-condition (v + 1) >> 2 is now incorporated into previous
  // add and right-shift commands.  Only 2 store instructions needed
  // because we are using the fact that 1/3 are stored just after 0/2.
  storeu_output(&in0, output + 0 * 4);
  storeu_output(&in1, output + 2 * 4);
}

void FDCT4x4_2D_LP(const int16_t *input, int16_t *output, int stride) {
  __m128i in0, in1;
  FDCT4x4_2D_HELPER(input, stride, &in0, &in1);
  _mm_storeu_si128((__m128i *)(output + 0 * 4), in0);
  _mm_storeu_si128((__m128i *)(output + 2 * 4), in1);
}

#undef ADD_EPI16
#undef SUB_EPI16
