#include <tmmintrin.h>
#include <immintrin.h>
#include "./av1_rtcd.h"
#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "av1/common/daala_tx.h"
#include "av1/encoder/daala_fwd_txfm.h"
#include "av1/common/x86/daala_tx_avx2.h"

#if CONFIG_DAALA_TX

/* Loads a 4x4 buffer of 16-bit values into 2 SSE registers. */
OD_SIMD_INLINE void od_load_buffer_compact_4x4_epi16(__m128i *q0, __m128i *q1,
                                                     const int16_t *in) {
  *q0 = _mm_loadu_si128((const __m128i *)in + 0);
  *q1 = _mm_loadu_si128((const __m128i *)in + 1);
}

OD_SIMD_INLINE void od_store_halves_buffer_4x4_epi16(int16_t *out, __m128i q0,
                                                     __m128i q1, __m128i q2,
                                                     __m128i q3) {
  q0 = _mm_unpacklo_epi64(q0, q1);
  q1 = _mm_unpacklo_epi64(q2, q3);
  _mm_storeu_si128((__m128i *)out + 0, q0);
  _mm_storeu_si128((__m128i *)out + 1, q1);
}

OD_SIMD_INLINE void od_store_buffer_4x8_epi32(tran_low_t *out, __m256i q0,
                                              __m256i q1, __m256i q2,
                                              __m256i q3) {
  _mm256_storeu_si256((__m256i *)out + 0, q0);
  _mm256_storeu_si256((__m256i *)out + 1, q1);
  _mm256_storeu_si256((__m256i *)out + 2, q2);
  _mm256_storeu_si256((__m256i *)out + 3, q3);
}

OD_SIMD_INLINE void od_transpose_good_name_here_4x4(__m128i *q0, __m128i *q1,
                                                    __m128i *q2, __m128i *q3) {
  __m128i a;
  __m128i b;
  /* Input:
     q0: q31 q21 q11 q01 | q30 q20 q10 q00
     q1: q33 q23 q13 q03 | q32 q22 q12 q02 */
  /* a: q32 q30 q22 q20 q12 q10 q02 q00 */
  a = _mm_unpacklo_epi16(*q0, *q2);
  /* b: q33 q31 q23 q21 q13 q11 q03 q01 */
  b = _mm_unpackhi_epi16(*q0, *q2);
  /* q0: q13 q12 q11 q10 q03 q02 q01 q00 */
  *q0 = _mm_unpacklo_epi16(a, b);
  /* q2: q33 q32 q31 q30 q23 q22 q21 q20 */
  *q2 = _mm_unpackhi_epi16(a, b);
  *q1 = _mm_unpackhi_epi64(*q0, *q0);
  *q3 = _mm_unpackhi_epi64(*q2, *q2);
}

OD_SIMD_INLINE void od_transpose_good_name_here_4x8(__m128i *q0, __m128i *q1,
                                                    __m128i *q2, __m128i *q3) {
  __m128i a;
  __m128i b;
  __m128i c;
  __m128i d;
  __m128i e;
  __m128i f;
  __m128i g;
  __m128i h;
  /* Input:
     q0: q31 q21 q11 q01 q30 q20 q10 q00
     q1: q33 q23 q13 q03 q32 q22 q12 q02
     q2: q35 q25 q15 q05 q34 q24 q14 q04
     q3: q37 q27 q17 q07 q36 q26 q16 q06 */
  /* a: q32 q30 q22 q20 q12 q10 q02 q00 */
  a = _mm_unpacklo_epi16(*q0, *q1);
  /* b: q33 q31 q23 q21 q13 q11 q03 q01 */
  b = _mm_unpackhi_epi16(*q0, *q1);
  /* a: q36 q34 q26 q24 q16 q14 q06 q04 */
  c = _mm_unpacklo_epi16(*q2, *q3);
  /* b: q37 q35 q27 q25 q17 q15 q07 q05 */
  d = _mm_unpackhi_epi16(*q2, *q3);
  /* e: q13 q12 q11 q10 q03 q02 q01 q00 */
  e = _mm_unpacklo_epi16(a, b);
  /* f: q33 q32 q31 q30 q23 q22 q21 q20 */
  f = _mm_unpackhi_epi16(a, b);
  /* g: q17 q16 q15 q14 q07 q06 q05 q04 */
  g = _mm_unpacklo_epi16(c, d);
  /* h: q37 q36 q35 q34 q27 q26 q25 q24 */
  h = _mm_unpackhi_epi16(c, d);
  *q0 = _mm_unpacklo_epi64(e, g);
  *q1 = _mm_unpackhi_epi64(e, g);
  *q2 = _mm_unpacklo_epi64(f, h);
  *q3 = _mm_unpackhi_epi64(f, h);
}

OD_SIMD_INLINE void od_transpose_good_name_here_8x4(__m128i *q0, __m128i *q1,
                                                    __m128i *q2, __m128i *q3,
                                                    __m128i *q4, __m128i *q5,
                                                    __m128i *q6, __m128i *q7) {
  od_transpose8x4(q0, q2, q4, q6);
  /* q0: q31 q21 q11 q01 | q30 q20 q10 q00
     q1: q33 q23 q13 q03 | q32 q22 q12 q02
     q2: q35 q25 q15 q05 | q34 q24 q14 q04
     q3: q37 q27 q17 q07 | q36 q26 q16 q06
  */
  *q1 = _mm_unpackhi_epi64(*q0, *q0);
  *q3 = _mm_unpackhi_epi64(*q2, *q2);
  *q5 = _mm_unpackhi_epi64(*q4, *q4);
  *q7 = _mm_unpackhi_epi64(*q6, *q6);
}

OD_SIMD_INLINE void od_transpose_16x8_unpack(__m256i *s0, __m256i *s1,
                                             __m256i *s2, __m256i *s3,
                                             __m256i *s4, __m256i *s5,
                                             __m256i *s6, __m256i *s7,
                                             __m256i *s8, __m256i *s9,
                                             __m256i *sa, __m256i *sb,
                                             __m256i *sc, __m256i *sd,
                                             __m256i *se, __m256i *sf) {
  const __m256i zero = _mm256_setzero_si256();
  __m256i a;
  __m256i b;
  __m256i c;
  __m256i d;
  __m256i e;
  __m256i f;
  __m256i g;
  __m256i h;
  __m256i s;
  __m256i t;
  __m256i u;
  __m256i v;
  __m256i w;
  __m256i x;
  __m256i y;
  __m256i z;
  __m256i sign0;
  __m256i sign1;
  __m256i sign2;
  __m256i sign3;
  __m256i sign4;
  __m256i sign5;
  __m256i sign6;
  __m256i sign7;
  /* Input:
     s0: rF0 rE0 rD0 rC0 rB0 rA0 r90 r80 r70 r60 r50 r40 r30 r20 r10 r00
     s1: rF1 rE1 rD1 rC1 rB1 rA1 r91 r81 r71 r61 r51 r41 r31 r21 r11 r01
     s2: rF2 rE2 rD2 rC2 rB2 rA2 r92 r82 r72 r62 r52 r42 r32 r22 r12 r02
     s3: rF3 rE3 rD3 rC3 rB3 rA3 r93 r83 r73 r63 r53 r43 r33 r23 r13 r03
     s4: rF4 rE4 rD4 rC4 rB4 rA4 r94 r84 r74 r64 r54 r44 r34 r24 r14 r04
     s5: rF5 rE5 rD5 rC5 rB5 rA5 r95 r85 r75 r65 r55 r45 r35 r25 r15 r05
     s6: rF6 rE6 rD6 rC6 rB6 rA6 r96 r86 r76 r66 r56 r46 r36 r26 r16 r06
     s7: rF7 rE7 rD7 rC7 rB7 rA7 r97 r87 r77 r67 r57 r47 r37 r27 r17 r07
  */
  /*     rB1 rB0 rA1 rA0 r91 r90 r81 r80|r31 r30 r21 r20 r11 r10 r01 r00 */
  a = _mm256_unpacklo_epi16(*s0, *s1);
  /*     rF1 rF0 rE1 rE0 rD1 rD0 rC1 rC0|r71 r70 r61 r60 r51 r50 r41 r40 */
  b = _mm256_unpackhi_epi16(*s0, *s1);
  /*     rB3 rB2 rA3 rA2 r93 r92 r83 r82|r33 r32 r23 r22 r13 r12 r03 r02 */
  c = _mm256_unpacklo_epi16(*s2, *s3);
  /*     rF3 rF2 rE3 rE2 rD3 rD2 rC3 rC2|r73 r72 r63 r62 r53 r52 r43 r42 */
  d = _mm256_unpackhi_epi16(*s2, *s3);
  /*     rB5 rB4 rA5 rA4 r95 r94 r85 r84|r35 r34 r25 r24 r15 r14 r05 r04 */
  e = _mm256_unpacklo_epi16(*s4, *s5);
  /*     rF5 rF4 rE5 rE4 rD5 rD4 rC5 rC4|r75 r74 r65 r64 r55 r54 r45 r44 */
  f = _mm256_unpackhi_epi16(*s4, *s5);
  /*     rB7 rB6 rA7 rA6 r97 r96 r87 r86|r37 r36 r27 r26 r17 r16 r07 r06 */
  g = _mm256_unpacklo_epi16(*s6, *s7);
  /*     rF7 rF6 rE7 rE6 rD7 rD6 rC7 rC6|r77 r76 r67 r66 r57 r56 r47 r46 */
  h = _mm256_unpackhi_epi16(*s6, *s7);

  /*     r93 r92 r91 r90 r83 r82 r81 r80|r13 r12 r11 r10 r03 r02 r01 r00 */
  s = _mm256_unpacklo_epi32(a, c);
  /*     rB3 rB2 rB1 rB0 rA3 rA2 rA1 rA0|r33 r32 r31 r30 r23 r22 r21 r20 */
  t = _mm256_unpackhi_epi32(a, c);
  /*     rD3 rD2 rD1 rD0 rC3 rC2 rC1 rC0|r53 r52 r51 r50 r43 r42 r41 r40 */
  u = _mm256_unpacklo_epi32(b, d);
  /*     rF3 rF2 rF1 rF0 rE3 rE2 rE1 rE0|r73 r72 r71 r70 r63 r62 r61 r60 */
  v = _mm256_unpackhi_epi32(b, d);
  /*     r97 r96 r95 r94 r87 r86 r85 r84|r17 r16 r15 r14 r07 r06 r05 r04 */
  w = _mm256_unpacklo_epi32(e, g);
  /*     rB7 rB6 rB5 rB4 rA7 rA6 rA5 rA4|r37 r36 r35 r34 r27 r26 r25 r24 */
  x = _mm256_unpackhi_epi32(e, g);
  /*     rD7 rD6 rD5 rD4 rE7 rE6 rE5 rE4|r57 r56 r55 r54 r47 r46 r45 r44 */
  y = _mm256_unpacklo_epi32(f, h);
  /*     rF7 rF6 rF5 rF4 rE7 rE6 rE5 rE4|r77 r76 r75 r74 r67 r66 r65 r64 */
  z = _mm256_unpackhi_epi32(f, h);

  /*     r17 r16 r15 r14 r07 r06 r05 r04|r13 r12 r11 r10 r03 r02 r01 r00 */
  a = _mm256_permute2x128_si256(s, w, 0 | (2 << 4));
  b = _mm256_permute2x128_si256(s, w, 1 | (3 << 4));
  c = _mm256_permute2x128_si256(t, x, 0 | (2 << 4));
  d = _mm256_permute2x128_si256(t, x, 1 | (3 << 4));
  e = _mm256_permute2x128_si256(u, y, 0 | (2 << 4));
  f = _mm256_permute2x128_si256(u, y, 1 | (3 << 4));
  g = _mm256_permute2x128_si256(v, z, 0 | (2 << 4));
  h = _mm256_permute2x128_si256(v, z, 1 | (3 << 4));

  sign0 = _mm256_cmpgt_epi16(zero, a);
  sign1 = _mm256_cmpgt_epi16(zero, b);
  sign2 = _mm256_cmpgt_epi16(zero, c);
  sign3 = _mm256_cmpgt_epi16(zero, d);
  sign4 = _mm256_cmpgt_epi16(zero, e);
  sign5 = _mm256_cmpgt_epi16(zero, f);
  sign6 = _mm256_cmpgt_epi16(zero, g);
  sign7 = _mm256_cmpgt_epi16(zero, h);
  *s0 = _mm256_unpacklo_epi16(a, sign0);
  *s1 = _mm256_unpackhi_epi16(a, sign0);
  *s8 = _mm256_unpacklo_epi16(b, sign1);
  *s9 = _mm256_unpackhi_epi16(b, sign1);
  *s2 = _mm256_unpacklo_epi16(c, sign2);
  *s3 = _mm256_unpackhi_epi16(c, sign2);
  *sa = _mm256_unpacklo_epi16(d, sign3);
  *sb = _mm256_unpackhi_epi16(d, sign3);
  *s4 = _mm256_unpacklo_epi16(e, sign4);
  *s5 = _mm256_unpackhi_epi16(e, sign4);
  *sc = _mm256_unpacklo_epi16(f, sign5);
  *sd = _mm256_unpackhi_epi16(f, sign5);
  *s6 = _mm256_unpacklo_epi16(g, sign6);
  *s7 = _mm256_unpackhi_epi16(g, sign6);
  *se = _mm256_unpacklo_epi16(h, sign7);
  *sf = _mm256_unpackhi_epi16(h, sign7);
}

OD_SIMD_INLINE void od_8x2_cvtepi16_epi32(__m128i *q0, __m128i *q1,
                                          __m128i *extra) {
  const __m128i zero = _mm_setzero_si128();
  __m128i sign0;
  __m128i sign1;
  sign0 = _mm_cmpgt_epi16(zero, *q0);
  sign1 = _mm_cmpgt_epi16(zero, *q1);
  extra[0] = _mm_unpackhi_epi16(*q0, sign0);
  extra[1] = _mm_unpackhi_epi16(*q1, sign1);
  *q0 = _mm_unpacklo_epi16(*q0, sign0);
  *q1 = _mm_unpacklo_epi16(*q1, sign1);
}

static void od_store_buffer_8x4_epi16(int16_t *out, int output_stride,
                                      __m128i q0, __m128i q1, __m128i q2,
                                      __m128i q3) {
  _mm_storeu_si128((__m128i *)(out + 0 * output_stride), q0);
  _mm_storeu_si128((__m128i *)(out + 1 * output_stride), q1);
  _mm_storeu_si128((__m128i *)(out + 2 * output_stride), q2);
  _mm_storeu_si128((__m128i *)(out + 3 * output_stride), q3);
}

static void od_col_fidtx_avx2(int16_t *out, int rows, int cols,
                              const int16_t *in, int input_stride, int bd) {
  __m128i q0;
  __m128i q1;
  __m128i q2;
  __m128i q3;
  if (cols <= 4) {
    int r;
    for (r = 0; r < rows; r += 4) {
      od_load_buffer_hbd_4x4_epi16(&q0, &q1, &q2, &q3, in + r * input_stride,
                                   input_stride, bd);
      od_store_halves_buffer_4x4_epi16(out + r * cols, q0, q1, q2, q3);
    }
  } else {
    int r;
    int c;
    for (r = 0; r < rows; r += 4) {
      for (c = 0; c < cols; c += 8) {
        od_load_buffer_hbd_8x4_epi16(
            &q0, &q1, &q2, &q3, in + r * input_stride + c, input_stride, bd);
        od_store_buffer_8x4_epi16(out + r * cols + c, cols, q0, q1, q2, q3);
      }
    }
  }
}

static void od_row_fidtx_avx2(tran_low_t *output_coeffs, int coeffs,
                              const int16_t *in) {
  int c;
  /* The number of rows and number of columns are both multiples of 4, so the
     total number of coefficients should be a multiple of 16. */
  assert(!(coeffs & 0xF));
  /* TODO(any): Use AVX2 for larger block sizes. */
  for (c = 0; c < coeffs; c += 16) {
    __m128i ext[2];
    __m128i q0;
    __m128i q1;
    od_load_buffer_compact_4x4_epi16(&q0, &q1, in + c);
    od_8x2_cvtepi16_epi32(&q0, &q1, ext);
    od_store_buffer_4x4_epi32(output_coeffs + c, q0, ext[0], q1, ext[1]);
  }
}

typedef void (*od_tx4_kernel8_epi16)(__m128i *q0, __m128i *q2, __m128i *q1,
                                     __m128i *q3);

static void od_col_tx4_avx2(int16_t *out, int cols, const int16_t *in,
                            int input_stride, int bd,
                            od_tx4_kernel8_epi16 kernel8) {
  __m128i q0;
  __m128i q1;
  __m128i q2;
  __m128i q3;
  if (cols <= 4) {
    od_load_buffer_hbd_4x4_epi16(&q0, &q1, &q2, &q3, in, input_stride, bd);
    kernel8(&q0, &q1, &q2, &q3);
    od_store_halves_buffer_4x4_epi16(out, q0, q2, q1, q3);
  } else {
    int c;
    for (c = 0; c < cols; c += 8) {
      od_load_buffer_hbd_8x4_epi16(
          &q0, &q1, &q2, &q3, in + c, input_stride, bd);
      kernel8(&q0, &q1, &q2, &q3);
      od_store_buffer_8x4_epi16(out + c, cols, q0, q2, q1, q3);
    }
  }
}

static void od_row_tx4_avx2(tran_low_t *output_coeffs, int rows,
                            const int16_t *in, od_tx4_kernel8_epi16 kernel8) {
  __m128i q0;
  __m128i q1;
  __m128i q2;
  __m128i q3;
  if (rows <= 4) {
    od_load_buffer_compact_4x4_epi16(&q0, &q2, in);
    od_transpose_good_name_here_4x4(&q0, &q1, &q2, &q3);
    kernel8(&q0, &q1, &q2, &q3);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose_unpack4x4(&q0, &q2, &q1, &q3);
    od_store_buffer_4x4_epi32(output_coeffs, q0, q2, q1, q3);
  } else {
    /* Higher row counts require 32-bit precision. */
    assert(rows <= 16);
    int r;
    for (r = 0; r < rows; r += 8) {
      __m128i ext[4];
      od_load_buffer_8x4_epi16(&q0, &q1, &q2, &q3, in + r * 4, 8);
      od_transpose_good_name_here_4x8(&q0, &q1, &q2, &q3);
      kernel8(&q0, &q1, &q2, &q3);
      /*TODO(any): Merge this transpose with coefficient scanning.*/
      od_transpose8x4(&q0, &q2, &q1, &q3);
      od_8x4_cvtepi16_epi32(&q0, &q1, &q2, &q3, ext);
      od_store_buffer_8x4_epi32(output_coeffs + r * 4, q0, ext[0], q2, ext[2], q1, ext[1],
                                q3, ext[3]);
    }
  }
}

static void od_col_fdct4_avx2(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd) {
  od_col_tx4_avx2(out, cols, in, input_stride, bd, od_fdct_4_kernel8_epi16);
}

static void od_row_fdct4_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx4_avx2(output_coeffs, rows, in, od_fdct_4_kernel8_epi16);
}

static void od_col_fdst4_avx2(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd) {
  od_col_tx4_avx2(out, cols, in, input_stride, bd, od_fdst_vii_4_kernel8_epi16);
}

static void od_row_fdst4_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx4_avx2(output_coeffs, rows, in, od_fdst_vii_4_kernel8_epi16);
}

static void od_col_flip_fdst4_avx2(int16_t *out, int cols, const int16_t *in,
                                   int input_stride, int bd) {
  od_col_tx4_avx2(out, cols, in, input_stride, bd,
                  od_flip_fdst_vii_4_kernel8_epi16);
}

static void od_row_flip_fdst4_avx2(tran_low_t *output_coeffs, int rows,
                                   const int16_t *in) {
  od_row_tx4_avx2(output_coeffs, rows, in, od_flip_fdst_vii_4_kernel8_epi16);
}

static void od_col_fidtx4_avx2(int16_t *out, int cols, const int16_t *in,
                               int input_stride, int bd) {
  od_col_fidtx_avx2(out, 4, cols, in, input_stride, bd);
}

static void od_row_fidtx4_avx2(tran_low_t *output_coeffs, int rows,
                               const int16_t *in) {
  od_row_fidtx_avx2(output_coeffs, rows * 4, in);
}

typedef void (*od_tx8_kernel8_epi16)(__m128i *r0, __m128i *r4, __m128i *r2,
                                     __m128i *r6, __m128i *r1, __m128i *r5,
                                     __m128i *r3, __m128i *r7);

static void od_col_tx8_avx2(int16_t *out, int cols, const int16_t *in,
                            int input_stride, int bd,
                            od_tx8_kernel8_epi16 kernel8) {
  __m128i r0;
  __m128i r1;
  __m128i r2;
  __m128i r3;
  __m128i r4;
  __m128i r5;
  __m128i r6;
  __m128i r7;
  if (cols <= 4) {
    od_load_buffer_hbd_4x4_epi16(&r0, &r1, &r2, &r3, in, input_stride, bd);
    od_load_buffer_hbd_4x4_epi16(&r4, &r5, &r6, &r7, in + 4 * input_stride,
                                 input_stride, bd);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    od_store_halves_buffer_4x4_epi16(out, r0, r4, r2, r6);
    od_store_halves_buffer_4x4_epi16(out + 16, r1, r5, r3, r7);
  } else {
    int c;
    for (c = 0; c < cols; c += 8) {
      od_load_buffer_hbd_8x4_epi16(
          &r0, &r1, &r2, &r3, in + c, input_stride, bd);
      od_load_buffer_hbd_8x4_epi16(
          &r4, &r5, &r6, &r7, in + 4 * input_stride + c, input_stride, bd);
      kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
      od_store_buffer_8x4_epi16(out + c, cols, r0, r4, r2, r6);
      od_store_buffer_8x4_epi16(out + 4 * cols + c, cols, r1, r5, r3, r7);
    }
  }
  //TODO: implement x16
}

static void od_row_tx8_avx2(tran_low_t *output_coeffs, int rows,
                            const int16_t *in, od_tx8_kernel8_epi16 kernel8) {
  __m128i r0;
  __m128i r1;
  __m128i r2;
  __m128i r3;
  __m128i r4;
  __m128i r5;
  __m128i r6;
  __m128i r7;
  if (rows <= 4) {
    od_load_buffer_8x4_epi16(&r0, &r2, &r4, &r6, in, 8);
    od_transpose_good_name_here_8x4(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose_unpack4x4(&r0, &r4, &r2, &r6);
    od_transpose_unpack4x4(&r1, &r5, &r3, &r7);
    od_store_buffer_4x4_epi32(output_coeffs, r0, r1, r4, r5);
    od_store_buffer_4x4_epi32(output_coeffs + 16, r2, r3, r6, r7);
  } else {
    int r;
    for (r = 0; r < rows; r += 8) {
      __m128i ext[8];
      od_load_buffer_8x4_epi16(&r0, &r1, &r2, &r3, in + r * 8, 8);
      od_load_buffer_8x4_epi16(&r4, &r5, &r6, &r7, in + 32 + r * 8, 8);
      od_transpose8x8_epi16(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
      kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
      od_transpose8x8_epi16(&r0, &r4, &r2, &r6, &r1, &r5, &r3, &r7);
      /*TODO(any): Merge this transpose with coefficient scanning.*/
      od_8x4_cvtepi16_epi32(&r0, &r1, &r2, &r3, ext);
      od_8x4_cvtepi16_epi32(&r4, &r5, &r6, &r7, ext + 4);
      od_store_buffer_8x4_epi32(output_coeffs + r * 8, r0, ext[0], r4, ext[4], r2, ext[2],
                                r6, ext[6]);
      od_store_buffer_8x4_epi32(output_coeffs + r * 8 + 32, r1, ext[1], r5, ext[5], r3,
                                ext[3], r7, ext[7]);
    }
  }
  // TODO: Must implement x16 because of overflow
}

static void od_row_fdct8_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx8_avx2(output_coeffs, rows, in, od_fdct_8_kernel8_epi16);
}

static void od_col_fdct8_avx2(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd) {
  od_col_tx8_avx2(out, cols, in, input_stride, bd, od_fdct_8_kernel8_epi16);
}

static void od_col_fdst8_avx2(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd) {
  od_col_tx8_avx2(out, cols, in, input_stride, bd, od_fdst_8_kernel8_epi16);
}

static void od_row_fdst8_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx8_avx2(output_coeffs, rows, in, od_fdst_8_kernel8_epi16);
}

static void od_col_flip_fdst8_avx2(int16_t *out, int cols, const int16_t *in,
                                   int input_stride, int bd) {
  od_col_tx8_avx2(out, cols, in, input_stride, bd,
                  od_flip_fdst_8_kernel8_epi16);
}

static void od_row_flip_fdst8_avx2(tran_low_t *output_coeffs, int rows,
                                   const int16_t *in) {
  od_row_tx8_avx2(output_coeffs, rows, in, od_flip_fdst_8_kernel8_epi16);
}

static void od_col_fidtx8_avx2(int16_t *out, int cols, const int16_t *in,
                               int input_stride, int bd) {
  od_col_fidtx_avx2(out, 8, cols, in, input_stride, bd);
}

static void od_row_fidtx8_avx2(tran_low_t *output_coeffs, int rows,
                               const int16_t *in) {
  od_row_fidtx_avx2(output_coeffs, rows * 8, in);
}

typedef void (*od_tx16_kernel8_epi16)(__m128i *s0, __m128i *s1, __m128i *s2,
                                      __m128i *s3, __m128i *s4, __m128i *s5,
                                      __m128i *s6, __m128i *s7, __m128i *s8,
                                      __m128i *s9, __m128i *sa, __m128i *sb,
                                      __m128i *sc, __m128i *sd, __m128i *se,
                                      __m128i *sf);

typedef void (*od_tx16_mm256_kernel)(__m256i *s0, __m256i *s1, __m256i *s2,
                                     __m256i *s3, __m256i *s4, __m256i *s5,
                                     __m256i *s6, __m256i *s7, __m256i *s8,
                                     __m256i *s9, __m256i *sa, __m256i *sb,
                                     __m256i *sc, __m256i *sd, __m256i *se,
                                     __m256i *sf);

static void od_row_tx16_avx2(tran_low_t *output_coeffs, int rows,
                             const int16_t *in,
                             od_tx16_kernel8_epi16 kernel8_epi16,
                             od_tx16_mm256_kernel kernel8_epi32) {
  if (rows <= 4) {
    __m128i s0;
    __m128i s1;
    __m128i s2;
    __m128i s3;
    __m128i s4;
    __m128i s5;
    __m128i s6;
    __m128i s7;
    __m128i s8;
    __m128i s9;
    __m128i sa;
    __m128i sb;
    __m128i sc;
    __m128i sd;
    __m128i se;
    __m128i sf;
    od_load_buffer_8x4_epi16(&s0, &s8, &s2, &sa, in, 8);
    od_load_buffer_8x4_epi16(&s4, &sc, &s6, &se, in + 32, 8);
    od_transpose_good_name_here_8x4(&s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7);
    od_transpose_good_name_here_8x4(&s8, &s9, &sa, &sb, &sc, &sd, &se, &sf);
    kernel8_epi16(&s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &sa, &sb,
            &sc, &sd, &se, &sf);
    od_transpose_unpack4x4(&s0, &s8, &s4, &sc);
    od_transpose_unpack4x4(&s2, &sa, &s6, &se);
    od_transpose_unpack4x4(&s1, &s9, &s5, &sd);
    od_transpose_unpack4x4(&s3, &sb, &s7, &sf);
    od_store_buffer_4x4_epi32(output_coeffs, s0, s2, s1, s3);
    od_store_buffer_4x4_epi32(output_coeffs + 16, s8, sa, s9, sb);
    od_store_buffer_4x4_epi32(output_coeffs + 32, s4, s6, s5, s7);
    od_store_buffer_4x4_epi32(output_coeffs + 48, sc, se, sd, sf);
    return;
  }
  {
    int r;
    for (r = 0; r < rows; r += 8) {
      __m256i ss0;
      __m256i ss1;
      __m256i ss2;
      __m256i ss3;
      __m256i ss4;
      __m256i ss5;
      __m256i ss6;
      __m256i ss7;
      __m256i ss8;
      __m256i ss9;
      __m256i ssa;
      __m256i ssb;
      __m256i ssc;
      __m256i ssd;
      __m256i sse;
      __m256i ssf;
      od_load_buffer_16x4_epi16(&ss0, &ss1, &ss2, &ss3, in + r * 16, 16);
      od_load_buffer_16x4_epi16(&ss4, &ss5, &ss6, &ss7, in + r * 16 + 64, 16);
      od_transpose_16x8_unpack(&ss0, &ss1, &ss2, &ss3, &ss4, &ss5, &ss6, &ss7, &ss8, &ss9,
                    &ssa, &ssb, &ssc, &ssd, &sse, &ssf);
      kernel8_epi32(&ss0, &ss1, &ss2, &ss3, &ss4, &ss5, &ss6, &ss7, &ss8, &ss9,
                    &ssa, &ssb, &ssc, &ssd, &sse, &ssf);
      /*TODO(any): Merge this transpose with coefficient scanning.*/
      od_transpose8x8_epi32(&ss0, &ss8, &ss4, &ssc, &ss2, &ssa, &ss6, &sse);
      od_transpose8x8_epi32(&ss1, &ss9, &ss5, &ssd, &ss3, &ssb, &ss7, &ssf);
      od_store_buffer_4x8_epi32(output_coeffs + r * 16, ss0, ss1, ss8, ss9);
      od_store_buffer_4x8_epi32(output_coeffs + r * 16 + 32, ss4, ss5, ssc, ssd);
      od_store_buffer_4x8_epi32(output_coeffs + r * 16 + 64, ss2, ss3, ssa, ssb);
      od_store_buffer_4x8_epi32(output_coeffs + r * 16 + 96, ss6, ss7, sse, ssf);
    }
  }
}

static void od_row_fdct16_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx16_avx2(output_coeffs, rows, in, od_fdct_16_kernel8_epi16, od_fdct_16_kernel8_epi32);
}

static void od_row_fdst16_avx2(tran_low_t *output_coeffs, int rows,
                              const int16_t *in) {
  od_row_tx16_avx2(output_coeffs, rows, in, od_fdst_16_kernel8_epi16, od_fdst_16_kernel8_epi32);
}

static void od_row_fidtx16_avx2(tran_low_t *output_coeffs, int rows,
                                const int16_t *in) {
  od_row_fidtx_avx2(output_coeffs, rows * 16, in);
}

static void od_col_fidtx16_avx2(int16_t *out, int cols, const int16_t *in,
                               int input_stride, int bd) {
  od_col_fidtx_avx2(out, 16, cols, in, input_stride, bd);
}

typedef void (*daala_col_ftx)(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd);
typedef void (*daala_row_ftx)(tran_low_t *output_coeffs, int rows,
                              const int16_t *in);

static const daala_col_ftx TX_COL_MAP[TX_SIZES][TX_TYPES] = {
  // 4-point transforms
  { od_col_fdct4_avx2, od_col_fdst4_avx2, od_col_flip_fdst4_avx2,
    od_col_fidtx4_avx2 },
  // 8-point transforms
  { od_col_fdct8_avx2, od_col_fdst8_avx2, od_col_flip_fdst8_avx2,
    od_col_fidtx8_avx2 },
  // 16-point transforms
  { NULL, NULL, NULL, od_col_fidtx16_avx2 },
  // 32-point transforms
  { NULL, NULL, NULL, NULL },
#if CONFIG_TX64X64
  // 64-point transforms
  { NULL, NULL, NULL, NULL },
#endif
};

static const daala_row_ftx TX_ROW_MAP[TX_SIZES][TX_TYPES] = {
  // 4-point transforms
  { od_row_fdct4_avx2, od_row_fdst4_avx2, od_row_flip_fdst4_avx2,
    od_row_fidtx4_avx2 },
  // 8-point transforms
  { od_row_fdct8_avx2, od_row_fdst8_avx2, od_row_flip_fdst8_avx2,
    od_row_fidtx8_avx2 },
  // 16-point transforms
  { od_row_fdct16_avx2, od_row_fdst16_avx2, NULL,
    od_row_fidtx16_avx2 },
  // 32-point transforms
  { NULL, NULL, NULL, NULL },
#if CONFIG_TX64X64
  // 64-point transforms
  { NULL, NULL, NULL, NULL },
#endif
};

/* Define this to verify the SIMD against the C versions of the transforms.
   This is intended to be replaced by real unit tests in the future. */
#undef DAALA_TX_VERIFY_SIMD

void daala_fwd_txfm_avx2(const int16_t *input_pixels, tran_low_t *output_coeffs,
                         int input_stride, TxfmParam *txfm_param) {
  const TX_SIZE tx_size = txfm_param->tx_size;
  const TX_TYPE tx_type = txfm_param->tx_type;
  assert(tx_size <= TX_SIZES_ALL);
  assert(tx_type <= TX_TYPES);

  if (txfm_param->lossless) {
    daala_fwd_txfm_c(input_pixels, output_coeffs, input_stride, txfm_param);
  } else {
    // General TX case
    assert(sizeof(tran_low_t) == sizeof(od_coeff));
    assert(sizeof(tran_low_t) >= 4);

    // Hook into existing map translation infrastructure to select
    // appropriate TX functions
    const TX_SIZE col_idx = txsize_vert_map[tx_size];
    const TX_SIZE row_idx = txsize_horz_map[tx_size];
    assert(col_idx <= TX_SIZES);
    assert(row_idx <= TX_SIZES);
    assert(vtx_tab[tx_type] <= (int)TX_TYPES_1D);
    assert(htx_tab[tx_type] <= (int)TX_TYPES_1D);
    daala_col_ftx col_tx = TX_COL_MAP[col_idx][vtx_tab[tx_type]];
    daala_row_ftx row_tx = TX_ROW_MAP[row_idx][htx_tab[tx_type]];

    int16_t tmpsq[MAX_TX_SQUARE];
    if (col_tx == NULL || row_tx == NULL) {
      daala_fwd_txfm_c(input_pixels, output_coeffs, input_stride, txfm_param);
    } else {
      const int cols = tx_size_wide[tx_size];
      const int rows = tx_size_high[tx_size];
      // Forward-transform columns
      col_tx(tmpsq, cols, input_pixels, input_stride, txfm_param->bd);
      // Forward-transform columns and sum with destination
      row_tx(output_coeffs, rows, tmpsq);
#if defined(DAALA_TX_VERIFY_SIMD)
      tran_low_t out_check_buf[MAX_TX_SQUARE];
      daala_fwd_txfm_c(input_pixels, out_check_buf, input_stride, txfm_param);
      {
        int r;
        for (r = 0; r < rows; r++) {
          if (memcmp(out_check_buf + r * cols, output_coeffs + r * cols,
                     cols * sizeof(*out_check_buf))) {
            fprintf(stderr, "%s(%i): Forward %ix%i %i_%i TX SIMD mismatch.\n",
                    __FILE__, __LINE__, rows, cols, vtx_tab[tx_type],
                    htx_tab[tx_type]);
            assert(0);
            exit(EXIT_FAILURE);
          }
        }
      }
#endif
    }
  }
}

#endif
