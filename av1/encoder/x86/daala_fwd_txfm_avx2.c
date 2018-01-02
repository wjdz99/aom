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

OD_SIMD_INLINE void od_8x2_cvtepi16_epi32(__m128i *q0, __m128i *q1, __m128i *extra) {
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
                              const int16_t *in, int input_stride,
                              int bd) {
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
  }
  else {
    int r;
    int c;
    for (r = 0; r < rows; r += 4) {
      for (c = 0; c < cols; c += 8) {
        od_load_buffer_hbd_8x4_epi16(&q0, &q1, &q2, &q3, in + r * input_stride + c, input_stride, bd);
        od_store_buffer_8x4_epi16(out + r * cols + c,
                                  cols, q0, q1, q2, q3);
      }
    }
  }
}

static void od_row_fidtx_avx2(tran_low_t *output_coeffs, int coeffs, const int16_t *in) {
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
    /* Higher col count unimplemented.*/
    assert(cols <= 8);
    od_load_buffer_hbd_8x4_epi16(&q0, &q1, &q2, &q3, in, input_stride, bd);
    kernel8(&q0, &q1, &q2, &q3);
    od_store_buffer_4x8_epi16(out, q0, q2, q1, q3);
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
    __m128i ext[4];
    /* Higher row count unimplemented.*/
    assert(rows <= 8);
    od_load_buffer_8x4_epi16(&q0, &q1, &q2, &q3, in, rows);
    od_transpose_good_name_here_4x8(&q0, &q1, &q2, &q3);
    kernel8(&q0, &q1, &q2, &q3);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose8x4(&q0, &q2, &q1, &q3);
    od_8x4_cvtepi16_epi32(&q0, &q1, &q2, &q3, ext);
    od_store_buffer_8x4_epi32(output_coeffs, q0, ext[0], q2, ext[2], q1, ext[1],
                              q3, ext[3]);
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
    /* Higher col count unimplemented.*/
    assert(cols <= 8);
    od_load_buffer_hbd_8x4_epi16(&r0, &r1, &r2, &r3, in, input_stride, bd);
    od_load_buffer_hbd_8x4_epi16(&r4, &r5, &r6, &r7, in + 4 * input_stride,
                                 input_stride, bd);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    od_store_buffer_4x8_epi16(out, r0, r4, r2, r6);
    od_store_buffer_4x8_epi16(out + 32, r1, r5, r3, r7);
  }
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
    /* Higher row count unimplemented.*/
    assert(rows <= 8);
    __m128i ext[8];
    od_load_buffer_8x4_epi16(&r0, &r1, &r2, &r3, in, rows);
    od_load_buffer_8x4_epi16(&r4, &r5, &r6, &r7, in + 32, rows);
    od_transpose8x8_epi16(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose8x8_epi16(&r0, &r4, &r2, &r6, &r1, &r5, &r3, &r7);
    od_8x4_cvtepi16_epi32(&r0, &r1, &r2, &r3, ext);
    od_8x4_cvtepi16_epi32(&r4, &r5, &r6, &r7, ext + 4);
    od_store_buffer_8x4_epi32(output_coeffs, r0, ext[0], r4, ext[4], r2, ext[2],
                              r6, ext[6]);
    od_store_buffer_8x4_epi32(output_coeffs + 32, r1, ext[1], r5, ext[5], r3,
                              ext[3], r7, ext[7]);
  }
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

typedef void (*daala_col_ftx)(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd);
typedef void (*daala_row_ftx)(tran_low_t *output_coeffs, int rows,
                              const int16_t *in);

static const daala_col_ftx TX_COL_MAP[TX_SIZES][TX_TYPES] = {
  // 4-point transforms
  { od_col_fdct4_avx2, od_col_fdst4_avx2, od_col_flip_fdst4_avx2, od_col_fidtx4_avx2 },
  // 8-point transforms
  { od_col_fdct8_avx2, od_col_fdst8_avx2, od_col_flip_fdst8_avx2, od_col_fidtx8_avx2 },
  // 16-point transforms
  { NULL, NULL, NULL, NULL },
  // 32-point transforms
  { NULL, NULL, NULL, NULL },
#if CONFIG_TX64X64
  // 64-point transforms
  { NULL, NULL, NULL, NULL },
#endif
};

static const daala_row_ftx TX_ROW_MAP[TX_SIZES][TX_TYPES] = {
  // 4-point transforms
  { od_row_fdct4_avx2, od_row_fdst4_avx2, od_row_flip_fdst4_avx2, od_row_fidtx4_avx2 },
  // 8-point transforms
  { od_row_fdct8_avx2, od_row_fdst8_avx2, od_row_flip_fdst8_avx2, od_row_fidtx8_avx2 },
  // 16-point transforms
  { NULL, NULL, NULL, NULL },
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
