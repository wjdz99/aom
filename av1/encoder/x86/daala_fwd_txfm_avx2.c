#include <tmmintrin.h>
#include <immintrin.h>
#include "./av1_rtcd.h"
#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "av1/common/daala_tx.h"
#include "av1/encoder/daala_fwd_txfm.h"
#include "av1/common/x86/daala_tx_avx2.h"

#if CONFIG_DAALA_TX

typedef void (*od_tx4_kernel8_epi16)(__m128i *q0, __m128i *q1, __m128i *q2,
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
    od_store_buffer_4x4_epi16(out, q0, q2, q1, q3);
  } else {
    /* Higher col count unimplemented.*/
    assert(cols <= 8);
    od_load_buffer_hbd_8x4_epi16(&q0, &q1, &q2, &q3, in, input_stride, bd);
    kernel8(&q0, &q1, &q2, &q3);
    od_store_buffer_8x4_epi16(out, cols, q0, q2, q1, q3);
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
    od_transpose_compact_in_4x4(&q0, &q1, &q2, &q3);
    kernel8(&q0, &q1, &q2, &q3);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose_unpack4x4(&q0, &q2, &q1, &q3);
    od_store_buffer_4x4_epi32(output_coeffs, q0, q2, q1, q3);
  } else {
    __m256i i32[4];
    /* Higher row count unimplemented.*/
    assert(rows <= 8);
    od_load_buffer_8x4_epi16(&q0, &q1, &q2, &q3, in, 8);
    od_transpose_compact_in_4x8(&q0, &q1, &q2, &q3);
    kernel8(&q0, &q1, &q2, &q3);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose_compact_out_8x4(&q0, &q2, &q1, &q3);
    od_4x8_cvtepi16_epi32(q0, q1, q2, q3, i32);
    od_store_buffer_4x8_epi32(output_coeffs, i32[0], i32[2], i32[1], i32[3]);
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

typedef void (*od_tx8_kernel8_epi16)(__m128i *r0, __m128i *r1, __m128i *r2,
                                     __m128i *r3, __m128i *r4, __m128i *r5,
                                     __m128i *r6, __m128i *r7);

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
    od_store_buffer_4x4_epi16(out, r0, r4, r2, r6);
    od_store_buffer_4x4_epi16(out + 16, r1, r5, r3, r7);
  } else {
    /* Higher col count unimplemented.*/
    assert(cols <= 8);
    od_load_buffer_hbd_8x4_epi16(&r0, &r1, &r2, &r3, in, input_stride, bd);
    od_load_buffer_hbd_8x4_epi16(&r4, &r5, &r6, &r7, in + 4 * input_stride,
                                 input_stride, bd);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    od_store_buffer_8x4_epi16(out, cols, r0, r4, r2, r6);
    od_store_buffer_8x4_epi16(out + 4 * cols, cols, r1, r5, r3, r7);
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
    od_transpose_8x4(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose_unpack4x4(&r0, &r4, &r2, &r6);
    od_transpose_unpack4x4(&r1, &r5, &r3, &r7);
    od_store_buffer_4x4_epi32(output_coeffs, r0, r1, r4, r5);
    od_store_buffer_4x4_epi32(output_coeffs + 16, r2, r3, r6, r7);
  } else {
    /* Higher row count unimplemented.*/
    assert(rows <= 8);
    __m256i i32[8];
    od_load_buffer_8x4_epi16(&r0, &r1, &r2, &r3, in, 8);
    od_load_buffer_8x4_epi16(&r4, &r5, &r6, &r7, in + 32, 8);
    od_transpose8x8_epi16(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    kernel8(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
    /*TODO(any): Merge this transpose with coefficient scanning.*/
    od_transpose8x8_epi16(&r0, &r4, &r2, &r6, &r1, &r5, &r3, &r7);
    od_4x8_cvtepi16_epi32(r0, r1, r2, r3, i32);
    od_4x8_cvtepi16_epi32(r4, r5, r6, r7, i32 + 4);
    od_store_buffer_4x8_epi32(output_coeffs, i32[0], i32[4], i32[2], i32[6]);
    od_store_buffer_4x8_epi32(output_coeffs + 32, i32[1], i32[5], i32[3],
                              i32[7]);
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

typedef void (*daala_col_ftx)(int16_t *out, int cols, const int16_t *in,
                              int input_stride, int bd);
typedef void (*daala_row_ftx)(tran_low_t *output_coeffs, int rows,
                              const int16_t *in);

static const daala_col_ftx TX_COL_MAP[TX_SIZES][TX_TYPES] = {
  // 4-point transforms
  { od_col_fdct4_avx2, od_col_fdst4_avx2, od_col_flip_fdst4_avx2, NULL },
  // 8-point transforms
  { od_col_fdct8_avx2, od_col_fdst8_avx2, od_col_flip_fdst8_avx2, NULL },
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
  { od_row_fdct4_avx2, od_row_fdst4_avx2, od_row_flip_fdst4_avx2, NULL },
  // 8-point transforms
  { od_row_fdct8_avx2, od_row_fdst8_avx2, od_row_flip_fdst8_avx2, NULL },
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
