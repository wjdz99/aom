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

#include "aom_dsp/fwd_txfm.h"
#include "aom_dsp/txfm_common.h"
#include "aom_dsp/x86/txfm_common_sse2.h"

// TODO(jingning) The high bit-depth version needs re-work for performance.
// The current SSE2 implementation also causes cross reference to the static
// functions in the C implementation file.
#if DCT_HIGH_BIT_DEPTH
#define ADD_EPI16 _mm_adds_epi16
#define SUB_EPI16 _mm_subs_epi16
#if FDCT32x32_HIGH_PRECISION
void aom_fdct32x32_rows_c(const int16_t *intermediate, tran_low_t *out) {
  int i, j;
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j) temp_in[j] = intermediate[j * 32 + i];
    aom_fdct32(temp_in, temp_out, 0);
    for (j = 0; j < 32; ++j)
      out[j + i * 32] =
          (tran_low_t)((temp_out[j] + 1 + (temp_out[j] < 0)) >> 2);
  }
}
#define HIGH_FDCT32x32_2D_C aom_highbd_fdct32x32_c
#define HIGH_FDCT32x32_2D_ROWS_C aom_fdct32x32_rows_c
#else
void aom_fdct32x32_rd_rows_c(const int16_t *intermediate, tran_low_t *out) {
  int i, j;
  for (i = 0; i < 32; ++i) {
    tran_high_t temp_in[32], temp_out[32];
    for (j = 0; j < 32; ++j) temp_in[j] = intermediate[j * 32 + i];
    aom_fdct32(temp_in, temp_out, 1);
    for (j = 0; j < 32; ++j) out[j + i * 32] = (tran_low_t)temp_out[j];
  }
}
#define HIGH_FDCT32x32_2D_C aom_highbd_fdct32x32_rd_c
#define HIGH_FDCT32x32_2D_ROWS_C aom_fdct32x32_rd_rows_c
#endif  // FDCT32x32_HIGH_PRECISION
#else
#define ADD_EPI16 _mm_add_epi16
#define SUB_EPI16 _mm_sub_epi16
#endif  // DCT_HIGH_BIT_DEPTH

#undef ADD_EPI16
#undef SUB_EPI16
#undef HIGH_FDCT32x32_2D_C
#undef HIGH_FDCT32x32_2D_ROWS_C
