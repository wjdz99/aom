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

#include <assert.h>
#include "aom_dsp/txfm_common.h"
#include "config/aom_dsp_rtcd.h"

void aom_fdct4x4_c(const int16_t *input, tran_low_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows.
  // We need an intermediate buffer between passes.
  tran_low_t intermediate[4 * 4];
  const tran_low_t *in_low = NULL;
  tran_low_t *out = intermediate;
  // Do the two transform passes
  for (int pass = 0; pass < 2; ++pass) {
    tran_high_t in_high[4];  // canbe16
    tran_high_t step[4];     // canbe16
    tran_low_t temp[4];
    for (int i = 0; i < 4; ++i) {
      // Load inputs.
      if (pass == 0) {
        in_high[0] = input[0 * stride] * 16;
        in_high[1] = input[1 * stride] * 16;
        in_high[2] = input[2 * stride] * 16;
        in_high[3] = input[3 * stride] * 16;
        if (i == 0 && in_high[0]) {
          ++in_high[0];
        }
        ++input;  // Next column
      } else {
        assert(in_low != NULL);
        in_high[0] = in_low[0 * 4];
        in_high[1] = in_low[1 * 4];
        in_high[2] = in_low[2 * 4];
        in_high[3] = in_low[3 * 4];
        ++in_low;  // Next column (which is a transposed row)
      }
      // Transform.
      step[0] = in_high[0] + in_high[3];
      step[1] = in_high[1] + in_high[2];
      step[2] = in_high[1] - in_high[2];
      step[3] = in_high[0] - in_high[3];
      temp[0] = (tran_low_t)fdct_round_shift((step[0] + step[1]) * cospi_16_64);
      temp[2] = (tran_low_t)fdct_round_shift((step[0] - step[1]) * cospi_16_64);
      temp[1] = (tran_low_t)fdct_round_shift(step[2] * cospi_24_64 +
                                             step[3] * cospi_8_64);
      temp[3] = (tran_low_t)fdct_round_shift(-step[2] * cospi_8_64 +
                                             step[3] * cospi_24_64);
      // Only transpose the first pass.
      if (pass == 0) {
        out[0] = temp[0];
        out[1] = temp[1];
        out[2] = temp[2];
        out[3] = temp[3];
        out += 4;
      } else {
        out[0 * 4] = temp[0];
        out[1 * 4] = temp[1];
        out[2 * 4] = temp[2];
        out[3 * 4] = temp[3];
        ++out;
      }
    }
    // Setup in/out for next pass.
    in_low = intermediate;
    out = output;
  }

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      output[j + i * 4] = (output[j + i * 4] + 1) >> 2;
  }
}

void aom_fdct4x4_lp_c(const int16_t *input, int16_t *output, int stride) {
  // The 2D transform is done with two passes which are actually pretty
  // similar. In the first one, we transform the columns and transpose
  // the results. In the second one, we transform the rows.
  // We need an intermediate buffer between passes.
  int16_t intermediate[4 * 4];
  const int16_t *in_low = NULL;
  int16_t *out = intermediate;
  // Do the two transform passes
  for (int pass = 0; pass < 2; ++pass) {
    int32_t in_high[4];  // canbe16
    int32_t step[4];     // canbe16
    int16_t temp[4];
    for (int i = 0; i < 4; ++i) {
      // Load inputs.
      if (pass == 0) {
        in_high[0] = input[0 * stride] * 16;
        in_high[1] = input[1 * stride] * 16;
        in_high[2] = input[2 * stride] * 16;
        in_high[3] = input[3 * stride] * 16;
        ++input;
        if (i == 0 && in_high[0]) {
          ++in_high[0];
        }
      } else {
        assert(in_low != NULL);
        in_high[0] = in_low[0 * 4];
        in_high[1] = in_low[1 * 4];
        in_high[2] = in_low[2 * 4];
        in_high[3] = in_low[3 * 4];
        ++in_low;
      }
      // Transform.
      step[0] = in_high[0] + in_high[3];
      step[1] = in_high[1] + in_high[2];
      step[2] = in_high[1] - in_high[2];
      step[3] = in_high[0] - in_high[3];
      temp[0] = (int16_t)fdct_round_shift((step[0] + step[1]) * cospi_16_64);
      temp[2] = (int16_t)fdct_round_shift((step[0] - step[1]) * cospi_16_64);
      temp[1] = (int16_t)fdct_round_shift(step[2] * cospi_24_64 +
                                          step[3] * cospi_8_64);
      temp[3] = (int16_t)fdct_round_shift(-step[2] * cospi_8_64 +
                                          step[3] * cospi_24_64);
      // Only transpose the first pass.
      if (pass == 0) {
        out[0] = temp[0];
        out[1] = temp[1];
        out[2] = temp[2];
        out[3] = temp[3];
        out += 4;
      } else {
        out[0 * 4] = temp[0];
        out[1 * 4] = temp[1];
        out[2 * 4] = temp[2];
        out[3 * 4] = temp[3];
        ++out;
      }
    }
    // Setup in/out for next pass.
    in_low = intermediate;
    out = output;
  }

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      output[j + i * 4] = (output[j + i * 4] + 1) >> 2;
  }
}