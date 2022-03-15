/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>

#include "av1/encoder/tune_visual_masking.h"

#include "av1/encoder/encodeframe.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/rdopt_utils.h"
#include "av1/encoder/encoder_utils.h"
#include "av1/encoder/extend.h"
#include "av1/encoder/var_based_part.h"

// TODO
static float mul_add(float a, float b, float c) { return a * b + c; }

static float masking_sqrt(float v) {
  static const float kLogOffset = 26.481471032459346f;
  static const float kMul = 211.50759899638012f;
  const float mul_v = kMul * 1e8;
  const float offset_v = kLogOffset;
  return 0.25f * sqrt(mul_add(v, sqrt(mul_v), offset_v));
}

float compute_mask(const float out_val) {
  const float kBase = -0.74174993f;
  const float kMul4 = 3.2353257320940401f;
  const float kMul2 = 12.906028311180409f;
  const float kOffset2 = 305.04035728311436f;
  const float kMul3 = 5.0220313103171232f;
  const float kOffset3 = 2.1925739705298404f;
  const float kOffset4 = 0.25f * kOffset3;
  const float kMul0 = 0.74760422233706747f;
  const float k1 = 1.0f;

  const float v1 = out_val * kMul0 > 1e-3f ? out_val * kMul0 : 1e-3f;
  const float v2 = k1 / (v1 + kOffset2);
  const float v3 = k1 / mul_add(v1, v1, kOffset3);
  const float v4 = k1 / mul_add(v1, v1, kOffset4);
  return (kBase + mul_add(kMul4, v4, mul_add(kMul2, v2, kMul3 * v3))) * 1.5f;
}

float compute_masking(const float pixels[6][6]) {
  float masking = 0.0f;
  for (int iy = 0; iy < 4; iy++) {
    for (int ix = 0; ix < 4; ix++) {
      const float base = 0.25f * (pixels[iy][ix + 1] + pixels[iy + 2][ix + 1] +
                                  pixels[iy + 1][ix] + pixels[iy + 1][ix + 2]);
      float diff = pixels[iy + 1][ix + 1] - base;
      diff *= diff;
      diff = masking_sqrt(diff);
      masking += diff;
    }
  }
  return compute_mask(masking / 16);
}

void av1_set_mb_visual_masking_rdmult_scaling(AV1_COMP *cpi) {
  // Higher value = fewer bits.
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  uint8_t *y_buffer = cpi->source->y_buffer;
  const int y_stride = cpi->source->y_stride;
  const int block_size = BLOCK_4X4;

  const int num_mi_w = mi_size_wide[block_size];
  const int num_mi_h = mi_size_high[block_size];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (mi_params->mi_rows + num_mi_h - 1) / num_mi_h;
  double log_sum = 0.0;

  const int use_hbd = cpi->source->flags & YV12_FLAG_HIGHBITDEPTH;

  float block[6][6];  // 4x4 with 1 border pixels.

  float scale = use_hbd ? 1.0 / ((1 << xd->bd) - 1) : 1.0 / 255;

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;

      for (int iy = -1; iy < 5; iy++) {
        int rowin = (row << 2) + iy;
        rowin = rowin == -1 ? 0 : rowin;
        rowin = rowin == cpi->source->y_height ? rowin - 1 : rowin;
        for (int ix = -1; ix < 5; ix++) {
          int colin = (col << 2) + ix;
          colin = colin == -1 ? 0 : colin;
          colin = colin == cpi->source->y_width ? colin - 1 : colin;
          const uint8_t *ptr = y_buffer + rowin * y_stride + colin;
          if (use_hbd) {
            block[1 + iy][1 + ix] = *CONVERT_TO_SHORTPTR(ptr) * scale;
          } else {
            block[1 + iy][1 + ix] = *ptr * scale;
          }
        }
      }

      float log_masking = -compute_masking(block);
      cpi->visual_masking_rdmult_scaling_factors[index] = exp(log_masking);
      log_sum += log_masking;
    }
  }
  log_sum = exp(log_sum / (double)(num_rows * num_cols));

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      cpi->visual_masking_rdmult_scaling_factors[index] /= log_sum;
    }
  }
}

void av1_set_visual_masking_rdmult(const AV1_COMP *const cpi, int *errorperbit,
                                   const BLOCK_SIZE bsize, const int mi_row,
                                   const int mi_col, int *const rdmult) {
  const AV1_COMMON *const cm = &cpi->common;

  const int bsize_base = BLOCK_4X4;
  const int num_mi_w = mi_size_wide[bsize_base];
  const int num_mi_h = mi_size_high[bsize_base];
  const int num_cols = (cm->mi_params.mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (cm->mi_params.mi_rows + num_mi_h - 1) / num_mi_h;
  const int num_bcols = (mi_size_wide[bsize] + num_mi_w - 1) / num_mi_w;
  const int num_brows = (mi_size_high[bsize] + num_mi_h - 1) / num_mi_h;
  int row, col;
  double num_of_mi = 0.0;
  double geom_mean_of_scale = 0.0;

  assert(cpi->oxcf.tune_cfg.tuning == AOM_TUNE_VISUAL_MASKING);

  for (row = mi_row / num_mi_w;
       row < num_rows && row < mi_row / num_mi_w + num_brows; ++row) {
    for (col = mi_col / num_mi_h;
         col < num_cols && col < mi_col / num_mi_h + num_bcols; ++col) {
      const int index = row * num_cols + col;
      geom_mean_of_scale +=
          log(cpi->visual_masking_rdmult_scaling_factors[index]);
      num_of_mi += 1.0;
    }
  }
  geom_mean_of_scale = exp(geom_mean_of_scale / num_of_mi);

  *rdmult = (int)((double)(*rdmult) * geom_mean_of_scale + 0.5);
  *rdmult = AOMMAX(*rdmult, 0);
  av1_set_error_per_bit(errorperbit, *rdmult);
}
