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

#ifndef AOM_AV1_ENCODER_RDOPT_UTILS_H_
#define AOM_AV1_ENCODER_RDOPT_UTILS_H_

#include "aom/aom_integer.h"
#include "av1/encoder/block.h"
#include "av1/common/blockd.h"
#include "av1/encoder/encoder.h"
#include "config/aom_dsp_rtcd.h"

#ifdef __cplusplus
extern "C" {
#endif

static const THR_MODES intra_to_mode_idx[INTRA_MODE_NUM] = {
  THR_DC,         // DC_PRED,
  THR_V_PRED,     // V_PRED,
  THR_H_PRED,     // H_PRED,
  THR_D45_PRED,   // D45_PRED,
  THR_D135_PRED,  // D135_PRED,
  THR_D113_PRED,  // D113_PRED,
  THR_D157_PRED,  // D157_PRED,
  THR_D203_PRED,  // D203_PRED,
  THR_D67_PRED,   // D67_PRED,
  THR_SMOOTH,     // SMOOTH_PRED,
  THR_SMOOTH_V,   // SMOOTH_V_PRED,
  THR_SMOOTH_H,   // SMOOTH_H_PRED,
  THR_PAETH,      // PAETH_PRED,
};

/* clang-format off */
    static const THR_MODES single_inter_to_mode_idx[SINGLE_INTER_MODE_NUM]
        [REF_FRAMES] = {
        // NEARESTMV,
            { THR_INVALID, THR_NEARESTMV, THR_NEARESTL2, THR_NEARESTL3,
            THR_NEARESTG, THR_NEARESTB, THR_NEARESTA2, THR_NEARESTA, },
            // NEARMV,
            { THR_INVALID, THR_NEARMV, THR_NEARL2, THR_NEARL3,
            THR_NEARG, THR_NEARB, THR_NEARA2, THR_NEARA, },
            // GLOBALMV,
            { THR_INVALID, THR_GLOBALMV, THR_GLOBALL2, THR_GLOBALL3,
            THR_GLOBALG, THR_GLOBALB, THR_GLOBALA2, THR_GLOBALA, },
            // NEWMV,
            { THR_INVALID, THR_NEWMV, THR_NEWL2, THR_NEWL3,
            THR_NEWG, THR_NEWB, THR_NEWA2, THR_NEWA, },
    };
/* clang-format on */

/* clang-format off */
static const THR_MODES comp_inter_to_mode_idx[COMP_INTER_MODE_NUM][REF_FRAMES]
    [REF_FRAMES] = {
    // NEAREST_NEARESTMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEAREST_NEARESTLL2, THR_COMP_NEAREST_NEARESTLL3,
            THR_COMP_NEAREST_NEARESTLG, THR_COMP_NEAREST_NEARESTLB,
            THR_COMP_NEAREST_NEARESTLA2, THR_COMP_NEAREST_NEARESTLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEARESTL2B,
            THR_COMP_NEAREST_NEARESTL2A2, THR_COMP_NEAREST_NEARESTL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEARESTL3B,
            THR_COMP_NEAREST_NEARESTL3A2, THR_COMP_NEAREST_NEARESTL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEARESTGB,
            THR_COMP_NEAREST_NEARESTGA2, THR_COMP_NEAREST_NEARESTGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEARESTBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEAR_NEARMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEAR_NEARLL2, THR_COMP_NEAR_NEARLL3,
            THR_COMP_NEAR_NEARLG, THR_COMP_NEAR_NEARLB,
            THR_COMP_NEAR_NEARLA2, THR_COMP_NEAR_NEARLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEARL2B,
            THR_COMP_NEAR_NEARL2A2, THR_COMP_NEAR_NEARL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEARL3B,
            THR_COMP_NEAR_NEARL3A2, THR_COMP_NEAR_NEARL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEARGB,
            THR_COMP_NEAR_NEARGA2, THR_COMP_NEAR_NEARGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEARBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEAREST_NEWMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEAREST_NEWLL2, THR_COMP_NEAREST_NEWLL3,
            THR_COMP_NEAREST_NEWLG, THR_COMP_NEAREST_NEWLB,
            THR_COMP_NEAREST_NEWLA2, THR_COMP_NEAREST_NEWLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEWL2B,
            THR_COMP_NEAREST_NEWL2A2, THR_COMP_NEAREST_NEWL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEWL3B,
            THR_COMP_NEAREST_NEWL3A2, THR_COMP_NEAREST_NEWL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEWGB,
            THR_COMP_NEAREST_NEWGA2, THR_COMP_NEAREST_NEWGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAREST_NEWBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEW_NEARESTMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEW_NEARESTLL2, THR_COMP_NEW_NEARESTLL3,
            THR_COMP_NEW_NEARESTLG, THR_COMP_NEW_NEARESTLB,
            THR_COMP_NEW_NEARESTLA2, THR_COMP_NEW_NEARESTLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARESTL2B,
            THR_COMP_NEW_NEARESTL2A2, THR_COMP_NEW_NEARESTL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARESTL3B,
            THR_COMP_NEW_NEARESTL3A2, THR_COMP_NEW_NEARESTL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARESTGB,
            THR_COMP_NEW_NEARESTGA2, THR_COMP_NEW_NEARESTGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARESTBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEAR_NEWMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEAR_NEWLL2, THR_COMP_NEAR_NEWLL3,
            THR_COMP_NEAR_NEWLG, THR_COMP_NEAR_NEWLB,
            THR_COMP_NEAR_NEWLA2, THR_COMP_NEAR_NEWLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEWL2B,
            THR_COMP_NEAR_NEWL2A2, THR_COMP_NEAR_NEWL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEWL3B,
            THR_COMP_NEAR_NEWL3A2, THR_COMP_NEAR_NEWL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEWGB,
            THR_COMP_NEAR_NEWGA2, THR_COMP_NEAR_NEWGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEAR_NEWBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEW_NEARMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEW_NEARLL2, THR_COMP_NEW_NEARLL3,
            THR_COMP_NEW_NEARLG, THR_COMP_NEW_NEARLB,
            THR_COMP_NEW_NEARLA2, THR_COMP_NEW_NEARLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARL2B,
            THR_COMP_NEW_NEARL2A2, THR_COMP_NEW_NEARL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARL3B,
            THR_COMP_NEW_NEARL3A2, THR_COMP_NEW_NEARL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARGB,
            THR_COMP_NEW_NEARGA2, THR_COMP_NEW_NEARGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEARBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // GLOBAL_GLOBALMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_GLOBAL_GLOBALLL2, THR_COMP_GLOBAL_GLOBALLL3,
            THR_COMP_GLOBAL_GLOBALLG, THR_COMP_GLOBAL_GLOBALLB,
            THR_COMP_GLOBAL_GLOBALLA2, THR_COMP_GLOBAL_GLOBALLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_GLOBAL_GLOBALL2B,
            THR_COMP_GLOBAL_GLOBALL2A2, THR_COMP_GLOBAL_GLOBALL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_GLOBAL_GLOBALL3B,
            THR_COMP_GLOBAL_GLOBALL3A2, THR_COMP_GLOBAL_GLOBALL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_GLOBAL_GLOBALGB,
            THR_COMP_GLOBAL_GLOBALGA2, THR_COMP_GLOBAL_GLOBALGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_GLOBAL_GLOBALBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
        // NEW_NEWMV,
        {
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID,
            THR_COMP_NEW_NEWLL2, THR_COMP_NEW_NEWLL3,
            THR_COMP_NEW_NEWLG, THR_COMP_NEW_NEWLB,
            THR_COMP_NEW_NEWLA2, THR_COMP_NEW_NEWLA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEWL2B,
            THR_COMP_NEW_NEWL2A2, THR_COMP_NEW_NEWL2A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEWL3B,
            THR_COMP_NEW_NEWL3A2, THR_COMP_NEW_NEWL3A, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEWGB,
            THR_COMP_NEW_NEWGA2, THR_COMP_NEW_NEWGA, },
            { THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID,
            THR_INVALID, THR_COMP_NEW_NEWBA, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
            { THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID, THR_INVALID,
            THR_INVALID, THR_INVALID, THR_INVALID, },
        },
};

/* clang-format on */
// Calculate rd threshold based on ref best rd and relevant scaling factors
static INLINE int64_t get_rd_thresh_from_best_rd(int64_t ref_best_rd,
                                                 int mul_factor,
                                                 int div_factor) {
  int64_t rd_thresh = ref_best_rd;
  if (div_factor != 0) {
    rd_thresh = ref_best_rd < (div_factor * (INT64_MAX / mul_factor))
                    ? ((ref_best_rd / div_factor) * mul_factor)
                    : INT64_MAX;
  }
  return rd_thresh;
}

static THR_MODES get_prediction_mode_idx(PREDICTION_MODE this_mode,
                                         MV_REFERENCE_FRAME ref_frame,
                                         MV_REFERENCE_FRAME second_ref_frame) {
  if (this_mode < INTRA_MODE_END) {
    assert(ref_frame == INTRA_FRAME);
    assert(second_ref_frame == NONE_FRAME);
    return intra_to_mode_idx[this_mode - INTRA_MODE_START];
  }
  if (this_mode >= SINGLE_INTER_MODE_START &&
      this_mode < SINGLE_INTER_MODE_END) {
    assert((ref_frame > INTRA_FRAME) && (ref_frame <= ALTREF_FRAME));
    return single_inter_to_mode_idx[this_mode - SINGLE_INTER_MODE_START]
                                   [ref_frame];
  }
  if (this_mode >= COMP_INTER_MODE_START && this_mode < COMP_INTER_MODE_END) {
    assert((ref_frame > INTRA_FRAME) && (ref_frame <= ALTREF_FRAME));
    assert((second_ref_frame > INTRA_FRAME) &&
           (second_ref_frame <= ALTREF_FRAME));
    return comp_inter_to_mode_idx[this_mode - COMP_INTER_MODE_START][ref_frame]
                                 [second_ref_frame];
  }
  assert(0);
  return THR_INVALID;
}

static int inter_mode_data_block_idx(BLOCK_SIZE bsize) {
  if (bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
      bsize == BLOCK_4X16 || bsize == BLOCK_16X4) {
    return -1;
  }
  return 1;
}

static AOM_INLINE double get_sad_norm(const int16_t *diff, int stride, int w,
                                      int h) {
  double sum = 0.0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      sum += abs(diff[j * stride + i]);
    }
  }
  assert(w > 0 && h > 0);
  return sum / (w * h);
}

static AOM_INLINE double get_sse_norm(const int16_t *diff, int stride, int w,
                                      int h) {
  double sum = 0.0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      const int err = diff[j * stride + i];
      sum += err * err;
    }
  }
  assert(w > 0 && h > 0);
  return sum / (w * h);
}

static AOM_INLINE void get_2x2_normalized_sses_and_sads(
    const AV1_COMP *const cpi, BLOCK_SIZE tx_bsize, const uint8_t *const src,
    int src_stride, const uint8_t *const dst, int dst_stride,
    const int16_t *const src_diff, int diff_stride, double *const sse_norm_arr,
    double *const sad_norm_arr) {
  const BLOCK_SIZE tx_bsize_half =
      get_partition_subsize(tx_bsize, PARTITION_SPLIT);
  if (tx_bsize_half == BLOCK_INVALID) {  // manually calculate stats
    const int half_width = block_size_wide[tx_bsize] / 2;
    const int half_height = block_size_high[tx_bsize] / 2;
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 2; ++col) {
        const int16_t *const this_src_diff =
            src_diff + row * half_height * diff_stride + col * half_width;
        if (sse_norm_arr) {
          sse_norm_arr[row * 2 + col] =
              get_sse_norm(this_src_diff, diff_stride, half_width, half_height);
        }
        if (sad_norm_arr) {
          sad_norm_arr[row * 2 + col] =
              get_sad_norm(this_src_diff, diff_stride, half_width, half_height);
        }
      }
    }
  } else {  // use function pointers to calculate stats
    const int half_width = block_size_wide[tx_bsize_half];
    const int half_height = block_size_high[tx_bsize_half];
    const int num_samples_half = half_width * half_height;
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 2; ++col) {
        const uint8_t *const this_src =
            src + row * half_height * src_stride + col * half_width;
        const uint8_t *const this_dst =
            dst + row * half_height * dst_stride + col * half_width;

        if (sse_norm_arr) {
          unsigned int this_sse;
          cpi->fn_ptr[tx_bsize_half].vf(this_src, src_stride, this_dst,
                                        dst_stride, &this_sse);
          sse_norm_arr[row * 2 + col] = (double)this_sse / num_samples_half;
        }

        if (sad_norm_arr) {
          const unsigned int this_sad = cpi->fn_ptr[tx_bsize_half].sdf(
              this_src, src_stride, this_dst, dst_stride);
          sad_norm_arr[row * 2 + col] = (double)this_sad / num_samples_half;
        }
      }
    }
  }
}

static int64_t get_sse(const AV1_COMP *cpi, const MACROBLOCK *x) {
  const AV1_COMMON *cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  const MACROBLOCKD *xd = &x->e_mbd;
  const MB_MODE_INFO *mbmi = xd->mi[0];
  int64_t total_sse = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const struct macroblock_plane *const p = &x->plane[plane];
    const struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE bs = get_plane_block_size(mbmi->sb_type, pd->subsampling_x,
                                               pd->subsampling_y);
    unsigned int sse;

    if (x->skip_chroma_rd && plane) continue;

    cpi->fn_ptr[bs].vf(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
                       &sse);
    total_sse += sse;
  }
  total_sse <<= 4;
  return total_sse;
}

static AOM_INLINE int64_t av1_block_error_c(const tran_low_t *coeff,
                                            const tran_low_t *dqcoeff,
                                            intptr_t block_size, int64_t *ssz) {
  int i;
  int64_t error = 0, sqcoeff = 0;

  for (i = 0; i < block_size; i++) {
    const int diff = coeff[i] - dqcoeff[i];
    error += diff * diff;
    sqcoeff += coeff[i] * coeff[i];
  }

  *ssz = sqcoeff;
  return error;
}

#if CONFIG_AV1_HIGHBITDEPTH
static int64_t av1_highbd_block_error_c(const tran_low_t *coeff,
                                        const tran_low_t *dqcoeff,
                                        intptr_t block_size, int64_t *ssz,
                                        int bd) {
  int i;
  int64_t error = 0, sqcoeff = 0;
  int shift = 2 * (bd - 8);
  int rounding = shift > 0 ? 1 << (shift - 1) : 0;

  for (i = 0; i < block_size; i++) {
    const int64_t diff = coeff[i] - dqcoeff[i];
    error += diff * diff;
    sqcoeff += (int64_t)coeff[i] * (int64_t)coeff[i];
  }
  assert(error >= 0 && sqcoeff >= 0);
  error = (error + rounding) >> shift;
  sqcoeff = (sqcoeff + rounding) >> shift;

  *ssz = sqcoeff;
  return error;
}
#endif

static AOM_INLINE double get_highbd_diff_mean(const uint8_t *src8,
                                              int src_stride,
                                              const uint8_t *dst8,
                                              int dst_stride, int w, int h) {
  const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  const uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
  double sum = 0.0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      const int diff = src[j * src_stride + i] - dst[j * dst_stride + i];
      sum += diff;
    }
  }
  assert(w > 0 && h > 0);
  return sum / (w * h);
}

static AOM_INLINE double get_diff_mean(const uint8_t *src, int src_stride,
                                       const uint8_t *dst, int dst_stride,
                                       int w, int h) {
  double sum = 0.0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      const int diff = src[j * src_stride + i] - dst[j * dst_stride + i];
      sum += diff;
    }
  }
  assert(w > 0 && h > 0);
  return sum / (w * h);
}

static int64_t calculate_sse(MACROBLOCKD *const xd,
                             const struct macroblock_plane *p,
                             struct macroblockd_plane *pd, const int bw,
                             const int bh) {
  int64_t sse = 0;
  const int shift = xd->bd - 8;
#if CONFIG_AV1_HIGHBITDEPTH
  if (is_cur_buf_hbd(xd)) {
    sse = aom_highbd_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
                         bw, bh);
  } else {
    sse =
        aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
  }
#else
  sse = aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
#endif
  sse = ROUND_POWER_OF_TWO(sse, shift * 2);
  return sse;
}

#define LEFT_TOP_MARGIN ((AOM_BORDER_IN_PIXELS - AOM_INTERP_EXTEND) << 3)
#define RIGHT_BOTTOM_MARGIN ((AOM_BORDER_IN_PIXELS - AOM_INTERP_EXTEND) << 3)

// TODO(jingning): this mv clamping function should be block size dependent.
static INLINE void clamp_mv2(MV *mv, const MACROBLOCKD *xd) {
  clamp_mv(mv, xd->mb_to_left_edge - LEFT_TOP_MARGIN,
           xd->mb_to_right_edge + RIGHT_BOTTOM_MARGIN,
           xd->mb_to_top_edge - LEFT_TOP_MARGIN,
           xd->mb_to_bottom_edge + RIGHT_BOTTOM_MARGIN);
}

#define BINS 32
static const float intra_hog_model_bias[DIRECTIONAL_MODES] = {
  0.450578f,  0.695518f,  -0.717944f, -0.639894f,
  -0.602019f, -0.453454f, 0.055857f,  -0.465480f,
};

static const float intra_hog_model_weights[BINS * DIRECTIONAL_MODES] = {
  -3.076402f, -3.757063f, -3.275266f, -3.180665f, -3.452105f, -3.216593f,
  -2.871212f, -3.134296f, -1.822324f, -2.401411f, -1.541016f, -1.195322f,
  -0.434156f, 0.322868f,  2.260546f,  3.368715f,  3.989290f,  3.308487f,
  2.277893f,  0.923793f,  0.026412f,  -0.385174f, -0.718622f, -1.408867f,
  -1.050558f, -2.323941f, -2.225827f, -2.585453f, -3.054283f, -2.875087f,
  -2.985709f, -3.447155f, 3.758139f,  3.204353f,  2.170998f,  0.826587f,
  -0.269665f, -0.702068f, -1.085776f, -2.175249f, -1.623180f, -2.975142f,
  -2.779629f, -3.190799f, -3.521900f, -3.375480f, -3.319355f, -3.897389f,
  -3.172334f, -3.594528f, -2.879132f, -2.547777f, -2.921023f, -2.281844f,
  -1.818988f, -2.041771f, -0.618268f, -1.396458f, -0.567153f, -0.285868f,
  -0.088058f, 0.753494f,  2.092413f,  3.215266f,  -3.300277f, -2.748658f,
  -2.315784f, -2.423671f, -2.257283f, -2.269583f, -2.196660f, -2.301076f,
  -2.646516f, -2.271319f, -2.254366f, -2.300102f, -2.217960f, -2.473300f,
  -2.116866f, -2.528246f, -3.314712f, -1.701010f, -0.589040f, -0.088077f,
  0.813112f,  1.702213f,  2.653045f,  3.351749f,  3.243554f,  3.199409f,
  2.437856f,  1.468854f,  0.533039f,  -0.099065f, -0.622643f, -2.200732f,
  -4.228861f, -2.875263f, -1.273956f, -0.433280f, 0.803771f,  1.975043f,
  3.179528f,  3.939064f,  3.454379f,  3.689386f,  3.116411f,  1.970991f,
  0.798406f,  -0.628514f, -1.252546f, -2.825176f, -4.090178f, -3.777448f,
  -3.227314f, -3.479403f, -3.320569f, -3.159372f, -2.729202f, -2.722341f,
  -3.054913f, -2.742923f, -2.612703f, -2.662632f, -2.907314f, -3.117794f,
  -3.102660f, -3.970972f, -4.891357f, -3.935582f, -3.347758f, -2.721924f,
  -2.219011f, -1.702391f, -0.866529f, -0.153743f, 0.107733f,  1.416882f,
  2.572884f,  3.607755f,  3.974820f,  3.997783f,  2.970459f,  0.791687f,
  -1.478921f, -1.228154f, -1.216955f, -1.765932f, -1.951003f, -1.985301f,
  -1.975881f, -1.985593f, -2.422371f, -2.419978f, -2.531288f, -2.951853f,
  -3.071380f, -3.277027f, -3.373539f, -4.462010f, -0.967888f, 0.805524f,
  2.794130f,  3.685984f,  3.745195f,  3.252444f,  2.316108f,  1.399146f,
  -0.136519f, -0.162811f, -1.004357f, -1.667911f, -1.964662f, -2.937579f,
  -3.019533f, -3.942766f, -5.102767f, -3.882073f, -3.532027f, -3.451956f,
  -2.944015f, -2.643064f, -2.529872f, -2.077290f, -2.809965f, -1.803734f,
  -1.783593f, -1.662585f, -1.415484f, -1.392673f, -0.788794f, -1.204819f,
  -1.998864f, -1.182102f, -0.892110f, -1.317415f, -1.359112f, -1.522867f,
  -1.468552f, -1.779072f, -2.332959f, -2.160346f, -2.329387f, -2.631259f,
  -2.744936f, -3.052494f, -2.787363f, -3.442548f, -4.245075f, -3.032172f,
  -2.061609f, -1.768116f, -1.286072f, -0.706587f, -0.192413f, 0.386938f,
  0.716997f,  1.481393f,  2.216702f,  2.737986f,  3.109809f,  3.226084f,
  2.490098f,  -0.095827f, -3.864816f, -3.507248f, -3.128925f, -2.908251f,
  -2.883836f, -2.881411f, -2.524377f, -2.624478f, -2.399573f, -2.367718f,
  -1.918255f, -1.926277f, -1.694584f, -1.723790f, -0.966491f, -1.183115f,
  -1.430687f, 0.872896f,  2.766550f,  3.610080f,  3.578041f,  3.334928f,
  2.586680f,  1.895721f,  1.122195f,  0.488519f,  -0.140689f, -0.799076f,
  -1.222860f, -1.502437f, -1.900969f, -3.206816f,
};

static void generate_hog(const uint8_t *src, int stride, int rows, int cols,
                         float *hist) {
  const float step = (float)PI / BINS;
  float total = 0.1f;
  src += stride;
  for (int r = 1; r < rows - 1; ++r) {
    for (int c = 1; c < cols - 1; ++c) {
      const uint8_t *above = &src[c - stride];
      const uint8_t *below = &src[c + stride];
      const uint8_t *left = &src[c - 1];
      const uint8_t *right = &src[c + 1];
      // Calculate gradient using Sobel fitlers.
      const int dx = (right[-stride] + 2 * right[0] + right[stride]) -
                     (left[-stride] + 2 * left[0] + left[stride]);
      const int dy = (below[-1] + 2 * below[0] + below[1]) -
                     (above[-1] + 2 * above[0] + above[1]);
      if (dx == 0 && dy == 0) continue;
      const int temp = abs(dx) + abs(dy);
      if (!temp) continue;
      total += temp;
      if (dx == 0) {
        hist[0] += temp / 2;
        hist[BINS - 1] += temp / 2;
      } else {
        const float angle = atanf(dy * 1.0f / dx);
        int idx = (int)roundf(angle / step) + BINS / 2;
        idx = AOMMIN(idx, BINS - 1);
        idx = AOMMAX(idx, 0);
        hist[idx] += temp;
      }
    }
    src += stride;
  }

  for (int i = 0; i < BINS; ++i) hist[i] /= total;
}

static void generate_hog_hbd(const uint8_t *src8, int stride, int rows,
                             int cols, float *hist) {
  const float step = (float)PI / BINS;
  float total = 0.1f;
  uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  src += stride;
  for (int r = 1; r < rows - 1; ++r) {
    for (int c = 1; c < cols - 1; ++c) {
      const uint16_t *above = &src[c - stride];
      const uint16_t *below = &src[c + stride];
      const uint16_t *left = &src[c - 1];
      const uint16_t *right = &src[c + 1];
      // Calculate gradient using Sobel fitlers.
      const int dx = (right[-stride] + 2 * right[0] + right[stride]) -
                     (left[-stride] + 2 * left[0] + left[stride]);
      const int dy = (below[-1] + 2 * below[0] + below[1]) -
                     (above[-1] + 2 * above[0] + above[1]);
      if (dx == 0 && dy == 0) continue;
      const int temp = abs(dx) + abs(dy);
      if (!temp) continue;
      total += temp;
      if (dx == 0) {
        hist[0] += temp / 2;
        hist[BINS - 1] += temp / 2;
      } else {
        const float angle = atanf(dy * 1.0f / dx);
        int idx = (int)roundf(angle / step) + BINS / 2;
        idx = AOMMIN(idx, BINS - 1);
        idx = AOMMAX(idx, 0);
        hist[idx] += temp;
      }
    }
    src += stride;
  }

  for (int i = 0; i < BINS; ++i) hist[i] /= total;
}

static INLINE void restore_dst_buf(MACROBLOCKD *xd, const BUFFER_SET dst,
                                   const int num_planes) {
  for (int i = 0; i < num_planes; i++) {
    xd->plane[i].dst.buf = dst.plane[i];
    xd->plane[i].dst.stride = dst.stride[i];
  }
}

static INLINE void swap_dst_buf(MACROBLOCKD *xd, const BUFFER_SET *dst_bufs[2],
                                int num_planes) {
  const BUFFER_SET *buf0 = dst_bufs[0];
  dst_bufs[0] = dst_bufs[1];
  dst_bufs[1] = buf0;
  restore_dst_buf(xd, *dst_bufs[0], num_planes);
}

// Represents a set of integers, from 0 to sizeof(int) * 8, as bits in
// an integer. 0 for the i-th bit means that integer is excluded, 1 means
// it is included.
static INLINE void mask_set_bit(int *mask, int index) { *mask |= (1 << index); }

static INLINE bool mask_check_bit(int mask, int index) {
  return (mask >> index) & 0x1;
}

// Updates best_mv for masked compound types
static INLINE void update_mask_best_mv(const MB_MODE_INFO *const mbmi,
                                       int_mv *best_mv, int_mv *cur_mv,
                                       const COMPOUND_TYPE cur_type,
                                       int *best_tmp_rate_mv, int tmp_rate_mv,
                                       const SPEED_FEATURES *const sf) {
  if (cur_type == COMPOUND_WEDGE ||
      (sf->inter_sf.enable_interinter_diffwtd_newmv_search &&
       cur_type == COMPOUND_DIFFWTD)) {
    *best_tmp_rate_mv = tmp_rate_mv;
    best_mv[0].as_int = mbmi->mv[0].as_int;
    best_mv[1].as_int = mbmi->mv[1].as_int;
  } else {
    best_mv[0].as_int = cur_mv[0].as_int;
    best_mv[1].as_int = cur_mv[1].as_int;
  }
}

static INLINE int bsize_to_num_blk(BLOCK_SIZE bsize) {
  int num_blk = 1 << (num_pels_log2_lookup[bsize] - 2 * MI_SIZE_LOG2);
  return num_blk;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_RDOPT_UTILS_H_
