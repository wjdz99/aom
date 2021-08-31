/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdio.h>
#include "av1/common/common_data.h"
#include <float.h>
#include "av1/common/enums.h"
#include "av1/common/idct.h"

#include "av1/common/reconinter.h"
#include "av1/encoder/allintra_vis.h"
#include "av1/encoder/hybrid_fwd_txfm.h"
#include "av1/encoder/rdopt_utils.h"

// Process the wiener variance in 16x16 block basis.
static int qsort_comp(const void *elem1, const void *elem2) {
  int a = *((const int *)elem1);
  int b = *((const int *)elem2);
  if (a > b) return 1;
  if (a < b) return -1;
  return 0;
}

void av1_init_mb_wiener_var_buffer(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;

  cpi->weber_bsize = BLOCK_8X8;

  if (cpi->mb_weber_stats) return;

  CHECK_MEM_ERROR(cm, cpi->mb_weber_stats,
                  aom_calloc(cpi->frame_info.mi_rows * cpi->frame_info.mi_cols,
                             sizeof(*cpi->mb_weber_stats)));
}

static int64_t get_satd(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                        int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  const int mi_step = mi_size_wide[cpi->weber_bsize];
  int mb_stride = cpi->frame_info.mi_cols;
  int mb_count = 0;
  int64_t satd = 0;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;

      satd += cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)]
                  .satd;
      ++mb_count;
    }
  }

  if (mb_count) satd = (int)(satd / mb_count);
  satd = AOMMAX(1, satd);

  return (int)satd;
}

static int64_t get_sse(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                       int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  const int mi_step = mi_size_wide[cpi->weber_bsize];
  int mb_stride = cpi->frame_info.mi_cols;
  int mb_count = 0;
  int64_t distortion = 0;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;

      distortion +=
          cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)]
              .distortion;
      ++mb_count;
    }
  }

  if (mb_count) distortion = (int)(distortion / mb_count);
  distortion = AOMMAX(1, distortion);

  return (int)distortion;
}

static double get_max_scale(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                            int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];
  const int mi_step = mi_size_wide[cpi->weber_bsize];
  int mb_stride = cpi->frame_info.mi_cols;
  double min_max_scale = 10.0;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;
      WeberStats *weber_stats =
          &cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)];
      if (weber_stats->max_scale < 1.0) continue;
      if (weber_stats->max_scale < min_max_scale)
        min_max_scale = weber_stats->max_scale;
    }
  }
  return min_max_scale;
}

static int get_window_wiener_var(AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                 int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  const int mi_step = mi_size_wide[cpi->weber_bsize];
  int sb_wiener_var = 0;
  int mb_stride = cpi->frame_info.mi_cols;
  int mb_count = 0;
  double base_num = 1;
  double base_den = 1;
  double base_reg = 1;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;

      WeberStats *weber_stats =
          &cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)];

      base_num += ((double)weber_stats->distortion) *
                  sqrt((double)weber_stats->src_variance) *
                  weber_stats->rec_pix_max;

      base_den += fabs(
          weber_stats->rec_pix_max * sqrt((double)weber_stats->src_variance) -
          weber_stats->src_pix_max * sqrt((double)weber_stats->rec_variance));

      base_reg += sqrt((double)weber_stats->distortion) *
                  sqrt((double)weber_stats->src_pix_max) * 0.1;
      ++mb_count;
    }
  }

  sb_wiener_var = (int)((base_num + base_reg) / (base_den + base_reg));
  sb_wiener_var = AOMMAX(1, sb_wiener_var);

  return (int)sb_wiener_var;
}

static int get_var_perceptual_ai(AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                 int mi_row, int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  int sb_wiener_var = get_window_wiener_var(cpi, bsize, mi_row, mi_col);

  if (mi_row >= (mi_high / 2)) {
    sb_wiener_var =
        AOMMIN(sb_wiener_var,
               get_window_wiener_var(cpi, bsize, mi_row - mi_high / 2, mi_col));
  }
  if (mi_row <= (cm->mi_params.mi_rows - mi_high - (mi_high / 2))) {
    sb_wiener_var =
        AOMMIN(sb_wiener_var,
               get_window_wiener_var(cpi, bsize, mi_row + mi_high / 2, mi_col));
  }
  if (mi_col >= (mi_wide / 2)) {
    sb_wiener_var =
        AOMMIN(sb_wiener_var,
               get_window_wiener_var(cpi, bsize, mi_row, mi_col - mi_wide / 2));
  }
  if (mi_col <= (cm->mi_params.mi_cols - mi_wide - (mi_wide / 2))) {
    sb_wiener_var =
        AOMMIN(sb_wiener_var,
               get_window_wiener_var(cpi, bsize, mi_row, mi_col + mi_wide / 2));
  }

  return sb_wiener_var;
}

// Copy from tpl_model.c
static INLINE double exp_bounded(double v) {
  // When v > 700 or <-700, the exp function will be close to overflow
  // For details, see the "Notes" in the following link.
  // https://en.cppreference.com/w/c/numeric/math/exp
  if (v > 700) {
    return DBL_MAX;
  } else if (v < -700) {
    return 0;
  }
  return exp(v);
}

// Copy from tpl_model.c
#define TPL_EPSILON 0.0000001
double estimate_coeff_entropy(double q_step, double b, double zero_bin_ratio,
                              int qcoeff) {
  b = AOMMAX(b, TPL_EPSILON);
  int abs_qcoeff = abs(qcoeff);
  double z0 = fmax(exp_bounded(-zero_bin_ratio / 2 * q_step / b), TPL_EPSILON);
  if (abs_qcoeff == 0) {
    double r = -log2(1 - z0);
    return r;
  } else {
    double z = fmax(exp_bounded(-q_step / b), TPL_EPSILON);
    double r = 1 - log2(z0) - log2(1 - z) - (abs_qcoeff - 1) * log2(z);
    return r;
  }
}

// Copy from tpl_model.c
static int rate_estimator(const tran_low_t *qcoeff, int eob, TX_SIZE tx_size) {
  const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];

  assert((1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]]) >= eob);
  int rate_cost = 1;

  for (int idx = 0; idx < eob; ++idx) {
    int abs_level = abs(qcoeff[scan_order->scan[idx]]);
    rate_cost += (int)(log(abs_level + 1.0) / log(2.0)) + 1;
  }

  return (rate_cost << AV1_PROB_COST_SHIFT);
}

// Copy from tpl_model.c
static double estimate_txfm_block_entropy(int q_index,
                                          const double *abs_coeff_mean,
                                          int *qcoeff_arr, int coeff_num) {
  double zero_bin_ratio = 2;
  double ac_q_step = av1_ac_quant_QTX(q_index, 0, AOM_BITS_8) / 4.;
  double est_rate = 0;
  // ac coeff
  for (int i = 1; i < coeff_num; ++i) {
    est_rate += estimate_coeff_entropy(ac_q_step, abs_coeff_mean[i],
                                       zero_bin_ratio, qcoeff_arr[i]);
  }
  return est_rate;
}

static void adjust_qindex(int *base_qindex, const int last_qindex,
                          const double avg_mb_entropy, const double mb_entropy,
                          double *upper_mb_entropy, double *lower_mb_entropy,
                          double *orig_entropy, const double avg_var,
                          const double this_mb_var, const double avg_saliency,
                          const double block_saliency, const double avg_mscn,
                          const double block_mscn, const int iter,
                          const int sb_row, const int sb_col) {
  (void)sb_row;
  (void)sb_col;
  // If a block has a large variance, we expect the block's entropy is above
  // average and it tolerates higher entropy.
  const double var_ratio = sqrt(this_mb_var / avg_var);

  // If a block has a high saliency score, we expect to use more bits/entropy.
  const double saliency_ratio = block_saliency / avg_saliency;

  // If a block has a high mscn score, we expect the block has more fine
  // structure or context, which means it tolerates higher entropy.
  const double mscn_ratio = block_mscn / avg_mscn;

  // const double weights[] = { 1, 1, 1 };
  // const double ratio = weights[0] * var_ratio + weights[1] * (1 -
  // saliency_ratio) + weights[2] * mscn_ratio;

  if (iter == 0) return;
  if (iter == 1) {
    *orig_entropy = mb_entropy;
    // Adjust upper and lower mb entropy limit according to variance
    if (*orig_entropy > *upper_mb_entropy) {
      const double complexity = (var_ratio + mscn_ratio) / 2;
      /*
      if (sb_row == 11 && sb_col == 11) {
        printf("orig upper_mb_entropy %.0f, ", *upper_mb_entropy);
      }
      */
      *upper_mb_entropy *= (1.5 * complexity);
      /*
      if (sb_row == 11 && sb_col == 11) {
        printf("complexity %.2f, upper_mb_entropy %.0f\n", complexity,
      *upper_mb_entropy);
      }
      */
      *upper_mb_entropy = AOMMIN(*upper_mb_entropy, 2.5 * avg_mb_entropy);
    }
    if (*orig_entropy < *lower_mb_entropy) {
      const double complexity = (var_ratio + mscn_ratio) / 2;
      *lower_mb_entropy *= (1.5 * complexity);
      *lower_mb_entropy = AOMMAX(*lower_mb_entropy, 0.25 * avg_mb_entropy);
    }
  }

  // const double entropy_ratio = mb_entropy / avg_mb_entropy;

  if (*orig_entropy > avg_mb_entropy) {
    double delta_qindex = 10;
    if (mb_entropy > avg_mb_entropy) {
      if (var_ratio > 1 && mb_entropy > *upper_mb_entropy) {
        delta_qindex *= AOMMIN(var_ratio, 2);
      } else if (var_ratio < 1) {
        delta_qindex /= AOMMAX(var_ratio, 0.5);
      }
      if (saliency_ratio > 1) {
        delta_qindex /= AOMMIN(saliency_ratio, 4);
      } else {
        delta_qindex /= AOMMAX(saliency_ratio, 0.5);
      }
      if (mscn_ratio < 1) {
        delta_qindex /= AOMMAX(mscn_ratio, 0.5);
      }
      *base_qindex = AOMMIN(last_qindex + (int)delta_qindex, MAXQ);
    } else if (mb_entropy < *lower_mb_entropy) {
      *base_qindex = AOMMAX(last_qindex - (int)delta_qindex, 0);
    }
  } else {
    double delta_qindex = 10;
    if (mb_entropy < avg_mb_entropy) {
      if (var_ratio > 1) {
        delta_qindex *= var_ratio;
      } else {
        delta_qindex /= AOMMAX(var_ratio, 0.5);
      }
      if (saliency_ratio > 1) {
        delta_qindex *= saliency_ratio;
      } else {
        delta_qindex *= AOMMAX(saliency_ratio, 0.25);
      }
      if (mscn_ratio > 1) {
        delta_qindex *= AOMMIN(mscn_ratio, 3);
      }
      *base_qindex = AOMMAX(last_qindex - (int)delta_qindex, 0);
    } else if (mb_entropy > *upper_mb_entropy) {
      *base_qindex = AOMMIN(last_qindex + (int)delta_qindex, MAXQ);
    }
  }
  /*
  if (sb_row == 11 && sb_col == 11) {
    printf(
        "iter %d, var_ratio %.2f, saliency_ratio %.2f, mscn_ratio %.2f, "
        "orig_entropy %.2f, avg_mb_entropy %.2f, mb_entropy %.2f, "
        "base_qindex %d, last_qindex %d\n",
        iter, var_ratio, saliency_ratio, mscn_ratio, *orig_entropy,
        avg_mb_entropy, mb_entropy, *base_qindex, last_qindex);
  }
  */

  /*
  // old strategy based on entropy only.
  if (iter == 0) return;
  if (iter == 1) *orig_entropy = mb_entropy;

  if (*orig_entropy > avg_mb_entropy * 2) {
    // Handle cases where some superblocks are much more complex
    const int delta_qindex = 10;
    *base_qindex = AOMMIN(last_qindex + delta_qindex, MAXQ);
  } else if (*orig_entropy < avg_mb_entropy / 4) {
    // Handle cases where some superblocks are much easier
    int delta_qindex = 15;
    if (mb_entropy < avg_mb_entropy / 2) delta_qindex = 10;
    *base_qindex = AOMMAX(0, last_qindex - delta_qindex);
  } else {
    const int delta_qindex = 5;
    if (mb_entropy > upper_mb_entropy) {
      *base_qindex = AOMMIN(last_qindex + delta_qindex, MAXQ);
    } else if (mb_entropy < lower_mb_entropy) {
      *base_qindex = AOMMAX(0, last_qindex - delta_qindex);
    }
  }
  */
}

static void entropy_based_sb_q(
    AV1_COMP *cpi, double *abs_coeff_mean, const double sum_frame_entropy,
    const double sum_frame_rate, const double avg_var,
    const double *const saliency_map, const double avg_saliency,
    const double *const mscn_map, const double avg_mscn) {
  (void)sum_frame_rate;
  AV1_COMMON *const cm = &cpi->common;
  const SequenceHeader *const seq_params = cm->seq_params;
  if (aom_realloc_frame_buffer(
          &cm->cur_frame->buf, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, cm->features.byte_alignment, NULL, NULL,
          NULL, cpi->oxcf.tool_cfg.enable_global_motion))
    aom_internal_error(cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");
  cm->quant_params.base_qindex = cpi->oxcf.rc_cfg.cq_level;
  av1_frame_init_quantizer(cpi);
  DECLARE_ALIGNED(32, int16_t, src_diff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, coeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, qcoeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, dqcoeff[32 * 32]);
  uint8_t *buffer = cpi->source->y_buffer;
  int buf_stride = cpi->source->y_stride;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO mbmi;
  memset(&mbmi, 0, sizeof(mbmi));
  MB_MODE_INFO *mbmi_ptr = &mbmi;
  xd->mi = &mbmi_ptr;
  xd->cur_buf = cpi->source;
  const TX_SIZE tx_size = TX_16X16;
  const int block_size = tx_size_wide[tx_size];
  const BitDepthInfo bd_info = get_bit_depth_info(xd);

  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int sb_width = mi_size_wide[sb_size];
  const int sb_height = mi_size_high[sb_size];
  const BLOCK_SIZE bsize = BLOCK_16X16;
  const int mb_step = mi_size_wide[bsize];
  const int num_col_sbs = (mi_params->mi_cols + sb_width - 1) / sb_width;
  const int num_row_sbs = (mi_params->mi_rows + sb_height - 1) / sb_height;
  const int orig_base_qindex = cpi->common.quant_params.base_qindex;
  const int num_mbs = cpi->frame_info.mb_rows * cpi->frame_info.mb_cols;
  const double avg_mb_entropy = sum_frame_entropy / num_mbs;
  const int max_iter = 5;

  // Iterate though each super block
  for (int sb_row = 0; sb_row < num_row_sbs; ++sb_row) {
    for (int sb_col = 0; sb_col < num_col_sbs; ++sb_col) {
      const int sb_idx = sb_row * num_col_sbs + sb_col;
      double upper_mb_entropy = 1.25 * avg_mb_entropy;
      double lower_mb_entropy = 0.75 * avg_mb_entropy;
      int base_qindex = orig_base_qindex;
      double mb_entropy = 0;
      int iter = 0;
      int last_qindex = base_qindex;
      double orig_entropy = 0;
      const double this_mb_var = cpi->mb_variance[sb_idx];
      double block_saliency = 0;
      int pixel_count = 0;
      double block_mscn = 0;
      for (int r = 0; r < block_size_high[sb_size]; ++r) {
        for (int c = 0; c < block_size_wide[sb_size]; ++c) {
          const int row = sb_row * block_size_high[sb_size] + r;
          const int col = sb_col * block_size_wide[sb_size] + c;
          if (row >= cpi->frame_info.frame_height ||
              col >= cpi->frame_info.frame_width) {
            continue;
          }
          block_saliency +=
              saliency_map[row * cpi->frame_info.frame_width + col];
          block_mscn += fabs(mscn_map[row * cpi->frame_info.frame_width + col]);
          ++pixel_count;
        }
      }
      block_saliency /= pixel_count;
      block_mscn /= pixel_count;

      do {
        // Adjust qindex based on the last qindex and corresponding sb entropy
        adjust_qindex(&base_qindex, last_qindex, avg_mb_entropy, mb_entropy,
                      &upper_mb_entropy, &lower_mb_entropy, &orig_entropy,
                      avg_var, this_mb_var, avg_saliency, block_saliency,
                      avg_mscn, block_mscn, iter, sb_row, sb_col);
        last_qindex = base_qindex;
        cm->quant_params.base_qindex = base_qindex;
        av1_frame_init_quantizer(cpi);
        ++iter;
        mb_entropy = 0;
        int count = 0;

        // Iterater though each mb block (16x16) inside a super block
        for (int mb_row = 0; mb_row < sb_height / mb_step; ++mb_row) {
          for (int mb_col = 0; mb_col < sb_width / mb_step; ++mb_col) {
            const PREDICTION_MODE best_mode =
                cpi->est_best_mode[mb_row * cpi->frame_info.mb_cols + mb_col];
            const int mi_row = sb_row * sb_height + mb_row * mb_step;
            const int mi_col = sb_col * sb_width + mb_col * mb_step;
            if (mi_row >= mi_params->mi_rows || mi_col >= mi_params->mi_cols)
              continue;
            xd->up_available = mi_row > 0;
            xd->left_available = mi_col > 0;
            const int mi_width = mi_size_wide[bsize];
            const int mi_height = mi_size_high[bsize];
            set_mode_info_offsets(&cpi->common.mi_params, &cpi->mbmi_ext_info,
                                  x, xd, mi_row, mi_col);
            set_mi_row_col(xd, &xd->tile, mi_row, mi_height, mi_col, mi_width,
                           cm->mi_params.mi_rows, cm->mi_params.mi_cols);
            set_plane_n4(xd, mi_size_wide[bsize], mi_size_high[bsize],
                         av1_num_planes(cm));
            xd->mi[0]->bsize = bsize;
            xd->mi[0]->motion_mode = SIMPLE_TRANSLATION;
            av1_setup_dst_planes(xd->plane, bsize, &cm->cur_frame->buf, mi_row,
                                 mi_col, 0, av1_num_planes(cm));
            int dst_buffer_stride = xd->plane[0].dst.stride;
            uint8_t *dst_buffer = xd->plane[0].dst.buf;
            uint8_t *mb_buffer =
                buffer + mi_row * MI_SIZE * buf_stride + mi_col * MI_SIZE;
            av1_predict_intra_block(xd, cm->seq_params->sb_size,
                                    cm->seq_params->enable_intra_edge_filter,
                                    block_size, block_size, tx_size, best_mode,
                                    0, 0, FILTER_INTRA_MODES, dst_buffer,
                                    dst_buffer_stride, dst_buffer,
                                    dst_buffer_stride, 0, 0, 0);
            av1_subtract_block(bd_info, block_size, block_size, src_diff,
                               block_size, mb_buffer, buf_stride, dst_buffer,
                               dst_buffer_stride);
            av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);
            const struct macroblock_plane *const p = &x->plane[0];
            uint16_t eob;
            const SCAN_ORDER *const scan_order =
                &av1_scan_orders[tx_size][DCT_DCT];
            QUANT_PARAM quant_param;
            int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
            av1_setup_quant(tx_size, 0, AV1_XFORM_QUANT_FP, 0, &quant_param);
#if CONFIG_AV1_HIGHBITDEPTH
            if (is_cur_buf_hbd(xd)) {
              av1_highbd_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff,
                                            &eob, scan_order, &quant_param);
            } else {
              av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                                     scan_order, &quant_param);
            }
#else
            av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                                   scan_order, &quant_param);
#endif  // CONFIG_AV1_HIGHBITDEPTH

            const double entropy = estimate_txfm_block_entropy(
                base_qindex, abs_coeff_mean, qcoeff, pix_num);
            mb_entropy += entropy;
            ++count;
          }
        }
        mb_entropy /= count;
        /*
        if (sb_row == 11 && sb_col == 11) {
          printf("iter %d, base_qindex %d, avg_mb_entropy %.0f, mb_entropy %.0f,
        upper %.0f, lower %.0f\n", iter, base_qindex, avg_mb_entropy,
        mb_entropy, upper_mb_entropy, lower_mb_entropy);
        }
        */
      } while (iter < max_iter && (mb_entropy > upper_mb_entropy ||
                                   mb_entropy < lower_mb_entropy));
      if (iter == 1) orig_entropy = mb_entropy;
      /*
      FILE *pfile = fopen("entropy.stat", "a");
      //fprintf(pfile, "(%d, %d), iter %d, avg_mb_entropy %.0f, orig_entropy
      %.0f, entropy %.0f, ratio %.2f, base_qindex %d\n",
      //        sb_row, sb_col, iter, avg_mb_entropy, orig_entropy, mb_entropy,
      mb_entropy / orig_entropy, base_qindex);
      //fprintf(pfile, "%.4f\n", orig_entropy / avg_mb_entropy);
      fclose(pfile);
      */

      // printf("sb (%d, %d), iter %d\n", sb_row, sb_col, iter);
      cpi->sb_qindex[sb_idx] = base_qindex;
    }
  }
}

void av1_set_mb_wiener_variance(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  uint8_t *buffer = cpi->source->y_buffer;
  int buf_stride = cpi->source->y_stride;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO mbmi;
  memset(&mbmi, 0, sizeof(mbmi));
  MB_MODE_INFO *mbmi_ptr = &mbmi;
  xd->mi = &mbmi_ptr;
  xd->cur_buf = cpi->source;

  const SequenceHeader *const seq_params = cm->seq_params;
  if (aom_realloc_frame_buffer(
          &cm->cur_frame->buf, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, cm->features.byte_alignment, NULL, NULL,
          NULL, cpi->oxcf.tool_cfg.enable_global_motion))
    aom_internal_error(cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");

  cm->quant_params.base_qindex = cpi->oxcf.rc_cfg.cq_level;
  av1_frame_init_quantizer(cpi);

  DECLARE_ALIGNED(32, int16_t, src_diff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, coeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, qcoeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, dqcoeff[32 * 32]);

  int mi_row, mi_col;

  BLOCK_SIZE bsize = cpi->weber_bsize;
  const TX_SIZE tx_size = max_txsize_lookup[bsize];
  const int block_size = tx_size_wide[tx_size];
  const int coeff_count = block_size * block_size;

  const BitDepthInfo bd_info = get_bit_depth_info(xd);
  cpi->norm_wiener_variance = 0;
  int mb_step = mi_size_wide[bsize];

  for (mi_row = 0; mi_row < cpi->frame_info.mi_rows; mi_row += mb_step) {
    for (mi_col = 0; mi_col < cpi->frame_info.mi_cols; mi_col += mb_step) {
      PREDICTION_MODE best_mode = DC_PRED;
      int best_intra_cost = INT_MAX;

      xd->up_available = mi_row > 0;
      xd->left_available = mi_col > 0;

      const int mi_width = mi_size_wide[bsize];
      const int mi_height = mi_size_high[bsize];
      set_mode_info_offsets(&cpi->common.mi_params, &cpi->mbmi_ext_info, x, xd,
                            mi_row, mi_col);
      set_mi_row_col(xd, &xd->tile, mi_row, mi_height, mi_col, mi_width,
                     cm->mi_params.mi_rows, cm->mi_params.mi_cols);
      set_plane_n4(xd, mi_size_wide[bsize], mi_size_high[bsize],
                   av1_num_planes(cm));
      xd->mi[0]->bsize = bsize;
      xd->mi[0]->motion_mode = SIMPLE_TRANSLATION;

      av1_setup_dst_planes(xd->plane, bsize, &cm->cur_frame->buf, mi_row,
                           mi_col, 0, av1_num_planes(cm));

      int dst_buffer_stride = xd->plane[0].dst.stride;
      uint8_t *dst_buffer = xd->plane[0].dst.buf;
      uint8_t *mb_buffer =
          buffer + mi_row * MI_SIZE * buf_stride + mi_col * MI_SIZE;

      for (PREDICTION_MODE mode = INTRA_MODE_START; mode < INTRA_MODE_END;
           ++mode) {
        av1_predict_intra_block(
            xd, cm->seq_params->sb_size,
            cm->seq_params->enable_intra_edge_filter, block_size, block_size,
            tx_size, mode, 0, 0, FILTER_INTRA_MODES, dst_buffer,
            dst_buffer_stride, dst_buffer, dst_buffer_stride, 0, 0, 0);

        av1_subtract_block(bd_info, block_size, block_size, src_diff,
                           block_size, mb_buffer, buf_stride, dst_buffer,
                           dst_buffer_stride);
        av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);
        int intra_cost = aom_satd(coeff, coeff_count);
        if (intra_cost < best_intra_cost) {
          best_intra_cost = intra_cost;
          best_mode = mode;
        }
      }

      int idx;
      av1_predict_intra_block(xd, cm->seq_params->sb_size,
                              cm->seq_params->enable_intra_edge_filter,
                              block_size, block_size, tx_size, best_mode, 0, 0,
                              FILTER_INTRA_MODES, dst_buffer, dst_buffer_stride,
                              dst_buffer, dst_buffer_stride, 0, 0, 0);
      av1_subtract_block(bd_info, block_size, block_size, src_diff, block_size,
                         mb_buffer, buf_stride, dst_buffer, dst_buffer_stride);
      av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);

      const struct macroblock_plane *const p = &x->plane[0];
      uint16_t eob;
      const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];
      QUANT_PARAM quant_param;
      int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
      av1_setup_quant(tx_size, 0, AV1_XFORM_QUANT_FP, 0, &quant_param);
#if CONFIG_AV1_HIGHBITDEPTH
      if (is_cur_buf_hbd(xd)) {
        av1_highbd_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                                      scan_order, &quant_param);
      } else {
        av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                               scan_order, &quant_param);
      }
#else
      av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                             scan_order, &quant_param);
#endif  // CONFIG_AV1_HIGHBITDEPTH
      av1_inverse_transform_block(xd, dqcoeff, 0, DCT_DCT, tx_size, dst_buffer,
                                  dst_buffer_stride, eob, 0);
      WeberStats *weber_stats =
          &cpi->mb_weber_stats[(mi_row / mb_step) * cpi->frame_info.mi_cols +
                               (mi_col / mb_step)];

      weber_stats->rec_pix_max = 1;
      weber_stats->rec_variance = 0;
      weber_stats->src_pix_max = 1;
      weber_stats->src_variance = 0;
      weber_stats->distortion = 0;

      int64_t src_mean = 0;
      int64_t rec_mean = 0;
      int64_t dist_mean = 0;

      for (int pix_row = 0; pix_row < block_size; ++pix_row) {
        for (int pix_col = 0; pix_col < block_size; ++pix_col) {
          int src_pix, rec_pix;
#if CONFIG_AV1_HIGHBITDEPTH
          if (is_cur_buf_hbd(xd)) {
            uint16_t *src = CONVERT_TO_SHORTPTR(mb_buffer);
            uint16_t *rec = CONVERT_TO_SHORTPTR(dst_buffer);
            src_pix = src[pix_row * buf_stride + pix_col];
            rec_pix = rec[pix_row * dst_buffer_stride + pix_col];
          } else {
            src_pix = mb_buffer[pix_row * buf_stride + pix_col];
            rec_pix = dst_buffer[pix_row * dst_buffer_stride + pix_col];
          }
#else
          src_pix = mb_buffer[pix_row * buf_stride + pix_col];
          rec_pix = dst_buffer[pix_row * dst_buffer_stride + pix_col];
#endif
          src_mean += src_pix;
          rec_mean += rec_pix;
          dist_mean += src_pix - rec_pix;
          weber_stats->src_variance += src_pix * src_pix;
          weber_stats->rec_variance += rec_pix * rec_pix;
          weber_stats->src_pix_max = AOMMAX(weber_stats->src_pix_max, src_pix);
          weber_stats->rec_pix_max = AOMMAX(weber_stats->rec_pix_max, rec_pix);
          weber_stats->distortion += (src_pix - rec_pix) * (src_pix - rec_pix);
        }
      }

      weber_stats->src_variance -= (src_mean * src_mean) / pix_num;
      weber_stats->rec_variance -= (rec_mean * rec_mean) / pix_num;
      weber_stats->distortion -= (dist_mean * dist_mean) / pix_num;
      weber_stats->satd = best_intra_cost;

      double reg = sqrt((double)weber_stats->distortion) *
                   sqrt((double)weber_stats->src_pix_max) * 0.1;
      double alpha_den = fabs(weber_stats->rec_pix_max *
                                  sqrt((double)weber_stats->src_variance) -
                              weber_stats->src_pix_max *
                                  sqrt((double)weber_stats->rec_variance)) +
                         reg;
      double alpha_num = ((double)weber_stats->distortion) *
                             sqrt((double)weber_stats->src_variance) *
                             weber_stats->rec_pix_max +
                         reg;

      weber_stats->alpha = AOMMAX(alpha_num, 1.0) / AOMMAX(alpha_den, 1.0);

      qcoeff[0] = 0;
      for (idx = 1; idx < coeff_count; ++idx) qcoeff[idx] = abs(qcoeff[idx]);
      qsort(qcoeff, coeff_count, sizeof(*coeff), qsort_comp);

      weber_stats->max_scale = (double)qcoeff[coeff_count - 1];
    }
  }

  int sb_step = mi_size_wide[cm->seq_params->sb_size];
  double sb_wiener_log = 0;
  double sb_count = 0;

  for (mi_row = 0; mi_row < cm->mi_params.mi_rows; mi_row += sb_step) {
    for (mi_col = 0; mi_col < cm->mi_params.mi_cols; mi_col += sb_step) {
      int sb_wiener_var =
          get_var_perceptual_ai(cpi, cm->seq_params->sb_size, mi_row, mi_col);
      int64_t satd = get_satd(cpi, cm->seq_params->sb_size, mi_row, mi_col);
      int64_t sse = get_sse(cpi, cm->seq_params->sb_size, mi_row, mi_col);
      double scaled_satd = (double)satd / sqrt((double)sse);
      sb_wiener_log += scaled_satd * log(sb_wiener_var);
      sb_count += scaled_satd;
    }
  }

  if (sb_count > 0)
    cpi->norm_wiener_variance = (int64_t)(exp(sb_wiener_log / sb_count));
  cpi->norm_wiener_variance = AOMMAX(1, cpi->norm_wiener_variance);

  for (int its_cnt = 0; its_cnt < 2; ++its_cnt) {
    sb_wiener_log = 0;
    sb_count = 0;
    for (mi_row = 0; mi_row < cm->mi_params.mi_rows; mi_row += sb_step) {
      for (mi_col = 0; mi_col < cm->mi_params.mi_cols; mi_col += sb_step) {
        int sb_wiener_var =
            get_var_perceptual_ai(cpi, cm->seq_params->sb_size, mi_row, mi_col);

        double beta = (double)cpi->norm_wiener_variance / sb_wiener_var;
        double min_max_scale =
            AOMMAX(1.0, get_max_scale(cpi, bsize, mi_row, mi_col));
        beta = 1.0 / AOMMIN(1.0 / beta, min_max_scale);
        beta = AOMMIN(beta, 4);
        beta = AOMMAX(beta, 0.25);

        sb_wiener_var = (int)(cpi->norm_wiener_variance / beta);

        int64_t satd = get_satd(cpi, cm->seq_params->sb_size, mi_row, mi_col);
        int64_t sse = get_sse(cpi, cm->seq_params->sb_size, mi_row, mi_col);
        double scaled_satd = (double)satd / sqrt((double)sse);
        sb_wiener_log += scaled_satd * log(sb_wiener_var);
        sb_count += scaled_satd;
      }
    }

    if (sb_count > 0)
      cpi->norm_wiener_variance = (int64_t)(exp(sb_wiener_log / sb_count));
    cpi->norm_wiener_variance = AOMMAX(1, cpi->norm_wiener_variance);
  }

  aom_free_frame_buffer(&cm->cur_frame->buf);
}

int av1_get_sbq_perceptual_ai(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                              int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int base_qindex = cm->quant_params.base_qindex;
  int sb_wiener_var = get_var_perceptual_ai(cpi, bsize, mi_row, mi_col);
  int offset = 0;
  double beta = (double)cpi->norm_wiener_variance / sb_wiener_var;
  double min_max_scale = AOMMAX(1.0, get_max_scale(cpi, bsize, mi_row, mi_col));
  beta = 1.0 / AOMMIN(1.0 / beta, min_max_scale);

  // Cap beta such that the delta q value is not much far away from the base q.
  beta = AOMMIN(beta, 4);
  beta = AOMMAX(beta, 0.25);
  offset = av1_get_deltaq_offset(cm->seq_params->bit_depth, base_qindex, beta);
  const DeltaQInfo *const delta_q_info = &cm->delta_q_info;
  offset = AOMMIN(offset, delta_q_info->delta_q_res * 20 - 1);
  offset = AOMMAX(offset, -delta_q_info->delta_q_res * 20 + 1);
  int qindex = cm->quant_params.base_qindex + offset;
  qindex = AOMMIN(qindex, MAXQ);
  qindex = AOMMAX(qindex, MINQ);
  if (base_qindex > MINQ) qindex = AOMMAX(qindex, MINQ + 1);

  return qindex;
}

void av1_init_mb_ur_var_buffer(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;

  if (cpi->mb_delta_q) return;

  CHECK_MEM_ERROR(cm, cpi->mb_delta_q,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->mb_delta_q)));
}

void av1_set_mb_ur_variance(AV1_COMP *cpi) {
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  uint8_t *y_buffer = cpi->source->y_buffer;
  const int y_stride = cpi->source->y_stride;
  const int block_size = cpi->common.seq_params->sb_size;

  const int num_mi_w = mi_size_wide[block_size];
  const int num_mi_h = mi_size_high[block_size];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (mi_params->mi_rows + num_mi_h - 1) / num_mi_h;
  const int use_hbd = cpi->source->flags & YV12_FLAG_HIGHBITDEPTH;

  double a = -23.06 * 4.0, b = 0.004065, c = 30.516 * 4.0;
  int delta_q_avg = 0;
  // Loop through each SB block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double var = 0.0, num_of_var = 0.0;
      const int index = row * num_cols + col;
      double sum_var = 0.0;

      // Loop through each 8x8 block.
      for (int mi_row = row * num_mi_h;
           mi_row < mi_params->mi_rows && mi_row < (row + 1) * num_mi_h;
           mi_row += 2) {
        for (int mi_col = col * num_mi_w;
             mi_col < mi_params->mi_cols && mi_col < (col + 1) * num_mi_w;
             mi_col += 2) {
          struct buf_2d buf;
          const int row_offset_y = mi_row << 2;
          const int col_offset_y = mi_col << 2;

          buf.buf = y_buffer + row_offset_y * y_stride + col_offset_y;
          buf.stride = y_stride;

          double block_variance;
          if (use_hbd) {
            block_variance = av1_high_get_sby_perpixel_variance(
                cpi, &buf, BLOCK_8X8, xd->bd);
          } else {
            block_variance =
                av1_get_sby_perpixel_variance(cpi, &buf, BLOCK_8X8);
          }

          block_variance = block_variance < 1.0 ? 1.0 : block_variance;
          var += log(block_variance);
          sum_var += block_variance;
          num_of_var += 1.0;
        }
      }
      var = exp(var / num_of_var);
      cpi->mb_variance[index] = (int)(sum_var / num_of_var);
      cpi->mb_delta_q[index] = (int)(a * exp(-b * var) + c + 0.5);
      delta_q_avg += cpi->mb_delta_q[index];
    }
  }

  delta_q_avg = (int)((double)delta_q_avg / (num_rows * num_cols) + 0.5);

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      cpi->mb_delta_q[index] -= delta_q_avg;
    }
  }
}

int av1_get_sbq_user_rating_based(AV1_COMP *const cpi, int mi_row, int mi_col) {
  const BLOCK_SIZE bsize = cpi->common.seq_params->sb_size;
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  AV1_COMMON *const cm = &cpi->common;
  const int base_qindex = cm->quant_params.base_qindex;
  if (base_qindex == MINQ || base_qindex == MAXQ) return base_qindex;

  const int num_mi_w = mi_size_wide[bsize];
  const int num_mi_h = mi_size_high[bsize];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int index = (mi_row / num_mi_h) * num_cols + (mi_col / num_mi_w);
  const int delta_q = cpi->mb_delta_q[index];

  int qindex = base_qindex + delta_q;
  qindex = AOMMIN(qindex, MAXQ);
  qindex = AOMMAX(qindex, MINQ + 1);

  return qindex;
}

void av1_init_sb_qindex_buffer(AV1_COMP *cpi) {
  if (cpi->sb_qindex) return;
  if (cpi->mb_variance) return;

  AV1_COMMON *cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int sb_width = mi_size_wide[sb_size];
  const int sb_height = mi_size_high[sb_size];
  const int num_col_sbs = (mi_params->mi_cols + sb_width - 1) / sb_width;
  const int num_row_sbs = (mi_params->mi_rows + sb_height - 1) / sb_height;
  const int num_sbs = num_col_sbs * num_row_sbs;

  CHECK_MEM_ERROR(cm, cpi->sb_qindex,
                  aom_calloc(num_sbs, sizeof(*cpi->sb_qindex)));

  if (cpi->est_best_mode) return;

  CHECK_MEM_ERROR(cm, cpi->est_best_mode,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->est_best_mode)));

  CHECK_MEM_ERROR(cm, cpi->mb_variance,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->mb_variance)));
}

static double compute_avg_variance(AV1_COMP *cpi) {
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const int block_size = cpi->common.seq_params->sb_size;
  const int num_mi_w = mi_size_wide[block_size];
  const int num_mi_h = mi_size_high[block_size];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (mi_params->mi_rows + num_mi_h - 1) / num_mi_h;

  double sum_var = 0;
  int count = 0;
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      // Note: cpi->mb_variance is the averaged var of a super block
      // in the unit of 8x8.
      int num_8x8 = 0;
      for (int mi_row = row * num_mi_h;
           mi_row < mi_params->mi_rows && mi_row < (row + 1) * num_mi_h;
           mi_row += 2) {
        for (int mi_col = col * num_mi_w;
             mi_col < mi_params->mi_cols && mi_col < (col + 1) * num_mi_w;
             mi_col += 2) {
          ++num_8x8;
        }
      }
      sum_var += cpi->mb_variance[index] * num_8x8;
      count += num_8x8;
    }
  }
  const double avg_var = sum_var / count;
  return avg_var;
}

static void median_filt_frame(const uint8_t *buffer, const int buf_stride,
                              const int frame_width, const int frame_height,
                              uint8_t *output_buffer) {
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      int vals[25];
      int count = 0;
      // Get neighbor 5x5 pixels
      for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
          const int r = row + dy;
          const int c = col + dx;
          if (r >= 0 && r < frame_height && c >= 0 && c < frame_width) {
            vals[count] = buffer[r * buf_stride + c];
            ++count;
          }
        }
      }
      // bubble sort
      for (int i = 0; i < count; ++i) {
        for (int j = i; j < count; ++j) {
          if (vals[j] > vals[i]) {
            const int tmp = vals[j];
            vals[j] = vals[i];
            vals[i] = tmp;
          }
        }
      }
      output_buffer[row * buf_stride + col] = vals[count / 2];
    }
  }
}

static int calc_histogram(const uint8_t *buffer, const int buf_stride,
                          const int frame_width, const int frame_height,
                          int *histogram) {
  int num_unique_colors = 0;
  const uint8_t *buf = buffer;
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      const int val = buf[col];
      if (histogram[val] == 0) ++num_unique_colors;
      ++histogram[val];
    }
    buf += buf_stride;
  }
  return num_unique_colors;
}

static int saliency_dist(int x, int y) {
  return abs(x - y);
  // return (x - y) * (x - y);
}

static void saliency_smoothing(double saliency[256], int window) {
  double orig_saliency[256];
  for (int i = 0; i < 256; ++i) orig_saliency[i] = saliency[i];
  const int half_win = window / 2;
  for (int i = 0; i < 256; ++i) {
    if (orig_saliency[i] == 0) continue;
    int sum_dist = 0;
    int count = 0;
    for (int j = -half_win; j <= half_win; ++j) {
      if (i + j >= 0 && i + j < 256 && orig_saliency[i + j] > 0) {
        sum_dist += saliency_dist(i, i + j);
        ++count;
      }
    }
    saliency[i] = 0;
    for (int j = -half_win; j <= half_win; ++j) {
      if (i + j >= 0 && i + j < 256 && orig_saliency[i + j] > 0) {
        saliency[i] +=
            (sum_dist - saliency_dist(i, i + j)) * orig_saliency[i + j];
      }
    }
    saliency[i] /= ((count - 1) * sum_dist);
  }
}

static double build_saliency_map(AV1_COMP *cpi, double *saliency_map) {
  const uint8_t *buffer = cpi->source->y_buffer;
  const int buf_stride = cpi->source->y_stride;
  const int frame_width = cpi->frame_info.frame_width;
  const int frame_height = cpi->frame_info.frame_height;
  uint8_t *filtered_buffer = cpi->ppi->filtered_buffer.y_buffer;

  // Filter the source
  median_filt_frame(buffer, buf_stride, frame_width, frame_height,
                    filtered_buffer);

  // Compute Y channel color histogram
  // Here we assume it is 8-bit, therefore 256 colors.
  int histogram[256] = { 0 };
  double saliency[256] = { 0 };
  calc_histogram(filtered_buffer, buf_stride, frame_width, frame_height,
                 histogram);

  // Calculate saliency
  for (int i = 0; i < 256; ++i) {
    if (histogram[i] == 0) continue;
    for (int j = 0; j < 256; ++j) {
      if (i == j || histogram[j] == 0) continue;
      saliency[i] += saliency_dist(i, j) * histogram[j];
    }
  }

  // Normalize
  double max_saliency = 0;
  for (int i = 0; i < 256; ++i) {
    max_saliency = AOMMAX(max_saliency, saliency[i]);
  }
  for (int i = 0; i < 256; ++i) {
    if (saliency[i] > 0) {
      saliency[i] /= max_saliency;
    }
  }

  // Color space smoothing
  saliency_smoothing(saliency, /*window=*/8);

  // Assign saliency map
  double sum_saliency = 0;
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      const int val = filtered_buffer[row * buf_stride + col];
      saliency_map[row * frame_width + col] = saliency[val];
      sum_saliency += saliency[val];
    }
  }
  const double avg_saliency = sum_saliency / (frame_width * frame_height);
  return avg_saliency;

  /*
  // Write to file
  FILE *pfile = fopen("saliency_map.csv", "w");
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      fprintf(pfile, "%d", saliency_map[row * frame_width + col]);
      if (col < frame_width - 1) fprintf(pfile, ",");
    }
    fprintf(pfile, "\n");
  }
  fclose(pfile);

  pfile = fopen("histogram.csv", "w");
  for (int i = 0; i < 256; ++i) {
    fprintf(pfile, "%d", histogram[i]);
    if (i < 255) fprintf(pfile, ",");
  }
  fclose(pfile);
  */

  // printf("saliency map generated.\n");
}

static double build_mscn_map(AV1_COMP *cpi, double *mscn_map) {
  const uint8_t *buffer = cpi->source->y_buffer;
  const int buf_stride = cpi->source->y_stride;
  const int frame_width = cpi->frame_info.frame_width;
  const int frame_height = cpi->frame_info.frame_height;
  const int half_win = 3;

  /*
  // Generate mscn map
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      double sum = 0;
      double sum_square = 0;
      int count = 0;
      for (int dy = -half_win; dy <= half_win; ++dy) {
        for (int dx = -half_win; dx <= half_win; ++dx) {
          if (row + dy < 0 || row + dy >= frame_height ||
              col + dx < 0 || col + dx >= frame_width) {
            continue;
          }
          const int pix = buffer[(row + dy) * buf_stride + col + dx];
          sum += pix;
          sum_square += pix * pix;
          ++count;
        }
      }
      const double mean = sum / count;
      const double sigma = sqrt(sum_square / count - mean * mean);
      mscn_map[row * frame_width + col] =
          (buffer[row * buf_stride + col] - mean) / (sigma + 1.0);
    }
  }
  */

  // h = round(fspecial('gaussian', 7, 3.0) * 1000)
  const int gauss_kernel[] = { 11, 15, 18, 19, 18, 15, 11, 15, 20, 23,
                               25, 23, 20, 15, 18, 23, 27, 29, 27, 23,
                               18, 19, 25, 29, 31, 29, 25, 19, 18, 23,
                               27, 29, 27, 23, 18, 15, 20, 23, 25, 23,
                               20, 15, 11, 15, 18, 19, 18, 15, 11 };
  // const int kernel_sum = 1003;

  // Generate mscn map with Gaussian kernel weights.
  double *mean_map = aom_calloc(frame_width * frame_height, sizeof(*mean_map));
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      double weighted_sum = 0;
      int count = 0;
      for (int dy = -half_win; dy <= half_win; ++dy) {
        for (int dx = -half_win; dx <= half_win; ++dx) {
          if (row + dy < 0 || row + dy >= frame_height || col + dx < 0 ||
              col + dx >= frame_width) {
            continue;
          }
          const int pix = buffer[(row + dy) * buf_stride + col + dx];
          weighted_sum +=
              pix * gauss_kernel[(dy + half_win) * (2 * half_win + 1) +
                                 (dx + half_win)];
          count += gauss_kernel[(dy + half_win) * (2 * half_win + 1) +
                                (dx + half_win)];
        }
      }
      const double weighted_mean = weighted_sum / count;
      mean_map[row * frame_width + col] = weighted_mean;
    }
  }
  double sum_mscn = 0;
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      double weighted_sum = 0;
      double count = 0;
      const double mean = mean_map[row * frame_width + col];
      for (int dy = -half_win; dy <= half_win; ++dy) {
        for (int dx = -half_win; dx <= half_win; ++dx) {
          if (row + dy < 0 || row + dy >= frame_height || col + dx < 0 ||
              col + dx >= frame_width) {
            continue;
          }
          const int pix = buffer[(row + dy) * buf_stride + col + dx];
          const double weight =
              gauss_kernel[(dy + half_win) * (2 * half_win + 1) +
                           (dx + half_win)];
          weighted_sum += weight * (pix - mean) * (pix - mean);
          count += weight;
        }
      }
      const double sigma = sqrt(weighted_sum / count);
      mscn_map[row * frame_width + col] =
          (buffer[row * buf_stride + col] - mean) / (sigma + 1.0);
      sum_mscn += fabs(mscn_map[row * frame_width + col]);
    }
  }
  aom_free(mean_map);
  const double avg_mscn = sum_mscn / (frame_width * frame_height);
  return avg_mscn;

  /*
  // Write to file
  FILE *pfile = fopen("mscn_map.csv", "w");
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      fprintf(pfile, "%.4f", mscn_map[row * frame_width + col]);
      if (col < frame_width - 1) fprintf(pfile, ",");
    }
    fprintf(pfile, "\n");
  }
  fclose(pfile);

  printf("mscn map generated.\n");
  */
}

void av1_set_sb_qindex(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  uint8_t *buffer = cpi->source->y_buffer;
  int buf_stride = cpi->source->y_stride;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  MB_MODE_INFO mbmi;
  memset(&mbmi, 0, sizeof(mbmi));
  MB_MODE_INFO *mbmi_ptr = &mbmi;
  xd->mi = &mbmi_ptr;
  xd->cur_buf = cpi->source;

  const SequenceHeader *const seq_params = cm->seq_params;
  if (aom_realloc_frame_buffer(
          &cm->cur_frame->buf, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, cm->features.byte_alignment, NULL, NULL,
          NULL, cpi->oxcf.tool_cfg.enable_global_motion))
    aom_internal_error(cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate frame buffer");

  cm->quant_params.base_qindex = cpi->oxcf.rc_cfg.cq_level;
  av1_frame_init_quantizer(cpi);

  DECLARE_ALIGNED(32, int16_t, src_diff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, coeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, qcoeff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, dqcoeff[32 * 32]);

  int mb_row, mb_col, count = 0;
  const TX_SIZE tx_size = TX_16X16;
  const int block_size = tx_size_wide[tx_size];
  const int coeff_count = block_size * block_size;

  const BitDepthInfo bd_info = get_bit_depth_info(xd);
  cpi->norm_wiener_variance = 0;

  int mb_step = mi_size_wide[BLOCK_16X16];
  BLOCK_SIZE bsize = BLOCK_16X16;
  double abs_coeff_mean[16 * 16] = { 0 };
  const int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];

  // --------------- First compute mean qcoeff and store best mode -----------
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      PREDICTION_MODE best_mode = DC_PRED;
      int best_intra_cost = INT_MAX;
      int mi_row = mb_row * mb_step;
      int mi_col = mb_col * mb_step;
      xd->up_available = mi_row > 0;
      xd->left_available = mi_col > 0;
      const int mi_width = mi_size_wide[bsize];
      const int mi_height = mi_size_high[bsize];
      set_mode_info_offsets(&cpi->common.mi_params, &cpi->mbmi_ext_info, x, xd,
                            mi_row, mi_col);
      set_mi_row_col(xd, &xd->tile, mi_row, mi_height, mi_col, mi_width,
                     cm->mi_params.mi_rows, cm->mi_params.mi_cols);
      set_plane_n4(xd, mi_size_wide[bsize], mi_size_high[bsize],
                   av1_num_planes(cm));
      xd->mi[0]->bsize = bsize;
      xd->mi[0]->motion_mode = SIMPLE_TRANSLATION;
      av1_setup_dst_planes(xd->plane, bsize, &cm->cur_frame->buf, mi_row,
                           mi_col, 0, av1_num_planes(cm));
      int dst_buffer_stride = xd->plane[0].dst.stride;
      uint8_t *dst_buffer = xd->plane[0].dst.buf;
      uint8_t *mb_buffer =
          buffer + mi_row * MI_SIZE * buf_stride + mi_col * MI_SIZE;
      for (PREDICTION_MODE mode = INTRA_MODE_START; mode < INTRA_MODE_END;
           ++mode) {
        av1_predict_intra_block(
            xd, cm->seq_params->sb_size,
            cm->seq_params->enable_intra_edge_filter, block_size, block_size,
            tx_size, mode, 0, 0, FILTER_INTRA_MODES, dst_buffer,
            dst_buffer_stride, dst_buffer, dst_buffer_stride, 0, 0, 0);

        av1_subtract_block(bd_info, block_size, block_size, src_diff,
                           block_size, mb_buffer, buf_stride, dst_buffer,
                           dst_buffer_stride);
        av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);
        int intra_cost = aom_satd(coeff, coeff_count);
        if (intra_cost < best_intra_cost) {
          best_intra_cost = intra_cost;
          best_mode = mode;
        }
      }

      av1_predict_intra_block(xd, cm->seq_params->sb_size,
                              cm->seq_params->enable_intra_edge_filter,
                              block_size, block_size, tx_size, best_mode, 0, 0,
                              FILTER_INTRA_MODES, dst_buffer, dst_buffer_stride,
                              dst_buffer, dst_buffer_stride, 0, 0, 0);
      av1_subtract_block(bd_info, block_size, block_size, src_diff, block_size,
                         mb_buffer, buf_stride, dst_buffer, dst_buffer_stride);
      av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);

      const struct macroblock_plane *const p = &x->plane[0];
      uint16_t eob;
      const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];
      QUANT_PARAM quant_param;
      av1_setup_quant(tx_size, 0, AV1_XFORM_QUANT_FP, 0, &quant_param);
#if CONFIG_AV1_HIGHBITDEPTH
      if (is_cur_buf_hbd(xd)) {
        av1_highbd_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                                      scan_order, &quant_param);
      } else {
        av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                               scan_order, &quant_param);
      }
#else
      av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                             scan_order, &quant_param);
#endif  // CONFIG_AV1_HIGHBITDEPTH
      av1_inverse_transform_block(xd, dqcoeff, 0, DCT_DCT, tx_size, dst_buffer,
                                  dst_buffer_stride, eob, 0);

      for (int i = 1; i < pix_num; ++i) {
        abs_coeff_mean[i] += abs(qcoeff[i]);
      }
      cpi->est_best_mode[mb_row * cpi->frame_info.mb_cols + mb_col] = best_mode;
      ++count;
    }
  }

  for (int i = 1; i < pix_num; ++i) {
    abs_coeff_mean[i] /= count;
  }
  count = 0;
  double sum_frame_entropy = 0;
  double sum_frame_rate = 0;
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      const PREDICTION_MODE best_mode =
          cpi->est_best_mode[mb_row * cpi->frame_info.mb_cols + mb_col];
      int mi_row = mb_row * mb_step;
      int mi_col = mb_col * mb_step;
      xd->up_available = mi_row > 0;
      xd->left_available = mi_col > 0;
      const int mi_width = mi_size_wide[bsize];
      const int mi_height = mi_size_high[bsize];
      set_mode_info_offsets(&cpi->common.mi_params, &cpi->mbmi_ext_info, x, xd,
                            mi_row, mi_col);
      set_mi_row_col(xd, &xd->tile, mi_row, mi_height, mi_col, mi_width,
                     cm->mi_params.mi_rows, cm->mi_params.mi_cols);
      set_plane_n4(xd, mi_size_wide[bsize], mi_size_high[bsize],
                   av1_num_planes(cm));
      xd->mi[0]->bsize = bsize;
      xd->mi[0]->motion_mode = SIMPLE_TRANSLATION;
      av1_setup_dst_planes(xd->plane, bsize, &cm->cur_frame->buf, mi_row,
                           mi_col, 0, av1_num_planes(cm));
      int dst_buffer_stride = xd->plane[0].dst.stride;
      uint8_t *dst_buffer = xd->plane[0].dst.buf;
      uint8_t *mb_buffer =
          buffer + mi_row * MI_SIZE * buf_stride + mi_col * MI_SIZE;

      av1_predict_intra_block(xd, cm->seq_params->sb_size,
                              cm->seq_params->enable_intra_edge_filter,
                              block_size, block_size, tx_size, best_mode, 0, 0,
                              FILTER_INTRA_MODES, dst_buffer, dst_buffer_stride,
                              dst_buffer, dst_buffer_stride, 0, 0, 0);
      av1_subtract_block(bd_info, block_size, block_size, src_diff, block_size,
                         mb_buffer, buf_stride, dst_buffer, dst_buffer_stride);
      av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);

      const struct macroblock_plane *const p = &x->plane[0];
      uint16_t eob;
      const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];
      QUANT_PARAM quant_param;
      // int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
      av1_setup_quant(tx_size, 0, AV1_XFORM_QUANT_FP, 0, &quant_param);
#if CONFIG_AV1_HIGHBITDEPTH
      if (is_cur_buf_hbd(xd)) {
        av1_highbd_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                                      scan_order, &quant_param);
      } else {
        av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                               scan_order, &quant_param);
      }
#else
      av1_quantize_fp_facade(coeff, pix_num, p, qcoeff, dqcoeff, &eob,
                             scan_order, &quant_param);
#endif  // CONFIG_AV1_HIGHBITDEPTH
      av1_inverse_transform_block(xd, dqcoeff, 0, DCT_DCT, tx_size, dst_buffer,
                                  dst_buffer_stride, eob, 0);

      const double entropy = estimate_txfm_block_entropy(
          cm->quant_params.base_qindex, abs_coeff_mean, qcoeff, pix_num);
      const int rate = rate_estimator(qcoeff, eob, tx_size);
      sum_frame_entropy += entropy;
      sum_frame_rate += rate;
      /*
      {
        FILE *pfile = fopen("entropy_variance.stat", "a");
        fprintf(pfile, "%.3f,%d,", entropy, rate);
        fclose(pfile);
      }
      */
      ++count;
    }
  }

  av1_set_mb_ur_variance(cpi);
  const int frame_width = cpi->frame_info.frame_width;
  const int frame_height = cpi->frame_info.frame_height;
  const double avg_var = compute_avg_variance(cpi);
  double *saliency_map =
      aom_calloc(frame_width * frame_height, sizeof(*saliency_map));
  const double avg_saliency = build_saliency_map(cpi, saliency_map);
  double *mscn_map = aom_calloc(frame_width * frame_height, sizeof(*mscn_map));
  const double avg_mscn = build_mscn_map(cpi, mscn_map);

  // determine qindex for each superblock based on even entropy distribution
  entropy_based_sb_q(cpi, abs_coeff_mean, sum_frame_entropy, sum_frame_rate,
                     avg_var, saliency_map, avg_saliency, mscn_map, avg_mscn);

  aom_free_frame_buffer(&cm->cur_frame->buf);
  aom_free(saliency_map);
  aom_free(mscn_map);
}

int av1_get_sb_qindex(AV1_COMP *const cpi, int mi_row, int mi_col) {
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int sb_row = mi_row / mi_size_high[sb_size];
  const int sb_col = mi_col / mi_size_wide[sb_size];
  const int num_col_sbs =
      (mi_params->mi_cols + mi_size_wide[sb_size] - 1) / mi_size_wide[sb_size];
  const int sb_idx = sb_row * num_col_sbs + sb_col;
  return cpi->sb_qindex[sb_idx];
}
