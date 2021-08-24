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

  if (cpi->mb_weber_stats) return;

  CHECK_MEM_ERROR(cm, cpi->mb_weber_stats,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->mb_weber_stats)));

  // TODO(chengchen):
  // The number allocated here is more than enough
  CHECK_MEM_ERROR(cm, cpi->sb_qindex,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->sb_qindex)));
}

static int64_t get_satd(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                        int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int mi_wide = mi_size_wide[bsize];
  const int mi_high = mi_size_high[bsize];

  const int mi_step = mi_size_wide[BLOCK_16X16];
  int mb_stride = cpi->frame_info.mb_cols;
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

  const int mi_step = mi_size_wide[BLOCK_16X16];
  int mb_stride = cpi->frame_info.mb_cols;
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
  const int mi_step = mi_size_wide[BLOCK_16X16];
  int mb_stride = cpi->frame_info.mb_cols;
  double min_max_scale = 10.0;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;
      WeberStats *weber_stats =
          &cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)];
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

  const int mi_step = mi_size_wide[BLOCK_16X16];
  int sb_wiener_var = 0;
  int mb_stride = cpi->frame_info.mb_cols;
  int mb_count = 0;
  int64_t mb_wiener_sum = 0;
  double base_num = 1;
  double base_den = 1;
  double base_reg = 1;

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;

      WeberStats *weber_stats =
          &cpi->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)];
      mb_wiener_sum += (int)(cpi->mb_weber_stats[(row / mi_step) * mb_stride +
                                                 (col / mi_step)]
                                 .alpha *
                             10000);

      base_num += ((double)weber_stats->distortion) *
                  // weber_stats->entropy *
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
                          const double avg_sb_entropy, const double sb_entropy,
                          double upper_sb_entropy, double lower_sb_entropy,
                          double *orig_entropy, const int iter) {
  (void)avg_sb_entropy;
  if (iter == 0) return;
  if (iter == 1) *orig_entropy = sb_entropy;

  // Handle cases when some super blocks are much more complex than the average
  if (*orig_entropy > avg_sb_entropy * 5) {
    const int delta_qindex = 10;
    if (sb_entropy > avg_sb_entropy * 2) {
      *base_qindex = AOMMIN(last_qindex + delta_qindex, MAXQ);
    }
  } else if (*orig_entropy < avg_sb_entropy / 10) {
    const int delta_qindex = 15;
    if (sb_entropy < avg_sb_entropy / 2) {
      *base_qindex = AOMMAX(0, last_qindex - delta_qindex);
    }
  } else {
    const int delta_qindex = 4;
    if (sb_entropy > upper_sb_entropy) {
      *base_qindex = AOMMIN(last_qindex + delta_qindex, MAXQ);
    } else if (sb_entropy < lower_sb_entropy) {
      *base_qindex = AOMMAX(0, last_qindex - delta_qindex);
    }
  }
}

static void entropy_based_sb_q(AV1_COMP *cpi, double *abs_coeff_mean,
                               const double sum_frame_entropy,
                               const double sum_frame_rate) {
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
  const int coeff_count = block_size * block_size;
  const BitDepthInfo bd_info = get_bit_depth_info(xd);

  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int sb_width = mi_size_wide[sb_size];
  const int sb_height = mi_size_high[sb_size];
  const BLOCK_SIZE bsize = BLOCK_16X16;
  const int mb_step = mi_size_wide[bsize];
  const int num_col_sbs = (mi_params->mi_cols + sb_width - 1) / sb_width;
  const int num_row_sbs = (mi_params->mi_rows + sb_height - 1) / sb_height;
  const int num_sbs = num_col_sbs * num_row_sbs;
  const int orig_base_qindex = cpi->common.quant_params.base_qindex;
  const double avg_sb_entropy = sum_frame_entropy / num_sbs;
  const double upper_sb_entropy = 2.5 * avg_sb_entropy;
  const double lower_sb_entropy = 0.75 * avg_sb_entropy;
  const int max_iter = 5;
  int base_qindex = orig_base_qindex;

  // Iterate though each super block
  for (int sb_row = 0; sb_row < num_row_sbs; ++sb_row) {
    for (int sb_col = 0; sb_col < num_col_sbs; ++sb_col) {
      base_qindex = orig_base_qindex;
      double sb_entropy = 0;
      int iter = 0;
      int last_qindex = base_qindex;
      double orig_entropy = 0;
      do {
        // Adjust qindex based on the last qindex and corresponding sb entropy
        adjust_qindex(&base_qindex, last_qindex, avg_sb_entropy, sb_entropy,
                      upper_sb_entropy, lower_sb_entropy, &orig_entropy, iter);
        last_qindex = base_qindex;
        cm->quant_params.base_qindex = base_qindex;
        av1_frame_init_quantizer(cpi);
        ++iter;
        sb_entropy = 0;

        // Iterater though each mb block (16x16) inside a super block
        for (int mb_row = 0; mb_row < sb_height / mb_step; ++mb_row) {
          for (int mb_col = 0; mb_col < sb_width / mb_step; ++mb_col) {
            PREDICTION_MODE best_mode = DC_PRED;
            int best_intra_cost = INT_MAX;
            const int mi_row = sb_row * sb_height + mb_row * mb_step;
            const int mi_col = sb_col * sb_width + mb_col * mb_step;
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
            for (PREDICTION_MODE mode = INTRA_MODE_START; mode < INTRA_MODE_END;
                 ++mode) {
              av1_predict_intra_block(xd, cm->seq_params->sb_size,
                                      cm->seq_params->enable_intra_edge_filter,
                                      block_size, block_size, tx_size, mode, 0,
                                      0, FILTER_INTRA_MODES, dst_buffer,
                                      dst_buffer_stride, dst_buffer,
                                      dst_buffer_stride, 0, 0, 0);
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
            sb_entropy += entropy;
          }
        }
        /*
        if (sb_row == 1 && sb_col == 12) {
          printf("iter %d, base_qindex %d, avg_sb_entropy %.0f, sb_entropy
        %.0f\n", iter, base_qindex, avg_sb_entropy, sb_entropy);
        }
        */
      } while (iter < max_iter && (sb_entropy > upper_sb_entropy ||
                                   sb_entropy < lower_sb_entropy));
      /*
      if (iter == 1) orig_entropy = sb_entropy;
      FILE *pfile = fopen("entropy.stat", "a");
      fprintf(pfile, "(%d, %d), avg_sb_entropy %.0f, orig_entropy %.0f, entropy
      %.0f, ratio %.2f, base_qindex %d\n", sb_row, sb_col, avg_sb_entropy,
      orig_entropy, sb_entropy, sb_entropy / orig_entropy, base_qindex);
      fclose(pfile);
      */

      const int sb_idx = sb_row * num_col_sbs + sb_col;
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

  int mb_row, mb_col, count = 0;
  const TX_SIZE tx_size = TX_16X16;
  const int block_size = tx_size_wide[tx_size];
  const int coeff_count = block_size * block_size;

  const BitDepthInfo bd_info = get_bit_depth_info(xd);
  cpi->norm_wiener_variance = 0;

  int mb_step = mi_size_wide[BLOCK_16X16];
  BLOCK_SIZE bsize = BLOCK_16X16;
  double abs_coeff_mean[16 * 16] = { 0 };

  // --------------- First compute mean qcoeff -----------------------
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

      int idx;
      int16_t median_val = 0;
      int64_t wiener_variance = 0;
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
          &cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col];

      {
        for (int i = 1; i < pix_num; ++i) {
          abs_coeff_mean[i] += abs(qcoeff[i]);
        }
      }

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

      coeff[0] = 0;
      for (idx = 1; idx < coeff_count; ++idx) coeff[idx] = abs(coeff[idx]);
      qsort(coeff, coeff_count - 1, sizeof(*coeff), qsort_comp);

      // Noise level estimation
      median_val = coeff[coeff_count / 2];

      // Wiener filter
      for (idx = 1; idx < coeff_count; ++idx) {
        int64_t sqr_coeff = (int64_t)coeff[idx] * coeff[idx];
        int64_t tmp_coeff = (int64_t)coeff[idx];
        if (median_val) {
          tmp_coeff = (sqr_coeff * coeff[idx]) /
                      (sqr_coeff + (int64_t)median_val * median_val);
        }
        wiener_variance += tmp_coeff * tmp_coeff;
      }
      cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
          .mb_wiener_variance = wiener_variance / coeff_count;
      ++count;
    }
  }

  // --------------- Then repeat -----------------------
  const int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
  for (int i = 1; i < pix_num; ++i) {
    abs_coeff_mean[i] /= count;
  }
  count = 0;
  double sum_frame_entropy = 0;
  double sum_frame_rate = 0;
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

      int idx;
      int16_t median_val = 0;
      int64_t wiener_variance = 0;
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
      WeberStats *weber_stats =
          &cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col];

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
      weber_stats->entropy = entropy;

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

      coeff[0] = 0;
      for (idx = 1; idx < coeff_count; ++idx) coeff[idx] = abs(coeff[idx]);
      qsort(coeff, coeff_count - 1, sizeof(*coeff), qsort_comp);

      // Noise level estimation
      median_val = coeff[coeff_count / 2];

      // Wiener filter
      for (idx = 1; idx < coeff_count; ++idx) {
        int64_t sqr_coeff = (int64_t)coeff[idx] * coeff[idx];
        int64_t tmp_coeff = (int64_t)coeff[idx];
        if (median_val) {
          tmp_coeff = (sqr_coeff * coeff[idx]) /
                      (sqr_coeff + (int64_t)median_val * median_val);
        }
        wiener_variance += tmp_coeff * tmp_coeff;
      }
      cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
          .mb_wiener_variance = wiener_variance / coeff_count;
      // printf("wiener variance %ld\n",
      //        cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols +
      //        mb_col].mb_wiener_variance);
      /*
      FILE *pfile = fopen("entropy_variance.stat", "a");
      fprintf(pfile, "%ld,%ld\n", weber_stats->distortion,
              cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols +
      mb_col].mb_wiener_variance); fclose(pfile);
      */
      ++count;
    }
  }
  /*
  FILE *pfile = fopen("mask.csv", "w");
  fprintf(pfile, "%d,%d,%ld\n", cpi->frame_info.frame_height,
          cpi->frame_info.frame_width, cpi->norm_wiener_variance);
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      fprintf(pfile, "%.0f,",
              cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols +
  mb_col].entropy);
    }
    fprintf(pfile, "\n");
  }
  fclose(pfile);
  */
  // determine qindex for each superblock based on even entropy distribution
  entropy_based_sb_q(cpi, abs_coeff_mean, sum_frame_entropy, sum_frame_rate);

  int sb_step = mi_size_wide[cm->seq_params->sb_size];
  double sb_wiener_log = 0;
  double sb_count = 0;

  for (int mi_row = 0; mi_row < cm->mi_params.mi_rows; mi_row += sb_step) {
    for (int mi_col = 0; mi_col < cm->mi_params.mi_cols; mi_col += sb_step) {
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

  sb_wiener_log = 0;
  sb_count = 0;
  for (int mi_row = 0; mi_row < cm->mi_params.mi_rows; mi_row += sb_step) {
    for (int mi_col = 0; mi_col < cm->mi_params.mi_cols; mi_col += sb_step) {
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

  aom_free_frame_buffer(&cm->cur_frame->buf);
}

static int get_sb_qindex(AV1_COMP *const cpi, int mi_row, int mi_col) {
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const BLOCK_SIZE sb_size = cpi->common.seq_params->sb_size;
  const int sb_row = mi_row / mi_size_high[sb_size];
  const int sb_col = mi_col / mi_size_wide[sb_size];
  const int num_col_sbs =
      (mi_params->mi_cols + mi_size_wide[sb_size] - 1) / mi_size_wide[sb_size];
  const int sb_idx = sb_row * num_col_sbs + sb_col;
  return cpi->sb_qindex[sb_idx];
}

int av1_get_sbq_perceptual_ai(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                              int mi_col) {
  return get_sb_qindex(cpi, mi_row, mi_col);

  AV1_COMMON *const cm = &cpi->common;
  const int base_qindex = cm->quant_params.base_qindex;
  int sb_wiener_var = get_var_perceptual_ai(cpi, bsize, mi_row, mi_col);
  int offset = 0;
  double beta = (double)cpi->norm_wiener_variance / sb_wiener_var;
  // printf("orig beta %.3f", beta);
  double min_max_scale = AOMMAX(1.0, get_max_scale(cpi, bsize, mi_row, mi_col));
  beta = 1.0 / AOMMIN(1.0 / beta, min_max_scale);
  // printf(", new beta %.3f, min_max_scale %.3f, ", beta, min_max_scale);

  // Cap beta such that the delta q value is not much far away from the base q.
  beta = AOMMIN(beta, 4);
  beta = AOMMAX(beta, 0.25);
  offset = av1_get_deltaq_offset(cm->seq_params->bit_depth, base_qindex, beta);
  // printf("base_qindex %d, offset %d \n", base_qindex, offset);
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

  if (cpi->mb_variance) return;

  CHECK_MEM_ERROR(cm, cpi->mb_variance,
                  aom_calloc(cpi->frame_info.mb_rows * cpi->frame_info.mb_cols,
                             sizeof(*cpi->mb_variance)));
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

  // Loop through each SB block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double var = 0.0, num_of_var = 0.0;
      const int index = row * num_cols + col;

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

          if (use_hbd) {
            var += av1_high_get_sby_perpixel_variance(cpi, &buf, BLOCK_8X8,
                                                      xd->bd);
          } else {
            var += av1_get_sby_perpixel_variance(cpi, &buf, BLOCK_8X8);
          }

          num_of_var += 1.0;
        }
      }
      var = var / num_of_var;
      cpi->mb_variance[index] = var;
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
  const double var = cpi->mb_variance[index];

  const int beta = 80;
  double a = -23.5 * 4.0, b = 0.00198, c = 30.65 * 4.0;
  if (base_qindex <= beta) {
    const double alpha = (double)base_qindex / (double)beta;
    a *= alpha;
    c *= alpha;
  } else {
    const double alpha = (double)(base_qindex - beta) / (double)(MAXQ - beta);
    a = a - a * alpha;
    c = c + ((double)MAXQ - c) * alpha;
  }

  int qindex = (int)(a * exp(-b * var) + c + 0.5);
  qindex = AOMMIN(qindex, MAXQ);
  qindex = AOMMAX(qindex, MINQ + 1);

  return qindex;
}
