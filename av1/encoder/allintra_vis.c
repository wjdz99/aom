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

#include "av1/encoder/allintra_vis.h"

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

static void build_saliency_map(AV1_COMP *cpi) {
  const uint8_t *buffer = cpi->source->y_buffer;
  const int buf_stride = cpi->source->y_stride;
  const int frame_width = cpi->frame_info.frame_width;
  const int frame_height = cpi->frame_info.frame_height;
  uint8_t *filtered_buffer = cpi->ppi->filtered_buffer.y_buffer;
  uint8_t *saliency_map = cpi->ppi->saliency_buffer.y_buffer;

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
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      const int val = filtered_buffer[row * buf_stride + col];
      saliency_map[row * buf_stride + col] = (uint8_t)(saliency[val] * 255);
    }
  }

  /*
  // Write to file
  FILE *pfile = fopen("saliency_map.csv", "w");
  for (int row = 0; row < frame_height; ++row) {
    for (int col = 0; col < frame_width; ++col) {
      fprintf(pfile, "%d", saliency_map[row * buf_stride + col]);
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

  for (int row = mi_row; row < mi_row + mi_high; row += mi_step) {
    for (int col = mi_col; col < mi_col + mi_wide; col += mi_step) {
      if (row >= cm->mi_params.mi_rows || col >= cm->mi_params.mi_cols)
        continue;

      mb_wiener_sum +=
          (int)cpi
              ->mb_weber_stats[(row / mi_step) * mb_stride + (col / mi_step)]
              .mb_wiener_variance;
      ++mb_count;
    }
  }

  if (mb_count) sb_wiener_var = (int)(mb_wiener_sum / mb_count);
  sb_wiener_var = AOMMAX(1, sb_wiener_var);

  return (int)(sqrt(sb_wiener_var));
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

  build_saliency_map(cpi);
  uint8_t *saliency_map = cpi->ppi->saliency_buffer.y_buffer;

  union {
#if CONFIG_AV1_HIGHBITDEPTH
    DECLARE_ALIGNED(32, uint16_t, zero_pred16[32 * 32]);
#endif
    DECLARE_ALIGNED(32, uint8_t, zero_pred8[32 * 32]);
  } pred_buffer_mem;
  uint8_t *pred_buf;

  DECLARE_ALIGNED(32, int16_t, src_diff[32 * 32]);
  DECLARE_ALIGNED(32, tran_low_t, coeff[32 * 32]);

  int mb_row, mb_col, count = 0;
  const TX_SIZE tx_size = TX_16X16;
  const int block_size = tx_size_wide[tx_size];
  const int coeff_count = block_size * block_size;

#if CONFIG_AV1_HIGHBITDEPTH
  xd->cur_buf = cpi->source;
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    pred_buf = CONVERT_TO_BYTEPTR(pred_buffer_mem.zero_pred16);
    memset(pred_buffer_mem.zero_pred16, 0,
           sizeof(*pred_buffer_mem.zero_pred16) * coeff_count);
  } else {
    pred_buf = pred_buffer_mem.zero_pred8;
    memset(pred_buffer_mem.zero_pred8, 0,
           sizeof(*pred_buffer_mem.zero_pred8) * coeff_count);
  }
#else
  pred_buf = pred_buffer_mem.zero_pred8;
  memset(pred_buffer_mem.zero_pred8, 0,
         sizeof(*pred_buffer_mem.zero_pred8) * coeff_count);
#endif

  const BitDepthInfo bd_info = get_bit_depth_info(xd);
  cpi->norm_wiener_variance = 0;

  int mb_step = mi_size_wide[BLOCK_16X16];

  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      PREDICTION_MODE best_mode = DC_PRED;
      int best_intra_cost = INT_MAX;

      int mi_row = mb_row * mb_step;
      int mi_col = mb_col * mb_step;

      xd->up_available = mi_row > 0;
      xd->left_available = mi_col > 0;

      int dst_mb_offset = mi_row * MI_SIZE * buf_stride + mi_col * MI_SIZE;
      uint8_t *dst_buffer = xd->cur_buf->y_buffer + dst_mb_offset;

      for (PREDICTION_MODE mode = INTRA_MODE_START; mode < INTRA_MODE_END;
           ++mode) {
        av1_predict_intra_block(xd, cm->seq_params->sb_size,
                                cm->seq_params->enable_intra_edge_filter,
                                block_size, block_size, tx_size, mode, 0, 0,
                                FILTER_INTRA_MODES, dst_buffer, buf_stride,
                                pred_buf, block_size, 0, 0, 0);

        av1_subtract_block(bd_info, block_size, block_size, src_diff,
                           block_size, dst_buffer, buf_stride, pred_buf,
                           block_size);
        av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);
        int intra_cost = aom_satd(coeff, coeff_count);
        if (intra_cost < best_intra_cost) {
          best_intra_cost = intra_cost;
          best_mode = mode;
        }
      }

      int idx;
      int16_t median_val = 0;
      uint8_t *mb_buffer =
          buffer + mb_row * block_size * buf_stride + mb_col * block_size;
      int64_t wiener_variance = 0;
      av1_predict_intra_block(
          xd, cm->seq_params->sb_size, cm->seq_params->enable_intra_edge_filter,
          block_size, block_size, tx_size, best_mode, 0, 0, FILTER_INTRA_MODES,
          dst_buffer, buf_stride, pred_buf, block_size, 0, 0, 0);
      av1_subtract_block(bd_info, block_size, block_size, src_diff, block_size,
                         mb_buffer, buf_stride, pred_buf, block_size);
      av1_quick_txfm(0, tx_size, bd_info, src_diff, block_size, coeff);
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

  // Take saliency into consideration as weights for wiener variance
  int64_t max_wiener_var = 0;
  int64_t min_wiener_var = INT64_MAX;
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      max_wiener_var =
          AOMMAX(max_wiener_var,
                 cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
                     .mb_wiener_variance);
      min_wiener_var =
          AOMMIN(min_wiener_var,
                 cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
                     .mb_wiener_variance);
    }
  }
  /*
  printf("-----------------------\n");
  printf("max wiener var: %ld, min wiener var: %ld \n",
         max_wiener_var, min_wiener_var);
  */
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      double block_saliency = 0;
      int pixel_count = 0;
      for (int r = 0; r < block_size; ++r) {
        for (int c = 0; c < block_size; ++c) {
          const int row = mb_row * block_size + r;
          const int col = mb_col * block_size + c;
          if (row >= cpi->frame_info.frame_height ||
              col >= cpi->frame_info.frame_width) {
            continue;
          }
          block_saliency += saliency_map[row * buf_stride + col];
          ++pixel_count;
        }
      }
      block_saliency /= 256;
      block_saliency /= pixel_count;
      /*
      printf("saliency %.3f, wiener variance %ld ", block_saliency,
             cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
                 .mb_wiener_variance);
      */
      const double weighted_wiener_var =
          pow(cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
                  .mb_wiener_variance,
              0.1) /
          (block_saliency * block_saliency);
      //(double)max_wiener_var;
      cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
          .mb_wiener_variance = (int64_t)(weighted_wiener_var);
      /*
      printf("after %ld \n",
             cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols + mb_col]
                 .mb_wiener_variance);
      */
    }
  }

  /*
  FILE *pfile = fopen("mask.csv", "w");
  fprintf(pfile, "%d,%d,%ld\n", cpi->frame_info.frame_height,
          cpi->frame_info.frame_width, max_wiener_var);
  for (mb_row = 0; mb_row < cpi->frame_info.mb_rows; ++mb_row) {
    for (mb_col = 0; mb_col < cpi->frame_info.mb_cols; ++mb_col) {
      fprintf(pfile, "%ld,",
              cpi->mb_weber_stats[mb_row * cpi->frame_info.mb_cols +
  mb_col].mb_wiener_variance);
    }
    fprintf(pfile, "\n");
  }
  fclose(pfile);
  */

  int sb_step = mi_size_wide[cm->seq_params->sb_size];
  double sb_wiener_log = 0;
  int sb_count = 0;

  for (int mi_row = 0; mi_row < cm->mi_params.mi_rows; mi_row += sb_step) {
    for (int mi_col = 0; mi_col < cm->mi_params.mi_cols; mi_col += sb_step) {
      int sb_wiener_var =
          get_var_perceptual_ai(cpi, cm->seq_params->sb_size, mi_row, mi_col);
      sb_wiener_log += log(sb_wiener_var);
      ++sb_count;
    }
  }

  if (sb_count)
    cpi->norm_wiener_variance = (int64_t)(exp(sb_wiener_log / sb_count));
  // printf("norm_wiener_variance %ld\n", cpi->norm_wiener_variance);
  cpi->norm_wiener_variance = AOMMAX(1, cpi->norm_wiener_variance);
}

int av1_get_sbq_perceptual_ai(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                              int mi_col) {
  AV1_COMMON *const cm = &cpi->common;
  const int base_qindex = cm->quant_params.base_qindex;
  int sb_wiener_var = get_var_perceptual_ai(cpi, bsize, mi_row, mi_col);
  int offset = 0;
  double beta = (double)cpi->norm_wiener_variance / sb_wiener_var;
  /*
  // Cap beta such that the delta q value is not much far away from the base q.
  beta = AOMMIN(beta, 4);
  beta = AOMMAX(beta, 0.25);
  */
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
