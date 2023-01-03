/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "config/aom_dsp_rtcd.h"

#include "av1/common/blockd.h"
#include "av1/common/reconinter.h"

#include "av1/encoder/encodemv.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/model_rd.h"
#include "av1/encoder/nonrd_opt.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/nonrd_utils.h"

static const SCAN_ORDER av1_fast_idtx_scan_order_16x16 = {
  av1_fast_idtx_scan_16x16, av1_fast_idtx_iscan_16x16
};

static INLINE void set_force_skip_flag(const AV1_COMP *const cpi,
                                       MACROBLOCK *const x, unsigned int sse,
                                       int *force_skip) {
  if (x->txfm_search_params.tx_mode_search_type == TX_MODE_SELECT &&
      cpi->sf.rt_sf.tx_size_level_based_on_qstep &&
      cpi->sf.rt_sf.tx_size_level_based_on_qstep >= 2) {
    const int qstep = x->plane[AOM_PLANE_Y].dequant_QTX[1] >> (x->e_mbd.bd - 5);
    const unsigned int qstep_sq = qstep * qstep;
    // If the sse is low for low source variance blocks, mark those as
    // transform skip.
    // Note: Though qstep_sq is based on ac qstep, the threshold is kept
    // low so that reliable early estimate of tx skip can be obtained
    // through its comparison with sse.
    if (sse < qstep_sq && x->source_variance < qstep_sq &&
        x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] == 0 &&
        x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] == 0)
      *force_skip = 1;
  }
}

#define CAP_TX_SIZE_FOR_BSIZE_GT32(tx_mode_search_type, bsize) \
  (((tx_mode_search_type) != ONLY_4X4 && (bsize) > BLOCK_32X32) ? true : false)
#define TX_SIZE_FOR_BSIZE_GT32 (TX_16X16)

static TX_SIZE calculate_tx_size(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                 MACROBLOCK *const x, unsigned int var,
                                 unsigned int sse, int *force_skip) {
  MACROBLOCKD *const xd = &x->e_mbd;
  TX_SIZE tx_size;
  const TxfmSearchParams *txfm_params = &x->txfm_search_params;
  if (txfm_params->tx_mode_search_type == TX_MODE_SELECT) {
    int multiplier = 8;
    unsigned int var_thresh = 0;
    unsigned int is_high_var = 1;
    // Use quantizer based thresholds to determine transform size.
    if (cpi->sf.rt_sf.tx_size_level_based_on_qstep) {
      const int qband = x->qindex >> (QINDEX_BITS - 2);
      const int mult[4] = { 8, 7, 6, 5 };
      assert(qband < 4);
      multiplier = mult[qband];
      const int qstep = x->plane[AOM_PLANE_Y].dequant_QTX[1] >> (xd->bd - 5);
      const unsigned int qstep_sq = qstep * qstep;
      var_thresh = qstep_sq * 2;
      if (cpi->sf.rt_sf.tx_size_level_based_on_qstep >= 2) {
        // If the sse is low for low source variance blocks, mark those as
        // transform skip.
        // Note: Though qstep_sq is based on ac qstep, the threshold is kept
        // low so that reliable early estimate of tx skip can be obtained
        // through its comparison with sse.
        if (sse < qstep_sq && x->source_variance < qstep_sq &&
            x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] == 0 &&
            x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] == 0)
          *force_skip = 1;
        // Further lower transform size based on aq mode only if residual
        // variance is high.
        is_high_var = (var >= var_thresh);
      }
    }
    // Choose larger transform size for blocks where dc component is dominant or
    // the ac component is low.
    if (sse > ((var * multiplier) >> 2) || (var < var_thresh))
      tx_size =
          AOMMIN(max_txsize_lookup[bsize],
                 tx_mode_to_biggest_tx_size[txfm_params->tx_mode_search_type]);
    else
      tx_size = TX_8X8;

    if (cpi->oxcf.q_cfg.aq_mode == CYCLIC_REFRESH_AQ &&
        cyclic_refresh_segment_id_boosted(xd->mi[0]->segment_id) && is_high_var)
      tx_size = TX_8X8;
    else if (tx_size > TX_16X16)
      tx_size = TX_16X16;
  } else {
    tx_size =
        AOMMIN(max_txsize_lookup[bsize],
               tx_mode_to_biggest_tx_size[txfm_params->tx_mode_search_type]);
  }

  if (CAP_TX_SIZE_FOR_BSIZE_GT32(txfm_params->tx_mode_search_type, bsize))
    tx_size = TX_SIZE_FOR_BSIZE_GT32;

  return AOMMIN(tx_size, TX_16X16);
}

static void block_variance(const uint8_t *src, int src_stride,
                           const uint8_t *ref, int ref_stride, int w, int h,
                           unsigned int *sse, int *sum, int block_size,
                           uint32_t *sse8x8, int *sum8x8, uint32_t *var8x8) {
  int k = 0;
  *sse = 0;
  *sum = 0;

  // This function is called for block sizes >= BLOCK_32x32. As per the design
  // the aom_get_var_sse_sum_8x8_quad() processes four 8x8 blocks (in a 8x32)
  // per call. Hence the width and height of the block need to be at least 8 and
  // 32 samples respectively.
  assert(w >= 32);
  assert(h >= 8);
  for (int row = 0; row < h; row += block_size) {
    for (int col = 0; col < w; col += 32) {
      aom_get_var_sse_sum_8x8_quad(src + src_stride * row + col, src_stride,
                                   ref + ref_stride * row + col, ref_stride,
                                   &sse8x8[k], &sum8x8[k], sse, sum,
                                   &var8x8[k]);
      k += 4;
    }
  }
}

static void block_variance_16x16_dual(const uint8_t *src, int src_stride,
                                      const uint8_t *ref, int ref_stride, int w,
                                      int h, unsigned int *sse, int *sum,
                                      int block_size, uint32_t *sse16x16,
                                      uint32_t *var16x16) {
  int k = 0;
  *sse = 0;
  *sum = 0;
  // This function is called for block sizes >= BLOCK_32x32. As per the design
  // the aom_get_var_sse_sum_16x16_dual() processes four 16x16 blocks (in a
  // 16x32) per call. Hence the width and height of the block need to be at
  // least 16 and 32 samples respectively.
  assert(w >= 32);
  assert(h >= 16);
  for (int row = 0; row < h; row += block_size) {
    for (int col = 0; col < w; col += 32) {
      aom_get_var_sse_sum_16x16_dual(src + src_stride * row + col, src_stride,
                                     ref + ref_stride * row + col, ref_stride,
                                     &sse16x16[k], sse, sum, &var16x16[k]);
      k += 2;
    }
  }
}

static void calculate_variance(int bw, int bh, TX_SIZE tx_size,
                               unsigned int *sse_i, int *sum_i,
                               unsigned int *var_o, unsigned int *sse_o,
                               int *sum_o) {
  const BLOCK_SIZE unit_size = txsize_to_bsize[tx_size];
  const int nw = 1 << (bw - b_width_log2_lookup[unit_size]);
  const int nh = 1 << (bh - b_height_log2_lookup[unit_size]);
  int row, col, k = 0;

  for (row = 0; row < nh; row += 2) {
    for (col = 0; col < nw; col += 2) {
      sse_o[k] = sse_i[row * nw + col] + sse_i[row * nw + col + 1] +
                 sse_i[(row + 1) * nw + col] + sse_i[(row + 1) * nw + col + 1];
      sum_o[k] = sum_i[row * nw + col] + sum_i[row * nw + col + 1] +
                 sum_i[(row + 1) * nw + col] + sum_i[(row + 1) * nw + col + 1];
      var_o[k] = sse_o[k] - (uint32_t)(((int64_t)sum_o[k] * sum_o[k]) >>
                                       (b_width_log2_lookup[unit_size] +
                                        b_height_log2_lookup[unit_size] + 6));
      k++;
    }
  }
}

// Adjust the ac_thr according to speed, width, height and normalized sum
static int ac_thr_factor(int speed, int width, int height, int norm_sum) {
  if (speed >= 8 && norm_sum < 5) {
    if (width <= 640 && height <= 480)
      return 4;
    else
      return 2;
  }
  return 1;
}

// Sets early_term flag based on chroma planes prediction
static INLINE void set_early_term_based_on_uv_plane(
    AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize, MACROBLOCKD *xd, int mi_row,
    int mi_col, int *early_term, int num_blk, const unsigned int *sse_tx,
    const unsigned int *var_tx, int sum, unsigned int var, unsigned int sse) {
  AV1_COMMON *const cm = &cpi->common;
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  const uint32_t dc_quant = p->dequant_QTX[0];
  const uint32_t ac_quant = p->dequant_QTX[1];
  const int64_t dc_thr = dc_quant * dc_quant >> 6;
  int64_t ac_thr = ac_quant * ac_quant >> 6;
  const int bw = b_width_log2_lookup[bsize];
  const int bh = b_height_log2_lookup[bsize];
  int ac_test = 1;
  int dc_test = 1;
  const int norm_sum = abs(sum) >> (bw + bh);

#if CONFIG_AV1_TEMPORAL_DENOISING
  if (cpi->oxcf.noise_sensitivity > 0 && denoise_svc(cpi) &&
      cpi->oxcf.speed > 5)
    ac_thr = av1_scale_acskip_thresh(ac_thr, cpi->denoiser.denoising_level,
                                     norm_sum, cpi->svc.temporal_layer_id);
  else
    ac_thr *= ac_thr_factor(cpi->oxcf.speed, cm->width, cm->height, norm_sum);
#else
  ac_thr *= ac_thr_factor(cpi->oxcf.speed, cm->width, cm->height, norm_sum);

#endif

  for (int k = 0; k < num_blk; k++) {
    // Check if all ac coefficients can be quantized to zero.
    if (!(var_tx[k] < ac_thr || var == 0)) {
      ac_test = 0;
      break;
    }
    // Check if dc coefficient can be quantized to zero.
    if (!(sse_tx[k] - var_tx[k] < dc_thr || sse == var)) {
      dc_test = 0;
      break;
    }
  }

  // Check if chroma can be skipped based on ac and dc test flags.
  if (ac_test && dc_test) {
    int skip_uv[2] = { 0 };
    unsigned int var_uv[2];
    unsigned int sse_uv[2];
    // Transform skipping test in UV planes.
    for (int plane = AOM_PLANE_U; plane <= AOM_PLANE_V; plane++) {
      int j = plane - 1;
      skip_uv[j] = 1;
      if (x->color_sensitivity[COLOR_SENS_IDX(plane)]) {
        skip_uv[j] = 0;
        struct macroblock_plane *const puv = &x->plane[plane];
        struct macroblockd_plane *const puvd = &xd->plane[plane];
        const BLOCK_SIZE uv_bsize = get_plane_block_size(
            bsize, puvd->subsampling_x, puvd->subsampling_y);
        // Adjust these thresholds for UV.
        const int64_t uv_dc_thr =
            (puv->dequant_QTX[0] * puv->dequant_QTX[0]) >> 3;
        const int64_t uv_ac_thr =
            (puv->dequant_QTX[1] * puv->dequant_QTX[1]) >> 3;
        av1_enc_build_inter_predictor(cm, xd, mi_row, mi_col, NULL, bsize,
                                      plane, plane);
        var_uv[j] = cpi->ppi->fn_ptr[uv_bsize].vf(puv->src.buf, puv->src.stride,
                                                  puvd->dst.buf,
                                                  puvd->dst.stride, &sse_uv[j]);
        if ((var_uv[j] < uv_ac_thr || var_uv[j] == 0) &&
            (sse_uv[j] - var_uv[j] < uv_dc_thr || sse_uv[j] == var_uv[j]))
          skip_uv[j] = 1;
        else
          break;
      }
    }
    if (skip_uv[0] & skip_uv[1]) {
      *early_term = 1;
    }
  }
}

static INLINE void calc_rate_dist_block_param(AV1_COMP *cpi, MACROBLOCK *x,
                                              RD_STATS *rd_stats,
                                              int calculate_rd, int *early_term,
                                              BLOCK_SIZE bsize,
                                              unsigned int sse) {
  if (calculate_rd) {
    if (!*early_term) {
      const int bw = block_size_wide[bsize];
      const int bh = block_size_high[bsize];

      model_rd_with_curvfit(cpi, x, bsize, AOM_PLANE_Y, rd_stats->sse, bw * bh,
                            &rd_stats->rate, &rd_stats->dist);
    }

    if (*early_term) {
      rd_stats->rate = 0;
      rd_stats->dist = sse << 4;
    }
  }
}

void av1_model_skip_for_sb_y_large_64(AV1_COMP *cpi, BLOCK_SIZE bsize,
                                      int mi_row, int mi_col, MACROBLOCK *x,
                                      MACROBLOCKD *xd, RD_STATS *rd_stats,
                                      int *early_term, int calculate_rd,
                                      int64_t best_sse,
                                      unsigned int *var_output,
                                      unsigned int var_prune_threshold) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  unsigned int sse;
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  int test_skip = 1;
  unsigned int var;
  int sum;
  const int bw = b_width_log2_lookup[bsize];
  const int bh = b_height_log2_lookup[bsize];
  unsigned int sse16x16[64] = { 0 };
  unsigned int var16x16[64] = { 0 };
  assert(xd->mi[0]->tx_size == TX_16X16);
  assert(bsize > BLOCK_32X32);

  // Calculate variance for whole partition, and also save 16x16 blocks'
  // variance to be used in following transform skipping test.
  block_variance_16x16_dual(p->src.buf, p->src.stride, pd->dst.buf,
                            pd->dst.stride, 4 << bw, 4 << bh, &sse, &sum, 16,
                            sse16x16, var16x16);

  var = sse - (unsigned int)(((int64_t)sum * sum) >> (bw + bh + 4));
  if (var_output) {
    *var_output = var;
    if (*var_output > var_prune_threshold) {
      return;
    }
  }

  rd_stats->sse = sse;
  // Skipping test
  *early_term = 0;
  set_force_skip_flag(cpi, x, sse, early_term);
  // The code below for setting skip flag assumes transform size of at least
  // 8x8, so force this lower limit on transform.
  MB_MODE_INFO *const mi = xd->mi[0];
  if (!calculate_rd && cpi->sf.rt_sf.sse_early_term_inter_search &&
      early_term_inter_search_with_sse(
          cpi->sf.rt_sf.sse_early_term_inter_search, bsize, sse, best_sse,
          mi->mode))
    test_skip = 0;

  if (*early_term) test_skip = 0;

  // Evaluate if the partition block is a skippable block in Y plane.
  if (test_skip) {
    const unsigned int *sse_tx = sse16x16;
    const unsigned int *var_tx = var16x16;
    const unsigned int num_block = (1 << (bw + bh - 2)) >> 2;
    set_early_term_based_on_uv_plane(cpi, x, bsize, xd, mi_row, mi_col,
                                     early_term, num_block, sse_tx, var_tx, sum,
                                     var, sse);
  }
  calc_rate_dist_block_param(cpi, x, rd_stats, calculate_rd, early_term, bsize,
                             sse);
}

void av1_model_skip_for_sb_y_large(AV1_COMP *cpi, BLOCK_SIZE bsize, int mi_row,
                                   int mi_col, MACROBLOCK *x, MACROBLOCKD *xd,
                                   RD_STATS *rd_stats, int *early_term,
                                   int calculate_rd, int64_t best_sse,
                                   unsigned int *var_output,
                                   unsigned int var_prune_threshold) {
  if (x->force_zeromv_skip_for_blk) {
    *early_term = 1;
    rd_stats->rate = 0;
    rd_stats->dist = 0;
    rd_stats->sse = 0;
    return;
  }

  // For block sizes greater than 32x32, the transform size is always 16x16.
  // This function avoids calling calculate_variance() for tx_size 16x16 cases
  // by directly populating variance at tx_size level from
  // block_variance_16x16_dual() function.
  const TxfmSearchParams *txfm_params = &x->txfm_search_params;
  if (CAP_TX_SIZE_FOR_BSIZE_GT32(txfm_params->tx_mode_search_type, bsize)) {
    xd->mi[0]->tx_size = TX_SIZE_FOR_BSIZE_GT32;
    av1_model_skip_for_sb_y_large_64(cpi, bsize, mi_row, mi_col, x, xd,
                                     rd_stats, early_term, calculate_rd,
                                     best_sse, var_output, var_prune_threshold);
    return;
  }

  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  unsigned int sse;
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  int test_skip = 1;
  unsigned int var;
  int sum;

  const int bw = b_width_log2_lookup[bsize];
  const int bh = b_height_log2_lookup[bsize];
  unsigned int sse8x8[256] = { 0 };
  int sum8x8[256] = { 0 };
  unsigned int var8x8[256] = { 0 };
  TX_SIZE tx_size;

  // Calculate variance for whole partition, and also save 8x8 blocks' variance
  // to be used in following transform skipping test.
  block_variance(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
                 4 << bw, 4 << bh, &sse, &sum, 8, sse8x8, sum8x8, var8x8);
  var = sse - (unsigned int)(((int64_t)sum * sum) >> (bw + bh + 4));
  if (var_output) {
    *var_output = var;
    if (*var_output > var_prune_threshold) {
      return;
    }
  }

  rd_stats->sse = sse;
  // Skipping test
  *early_term = 0;
  tx_size = calculate_tx_size(cpi, bsize, x, var, sse, early_term);
  assert(tx_size <= TX_16X16);
  // The code below for setting skip flag assumes transform size of at least
  // 8x8, so force this lower limit on transform.
  if (tx_size < TX_8X8) tx_size = TX_8X8;
  xd->mi[0]->tx_size = tx_size;

  MB_MODE_INFO *const mi = xd->mi[0];
  if (!calculate_rd && cpi->sf.rt_sf.sse_early_term_inter_search &&
      early_term_inter_search_with_sse(
          cpi->sf.rt_sf.sse_early_term_inter_search, bsize, sse, best_sse,
          mi->mode))
    test_skip = 0;

  if (*early_term) test_skip = 0;

  // Evaluate if the partition block is a skippable block in Y plane.
  if (test_skip) {
    unsigned int sse16x16[64] = { 0 };
    int sum16x16[64] = { 0 };
    unsigned int var16x16[64] = { 0 };
    const unsigned int *sse_tx = sse8x8;
    const unsigned int *var_tx = var8x8;
    unsigned int num_blks = 1 << (bw + bh - 2);

    if (tx_size >= TX_16X16) {
      calculate_variance(bw, bh, TX_8X8, sse8x8, sum8x8, var16x16, sse16x16,
                         sum16x16);
      sse_tx = sse16x16;
      var_tx = var16x16;
      num_blks = num_blks >> 2;
    }
    set_early_term_based_on_uv_plane(cpi, x, bsize, xd, mi_row, mi_col,
                                     early_term, num_blks, sse_tx, var_tx, sum,
                                     var, sse);
  }
  calc_rate_dist_block_param(cpi, x, rd_stats, calculate_rd, early_term, bsize,
                             sse);
}

void av1_model_rd_for_sb_y(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                           MACROBLOCK *x, MACROBLOCKD *xd, RD_STATS *rd_stats,
                           unsigned int *var_out, int calculate_rd,
                           int *early_term) {
  if (x->force_zeromv_skip_for_blk && early_term != NULL) {
    *early_term = 1;
    rd_stats->rate = 0;
    rd_stats->dist = 0;
    rd_stats->sse = 0;
  }

  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  const int ref = xd->mi[0]->ref_frame[0];

  assert(bsize < BLOCK_SIZES_ALL);

  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  unsigned int sse;
  int rate;
  int64_t dist;

  unsigned int var = cpi->ppi->fn_ptr[bsize].vf(
      p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, &sse);
  int force_skip = 0;
  xd->mi[0]->tx_size = calculate_tx_size(cpi, bsize, x, var, sse, &force_skip);
  if (var_out) {
    *var_out = var;
  }

  if (calculate_rd && (!force_skip || ref == INTRA_FRAME)) {
    const int bwide = block_size_wide[bsize];
    const int bhigh = block_size_high[bsize];
    model_rd_with_curvfit(cpi, x, bsize, AOM_PLANE_Y, sse, bwide * bhigh, &rate,
                          &dist);
  } else {
    rate = INT_MAX;  // this will be overwritten later with av1_block_yrd
    dist = INT_MAX;
  }
  rd_stats->sse = sse;
  x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

  if (force_skip && ref > INTRA_FRAME) {
    rate = 0;
    dist = (int64_t)sse << 4;
  }

  assert(rate >= 0);

  rd_stats->skip_txfm = (rate == 0);
  rate = AOMMIN(rate, INT_MAX);
  rd_stats->rate = rate;
  rd_stats->dist = dist;
}

static INLINE void aom_process_hadamard_lp_8x16(MACROBLOCK *x,
                                                int max_blocks_high,
                                                int max_blocks_wide,
                                                int num_4x4_w, int step,
                                                int block_step) {
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  const int bw = 4 * num_4x4_w;
  const int num_4x4 = AOMMIN(num_4x4_w, max_blocks_wide);
  int block = 0;

  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0; c < num_4x4; c += 2 * block_step) {
      const int16_t *src_diff = &p->src_diff[(r * bw + c) << 2];
      int16_t *low_coeff = (int16_t *)p->coeff + BLOCK_OFFSET(block);
      aom_hadamard_lp_8x8_dual(src_diff, (ptrdiff_t)bw, low_coeff);
      block += 2 * step;
    }
  }
}

#define DECLARE_BLOCK_YRD_BUFFERS()                      \
  DECLARE_ALIGNED(64, tran_low_t, dqcoeff_buf[16 * 16]); \
  DECLARE_ALIGNED(64, tran_low_t, qcoeff_buf[16 * 16]);  \
  DECLARE_ALIGNED(64, tran_low_t, coeff_buf[16 * 16]);   \
  uint16_t eob[1];

#define DECLARE_BLOCK_YRD_VARS()                                          \
  /* When is_tx_8x8_dual_applicable is true, we compute the txfm for the  \
   * entire bsize and write macroblock_plane::coeff. So low_coeff is kept \
   * as a non-const so we can reassign it to macroblock_plane::coeff. */  \
  int16_t *low_coeff = (int16_t *)coeff_buf;                              \
  int16_t *const low_qcoeff = (int16_t *)qcoeff_buf;                      \
  int16_t *const low_dqcoeff = (int16_t *)dqcoeff_buf;                    \
  const int diff_stride = bw;

#define DECLARE_LOOP_VARS_BLOCK_YRD() \
  const int16_t *src_diff = &p->src_diff[(r * diff_stride + c) << 2];

#if CONFIG_AV1_HIGHBITDEPTH
#define DECLARE_BLOCK_YRD_HBD_VARS()     \
  tran_low_t *const coeff = coeff_buf;   \
  tran_low_t *const qcoeff = qcoeff_buf; \
  tran_low_t *const dqcoeff = dqcoeff_buf;

static AOM_FORCE_INLINE void update_yrd_loop_vars_hbd(
    MACROBLOCK *x, int *skippable, int step, int ncoeffs,
    tran_low_t *const coeff, tran_low_t *const qcoeff,
    tran_low_t *const dqcoeff, RD_STATS *this_rdc, int *eob_cost,
    int tx_blk_id) {
  const int is_txfm_skip = (ncoeffs == 0);
  *skippable &= is_txfm_skip;
  x->txfm_search_info.blk_skip[tx_blk_id] = is_txfm_skip;
  *eob_cost += get_msb(ncoeffs + 1);

  int64_t dummy;
  if (ncoeffs == 1)
    this_rdc->rate += (int)abs(qcoeff[0]);
  else if (ncoeffs > 1)
    this_rdc->rate += aom_satd(qcoeff, step << 4);

  this_rdc->dist += av1_block_error(coeff, dqcoeff, step << 4, &dummy) >> 2;
}
#endif

static AOM_FORCE_INLINE void update_yrd_loop_vars(
    MACROBLOCK *x, int *skippable, int step, int ncoeffs,
    int16_t *const low_coeff, int16_t *const low_qcoeff,
    int16_t *const low_dqcoeff, RD_STATS *this_rdc, int *eob_cost,
    int tx_blk_id) {
  const int is_txfm_skip = (ncoeffs == 0);
  *skippable &= is_txfm_skip;
  x->txfm_search_info.blk_skip[tx_blk_id] = is_txfm_skip;
  *eob_cost += get_msb(ncoeffs + 1);
  if (ncoeffs == 1)
    this_rdc->rate += (int)abs(low_qcoeff[0]);
  else if (ncoeffs > 1)
    this_rdc->rate += aom_satd_lp(low_qcoeff, step << 4);

  this_rdc->dist += av1_block_error_lp(low_coeff, low_dqcoeff, step << 4) >> 2;
}

/*!\brief Calculates RD Cost using Hadamard transform.
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Calculates RD Cost using Hadamard transform. For low bit depth this function
 * uses low-precision set of functions (16-bit) and 32 bit for high bit depth
 * \param[in]    x              Pointer to structure holding all the data for
                                the current macroblock
 * \param[in]    this_rdc       Pointer to calculated RD Cost
 * \param[in]    skippable      Pointer to a flag indicating possible tx skip
 * \param[in]    bsize          Current block size
 * \param[in]    tx_size        Transform size
 * \param[in]    is_inter_mode  Flag to indicate inter mode
 *
 * \remark Nothing is returned. Instead, calculated RD cost is placed to
 * \c this_rdc. \c skippable flag is set if there is no non-zero quantized
 * coefficients for Hadamard transform
 */
void av1_block_yrd(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                   BLOCK_SIZE bsize, TX_SIZE tx_size, int is_inter_mode) {
  MACROBLOCKD *xd = &x->e_mbd;
  const struct macroblockd_plane *pd = &xd->plane[AOM_PLANE_Y];
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  assert(bsize < BLOCK_SIZES_ALL);
  const int num_4x4_w = mi_size_wide[bsize];
  const int num_4x4_h = mi_size_high[bsize];
  const int step = 1 << (tx_size << 1);
  const int block_step = (1 << tx_size);
  const int row_step = step * num_4x4_w >> tx_size;
  int block = 0;
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >> 5);
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 : xd->mb_to_bottom_edge >> 5);
  int eob_cost = 0;
  const int bw = 4 * num_4x4_w;
  const int bh = 4 * num_4x4_h;
  const int use_hbd = is_cur_buf_hbd(xd);
  int num_blk_skip_w = num_4x4_w;
  int sh_blk_skip = 0;
  if (is_inter_mode) {
    num_blk_skip_w = num_4x4_w >> 1;
    sh_blk_skip = 1;
  }

#if CONFIG_AV1_HIGHBITDEPTH
  if (use_hbd) {
    aom_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride);
  } else {
    aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                       pd->dst.buf, pd->dst.stride);
  }
#else
  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
#endif

  // Keep the intermediate value on the stack here. Writing directly to
  // skippable causes speed regression due to load-and-store issues in
  // update_yrd_loop_vars.
  int temp_skippable = 1;
  this_rdc->dist = 0;
  this_rdc->rate = 0;
  // For block sizes 8x16 or above, Hadamard txfm of two adjacent 8x8 blocks
  // can be done per function call. Hence the call of Hadamard txfm is
  // abstracted here for the specified cases.
  int is_tx_8x8_dual_applicable =
      (tx_size == TX_8X8 && block_size_wide[bsize] >= 16 &&
       block_size_high[bsize] >= 8);

#if CONFIG_AV1_HIGHBITDEPTH
  // As of now, dual implementation of hadamard txfm is available for low
  // bitdepth.
  if (use_hbd) is_tx_8x8_dual_applicable = 0;
#endif

  if (is_tx_8x8_dual_applicable) {
    aom_process_hadamard_lp_8x16(x, max_blocks_high, max_blocks_wide, num_4x4_w,
                                 step, block_step);
  }

  const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];
  DECLARE_BLOCK_YRD_BUFFERS()
  DECLARE_BLOCK_YRD_VARS()
#if CONFIG_AV1_HIGHBITDEPTH
  DECLARE_BLOCK_YRD_HBD_VARS()
#else
  (void)use_hbd;
#endif

  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0, s = 0; c < max_blocks_wide; c += block_step, s += step) {
      DECLARE_LOOP_VARS_BLOCK_YRD()

      switch (tx_size) {
#if CONFIG_AV1_HIGHBITDEPTH
        case TX_16X16:
          if (use_hbd) {
            aom_hadamard_16x16(src_diff, diff_stride, coeff);
            av1_quantize_fp(coeff, 16 * 16, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob,
                            // default_scan_fp_16x16_transpose and
                            // av1_default_iscan_fp_16x16_transpose have to be
                            // used together.
                            default_scan_fp_16x16_transpose,
                            av1_default_iscan_fp_16x16_transpose);
          } else {
            aom_hadamard_lp_16x16(src_diff, diff_stride, low_coeff);
            av1_quantize_lp(low_coeff, 16 * 16, p->round_fp_QTX,
                            p->quant_fp_QTX, low_qcoeff, low_dqcoeff,
                            p->dequant_QTX, eob,
                            // default_scan_lp_16x16_transpose and
                            // av1_default_iscan_lp_16x16_transpose have to be
                            // used together.
                            default_scan_lp_16x16_transpose,
                            av1_default_iscan_lp_16x16_transpose);
          }
          break;
        case TX_8X8:
          if (use_hbd) {
            aom_hadamard_8x8(src_diff, diff_stride, coeff);
            av1_quantize_fp(
                coeff, 8 * 8, p->zbin_QTX, p->round_fp_QTX, p->quant_fp_QTX,
                p->quant_shift_QTX, qcoeff, dqcoeff, p->dequant_QTX, eob,
                default_scan_8x8_transpose, av1_default_iscan_8x8_transpose);
          } else {
            if (is_tx_8x8_dual_applicable) {
              // The coeffs are pre-computed for the whole block, so re-assign
              // low_coeff to the appropriate location.
              const int block_offset = BLOCK_OFFSET(block + s);
              low_coeff = (int16_t *)p->coeff + block_offset;
            } else {
              aom_hadamard_lp_8x8(src_diff, diff_stride, low_coeff);
            }
            av1_quantize_lp(
                low_coeff, 8 * 8, p->round_fp_QTX, p->quant_fp_QTX, low_qcoeff,
                low_dqcoeff, p->dequant_QTX, eob,
                // default_scan_8x8_transpose and
                // av1_default_iscan_8x8_transpose have to be used together.
                default_scan_8x8_transpose, av1_default_iscan_8x8_transpose);
          }
          break;
        default:
          assert(tx_size == TX_4X4);
          // In tx_size=4x4 case, aom_fdct4x4 and aom_fdct4x4_lp generate
          // normal coefficients order, so we don't need to change the scan
          // order here.
          if (use_hbd) {
            aom_fdct4x4(src_diff, coeff, diff_stride);
            av1_quantize_fp(coeff, 4 * 4, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob, scan_order->scan,
                            scan_order->iscan);
          } else {
            aom_fdct4x4_lp(src_diff, low_coeff, diff_stride);
            av1_quantize_lp(low_coeff, 4 * 4, p->round_fp_QTX, p->quant_fp_QTX,
                            low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                            scan_order->scan, scan_order->iscan);
          }
          break;
#else
        case TX_16X16:
          aom_hadamard_lp_16x16(src_diff, diff_stride, low_coeff);
          av1_quantize_lp(low_coeff, 16 * 16, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          default_scan_lp_16x16_transpose,
                          av1_default_iscan_lp_16x16_transpose);
          break;
        case TX_8X8:
          if (is_tx_8x8_dual_applicable) {
            // The coeffs are pre-computed for the whole block, so re-assign
            // low_coeff to the appropriate location.
            const int block_offset = BLOCK_OFFSET(block + s);
            low_coeff = (int16_t *)p->coeff + block_offset;
          } else {
            aom_hadamard_lp_8x8(src_diff, diff_stride, low_coeff);
          }
          av1_quantize_lp(low_coeff, 8 * 8, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          default_scan_8x8_transpose,
                          av1_default_iscan_8x8_transpose);
          break;
        default:
          aom_fdct4x4_lp(src_diff, low_coeff, diff_stride);
          av1_quantize_lp(low_coeff, 4 * 4, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          scan_order->scan, scan_order->iscan);
          break;
#endif
      }
      assert(*eob <= 1024);
#if CONFIG_AV1_HIGHBITDEPTH
      if (use_hbd)
        update_yrd_loop_vars_hbd(x, &temp_skippable, step, *eob, coeff, qcoeff,
                                 dqcoeff, this_rdc, &eob_cost,
                                 (r * num_blk_skip_w + c) >> sh_blk_skip);
      else
#endif
        update_yrd_loop_vars(x, &temp_skippable, step, *eob, low_coeff,
                             low_qcoeff, low_dqcoeff, this_rdc, &eob_cost,
                             (r * num_blk_skip_w + c) >> sh_blk_skip);
    }
    block += row_step;
  }

  this_rdc->skip_txfm = *skippable = temp_skippable;
  if (this_rdc->sse < INT64_MAX) {
    this_rdc->sse = (this_rdc->sse << 6) >> 2;
    if (temp_skippable) {
      this_rdc->dist = 0;
      this_rdc->dist = this_rdc->sse;
      return;
    }
  }

  // If skippable is set, rate gets clobbered later.
  this_rdc->rate <<= (2 + AV1_PROB_COST_SHIFT);
  this_rdc->rate += (eob_cost << AV1_PROB_COST_SHIFT);
}

// Explicitly enumerate the cases so the compiler can generate SIMD for the
// function. According to the disassembler, gcc generates SSE codes for each of
// the possible block sizes. The hottest case is tx_width 16, which takes up
// about 8% of the self cycle of av1_nonrd_pick_inter_mode_sb. Since
// av1_nonrd_pick_inter_mode_sb takes up about 3% of total encoding time, the
// potential room of improvement for writing AVX2 optimization is only 3% * 8% =
// 0.24% of total encoding time.
static AOM_INLINE void scale_square_buf_vals(int16_t *dst, int tx_width,
                                             const int16_t *src,
                                             int src_stride) {
#define DO_SCALING                                                   \
  do {                                                               \
    for (int idy = 0; idy < tx_width; ++idy) {                       \
      for (int idx = 0; idx < tx_width; ++idx) {                     \
        dst[idy * tx_width + idx] = src[idy * src_stride + idx] * 8; \
      }                                                              \
    }                                                                \
  } while (0)

  if (tx_width == 4) {
    DO_SCALING;
  } else if (tx_width == 8) {
    DO_SCALING;
  } else if (tx_width == 16) {
    DO_SCALING;
  } else {
    assert(0);
  }

#undef DO_SCALING
}

/*!\brief Calculates RD Cost when the block uses Identity transform.
 * Note that thie function is only for low bit depth encoding, since it
 * is called in real-time mode for now, which sets high bit depth to 0:
 * -DCONFIG_AV1_HIGHBITDEPTH=0
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Calculates RD Cost. For low bit depth this function
 * uses low-precision set of functions (16-bit) and 32 bit for high bit depth
 * \param[in]    x              Pointer to structure holding all the data for
                                the current macroblock
 * \param[in]    this_rdc       Pointer to calculated RD Cost
 * \param[in]    skippable      Pointer to a flag indicating possible tx skip
 * \param[in]    bsize          Current block size
 * \param[in]    tx_size        Transform size
 *
 * \remark Nothing is returned. Instead, calculated RD cost is placed to
 * \c this_rdc. \c skippable flag is set if all coefficients are zero.
 */
void av1_block_yrd_idtx(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                        BLOCK_SIZE bsize, TX_SIZE tx_size) {
  MACROBLOCKD *xd = &x->e_mbd;
  const struct macroblockd_plane *pd = &xd->plane[AOM_PLANE_Y];
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  assert(bsize < BLOCK_SIZES_ALL);
  const int num_4x4_w = mi_size_wide[bsize];
  const int num_4x4_h = mi_size_high[bsize];
  const int step = 1 << (tx_size << 1);
  const int block_step = (1 << tx_size);
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >> 5);
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 : xd->mb_to_bottom_edge >> 5);
  int eob_cost = 0;
  const int bw = 4 * num_4x4_w;
  const int bh = 4 * num_4x4_h;
  const int num_blk_skip_w = num_4x4_w >> 1;
  const int sh_blk_skip = 1;
  // Keep the intermediate value on the stack here. Writing directly to
  // skippable causes speed regression due to load-and-store issues in
  // update_yrd_loop_vars.
  int temp_skippable = 1;
  int tx_wd = 0;
  const SCAN_ORDER *scan_order = NULL;
  switch (tx_size) {
    case TX_64X64:
      assert(0);  // Not implemented
      break;
    case TX_32X32:
      assert(0);  // Not used
      break;
    case TX_16X16:
      scan_order = &av1_fast_idtx_scan_order_16x16;
      tx_wd = 16;
      break;
    case TX_8X8:
      scan_order = &av1_fast_idtx_scan_order_8x8;
      tx_wd = 8;
      break;
    default:
      assert(tx_size == TX_4X4);
      scan_order = &av1_fast_idtx_scan_order_4x4;
      tx_wd = 4;
      break;
  }
  assert(scan_order != NULL);

  this_rdc->dist = 0;
  this_rdc->rate = 0;
  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  DECLARE_BLOCK_YRD_BUFFERS()
  DECLARE_BLOCK_YRD_VARS()
  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0, s = 0; c < max_blocks_wide; c += block_step, s += step) {
      DECLARE_LOOP_VARS_BLOCK_YRD()
      scale_square_buf_vals(low_coeff, tx_wd, src_diff, diff_stride);
      av1_quantize_lp(low_coeff, tx_wd * tx_wd, p->round_fp_QTX,
                      p->quant_fp_QTX, low_qcoeff, low_dqcoeff, p->dequant_QTX,
                      eob, scan_order->scan, scan_order->iscan);
      assert(*eob <= 1024);
      update_yrd_loop_vars(x, &temp_skippable, step, *eob, low_coeff,
                           low_qcoeff, low_dqcoeff, this_rdc, &eob_cost,
                           (r * num_blk_skip_w + c) >> sh_blk_skip);
    }
  }
  this_rdc->skip_txfm = *skippable = temp_skippable;
  if (this_rdc->sse < INT64_MAX) {
    this_rdc->sse = (this_rdc->sse << 6) >> 2;
    if (temp_skippable) {
      this_rdc->dist = 0;
      this_rdc->dist = this_rdc->sse;
      return;
    }
  }
  // If skippable is set, rate gets clobbered later.
  this_rdc->rate <<= (2 + AV1_PROB_COST_SHIFT);
  this_rdc->rate += (eob_cost << AV1_PROB_COST_SHIFT);
}

#if CONFIG_INTERNAL_STATS
void av1_store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                              int mode_index) {
#else
void av1_store_coding_context(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
#endif  // CONFIG_INTERNAL_STATS
  MACROBLOCKD *const xd = &x->e_mbd;
  TxfmSearchInfo *txfm_info = &x->txfm_search_info;

  // Take a snapshot of the coding context so it can be
  // restored if we decide to encode this way
  ctx->rd_stats.skip_txfm = txfm_info->skip_txfm;

  ctx->skippable = txfm_info->skip_txfm;
#if CONFIG_INTERNAL_STATS
  ctx->best_mode_index = mode_index;
#endif  // CONFIG_INTERNAL_STATS
  ctx->mic = *xd->mi[0];
  ctx->skippable = txfm_info->skip_txfm;
  av1_copy_mbmi_ext_to_mbmi_ext_frame(&ctx->mbmi_ext_best, &x->mbmi_ext,
                                      av1_ref_frame_type(xd->mi[0]->ref_frame));
}

int64_t av1_model_rd_for_sb_uv(AV1_COMP *cpi, BLOCK_SIZE plane_bsize,
                               MACROBLOCK *x, MACROBLOCKD *xd,
                               RD_STATS *this_rdc, int start_plane,
                               int stop_plane) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  unsigned int sse;
  int rate;
  int64_t dist;
  int plane;
  int64_t tot_sse = 0;

  this_rdc->rate = 0;
  this_rdc->dist = 0;
  this_rdc->skip_txfm = 0;

  for (plane = start_plane; plane <= stop_plane; ++plane) {
    struct macroblock_plane *const p = &x->plane[plane];
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const uint32_t dc_quant = p->dequant_QTX[0];
    const uint32_t ac_quant = p->dequant_QTX[1];
    const BLOCK_SIZE bs = plane_bsize;
    unsigned int var;
    if (!x->color_sensitivity[COLOR_SENS_IDX(plane)]) continue;

    var = cpi->ppi->fn_ptr[bs].vf(p->src.buf, p->src.stride, pd->dst.buf,
                                  pd->dst.stride, &sse);
    assert(sse >= var);
    tot_sse += sse;

    av1_model_rd_from_var_lapndz(sse - var, num_pels_log2_lookup[bs],
                                 dc_quant >> 3, &rate, &dist);

    this_rdc->rate += rate >> 1;
    this_rdc->dist += dist << 3;

    av1_model_rd_from_var_lapndz(var, num_pels_log2_lookup[bs], ac_quant >> 3,
                                 &rate, &dist);

    this_rdc->rate += rate;
    this_rdc->dist += dist << 4;
  }

  if (this_rdc->rate == 0) {
    this_rdc->skip_txfm = 1;
  }

  if (RDCOST(x->rdmult, this_rdc->rate, this_rdc->dist) >=
      RDCOST(x->rdmult, 0, tot_sse << 4)) {
    this_rdc->rate = 0;
    this_rdc->dist = tot_sse << 4;
    this_rdc->skip_txfm = 1;
  }

  return tot_sse;
}

void av1_set_color_sensitivity(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                               int y_sad, unsigned int source_variance,
                               struct buf_2d yv12_mb[MAX_MB_PLANE]) {
  const int subsampling_x = cpi->common.seq_params->subsampling_x;
  const int subsampling_y = cpi->common.seq_params->subsampling_y;
  int factor = (bsize >= BLOCK_32X32) ? 2 : 3;
  int shift = 3;
  if (cpi->oxcf.tune_cfg.content == AOM_CONTENT_SCREEN &&
      cpi->rc.high_source_sad) {
    factor = 1;
    shift = 6;
  }
  NOISE_LEVEL noise_level = kLow;
  int norm_sad =
      y_sad >> (b_width_log2_lookup[bsize] + b_height_log2_lookup[bsize]);
  unsigned int thresh_spatial = (cpi->common.width > 1920) ? 5000 : 1000;
  // If the spatial source variance is high and the normalized y_sad
  // is low, then y-channel is likely good for mode estimation, so keep
  // color_sensitivity off. For low noise content for now, since there is
  // some bdrate regression for noisy color clip.
  if (cpi->noise_estimate.enabled)
    noise_level = av1_noise_estimate_extract_level(&cpi->noise_estimate);
  if (noise_level == kLow && source_variance > thresh_spatial &&
      cpi->oxcf.tune_cfg.content != AOM_CONTENT_SCREEN && norm_sad < 50) {
    x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] = 0;
    x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] = 0;
    return;
  }
  const int num_planes = av1_num_planes(&cpi->common);

  for (int plane = AOM_PLANE_U; plane < num_planes; ++plane) {
    if (x->color_sensitivity[COLOR_SENS_IDX(plane)] == 2 ||
        source_variance < 50) {
      struct macroblock_plane *const p = &x->plane[plane];
      const BLOCK_SIZE bs =
          get_plane_block_size(bsize, subsampling_x, subsampling_y);

      const int uv_sad = cpi->ppi->fn_ptr[bs].sdf(
          p->src.buf, p->src.stride, yv12_mb[plane].buf, yv12_mb[plane].stride);

      const int norm_uv_sad =
          uv_sad >> (b_width_log2_lookup[bs] + b_height_log2_lookup[bs]);
      x->color_sensitivity[COLOR_SENS_IDX(plane)] =
          uv_sad > (factor * (y_sad >> shift)) && norm_uv_sad > 40;
      if (source_variance < 50 && norm_uv_sad > 100)
        x->color_sensitivity[COLOR_SENS_IDX(plane)] = 1;
    }
  }
}
