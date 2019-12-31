/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "aom_ports/system_state.h"
#include "av1/encoder/model_rd.h"
#include "av1/encoder/rdopt_utils.h"
#include "config/aom_dsp_rtcd.h"
#include "av1/encoder/pustats.h"

static double get_sse_norm(const int16_t *diff, int stride, int w, int h) {
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

static double get_sad_norm(const int16_t *diff, int stride, int w, int h) {
  double sum = 0.0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      sum += abs(diff[j * stride + i]);
    }
  }
  assert(w > 0 && h > 0);
  return sum / (w * h);
}

void get_2x2_normalized_sses_and_sads(
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

static double get_highbd_diff_mean(const uint8_t *src8, int src_stride,
                                   const uint8_t *dst8, int dst_stride, int w,
                                   int h) {
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

static double get_diff_mean(const uint8_t *src, int src_stride,
                            const uint8_t *dst, int dst_stride, int w, int h) {
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

void model_rd_from_sse(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                       BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                       int num_samples, int *rate, int64_t *dist) {
  (void)num_samples;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const int dequant_shift = (is_cur_buf_hbd(xd)) ? xd->bd - 5 : 3;

  // Fast approximate the modelling function.
  if (cpi->sf.rd_sf.simple_model_rd_from_var) {
    const int64_t square_error = sse;
    int quantizer = p->dequant_QTX[1] >> dequant_shift;
    if (quantizer < 120)
      *rate = (int)AOMMIN(
          (square_error * (280 - quantizer)) >> (16 - AV1_PROB_COST_SHIFT),
          INT_MAX);
    else
      *rate = 0;
    assert(*rate >= 0);
    *dist = (square_error * quantizer) >> 8;
  } else {
    av1_model_rd_from_var_lapndz(sse, num_pels_log2_lookup[plane_bsize],
                                 p->dequant_QTX[1] >> dequant_shift, rate,
                                 dist);
  }
  *dist <<= 4;
}

void model_rd_for_sb(const AV1_COMP *const cpi, BLOCK_SIZE bsize, MACROBLOCK *x,
                     MACROBLOCKD *xd, int plane_from, int plane_to,
                     int *out_rate_sum, int64_t *out_dist_sum,
                     int *skip_txfm_sb, int64_t *skip_sse_sb, int *plane_rate,
                     int64_t *plane_sse, int64_t *plane_dist) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  int plane;
  const int ref = xd->mi[0]->ref_frame[0];

  int64_t rate_sum = 0;
  int64_t dist_sum = 0;
  int64_t total_sse = 0;

  assert(bsize < BLOCK_SIZES_ALL);

  for (plane = plane_from; plane <= plane_to; ++plane) {
    struct macroblock_plane *const p = &x->plane[plane];
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    assert(plane_bsize < BLOCK_SIZES_ALL);
    const int bw = block_size_wide[plane_bsize];
    const int bh = block_size_high[plane_bsize];
    int64_t sse;
    int rate;
    int64_t dist;

    if (x->skip_chroma_rd && plane) continue;

    sse = calculate_sse(xd, p, pd, bw, bh);

    model_rd_from_sse(cpi, x, plane_bsize, plane, sse, bw * bh, &rate, &dist);

    if (plane == 0) x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

    total_sse += sse;
    rate_sum += rate;
    dist_sum += dist;
    if (plane_rate) plane_rate[plane] = rate;
    if (plane_sse) plane_sse[plane] = sse;
    if (plane_dist) plane_dist[plane] = dist;
    assert(rate_sum >= 0);
  }

  if (skip_txfm_sb) *skip_txfm_sb = total_sse == 0;
  if (skip_sse_sb) *skip_sse_sb = total_sse << 4;
  rate_sum = AOMMIN(rate_sum, INT_MAX);
  *out_rate_sum = (int)rate_sum;
  *out_dist_sum = dist_sum;
}

// Fits a curve for rate and distortion using as feature:
// log2(sse_norm/qstep^2)
void model_rd_with_curvfit(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                           BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                           int num_samples, int *rate, int64_t *dist) {
  (void)cpi;
  (void)plane_bsize;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const int dequant_shift = (is_cur_buf_hbd(xd)) ? xd->bd - 5 : 3;
  const int qstep = AOMMAX(p->dequant_QTX[1] >> dequant_shift, 1);

  if (sse == 0) {
    if (rate) *rate = 0;
    if (dist) *dist = 0;
    return;
  }
  aom_clear_system_state();
  const double sse_norm = (double)sse / num_samples;
  const double qstepsqr = (double)qstep * qstep;
  const double xqr = log2(sse_norm / qstepsqr);
  double rate_f, dist_by_sse_norm_f;
  av1_model_rd_curvfit(plane_bsize, sse_norm, xqr, &rate_f,
                       &dist_by_sse_norm_f);

  const double dist_f = dist_by_sse_norm_f * sse_norm;
  int rate_i = (int)(AOMMAX(0.0, rate_f * num_samples) + 0.5);
  int64_t dist_i = (int64_t)(AOMMAX(0.0, dist_f * num_samples) + 0.5);
  aom_clear_system_state();

  // Check if skip is better
  if (rate_i == 0) {
    dist_i = sse << 4;
  } else if (RDCOST(x->rdmult, rate_i, dist_i) >=
             RDCOST(x->rdmult, 0, sse << 4)) {
    rate_i = 0;
    dist_i = sse << 4;
  }

  if (rate) *rate = rate_i;
  if (dist) *dist = dist_i;
}

void model_rd_for_sb_with_curvfit(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                  MACROBLOCK *x, MACROBLOCKD *xd,
                                  int plane_from, int plane_to,
                                  int *out_rate_sum, int64_t *out_dist_sum,
                                  int *skip_txfm_sb, int64_t *skip_sse_sb,
                                  int *plane_rate, int64_t *plane_sse,
                                  int64_t *plane_dist) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  const int ref = xd->mi[0]->ref_frame[0];

  int64_t rate_sum = 0;
  int64_t dist_sum = 0;
  int64_t total_sse = 0;

  for (int plane = plane_from; plane <= plane_to; ++plane) {
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    int64_t dist, sse;
    int rate;

    if (x->skip_chroma_rd && plane) continue;

    int bw, bh;
    const struct macroblock_plane *const p = &x->plane[plane];
    get_txb_dimensions(xd, plane, plane_bsize, 0, 0, plane_bsize, NULL, NULL,
                       &bw, &bh);

    sse = calculate_sse(xd, p, pd, bw, bh);
    model_rd_with_curvfit(cpi, x, plane_bsize, plane, sse, bw * bh, &rate,
                          &dist);

    if (plane == 0) x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

    total_sse += sse;
    rate_sum += rate;
    dist_sum += dist;

    if (plane_rate) plane_rate[plane] = rate;
    if (plane_sse) plane_sse[plane] = sse;
    if (plane_dist) plane_dist[plane] = dist;
  }

  if (skip_txfm_sb) *skip_txfm_sb = rate_sum == 0;
  if (skip_sse_sb) *skip_sse_sb = total_sse << 4;
  *out_rate_sum = (int)rate_sum;
  *out_dist_sum = dist_sum;
}

// Fits a surface for rate and distortion using as features:
// log2(sse_norm + 1) and log2(sse_norm/qstep^2)
void model_rd_with_surffit(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                           BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                           int num_samples, int *rate, int64_t *dist) {
  (void)cpi;
  (void)plane_bsize;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const int dequant_shift = (is_cur_buf_hbd(xd)) ? xd->bd - 5 : 3;
  const int qstep = AOMMAX(p->dequant_QTX[1] >> dequant_shift, 1);
  if (sse == 0) {
    if (rate) *rate = 0;
    if (dist) *dist = 0;
    return;
  }
  aom_clear_system_state();
  const double sse_norm = (double)sse / num_samples;
  const double qstepsqr = (double)qstep * qstep;
  const double xm = log(sse_norm + 1.0) / log(2.0);
  const double yl = log(sse_norm / qstepsqr) / log(2.0);
  double rate_f, dist_by_sse_norm_f;

  av1_model_rd_surffit(plane_bsize, sse_norm, xm, yl, &rate_f,
                       &dist_by_sse_norm_f);

  const double dist_f = dist_by_sse_norm_f * sse_norm;
  int rate_i = (int)(AOMMAX(0.0, rate_f * num_samples) + 0.5);
  int64_t dist_i = (int64_t)(AOMMAX(0.0, dist_f * num_samples) + 0.5);
  aom_clear_system_state();

  // Check if skip is better
  if (rate_i == 0) {
    dist_i = sse << 4;
  } else if (RDCOST(x->rdmult, rate_i, dist_i) >=
             RDCOST(x->rdmult, 0, sse << 4)) {
    rate_i = 0;
    dist_i = sse << 4;
  }

  if (rate) *rate = rate_i;
  if (dist) *dist = dist_i;
}

void model_rd_for_sb_with_surffit(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                  MACROBLOCK *x, MACROBLOCKD *xd,
                                  int plane_from, int plane_to,
                                  int *out_rate_sum, int64_t *out_dist_sum,
                                  int *skip_txfm_sb, int64_t *skip_sse_sb,
                                  int *plane_rate, int64_t *plane_sse,
                                  int64_t *plane_dist) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  const int ref = xd->mi[0]->ref_frame[0];

  int64_t rate_sum = 0;
  int64_t dist_sum = 0;
  int64_t total_sse = 0;

  for (int plane = plane_from; plane <= plane_to; ++plane) {
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    int64_t dist, sse;
    int rate;

    if (x->skip_chroma_rd && plane) continue;

    int bw, bh;
    const struct macroblock_plane *const p = &x->plane[plane];
    get_txb_dimensions(xd, plane, plane_bsize, 0, 0, plane_bsize, NULL, NULL,
                       &bw, &bh);
    sse = calculate_sse(xd, p, pd, bw, bh);

    model_rd_with_surffit(cpi, x, plane_bsize, plane, sse, bw * bh, &rate,
                          &dist);

    if (plane == 0) x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

    total_sse += sse;
    rate_sum += rate;
    dist_sum += dist;

    if (plane_rate) plane_rate[plane] = rate;
    if (plane_sse) plane_sse[plane] = sse;
    if (plane_dist) plane_dist[plane] = dist;
  }

  if (skip_txfm_sb) *skip_txfm_sb = rate_sum == 0;
  if (skip_sse_sb) *skip_sse_sb = total_sse << 4;
  *out_rate_sum = (int)rate_sum;
  *out_dist_sum = dist_sum;
}

void model_rd_with_dnn(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                       BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                       int num_samples, int *rate, int64_t *dist) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const struct macroblock_plane *const p = &x->plane[plane];
  const int log_numpels = num_pels_log2_lookup[plane_bsize];

  const int dequant_shift = (is_cur_buf_hbd(xd)) ? xd->bd - 5 : 3;
  const int q_step = AOMMAX(p->dequant_QTX[1] >> dequant_shift, 1);

  int bw, bh;
  get_txb_dimensions(xd, plane, plane_bsize, 0, 0, plane_bsize, NULL, NULL, &bw,
                     &bh);
  const int src_stride = p->src.stride;
  const uint8_t *const src = p->src.buf;
  const int dst_stride = pd->dst.stride;
  const uint8_t *const dst = pd->dst.buf;
  const int16_t *const src_diff = p->src_diff;
  const int diff_stride = block_size_wide[plane_bsize];
  const int shift = (xd->bd - 8);

  if (sse == 0) {
    if (rate) *rate = 0;
    if (dist) *dist = 0;
    return;
  }
  if (plane) {
    int model_rate;
    int64_t model_dist;
    model_rd_with_curvfit(cpi, x, plane_bsize, plane, sse, num_samples,
                          &model_rate, &model_dist);
    if (rate) *rate = model_rate;
    if (dist) *dist = model_dist;
    return;
  }

  aom_clear_system_state();
  const double sse_norm = (double)sse / num_samples;

  double sse_norm_arr[4];
  get_2x2_normalized_sses_and_sads(cpi, plane_bsize, src, src_stride, dst,
                                   dst_stride, src_diff, diff_stride,
                                   sse_norm_arr, NULL);
  double mean;
  if (is_cur_buf_hbd(xd)) {
    mean = get_highbd_diff_mean(src, src_stride, dst, dst_stride, bw, bh);
  } else {
    mean = get_diff_mean(src, src_stride, dst, dst_stride, bw, bh);
  }
  if (shift) {
    for (int k = 0; k < 4; ++k) sse_norm_arr[k] /= (1 << (2 * shift));
    mean /= (1 << shift);
  }
  double sse_norm_sum = 0.0, sse_frac_arr[3];
  for (int k = 0; k < 4; ++k) sse_norm_sum += sse_norm_arr[k];
  for (int k = 0; k < 3; ++k)
    sse_frac_arr[k] =
        sse_norm_sum > 0.0 ? sse_norm_arr[k] / sse_norm_sum : 0.25;
  const double q_sqr = (double)(q_step * q_step);
  const double q_sqr_by_sse_norm = q_sqr / (sse_norm + 1.0);
  const double mean_sqr_by_sse_norm = mean * mean / (sse_norm + 1.0);
  float hor_corr, vert_corr;
  av1_get_horver_correlation_full(src_diff, diff_stride, bw, bh, &hor_corr,
                                  &vert_corr);

  float features[NUM_FEATURES_PUSTATS];
  features[0] = (float)hor_corr;
  features[1] = (float)log_numpels;
  features[2] = (float)mean_sqr_by_sse_norm;
  features[3] = (float)q_sqr_by_sse_norm;
  features[4] = (float)sse_frac_arr[0];
  features[5] = (float)sse_frac_arr[1];
  features[6] = (float)sse_frac_arr[2];
  features[7] = (float)vert_corr;

  float rate_f, dist_by_sse_norm_f;
  av1_nn_predict(features, &av1_pustats_dist_nnconfig, 1, &dist_by_sse_norm_f);
  av1_nn_predict(features, &av1_pustats_rate_nnconfig, 1, &rate_f);
  aom_clear_system_state();
  const float dist_f = (float)((double)dist_by_sse_norm_f * (1.0 + sse_norm));
  int rate_i = (int)(AOMMAX(0.0, rate_f * num_samples) + 0.5);
  int64_t dist_i = (int64_t)(AOMMAX(0.0, dist_f * num_samples) + 0.5);

  // Check if skip is better
  if (rate_i == 0) {
    dist_i = sse << 4;
  } else if (RDCOST(x->rdmult, rate_i, dist_i) >=
             RDCOST(x->rdmult, 0, sse << 4)) {
    rate_i = 0;
    dist_i = sse << 4;
  }

  if (rate) *rate = rate_i;
  if (dist) *dist = dist_i;
  return;
}

void model_rd_for_sb_with_dnn(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                              MACROBLOCK *x, MACROBLOCKD *xd, int plane_from,
                              int plane_to, int *out_rate_sum,
                              int64_t *out_dist_sum, int *skip_txfm_sb,
                              int64_t *skip_sse_sb, int *plane_rate,
                              int64_t *plane_sse, int64_t *plane_dist) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  const int ref = xd->mi[0]->ref_frame[0];

  int64_t rate_sum = 0;
  int64_t dist_sum = 0;
  int64_t total_sse = 0;

  for (int plane = plane_from; plane <= plane_to; ++plane) {
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    int64_t dist, sse;
    int rate;

    if (x->skip_chroma_rd && plane) continue;

    const struct macroblock_plane *const p = &x->plane[plane];
    int bw, bh;
    get_txb_dimensions(xd, plane, plane_bsize, 0, 0, plane_bsize, NULL, NULL,
                       &bw, &bh);
    sse = calculate_sse(xd, p, pd, bw, bh);

    model_rd_with_dnn(cpi, x, plane_bsize, plane, sse, bw * bh, &rate, &dist);

    if (plane == 0) x->pred_sse[ref] = (unsigned int)AOMMIN(sse, UINT_MAX);

    total_sse += sse;
    rate_sum += rate;
    dist_sum += dist;

    if (plane_rate) plane_rate[plane] = rate;
    if (plane_sse) plane_sse[plane] = sse;
    if (plane_dist) plane_dist[plane] = dist;
  }

  if (skip_txfm_sb) *skip_txfm_sb = rate_sum == 0;
  if (skip_sse_sb) *skip_sse_sb = total_sse << 4;
  *out_rate_sum = (int)rate_sum;
  *out_dist_sum = dist_sum;
}
