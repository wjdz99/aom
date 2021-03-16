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

#include <math.h>

#include "av1/encoder/tune_butteraugli.h"

#include "aom_dsp/butteraugli.h"
#include "aom_ports/system_state.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/extend.h"

static const int resize_factor = 2;

struct block_rd_info {
  double alpha;
  double var;
  double mse;
  double dbutteraugli;
};

static void update_scaling_factors(const int blk_count,
                                   const struct block_rd_info *const rds,
                                   double *const new_scaling_facotrs,
                                   double k) {
  double log_sum = 0.0;
  double non_zero_blk_count = 0.0;
  for (int i = 0; i < blk_count; ++i) {
    const float eps = 0.1f;
    double weight;
    if (rds[i].dbutteraugli < eps || rds[i].mse < eps) {
      weight = -1.0;
    } else {
      weight = rds[i].mse / rds[i].dbutteraugli;
      weight = AOMMIN(weight, 5.0);
      weight = AOMMAX(weight, 0.2);
      weight += k;
      log_sum += log(weight);
      non_zero_blk_count += 1.0;
    }
    new_scaling_facotrs[i] = weight;
  }

  // Geometric average of the weights.
  log_sum = exp(log_sum / non_zero_blk_count);

  for (int i = 0; i < blk_count; ++i) {
    if (new_scaling_facotrs[i] > 0.0) {
      new_scaling_facotrs[i] /= log_sum;
      new_scaling_facotrs[i] = AOMMIN(new_scaling_facotrs[i], 2.5);
      new_scaling_facotrs[i] = AOMMAX(new_scaling_facotrs[i], 0.4);
    } else {
      new_scaling_facotrs[i] = 1.0;
    }
  }
}

static double solve_baseline_frame_rate(const int blk_count,
                                        const struct block_rd_info *const rds,
                                        const double orig_rdmult,
                                        const double target_distortion) {
  double rdmult_l = orig_rdmult, rdmult_r = orig_rdmult * 4.0;
  int iter = 20;
  double distortion = 0.0, rdmult = 0.0;
  while (iter > 0 && fabs(distortion - target_distortion) > blk_count * 0.001) {
    distortion = 0.0;
    rdmult = (rdmult_l + rdmult_r) / 2.0;
    for (int i = 0; i < blk_count; ++i) distortion += rdmult * rds[i].alpha;
    if (distortion > target_distortion) rdmult_r = rdmult;
    if (distortion < target_distortion) rdmult_l = rdmult;
    iter--;
  }

  const float eps = 0.1f;
  double frame_rate = 0.0;
  for (int i = 0; i < blk_count; ++i) {
    if (rds[i].dbutteraugli > eps && rds[i].mse > eps && rds[i].var > eps) {
      const double d = rdmult * rds[i].alpha;
      const double r = rds[i].alpha * log(rds[i].var / d);
      frame_rate += AOMMAX(0.0, r);
    }
  }
  return frame_rate;
}

static double find_optimal_k(AV1_COMMON *const cm, const int blk_count,
                             const double distortion_limit,
                             const struct block_rd_info *rds,
                             const double orig_rdmult,
                             const double baseline_distortion) {
  const float eps = 0.1f;
  const double step_size = 0.1;
  const double max_k = 3.0;
  double k = 0.1, baseline_frame_rate, frame_rate, *new_scaling_factors;
  CHECK_MEM_ERROR(cm, new_scaling_factors,
                  aom_malloc(blk_count * sizeof(*new_scaling_factors)));
  double frame_distortion;

  do {
    frame_distortion = 0.0;
    frame_rate = 0.0;
    k += step_size;
    update_scaling_factors(blk_count, rds, new_scaling_factors, k);
    // Solve frame_rate and distortion with scaling_factors.
    for (int i = 0; i < blk_count; ++i) {
      const double rdmult = new_scaling_factors[i] * orig_rdmult;
      const double d = rdmult * rds[i].alpha;
      frame_distortion += d;
      if (rds[i].dbutteraugli > eps && rds[i].mse > eps && rds[i].var > eps) {
        const double r = rds[i].alpha * log(rds[i].var / d);
        frame_rate += r;
      }
    }

    baseline_frame_rate = solve_baseline_frame_rate(blk_count, rds, orig_rdmult,
                                                    frame_distortion);
  } while (k < max_k &&
           frame_rate > distortion_limit * baseline_frame_rate);
  return k;
}

static void set_mb_butteraugli_rdmult_scaling(AV1_COMP *cpi,
                                              const YV12_BUFFER_CONFIG *source,
                                              const YV12_BUFFER_CONFIG *recon) {
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int bit_depth = cpi->td.mb.e_mbd.bd;
  const int width = source->y_crop_width;
  const int height = source->y_crop_height;

  float *diffmap;
  CHECK_MEM_ERROR(cm, diffmap, aom_malloc(width * height * sizeof(*diffmap)));
  if (!aom_calc_butteraugli(source, recon, bit_depth, diffmap)) {
    aom_internal_error(&cm->error, AOM_CODEC_ERROR,
                       "Failed to calculate Butteraugli distances.");
  }

  const int num_mi_w = mi_size_wide[butteraugli_rdo_bsize] / resize_factor;
  const int num_mi_h = mi_size_high[butteraugli_rdo_bsize] / resize_factor;
  const int num_cols =
      (mi_params->mi_cols / resize_factor + num_mi_w - 1) / num_mi_w;
  const int num_rows =
      (mi_params->mi_rows / resize_factor + num_mi_h - 1) / num_mi_h;
  const int block_w = num_mi_w << 2;
  const int block_h = num_mi_h << 2;

  const double rdmult = (double)cpi->rd.RDMULT;
  struct block_rd_info *rds;
  CHECK_MEM_ERROR(cm, rds, aom_malloc(num_rows * num_cols * sizeof(*rds)));
  double *const scaling_factors = cpi->butteraugli_info.rdmult_scaling_factors;
  double baseline_distortion = 0.0;

  // for (int i = 0; i < num_rows * num_cols; ++i) {
  //   diffmap[i] = (float)i / (num_rows * num_cols / 2) + 0.01f;
  // }

  // Loop through each block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      const int y_start = row * block_h;
      const int x_start = col * block_w;
      float dbutteraugli = 0.0f, dmse = 0.0f, src_px_d = 0.0f, src_px_sd = 0.0f;
      float px_count = 0.0f;

      // Loop through each pixel.
      for (int y = y_start; y < y_start + block_h && y < height; y++) {
        for (int x = x_start; x < x_start + block_w && x < width; x++) {
          dbutteraugli += powf(diffmap[y * width + x], 12.0f);
          const float src_y = source->y_buffer[y * source->y_stride + x];
          const float px_diff =
              src_y - (float)recon->y_buffer[y * recon->y_stride + x];
          dmse += px_diff * px_diff;
          src_px_d += src_y;
          src_px_sd += src_y * src_y;
          px_count += 1.0f;
        }
      }
      for (int y = y_start / 2; y < y_start / 2 + block_h / 2 && y < height / 2;
           y++) {
        for (int x = x_start / 2;
             x < x_start / 2 + block_w / 2 && x < width / 2; x++) {
          const int src_px_index = y * source->uv_stride + x;
          const int recon_px_index = y * recon->uv_stride + x;
          const float src_u = source->u_buffer[src_px_index];
          const float src_v = source->v_buffer[src_px_index];
          const float px_diff_u =
              src_u - (float)recon->u_buffer[recon_px_index];
          const float px_diff_v =
              src_v - (float)recon->v_buffer[recon_px_index];
          dmse += px_diff_u * px_diff_u + px_diff_v * px_diff_v;
          src_px_d += src_u + src_v;
          src_px_sd += src_u * src_u + src_v * src_v;
          px_count += 2.0f;
        }
      }

      baseline_distortion += dmse;
      dbutteraugli = powf(dbutteraugli, 1.0f / 12.0f);
      rds[index].alpha = dmse / (double)rdmult;
      dmse = dmse / px_count;
      rds[index].var = src_px_sd / px_count - powf(src_px_d / px_count, 2.0f);
      rds[index].mse = dmse;
      rds[index].dbutteraugli = dbutteraugli;
    }
  }

  const double distortion_limit = 1.5;
  // 'K' is used to balance the rate-distortion distribution between PSNR
  // and Butteraugli.
  const double k = find_optimal_k(cm, num_rows * num_cols, distortion_limit,
                                  rds, rdmult, baseline_distortion);

  update_scaling_factors(num_rows * num_cols, rds, scaling_factors, k);

  aom_free(diffmap);
  aom_free(rds);
}

void av1_set_butteraugli_rdmult(const AV1_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE bsize, int mi_row, int mi_col,
                                int *rdmult) {
  assert(cpi->oxcf.tune_cfg.tuning == AOM_TUNE_BUTTERAUGLI);
  if (!cpi->butteraugli_info.recon_set) {
    return;
  }
  const AV1_COMMON *const cm = &cpi->common;

  const int num_mi_w = mi_size_wide[butteraugli_rdo_bsize];
  const int num_mi_h = mi_size_high[butteraugli_rdo_bsize];
  const int num_cols = (cm->mi_params.mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (cm->mi_params.mi_rows + num_mi_h - 1) / num_mi_h;
  const int num_bcols = (mi_size_wide[bsize] + num_mi_w - 1) / num_mi_w;
  const int num_brows = (mi_size_high[bsize] + num_mi_h - 1) / num_mi_h;
  double num_of_mi = 0.0;
  double geom_mean_of_scale = 0.0;

  aom_clear_system_state();
  for (int row = mi_row / num_mi_w;
       row < num_rows && row < mi_row / num_mi_w + num_brows; ++row) {
    for (int col = mi_col / num_mi_h;
         col < num_cols && col < mi_col / num_mi_h + num_bcols; ++col) {
      const int index = row * num_cols + col;
      geom_mean_of_scale +=
          log(cpi->butteraugli_info.rdmult_scaling_factors[index]);
      num_of_mi += 1.0;
    }
  }
  geom_mean_of_scale = exp(geom_mean_of_scale / num_of_mi);

  *rdmult = (int)((double)(*rdmult) * geom_mean_of_scale + 0.5);
  *rdmult = AOMMAX(*rdmult, 0);
  av1_set_error_per_bit(&x->errorperbit, *rdmult);
  aom_clear_system_state();
}

static void copy_plane(const uint8_t *src, int src_stride, uint8_t *dst,
                       int dst_stride, int w, int h) {
  for (int row = 0; row < h; row++) {
    memcpy(dst, src, w);
    src += src_stride;
    dst += dst_stride;
  }
}

static void copy_img(const YV12_BUFFER_CONFIG *src, YV12_BUFFER_CONFIG *dst,
                     int width, int height) {
  copy_plane(src->y_buffer, src->y_stride, dst->y_buffer, dst->y_stride, width,
             height);
  copy_plane(src->u_buffer, src->uv_stride, dst->u_buffer, dst->uv_stride,
             width / 2, height / 2);
  copy_plane(src->v_buffer, src->uv_stride, dst->v_buffer, dst->uv_stride,
             width / 2, height / 2);
}

static void zero_plane(uint8_t *dst, int dst_stride, int h) {
  for (int row = 0; row < h; row++) {
    memset(dst, 0, dst_stride);
    dst += dst_stride;
  }
}

static void zero_img(YV12_BUFFER_CONFIG *dst) {
  zero_plane(dst->y_buffer, dst->y_stride, dst->y_height);
  zero_plane(dst->u_buffer, dst->uv_stride, dst->uv_height);
  zero_plane(dst->v_buffer, dst->uv_stride, dst->uv_height);
}

void av1_setup_butteraugli_source(AV1_COMP *cpi) {
  aom_clear_system_state();
  YV12_BUFFER_CONFIG *const dst = &cpi->butteraugli_info.source;
  AV1_COMMON *const cm = &cpi->common;
  const int width = cpi->source->y_crop_width;
  const int height = cpi->source->y_crop_height;
  const int bit_depth = cpi->td.mb.e_mbd.bd;
  if (dst->buffer_alloc_sz == 0) {
    aom_alloc_frame_buffer(
        dst, width, height, 1, 1, cm->seq_params.use_highbitdepth,
        cpi->oxcf.border_in_pixels, cm->features.byte_alignment);
  }
  av1_copy_and_extend_frame(cpi->source, dst);

  YV12_BUFFER_CONFIG *const resized_dst = &cpi->butteraugli_info.resized_source;
  if (resized_dst->buffer_alloc_sz == 0) {
    aom_alloc_frame_buffer(
        resized_dst, width / resize_factor, height / resize_factor, 1, 1,
        cm->seq_params.use_highbitdepth, cpi->oxcf.border_in_pixels,
        cm->features.byte_alignment);
  }
  av1_resize_and_extend_frame_nonnormative(cpi->source, resized_dst, bit_depth,
                                           av1_num_planes(cm));

  zero_img(cpi->source);
  copy_img(resized_dst, cpi->source, width / resize_factor,
           height / resize_factor);
  aom_clear_system_state();
}

void av1_restore_butteraugli_source(AV1_COMP *cpi) {
  aom_clear_system_state();
  av1_copy_and_extend_frame(&cpi->butteraugli_info.source, cpi->source);
  AV1_COMMON *const cm = &cpi->common;
  const int width = cpi->source->y_crop_width;
  const int height = cpi->source->y_crop_height;

  YV12_BUFFER_CONFIG resized_recon;
  memset(&resized_recon, 0, sizeof(resized_recon));
  aom_alloc_frame_buffer(
      &resized_recon, width / resize_factor, height / resize_factor, 1, 1,
      cm->seq_params.use_highbitdepth, cpi->oxcf.border_in_pixels,
      cm->features.byte_alignment);
  copy_img(&cpi->common.cur_frame->buf, &resized_recon, width / resize_factor,
           height / resize_factor);

  set_mb_butteraugli_rdmult_scaling(cpi, &cpi->butteraugli_info.resized_source,
                                    &resized_recon);
  cpi->butteraugli_info.recon_set = true;
  aom_free_frame_buffer(&resized_recon);
  aom_clear_system_state();
}
