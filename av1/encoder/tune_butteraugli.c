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

void av1_setup_butteraugli_recon(AV1_COMP *cpi,
                                 const YV12_BUFFER_CONFIG *recon) {
  YV12_BUFFER_CONFIG *const dst = &cpi->butteraugli_info.recon;
  AV1_COMMON *const cm = &cpi->common;
  const int width = recon->y_crop_width;
  const int height = recon->y_crop_height;
  if (dst->buffer_alloc_sz == 0) {
    aom_alloc_frame_buffer(
        dst, width, height, 1, 1, cm->seq_params.use_highbitdepth,
        cpi->oxcf.border_in_pixels, cm->features.byte_alignment);
  }
  av1_copy_and_extend_frame(recon, dst);
  cpi->butteraugli_info.recon_set = true;
}

void av1_set_mb_butteraugli_rdmult_scaling(AV1_COMP *cpi) {
  if (!cpi->butteraugli_info.recon_set) {
    return;
  }
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  YV12_BUFFER_CONFIG *source = cpi->source;
  const int bit_depth = cpi->td.mb.e_mbd.bd;
  const int resize_factor = 2;
  const int width = source->y_crop_width / resize_factor;
  const int height = source->y_crop_height / resize_factor;

  aom_clear_system_state();
  YV12_BUFFER_CONFIG resized_source, resized_recon;
  memset(&resized_source, 0, sizeof(resized_source));
  memset(&resized_recon, 0, sizeof(resized_recon));
  aom_alloc_frame_buffer(
      &resized_source, width, height, 1, 1, cm->seq_params.use_highbitdepth,
      cpi->oxcf.border_in_pixels, cm->features.byte_alignment);
  aom_alloc_frame_buffer(
      &resized_recon, width, height, 1, 1, cm->seq_params.use_highbitdepth,
      cpi->oxcf.border_in_pixels, cm->features.byte_alignment);
  av1_resize_and_extend_frame_nonnormative(cpi->source, &resized_source,
                                           bit_depth, av1_num_planes(cm));
  av1_resize_and_extend_frame_nonnormative(&cpi->butteraugli_info.recon,
                                           &resized_recon, bit_depth,
                                           av1_num_planes(cm));

  float *diffmap;
  CHECK_MEM_ERROR(cm, diffmap, aom_malloc(width * height * sizeof(*diffmap)));
  if (!aom_calc_butteraugli(&resized_source, &resized_recon, bit_depth,
                            diffmap)) {
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
  double log_sum = 0.0;
  double blk_count = 0.0;

  // Loop through each block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      const int y_start = row * block_h;
      const int x_start = col * block_w;
      float dbutteraugli = 0.0f;
      float dmse = 0.0f;

      // Loop through each pixel.
      for (int y = y_start; y < y_start + block_h && y < height; y++) {
        for (int x = x_start; x < x_start + block_w && x < width; x++) {
          dbutteraugli += powf(diffmap[y * width + x], 6.0f);
          float px_diff =
              resized_source.y_buffer[y * resized_source.y_stride + x] -
              resized_recon.y_buffer[y * resized_recon.y_stride + x];
          dmse += px_diff * px_diff;
        }
      }
      for (int y = y_start; y < y_start + block_h && y < height; y += 2) {
        for (int x = x_start; x < x_start + block_w && x < width; x += 2) {
          const int src_px_index = y / 2 * resized_source.uv_stride + x / 2;
          const int recon_px_index = y / 2 * resized_recon.uv_stride + x / 2;
          const float px_diff_u =
              (float)(resized_source.u_buffer[src_px_index] -
                      resized_recon.u_buffer[recon_px_index]);
          const float px_diff_v =
              (float)(resized_source.v_buffer[src_px_index] -
                      resized_recon.v_buffer[recon_px_index]);
          dmse += px_diff_u * px_diff_u + px_diff_v * px_diff_v;
        }
      }

      dbutteraugli = powf(dbutteraugli, 1.0f / 6.0f);
      dmse = dmse / (2.0f * (float)block_w * (float)block_h);
      // 'K' is used to balance the rate-distortion distribution between PSNR
      // and Butteraugli.
      const double K = 0.1;
      const float eps = 0.01f;
      double weight;
      if (dbutteraugli < eps || dmse < eps) {
        weight = -1.0;
      } else {
        blk_count += 1.0;
        weight = dmse / dbutteraugli;
        weight = AOMMIN(weight, 3.0);
        weight += K;
        log_sum += log(weight);
      }
      cpi->butteraugli_info.rdmult_scaling_factors[index] = weight;
    }
  }
  // Geometric average of the weights.
  log_sum = exp(log_sum / blk_count);

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      double *weight = &cpi->butteraugli_info.rdmult_scaling_factors[index];
      if (*weight <= 0.0) {
        *weight = 1.0;
      } else {
        *weight /= log_sum;
      }
    }
  }

  aom_free(diffmap);
  aom_free_frame_buffer(&resized_source);
  aom_free_frame_buffer(&resized_recon);
  aom_clear_system_state();
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
