/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <immintrin.h>
#include <math.h>

#include "aom_dsp/aom_dsp_common.h"
#include "av1/common/av1_common_int.h"
#include "av1/encoder/cnn.h"

// This mask rearranges source pixels in the order shown below.
// shuffle_src_layer0[0][8]: applied on source pixels 0 to 7.
// shuffle_src_layer0[1][8]: applied on source pixels 7 to 14.
// This shuffling is needed to process 3 5x5 blocks which needs
// source pixels in the following order.
// 1st 5x5 block: source pixels needed are 0 to 4,
// 2nd 5x5 block: source pixels needed are 4 to 8,
// 3rd 5x5 block: source pixels needed are 8 to 12.
const uint32_t shuffle_src_layer0[2][8] = { { 0, 1, 2, 3, 4, 4, 5, 6 },
                                            { 0, 1, 1, 2, 3, 4, 5, 0 } };

// This mask rearranges weights to match source pixels order.
const uint32_t shuffle_weight_layer0[2][8] = { { 0, 1, 2, 3, 4, 0, 1, 2 },
                                               { 3, 4, 0, 1, 2, 3, 4, 0 } };

// Load weights needed for layer 0 (for 5x5 block processing),
// and fill the registers appropriately to match source pixel mapping.
#define PREPARE_WEIGHTS_FOR_5X5_CONV()                              \
  {                                                                 \
    for (int col = 0; col < 5; ++col) {                             \
      for (int row = 0; row < 5; ++row) {                           \
        weight[col][row] = layer_config->weights[off];              \
        off += cstep;                                               \
      }                                                             \
    }                                                               \
                                                                    \
    shuffle_weight[0] = _mm256_loadu_ps(&weight[0][0]);             \
    shuffle_weight[1] = _mm256_loadu_ps(&weight[1][0]);             \
    shuffle_weight[2] = _mm256_loadu_ps(&weight[2][0]);             \
    shuffle_weight[3] = _mm256_loadu_ps(&weight[3][0]);             \
    shuffle_weight[4] = _mm256_loadu_ps(&weight[4][0]);             \
                                                                    \
    shuffle_weight[0] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[0], weight_mask_0); \
    shuffle_weight[1] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[1], weight_mask_0); \
    shuffle_weight[2] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[2], weight_mask_0); \
    shuffle_weight[3] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[3], weight_mask_0); \
    shuffle_weight[4] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[4], weight_mask_0); \
    shuffle_weight[5] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[0], weight_mask_1); \
    shuffle_weight[6] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[1], weight_mask_1); \
    shuffle_weight[7] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[2], weight_mask_1); \
    shuffle_weight[8] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[3], weight_mask_1); \
    shuffle_weight[9] =                                             \
        _mm256_permutevar8x32_ps(shuffle_weight[4], weight_mask_1); \
  }

// For each row, loads source pixels 0 to 7(load_src_0), 7 to 14(load_src_1) and
// arranges them appropritely to prcoess 3 blocks.
#define PERFORM_CONVOLVE_FOR_3_5X5_BLOCKS()                            \
  {                                                                    \
    for (int row = 0; row < 5; row++) {                                \
      load_src_0 = _mm256_loadu_ps(input_ptr);                         \
      load_src_1 = _mm256_loadu_ps(input_ptr + 7);                     \
      load_src_0 = _mm256_permutevar8x32_ps(load_src_0, block0_1);     \
      load_src_1 = _mm256_permutevar8x32_ps(load_src_1, block1_2);     \
      load_src_0 = _mm256_mul_ps(load_src_0, shuffle_weight[0 + row]); \
      load_src_1 = _mm256_mul_ps(load_src_1, shuffle_weight[5 + row]); \
      accum_src_0 = _mm256_add_ps(load_src_0, accum_src_0);            \
      accum_src_1 = _mm256_add_ps(load_src_1, accum_src_1);            \
      input_ptr += in_stride;                                          \
    }                                                                  \
  }

// Do convolution of one 5x5 block.
#define PERFORM_CONVOLVE_FOR_1_5X5_BLOCK(blk0, w, accum0, in_stride)     \
  {                                                                      \
    __m128 load_src[5];                                                  \
    load_src[0] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[0][4];                      \
    input_ptr += in_stride;                                              \
    load_src[1] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[1][4];                      \
    input_ptr += in_stride;                                              \
    load_src[2] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[2][4];                      \
    input_ptr += in_stride;                                              \
    load_src[3] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[3][4];                      \
    input_ptr += in_stride;                                              \
    load_src[4] = _mm_loadu_ps(input_ptr);                               \
    last_column_sum += input_ptr[4] * weight[4][4];                      \
                                                                         \
    load_src[0] = _mm_mul_ps(load_src[0], _mm256_castps256_ps128(w[0])); \
    load_src[1] = _mm_mul_ps(load_src[1], _mm256_castps256_ps128(w[1])); \
    load_src[2] = _mm_mul_ps(load_src[2], _mm256_castps256_ps128(w[2])); \
    load_src[3] = _mm_mul_ps(load_src[3], _mm256_castps256_ps128(w[3])); \
    load_src[4] = _mm_mul_ps(load_src[4], _mm256_castps256_ps128(w[4])); \
                                                                         \
    accum0 = _mm_add_ps(load_src[0], accum0);                            \
    load_src[1] = _mm_add_ps(load_src[1], load_src[2]);                  \
    load_src[3] = _mm_add_ps(load_src[3], load_src[4]);                  \
    load_src[1] = _mm_add_ps(load_src[1], load_src[3]);                  \
    accum0 = _mm_add_ps(accum0, load_src[1]);                            \
  }

// Perform convolution when filter_width and filter_height are equal to 5.
static void no_maxpool_padding_valid_convolve_5x5_avx2(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int channel_step, const int cstep) {
  assert(layer_config->filter_width == 5 && layer_config->filter_height == 5);

  // Load shuffle buffers needed for source.
  const __m256i block0_1 = _mm256_loadu_si256((__m256i *)shuffle_src_layer0[0]);
  const __m256i block1_2 = _mm256_loadu_si256((__m256i *)shuffle_src_layer0[1]);

  // Load shuffle buffers needed for weight.
  const __m256i weight_mask_0 =
      _mm256_loadu_si256((__m256i *)shuffle_weight_layer0[0]);
  const __m256i weight_mask_1 =
      _mm256_loadu_si256((__m256i *)shuffle_weight_layer0[1]);

  // Width needed to process 3 5x5 blocks.
  const int process_width = layer_config->skip_width * 3;
  for (int i = start_idx; i < layer_config->out_channels; i += channel_step) {
    const float out_ch_bias = layer_config->bias[i];
    for (int k = 0; k < layer_config->in_channels; ++k) {
      __m256 shuffle_weight[10];
      float weight[5][8] = { { 0 } };
      int off = k * layer_config->out_channels + i;

      // In layer 0, the convolution process happens at 5x5.
      // The weights needed for 5x5 block are same across the in-channels,
      // which is why the load of weights happens once for each in-channel.
      PREPARE_WEIGHTS_FOR_5X5_CONV();

      for (int h = 0, u = 0; h < in_height - layer_config->filter_height + 1;
           h += layer_config->skip_height, ++u) {
        const int out_h = u * out_stride;
        int v = 0;
        int w = 0;
        int rem_width = in_width;
        // Processing 3 5x5 blocks at a time, if sufficient width is present.
        while (rem_width > process_width) {
          __m256 load_src_0, load_src_1;
          __m256 accum_src_0 = _mm256_setzero_ps();
          __m256 accum_src_1 = _mm256_setzero_ps();
          const float *input_ptr = &input[k][h * in_stride + w];
          PERFORM_CONVOLVE_FOR_3_5X5_BLOCKS();

          // Accumulate across column.
          __m256 accum = _mm256_hadd_ps(accum_src_0, accum_src_1);
          __m128 tmp_reg_0 = _mm256_extractf128_ps(accum_src_0, 1);
          __m128 tmp_reg_1 = _mm256_extractf128_ps(accum_src_1, 1);

          __m128 accum_l = _mm256_castps256_ps128(accum);
          __m128 accum_h = _mm256_extractf128_ps(accum, 1);

          __m128 tmp_reg_2 = _mm_add_ps(accum_l, tmp_reg_0);
          __m128 tmp_reg_3 = _mm_add_ps(tmp_reg_0, accum_h);
          __m128 tmp_reg_4 = _mm_add_ps(tmp_reg_1, accum_h);

          // 1st 5x5 block output.
          output[i][out_h + v] =
              out_ch_bias + _mm_cvtss_f32(tmp_reg_2) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 1));

          // 2nd 5x5 block output.
          output[i][out_h + v + 1] =
              out_ch_bias +
              _mm_cvtss_f32(_mm_shuffle_ps(tmp_reg_3, tmp_reg_3, 1)) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 2));

          // 3rd 5x5 block output.
          output[i][out_h + v + 2] =
              out_ch_bias + _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 3)) +
              _mm_cvtss_f32(_mm_shuffle_ps(tmp_reg_4, tmp_reg_4, 2));

          v += 3;
          w += process_width;
          rem_width -= process_width;
        }

        // Process remaining blocks as single 5x5 block at a time.
        while (rem_width > layer_config->skip_width) {
          float last_column_sum = 0;
          __m128 accum = _mm_setzero_ps();
          const float *input_ptr = &input[k][h * in_stride + w];
          PERFORM_CONVOLVE_FOR_1_5X5_BLOCK(block0_1, shuffle_weight, accum,
                                           in_stride);

          // Accumulate across column.
          accum = _mm_hadd_ps(accum, accum);
          output[i][out_h + v] = out_ch_bias + last_column_sum +
                                 _mm_cvtss_f32(accum) +
                                 _mm_cvtss_f32(_mm_shuffle_ps(accum, accum, 1));

          v += 1;
          w += layer_config->skip_width;
          rem_width -= layer_config->skip_width;
        }
      }
    }
  }
}

// AVX2 variant of av1_cnn_convolve_c().
void av1_cnn_convolve_avx2(const float **input, int in_width, int in_height,
                           int in_stride, const CNN_LAYER_CONFIG *layer_config,
                           float **output, int out_stride, int start_idx,
                           int step) {
  assert(!layer_config->deconvolve);
  const int cstep = layer_config->in_channels * layer_config->out_channels;
  const int filter_height_half = layer_config->filter_height >> 1;
  const int filter_width_half = layer_config->filter_width >> 1;
  const int channel_step = AOMMAX(step, 1);

  if (layer_config->maxpool &&
      (layer_config->skip_height > 1 || layer_config->skip_width > 1)) {
    switch (layer_config->pad) {
      case PADDING_SAME_ZERO:
        maxpool_padding_zero(input, in_width, in_height, in_stride,
                             layer_config, output, out_stride, cstep,
                             filter_width_half, filter_height_half);
        break;
      case PADDING_SAME_REPLICATE:
        maxpool_padding_replicate(input, in_width, in_height, in_stride,
                                  layer_config, output, out_stride, cstep,
                                  filter_width_half, filter_height_half);
        break;
      case PADDING_VALID:
        maxpool_padding_valid(input, in_width, in_height, in_stride,
                              layer_config, output, out_stride, cstep);
        break;
      default: assert(0 && "Unknown padding type");
    }
  } else {
    // Results in element-wise matrix multiplication.
    if (layer_config->filter_height == 1 && layer_config->filter_width == 1) {
      element_wise_convolve(input, in_width, in_height, in_stride, layer_config,
                            output, out_stride, start_idx, step);
      return;
    }
    const int ii_shift =
        filter_height_half - (layer_config->filter_height - 1) % 2;
    const int jj_shift =
        filter_width_half - (layer_config->filter_width - 1) % 2;
    switch (layer_config->pad) {
      case PADDING_SAME_ZERO: {
        no_maxpool_padding_zero(input, in_width, in_height, in_stride,
                                layer_config, output, out_stride, start_idx,
                                cstep, filter_width_half, filter_height_half,
                                ii_shift, jj_shift, channel_step);
        break;
      }
      case PADDING_SAME_REPLICATE: {
        no_maxpool_padding_replicate(
            input, in_width, in_height, in_stride, layer_config, output,
            out_stride, start_idx, cstep, ii_shift, jj_shift, channel_step);
        break;
      }
      case PADDING_VALID: {
        if (layer_config->filter_width == 5 &&
            layer_config->filter_height == 5) {
          no_maxpool_padding_valid_convolve_5x5_avx2(
              input, in_width, in_height, in_stride, layer_config, output,
              out_stride, start_idx, channel_step, cstep);
        } else {
          no_maxpool_padding_valid(input, in_width, in_height, in_stride,
                                   layer_config, output, out_stride, start_idx,
                                   cstep, channel_step);
        }
        break;
      }
      default: assert(0 && "Unknown padding type");
    }
  }
}
