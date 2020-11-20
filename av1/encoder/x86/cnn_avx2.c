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

#define CLAMPINDEX(a, hi) ((a) < 0 ? 0 : ((a) >= (hi) ? ((hi)-1) : (a)))

// This mask rearranges source pixels in the order shown below.
// shuffle_src_layer0[0][8]: applied on source pixels 0 to 7.
// shuffle_src_layer0[1][8]: applied on source pixels 7 to 14.
// This shuffling is needed to process 3 5x5 blocks which needs
// source pixels in the following order.
// 1st 5x5 block: source pixels needed are 0 to 4,
// 2nd 5x5 block: source pixels needed are 4 to 8,
// 3rd 5x5 block: source pixels needed are 8 to 12.
const int shuffle_src_layer0[2][8] = { { 0, 1, 2, 3, 4, 4, 5, 6 },
                                       { 0, 1, 1, 2, 3, 4, 5, 0 } };

// This mask rearranges weights to match source pixels order.
const int shuffle_weight_layer0[2][8] = { { 0, 1, 2, 3, 4, 0, 1, 2 },
                                          { 3, 4, 0, 1, 2, 3, 4, 0 } };

// Load weights needed for layer 0 (for 5x5 block processing),
// and fill the registers appropriately to match source pixel mapping.
static INLINE void set_5x5_weights(__m256 load_weight[], __m256i weight_mask_0,
                                   __m256i weight_mask_1,
                                   const CNN_LAYER_CONFIG *layer_config,
                                   int off, int cstep) {
  float weight[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      weight[i][j] = layer_config->weights[off];
      off += cstep;
    }
  }

  load_weight[0] = _mm256_loadu_ps(&weight[0][0]);
  load_weight[1] = _mm256_loadu_ps(&weight[1][0]);
  load_weight[2] = _mm256_loadu_ps(&weight[2][0]);
  load_weight[3] = _mm256_loadu_ps(&weight[3][0]);
  load_weight[4] = _mm256_loadu_ps(&weight[4][0]);
  // set weight according to block mapping.
  load_weight[0] = _mm256_permutevar8x32_ps(load_weight[0], weight_mask_0);
  load_weight[1] = _mm256_permutevar8x32_ps(load_weight[1], weight_mask_0);
  load_weight[2] = _mm256_permutevar8x32_ps(load_weight[2], weight_mask_0);
  load_weight[3] = _mm256_permutevar8x32_ps(load_weight[3], weight_mask_0);
  load_weight[4] = _mm256_permutevar8x32_ps(load_weight[4], weight_mask_0);
  load_weight[5] = _mm256_permutevar8x32_ps(load_weight[0], weight_mask_1);
  load_weight[6] = _mm256_permutevar8x32_ps(load_weight[1], weight_mask_1);
  load_weight[7] = _mm256_permutevar8x32_ps(load_weight[2], weight_mask_1);
  load_weight[8] = _mm256_permutevar8x32_ps(load_weight[3], weight_mask_1);
  load_weight[9] = _mm256_permutevar8x32_ps(load_weight[4], weight_mask_1);
}

// For each row, loads source pixels 0 to 7(load_src_0), 7 to 14(load_src_1) and
// arranges them appropritely to prcoess 3 blocks.
#define DO_3_5X5_BLOCK_CONVOLVE()                                      \
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
#define DO_1_5X5_BLOCK_CONVOLVE(blk0, w, accum0, input, in_stride) \
  {                                                                \
    __m256 load_src[5];                                            \
    load_src[0] = _mm256_loadu_ps(input);                          \
    input += in_stride;                                            \
    load_src[1] = _mm256_loadu_ps(input);                          \
    input_ptr += in_stride;                                        \
    load_src[2] = _mm256_loadu_ps(input);                          \
    input_ptr += in_stride;                                        \
    load_src[3] = _mm256_loadu_ps(input);                          \
    input_ptr += in_stride;                                        \
    load_src[4] = _mm256_loadu_ps(input);                          \
                                                                   \
    load_src[0] = _mm256_permutevar8x32_ps(load_src[0], blk0);     \
    load_src[1] = _mm256_permutevar8x32_ps(load_src[1], blk0);     \
    load_src[2] = _mm256_permutevar8x32_ps(load_src[2], blk0);     \
    load_src[3] = _mm256_permutevar8x32_ps(load_src[3], blk0);     \
    load_src[4] = _mm256_permutevar8x32_ps(load_src[4], blk0);     \
                                                                   \
    load_src[0] = _mm256_mul_ps(load_src[0], w[0]);                \
    load_src[1] = _mm256_mul_ps(load_src[1], w[1]);                \
    load_src[2] = _mm256_mul_ps(load_src[2], w[2]);                \
    load_src[3] = _mm256_mul_ps(load_src[3], w[3]);                \
    load_src[4] = _mm256_mul_ps(load_src[4], w[4]);                \
                                                                   \
    accum0 = _mm256_add_ps(load_src[0], accum0);                   \
    load_src[1] = _mm256_add_ps(load_src[1], load_src[2]);         \
    load_src[3] = _mm256_add_ps(load_src[3], load_src[4]);         \
    load_src[1] = _mm256_add_ps(load_src[1], load_src[3]);         \
    accum0 = _mm256_add_ps(accum0, load_src[1]);                   \
  }

// Perform convolution when filter_width and filter_height are equal to 5.
void cnn_convlove_5x5_avx2(const float **input, int in_width, int in_height,
                           int in_stride, const CNN_LAYER_CONFIG *layer_config,
                           float **output, int out_stride, int start_idx,
                           const int channel_step, const int cstep) {
  const __m256i block0_1 = _mm256_loadu_si256((__m256i *)shuffle_src_layer0[0]);
  const __m256i block1_2 = _mm256_loadu_si256((__m256i *)shuffle_src_layer0[1]);

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
      // In layer 0, the convolution process happens at 5x5.
      // The weights needed for 5x5 block are same across the in-channels,
      // which is why the load of weights happens once for each in-channel.
      set_5x5_weights(shuffle_weight, weight_mask_0, weight_mask_1,
                      layer_config, k * layer_config->out_channels + i, cstep);

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
          DO_3_5X5_BLOCK_CONVOLVE();

          // Accumulate across column
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

        // Process remaining blocks as single 5x5 block.
        while (rem_width > layer_config->skip_width) {
          __m256 accum = _mm256_setzero_ps();
          const float *input_ptr = &input[k][h * in_stride + w];
          DO_1_5X5_BLOCK_CONVOLVE(block0_1, shuffle_weight, accum, input_ptr,
                                  in_stride);

          // Accumulate across column
          __m128 accum_l = _mm256_castps256_ps128(accum);
          __m128 accum_h = _mm256_extractf128_ps(accum, 1);
          accum_l = _mm_hadd_ps(accum_l, accum_h);
          output[i][out_h + v] =
              out_ch_bias + _mm_cvtss_f32(accum_l) +
              _mm_cvtss_f32(_mm_shuffle_ps(accum_l, accum_l, 1)) +
              _mm_cvtss_f32(accum_h);

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
        for (int i = 0; i < layer_config->out_channels; ++i) {
          for (int h = 0, u = 0; h < in_height;
               h += layer_config->skip_height, ++u) {
            for (int w = 0, v = 0; w < in_width;
                 w += layer_config->skip_width, ++v) {
              for (int hh = h;
                   hh < AOMMIN(in_height, h + layer_config->skip_height);
                   ++hh) {
                for (int ww = w;
                     ww < AOMMIN(in_width, w + layer_config->skip_width);
                     ++ww) {
                  float sum = layer_config->bias[i];
                  for (int k = 0; k < layer_config->in_channels; ++k) {
                    int off = k * layer_config->out_channels + i;
                    for (int l = 0; l < layer_config->filter_height; ++l) {
                      const int ii = hh + l - filter_height_half;
                      for (int m = 0; m < layer_config->filter_width;
                           ++m, off += cstep) {
                        const int jj = ww + m - filter_width_half;
                        if (ii < 0 || ii >= in_height || jj < 0 ||
                            jj >= in_width)
                          continue;
                        sum += layer_config->weights[off] *
                               input[k][ii * in_stride + jj];
                      }
                    }
                  }
                  const float a = sum;
                  if (h == hh && w == ww)
                    output[i][u * out_stride + v] = a;
                  else
                    output[i][u * out_stride + v] =
                        AOMMAX(output[i][u * out_stride + v], a);
                }
              }
            }
          }
        }
        break;
      case PADDING_SAME_REPLICATE:
        for (int i = 0; i < layer_config->out_channels; ++i) {
          for (int h = 0, u = 0; h < in_height;
               h += layer_config->skip_height, ++u) {
            for (int w = 0, v = 0; w < in_width;
                 w += layer_config->skip_width, ++v) {
              for (int hh = h;
                   hh < AOMMIN(in_height, h + layer_config->skip_height);
                   ++hh) {
                for (int ww = w;
                     ww < AOMMIN(in_width, w + layer_config->skip_width);
                     ++ww) {
                  float sum = layer_config->bias[i];
                  for (int k = 0; k < layer_config->in_channels; ++k) {
                    int off = k * layer_config->out_channels + i;
                    for (int l = 0; l < layer_config->filter_height; ++l) {
                      const int ii =
                          CLAMPINDEX(hh + l - filter_height_half, in_height);
                      for (int m = 0; m < layer_config->filter_width;
                           ++m, off += cstep) {
                        const int jj =
                            CLAMPINDEX(ww + m - filter_width_half, in_width);
                        assert(ii >= 0 && ii < in_height && jj >= 0 &&
                               jj < in_width);
                        sum += layer_config->weights[off] *
                               input[k][ii * in_stride + jj];
                      }
                    }
                  }
                  const float a = sum;
                  if (h == hh && w == ww)
                    output[i][u * out_stride + v] = a;
                  else
                    output[i][u * out_stride + v] =
                        AOMMAX(output[i][u * out_stride + v], a);
                }
              }
            }
          }
        }
        break;
      case PADDING_VALID:
        for (int i = 0; i < layer_config->out_channels; ++i) {
          for (int h = 0, u = 0;
               h < in_height - layer_config->filter_height + 1;
               h += layer_config->skip_height, ++u) {
            for (int w = 0, v = 0;
                 w < in_width - layer_config->filter_width + 1;
                 w += layer_config->skip_width, ++v) {
              for (int hh = h;
                   hh < AOMMIN(in_height, h + layer_config->skip_height);
                   ++hh) {
                for (int ww = w;
                     ww < AOMMIN(in_width, w + layer_config->skip_width);
                     ++ww) {
                  float sum = layer_config->bias[i];
                  for (int k = 0; k < layer_config->in_channels; ++k) {
                    int off = k * layer_config->out_channels + i;
                    for (int l = 0; l < layer_config->filter_height; ++l) {
                      const int ii = hh + l;
                      for (int m = 0; m < layer_config->filter_width;
                           ++m, off += cstep) {
                        const int jj = ww + m;
                        assert(ii >= 0 && ii < in_height && jj >= 0 &&
                               jj < in_width);
                        sum += layer_config->weights[off] *
                               input[k][ii * in_stride + jj];
                      }
                    }
                  }
                  const float a = sum;
                  if (h == hh && w == ww)
                    output[i][u * out_stride + v] = a;
                  else
                    output[i][u * out_stride + v] =
                        AOMMAX(output[i][u * out_stride + v], a);
                }
              }
            }
          }
        }
        break;
      default: assert(0 && "Unknown padding type");
    }
  } else {
    // Results in element-wise matrix multiplication.
    if (layer_config->filter_height == 1 && layer_config->filter_width == 1) {
      const int start_h = get_start_shift_convolve(
          in_height, layer_config->filter_height, layer_config->skip_height);
      const int start_w =
          get_start_shift_convolve(in_width, layer_config->filter_width,
                                   layer_config->skip_width) +
          start_idx * layer_config->skip_width;
      const int out_w_step = AOMMAX(step, 1);
      const int in_w_step = layer_config->skip_width * out_w_step;
      for (int i = 0; i < layer_config->out_channels; ++i) {
        for (int h = start_h, u = 0; h < in_height;
             h += layer_config->skip_height, ++u) {
          const int in_h = h * in_stride;
          const int out_h = u * out_stride + start_idx;
          for (int w = start_w, out_index = out_h; w < in_width;
               w += in_w_step, out_index += out_w_step) {
            float sum = layer_config->bias[i];
            for (int k = 0; k < layer_config->in_channels; ++k) {
              sum += layer_config->weights[k * layer_config->out_channels + i] *
                     input[k][in_h + w];
            }
            output[i][out_index] = sum;
          }
        }
      }
      return;
    }
    const int ii_shift =
        filter_height_half - (layer_config->filter_height - 1) % 2;
    const int jj_shift =
        filter_width_half - (layer_config->filter_width - 1) % 2;
    switch (layer_config->pad) {
      case PADDING_SAME_ZERO: {
        const int start_h = get_start_shift_convolve(
            in_height, layer_config->filter_height, layer_config->skip_height);
        const int start_w = get_start_shift_convolve(
            in_width, layer_config->filter_width, layer_config->skip_width);
        const int end_ii_shift = filter_height_half + 1;
        const int end_jj_shift = filter_width_half + 1;
        // *_filter_margin stores the number of pixels along a dimension in the
        // intersection of the complement of the image in the extended image
        // and the filter.
        const int top_filter_margin = layer_config->filter_width * ii_shift;
        const int right_filter_margin = end_jj_shift - in_width;
        for (int i = start_idx; i < layer_config->out_channels;
             i += channel_step) {
          for (int h = start_h, u = 0; h < in_height;
               h += layer_config->skip_height, ++u) {
            const int out_h = u * out_stride;
            const int top_cstep =
                AOMMAX(0, top_filter_margin - h * layer_config->filter_width) *
                    cstep +
                i;
            const int start_ii = AOMMAX(0, h - ii_shift);
            const int end_ii = AOMMIN(in_height, h + end_ii_shift);
            for (int w = start_w, out_index = out_h; w < in_width;
                 w += layer_config->skip_width, ++out_index) {
              const int left_cstep = AOMMAX(0, jj_shift - w) * cstep;
              const int right_cstep =
                  AOMMAX(0, right_filter_margin + w) * cstep;
              const int start_jj = AOMMAX(0, w - jj_shift);
              const int end_jj = AOMMIN(in_width, w + end_jj_shift);
              float sum = layer_config->bias[i];
              for (int k = 0; k < layer_config->in_channels; ++k) {
                int off = k * layer_config->out_channels + top_cstep;
                for (int ii = start_ii; ii < end_ii; ++ii) {
                  off += left_cstep;
                  for (int jj = start_jj; jj < end_jj; ++jj, off += cstep) {
                    sum += layer_config->weights[off] *
                           input[k][ii * in_stride + jj];
                  }
                  off += right_cstep;
                }
              }
              output[i][out_index] = sum;
            }
          }
        }
        break;
      }
      case PADDING_SAME_REPLICATE: {
        // h and w are shifted to an offset coordinate system to reduce in-loop
        // computation.
        const int start_h =
            get_start_shift_convolve(in_height, layer_config->filter_height,
                                     layer_config->skip_height) -
            ii_shift;
        const int start_w =
            get_start_shift_convolve(in_width, layer_config->filter_width,
                                     layer_config->skip_width) -
            jj_shift;
        const int end_h = in_height - ii_shift;
        const int end_w = in_width - jj_shift;
        for (int i = start_idx; i < layer_config->out_channels;
             i += channel_step) {
          for (int h = start_h, u = 0; h < end_h;
               h += layer_config->skip_height, ++u) {
            const int out_h = u * out_stride;
            const int upper_ii_index = layer_config->filter_height + h;
            for (int w = start_w, out_index = out_h; w < end_w;
                 w += layer_config->skip_width, ++out_index) {
              const int upper_jj_index = layer_config->filter_width + w;
              float sum = layer_config->bias[i];
              for (int k = 0; k < layer_config->in_channels; ++k) {
                int off = k * layer_config->out_channels + i;
                for (int ii = h; ii < upper_ii_index; ++ii) {
                  const int clamped_ii = CLAMPINDEX(ii, in_height);
                  for (int jj = w; jj < upper_jj_index; ++jj) {
                    const int clamped_jj = CLAMPINDEX(jj, in_width);
                    assert(clamped_ii >= 0 && clamped_ii < in_height &&
                           clamped_jj >= 0 && clamped_jj < in_width);
                    sum += layer_config->weights[off] *
                           input[k][clamped_ii * in_stride + clamped_jj];
                    off += cstep;
                  }
                }
              }
              output[i][out_index] = sum;
            }
          }
        }
        break;
      }
      case PADDING_VALID: {
        if (layer_config->filter_width == 5 &&
            layer_config->filter_height == 5) {
          cnn_convlove_5x5_avx2(input, in_width, in_height, in_stride,
                                layer_config, output, out_stride, start_idx,
                                channel_step, cstep);
        } else {
          for (int i = start_idx; i < layer_config->out_channels;
               i += channel_step) {
            float out[16][16];
            for (int h = 0; h < 16; ++h) {
              for (int w = 0; w < 16; ++w) {
                out[h][w] = layer_config->bias[i];
              }
            }
            for (int k = 0; k < layer_config->in_channels; ++k) {
              float weight[5][5];
              int off = k * layer_config->out_channels + i;
              for (int weight_h = 0; weight_h < layer_config->filter_height;
                   weight_h++) {
                for (int weight_w = 0; weight_w < layer_config->filter_width;
                     weight_w++) {
                  weight[weight_h][weight_w] = layer_config->weights[off];
                  off += cstep;
                }
              }

              for (int h = 0, u = 0;
                   h < in_height - layer_config->filter_height + 1;
                   h += layer_config->skip_height, ++u) {
                const int upper_ii_index = layer_config->filter_height + h;
                for (int w = 0, v = 0;
                     w < in_width - layer_config->filter_width + 1;
                     w += layer_config->skip_width, ++v) {
                  const int upper_jj_index = layer_config->filter_width + w;
                  // float sum = layer_config->bias[i];
                  for (int ii = h, i_weight = 0; ii < upper_ii_index;
                       ++ii, ++i_weight) {
                    for (int jj = w, j_weight = 0; jj < upper_jj_index;
                         ++jj, ++j_weight) {
                      assert(ii >= 0 && ii < in_height && jj >= 0 &&
                             jj < in_width);
                      out[u][v] += weight[i_weight][j_weight] *
                                   input[k][ii * in_stride + jj];
                    }
                  }
                }
              }
            }

            for (int h = 0, u = 0;
                 h < in_height - layer_config->filter_height + 1;
                 h += layer_config->skip_height, ++u) {
              const int out_h = u * out_stride;
              for (int w = 0, v = 0, out_index = out_h;
                   w < in_width - layer_config->filter_width + 1;
                   w += layer_config->skip_width, ++v, ++out_index) {
                output[i][out_index] = out[u][v];
              }
            }
          }
        }
        break;
      }
      default: assert(0 && "Unknown padding type");
    }
  }
}
