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

#ifndef AOM_AV1_ENCODER_CNN_INTERNAL_H_
#define AOM_AV1_ENCODER_CNN_INTERNAL_H_

#ifdef __cplusplus
extern "C" {
#endif

// CNNConvolve get start shift value.
static INLINE int get_start_shift_convolve(int width, int filt_width,
                                           int stride) {
  const int mod = (width % stride);
  const int filt_off = (filt_width - 1) / 2;
  const int dif = (mod ? mod - 1 : stride - 1);
  return AOMMIN((dif + (filt_width % 2)) / 2, filt_off);
}

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_SAME_ZERO.
static INLINE void cnn_maxpool_padding_zero(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep, const int filter_width_half,
    const int filter_height_half) {
  for (int i = 0; i < layer_config->out_channels; ++i) {
    for (int h = 0, u = 0; h < in_height; h += layer_config->skip_height, ++u) {
      for (int w = 0, v = 0; w < in_width; w += layer_config->skip_width, ++v) {
        for (int hh = h; hh < AOMMIN(in_height, h + layer_config->skip_height);
             ++hh) {
          for (int ww = w; ww < AOMMIN(in_width, w + layer_config->skip_width);
               ++ww) {
            float sum = layer_config->bias[i];
            for (int k = 0; k < layer_config->in_channels; ++k) {
              int off = k * layer_config->out_channels + i;
              for (int l = 0; l < layer_config->filter_height; ++l) {
                const int ii = hh + l - filter_height_half;
                for (int m = 0; m < layer_config->filter_width;
                     ++m, off += cstep) {
                  const int jj = ww + m - filter_width_half;
                  if (ii < 0 || ii >= in_height || jj < 0 || jj >= in_width)
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
}

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_SAME_REPLICATE.
static INLINE void cnn_maxpool_padding_replicate(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep, const int filter_width_half,
    const int filter_height_half) {
  for (int i = 0; i < layer_config->out_channels; ++i) {
    for (int h = 0, u = 0; h < in_height; h += layer_config->skip_height, ++u) {
      for (int w = 0, v = 0; w < in_width; w += layer_config->skip_width, ++v) {
        for (int hh = h; hh < AOMMIN(in_height, h + layer_config->skip_height);
             ++hh) {
          for (int ww = w; ww < AOMMIN(in_width, w + layer_config->skip_width);
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
                  assert(ii >= 0 && ii < in_height && jj >= 0 && jj < in_width);
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
}

// CNNConvolve specific to maxpool set as 1, either skip_width or skip_height
// greater than 1 and padding equal to PADDING_VALID.
static INLINE void cnn_maxpool_padding_valid(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    const int cstep) {
  for (int i = 0; i < layer_config->out_channels; ++i) {
    for (int h = 0, u = 0; h < in_height - layer_config->filter_height + 1;
         h += layer_config->skip_height, ++u) {
      for (int w = 0, v = 0; w < in_width - layer_config->filter_width + 1;
           w += layer_config->skip_width, ++v) {
        for (int hh = h; hh < AOMMIN(in_height, h + layer_config->skip_height);
             ++hh) {
          for (int ww = w; ww < AOMMIN(in_width, w + layer_config->skip_width);
               ++ww) {
            float sum = layer_config->bias[i];
            for (int k = 0; k < layer_config->in_channels; ++k) {
              int off = k * layer_config->out_channels + i;
              for (int l = 0; l < layer_config->filter_height; ++l) {
                const int ii = hh + l;
                for (int m = 0; m < layer_config->filter_width;
                     ++m, off += cstep) {
                  const int jj = ww + m;
                  assert(ii >= 0 && ii < in_height && jj >= 0 && jj < in_width);
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
}

// CNNConvolve specific to maxpool set as 0 with filter_height and filter_width
// equal to 1.
static INLINE void cnn_element_wise_convolve(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, int step) {
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
}

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_SAME_ZERO.
static INLINE void cnn_no_maxpool_padding_zero(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int filter_width_half,
    const int filter_height_half, const int ii_shift, const int jj_shift,
    const int channel_step) {
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
  for (int i = start_idx; i < layer_config->out_channels; i += channel_step) {
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
        const int right_cstep = AOMMAX(0, right_filter_margin + w) * cstep;
        const int start_jj = AOMMAX(0, w - jj_shift);
        const int end_jj = AOMMIN(in_width, w + end_jj_shift);
        float sum = layer_config->bias[i];
        for (int k = 0; k < layer_config->in_channels; ++k) {
          int off = k * layer_config->out_channels + top_cstep;
          for (int ii = start_ii; ii < end_ii; ++ii) {
            off += left_cstep;
            for (int jj = start_jj; jj < end_jj; ++jj, off += cstep) {
              sum += layer_config->weights[off] * input[k][ii * in_stride + jj];
            }
            off += right_cstep;
          }
        }
        output[i][out_index] = sum;
      }
    }
  }
}

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_SAME_REPLICATE.
static INLINE void cnn_no_maxpool_padding_replicate(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int ii_shift, const int jj_shift,
    const int channel_step) {
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
  for (int i = start_idx; i < layer_config->out_channels; i += channel_step) {
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
}

// CNNConvolve specific to maxpool set as 0 and padding equal to
// PADDING_VALID.
static INLINE void cnn_no_maxpool_padding_valid(
    const float **input, int in_width, int in_height, int in_stride,
    const CNN_LAYER_CONFIG *const layer_config, float **output, int out_stride,
    int start_idx, const int cstep, const int channel_step) {
  for (int i = start_idx; i < layer_config->out_channels; i += channel_step) {
    for (int h = 0, u = 0; h < in_height - layer_config->filter_height + 1;
         h += layer_config->skip_height, ++u) {
      const int out_h = u * out_stride;
      const int upper_ii_index = layer_config->filter_height + h;
      for (int w = 0, out_index = out_h;
           w < in_width - layer_config->filter_width + 1;
           w += layer_config->skip_width, ++out_index) {
        const int upper_jj_index = layer_config->filter_width + w;
        float sum = layer_config->bias[i];
        for (int k = 0; k < layer_config->in_channels; ++k) {
          int off = k * layer_config->out_channels + i;
          for (int ii = h; ii < upper_ii_index; ++ii) {
            for (int jj = w; jj < upper_jj_index; ++jj) {
              assert(ii >= 0 && ii < in_height && jj >= 0 && jj < in_width);
              sum += layer_config->weights[off] * input[k][ii * in_stride + jj];
              off += cstep;
            }
          }
        }
        output[i][out_index] = sum;
      }
    }
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_CNN_INTERNAL_H_
