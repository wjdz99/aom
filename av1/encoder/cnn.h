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

#ifndef AOM_AV1_COMMON_CNN_H_
#define AOM_AV1_COMMON_CNN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

#include "aom_util/aom_thread.h"
#include "config/av1_rtcd.h"

struct AV1Common;

#define CNN_MAX_HIDDEN_LAYERS 64
#define CNN_MAX_LAYERS (CNN_MAX_HIDDEN_LAYERS + 1)
#define CNN_MAX_CHANNELS 256
#define CNN_MAX_BRANCHES 4
#define CNN_MAX_THREADS 32

#define NO_BRANCH_CONFIG \
  { 0, 0, 0 }
#define NO_BN_PARAMS \
  { NULL, NULL, NULL, NULL }

#define CLAMPINDEX(a, hi) ((a) < 0 ? 0 : ((a) >= (hi) ? ((hi)-1) : (a)))

enum {
  PADDING_SAME_ZERO,       // tensorflow's SAME padding with pixels outside
                           // the image area assumed to be 0 (default)
  PADDING_SAME_REPLICATE,  // tensorflow's SAME padding with pixels outside
                           // the image area replicated from closest edge
  PADDING_VALID            // tensorflow's VALID padding
} UENUM1BYTE(PADDING_TYPE);

// enum { NONE, RELU, SOFTSIGN } UENUM1BYTE(ACTIVATION);

// Times when input tensor may be copied to branches given in input_to_branches.
// BRANCH_NO_COPY: doesn't copy any tensor.
// BRANCH_INPUT: copies the input tensor to branches.
// BRANCH_OUTPUT: copies the convolved tensor to branches.
// BRANCH_COMBINED: copies the combined (after convolving and branch combining)
//   tensor. If no combinations happen at this layer, then this option
//   has the same effect as COPY_OUTPUT.
enum {
  BRANCH_NO_COPY,
  BRANCH_INPUT,
  BRANCH_OUTPUT,
  BRANCH_COMBINED
} UENUM1BYTE(BRANCH_COPY);

// Types of combining branches with output of current layer:
// BRANCH_NOC: no branch combining
// BRANCH_ADD: Add previously stored branch tensor to output of layer
// BRANCH_CAT: Concatenate branch tensor to output of layer
enum { BRANCH_NOC, BRANCH_ADD, BRANCH_CAT } UENUM1BYTE(BRANCH_COMBINE);

// The parameters used to scale each channel in batch
// normalization. The processing in done on a per-channel basis.
// e.g. bn_mean[c] is the mean for all pixels in channel c. This
// is always applied after activation. The output is given by
// out[c,i,j] = norm[c,i,j] * bn_gamma[c] + bn_beta[c] where
// norm[c,i,j] = (in[c,i,j] - bn_mean[c]) / bn_std[c]
// here we assume that the effect of variance_epsilon is already
// taken into account when bn_std is calculated. The pointers
// needs to be either all zero or all valid. If all zero, then
// batchnorm is disabled, else batchnorm is applied.
struct CNN_BATCHNORM_PARAMS {
  const float *bn_gamma;
  const float *bn_beta;
  const float *bn_mean;
  const float *bn_std;
};

struct CNN_BRANCH_CONFIG {
  int input_to_branches;  // If nonzero, copy the active tensor to the current
  // layer and store for future use in branches
  // specified in the field as a binary mask. For
  // example, if input_to_branch = 0x06, it means the
  // input tensor to the current branch is copied to
  // branches 1 and 2 (where 0 represents the primary
  // branch). One restriction is that the mask
  // cannot indicate copying to the current branch.
  // If greater than 0, only copies the channels up
  // to the given index.
  int channels_to_copy;  // Within the layer, input a copy of active
  // tensor to branches given in input_to_branches.
  int branches_to_combine;  // mask of branches to combine with output of
  // current layer, if
  // branch_combine_type != BRANCH_NOC
  // For example, if branches_to_combine = 0x0A,
  // it means that braches 1 and 3 are combined
  // with the current branch.
};

struct CNN_LAYER_CONFIG {
  int in_channels;
  int filter_width;
  int filter_height;
  int out_channels;
  int skip_width;
  int skip_height;
  int maxpool;            // whether to use maxpool or not (only effective when
                          // skip width or skip_height are > 1)
  const float *weights;   // array of length filter_height x filter_width x
                          // in_channels x out_channels where the inner-most
                          // scan is out_channels and the outer most scan is
                          // filter_height.
  const float *bias;      // array of length out_channels
  PADDING_TYPE pad;       // padding type
  ACTIVATION activation;  // the activation function to use after convolution
  int deconvolve;         // whether this is a deconvolution layer.
                          // 0: If skip_width or skip_height are > 1, then we
                          // reduce resolution
                          // 1: If skip_width or skip_height are > 1, then we
                          // increase resolution
  int branch;             // branch index in [0, CNN_MAX_BRANCHES - 1], where
                          // 0 refers to the primary branch.
  BRANCH_COPY branch_copy_type;
  BRANCH_COMBINE branch_combine_type;
  struct CNN_BRANCH_CONFIG branch_config;
  struct CNN_BATCHNORM_PARAMS
      bn_params;   // A struct that contains the parameters
                   // used for batch normalization.
  int output_num;  // The output buffer idx to which the layer output is
                   // written. Set to -1 to disable writing it to the output. In
                   // the case that branch_combine_type is BRANCH_CAT, all
                   // concatenated channels will be written to output. In the
                   // case of BRANCH_ADD, the output will be the result of
                   // summation.
};

struct CNN_CONFIG {
  int num_layers;  // number of CNN layers ( = number of hidden layers + 1)
  int is_residue;  // whether the output activation is a residue
  int ext_width, ext_height;  // extension horizontally and vertically
  int strict_bounds;          // whether the input bounds are strict or not.
                              // If strict, the extension area is filled by
                              // replication; if not strict, image data is
                              // assumed available beyond the bounds.
  CNN_LAYER_CONFIG layer_config[CNN_MAX_LAYERS];
};

struct CNN_THREAD_DATA {
  int num_workers;
  AVxWorker *workers;
};

struct CNN_MULTI_OUT {
  int num_outputs;
  const int *output_channels;
  const int *output_strides;
  float **output_buffer;
};

// Function to return size of output
void av1_find_cnn_output_size(int in_width, int in_height,
                              const CNN_CONFIG *cnn_config, int *out_width,
                              int *out_height, int *out_channels);

static INLINE void find_layer_output_size(int in_width, int in_height,
                                          const CNN_LAYER_CONFIG *layer_config,
                                          int *out_width, int *out_height) {
  if (!layer_config->deconvolve) {
    switch (layer_config->pad) {
      case PADDING_SAME_ZERO:
      case PADDING_SAME_REPLICATE:
        *out_width = (in_width + layer_config->skip_width - 1) /
                     layer_config->skip_width;
        *out_height = (in_height + layer_config->skip_height - 1) /
                      layer_config->skip_height;
        break;
      case PADDING_VALID:
        *out_width =
            (in_width - layer_config->filter_width + layer_config->skip_width) /
            layer_config->skip_width;
        *out_height = (in_height - layer_config->filter_height +
                       layer_config->skip_height) /
                      layer_config->skip_height;
        break;
      default: assert(0 && "Unknown padding type");
    }
  } else {
    switch (layer_config->pad) {
      case PADDING_SAME_ZERO:
      case PADDING_SAME_REPLICATE:
        *out_width = in_width * layer_config->skip_width;
        *out_height = in_height * layer_config->skip_height;
        break;
      case PADDING_VALID:
        *out_width = (in_width - 1) * layer_config->skip_width +
                     layer_config->filter_width;
        *out_height = (in_height - 1) * layer_config->skip_height +
                      layer_config->filter_height;
        break;
      default: assert(0 && "Unknown padding type");
    }
  }
}

// Prediction functions from set of input image buffers. This function supports
// CNN with multiple outputs.
void av1_cnn_predict_img_multi_out(uint8_t **dgd, int width, int height,
                                   int stride, const CNN_CONFIG *cnn_config,
                                   const CNN_THREAD_DATA *thread_data,
                                   struct CNN_MULTI_OUT *output);
void av1_cnn_predict_img_multi_out_highbd(uint16_t **dgd, int width, int height,
                                          int stride,
                                          const CNN_CONFIG *cnn_config,
                                          const CNN_THREAD_DATA *thread_data,
                                          int bit_depth, CNN_MULTI_OUT *output);

// Prediction functions from set of input image buffers. This function only
// supports a single output.
void av1_cnn_predict_img(uint8_t **dgd, int width, int height, int stride,
                         const CNN_CONFIG *cnn_config,
                         const CNN_THREAD_DATA *thread_data, float **output,
                         int out_stride);
void av1_cnn_predict_img_highbd(uint16_t **dgd, int width, int height,
                                int stride, const CNN_CONFIG *cnn_config,
                                const CNN_THREAD_DATA *thread_data,
                                int bit_depth, float **output, int out_stride);

// CNNConvolve get start shift value.
static INLINE int get_start_shift_convolve(int width, int filt_width,
                                           int stride) {
  const int mod = (width % stride);
  const int filt_off = (filt_width - 1) / 2;
  const int dif = (mod ? mod - 1 : stride - 1);
  return AOMMIN((dif + (filt_width % 2)) / 2, filt_off);
}

// CNNConvolve specific to maxpool set as 1 and padding equal to
// PADDING_SAME_ZERO.
static INLINE void maxpool_padding_zero(
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

// CNNConvolve specific to maxpool set as 1 and padding equal to
// PADDING_SAME_REPLICATE.
static INLINE void maxpool_padding_replicate(
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

// CNNConvolve specific to maxpool set as 1 and padding equal to
// PADDING_VALID.
static INLINE void maxpool_padding_valid(
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

// CNNConvolve specific to maxpool set as 0 with filter_height & filter_width
// equals 1.
static INLINE void element_wise_convolve(
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
static INLINE void no_maxpool_padding_zero(
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
static INLINE void no_maxpool_padding_replicate(
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
static INLINE void no_maxpool_padding_valid(
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

#endif  // AOM_AV1_COMMON_CNN_H_
