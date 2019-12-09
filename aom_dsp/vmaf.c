/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <stdlib.h>
#include "vmaf.h"
#include "aom_mem/aom_mem.h"

#if CONFIG_TUNE_VMAF

typedef struct framedata {
  const YV12_BUFFER_CONFIG *source;
  const YV12_BUFFER_CONFIG *distorted;
  int frame_set;
} framedata;

int read_frame_8bd(float *ref_data, float *main_data, float *temp_data,
                   int stride, void *user_data) {
  framedata *frames = (framedata *)user_data;

  if (!frames->frame_set) {
    const int width = frames->source->y_width;
    const int height = frames->source->y_height;
    assert(width == frames->distorted->y_width);
    assert(height == frames->distorted->y_height);
    uint8_t *ref_ptr = frames->source->y_buffer;
    uint8_t *main_ptr = frames->distorted->y_buffer;

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        ref_data[col] = (float)ref_ptr[col];
      }
      ref_ptr += frames->source->y_stride;
      ref_data += stride / sizeof(*ref_data);
    }

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        main_data[col] = (float)main_ptr[col];
      }
      main_ptr += frames->distorted->y_stride;
      main_data += stride / sizeof(*main_data);
    }
    frames->frame_set = 1;
    return 0;
  }

  (void)temp_data;
  return 2;
}

double aom_calc_vmaf(const char *model_path, const YV12_BUFFER_CONFIG *source,
                     const YV12_BUFFER_CONFIG *distorted) {
  const int width = source->y_width;
  const int height = source->y_height;
  struct framedata frames = { source, distorted, 0 };
  double vmaf_score;
  int (*read_frame)(float *reference_data, float *distorted_data,
                    float *temp_data, int stride, void *s);
  read_frame = read_frame_8bd;
  int ret = compute_vmaf(
      &vmaf_score, (char *)"yuv420p", width, height, read_frame,
      /* user_data */ &frames, (char *)model_path,
      /* log_path */ NULL, /* log_fmt */ NULL, /* disable_clip */ 0,
      /* disable_avx */ 0, /* enable_transform */ 0,
      /* phone_model */ 0, /* do_psnr */ 0, /* do_ssim */ 0,
      /* do_ms_ssim */ 0, /* pool_method */ NULL, /* thread */ 1,
      /* subsample */ 1, /* enable_conf_interval */ 0);

  if (ret != 0) exit(EXIT_FAILURE);
  return vmaf_score;
}

#endif  // CONFIG_TUNE_VMAF
