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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif
#include <libvmaf.rc.h>

#include "aom_dsp/blend.h"
#include "aom_dsp/vmaf_neg.h"
#include "aom_ports/system_state.h"

static void vmaf_fatal_error(const char *message) {
  fprintf(stderr, "Fatal error: %s\n", message);
  exit(EXIT_FAILURE);
}

void aom_init_vmaf_rc(VmafContext **vmaf_context, VmafModel **vmaf_model,
                      const char *model_path) {
  VmafConfiguration cfg;
  cfg.log_level = VMAF_LOG_LEVEL_INFO;
  cfg.n_threads = 0;
  cfg.n_subsample = 1;
  cfg.cpumask = 0;

  VmafModelConfig model_cfg;
  model_cfg.flags = VMAF_MODEL_FLAGS_DEFAULT;
  model_cfg.name = "vmaf";
  model_cfg.path = (char *)model_path;

  if (vmaf_init(vmaf_context, cfg)) {
    vmaf_fatal_error("Failed to init VMAF context.");
  }

  if (vmaf_model_load_from_path(vmaf_model, &model_cfg)) {
    vmaf_fatal_error("Failed to load VMAF model.");
  }

  if (vmaf_use_features_from_model(*vmaf_context, *vmaf_model)) {
    vmaf_fatal_error("Failed to load feature extractors from VMAF model.");
  }
}

void aom_close_vmaf_rc(VmafContext *vmaf_context, VmafModel *vmaf_model) {
  vmaf_model_destroy(vmaf_model);
  if (vmaf_close(vmaf_context)) {
    vmaf_fatal_error("Failed to close VMAF context.");
  }
}

void picture_copy(const int bit_depth, const YV12_BUFFER_CONFIG *src,
                  VmafPicture *dst) {
  const int width = src->y_width;
  const int height = src->y_height;
  float *dst_ptr = (float *)dst->data[0];

  if (bit_depth > 8) {
    const float scale_factor = 1.0f / (float)(1 << (bit_depth - 8));
    uint16_t *src_ptr = CONVERT_TO_SHORTPTR(src->y_buffer);

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        dst_ptr[col] = scale_factor * (float)src_ptr[col];
      }
      src_ptr += src->y_stride;
      dst_ptr += dst->stride[0];
    }
  } else {
    uint8_t *src_ptr = src->y_buffer;

    for (int row = 0; row < height; ++row) {
      for (int col = 0; col < width; ++col) {
        dst_ptr[col] = (float)src_ptr[col];
      }
      src_ptr += src->y_stride;
      dst_ptr += dst->stride[0];
    }
  }
}

void aom_calc_vmaf_rc(VmafContext *vmaf_context, VmafModel *vmaf_model,
                      const YV12_BUFFER_CONFIG *source,
                      const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                      double *vmaf) {
  aom_clear_system_state();
  VmafPicture ref, dist;
  printf("before vmaf alloc pictures\n");
  if (vmaf_picture_alloc(&ref, VMAF_PIX_FMT_YUV420P, /*bit depth=*/8,
                         source->y_width, source->y_height) ||
      vmaf_picture_alloc(&dist, VMAF_PIX_FMT_YUV420P, /*bit depth=*/8,
                         distorted->y_width, distorted->y_height)) {
    vmaf_fatal_error("Failed to alloc VMAF pictures.");
  }
  picture_copy(bit_depth, source, &ref);
  picture_copy(bit_depth, distorted, &dist);
  printf("after vmaf alloc pictures\n");
  if (vmaf_read_pictures(vmaf_context, &ref, &dist, /*picture index=*/0)) {
    vmaf_fatal_error("Failed to read VMAF pictures.");
  }
  printf("after vmaf read pictures\n");
  if (vmaf_score_at_index(vmaf_context, vmaf_model, vmaf,
                          /*picture index=*/0)) {
    vmaf_fatal_error("Failed to calculate VMAF score.");
  }
  aom_clear_system_state();
}
