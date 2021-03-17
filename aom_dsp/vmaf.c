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

#include "aom_dsp/vmaf.h"

#include <assert.h>
#include <libvmaf/libvmaf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

#include "aom_dsp/blend.h"
#include "aom_ports/system_state.h"

static void vmaf_fatal_error(const char *message) {
  fprintf(stderr, "Fatal error: %s\n", message);
  exit(EXIT_FAILURE);
}

void aom_init_vmaf_model(VmafModel **vmaf_model, const char *model_path) {
  if (*vmaf_model != NULL) return;
  VmafModelConfig model_cfg;
  model_cfg.flags = VMAF_MODEL_FLAG_DISABLE_CLIP;
  model_cfg.name = "vmaf";

  if (vmaf_model_load_from_path(vmaf_model, &model_cfg, model_path)) {
    vmaf_fatal_error("Failed to load VMAF model.");
  }
}

void aom_close_vmaf_model(VmafModel *vmaf_model) {
  vmaf_model_destroy(vmaf_model);
}

static void copy_picture(const int bit_depth, const YV12_BUFFER_CONFIG *src,
                         VmafPicture *dst) {
  const int width = src->y_width;
  const int height = src->y_height;

  if (bit_depth > 8) {
    uint16_t *src_ptr = CONVERT_TO_SHORTPTR(src->y_buffer);
    uint16_t *dst_ptr = dst->data[0];

    for (int row = 0; row < height; ++row) {
      memcpy(dst_ptr, src_ptr, width * sizeof(dst_ptr[0]));
      src_ptr += src->y_stride;
      dst_ptr += dst->stride[0] / 2;
    }
  } else {
    uint8_t *src_ptr = src->y_buffer;
    uint8_t *dst_ptr = (uint8_t *)dst->data[0];

    for (int row = 0; row < height; ++row) {
      memcpy(dst_ptr, src_ptr, width * sizeof(dst_ptr[0]));
      src_ptr += src->y_stride;
      dst_ptr += dst->stride[0];
    }
  }
}

void aom_init_vmaf_context(VmafContext **vmaf_context, VmafModel *vmaf_model,
                           bool cal_vmaf_neg) {
  VmafConfiguration cfg;
  cfg.log_level = VMAF_LOG_LEVEL_NONE;
  cfg.n_threads = 0;
  cfg.n_subsample = 0;
  cfg.cpumask = 0;

  if (vmaf_init(vmaf_context, cfg)) {
    vmaf_fatal_error("Failed to init VMAF context.");
  }

  if (vmaf_use_features_from_model(*vmaf_context, vmaf_model)) {
    vmaf_fatal_error("Failed to load feature extractors from VMAF model.");
  }

  if (cal_vmaf_neg) {
    VmafFeatureDictionary *vif_feature = NULL;
    vmaf_feature_dictionary_set(&vif_feature, "vif_enhn_gain_limit", "1.0");
    if (vmaf_use_feature(*vmaf_context, "float_vif", vif_feature)) {
      vmaf_fatal_error("Failed to use feature float_vif.");
    }

    VmafFeatureDictionary *adm_feature = NULL;
    vmaf_feature_dictionary_set(&adm_feature, "adm_enhn_gain_limit", "1.0");
    if (vmaf_use_feature(*vmaf_context, "float_adm", adm_feature)) {
      vmaf_fatal_error("Failed to use feature float_adm.");
    }
  }

  VmafFeatureDictionary *motion_force_zero = NULL;
  vmaf_feature_dictionary_set(&motion_force_zero, "motion_force_zero", "true");
  if (vmaf_use_feature(*vmaf_context, "float_motion", motion_force_zero)) {
    vmaf_fatal_error("Failed to use feature float_motion.");
  }
}

void aom_close_vmaf_context(VmafContext *vmaf_context) {
  if (vmaf_close(vmaf_context)) {
    vmaf_fatal_error("Failed to close VMAF context.");
  }
}

void aom_calc_vmaf_at_index(VmafContext *vmaf_context, VmafModel *vmaf_model,
                            const YV12_BUFFER_CONFIG *source,
                            const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                            int frame_index, double *vmaf) {
  VmafPicture ref, dist;
  if (vmaf_picture_alloc(&ref, VMAF_PIX_FMT_YUV420P, bit_depth, source->y_width,
                         source->y_height) ||
      vmaf_picture_alloc(&dist, VMAF_PIX_FMT_YUV420P, bit_depth,
                         source->y_width, source->y_height)) {
    vmaf_fatal_error("Failed to alloc VMAF pictures.");
  }
  copy_picture(bit_depth, source, &ref);
  copy_picture(bit_depth, distorted, &dist);
  if (vmaf_read_pictures(vmaf_context, &ref, &dist,
                         /*picture index=*/frame_index)) {
    vmaf_fatal_error("Failed to read VMAF pictures.");
  }

  vmaf_picture_unref(&ref);
  vmaf_picture_unref(&dist);

  vmaf_score_at_index(vmaf_context, vmaf_model, vmaf, frame_index);
}
