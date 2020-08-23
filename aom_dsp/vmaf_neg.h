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

#ifndef AOM_AOM_DSP_VMAF_H_
#define AOM_AOM_DSP_VMAF_H_

#include "aom_scale/yv12config.h"

typedef struct VmafContext VmafContext;
typedef struct VmafModel VmafModel;

typedef struct {
  // Stores the scaling factors for rdmult when tuning for VMAF.
  // rdmult_scaling_factors[row * num_cols + col] stores the scaling factors for
  // 64x64 block at (row, col).
  double *rdmult_scaling_factors;

  // Stores the luma sse of the last frame.
  double last_frame_ysse;

  // Stores the VMAF of the last frame.
  double last_frame_vmaf;

  // Stores the filter strength of the last frame.
  double last_frame_unsharp_amount;

  double best_unsharp_amount;

  int original_qindex;

  // VMAF context used in VMAF_RC caculations.
  VmafContext *vmaf_context;

  // VMAF model used in VMAF_RC caculations.
  VmafModel *vmaf_model;
} TuneVMAFInfo;

void aom_init_vmaf_rc(VmafContext **vmaf_context, VmafModel **vmaf_model,
                      const char *model_path);

void aom_calc_vmaf_rc(VmafContext *vmaf_context, VmafModel *vmaf_model,
                      const YV12_BUFFER_CONFIG *source,
                      const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                      double *vmaf);

void aom_close_vmaf_rc(VmafContext *vmaf_context, VmafModel *vmaf_model);

#endif  // AOM_AOM_DSP_VMAF_H_
