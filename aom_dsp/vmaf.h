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

#include <stdbool.h>
#include "aom_scale/yv12config.h"

typedef struct VmafModel VmafModel;

void aom_init_vmaf_model(VmafModel **vmaf_model, const char *model_path);
void aom_close_vmaf_model(VmafModel *vmaf_model);

void aom_calc_vmaf(VmafModel *vmaf_model, const YV12_BUFFER_CONFIG *source,
                   const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                   bool cal_vmaf_neg, double *vmaf);

#endif  // AOM_AOM_DSP_VMAF_H_
