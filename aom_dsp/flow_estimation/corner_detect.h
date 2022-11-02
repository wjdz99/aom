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

#ifndef AOM_FLOW_ESTIMATION_CORNER_DETECT_H_
#define AOM_FLOW_ESTIMATION_CORNER_DETECT_H_

#include "aom_scale/yv12config.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#ifdef __cplusplus
extern "C" {
#endif

void aom_find_corners_in_frame(YV12_BUFFER_CONFIG *frm, int bit_depth);

#ifdef __cplusplus
}
#endif

#endif  // AOM_FLOW_ESTIMATION_CORNER_DETECT_H_
