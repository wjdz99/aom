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
#ifndef AV1_ENCODER_PICKRST_H_
#define AV1_ENCODER_PICKRST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/encoder.h"

struct yv12_buffer_config;
struct AV1_COMP;

int64_t sse_restoration_frame(AV1_COMMON *const cm,
                              const YV12_BUFFER_CONFIG *src,
                              const YV12_BUFFER_CONFIG *dst,
                              int components_pattern);
uint64_t av1_pick_filter_restoration(const YV12_BUFFER_CONFIG *sd,
                                     AV1_COMP *cpi, LPF_PICK_METHOD method);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_PICKRST_H_
