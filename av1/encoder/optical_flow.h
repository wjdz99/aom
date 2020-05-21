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

#ifndef AOM_AV1_ENCODER_OPTICAL_FLOW_H_
#define AOM_AV1_ENCODER_OPTICAL_FLOW_H_

#ifdef __cplusplus
extern "C" {
#endif
typedef enum { LucasKanade = 0 } optflow_method;
typedef enum { None = 0, Smooth = 1, Median = 2 } mv_filter_type;
// default options for optical flow
#define OPFL_WINDOW_SIZE 25
#define OPFL_PYRAMID_LEVELS 3  // total levels (max is 5)

void optical_flow(const YV12_BUFFER_CONFIG *from_frame,
                  const YV12_BUFFER_CONFIG *to_frame, const int from_frame_idx,
                  const int to_frame_idx, const int bit_depth, int levels,
                  int window_size, const mv_filter_type mv_filter,
                  const optflow_method method, MV *mvs);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_OPTICAL_FLOW_H_
