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
#if CONFIG_OPTICAL_FLOW_API

#include <math.h>
#include <limits.h>

#include "config/aom_config.h"
#include "av1/common/av1_common_int.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/mathutils.h"
#include "av1/encoder/optical_flow.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/temporal_filter.h"
#include "aom_mem/aom_mem.h"

typedef struct LOCALMV {
  double row;
  double col;
} LOCALMV;

// Computes optical flow by applying algorithm at
// multiple pyramid levels of images (lower-resolution, smoothed images)
// This accounts for larger motions.
// Inputs:
//   from_frame Frame buffer.
//   to_frame: Frame buffer. MVs point from_frame -> to_frame.
//   from_frame_idx: Index of from_frame.
//   to_frame_idx: Index of to_frame. Return all zero MVs when idx are equal.
//   bit_depth:
//   levels: total pyramid levels. Must be in [1,5].
//   window_size: as used by specific algorithm.
//   mv_filter: None, Smooth, or Median.
//   method: LucasKanade,
//   mvs: pointer to MVs. Contains initialization, and modified
//   based on optical flow.
void optical_flow(const YV12_BUFFER_CONFIG *from_frame,
                  const YV12_BUFFER_CONFIG *to_frame, const int from_frame_idx,
                  const int to_frame_idx, const int bit_depth,
                  OPFL_PARAMS opfl_params, const mv_filter_type mv_filter,
                  const optflow_method method, MV *mvs) {
  int levels = opfl_params.pyramid_levels;
  if (levels < 1 || levels > 5) {
    printf("Pyramid levels out of bounds. Choose a value within [%d, %d].\n", 1,
           5);
    printf("Resetting to default value %d\n", OPFL_PYRAMID_LEVELS);
    levels = OPFL_PYRAMID_LEVELS;
  }
  const int frame_height = from_frame->y_crop_height;
  const int frame_width = from_frame->y_crop_width;
  LOCALMV *localmvs = malloc(frame_height * frame_width * sizeof(LOCALMV));
  if (from_frame_idx == to_frame_idx) {
    // immediately return all zero mvs when frame indices are equal
    for (int yy = 0; yy < frame_height; yy++) {
      for (int xx = 0; xx < frame_width; xx++) {
        MV mv = { .row = 0, .col = 0 };
        mvs[yy * frame_width + xx] = mv;
      }
    }
    return;
  }
  // Initialize double mvs based on input parameter mvs array
  for (int i = 0; i < frame_width * frame_height; i++) {
    MV mv = mvs[i];
    LOCALMV localmv = { .row = mv.row / 8.0, .col = mv.col / 8.0 };
    localmvs[i] = localmv;
  }
  // Apply optical flow algorithm

  // Update original mvs array
  for (int j = 0; j < frame_height; j++) {
    for (int i = 0; i < frame_width; i++) {
      int idx = j * frame_width + i;
      int new_x = localmvs[idx].row + i;
      int new_y = localmvs[idx].col + j;
      if ((fabs(localmvs[idx].row) >= 0.125 ||
           fabs(localmvs[idx].col) >= 0.125)) {
        // if mv points outside of frame (lost feature), keep old mv.
        if (new_x < frame_width && new_x >= 0 && new_y < frame_height &&
            new_y >= 0) {
          MV mv = { .row = round(8 * localmvs[idx].row),
                    .col = round(8 * localmvs[idx].col) };
          mvs[idx] = mv;
        }
      }
    }
  }
  free(localmvs);
}
#endif
