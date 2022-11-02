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

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "third_party/fastfeat/fast.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/flow_estimation/corner_detect.h"
#include "aom_dsp/flow_estimation/flow_estimation.h"
#include "aom_dsp/flow_estimation/pyramid.h"
#include "aom_mem/aom_mem.h"

// Fast_9 wrapper
#define FAST_BARRIER 18
void aom_find_corners_in_frame(YV12_BUFFER_CONFIG *frm, int bit_depth) {
  if (frm->corners) {
    // Already computed, no need to do it again
    return;
  }

  ImagePyramid *pyr = aom_compute_pyramid(frm, bit_depth, MAX_PYRAMID_LEVELS);
  PyramidLayer *layer = &pyr->layers[0];

  int num_points;
  xy *const frm_corners_xy =
      aom_fast9_detect_nonmax(layer->buffer, layer->width, layer->height,
                              layer->stride, FAST_BARRIER, &num_points);

  num_points = AOMMIN(num_points, MAX_CORNERS);

  if (num_points > 0 && frm_corners_xy) {
    frm->corners = aom_malloc(2 * num_points * sizeof(*frm->corners));
    memcpy(frm->corners, frm_corners_xy, sizeof(*frm_corners_xy) * num_points);
    frm->num_corners = num_points;
  } else {
    frm->corners = NULL;
    frm->num_corners = 0;
  }

  free(frm_corners_xy);
}
