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
#include "aom_mem/aom_mem.h"

#define FAST_BARRIER 18

size_t av1_get_corner_list_size() { return sizeof(CornerList); }

CornerList *av1_alloc_corner_list() {
  CornerList *corners = (CornerList *)aom_calloc(1, sizeof(CornerList));
  if (!corners) {
    return NULL;
  }

  corners->valid = false;
  pthread_mutex_init(&corners->mutex, NULL);
  return corners;
}

void compute_corner_list(const ImagePyramid *pyr, CornerList *corners) {
  const uint8_t *buf = pyr->layers[0].buffer;
  size_t width = pyr->layers[0].width;
  size_t height = pyr->layers[0].height;
  size_t stride = pyr->layers[0].stride;

  int num_corners;
  xy *const frm_corners_xy = aom_fast9_detect_nonmax(
      buf, width, height, stride, FAST_BARRIER, &num_corners);
  num_corners = AOMMIN(num_corners, MAX_CORNERS);
  if (num_corners > 0 && frm_corners_xy) {
    memcpy(corners->corners, frm_corners_xy,
           sizeof(*frm_corners_xy) * num_corners);
    corners->num_corners = num_corners;
  } else {
    corners->num_corners = 0;
  }
  free(frm_corners_xy);
}

void av1_compute_corner_list(const ImagePyramid *pyr, CornerList *corners) {
  assert(corners);

#if CONFIG_MULTITHREAD
  pthread_mutex_lock(&corners->mutex);
#endif  // CONFIG_MULTITHREAD

  if (!corners->valid) {
    compute_corner_list(pyr, corners);
    corners->valid = true;
  }

#if CONFIG_MULTITHREAD
  pthread_mutex_unlock(&corners->mutex);
#endif  // CONFIG_MULTITHREAD
}

void av1_invalidate_corner_list(CornerList *corners) {
  if (corners) {
#if CONFIG_MULTITHREAD
    pthread_mutex_lock(&corners->mutex);
#endif  // CONFIG_MULTITHREAD
    corners->valid = false;
#if CONFIG_MULTITHREAD
    pthread_mutex_unlock(&corners->mutex);
#endif  // CONFIG_MULTITHREAD
  }
}

void av1_free_corner_list(CornerList *corners) {
  if (corners) {
#if CONFIG_MULTITHREAD
    pthread_mutex_destroy(&corners->mutex);
#endif  // CONFIG_MULTITHREAD
    aom_free(corners);
  }
}
