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
#include <memory.h>
#include <math.h>

#include "config/aom_dsp_rtcd.h"

#include "aom_dsp/flow_estimation/corner_detect.h"
#include "aom_dsp/flow_estimation/corner_match.h"
#include "aom_dsp/flow_estimation/disflow.h"
#include "aom_dsp/flow_estimation/flow_estimation.h"
#include "aom_dsp/flow_estimation/ransac.h"
#include "aom_dsp/pyramid.h"
#include "aom_scale/yv12config.h"

// Compute SAD over a window that isn't necessarily 16px-aligned
static INLINE unsigned int generic_sad(const uint8_t *const ref, int ref_stride,
                                       const uint8_t *const dst, int dst_stride,
                                       int p_width, int p_height) {
  unsigned int sad = 0;
  for (int i = 0; i < p_height; ++i) {
    for (int j = 0; j < p_width; ++j) {
      sad += abs(dst[j + i * dst_stride] - ref[j + i * ref_stride]);
    }
  }
  return sad;
}

static int is_eligible_point(int pointx, int pointy, int width, int height) {
  return (pointx >= MATCH_SZ_BY2 && pointy >= MATCH_SZ_BY2 &&
          pointx + MATCH_SZ_BY2 < width && pointy + MATCH_SZ_BY2 < height);
}

static int is_eligible_distance(int point1x, int point1y, int point2x,
                                int point2y, int width, int height) {
  const int thresh = (width < height ? height : width) >> 4;
  return ((point1x - point2x) * (point1x - point2x) +
          (point1y - point2y) * (point1y - point2y)) <= thresh * thresh;
}

typedef struct {
  int x;
  int y;
  int best_match_idx;
  unsigned int best_match_sad;
} PointInfo;

static int determine_correspondence(const unsigned char *src,
                                    const int *src_corners, int num_src_corners,
                                    const unsigned char *ref,
                                    const int *ref_corners, int num_ref_corners,
                                    int width, int height, int src_stride,
                                    int ref_stride,
                                    Correspondence *correspondences) {
  PointInfo *src_point_info = NULL;
  PointInfo *ref_point_info = NULL;
  int num_correspondences = 0;

  src_point_info =
      (PointInfo *)aom_calloc(num_src_corners, sizeof(*src_point_info));
  if (!src_point_info) {
    goto finished;
  }

  ref_point_info =
      (PointInfo *)aom_calloc(num_ref_corners, sizeof(*ref_point_info));
  if (!ref_point_info) {
    goto finished;
  }

  // First pass (linear):
  // Filter corner lists and compute per-patch means and standard deviations,
  // for the src and ref frames independently
  int src_point_count = 0;
  for (int i = 0; i < num_src_corners; i++) {
    int src_x = src_corners[2 * i];
    int src_y = src_corners[2 * i + 1];
    if (!is_eligible_point(src_x, src_y, width, height)) continue;

    PointInfo *point = &src_point_info[src_point_count];
    point->x = src_x;
    point->y = src_y;
    point->best_match_sad = UINT_MAX;
    src_point_count++;
  }
  if (src_point_count == 0) {
    goto finished;
  }

  int ref_point_count = 0;
  for (int j = 0; j < num_ref_corners; j++) {
    int ref_x = ref_corners[2 * j];
    int ref_y = ref_corners[2 * j + 1];
    if (!is_eligible_point(ref_x, ref_y, width, height)) continue;

    PointInfo *point = &ref_point_info[ref_point_count];
    point->x = ref_x;
    point->y = ref_y;
    point->best_match_sad = UINT_MAX;
    ref_point_count++;
  }
  if (ref_point_count == 0) {
    goto finished;
  }

  // Second pass (quadratic):
  // For each pair of points, compute correlation, and use this to determine
  // the best match of each corner, in both directions
  for (int i = 0; i < src_point_count; ++i) {
    PointInfo *src_point = &src_point_info[i];
    for (int j = 0; j < ref_point_count; ++j) {
      PointInfo *ref_point = &ref_point_info[j];
      const int sx = src_point->x;
      const int sy = src_point->y;
      const int rx = ref_point->x;
      const int ry = ref_point->y;

      if (!is_eligible_distance(sx, sy, rx, ry, width, height)) continue;

      const uint8_t *src_tl =
          src + (sy - MATCH_SZ_BY2) * src_stride + (sx - MATCH_SZ_BY2);
      const uint8_t *ref_tl =
          ref + (ry - MATCH_SZ_BY2) * ref_stride + (rx - MATCH_SZ_BY2);

      // TODO(rachelbarker): Modify aom_sad16x16_sse2 to not require the
      // source pointer to be 16-byte aligned, then use aom_sad16x16() here
      const unsigned int match_sad = generic_sad(
          src_tl, src_stride, ref_tl, ref_stride, MATCH_SZ, MATCH_SZ);
      if (match_sad < src_point->best_match_sad) {
        src_point->best_match_idx = j;
        src_point->best_match_sad = match_sad;
      }
      if (match_sad < ref_point->best_match_sad) {
        ref_point->best_match_idx = i;
        ref_point->best_match_sad = match_sad;
      }
    }
  }

  // Third pass (linear):
  // Scan through source corners, generating a correspondence for each corner
  // iff ref_best_match[src_best_match[i]] == i
  // Then refine the generated correspondences using optical flow
  for (int i = 0; i < src_point_count; i++) {
    PointInfo *point = &src_point_info[i];

    // Skip points with no match candidates
    if (point->best_match_sad == UINT_MAX) continue;

    PointInfo *match_point = &ref_point_info[point->best_match_idx];
    if (match_point->best_match_idx == i) {
      // Refine match using optical flow and store
      const int sx = point->x;
      const int sy = point->y;
      const int rx = match_point->x;
      const int ry = match_point->y;
      double u = (double)(rx - sx);
      double v = (double)(ry - sy);

      const int patch_tl_x = sx - DISFLOW_PATCH_CENTER;
      const int patch_tl_y = sy - DISFLOW_PATCH_CENTER;
      aom_compute_flow_at_point(src, ref, patch_tl_x, patch_tl_y, width, height,
                                src_stride, &u, &v);

      Correspondence *correspondence = &correspondences[num_correspondences];
      correspondence->x = (double)sx;
      correspondence->y = (double)sy;
      correspondence->rx = (double)sx + u;
      correspondence->ry = (double)sy + v;
      num_correspondences++;
    }
  }

finished:
  aom_free(src_point_info);
  aom_free(ref_point_info);
  return num_correspondences;
}

bool av1_compute_global_motion_feature_match(
    TransformationType type, YV12_BUFFER_CONFIG *src, YV12_BUFFER_CONFIG *ref,
    int bit_depth, MotionModel *motion_models, int num_motion_models,
    bool *mem_alloc_failed) {
  int num_correspondences;
  Correspondence *correspondences;
  ImagePyramid *src_pyramid = src->y_pyramid;
  CornerList *src_corners = src->corners;
  ImagePyramid *ref_pyramid = ref->y_pyramid;
  CornerList *ref_corners = ref->corners;

  // Precompute information we will need about each frame
  if (!aom_compute_pyramid(src, bit_depth, src_pyramid)) {
    *mem_alloc_failed = true;
    return false;
  }
  if (!av1_compute_corner_list(src_pyramid, src_corners)) {
    *mem_alloc_failed = true;
    return false;
  }
  if (!aom_compute_pyramid(ref, bit_depth, ref_pyramid)) {
    *mem_alloc_failed = true;
    return false;
  }
  if (!av1_compute_corner_list(ref_pyramid, ref_corners)) {
    *mem_alloc_failed = true;
    return false;
  }

  const uint8_t *src_buffer = src_pyramid->layers[0].buffer;
  const int src_width = src_pyramid->layers[0].width;
  const int src_height = src_pyramid->layers[0].height;
  const int src_stride = src_pyramid->layers[0].stride;

  const uint8_t *ref_buffer = ref_pyramid->layers[0].buffer;
  assert(ref_pyramid->layers[0].width == src_width);
  assert(ref_pyramid->layers[0].height == src_height);
  const int ref_stride = ref_pyramid->layers[0].stride;

  // find correspondences between the two images
  correspondences = (Correspondence *)aom_malloc(src_corners->num_corners *
                                                 sizeof(*correspondences));
  if (!correspondences) {
    *mem_alloc_failed = true;
    return false;
  }
  num_correspondences = determine_correspondence(
      src_buffer, src_corners->corners, src_corners->num_corners, ref_buffer,
      ref_corners->corners, ref_corners->num_corners, src_width, src_height,
      src_stride, ref_stride, correspondences);

  bool result = ransac(correspondences, num_correspondences, type,
                       motion_models, num_motion_models, mem_alloc_failed);

  aom_free(correspondences);
  return result;
}
