/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 3-Clause Clear License
 * and the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear
 * License was not distributed with this source code in the LICENSE file, you
 * can obtain it at aomedia.org/license/software-license/bsd-3-c-c/.  If the
 * Alliance for Open Media Patent License 1.0 was not distributed with this
 * source code in the PATENTS file, you can obtain it at
 * aomedia.org/license/patent-license/.
 */

#include "aom_dsp/flow_estimation/pyramid.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/bitops.h"

// TODO(rachelbarker): Move needed code from av1/ to aom_dsp/
#include "av1/common/resize.h"

#include <assert.h>
#include <string.h>

// TODO(rachelbarker):
// If original image is stored in an 8-bit buffer, then we should be able
// to point the lowest pyramid level at the same buffer, rather than copying
static INLINE ImagePyramid *alloc_pyramid(int width, int height, int n_levels) {
  assert(n_levels <= MAX_PYRAMID_LEVELS);

  // Limit number of levels on small frames
  const int msb = get_msb(AOMMIN(width, height));
  const int max_levels = AOMMAX(msb - MIN_PYRAMID_SIZE_LOG2, 1);
  n_levels = AOMMIN(n_levels, max_levels);

  ImagePyramid *pyr = aom_malloc(sizeof(*pyr));
  if (!pyr) {
    return NULL;
  }

  pyr->n_levels = n_levels;

  // Compute sizes and offsets for each pyramid level
  // These are gathered up first, so that we can allocate all pyramid levels
  // in a single buffer
  size_t buffer_size = 0;
  size_t layer_offsets[MAX_PYRAMID_LEVELS];

  // Work out if we need to allocate a few extra bytes for alignment.
  // aom_memalign() will ensure that the start of the allocation is aligned
  // to a multiple of PYRAMID_ALIGNMENT. But we want the first image pixel
  // to be aligned, not the first byte of the allocation.
  //
  // In the loop below, we ensure that the stride of every image is a multiple
  // of PYRAMID_ALIGNMENT. Thus the allocated size of each pyramid level will
  // also be a multiple of PYRAMID_ALIGNMENT. Thus, as long as we can get the
  // first pixel in the first pyramid layer aligned properly, that will
  // automatically mean that the first pixel of every row of every layer is
  // properly aligned too.
  //
  // Thus all we need to consider is the first pixel in the first layer.
  // This is located at offset
  //   extra_bytes + level_stride * PYRAMID_PADDING + PYRAMID_PADDING
  // bytes into the buffer. Since level_stride is a multiple of
  // PYRAMID_ALIGNMENT, we can ignore that. So we need
  //   extra_bytes + PYRAMID_PADDING = multiple of PYRAMID_ALIGNMENT
  //
  // To solve this, we can round PYRAMID_PADDING up to the next multiple
  // of PYRAMID_ALIGNMENT, then subtract the orginal value to calculate
  // how many extra bytes are needed.
  size_t first_px_offset =
      (PYRAMID_PADDING + PYRAMID_ALIGNMENT - 1) & ~(PYRAMID_ALIGNMENT - 1);
  size_t extra_bytes = first_px_offset - PYRAMID_PADDING;
  buffer_size += extra_bytes;

  for (int level = 0; level < n_levels; level++) {
    PyramidLayer *layer = &pyr->layers[level];

    int level_width = width >> level;
    int level_height = height >> level;

    // Allocate padding for each layer
    int padded_width = level_width + 2 * PYRAMID_PADDING;
    int padded_height = level_height + 2 * PYRAMID_PADDING;

    // Align the layer stride to be a multiple of PYRAMID_ALIGNMENT
    // This ensures that, as long as the top-left pixel in this pyramid level is
    // properly aligned, then so will the leftmost pixel in every row of the
    // pyramid level.
    int level_stride =
        (padded_width + PYRAMID_ALIGNMENT - 1) & ~(PYRAMID_ALIGNMENT - 1);

    size_t level_alloc_start = buffer_size;
    size_t level_start =
        level_alloc_start + PYRAMID_PADDING * level_stride + PYRAMID_PADDING;

    buffer_size += level_stride * padded_height;

    layer_offsets[level] = level_start;
    layer->width = level_width;
    layer->height = level_height;
    layer->stride = level_stride;
  }

  pyr->buffer_start =
      aom_memalign(PYRAMID_ALIGNMENT, buffer_size * sizeof(*pyr->buffer_start));
  if (!pyr->buffer_start) {
    aom_free(pyr);
    return NULL;
  }

  // Fill in pointers for each level
  for (int level = 0; level < n_levels; level++) {
    PyramidLayer *layer = &pyr->layers[level];
    layer->buffer = pyr->buffer_start + layer_offsets[level];
  }

  return pyr;
}

// Fill the border region of a pyramid frame.
// This must be called after the main image area is filled out.
// `img_buf` should point to the first pixel in the image area,
// ie. it should be pyr->level_buffer + pyr->level_loc[level].
static INLINE void fill_border(uint8_t *img_buf, const int width,
                               const int height, const int stride) {
  // Fill left and right areas
  for (int row = 0; row < height; row++) {
    uint8_t *row_start = &img_buf[row * stride];
    uint8_t left_pixel = row_start[0];
    memset(row_start - PYRAMID_PADDING, left_pixel, PYRAMID_PADDING);
    uint8_t right_pixel = row_start[width - 1];
    memset(row_start + width, right_pixel, PYRAMID_PADDING);
  }

  // Fill top area
  for (int row = -PYRAMID_PADDING; row < 0; row++) {
    uint8_t *row_start = &img_buf[row * stride];
    memcpy(row_start - PYRAMID_PADDING, img_buf - PYRAMID_PADDING,
           width + 2 * PYRAMID_PADDING);
  }

  // Fill bottom area
  uint8_t *last_row_start = &img_buf[(height - 1) * stride];
  for (int row = height; row < height + PYRAMID_PADDING; row++) {
    uint8_t *row_start = &img_buf[row * stride];
    memcpy(row_start - PYRAMID_PADDING, last_row_start - PYRAMID_PADDING,
           width + 2 * PYRAMID_PADDING);
  }
}

// Compute coarse to fine pyramids for a frame
static INLINE void fill_pyramid(const YV12_BUFFER_CONFIG *frm, int bit_depth,
                                ImagePyramid *frm_pyr) {
  int n_levels = frm_pyr->n_levels;
  const int frm_width = frm->y_width;
  const int frm_height = frm->y_height;
  const int frm_stride = frm->y_stride;
  assert((frm_width >> n_levels) >= 0);
  assert((frm_height >> n_levels) >= 0);

  PyramidLayer *first_layer = &frm_pyr->layers[0];
  assert(first_layer->width == frm_width);
  assert(first_layer->height == frm_height);

  uint8_t *cur_buffer = first_layer->buffer;
  int cur_width = first_layer->width;
  int cur_height = first_layer->height;
  int cur_stride = first_layer->stride;

  // Fill out the initial pyramid level
  // For frames stored in 16-bit buffers, we need to downconvert to 8 bits.
  // For frames stored in 8-bit buffers, we can just copy each row
  if (frm->flags & YV12_FLAG_HIGHBITDEPTH) {
    uint16_t *frm_buffer = CONVERT_TO_SHORTPTR(frm->y_buffer);
    uint8_t *pyr_buffer = cur_buffer;
    for (int y = 0; y < frm_height; y++) {
      uint16_t *frm_row = frm_buffer + y * frm_stride;
      uint8_t *pyr_row = pyr_buffer + y * cur_stride;
      for (int x = 0; x < frm_width; x++) {
        pyr_row[x] = frm_row[x] >> (bit_depth - 8);
      }
    }
  } else {
    uint8_t *frm_buffer = frm->y_buffer;
    uint8_t *pyr_buffer = cur_buffer;
    for (int y = 0; y < frm_height; y++) {
      uint8_t *frm_row = frm_buffer + y * frm_stride;
      uint8_t *pyr_row = pyr_buffer + y * cur_stride;
      memcpy(pyr_row, frm_row, frm_width);
    }
  }

  fill_border(cur_buffer, cur_width, cur_height, cur_stride);

  // Start at the finest level and resize down to the coarsest level
  for (int level = 1; level < n_levels; ++level) {
    PyramidLayer *next_layer = &frm_pyr->layers[level];
    uint8_t *next_buffer = next_layer->buffer;
    int next_width = next_layer->width;
    int next_height = next_layer->height;
    int next_stride = next_layer->stride;

    // Compute the next pyramid level by downsampling the current level.
    //
    // We downsample by a factor of exactly 2, clipping the rightmost and
    // bottommost pixel off of the current level if needed. We do this for
    // two main reasons:
    //
    // 1) In the disflow code, when stepping from a higher pyramid level to a
    //    lower pyramid level, we need to not just interpolate the flow field
    //    but also to scale each flow vector by the upsampling ratio.
    //    So it is much more convenient if this ratio is simply 2.
    //
    // 2) Up/downsampling by a factor of 2 can be implemented much more
    //    efficiently than up/downsampling by a generic ratio.
    //    TODO(rachelbarker): Implement optimized downsample-by-2 function
    av1_resize_plane(cur_buffer, next_height << 1, next_width << 1, cur_stride,
                     next_buffer, next_height, next_width, next_stride);
    fill_border(next_buffer, next_width, next_height, next_stride);

    cur_buffer = next_buffer;
    cur_width = next_width;
    cur_height = next_height;
    cur_stride = next_stride;
  }
}

// Allocate and fill out a pyramid structure for a given frame
ImagePyramid *aom_compute_pyramid(struct yv12_buffer_config *frm, int bit_depth,
                                  int n_levels) {
  if (frm->y_pyramid) {
    // Already computed, no need to do it again
    return frm->y_pyramid;
  }

  frm->y_pyramid = alloc_pyramid(frm->y_width, frm->y_height, n_levels);
  fill_pyramid(frm, bit_depth, frm->y_pyramid);
  return frm->y_pyramid;
}

void aom_free_pyramid(ImagePyramid *pyr) {
  if (pyr) {
    aom_free(pyr->buffer_start);
    aom_free(pyr);
  }
}
