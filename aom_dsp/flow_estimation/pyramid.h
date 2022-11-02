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

#ifndef AOM_FLOW_ESTIMATION_PYRAMID_H_
#define AOM_FLOW_ESTIMATION_PYRAMID_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

#include "config/aom_config.h"

// Maximum number of pyramid levels
#if CONFIG_GM_USE_DISFLOW
// Disflow requires as many pyramid levels as possible, subject
// to the limit imposed by the frame size and MIN_PYRAMID_SIZE.
// At the same time, some allocations are sized according to
// MAX_PYRAMID_LEVELS, so we don't want to make it too big.
//
// With MIN_PYRAMID_SIZE == 8, 12 pyramid levels is enough for
// frames up to 32k on their shortest side.
#define MAX_PYRAMID_LEVELS 12
#else
// Feature based code only requires one pyramid level
#define MAX_PYRAMID_LEVELS 1
#endif  // CONFIG_GM_USE_DISFLOW

// Minimum dimensions of a downsampled image
#define MIN_PYRAMID_SIZE_LOG2 3
#define MIN_PYRAMID_SIZE (1 << MIN_PYRAMID_SIZE_LOG2)

// Size of border around each pyramid image, in pixels
// Similarly to the border around regular image buffers, this border is filled
// with copies of the outermost pixels of the frame, to allow for more efficient
// convolution code
// TODO(rachelbarker): How many pixels do we actually need here?
// I think we only need 9 for disflow, but how many for corner matching?
#define PYRAMID_PADDING 16

// Byte alignment of each line within the image pyramids.
// That is, the first pixel inside the image (ie, not in the border region),
// on each row of each pyramid level, is aligned to this byte alignment.
// This value must be a power of 2.
#define PYRAMID_ALIGNMENT 32

// Forward declare this struct rather than including "aom_scale/yv12config.h",
// so that that file can include this one without causing circular dependencies
struct yv12_buffer_config;

typedef struct {
  uint8_t *buffer;
  int width;
  int height;
  int stride;
} PyramidLayer;

// Struct for an image pyramid
// TODO(rachelbarker):
// Allocate "layers" array dynamically, and store the number of layers
// configured.
// Then, if we start with a small number of layers and later need more,
// we can easily reallocate
typedef struct {
  int n_levels;
  // Pointer to start of allocated buffer
  uint8_t *buffer_start;
  // Data for each level
  PyramidLayer layers[MAX_PYRAMID_LEVELS];
} ImagePyramid;

// Allocate and fill out a downsampling pyramid for a given frame.
//
// The top level (index 0) will always be an 8-bit copy of the input frame,
// regardless of the input bit depth. Additional levels are then downscaled
// by powers of 2.
//
// Note on n_levels:
// * For feature-based global motion, n_levels need only be 1,
//   which just constructs an 8-bit version of the input frame.
// * For disflow-based global motion, n_levels should equal
//   DISFLOW_PYRAMID_LEVELS
// * In any case, n_levels must be <= MAX_PYRAMID_LEVELS
//
// For small input frames, the number of levels actually constructed
// will be limited so that the smallest image is at least MIN_PYRAMID_SIZE
// pixels along each side.
//
// However, if the input frame has a side of length < MIN_PYRAMID_SIZE,
// we will still construct the top level.
ImagePyramid *aom_compute_pyramid(struct yv12_buffer_config *frm, int bit_depth,
                                  int n_levels);

void aom_free_pyramid(ImagePyramid *pyr);

#ifdef __cplusplus
}
#endif

#endif  // AOM_FLOW_ESTIMATION_PYRAMID_H_
