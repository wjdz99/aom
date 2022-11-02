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
#ifndef AOM_FLOW_ESTIMATION_DISFLOW_BASED_H_
#define AOM_FLOW_ESTIMATION_DISFLOW_BASED_H_

#include <stdbool.h>

#include "aom_dsp/flow_estimation/flow_estimation.h"
#include "aom_dsp/rect.h"
#include "aom_scale/yv12config.h"

#ifdef __cplusplus
extern "C" {
#endif

// Number of pyramid levels in disflow computation
#define DISFLOW_PYRAMID_LEVELS 12

// Size of square patches in the disflow dense grid
// Must be a power of 2
#define DISFLOW_PATCH_SIZE_LOG2 3
#define DISFLOW_PATCH_SIZE (1 << DISFLOW_PATCH_SIZE_LOG2)
// Center point of square patch
#define DISFLOW_PATCH_CENTER ((DISFLOW_PATCH_SIZE / 2) - 1)

// Warp error convergence threshold for disflow
// This is specified as a mean squared error per pixel, but is converted
// to a summed squared error per patch for efficiency
// TODO(rachelbarker): Work out a good value for DISFLOW_MSE_TR
// Its current value of 0.01 is effectively 0; something more on the order
// of 1 would probably be sensible.
#define DISFLOW_MSE_TR 0.01
#define DISFLOW_SSE_TR \
  ((int)(DISFLOW_MSE_TR * DISFLOW_PATCH_SIZE * DISFLOW_PATCH_SIZE + 0.5))
// Max number of iterations if warp convergence is not found
#define DISFLOW_MAX_ITR 10

FlowField *aom_alloc_flow_field(int width, int height, int stride);
void aom_free_flow_field(FlowField *flow);

FlowField *aom_compute_flow_field(YV12_BUFFER_CONFIG *frm,
                                  YV12_BUFFER_CONFIG *ref, int bit_depth);

bool aom_fit_global_model_to_flow_field(FlowField *flow,
                                        TransformationType type,
                                        YV12_BUFFER_CONFIG *frm, int bit_depth,
                                        MotionModel *params_by_motion,
                                        int num_motions);

bool aom_fit_local_model_to_flow_field(const FlowField *flow,
                                       const PixelRect *rect,
                                       TransformationType type, double *mat);

#ifdef __cplusplus
}
#endif

#endif  // AOM_FLOW_ESTIMATION_DISFLOW_BASED_H_
