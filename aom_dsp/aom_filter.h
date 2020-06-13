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

#ifndef AOM_AOM_DSP_AOM_FILTER_H_
#define AOM_AOM_DSP_AOM_FILTER_H_

#include "aom/aom_integer.h"
#include "av1/common/filter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BIL_SUBPEL_BITS 3
#define BIL_SUBPEL_SHIFTS (1 << BIL_SUBPEL_BITS)

// 2 tap bilinear filters
static const uint8_t bilinear_filters_2t[BIL_SUBPEL_SHIFTS][2] = {
  { 128, 0 }, { 112, 16 }, { 96, 32 }, { 80, 48 },
  { 64, 64 }, { 48, 80 },  { 32, 96 }, { 16, 112 },
};

static INLINE const int16_t *av1_get_interp_filter_kernel(
    const InterpFilter interp_filter, int subpel_search) {
  assert(subpel_search >= USE_2_TAPS);
  return (subpel_search == USE_2_TAPS)
             ? av1_interp_4tap[BILINEAR].filter_ptr
             : ((subpel_search == USE_4_TAPS)
                    ? av1_interp_4tap[interp_filter].filter_ptr
                    : av1_interp_filter_params_list[interp_filter].filter_ptr);
}

static INLINE const InterpFilterParams *av1_get_filter(int subpel_search) {
  assert(subpel_search >= USE_2_TAPS);

  switch (subpel_search) {
    case USE_2_TAPS: return &av1_interp_4tap[BILINEAR];
    case USE_4_TAPS: return &av1_interp_4tap[EIGHTTAP_REGULAR];
    case USE_8_TAPS: return &av1_interp_filter_params_list[EIGHTTAP_REGULAR];
    default: assert(0); return NULL;
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AOM_DSP_AOM_FILTER_H_
