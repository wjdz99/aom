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

#include <assert.h>

#include "av1/common/filter.h"

DECLARE_ALIGNED(256, static const int16_t,
                bilinear_filters[SUBPEL_SHIFTS][8]) = {
  { 0, 0, 0, 128, 0, 0, 0, 0 },  { 0, 0, 0, 120, 8, 0, 0, 0 },
  { 0, 0, 0, 112, 16, 0, 0, 0 }, { 0, 0, 0, 104, 24, 0, 0, 0 },
  { 0, 0, 0, 96, 32, 0, 0, 0 },  { 0, 0, 0, 88, 40, 0, 0, 0 },
  { 0, 0, 0, 80, 48, 0, 0, 0 },  { 0, 0, 0, 72, 56, 0, 0, 0 },
  { 0, 0, 0, 64, 64, 0, 0, 0 },  { 0, 0, 0, 56, 72, 0, 0, 0 },
  { 0, 0, 0, 48, 80, 0, 0, 0 },  { 0, 0, 0, 40, 88, 0, 0, 0 },
  { 0, 0, 0, 32, 96, 0, 0, 0 },  { 0, 0, 0, 24, 104, 0, 0, 0 },
  { 0, 0, 0, 16, 112, 0, 0, 0 }, { 0, 0, 0, 8, 120, 0, 0, 0 }
};

// Lagrangian interpolation filter
DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_8[SUBPEL_SHIFTS][8]) = {
  { 0, 0, 0, 128, 0, 0, 0, 0 },      { 0, 2, -6, 126, 8, -2, 0, 0 },
  { 0, 2, -10, 122, 18, -4, 0, 0 },  { 0, 2, -12, 116, 28, -8, 2, 0 },
  { 0, 2, -14, 110, 38, -10, 2, 0 }, { 0, 2, -14, 102, 48, -12, 2, 0 },
  { 0, 2, -16, 94, 58, -12, 2, 0 },  { 0, 2, -14, 84, 66, -12, 2, 0 },
  { 0, 2, -14, 76, 76, -14, 2, 0 },  { 0, 2, -12, 66, 84, -14, 2, 0 },
  { 0, 2, -12, 58, 94, -16, 2, 0 },  { 0, 2, -12, 48, 102, -14, 2, 0 },
  { 0, 2, -10, 38, 110, -14, 2, 0 }, { 0, 2, -8, 28, 116, -12, 2, 0 },
  { 0, 0, -4, 18, 122, -10, 2, 0 },  { 0, 0, -2, 8, 126, -6, 2, 0 }
};

// DCT based filter
DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_8sharp[SUBPEL_SHIFTS][8]) = {
  { 0, 0, 0, 128, 0, 0, 0, 0 },         { -2, 2, -6, 126, 8, -2, 2, 0 },
  { -2, 6, -12, 124, 16, -6, 4, -2 },   { -2, 8, -18, 120, 26, -10, 6, -2 },
  { -4, 10, -22, 116, 38, -14, 6, -2 }, { -4, 10, -22, 108, 48, -18, 8, -2 },
  { -4, 10, -24, 100, 60, -20, 8, -2 }, { -4, 10, -24, 90, 70, -22, 10, -2 },
  { -4, 12, -24, 80, 80, -24, 12, -4 }, { -2, 10, -22, 70, 90, -24, 10, -4 },
  { -2, 8, -20, 60, 100, -24, 10, -4 }, { -2, 8, -18, 48, 108, -22, 10, -4 },
  { -2, 6, -14, 38, 116, -22, 10, -4 }, { -2, 6, -10, 26, 120, -18, 8, -2 },
  { -2, 4, -6, 16, 124, -12, 6, -2 },   { 0, 2, -2, 8, 126, -6, 2, -2 }
};

#if CONFIG_EXT_INTERP
DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_10sharp[SUBPEL_SHIFTS][10]) = {
  // intfilt 0.77
  { 0, 0, 0, 0, 128, 0, 0, 0, 0, 0 },
  { 0, -1, 3, -6, 127, 8, -4, 2, -1, 0 },
  { 1, -2, 5, -12, 124, 18, -7, 3, -2, 0 },
  { 1, -3, 7, -17, 119, 28, -11, 5, -2, 1 },
  { 1, -4, 8, -20, 114, 38, -14, 7, -3, 1 },
  { 1, -4, 9, -22, 107, 49, -17, 8, -4, 1 },
  { 2, -5, 10, -24, 99, 59, -20, 9, -4, 2 },
  { 2, -5, 10, -24, 90, 70, -22, 10, -5, 2 },
  { 2, -5, 10, -23, 80, 80, -23, 10, -5, 2 },
  { 2, -5, 10, -22, 70, 90, -24, 10, -5, 2 },
  { 2, -4, 9, -20, 59, 99, -24, 10, -5, 2 },
  { 1, -4, 8, -17, 49, 107, -22, 9, -4, 1 },
  { 1, -3, 7, -14, 38, 114, -20, 8, -4, 1 },
  { 1, -2, 5, -11, 28, 119, -17, 7, -3, 1 },
  { 0, -2, 3, -7, 18, 124, -12, 5, -2, 1 },
  { 0, -1, 2, -4, 8, 127, -6, 3, -1, 0 },
};

DECLARE_ALIGNED(16, static const int16_t,
                sub_pel_filters_12sharp[SUBPEL_SHIFTS][12]) = {
  // intfilt 0.85
  { 0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0 },
  { 0, 1, -2, 3, -7, 127, 8, -4, 2, -1, 1, 0 },
  { -1, 2, -3, 6, -13, 124, 18, -8, 4, -2, 2, -1 },
  { -1, 3, -4, 8, -18, 120, 28, -12, 7, -4, 2, -1 },
  { -1, 3, -6, 10, -21, 115, 38, -15, 8, -5, 3, -1 },
  { -2, 4, -6, 12, -24, 108, 49, -18, 10, -6, 3, -2 },
  { -2, 4, -7, 13, -25, 100, 60, -21, 11, -7, 4, -2 },
  { -2, 4, -7, 13, -26, 91, 71, -24, 13, -7, 4, -2 },
  { -2, 4, -7, 13, -25, 81, 81, -25, 13, -7, 4, -2 },
  { -2, 4, -7, 13, -24, 71, 91, -26, 13, -7, 4, -2 },
  { -2, 4, -7, 11, -21, 60, 100, -25, 13, -7, 4, -2 },
  { -2, 3, -6, 10, -18, 49, 108, -24, 12, -6, 4, -2 },
  { -1, 3, -5, 8, -15, 38, 115, -21, 10, -6, 3, -1 },
  { -1, 2, -4, 7, -12, 28, 120, -18, 8, -4, 3, -1 },
  { -1, 2, -2, 4, -8, 18, 124, -13, 6, -3, 2, -1 },
  { 0, 1, -1, 2, -4, 8, 127, -7, 3, -2, 1, 0 },
};

DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_8smooth[SUBPEL_SHIFTS][8]) = {
  // freqmultiplier = 0.75
  { 0, 0, 0, 128, 0, 0, 0, 0 },     { 2, -10, 19, 95, 31, -11, 2, 0 },
  { 2, -9, 14, 94, 37, -12, 2, 0 }, { 2, -8, 9, 92, 43, -12, 1, 1 },
  { 2, -7, 5, 90, 49, -12, 1, 0 },  { 2, -5, 1, 86, 55, -12, 0, 1 },
  { 1, -4, -2, 82, 61, -11, 0, 1 }, { 1, -3, -5, 77, 67, -9, -1, 1 },
  { 1, -2, -7, 72, 72, -7, -2, 1 }, { 1, -1, -9, 67, 77, -5, -3, 1 },
  { 1, 0, -11, 61, 82, -2, -4, 1 }, { 1, 0, -12, 55, 86, 1, -5, 2 },
  { 0, 1, -12, 49, 90, 5, -7, 2 },  { 1, 1, -12, 43, 92, 9, -8, 2 },
  { 0, 2, -12, 37, 94, 14, -9, 2 }, { 0, 2, -11, 31, 95, 19, -10, 2 },
};

DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_8smooth2[SUBPEL_SHIFTS][8]) = {
  // freqmultiplier = 0.35
  { 0, 0, 0, 128, 0, 0, 0, 0 },     { -1, 8, 31, 47, 34, 10, 0, -1 },
  { -1, 7, 29, 46, 36, 12, 0, -1 }, { -1, 6, 28, 46, 37, 13, 0, -1 },
  { -1, 5, 26, 46, 38, 14, 1, -1 }, { -1, 4, 25, 45, 39, 16, 1, -1 },
  { -1, 4, 23, 44, 41, 17, 1, -1 }, { -1, 3, 21, 44, 42, 18, 2, -1 },
  { -1, 2, 20, 43, 43, 20, 2, -1 }, { -1, 2, 18, 42, 44, 21, 3, -1 },
  { -1, 1, 17, 41, 44, 23, 4, -1 }, { -1, 1, 16, 39, 45, 25, 4, -1 },
  { -1, 1, 14, 38, 46, 26, 5, -1 }, { -1, 0, 13, 37, 46, 28, 6, -1 },
  { -1, 0, 12, 36, 46, 29, 7, -1 }, { -1, 0, 10, 34, 47, 31, 8, -1 },
};

#else  // CONFIG_EXT_INTERP

// freqmultiplier = 0.5
DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_8smooth[SUBPEL_SHIFTS][8]) = {
  { 0, 0, 0, 128, 0, 0, 0, 0 },     { 0, 2, 28, 62, 34, 2, 0, 0 },
  { 0, 0, 26, 62, 36, 4, 0, 0 },    { 0, 0, 22, 62, 40, 4, 0, 0 },
  { 0, 0, 20, 60, 42, 6, 0, 0 },    { 0, 0, 18, 58, 44, 8, 0, 0 },
  { 0, 0, 16, 56, 46, 10, 0, 0 },   { 0, -2, 16, 54, 48, 12, 0, 0 },
  { 0, -2, 14, 52, 52, 14, -2, 0 }, { 0, 0, 12, 48, 54, 16, -2, 0 },
  { 0, 0, 10, 46, 56, 16, 0, 0 },   { 0, 0, 8, 44, 58, 18, 0, 0 },
  { 0, 0, 6, 42, 60, 20, 0, 0 },    { 0, 0, 4, 40, 62, 22, 0, 0 },
  { 0, 0, 4, 36, 62, 26, 0, 0 },    { 0, 0, 2, 34, 62, 28, 2, 0 }
};

#endif  // CONFIG_EXT_INTERP

const InterpKernel *av1_filter_kernels[4] = { sub_pel_filters_8,
                                              sub_pel_filters_8smooth,
                                              sub_pel_filters_8sharp,
                                              bilinear_filters };
#if CONFIG_EXT_INTERP
static const InterpFilterParams interp_filter_params_list[SWITCHABLE_FILTERS +
                                                          1] = {
  { (const int16_t *)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP },
  { (const int16_t *)sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SMOOTH },
  { (const int16_t *)sub_pel_filters_10sharp, 10, SUBPEL_SHIFTS,
    MULTITAP_SHARP },
  { (const int16_t *)sub_pel_filters_8smooth2, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SMOOTH2 },
  { (const int16_t *)sub_pel_filters_12sharp, 12, SUBPEL_SHIFTS,
    MULTITAP_SHARP2 },
  { (const int16_t *)bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR }
};
#else   // CONFIG_EXT_INTERP
static const InterpFilterParams interp_filter_params_list[SWITCHABLE_FILTERS +
                                                          1] = {
  { (const int16_t *)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP },
  { (const int16_t *)sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SMOOTH },
  { (const int16_t *)sub_pel_filters_8sharp, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SHARP },
  { (const int16_t *)bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR }
};
#endif  // CONFIG_EXT_INTERP

InterpFilterParams get_interp_filter_params(InterpFilter interp_filter) {
  InterpFilterParams params = interp_filter_params_list[interp_filter];
  assert(params.interp_filter == interp_filter);
  return params;
}
