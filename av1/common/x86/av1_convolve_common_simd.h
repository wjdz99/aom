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

#include "./aom_config.h"

typedef const int8_t (*SubpelFilterCoeffs)[16];
#if CONFIG_AOM_HIGHBITDEPTH
typedef const int16_t (*HbdSubpelFilterCoeffs)[8];
#endif

static INLINE SubpelFilterCoeffs
av1_get_subpel_filter_signal_dir(const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}

static INLINE SubpelFilterCoeffs
av1_get_subpel_filter_ver_signal_dir(const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_ver_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_ver_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_ver_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}

#if CONFIG_AOM_HIGHBITDEPTH
static INLINE HbdSubpelFilterCoeffs av1_hbd_get_subpel_filter_ver_signal_dir(
    const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSE4_1
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_highbd_ver_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_highbd_ver_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSE4_1
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_highbd_ver_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}
#endif
