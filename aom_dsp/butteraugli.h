/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_BUTTERAUGLI_H_
#define AOM_AOM_DSP_BUTTERAUGLI_H_

#include "aom_scale/yv12config.h"

// Please note that:
// 1. We have to convert YUV to RGB first and then convert RGB to XYB, and vice
// versa.
// 2. When converting between YUV and RGB, these two functions use the
// BT.709 matrix coefficients and assume the Y, U, V samples are limited range.
// 3. When converting between RGB and XYB, these two functions assume the RGB
// color space is sRGB.
void aom_yuv_to_xyb(const YV12_BUFFER_CONFIG *yuv, YV12_BUFFER_CONFIG *xyb);
void aom_xyb_to_yuv(const YV12_BUFFER_CONFIG *xyb, YV12_BUFFER_CONFIG *yuv);

int aom_calc_butteraugli(const YV12_BUFFER_CONFIG *source,
                         const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                         float *dist_map);

#endif  // AOM_AOM_DSP_BUTTERAUGLI_H_
