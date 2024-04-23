/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_SKIN_DETECTION_H_
#define AOM_AOM_DSP_SKIN_DETECTION_H_

#ifdef __cplusplus
extern "C" {
#endif

int aom_skin_pixel(const int y, const int cb, const int cr, int motion);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AOM_DSP_SKIN_DETECTION_H_
