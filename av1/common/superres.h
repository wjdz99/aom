/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_COMMON_SUPERRES_H_
#define AV1_COMMON_SUPERRES_H_

#include "./aom_config.h"

#if CONFIG_LOOP_RESTORATION
#if CONFIG_FRAME_SUPERRES

#ifdef __cplusplus
extern "C" {
#endif

void av1_frame_superres_pre_encode(AV1_COMMON *cm);
void av1_frame_superres_post_encode(AV1_COMMON *cm);
void av1_frame_superres_post_decode(AV1_COMMON *cm);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CONFIG_FRAME_SUPERRES
#endif  // CONFIG_LOOP_RESTORATION

#endif
