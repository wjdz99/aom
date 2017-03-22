/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 *
 */

#include "./aom_config.h"

#if CONFIG_LOOP_RESTORATION
#if CONFIG_FRAME_SUPERRES

// TODO(afergs): call before decode
void av1_frame_superres_pre_encode(AV1_COMMON *cm) {

  // TODO(afergs): decide whether to scale (return if not)

  // TODO(afergs): reserve buffers for smaller image

  // TODO(afergs): keep copy of high resolution source in cm somewhere

  // TODO(afergs): downscale

}

// TODO(afergs): call after encode (after encoder test decode?)
void av1_frame_superres_post_encode(AV1_COMMON *cm) {

  return; // TODO(afergs): [STAGE II] implement in second stage
          // no Wiener filter in stage 1

  // TODO(afergs): decode encoded image (done elsewhere?)

  // TODO(afergs): upscale decoded image

  // TODO(afergs): fit Wiener filter to images, save params to cm

}


// TODO(afergs): call after decode
void av1_frame_superres_post_decode(AV1_COMMON *cm) {

  // TODO(afergs): upscale decoded image

  // TODO(afergs): [STAGE II]

}

#endif  // CONFIG_FRAME_SUPERRES
#endif  // CONFIG_LOOP_RESTORATION
