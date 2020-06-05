/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_COMMON_CNN_WRAPPER_H_
#define AOM_AV1_COMMON_CNN_WRAPPER_H_

#include "av1/common/onyxc_int.h"
#include "av1/common/resize.h"
#include "av1/encoder/firstpass.h"

#ifdef __cplusplus
extern "C" {
#endif

// Minimum base_qindex needed to run cnn.
#define MIN_CNN_Q_INDEX 67

static INLINE int av1_use_cnn_encode(const AV1_COMMON *cm,
                                     const GF_GROUP *gf_group) {
  bool is_overlay_update = false;
  if (gf_group->update_type[gf_group->index] == OVERLAY_UPDATE ||
      gf_group->update_type[gf_group->index] == INTNL_OVERLAY_UPDATE) {
    is_overlay_update = true;
  }

  return ((cm->base_qindex > MIN_CNN_Q_INDEX) && !av1_superres_scaled(cm) &&
          !is_overlay_update);
}

static INLINE int av1_use_cnn(const AV1_COMMON *cm) {
  return ((cm->base_qindex > MIN_CNN_Q_INDEX) && !av1_superres_scaled(cm));
}

void av1_encode_restore_cnn(AV1_COMMON *cm, AVxWorker *workers,
                            int num_workers);

void av1_decode_restore_cnn(AV1_COMMON *cm, AVxWorker *workers,
                            int num_workers);

const CNN_CONFIG *av1_get_cnn_config_from_qindex(int qindex,
                                                 FRAME_TYPE frame_type);

#ifdef __cplusplus
}
#endif
#endif  // AOM_AV1_COMMON_CNN_WRAPPER_H_
