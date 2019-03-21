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

#include <stdio.h>

#include "av1/common/cnn.h"
#include "av1/common/onyxc_int.h"

#include "av1/models/intra_frame_model/trial_model.h"

void cnn_restoration(AV1_COMMON *cm) {
  // TODO(logangw): ADD INFRASTRUCTURE TO CHOOSE MODELS
  av1_restore_cnn_plane(cm, &model, AOM_PLANE_Y);
}
