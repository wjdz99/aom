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

#include "av1/common/cfl.h"

#include "./av1_rtcd.h"

void subtract_average_c(int16_t *pred_buf_q3, int width, int height,
                        int16_t avg_q3) {
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      pred_buf_q3[i] -= avg_q3;
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}
