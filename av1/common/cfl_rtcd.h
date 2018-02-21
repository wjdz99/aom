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

#ifndef AV1_COMMON_CFL_RTCD_H_
#define AV1_COMMON_CFL_RTCD_H_

#include "av1/common/enums.h"

typedef void (*cfl_subsample_lbd_fn)(const uint8_t *input, int input_stride,
                                     int16_t *output_q3, int width, int height);

typedef void (*cfl_subsample_hbd_fn)(const uint16_t *input, int input_stride,
                                     int16_t *output_q3, int width, int height);

typedef void (*cfl_subtract_average_fn)(int16_t *pred_buf_q3);

typedef void (*cfl_predict_lbd_fn)(const int16_t *pred_buf_q3, uint8_t *dst,
                                   int dst_stride, TX_SIZE tx_size,
                                   int alpha_q3);

typedef void (*cfl_predict_hbd_fn)(const int16_t *pred_buf_q3, uint16_t *dst,
                                   int dst_stride, TX_SIZE tx_size,
                                   int alpha_q3, int bd);

#endif  // AV1_COMMON_CFL_RTCD_H_
