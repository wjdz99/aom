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

#ifndef AOM_AV1_ENCODER_PARTITIONING_H_
#define AOM_AV1_ENCODER_PARTITIONING_H_

#include "av1/encoder/encoder.h"
#include "av1/encoder/tokenize.h"

void av1_encode_nonrd_sb(AV1_COMP *cpi, ThreadData *td, TileDataEnc *tile_data,
                         TokenExtra **tp, const int mi_row, const int mi_col,
                         const int seg_skip);

void av1_encode_rd_sb(AV1_COMP *cpi, ThreadData *td, TileDataEnc *tile_data,
                      TokenExtra **tp, const int mi_row, const int mi_col,
                      const int seg_skip);

#endif  // AOM_AV1_ENCODER_PARTITIONING_H_
