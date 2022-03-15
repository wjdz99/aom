/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_TUNE_VISUAL_MASKING_H_
#define AOM_AV1_ENCODER_TUNE_VISUAL_MASKING_H_

#include "aom_scale/yv12config.h"
#include "av1/common/enums.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/block.h"
#include "av1/encoder/encoder.h"

void av1_set_visual_masking_rdmult(const AV1_COMP *const cpi, int *errorperbit,
                                   const BLOCK_SIZE bsize, const int mi_row,
                                   const int mi_col, int *const rdmult);

void av1_set_mb_visual_masking_rdmult_scaling(AV1_COMP *cpi);

#endif  // AOM_AV1_ENCODER_TUNE_VISUAL_MASKING_H_
