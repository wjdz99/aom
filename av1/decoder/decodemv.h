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

#ifndef AV1_DECODER_DECODEMV_H_
#define AV1_DECODER_DECODEMV_H_

#include "aom_dsp/bitreader.h"

#include "av1/decoder/decoder.h"

#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
void av1_read_mode_info(AV1Decoder *const pbi, MACROBLOCKD *xd,
#if CONFIG_SUPERTX
                        int supertx_enabled,
#endif

                        int mi_row, int mi_col, aom_reader *r, int x_mis,
                        int y_mis);
=======
void av1_read_mode_info(AV1Decoder *const pbi, MACROBLOCKD *xd, int mi_row,
                        int mi_col, aom_reader *r, int x_mis, int y_mis);
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_DECODER_DECODEMV_H_
