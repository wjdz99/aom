/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AOM_AV1_ENCODER_AV1_SKIN_DETECTION_H_
#define AOM_AV1_ENCODER_AV1_SKIN_DETECTION_H_

#include "av1/common/blockd.h"
#include "aom_dsp/skin_detection.h"

#ifdef __cplusplus
extern "C" {
#endif

//#define OUTPUT_YUV_SKINMAP

struct AV1_COMP;

int av1_compute_skin_block(const uint8_t *y, const uint8_t *u, const uint8_t *v,
                           int stride, int strideuv, int bsize,
                           int consec_zeromv, int curr_motion_magn);

void av1_compute_skin_sb(struct AV1_COMP *const cpi, BLOCK_SIZE bsize,
                         int mi_row, int mi_col, int low_motion);

#ifdef OUTPUT_YUV_SKINMAP
// For viewing skin map on input source.
void av1_output_skin_map(struct AV1_COMP *const cpi, FILE *yuv_skinmap_file);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_AV1_SKIN_DETECTION_H_
