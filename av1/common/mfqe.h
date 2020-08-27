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

#ifndef AOM_AV1_COMMON_MFQE_H_
#define AOM_AV1_COMMON_MFQE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "config/av1_rtcd.h"

#include "av1/common/onyxc_int.h"

typedef struct y_buffer_config {
  uint8_t *buffer;
  int stride;
  int height;
  int width;
} Y_BUFFER_CONFIG;

typedef struct mv_mfqe {
  MV mv;
  int subpel_x_qn;
  int subpel_y_qn;
  uint8_t ref_index;
  uint8_t valid;
} MV_MFQE;

static const MV_MFQE kZeroMvMFQE = { { 0, 0 }, 0, 0, 0, 0 };

static const int16_t grid_search_rows[13] = {
  0, 0, 0, 0, 0, 1, 1, 1, -1, -1, -1, 2, -2,
};

static const int16_t grid_search_cols[13] = {
  -2, -1, 0, 1, 2, -1, 0, 1, -1, 0, 1, 0, 0,
};

static INLINE int cmpref(const void *a, const void *b) {
  RefCntBuffer *ref1 = *((RefCntBuffer **)a);
  RefCntBuffer *ref2 = *((RefCntBuffer **)b);
  return ref1->base_qindex - ref2->base_qindex;
}

// Actually apply In-Loop Multi-Frame Quality Enhancement to the tmp buffer,
// using the reference frames. Perform full-pixel motion search on 8x8 blocks,
// then perform finer-grained search to obtain subpel motion vectors. Finally,
// replace the blocks in current frame by interpolation.
void av1_apply_loop_mfqe(Y_BUFFER_CONFIG *tmp, RefCntBuffer *ref_frames[],
                         int block_size, int scale, int high_bd, int bd);

// Apply In-Loop Multi-Frame Quality Enhancement from the decoder side.
void av1_decode_restore_mfqe(AV1_COMMON *cm, int scale, int bsize);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_MFQE_H_
