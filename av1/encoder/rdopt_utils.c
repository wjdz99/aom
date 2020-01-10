/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/encoder/block.h"
#include "av1/encoder/rdopt_data_defs.h"

int inter_mode_data_block_idx(BLOCK_SIZE bsize) {
  if (bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
      bsize == BLOCK_4X16 || bsize == BLOCK_16X4) {
    return -1;
  }
  return 1;
}

THR_MODES get_prediction_mode_idx(PREDICTION_MODE this_mode,
                                  MV_REFERENCE_FRAME ref_frame,
                                  MV_REFERENCE_FRAME second_ref_frame) {
  if (this_mode < INTRA_MODE_END) {
    assert(ref_frame == INTRA_FRAME);
    assert(second_ref_frame == NONE_FRAME);
    return intra_to_mode_idx[this_mode - INTRA_MODE_START];
  }
  if (this_mode >= SINGLE_INTER_MODE_START &&
      this_mode < SINGLE_INTER_MODE_END) {
    assert((ref_frame > INTRA_FRAME) && (ref_frame <= ALTREF_FRAME));
    return single_inter_to_mode_idx[this_mode - SINGLE_INTER_MODE_START]
                                   [ref_frame];
  }
  if (this_mode >= COMP_INTER_MODE_START && this_mode < COMP_INTER_MODE_END) {
    assert((ref_frame > INTRA_FRAME) && (ref_frame <= ALTREF_FRAME));
    assert((second_ref_frame > INTRA_FRAME) &&
           (second_ref_frame <= ALTREF_FRAME));
    return comp_inter_to_mode_idx[this_mode - COMP_INTER_MODE_START][ref_frame]
                                 [second_ref_frame];
  }
  assert(0);
  return THR_INVALID;
}
