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

#include "av1/encoder/rdopt_utils.h"
#include "config/aom_dsp_rtcd.h"

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

int64_t calculate_sse(MACROBLOCKD *const xd, const struct macroblock_plane *p,
                      struct macroblockd_plane *pd, const int bw,
                      const int bh) {
  int64_t sse = 0;
  const int shift = xd->bd - 8;
#if CONFIG_AV1_HIGHBITDEPTH
  if (is_cur_buf_hbd(xd)) {
    sse = aom_highbd_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
                         bw, bh);
  } else {
    sse =
        aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
  }
#else
  sse = aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
#endif
  sse = ROUND_POWER_OF_TWO(sse, shift * 2);
  return sse;
}
