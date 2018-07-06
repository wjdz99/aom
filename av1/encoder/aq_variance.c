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

#include <math.h>

#include "aom_ports/mem.h"

#include "av1/encoder/aq_variance.h"

#include "av1/common/seg_common.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/segmentation.h"
#include "aom_ports/system_state.h"

static const double rate_ratio[MAX_SEGMENTS] = { 2.2, 1.7, 1.3, 1.0,
                                                 0.9, .8,  .7,  .6 };

#define SEGMENT_ID(i) i;

DECLARE_ALIGNED(16, static const uint8_t, av1_all_zeros[MAX_SB_SIZE]) = { 0 };
DECLARE_ALIGNED(16, static const uint16_t,
                av1_highbd_all_zeros[MAX_SB_SIZE]) = { 0 };

unsigned int av1_vaq_segment_id(int energy) { return SEGMENT_ID(energy); }

void av1_vaq_frame_setup(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;
  struct segmentation *seg = &cm->seg;
  int i;

  int resolution_change =
      cm->prev_frame && (cm->width != cm->prev_frame->width ||
                         cm->height != cm->prev_frame->height);
  int avg_energy = (int)(cpi->twopass.mb_av_energy - 2);
  double avg_ratio;
  if (avg_energy > 7) avg_energy = 7;
  if (avg_energy < 0) avg_energy = 0;
  avg_ratio = rate_ratio[avg_energy];

  if (resolution_change) {
    memset(cpi->segmentation_map, 0, cm->mi_rows * cm->mi_cols);
    av1_clearall_segfeatures(seg);
    aom_clear_system_state();
    av1_disable_segmentation(seg);
    return;
  }
  if (frame_is_intra_only(cm) || cm->error_resilient_mode ||
      cpi->refresh_alt_ref_frame ||
      (cpi->refresh_golden_frame && !cpi->rc.is_src_frame_alt_ref)) {
    cpi->vaq_refresh = 1;

    av1_enable_segmentation(seg);
    av1_clearall_segfeatures(seg);

    aom_clear_system_state();

    for (i = 0; i < MAX_SEGMENTS; ++i) {
      // Set up avg segment id to be 1.0 and adjust the other segments around
      // it.
      int qindex_delta =
          av1_compute_qdelta_by_rate(&cpi->rc, cm->frame_type, cm->base_qindex,
                                     rate_ratio[i] / avg_ratio, cm->bit_depth);

      // We don't allow qindex 0 in a segment if the base value is not 0.
      // Q index 0 (lossless) implies 4x4 encoding only and in AQ mode a segment
      // Q delta is sometimes applied without going back around the rd loop.
      // This could lead to an illegal combination of partition size and q.
      if ((cm->base_qindex != 0) && ((cm->base_qindex + qindex_delta) == 0)) {
        qindex_delta = -cm->base_qindex + 1;
      }

      av1_set_segdata(seg, i, SEG_LVL_ALT_Q, qindex_delta);
      av1_enable_segfeature(seg, i, SEG_LVL_ALT_Q);
    }
  }
}

static unsigned int block_lvariance(const AV1_COMP *const cpi, MACROBLOCK *x,
                                    BLOCK_SIZE bs) {
  // This functions returns a score for the blocks local variance as calculated
  // by: sum of the log of the (4x4 variances) of each subblock to the current
  // block (x,bs)
  // * 32 / number of pixels in the block_size.
  // This is used for segmentation because to avoid situations in which a large
  // block with a gentle gradient gets marked high variance even though each
  // subblock has a low variance.   This allows us to assign the same segment
  // number for the same sorts of area regardless of how the partitioning goes.

  MACROBLOCKD *xd = &x->e_mbd;
  double var = 0;
  unsigned int sse;
  int i, j;

  int right_overflow =
      (xd->mb_to_right_edge < 0) ? ((-xd->mb_to_right_edge) >> 3) : 0;
  int bottom_overflow =
      (xd->mb_to_bottom_edge < 0) ? ((-xd->mb_to_bottom_edge) >> 3) : 0;

  const int bw = MI_SIZE * mi_size_wide[bs] - right_overflow;
  const int bh = MI_SIZE * mi_size_high[bs] - bottom_overflow;

  for (i = 0; i < bh; i += 4) {
    for (j = 0; j < bw; j += 4) {
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        var +=
            log(1.0 + cpi->fn_ptr[BLOCK_4X4].vf(
                          x->plane[0].src.buf + i * x->plane[0].src.stride + j,
                          x->plane[0].src.stride,
                          CONVERT_TO_BYTEPTR(av1_highbd_all_zeros), 0, &sse) /
                          16);
      } else {
        var +=
            log(1.0 + cpi->fn_ptr[BLOCK_4X4].vf(
                          x->plane[0].src.buf + i * x->plane[0].src.stride + j,
                          x->plane[0].src.stride, av1_all_zeros, 0, &sse) /
                          16);
      }
    }
  }
  // Use average of 4x4 log variance. The range for 8 bit 0 - 9.704121561.
  var /= (bw / 4 * bh / 4);
  if (var > 7) var = 7;

  // Convert to integer.
  return (unsigned int)var;
}

double av1_log_block_var(const AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bs) {
  unsigned int var = block_lvariance(cpi, x, bs);
  aom_clear_system_state();
  return var;
}

#define DEFAULT_E_MIDPOINT 10.0
int av1_block_energy(const AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bs) {
  aom_clear_system_state();
  return block_lvariance(cpi, x, bs);
}
