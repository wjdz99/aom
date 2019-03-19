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

#ifndef AV1_ENCODER_PARTITION_STRATEGY_H_
#define AV1_ENCODER_PARTITION_STRATEGY_H_

#include "av1/encoder/encodeframe.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodemb.h"

// Performs a motion search in SIMPLE_TRANSLATION mode using reference frame
// ref. Note that this sets the offset of mbmi, so we will need to reset it
// after calling this function.
void av1_simple_motion_search(AV1_COMP *const cpi, MACROBLOCK *x, int mi_row,
                              int mi_col, BLOCK_SIZE bsize, int ref,
                              MV ref_mv_full, int num_planes, int use_subpixel);

// Performs a simple motion search to calculate the sse and var of the residue
void av1_simple_motion_sse_var(AV1_COMP *cpi, MACROBLOCK *x, int mi_row,
                               int mi_col, BLOCK_SIZE bsize,
                               const MV ref_mv_full, int use_subpixel,
                               unsigned int *sse, unsigned int *var);

// Performs a simple_motion_search with a single reference frame and extract
// the variance of residues. Then use the features to determine whether we want
// to go straight to splitting without trying PARTITION_NONE
void av1_simple_motion_search_based_split(
    AV1_COMP *const cpi, MACROBLOCK *x, int mi_row, int mi_col,
    BLOCK_SIZE bsize, int *partition_none_allowed, int *partition_horz_allowed,
    int *partition_vert_allowed, int *do_rectangular_split);

// Given a list of ref frames in refs, performs simple_motion_search on each of
// the refs and returns the ref with the smallest sse. Returns -1 if none of the
// ref in the list is available. Also stores the best sse and var in best_sse,
// best_var, respectively. If save_mv_code is -1, don't update mv_ref_fulls in
// pc_tree. If save_mv_code is between 0 and 3, update mv_ref_fulls under
// pc_tree->split[i]. If save_mv_code is 4, update mv_ref_fulls under pc_tree.
static int simple_motion_search_get_best_ref(
    AV1_COMP *const cpi, MACROBLOCK *x, PC_TREE *pc_tree, int mi_row,
    int mi_col, BLOCK_SIZE bsize, const int *const refs, int num_refs,
    int use_subpixel, int save_mv_code, unsigned int *best_sse,
    unsigned int *best_var) {
  // TODO(chiyotsai@google.com): The calculation of variance currently uses
  // bsize, so we might take area outside of the image into account. We need to
  // modify the SIMD functions to fix this later.
  const AV1_COMMON *const cm = &cpi->common;
  int best_ref = -1;

  if (mi_col >= cm->mi_cols || mi_row >= cm->mi_rows) {
    // If the whole block is outside of the image, set the var and sse to 0.
    *best_var = 0;
    *best_sse = 0;

    return best_ref;
  }

  // Otherwise do loop through the reference frames and find the one with the
  // minimum SSE
  const MACROBLOCKD *xd = &x->e_mbd;
  const MV *mv_ref_fulls = pc_tree->mv_ref_fulls;

  const int num_planes = 1;

  *best_sse = INT_MAX;

  for (int ref_idx = 0; ref_idx < num_refs; ref_idx++) {
    const int ref = refs[ref_idx];

    if (cpi->ref_frame_flags & av1_ref_frame_flag_list[ref]) {
      unsigned int curr_sse = 0, curr_var = 0;
      av1_simple_motion_search(cpi, x, mi_row, mi_col, bsize, ref,
                               mv_ref_fulls[ref], num_planes, use_subpixel);
      curr_var = cpi->fn_ptr[bsize].vf(
          x->plane[0].src.buf, x->plane[0].src.stride, xd->plane[0].dst.buf,
          xd->plane[0].dst.stride, &curr_sse);
      if (curr_sse < *best_sse) {
        *best_sse = curr_sse;
        *best_var = curr_var;
        best_ref = ref;
      }

      const int new_mv_row = x->best_mv.as_mv.row / 8;
      const int new_mv_col = x->best_mv.as_mv.col / 8;
      if (save_mv_code == 4) {
        pc_tree->mv_ref_fulls[ref].row = new_mv_row;
        pc_tree->mv_ref_fulls[ref].col = new_mv_col;
      } else if (save_mv_code >= 0 && save_mv_code < 4) {
        // Propagate the new motion vectors to a lower level
        pc_tree->split[save_mv_code]->mv_ref_fulls[ref].row = new_mv_row;
        pc_tree->split[save_mv_code]->mv_ref_fulls[ref].col = new_mv_col;
      } else {
        assert(save_mv_code == -1 &&
               "Unknown code in simple_motion_search_get_best_ref.");
      }
    }
  }

  return best_ref;
}
#endif  // AV1_ENCODER_PARTITION_STRATEGY_H_
