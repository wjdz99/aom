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

#include "av1/common/blockd.h"

PREDICTION_MODE av1_left_block_mode(const MODE_INFO *cur_mi,
                                    const MODE_INFO *left_mi, int b) {
  if (b == 0 || b == 2) {
    if (!left_mi || is_inter_block(&left_mi->mbmi)) return DC_PRED;

    return get_y_mode(left_mi, b + 1);
  } else {
    assert(b == 1 || b == 3);
    return cur_mi->bmi[b - 1].as_mode;
  }
}

PREDICTION_MODE av1_above_block_mode(const MODE_INFO *cur_mi,
                                     const MODE_INFO *above_mi, int b) {
  if (b == 0 || b == 1) {
    if (!above_mi || is_inter_block(&above_mi->mbmi)) return DC_PRED;

    return get_y_mode(above_mi, b + 2);
  } else {
    assert(b == 2 || b == 3);
    return cur_mi->bmi[b - 2].as_mode;
  }
}

#if CONFIG_COEF_INTERLEAVE
void av1_foreach_transformed_block_interleave(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize,
    foreach_transformed_block_visitor visit, void *arg) {
  const struct macroblockd_plane *const pd_y = &xd->plane[0];
  const struct macroblockd_plane *const pd_c = &xd->plane[1];
  const MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;

  const TX_SIZE tx_log2_y = mbmi->tx_size;
  const TX_SIZE tx_log2_c = get_uv_tx_size(mbmi, pd_c);
  const int tx_sz_y = (1 << tx_log2_y);
  const int tx_sz_c = (1 << tx_log2_c);

  const BLOCK_SIZE plane_bsize_y = get_plane_block_size(bsize, pd_y);
  const BLOCK_SIZE plane_bsize_c = get_plane_block_size(bsize, pd_c);

  const int num_4x4_w_y = num_4x4_blocks_wide_lookup[plane_bsize_y];
  const int num_4x4_w_c = num_4x4_blocks_wide_lookup[plane_bsize_c];
  const int num_4x4_h_y = num_4x4_blocks_high_lookup[plane_bsize_y];
  const int num_4x4_h_c = num_4x4_blocks_high_lookup[plane_bsize_c];

  const int step_y = 1 << (tx_log2_y << 1);
  const int step_c = 1 << (tx_log2_c << 1);

  const int max_4x4_w_y =
      get_max_4x4_size(num_4x4_w_y, xd->mb_to_right_edge, pd_y->subsampling_x);
  const int max_4x4_h_y =
      get_max_4x4_size(num_4x4_h_y, xd->mb_to_bottom_edge, pd_y->subsampling_y);

  const int extra_step_y = ((num_4x4_w_y - max_4x4_w_y) >> tx_log2_y) * step_y;

  const int max_4x4_w_c =
      get_max_4x4_size(num_4x4_w_c, xd->mb_to_right_edge, pd_c->subsampling_x);
  const int max_4x4_h_c =
      get_max_4x4_size(num_4x4_h_c, xd->mb_to_bottom_edge, pd_c->subsampling_y);

  const int extra_step_c = ((num_4x4_w_c - max_4x4_w_c) >> tx_log2_c) * step_c;

  // The max_4x4_w/h may be smaller than tx_sz under some corner cases,
  // i.e. when the SB is splitted by tile boundaries.
  const int tu_num_w_y = (max_4x4_w_y + tx_sz_y - 1) / tx_sz_y;
  const int tu_num_h_y = (max_4x4_h_y + tx_sz_y - 1) / tx_sz_y;
  const int tu_num_w_c = (max_4x4_w_c + tx_sz_c - 1) / tx_sz_c;
  const int tu_num_h_c = (max_4x4_h_c + tx_sz_c - 1) / tx_sz_c;
  const int tu_num_y = tu_num_w_y * tu_num_h_y;
  const int tu_num_c = tu_num_w_c * tu_num_h_c;

  int tu_idx_c = 0;
  int offset_y, row_y, col_y;
  int offset_c, row_c, col_c;

  for (row_y = 0; row_y < tu_num_h_y; row_y++) {
    for (col_y = 0; col_y < tu_num_w_y; col_y++) {
      // luma
      offset_y = (row_y * tu_num_w_y + col_y) * step_y + row_y * extra_step_y;
      visit(0, offset_y, row_y * tx_sz_y, col_y * tx_sz_y, plane_bsize_y,
            tx_log2_y, arg);
      // chroma
      if (tu_idx_c < tu_num_c) {
        row_c = (tu_idx_c / tu_num_w_c) * tx_sz_c;
        col_c = (tu_idx_c % tu_num_w_c) * tx_sz_c;
        offset_c = tu_idx_c * step_c + (tu_idx_c / tu_num_w_c) * extra_step_c;
        visit(1, offset_c, row_c, col_c, plane_bsize_c, tx_log2_c, arg);
        visit(2, offset_c, row_c, col_c, plane_bsize_c, tx_log2_c, arg);
        tu_idx_c++;
      }
    }
  }

  // In 422 case, it's possilbe that Chroma has more TUs than Luma
  while (tu_idx_c < tu_num_c) {
    row_c = (tu_idx_c / tu_num_w_c) * tx_sz_c;
    col_c = (tu_idx_c % tu_num_w_c) * tx_sz_c;
    offset_c = tu_idx_c * step_c + row_c * extra_step_c;
    visit(1, offset_c, row_c, col_c, plane_bsize_c, tx_log2_c, arg);
    visit(2, offset_c, row_c, col_c, plane_bsize_c, tx_log2_c, arg);
    tu_idx_c++;
  }
}
#endif // CONFIG_COEF_INTERLEAVE

void av1_foreach_transformed_block_in_plane(
    const MACROBLOCKD *const xd, BLOCK_SIZE bsize, int plane,
    foreach_transformed_block_visitor visit, void *arg) {
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  // block and transform sizes, in number of 4x4 blocks log 2 ("*_b")
  // 4x4=0, 8x8=2, 16x16=4, 32x32=6, 64x64=8
  // transform size varies per plane, look it up in a common way.
  const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int num_4x4_w = num_4x4_blocks_wide_lookup[plane_bsize];
  const int num_4x4_h = num_4x4_blocks_high_lookup[plane_bsize];
  const int step = 1 << (tx_size_1d_in_unit_log2[tx_size] * 2);
  int i = 0, r, c;

  // If mb_to_right_edge is < 0 we are in a situation in which
  // the current block size extends into the UMV and we won't
  // visit the sub blocks that are wholly within the UMV.
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >>
                                                       (5 + pd->subsampling_x));
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0
                       ? 0
                       : xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));
  const int extra_step =
      ((num_4x4_w - max_blocks_wide) >> tx_size_1d_in_unit_log2[tx_size]) *
      step;

  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (r = 0; r < max_blocks_high; r += tx_size_1d_in_unit[tx_size]) {
    // Skip visiting the sub blocks that are wholly within the UMV.
    for (c = 0; c < max_blocks_wide; c += tx_size_1d_in_unit[tx_size]) {
      visit(plane, i, r, c, plane_bsize, tx_size, arg);
      i += step;
    }
    i += extra_step;
  }
}

void av1_foreach_transformed_block(const MACROBLOCKD *const xd,
                                   BLOCK_SIZE bsize,
                                   foreach_transformed_block_visitor visit,
                                   void *arg) {
  int plane;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane)
    av1_foreach_transformed_block_in_plane(xd, bsize, plane, visit, arg);
}

void av1_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      TX_SIZE tx_size, int has_eob, int aoff, int loff) {
  ENTROPY_CONTEXT *const a = pd->above_context + aoff;
  ENTROPY_CONTEXT *const l = pd->left_context + loff;
  const int tx_size_in_blocks = tx_size_1d_in_unit[tx_size];

  // above
  if (has_eob && xd->mb_to_right_edge < 0) {
    int i;
    const int blocks_wide =
        pd->n4_w + (xd->mb_to_right_edge >> (5 + pd->subsampling_x));
    int above_contexts = tx_size_in_blocks;
    if (above_contexts + aoff > blocks_wide)
      above_contexts = blocks_wide - aoff;

    for (i = 0; i < above_contexts; ++i) a[i] = has_eob;
    for (i = above_contexts; i < tx_size_in_blocks; ++i) a[i] = 0;
  } else {
    memset(a, has_eob, sizeof(ENTROPY_CONTEXT) * tx_size_in_blocks);
  }

  // left
  if (has_eob && xd->mb_to_bottom_edge < 0) {
    int i;
    const int blocks_high =
        pd->n4_h + (xd->mb_to_bottom_edge >> (5 + pd->subsampling_y));
    int left_contexts = tx_size_in_blocks;
    if (left_contexts + loff > blocks_high) left_contexts = blocks_high - loff;

    for (i = 0; i < left_contexts; ++i) l[i] = has_eob;
    for (i = left_contexts; i < tx_size_in_blocks; ++i) l[i] = 0;
  } else {
    memset(l, has_eob, sizeof(ENTROPY_CONTEXT) * tx_size_in_blocks);
  }
}

void av1_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y) {
  int i;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    xd->plane[i].plane_type = i ? PLANE_TYPE_UV : PLANE_TYPE_Y;
    xd->plane[i].subsampling_x = i ? ss_x : 0;
    xd->plane[i].subsampling_y = i ? ss_y : 0;
  }
}
