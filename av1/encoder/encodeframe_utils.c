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

#include "aom_ports/system_state.h"

#include "av1/common/pred_common.h"
#include "av1/common/reconintra.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/encodeframe_utils.h"
#include "av1/encoder/partition_strategy.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/tpl_model.h"

#if !CONFIG_REALTIME_ONLY
void av1_save_context(const MACROBLOCK *x, RD_SEARCH_MACROBLOCK_CONTEXT *ctx,
                      int mi_row, int mi_col, BLOCK_SIZE bsize,
                      const int num_planes) {
  const MACROBLOCKD *xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide =
      block_size_wide[bsize] >> tx_size_wide_log2[0];
  const int num_4x4_blocks_high =
      block_size_high[bsize] >> tx_size_high_log2[0];
  int mi_width = mi_size_wide[bsize];
  int mi_height = mi_size_high[bsize];

  // buffer the above/left context information of the block in search.
  for (p = 0; p < num_planes; ++p) {
    int tx_col = mi_col;
    int tx_row = mi_row & MAX_MIB_MASK;
    memcpy(ctx->a + num_4x4_blocks_wide * p,
           xd->above_context[p] + (tx_col >> xd->plane[p].subsampling_x),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(ctx->l + num_4x4_blocks_high * p,
           xd->left_context[p] + (tx_row >> xd->plane[p].subsampling_y),
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(ctx->sa, xd->above_seg_context + mi_col,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(ctx->sl, xd->left_seg_context + (mi_row & MAX_MIB_MASK),
         sizeof(xd->left_seg_context[0]) * mi_height);
  memcpy(ctx->ta, xd->above_txfm_context,
         sizeof(*xd->above_txfm_context) * mi_width);
  memcpy(ctx->tl, xd->left_txfm_context,
         sizeof(*xd->left_txfm_context) * mi_height);
  ctx->p_ta = xd->above_txfm_context;
  ctx->p_tl = xd->left_txfm_context;
}

void av1_restore_context(const AV1_COMMON *cm, MACROBLOCK *x,
                         const RD_SEARCH_MACROBLOCK_CONTEXT *ctx, int mi_row,
                         int mi_col, BLOCK_SIZE bsize, const int num_planes) {
  MACROBLOCKD *xd = &x->e_mbd;
  int p;
  const int num_4x4_blocks_wide =
      block_size_wide[bsize] >> tx_size_wide_log2[0];
  const int num_4x4_blocks_high =
      block_size_high[bsize] >> tx_size_high_log2[0];
  int mi_width = mi_size_wide[bsize];
  int mi_height = mi_size_high[bsize];
  for (p = 0; p < num_planes; p++) {
    int tx_col = mi_col;
    int tx_row = mi_row & MAX_MIB_MASK;
    memcpy(xd->above_context[p] + (tx_col >> xd->plane[p].subsampling_x),
           ctx->a + num_4x4_blocks_wide * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_wide) >>
               xd->plane[p].subsampling_x);
    memcpy(xd->left_context[p] + (tx_row >> xd->plane[p].subsampling_y),
           ctx->l + num_4x4_blocks_high * p,
           (sizeof(ENTROPY_CONTEXT) * num_4x4_blocks_high) >>
               xd->plane[p].subsampling_y);
  }
  memcpy(xd->above_seg_context + mi_col, ctx->sa,
         sizeof(*xd->above_seg_context) * mi_width);
  memcpy(xd->left_seg_context + (mi_row & MAX_MIB_MASK), ctx->sl,
         sizeof(xd->left_seg_context[0]) * mi_height);
  xd->above_txfm_context = ctx->p_ta;
  xd->left_txfm_context = ctx->p_tl;
  memcpy(xd->above_txfm_context, ctx->ta,
         sizeof(*xd->above_txfm_context) * mi_width);
  memcpy(xd->left_txfm_context, ctx->tl,
         sizeof(*xd->left_txfm_context) * mi_height);

  av1_mark_block_as_not_coded(xd, mi_row, mi_col, bsize,
                              cm->seq_params.sb_size);
}

void av1_set_offsets_without_segment_id(const AV1_COMP *const cpi,
                                        const TileInfo *const tile,
                                        MACROBLOCK *const x, int mi_row,
                                        int mi_col, BLOCK_SIZE bsize,
                                        const CHROMA_REF_INFO *chr_ref_info) {
  const AV1_COMMON *const cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCKD *const xd = &x->e_mbd;
  assert(bsize < BLOCK_SIZES_ALL);
  const int mi_width = mi_size_wide[bsize];
  const int mi_height = mi_size_high[bsize];

  set_mode_info_offsets(cpi, x, xd, mi_row, mi_col);
  set_skip_context(xd, mi_row, mi_col, num_planes, chr_ref_info);
  xd->above_txfm_context = cm->above_txfm_context[tile->tile_row] + mi_col;
  xd->left_txfm_context =
      xd->left_txfm_context_buffer + (mi_row & MAX_MIB_MASK);

  // Set up destination pointers.
  av1_setup_dst_planes(xd->plane, &cm->cur_frame->buf, mi_row, mi_col, 0,
                       num_planes, chr_ref_info);

  // Set up limit values for MV components.
  // Mv beyond the range do not produce new/different prediction block.
  x->mv_limits.row_min =
      -(((mi_row + mi_height) * MI_SIZE) + AOM_INTERP_EXTEND);
  x->mv_limits.col_min = -(((mi_col + mi_width) * MI_SIZE) + AOM_INTERP_EXTEND);
  x->mv_limits.row_max = (cm->mi_rows - mi_row) * MI_SIZE + AOM_INTERP_EXTEND;
  x->mv_limits.col_max = (cm->mi_cols - mi_col) * MI_SIZE + AOM_INTERP_EXTEND;

  set_plane_n4(xd, bsize, num_planes, chr_ref_info);

  // Set up distance of MB to edge of frame in 1/8th pel units.
#if !CONFIG_EXT_RECUR_PARTITIONS
  assert(!(mi_col & (mi_width - 1)) && !(mi_row & (mi_height - 1)));
#endif  // !CONFIG_EXT_RECUR_PARTITIONS
  set_mi_row_col(xd, tile, mi_row, mi_height, mi_col, mi_width, cm->mi_rows,
                 cm->mi_cols, chr_ref_info);

  // Set up source buffers.
  av1_setup_src_planes(x, cpi->source, mi_row, mi_col, num_planes,
                       chr_ref_info);

  // required by av1_append_sub8x8_mvs_for_idx() and av1_find_best_ref_mvs()
  xd->tile = *tile;

  xd->cfl.mi_row = mi_row;
  xd->cfl.mi_col = mi_col;
}

void av1_enc_set_offsets(const AV1_COMP *const cpi, const TileInfo *const tile,
                         MACROBLOCK *const x, int mi_row, int mi_col,
                         BLOCK_SIZE bsize, CHROMA_REF_INFO *chr_ref_info) {
  const AV1_COMMON *const cm = &cpi->common;
  const struct segmentation *const seg = &cm->seg;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi;

  av1_set_offsets_without_segment_id(cpi, tile, x, mi_row, mi_col, bsize,
                                     chr_ref_info);

  // Setup segment ID.
  mbmi = xd->mi[0];
  mbmi->segment_id = 0;
  if (seg->enabled) {
    if (seg->enabled && !cpi->vaq_refresh) {
      const uint8_t *const map =
          seg->update_map ? cpi->segmentation_map : cm->last_frame_seg_map;
      mbmi->segment_id =
          map ? get_segment_id(cm, map, bsize, mi_row, mi_col) : 0;
    }
    av1_init_plane_quantizers(cpi, x, mbmi->segment_id);
  }
}

#endif  // !CONFIG_REALTIME_ONLY
