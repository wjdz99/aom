/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>

#include "av1/encoder/encoder.h"
#include "av1/encoder/ethread.h"
#include "av1/encoder/av1_multi_thread.h"

void av1_row_mt_mem_alloc(AV1_COMP *cpi, int max_sb_rows) {
  struct AV1Common *cm = &cpi->common;
  int tile_cols = cm->tiles.cols;
  int tile_rows = cm->tiles.rows;
  EncRowMT *enc_row_mt = &cpi->enc_row_mt;
  int tile_row, tile_col;

  CHECK_MEM_ERROR(cm, enc_row_mt->row_mt_sync,
                  aom_memalign(32, tile_cols * tile_rows *
                                       sizeof(*enc_row_mt->row_mt_sync)));

  // Allocate memory for row based multi-threading
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      int tile_index = tile_row * tile_cols + tile_col;
      TileDataEnc *this_tile = &cpi->tile_data[tile_index];
      AV1RowMTSync *row_mt_sync = &enc_row_mt->row_mt_sync[tile_index];
      av1_row_mt_sync_mem_alloc(row_mt_sync, cm, max_sb_rows);
      if (cpi->oxcf.cdf_update_mode)
        CHECK_MEM_ERROR(
            cm, row_mt_sync->row_ctx,
            (FRAME_CONTEXT *)aom_memalign(
                16,
                AOMMAX(1, (av1_get_sb_cols_in_tile(cm, this_tile->tile_info) -
                           1)) *
                    sizeof(*row_mt_sync->row_ctx)));
    }
  }

  enc_row_mt->allocated_tile_cols = tile_cols;
  enc_row_mt->allocated_tile_rows = tile_rows;
  enc_row_mt->allocated_sb_rows = max_sb_rows;
}

void av1_row_mt_mem_dealloc(AV1_COMP *cpi) {
  EncRowMT *enc_row_mt = &cpi->enc_row_mt;
  int tile_cols = enc_row_mt->allocated_tile_cols;
  int tile_rows = enc_row_mt->allocated_tile_rows;
  int tile_col, tile_row;

  // Free row based multi-threading sync memory
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      AV1RowMTSync *row_mt_sync =
          &enc_row_mt->row_mt_sync[tile_row * tile_cols + tile_col];
      av1_row_mt_sync_mem_dealloc(row_mt_sync);
      if (cpi->oxcf.cdf_update_mode) aom_free(row_mt_sync->row_ctx);
    }
  }
  enc_row_mt->allocated_sb_rows = 0;
  enc_row_mt->allocated_tile_cols = 0;
  enc_row_mt->allocated_tile_rows = 0;
}
