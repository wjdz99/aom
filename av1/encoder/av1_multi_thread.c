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

void av1_row_mt_mem_alloc(AV1_COMP *const cpi, int max_sb_rows) {
  struct AV1Common *cm = &cpi->common;
  AV1EncRowMT *const enc_row_mt = &cpi->enc_row_mt;
  const int tile_cols = cm->tiles.cols;
  const int tile_rows = cm->tiles.rows;
  int tile_col, tile_row;

  CHECK_MEM_ERROR(
      cm, enc_row_mt->sync,
      aom_memalign(32, tile_cols * tile_rows * sizeof(*enc_row_mt->sync)));

  // Allocate memory for row based multi-threading
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      int tile_index = tile_row * tile_cols + tile_col;
      TileInfo tile_info = cpi->tile_data[tile_index].tile_info;
      SyncData *const sync = &enc_row_mt->sync[tile_index];

      av1_row_mt_sync_mem_alloc(sync, cm, max_sb_rows);

      if (cpi->oxcf.cdf_update_mode) {
        const int sb_cols_in_tile = av1_get_sb_cols_in_tile(cm, tile_info);
        const int num_row_ctx = AOMMAX(1, (sb_cols_in_tile - 1));
        CHECK_MEM_ERROR(cm, sync->row_ctx,
                        (FRAME_CONTEXT *)aom_memalign(
                            16, num_row_ctx * sizeof(*sync->row_ctx)));
      }
    }
  }
  enc_row_mt->allocated_tile_cols = tile_cols;
  enc_row_mt->allocated_tile_rows = tile_rows;
  enc_row_mt->allocated_sb_rows = max_sb_rows;
}

void av1_row_mt_mem_dealloc(AV1EncRowMT *const enc_row_mt,
                            uint8_t cdf_update_mode) {
  const int tile_cols = enc_row_mt->allocated_tile_cols;
  const int tile_rows = enc_row_mt->allocated_tile_rows;
  int tile_col, tile_row;

  // Free row based multi-threading sync memory
  for (tile_row = 0; tile_row < tile_rows; tile_row++) {
    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
      SyncData *const sync = &enc_row_mt->sync[tile_row * tile_cols + tile_col];

      av1_row_mt_sync_mem_dealloc(sync);

      if (cdf_update_mode) aom_free(sync->row_ctx);
    }
  }
  enc_row_mt->allocated_sb_rows = 0;
  enc_row_mt->allocated_tile_cols = 0;
  enc_row_mt->allocated_tile_rows = 0;
}
