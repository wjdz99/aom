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

#include "av1/common/cfl.h"
#include "av1/common/idct.h"
#include "av1/common/pred_common.h"
#include "av1/common/reconintra.h"
#include "av1/encoder/hybrid_fwd_txfm.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/transform_search.h"
#include "av1/encoder/tx_prune_model_weights.h"

struct rdcost_block_args {
  const AV1_COMP *cpi;
  MACROBLOCK *x;
  ENTROPY_CONTEXT t_above[MAX_MIB_SIZE];
  ENTROPY_CONTEXT t_left[MAX_MIB_SIZE];
  RD_STATS rd_stats;
  int64_t this_rd;
  int64_t best_rd;
  int exit_early;
  int incomplete_exit;
  int use_fast_coef_costing;
  FAST_TX_SEARCH_MODE ftxs_mode;
  int skip_trellis;
};

typedef struct {
  int64_t rd;
  int txb_entropy_ctx;
  TX_TYPE tx_type;
} TxCandidateInfo;

typedef struct {
  int leaf;
  int8_t children[4];
} RD_RECORD_IDX_NODE;

// origin_threshold * 128 / 100
static const uint32_t skip_pred_threshold[3][BLOCK_SIZES_ALL] = {
  {
      64, 64, 64, 70, 60, 60, 68, 68, 68, 68, 68,
      68, 68, 68, 68, 68, 64, 64, 70, 70, 68, 68,
  },
  {
      88, 88, 88, 86, 87, 87, 68, 68, 68, 68, 68,
      68, 68, 68, 68, 68, 88, 88, 86, 86, 68, 68,
  },
  {
      90, 93, 93, 90, 93, 93, 74, 74, 74, 74, 74,
      74, 74, 74, 74, 74, 90, 90, 90, 90, 74, 74,
  },
};

// lookup table for predict_skip_flag
// int max_tx_size = max_txsize_rect_lookup[bsize];
// if (tx_size_high[max_tx_size] > 16 || tx_size_wide[max_tx_size] > 16)
//   max_tx_size = AOMMIN(max_txsize_lookup[bsize], TX_16X16);
static const TX_SIZE max_predict_sf_tx_size[BLOCK_SIZES_ALL] = {
  TX_4X4,   TX_4X8,   TX_8X4,   TX_8X8,   TX_8X16,  TX_16X8,
  TX_16X16, TX_16X16, TX_16X16, TX_16X16, TX_16X16, TX_16X16,
  TX_16X16, TX_16X16, TX_16X16, TX_16X16, TX_4X16,  TX_16X4,
  TX_8X8,   TX_8X8,   TX_16X16, TX_16X16,
};

static const RD_RECORD_IDX_NODE rd_record_tree_8x8[] = {
  { 1, { 0 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_8x16[] = {
  { 0, { 1, 2, -1, -1 } },
  { 1, { 0, 0, 0, 0 } },
  { 1, { 0, 0, 0, 0 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_16x8[] = {
  { 0, { 1, 2, -1, -1 } },
  { 1, { 0 } },
  { 1, { 0 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_16x16[] = {
  { 0, { 1, 2, 3, 4 } }, { 1, { 0 } }, { 1, { 0 } }, { 1, { 0 } }, { 1, { 0 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_1_2[] = {
  { 0, { 1, 2, -1, -1 } },
  { 0, { 3, 4, 5, 6 } },
  { 0, { 7, 8, 9, 10 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_2_1[] = {
  { 0, { 1, 2, -1, -1 } },
  { 0, { 3, 4, 7, 8 } },
  { 0, { 5, 6, 9, 10 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_sqr[] = {
  { 0, { 1, 2, 3, 4 } },     { 0, { 5, 6, 9, 10 } },    { 0, { 7, 8, 11, 12 } },
  { 0, { 13, 14, 17, 18 } }, { 0, { 15, 16, 19, 20 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_64x128[] = {
  { 0, { 2, 3, 4, 5 } },     { 0, { 6, 7, 8, 9 } },
  { 0, { 10, 11, 14, 15 } }, { 0, { 12, 13, 16, 17 } },
  { 0, { 18, 19, 22, 23 } }, { 0, { 20, 21, 24, 25 } },
  { 0, { 26, 27, 30, 31 } }, { 0, { 28, 29, 32, 33 } },
  { 0, { 34, 35, 38, 39 } }, { 0, { 36, 37, 40, 41 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_128x64[] = {
  { 0, { 2, 3, 6, 7 } },     { 0, { 4, 5, 8, 9 } },
  { 0, { 10, 11, 18, 19 } }, { 0, { 12, 13, 20, 21 } },
  { 0, { 14, 15, 22, 23 } }, { 0, { 16, 17, 24, 25 } },
  { 0, { 26, 27, 34, 35 } }, { 0, { 28, 29, 36, 37 } },
  { 0, { 30, 31, 38, 39 } }, { 0, { 32, 33, 40, 41 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_128x128[] = {
  { 0, { 4, 5, 8, 9 } },     { 0, { 6, 7, 10, 11 } },
  { 0, { 12, 13, 16, 17 } }, { 0, { 14, 15, 18, 19 } },
  { 0, { 20, 21, 28, 29 } }, { 0, { 22, 23, 30, 31 } },
  { 0, { 24, 25, 32, 33 } }, { 0, { 26, 27, 34, 35 } },
  { 0, { 36, 37, 44, 45 } }, { 0, { 38, 39, 46, 47 } },
  { 0, { 40, 41, 48, 49 } }, { 0, { 42, 43, 50, 51 } },
  { 0, { 52, 53, 60, 61 } }, { 0, { 54, 55, 62, 63 } },
  { 0, { 56, 57, 64, 65 } }, { 0, { 58, 59, 66, 67 } },
  { 0, { 68, 69, 76, 77 } }, { 0, { 70, 71, 78, 79 } },
  { 0, { 72, 73, 80, 81 } }, { 0, { 74, 75, 82, 83 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_1_4[] = {
  { 0, { 1, -1, 2, -1 } },
  { 0, { 3, 4, -1, -1 } },
  { 0, { 5, 6, -1, -1 } },
};

static const RD_RECORD_IDX_NODE rd_record_tree_4_1[] = {
  { 0, { 1, 2, -1, -1 } },
  { 0, { 3, 4, -1, -1 } },
  { 0, { 5, 6, -1, -1 } },
};

static const RD_RECORD_IDX_NODE *rd_record_tree[BLOCK_SIZES_ALL] = {
  NULL,                    // BLOCK_4X4
  NULL,                    // BLOCK_4X8
  NULL,                    // BLOCK_8X4
  rd_record_tree_8x8,      // BLOCK_8X8
  rd_record_tree_8x16,     // BLOCK_8X16
  rd_record_tree_16x8,     // BLOCK_16X8
  rd_record_tree_16x16,    // BLOCK_16X16
  rd_record_tree_1_2,      // BLOCK_16X32
  rd_record_tree_2_1,      // BLOCK_32X16
  rd_record_tree_sqr,      // BLOCK_32X32
  rd_record_tree_1_2,      // BLOCK_32X64
  rd_record_tree_2_1,      // BLOCK_64X32
  rd_record_tree_sqr,      // BLOCK_64X64
  rd_record_tree_64x128,   // BLOCK_64X128
  rd_record_tree_128x64,   // BLOCK_128X64
  rd_record_tree_128x128,  // BLOCK_128X128
  NULL,                    // BLOCK_4X16
  NULL,                    // BLOCK_16X4
  rd_record_tree_1_4,      // BLOCK_8X32
  rd_record_tree_4_1,      // BLOCK_32X8
  rd_record_tree_1_4,      // BLOCK_16X64
  rd_record_tree_4_1,      // BLOCK_64X16
};

static const int rd_record_tree_size[BLOCK_SIZES_ALL] = {
  0,                                                            // BLOCK_4X4
  0,                                                            // BLOCK_4X8
  0,                                                            // BLOCK_8X4
  sizeof(rd_record_tree_8x8) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_8X8
  sizeof(rd_record_tree_8x16) / sizeof(RD_RECORD_IDX_NODE),     // BLOCK_8X16
  sizeof(rd_record_tree_16x8) / sizeof(RD_RECORD_IDX_NODE),     // BLOCK_16X8
  sizeof(rd_record_tree_16x16) / sizeof(RD_RECORD_IDX_NODE),    // BLOCK_16X16
  sizeof(rd_record_tree_1_2) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_16X32
  sizeof(rd_record_tree_2_1) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_32X16
  sizeof(rd_record_tree_sqr) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_32X32
  sizeof(rd_record_tree_1_2) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_32X64
  sizeof(rd_record_tree_2_1) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_64X32
  sizeof(rd_record_tree_sqr) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_64X64
  sizeof(rd_record_tree_64x128) / sizeof(RD_RECORD_IDX_NODE),   // BLOCK_64X128
  sizeof(rd_record_tree_128x64) / sizeof(RD_RECORD_IDX_NODE),   // BLOCK_128X64
  sizeof(rd_record_tree_128x128) / sizeof(RD_RECORD_IDX_NODE),  // BLOCK_128X128
  0,                                                            // BLOCK_4X16
  0,                                                            // BLOCK_16X4
  sizeof(rd_record_tree_1_4) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_8X32
  sizeof(rd_record_tree_4_1) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_32X8
  sizeof(rd_record_tree_1_4) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_16X64
  sizeof(rd_record_tree_4_1) / sizeof(RD_RECORD_IDX_NODE),      // BLOCK_64X16
};

static INLINE void init_rd_record_tree(TXB_RD_INFO_NODE *tree,
                                       BLOCK_SIZE bsize) {
  const RD_RECORD_IDX_NODE *rd_record = rd_record_tree[bsize];
  const int size = rd_record_tree_size[bsize];
  for (int i = 0; i < size; ++i) {
    if (rd_record[i].leaf) {
      av1_zero(tree[i].children);
    } else {
      for (int j = 0; j < 4; ++j) {
        const int8_t idx = rd_record[i].children[j];
        tree[i].children[j] = idx > 0 ? &tree[idx] : NULL;
      }
    }
  }
}

static int find_tx_size_rd_info(TXB_RD_RECORD *cur_record,
                                const uint32_t hash) {
  // Linear search through the circular buffer to find matching hash.
  for (int i = cur_record->index_start - 1; i >= 0; i--) {
    if (cur_record->hash_vals[i] == hash) return i;
  }
  for (int i = cur_record->num - 1; i >= cur_record->index_start; i--) {
    if (cur_record->hash_vals[i] == hash) return i;
  }
  int index;
  // If not found - add new RD info into the buffer and return its index
  if (cur_record->num < TX_SIZE_RD_RECORD_BUFFER_LEN) {
    index = (cur_record->index_start + cur_record->num) %
            TX_SIZE_RD_RECORD_BUFFER_LEN;
    cur_record->num++;
  } else {
    index = cur_record->index_start;
    cur_record->index_start =
        (cur_record->index_start + 1) % TX_SIZE_RD_RECORD_BUFFER_LEN;
  }

  cur_record->hash_vals[index] = hash;
  av1_zero(cur_record->tx_rd_info[index]);
  return index;
}

// Go through all TX blocks that could be used in TX size search, compute
// residual hash values for them and find matching RD info that stores previous
// RD search results for these TX blocks. The idea is to prevent repeated
// rate/distortion computations that happen because of the combination of
// partition and TX size search. The resulting RD info records are returned in
// the form of a quadtree for easier access in actual TX size search.
static int find_tx_size_rd_records(MACROBLOCK *x, BLOCK_SIZE bsize,
                                   TXB_RD_INFO_NODE *dst_rd_info) {
  TXB_RD_RECORD *rd_records_table[4] = { x->txb_rd_record_8X8,
                                         x->txb_rd_record_16X16,
                                         x->txb_rd_record_32X32,
                                         x->txb_rd_record_64X64 };
  const TX_SIZE max_square_tx_size = max_txsize_lookup[bsize];
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];

  // Hashing is performed only for square TX sizes larger than TX_4X4
  if (max_square_tx_size < TX_8X8) return 0;
  const int diff_stride = bw;
  const struct macroblock_plane *const p = &x->plane[0];
  const int16_t *diff = &p->src_diff[0];
  init_rd_record_tree(dst_rd_info, bsize);
  // Coordinates of the top-left corner of current block within the superblock
  // measured in pixels:
  const int mi_row = x->e_mbd.mi_row;
  const int mi_col = x->e_mbd.mi_col;
  const int mi_row_in_sb = (mi_row % MAX_MIB_SIZE) << MI_SIZE_LOG2;
  const int mi_col_in_sb = (mi_col % MAX_MIB_SIZE) << MI_SIZE_LOG2;
  int cur_rd_info_idx = 0;
  int cur_tx_depth = 0;
  TX_SIZE cur_tx_size = max_txsize_rect_lookup[bsize];
  while (cur_tx_depth <= MAX_VARTX_DEPTH) {
    const int cur_tx_bw = tx_size_wide[cur_tx_size];
    const int cur_tx_bh = tx_size_high[cur_tx_size];
    if (cur_tx_bw < 8 || cur_tx_bh < 8) break;
    const TX_SIZE next_tx_size = sub_tx_size_map[cur_tx_size];
    const int tx_size_idx = cur_tx_size - TX_8X8;
    for (int row = 0; row < bh; row += cur_tx_bh) {
      for (int col = 0; col < bw; col += cur_tx_bw) {
        if (cur_tx_bw != cur_tx_bh) {
          // Use dummy nodes for all rectangular transforms within the
          // TX size search tree.
          dst_rd_info[cur_rd_info_idx].rd_info_array = NULL;
        } else {
          // Get spatial location of this TX block within the superblock
          // (measured in cur_tx_bsize units).
          const int row_in_sb = (mi_row_in_sb + row) / cur_tx_bh;
          const int col_in_sb = (mi_col_in_sb + col) / cur_tx_bw;

          int16_t hash_data[MAX_SB_SQUARE];
          int16_t *cur_hash_row = hash_data;
          const int16_t *cur_diff_row = diff + row * diff_stride + col;
          for (int i = 0; i < cur_tx_bh; i++) {
            memcpy(cur_hash_row, cur_diff_row, sizeof(*hash_data) * cur_tx_bw);
            cur_hash_row += cur_tx_bw;
            cur_diff_row += diff_stride;
          }
          const int hash = av1_get_crc32c_value(&x->mb_rd_record.crc_calculator,
                                                (uint8_t *)hash_data,
                                                2 * cur_tx_bw * cur_tx_bh);
          // Find corresponding RD info based on the hash value.
          const int record_idx =
              row_in_sb * (MAX_MIB_SIZE >> (tx_size_idx + 1)) + col_in_sb;
          TXB_RD_RECORD *records = &rd_records_table[tx_size_idx][record_idx];
          int idx = find_tx_size_rd_info(records, hash);
          dst_rd_info[cur_rd_info_idx].rd_info_array =
              &records->tx_rd_info[idx];
        }
        ++cur_rd_info_idx;
      }
    }
    cur_tx_size = next_tx_size;
    ++cur_tx_depth;
  }
  return 1;
}

int64_t txfm_yrd(const AV1_COMP *const cpi, MACROBLOCK *x, RD_STATS *rd_stats,
                 int64_t ref_best_rd, BLOCK_SIZE bs, TX_SIZE tx_size,
                 FAST_TX_SEARCH_MODE ftxs_mode, int skip_trellis) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  int64_t rd = INT64_MAX;
  const int skip_ctx = av1_get_skip_context(xd);
  int s0, s1;
  const int is_inter = is_inter_block(mbmi);
  const int tx_select = x->tx_mode_search_type == TX_MODE_SELECT &&
                        block_signals_txsize(mbmi->sb_type);
  int ctx = txfm_partition_context(
      xd->above_txfm_context, xd->left_txfm_context, mbmi->sb_type, tx_size);
  const int r_tx_size =
      is_inter ? x->txfm_partition_cost[ctx][0] : tx_size_cost(x, bs, tx_size);

  assert(IMPLIES(is_rect_tx(tx_size), is_rect_tx_allowed_bsize(bs)));

  s0 = x->skip_cost[skip_ctx][0];
  s1 = x->skip_cost[skip_ctx][1];

  int64_t skip_rd = INT64_MAX;
  int64_t this_rd = RDCOST(x->rdmult, s0 + r_tx_size * tx_select, 0);

  if (is_inter) skip_rd = RDCOST(x->rdmult, s1, 0);

  mbmi->tx_size = tx_size;
  txfm_rd_in_plane(
      x, cpi, rd_stats, ref_best_rd, AOMMIN(this_rd, skip_rd), AOM_PLANE_Y, bs,
      tx_size, cpi->sf.rd_sf.use_fast_coef_costing, ftxs_mode, skip_trellis);
  if (rd_stats->rate == INT_MAX) return INT64_MAX;

  // rdstats->rate should include all the rate except skip/non-skip cost as the
  // same is accounted in the caller functions after rd evaluation of all
  // planes. However the decisions should be done after considering the
  // skip/non-skip header cost
  if (rd_stats->skip && is_inter) {
    rd = RDCOST(x->rdmult, s1, rd_stats->sse);
  } else {
    // Intra blocks are always signalled as non-skip
    rd = RDCOST(x->rdmult, rd_stats->rate + s0 + r_tx_size * tx_select,
                rd_stats->dist);
    rd_stats->rate += r_tx_size * tx_select;
  }
  if (is_inter && !xd->lossless[xd->mi[0]->segment_id]) {
    int64_t temp_skip_rd = RDCOST(x->rdmult, s1, rd_stats->sse);
    if (temp_skip_rd <= rd) {
      rd = temp_skip_rd;
      rd_stats->rate = 0;
      rd_stats->dist = rd_stats->sse;
      rd_stats->skip = 1;
    }
  }

  return rd;
}

static AOM_INLINE void choose_tx_size_type_from_rd(const AV1_COMP *const cpi,
                                                   MACROBLOCK *x,
                                                   RD_STATS *rd_stats,
                                                   int64_t ref_best_rd,
                                                   BLOCK_SIZE bs) {
  av1_invalid_rd_stats(rd_stats);

  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const TX_SIZE max_rect_tx_size = max_txsize_rect_lookup[bs];
  const int tx_select = x->tx_mode_search_type == TX_MODE_SELECT;
  int start_tx;
  int depth, init_depth;

  if (tx_select) {
    start_tx = max_rect_tx_size;
    init_depth = get_search_init_depth(mi_size_wide[bs], mi_size_high[bs],
                                       is_inter_block(mbmi), &cpi->sf,
                                       x->tx_size_search_method);
  } else {
    const TX_SIZE chosen_tx_size =
        tx_size_from_tx_mode(bs, x->tx_mode_search_type);
    start_tx = chosen_tx_size;
    init_depth = MAX_TX_DEPTH;
  }

  uint8_t best_txk_type_map[MAX_MIB_SIZE * MAX_MIB_SIZE];
  uint8_t best_blk_skip[MAX_MIB_SIZE * MAX_MIB_SIZE];
  TX_SIZE best_tx_size = max_rect_tx_size;
  int64_t best_rd = INT64_MAX;
  const int n4 = bsize_to_num_blk(bs);
  x->rd_model = FULL_TXFM_RD;
  depth = init_depth;
  int64_t rd[MAX_TX_DEPTH + 1] = { INT64_MAX, INT64_MAX, INT64_MAX };
  for (int n = start_tx; depth <= MAX_TX_DEPTH;
       depth++, n = sub_tx_size_map[n]) {
#if CONFIG_DIST_8X8
    if (x->using_dist_8x8) {
      if (tx_size_wide[n] < 8 || tx_size_high[n] < 8) continue;
    }
#endif
    if (!cpi->oxcf.enable_tx64 && txsize_sqr_up_map[n] == TX_64X64) continue;

    RD_STATS this_rd_stats;
    rd[depth] =
        txfm_yrd(cpi, x, &this_rd_stats, ref_best_rd, bs, n, FTXS_NONE, 0);

    if (rd[depth] < best_rd) {
      av1_copy_array(best_blk_skip, x->blk_skip, n4);
      av1_copy_array(best_txk_type_map, xd->tx_type_map, n4);
      best_tx_size = n;
      best_rd = rd[depth];
      *rd_stats = this_rd_stats;
    }
    if (n == TX_4X4) break;
    // If we are searching three depths, prune the smallest size depending
    // on rd results for the first two depths for low contrast blocks.
    if (depth > init_depth && depth != MAX_TX_DEPTH &&
        x->source_variance < 256) {
      if (rd[depth - 1] != INT64_MAX && rd[depth] > rd[depth - 1]) break;
    }
  }

  if (rd_stats->rate != INT_MAX) {
    mbmi->tx_size = best_tx_size;
    av1_copy_array(xd->tx_type_map, best_txk_type_map, n4);
    av1_copy_array(x->blk_skip, best_blk_skip, n4);
  }
}

// Compute the pixel domain distortion from diff on all visible 4x4s in the
// transform block.
static INLINE int64_t pixel_diff_dist(const MACROBLOCK *x, int plane,
                                      int blk_row, int blk_col,
                                      const BLOCK_SIZE plane_bsize,
                                      const BLOCK_SIZE tx_bsize,
                                      unsigned int *block_mse_q8) {
  int visible_rows, visible_cols;
  const MACROBLOCKD *xd = &x->e_mbd;
  get_txb_dimensions(xd, plane, plane_bsize, blk_row, blk_col, tx_bsize, NULL,
                     NULL, &visible_cols, &visible_rows);
  const int diff_stride = block_size_wide[plane_bsize];
  const int16_t *diff = x->plane[plane].src_diff;
#if CONFIG_DIST_8X8
  int txb_height = block_size_high[tx_bsize];
  int txb_width = block_size_wide[tx_bsize];
  if (x->using_dist_8x8 && plane == 0) {
    const int src_stride = x->plane[plane].src.stride;
    const int src_idx = (blk_row * src_stride + blk_col) << MI_SIZE_LOG2;
    const int diff_idx = (blk_row * diff_stride + blk_col) << MI_SIZE_LOG2;
    const uint8_t *src = &x->plane[plane].src.buf[src_idx];
    return dist_8x8_diff(x, src, src_stride, diff + diff_idx, diff_stride,
                         txb_width, txb_height, visible_cols, visible_rows,
                         x->qindex);
  }
#endif
  diff += ((blk_row * diff_stride + blk_col) << MI_SIZE_LOG2);
  uint64_t sse =
      aom_sum_squares_2d_i16(diff, diff_stride, visible_cols, visible_rows);
  if (block_mse_q8 != NULL) {
    if (visible_cols > 0 && visible_rows > 0)
      *block_mse_q8 =
          (unsigned int)((256 * sse) / (visible_cols * visible_rows));
    else
      *block_mse_q8 = UINT_MAX;
  }
  return sse;
}

// Uses simple features on top of DCT coefficients to quickly predict
// whether optimal RD decision is to skip encoding the residual.
// The sse value is stored in dist.
static int predict_skip_flag(MACROBLOCK *x, BLOCK_SIZE bsize, int64_t *dist,
                             int reduced_tx_set) {
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  const MACROBLOCKD *xd = &x->e_mbd;
  const int16_t dc_q = av1_dc_quant_QTX(x->qindex, 0, xd->bd);

  *dist = pixel_diff_dist(x, 0, 0, 0, bsize, bsize, NULL);

  const int64_t mse = *dist / bw / bh;
  // Normalized quantizer takes the transform upscaling factor (8 for tx size
  // smaller than 32) into account.
  const int16_t normalized_dc_q = dc_q >> 3;
  const int64_t mse_thresh = (int64_t)normalized_dc_q * normalized_dc_q / 8;
  // For faster early skip decision, use dist to compare against threshold so
  // that quality risk is less for the skip=1 decision. Otherwise, use mse
  // since the fwd_txfm coeff checks will take care of quality
  // TODO(any): Use dist to return 0 when predict_skip_level is 1
  int64_t pred_err = (x->predict_skip_level >= 2) ? *dist : mse;
  // Predict not to skip when error is larger than threshold.
  if (pred_err > mse_thresh) return 0;
  // Return as skip otherwise for aggressive early skip
  else if (x->predict_skip_level >= 2)
    return 1;

  const int max_tx_size = max_predict_sf_tx_size[bsize];
  const int tx_h = tx_size_high[max_tx_size];
  const int tx_w = tx_size_wide[max_tx_size];
  DECLARE_ALIGNED(32, tran_low_t, coefs[32 * 32]);
  TxfmParam param;
  param.tx_type = DCT_DCT;
  param.tx_size = max_tx_size;
  param.bd = xd->bd;
  param.is_hbd = is_cur_buf_hbd(xd);
  param.lossless = 0;
  param.tx_set_type = av1_get_ext_tx_set_type(
      param.tx_size, is_inter_block(xd->mi[0]), reduced_tx_set);
  const int bd_idx = (xd->bd == 8) ? 0 : ((xd->bd == 10) ? 1 : 2);
  const uint32_t max_qcoef_thresh = skip_pred_threshold[bd_idx][bsize];
  const int16_t *src_diff = x->plane[0].src_diff;
  const int n_coeff = tx_w * tx_h;
  const int16_t ac_q = av1_ac_quant_QTX(x->qindex, 0, xd->bd);
  const uint32_t dc_thresh = max_qcoef_thresh * dc_q;
  const uint32_t ac_thresh = max_qcoef_thresh * ac_q;
  for (int row = 0; row < bh; row += tx_h) {
    for (int col = 0; col < bw; col += tx_w) {
      av1_fwd_txfm(src_diff + col, coefs, bw, &param);
      // Operating on TX domain, not pixels; we want the QTX quantizers
      const uint32_t dc_coef = (((uint32_t)abs(coefs[0])) << 7);
      if (dc_coef >= dc_thresh) return 0;
      for (int i = 1; i < n_coeff; ++i) {
        const uint32_t ac_coef = (((uint32_t)abs(coefs[i])) << 7);
        if (ac_coef >= ac_thresh) return 0;
      }
    }
    src_diff += tx_h * bw;
  }
  return 1;
}

// Used to set proper context for early termination with skip = 1.
static AOM_INLINE void set_skip_flag(MACROBLOCK *x, RD_STATS *rd_stats,
                                     int bsize, int64_t dist) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int n4 = bsize_to_num_blk(bsize);
  const TX_SIZE tx_size = max_txsize_rect_lookup[bsize];
  memset(xd->tx_type_map, DCT_DCT, sizeof(xd->tx_type_map[0]) * n4);
  memset(mbmi->inter_tx_size, tx_size, sizeof(mbmi->inter_tx_size));
  mbmi->tx_size = tx_size;
  for (int i = 0; i < n4; ++i) set_blk_skip(x, 0, i, 1);
  rd_stats->skip = 1;
  if (is_cur_buf_hbd(xd)) dist = ROUND_POWER_OF_TWO(dist, (xd->bd - 8) * 2);
  rd_stats->dist = rd_stats->sse = (dist << 4);
  // Though decision is to make the block as skip based on luma stats,
  // it is possible that block becomes non skip after chroma rd. In addition
  // intermediate non skip costs calculated by caller function will be
  // incorrect, if rate is set as  zero (i.e., if zero_blk_rate is not
  // accounted). Hence intermediate rate is populated to code the luma tx blks
  // as skip, the caller function based on final rd decision (i.e., skip vs
  // non-skip) sets the final rate accordingly. Here the rate populated
  // corresponds to coding all the tx blocks with zero_blk_rate (based on max tx
  // size possible) in the current block. Eg: For 128*128 block, rate would be
  // 4 * zero_blk_rate where zero_blk_rate corresponds to coding of one 64x64 tx
  // block as 'all zeros'
  ENTROPY_CONTEXT ctxa[MAX_MIB_SIZE];
  ENTROPY_CONTEXT ctxl[MAX_MIB_SIZE];
  av1_get_entropy_contexts(bsize, &xd->plane[0], ctxa, ctxl);
  ENTROPY_CONTEXT *ta = ctxa;
  ENTROPY_CONTEXT *tl = ctxl;
  const TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);
  TXB_CTX txb_ctx;
  get_txb_ctx(bsize, tx_size, 0, ta, tl, &txb_ctx);
  const int zero_blk_rate = x->coeff_costs[txs_ctx][PLANE_TYPE_Y]
                                .txb_skip_cost[txb_ctx.txb_skip_ctx][1];
  rd_stats->rate = zero_blk_rate *
                   (block_size_wide[bsize] >> tx_size_wide_log2[tx_size]) *
                   (block_size_high[bsize] >> tx_size_high_log2[tx_size]);
}

static AOM_INLINE void save_tx_rd_info(int n4, uint32_t hash,
                                       const MACROBLOCK *const x,
                                       const RD_STATS *const rd_stats,
                                       MB_RD_RECORD *tx_rd_record) {
  int index;
  if (tx_rd_record->num < RD_RECORD_BUFFER_LEN) {
    index =
        (tx_rd_record->index_start + tx_rd_record->num) % RD_RECORD_BUFFER_LEN;
    ++tx_rd_record->num;
  } else {
    index = tx_rd_record->index_start;
    tx_rd_record->index_start =
        (tx_rd_record->index_start + 1) % RD_RECORD_BUFFER_LEN;
  }
  MB_RD_INFO *const tx_rd_info = &tx_rd_record->tx_rd_info[index];
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  tx_rd_info->hash_value = hash;
  tx_rd_info->tx_size = mbmi->tx_size;
  memcpy(tx_rd_info->blk_skip, x->blk_skip,
         sizeof(tx_rd_info->blk_skip[0]) * n4);
  av1_copy(tx_rd_info->inter_tx_size, mbmi->inter_tx_size);
  av1_copy_array(tx_rd_info->tx_type_map, xd->tx_type_map, n4);
  tx_rd_info->rd_stats = *rd_stats;
}

static INLINE uint32_t get_block_residue_hash(MACROBLOCK *x, BLOCK_SIZE bsize) {
  const int rows = block_size_high[bsize];
  const int cols = block_size_wide[bsize];
  const int16_t *diff = x->plane[0].src_diff;
  const uint32_t hash = av1_get_crc32c_value(&x->mb_rd_record.crc_calculator,
                                             (uint8_t *)diff, 2 * rows * cols);
  return (hash << 5) + bsize;
}

static AOM_INLINE void inverse_transform_block_facade(MACROBLOCKD *xd,
                                                      int plane, int block,
                                                      int blk_row, int blk_col,
                                                      int eob,
                                                      int reduced_tx_set) {
  if (!eob) return;

  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *dqcoeff = pd->dqcoeff + BLOCK_OFFSET(block);
  const PLANE_TYPE plane_type = get_plane_type(plane);
  const TX_SIZE tx_size = av1_get_tx_size(plane, xd);
  const TX_TYPE tx_type = av1_get_tx_type(xd, plane_type, blk_row, blk_col,
                                          tx_size, reduced_tx_set);
  const int dst_stride = pd->dst.stride;
  uint8_t *dst = &pd->dst.buf[(blk_row * dst_stride + blk_col) << MI_SIZE_LOG2];
  av1_inverse_transform_block(xd, dqcoeff, plane, tx_type, tx_size, dst,
                              dst_stride, eob, reduced_tx_set);
}

static int find_tx_size_rd_info(TXB_RD_RECORD *cur_record, const uint32_t hash);

static uint32_t get_intra_txb_hash(MACROBLOCK *x, int plane, int blk_row,
                                   int blk_col, BLOCK_SIZE plane_bsize,
                                   TX_SIZE tx_size) {
  int16_t tmp_data[64 * 64];
  const int diff_stride = block_size_wide[plane_bsize];
  const int16_t *diff = x->plane[plane].src_diff;
  const int16_t *cur_diff_row = diff + 4 * blk_row * diff_stride + 4 * blk_col;
  const int txb_w = tx_size_wide[tx_size];
  const int txb_h = tx_size_high[tx_size];
  uint8_t *hash_data = (uint8_t *)cur_diff_row;
  if (txb_w != diff_stride) {
    int16_t *cur_hash_row = tmp_data;
    for (int i = 0; i < txb_h; i++) {
      memcpy(cur_hash_row, cur_diff_row, sizeof(*diff) * txb_w);
      cur_hash_row += txb_w;
      cur_diff_row += diff_stride;
    }
    hash_data = (uint8_t *)tmp_data;
  }
  CRC32C *crc = &x->mb_rd_record.crc_calculator;
  const uint32_t hash = av1_get_crc32c_value(crc, hash_data, 2 * txb_w * txb_h);
  return (hash << 5) + tx_size;
}

static INLINE void dist_block_tx_domain(MACROBLOCK *x, int plane, int block,
                                        TX_SIZE tx_size, int64_t *out_dist,
                                        int64_t *out_sse) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  // Transform domain distortion computation is more efficient as it does
  // not involve an inverse transform, but it is less accurate.
  const int buffer_length = av1_get_max_eob(tx_size);
  int64_t this_sse;
  // TX-domain results need to shift down to Q2/D10 to match pixel
  // domain distortion values which are in Q2^2
  int shift = (MAX_TX_SCALE - av1_get_tx_scale(tx_size)) * 2;
  const int block_offset = BLOCK_OFFSET(block);
  tran_low_t *const coeff = p->coeff + block_offset;
  tran_low_t *const dqcoeff = pd->dqcoeff + block_offset;
#if CONFIG_AV1_HIGHBITDEPTH
  if (is_cur_buf_hbd(xd))
    *out_dist = av1_highbd_block_error(coeff, dqcoeff, buffer_length, &this_sse,
                                       xd->bd);
  else
    *out_dist = av1_block_error(coeff, dqcoeff, buffer_length, &this_sse);
#else
  *out_dist = av1_block_error(coeff, dqcoeff, buffer_length, &this_sse);
#endif
  *out_dist = RIGHT_SIGNED_SHIFT(*out_dist, shift);
  *out_sse = RIGHT_SIGNED_SHIFT(this_sse, shift);
}

static unsigned pixel_dist_visible_only(
    const AV1_COMP *const cpi, const MACROBLOCK *x, const uint8_t *src,
    const int src_stride, const uint8_t *dst, const int dst_stride,
    const BLOCK_SIZE tx_bsize, int txb_rows, int txb_cols, int visible_rows,
    int visible_cols) {
  unsigned sse;

  if (txb_rows == visible_rows && txb_cols == visible_cols) {
    cpi->fn_ptr[tx_bsize].vf(src, src_stride, dst, dst_stride, &sse);
    return sse;
  }

#if CONFIG_AV1_HIGHBITDEPTH
  const MACROBLOCKD *xd = &x->e_mbd;
  if (is_cur_buf_hbd(xd)) {
    uint64_t sse64 = aom_highbd_sse_odd_size(src, src_stride, dst, dst_stride,
                                             visible_cols, visible_rows);
    return (unsigned int)ROUND_POWER_OF_TWO(sse64, (xd->bd - 8) * 2);
  }
#else
  (void)x;
#endif
  sse = aom_sse_odd_size(src, src_stride, dst, dst_stride, visible_cols,
                         visible_rows);
  return sse;
}

// Compute the pixel domain distortion from src and dst on all visible 4x4s in
// the
// transform block.
static unsigned pixel_dist(const AV1_COMP *const cpi, const MACROBLOCK *x,
                           int plane, const uint8_t *src, const int src_stride,
                           const uint8_t *dst, const int dst_stride,
                           int blk_row, int blk_col,
                           const BLOCK_SIZE plane_bsize,
                           const BLOCK_SIZE tx_bsize) {
  int txb_rows, txb_cols, visible_rows, visible_cols;
  const MACROBLOCKD *xd = &x->e_mbd;

  get_txb_dimensions(xd, plane, plane_bsize, blk_row, blk_col, tx_bsize,
                     &txb_cols, &txb_rows, &visible_cols, &visible_rows);
  assert(visible_rows > 0);
  assert(visible_cols > 0);

#if CONFIG_DIST_8X8
  if (x->using_dist_8x8 && plane == 0)
    return (unsigned)av1_dist_8x8(cpi, x, src, src_stride, dst, dst_stride,
                                  tx_bsize, txb_cols, txb_rows, visible_cols,
                                  visible_rows, x->qindex);
#endif  // CONFIG_DIST_8X8

  unsigned sse = pixel_dist_visible_only(cpi, x, src, src_stride, dst,
                                         dst_stride, tx_bsize, txb_rows,
                                         txb_cols, visible_rows, visible_cols);

  return sse;
}

static INLINE int64_t dist_block_px_domain(const AV1_COMP *cpi, MACROBLOCK *x,
                                           int plane, BLOCK_SIZE plane_bsize,
                                           int block, int blk_row, int blk_col,
                                           TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const uint16_t eob = p->eobs[block];
  const BLOCK_SIZE tx_bsize = txsize_to_bsize[tx_size];
  const int bsw = block_size_wide[tx_bsize];
  const int bsh = block_size_high[tx_bsize];
  const int src_stride = x->plane[plane].src.stride;
  const int dst_stride = xd->plane[plane].dst.stride;
  // Scale the transform block index to pixel unit.
  const int src_idx = (blk_row * src_stride + blk_col) << MI_SIZE_LOG2;
  const int dst_idx = (blk_row * dst_stride + blk_col) << MI_SIZE_LOG2;
  const uint8_t *src = &x->plane[plane].src.buf[src_idx];
  const uint8_t *dst = &xd->plane[plane].dst.buf[dst_idx];
  const tran_low_t *dqcoeff = pd->dqcoeff + BLOCK_OFFSET(block);

  assert(cpi != NULL);
  assert(tx_size_wide_log2[0] == tx_size_high_log2[0]);

  uint8_t *recon;
  DECLARE_ALIGNED(16, uint16_t, recon16[MAX_TX_SQUARE]);

#if CONFIG_AV1_HIGHBITDEPTH
  if (is_cur_buf_hbd(xd)) {
    recon = CONVERT_TO_BYTEPTR(recon16);
    av1_highbd_convolve_2d_copy_sr(CONVERT_TO_SHORTPTR(dst), dst_stride,
                                   CONVERT_TO_SHORTPTR(recon), MAX_TX_SIZE, bsw,
                                   bsh, NULL, NULL, 0, 0, NULL, xd->bd);
  } else {
    recon = (uint8_t *)recon16;
    av1_convolve_2d_copy_sr(dst, dst_stride, recon, MAX_TX_SIZE, bsw, bsh, NULL,
                            NULL, 0, 0, NULL);
  }
#else
  recon = (uint8_t *)recon16;
  av1_convolve_2d_copy_sr(dst, dst_stride, recon, MAX_TX_SIZE, bsw, bsh, NULL,
                          NULL, 0, 0, NULL);
#endif

  const PLANE_TYPE plane_type = get_plane_type(plane);
  TX_TYPE tx_type = av1_get_tx_type(xd, plane_type, blk_row, blk_col, tx_size,
                                    cpi->common.reduced_tx_set_used);
  av1_inverse_transform_block(xd, dqcoeff, plane, tx_type, tx_size, recon,
                              MAX_TX_SIZE, eob,
                              cpi->common.reduced_tx_set_used);

  return 16 * pixel_dist(cpi, x, plane, src, src_stride, recon, MAX_TX_SIZE,
                         blk_row, blk_col, plane_bsize, tx_bsize);
}

static INLINE int32_t find_mb_rd_info(const MB_RD_RECORD *const mb_rd_record,
                                      const int64_t ref_best_rd,
                                      const uint32_t hash) {
  int32_t match_index = -1;
  if (ref_best_rd != INT64_MAX) {
    for (int i = 0; i < mb_rd_record->num; ++i) {
      const int index = (mb_rd_record->index_start + i) % RD_RECORD_BUFFER_LEN;
      // If there is a match in the tx_rd_record, fetch the RD decision and
      // terminate early.
      if (mb_rd_record->tx_rd_info[index].hash_value == hash) {
        match_index = index;
        break;
      }
    }
  }
  return match_index;
}

static AOM_INLINE void fetch_tx_rd_info(int n4,
                                        const MB_RD_INFO *const tx_rd_info,
                                        RD_STATS *const rd_stats,
                                        MACROBLOCK *const x) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  mbmi->tx_size = tx_rd_info->tx_size;
  memcpy(x->blk_skip, tx_rd_info->blk_skip,
         sizeof(tx_rd_info->blk_skip[0]) * n4);
  av1_copy(mbmi->inter_tx_size, tx_rd_info->inter_tx_size);
  av1_copy_array(xd->tx_type_map, tx_rd_info->tx_type_map, n4);
  *rd_stats = tx_rd_info->rd_stats;
}

static AOM_INLINE void get_energy_distribution_finer(const int16_t *diff,
                                                     int stride, int bw, int bh,
                                                     float *hordist,
                                                     float *verdist) {
  // First compute downscaled block energy values (esq); downscale factors
  // are defined by w_shift and h_shift.
  unsigned int esq[256];
  const int w_shift = bw <= 8 ? 0 : 1;
  const int h_shift = bh <= 8 ? 0 : 1;
  const int esq_w = bw >> w_shift;
  const int esq_h = bh >> h_shift;
  const int esq_sz = esq_w * esq_h;
  int i, j;
  memset(esq, 0, esq_sz * sizeof(esq[0]));
  if (w_shift) {
    for (i = 0; i < bh; i++) {
      unsigned int *cur_esq_row = esq + (i >> h_shift) * esq_w;
      const int16_t *cur_diff_row = diff + i * stride;
      for (j = 0; j < bw; j += 2) {
        cur_esq_row[j >> 1] += (cur_diff_row[j] * cur_diff_row[j] +
                                cur_diff_row[j + 1] * cur_diff_row[j + 1]);
      }
    }
  } else {
    for (i = 0; i < bh; i++) {
      unsigned int *cur_esq_row = esq + (i >> h_shift) * esq_w;
      const int16_t *cur_diff_row = diff + i * stride;
      for (j = 0; j < bw; j++) {
        cur_esq_row[j] += cur_diff_row[j] * cur_diff_row[j];
      }
    }
  }

  uint64_t total = 0;
  for (i = 0; i < esq_sz; i++) total += esq[i];

  // Output hordist and verdist arrays are normalized 1D projections of esq
  if (total == 0) {
    float hor_val = 1.0f / esq_w;
    for (j = 0; j < esq_w - 1; j++) hordist[j] = hor_val;
    float ver_val = 1.0f / esq_h;
    for (i = 0; i < esq_h - 1; i++) verdist[i] = ver_val;
    return;
  }

  const float e_recip = 1.0f / (float)total;
  memset(hordist, 0, (esq_w - 1) * sizeof(hordist[0]));
  memset(verdist, 0, (esq_h - 1) * sizeof(verdist[0]));
  const unsigned int *cur_esq_row;
  for (i = 0; i < esq_h - 1; i++) {
    cur_esq_row = esq + i * esq_w;
    for (j = 0; j < esq_w - 1; j++) {
      hordist[j] += (float)cur_esq_row[j];
      verdist[i] += (float)cur_esq_row[j];
    }
    verdist[i] += (float)cur_esq_row[j];
  }
  cur_esq_row = esq + i * esq_w;
  for (j = 0; j < esq_w - 1; j++) hordist[j] += (float)cur_esq_row[j];

  for (j = 0; j < esq_w - 1; j++) hordist[j] *= e_recip;
  for (i = 0; i < esq_h - 1; i++) verdist[i] *= e_recip;
}

// These thresholds were calibrated to provide a certain number of TX types
// pruned by the model on average, i.e. selecting a threshold with index i
// will lead to pruning i+1 TX types on average
static const float *prune_2D_adaptive_thresholds[] = {
  // TX_4X4
  (float[]){ 0.00549f, 0.01306f, 0.02039f, 0.02747f, 0.03406f, 0.04065f,
             0.04724f, 0.05383f, 0.06067f, 0.06799f, 0.07605f, 0.08533f,
             0.09778f, 0.11780f },
  // TX_8X8
  (float[]){ 0.00037f, 0.00183f, 0.00525f, 0.01038f, 0.01697f, 0.02502f,
             0.03381f, 0.04333f, 0.05286f, 0.06287f, 0.07434f, 0.08850f,
             0.10803f, 0.14124f },
  // TX_16X16
  (float[]){ 0.01404f, 0.02000f, 0.04211f, 0.05164f, 0.05798f, 0.06335f,
             0.06897f, 0.07629f, 0.08875f, 0.11169f },
  // TX_32X32
  NULL,
  // TX_64X64
  NULL,
  // TX_4X8
  (float[]){ 0.00183f, 0.00745f, 0.01428f, 0.02185f, 0.02966f, 0.03723f,
             0.04456f, 0.05188f, 0.05920f, 0.06702f, 0.07605f, 0.08704f,
             0.10168f, 0.12585f },
  // TX_8X4
  (float[]){ 0.00085f, 0.00476f, 0.01135f, 0.01892f, 0.02698f, 0.03528f,
             0.04358f, 0.05164f, 0.05994f, 0.06848f, 0.07849f, 0.09021f,
             0.10583f, 0.13123f },
  // TX_8X16
  (float[]){ 0.00037f, 0.00232f, 0.00671f, 0.01257f, 0.01965f, 0.02722f,
             0.03552f, 0.04382f, 0.05237f, 0.06189f, 0.07336f, 0.08728f,
             0.10730f, 0.14221f },
  // TX_16X8
  (float[]){ 0.00061f, 0.00330f, 0.00818f, 0.01453f, 0.02185f, 0.02966f,
             0.03772f, 0.04578f, 0.05383f, 0.06262f, 0.07288f, 0.08582f,
             0.10339f, 0.13464f },
  // TX_16X32
  NULL,
  // TX_32X16
  NULL,
  // TX_32X64
  NULL,
  // TX_64X32
  NULL,
  // TX_4X16
  (float[]){ 0.00232f, 0.00671f, 0.01257f, 0.01941f, 0.02673f, 0.03430f,
             0.04211f, 0.04968f, 0.05750f, 0.06580f, 0.07507f, 0.08655f,
             0.10242f, 0.12878f },
  // TX_16X4
  (float[]){ 0.00110f, 0.00525f, 0.01208f, 0.01990f, 0.02795f, 0.03601f,
             0.04358f, 0.05115f, 0.05896f, 0.06702f, 0.07629f, 0.08752f,
             0.10217f, 0.12610f },
  // TX_8X32
  NULL,
  // TX_32X8
  NULL,
  // TX_16X64
  NULL,
  // TX_64X16
  NULL,
};

static float get_dev(float mean, double x2_sum, int num) {
  const float e_x2 = (float)(x2_sum / num);
  const float diff = e_x2 - mean * mean;
  const float dev = (diff > 0) ? sqrtf(diff) : 0;
  return dev;
}

// Feature used by the model to predict tx split: the mean and standard
// deviation values of the block and sub-blocks.
static AOM_INLINE void get_mean_dev_features(const int16_t *data, int stride,
                                             int bw, int bh, float *feature) {
  const int16_t *const data_ptr = &data[0];
  const int subh = (bh >= bw) ? (bh >> 1) : bh;
  const int subw = (bw >= bh) ? (bw >> 1) : bw;
  const int num = bw * bh;
  const int sub_num = subw * subh;
  int feature_idx = 2;
  int total_x_sum = 0;
  int64_t total_x2_sum = 0;
  int blk_idx = 0;
  double mean2_sum = 0.0f;
  float dev_sum = 0.0f;

  for (int row = 0; row < bh; row += subh) {
    for (int col = 0; col < bw; col += subw) {
      int x_sum;
      int64_t x2_sum;
      // TODO(any): Write a SIMD version. Clear registers.
      aom_get_blk_sse_sum(data_ptr + row * stride + col, stride, subw, subh,
                          &x_sum, &x2_sum);
      total_x_sum += x_sum;
      total_x2_sum += x2_sum;

      aom_clear_system_state();
      const float mean = (float)x_sum / sub_num;
      const float dev = get_dev(mean, (double)x2_sum, sub_num);
      feature[feature_idx++] = mean;
      feature[feature_idx++] = dev;
      mean2_sum += (double)(mean * mean);
      dev_sum += dev;
      blk_idx++;
    }
  }

  const float lvl0_mean = (float)total_x_sum / num;
  feature[0] = lvl0_mean;
  feature[1] = get_dev(lvl0_mean, (double)total_x2_sum, num);

  if (blk_idx > 1) {
    // Deviation of means.
    feature[feature_idx++] = get_dev(lvl0_mean, mean2_sum, blk_idx);
    // Mean of deviations.
    feature[feature_idx++] = dev_sum / blk_idx;
  }
}

static int ml_predict_tx_split(MACROBLOCK *x, BLOCK_SIZE bsize, int blk_row,
                               int blk_col, TX_SIZE tx_size) {
  const NN_CONFIG *nn_config = av1_tx_split_nnconfig_map[tx_size];
  if (!nn_config) return -1;

  const int diff_stride = block_size_wide[bsize];
  const int16_t *diff =
      x->plane[0].src_diff + 4 * blk_row * diff_stride + 4 * blk_col;
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];
  aom_clear_system_state();

  float features[64] = { 0.0f };
  get_mean_dev_features(diff, diff_stride, bw, bh, features);

  float score = 0.0f;
  av1_nn_predict(features, nn_config, 1, &score);
  aom_clear_system_state();

  int int_score = (int)(score * 10000);
  return clamp(int_score, -80000, 80000);
}

// Probablities are sorted in descending order.
static INLINE void sort_probability(float prob[], int txk[], int len) {
  int i, j, k;

  for (i = 1; i <= len - 1; ++i) {
    for (j = 0; j < i; ++j) {
      if (prob[j] < prob[i]) {
        float temp;
        int tempi;

        temp = prob[i];
        tempi = txk[i];

        for (k = i; k > j; k--) {
          prob[k] = prob[k - 1];
          txk[k] = txk[k - 1];
        }

        prob[j] = temp;
        txk[j] = tempi;
        break;
      }
    }
  }
}

static void prune_tx_2D(MACROBLOCK *x, BLOCK_SIZE bsize, TX_SIZE tx_size,
                        int blk_row, int blk_col, TxSetType tx_set_type,
                        TX_TYPE_PRUNE_MODE prune_mode, int *txk_map,
                        uint16_t *allowed_tx_mask) {
  int tx_type_table_2D[16] = {
    DCT_DCT,      DCT_ADST,      DCT_FLIPADST,      V_DCT,
    ADST_DCT,     ADST_ADST,     ADST_FLIPADST,     V_ADST,
    FLIPADST_DCT, FLIPADST_ADST, FLIPADST_FLIPADST, V_FLIPADST,
    H_DCT,        H_ADST,        H_FLIPADST,        IDTX
  };
  if (tx_set_type != EXT_TX_SET_ALL16 &&
      tx_set_type != EXT_TX_SET_DTT9_IDTX_1DDCT)
    return;
#if CONFIG_NN_V2
  NN_CONFIG_V2 *nn_config_hor = av1_tx_type_nnconfig_map_hor[tx_size];
  NN_CONFIG_V2 *nn_config_ver = av1_tx_type_nnconfig_map_ver[tx_size];
#else
  const NN_CONFIG *nn_config_hor = av1_tx_type_nnconfig_map_hor[tx_size];
  const NN_CONFIG *nn_config_ver = av1_tx_type_nnconfig_map_ver[tx_size];
#endif
  if (!nn_config_hor || !nn_config_ver) return;  // Model not established yet.

  aom_clear_system_state();
  float hfeatures[16], vfeatures[16];
  float hscores[4], vscores[4];
  float scores_2D_raw[16];
  float scores_2D[16];
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];
  const int hfeatures_num = bw <= 8 ? bw : bw / 2;
  const int vfeatures_num = bh <= 8 ? bh : bh / 2;
  assert(hfeatures_num <= 16);
  assert(vfeatures_num <= 16);

  const struct macroblock_plane *const p = &x->plane[0];
  const int diff_stride = block_size_wide[bsize];
  const int16_t *diff = p->src_diff + 4 * blk_row * diff_stride + 4 * blk_col;
  get_energy_distribution_finer(diff, diff_stride, bw, bh, hfeatures,
                                vfeatures);
  av1_get_horver_correlation_full(diff, diff_stride, bw, bh,
                                  &hfeatures[hfeatures_num - 1],
                                  &vfeatures[vfeatures_num - 1]);
  aom_clear_system_state();
#if CONFIG_NN_V2
  av1_nn_predict_v2(hfeatures, nn_config_hor, 0, hscores);
  av1_nn_predict_v2(vfeatures, nn_config_ver, 0, vscores);
#else
  av1_nn_predict(hfeatures, nn_config_hor, 1, hscores);
  av1_nn_predict(vfeatures, nn_config_ver, 1, vscores);
#endif
  aom_clear_system_state();

  for (int i = 0; i < 4; i++) {
    float *cur_scores_2D = scores_2D_raw + i * 4;
    cur_scores_2D[0] = vscores[i] * hscores[0];
    cur_scores_2D[1] = vscores[i] * hscores[1];
    cur_scores_2D[2] = vscores[i] * hscores[2];
    cur_scores_2D[3] = vscores[i] * hscores[3];
  }

  av1_nn_softmax(scores_2D_raw, scores_2D, 16);

  const int prune_aggr_table[4][2] = { { 4, 1 }, { 6, 3 }, { 9, 6 }, { 9, 6 } };
  int pruning_aggressiveness = 0;
  if (tx_set_type == EXT_TX_SET_ALL16) {
    pruning_aggressiveness =
        prune_aggr_table[prune_mode - PRUNE_2D_ACCURATE][0];
  } else if (tx_set_type == EXT_TX_SET_DTT9_IDTX_1DDCT) {
    pruning_aggressiveness =
        prune_aggr_table[prune_mode - PRUNE_2D_ACCURATE][1];
  }

  // Always keep the TX type with the highest score, prune all others with
  // score below score_thresh.
  int max_score_i = 0;
  float max_score = 0.0f;
  for (int i = 0; i < 16; i++) {
    if (scores_2D[i] > max_score &&
        (*allowed_tx_mask & (1 << tx_type_table_2D[i]))) {
      max_score = scores_2D[i];
      max_score_i = i;
    }
  }

  const float score_thresh =
      prune_2D_adaptive_thresholds[tx_size][pruning_aggressiveness];

  uint16_t allow_bitmask = 0;
  float sum_score = 0.0;
  // Calculate sum of allowed tx type score and Populate allow bit mask based
  // on score_thresh and allowed_tx_mask
  for (int tx_idx = 0; tx_idx < TX_TYPES; tx_idx++) {
    int allow_tx_type = *allowed_tx_mask & (1 << tx_type_table_2D[tx_idx]);
    if ((scores_2D[tx_idx] >= score_thresh && allow_tx_type) ||
        tx_idx == max_score_i) {
      // Set allow mask based on score_thresh and tx type with max score
      allow_bitmask |= (1 << tx_type_table_2D[tx_idx]);

      // Accumulate score of allowed tx type
      sum_score += scores_2D[tx_idx];
    }
  }
  // Sort tx type probability of all types
  sort_probability(scores_2D, tx_type_table_2D, TX_TYPES);

  // Enable more pruning based on tx type probability and number of allowed tx
  // types
  if (prune_mode == PRUNE_2D_AGGRESSIVE) {
    float temp_score = 0.0;
    float score_ratio = 0.0;
    int tx_idx, tx_count = 0;
    const float inv_sum_score = 100 / sum_score;
    // Get allowed tx types based on sorted probability score and tx count
    for (tx_idx = 0; tx_idx < TX_TYPES; tx_idx++) {
      // Skip the tx type which has more than 30% of cumulative
      // probability and allowed tx type count is more than 2
      if (score_ratio > 30.0 && tx_count >= 2) break;

      // Calculate cumulative probability of allowed tx types
      if (allow_bitmask & (1 << tx_type_table_2D[tx_idx])) {
        // Calculate cumulative probability
        temp_score += scores_2D[tx_idx];

        // Calculate percentage of cumulative probability of allowed tx type
        score_ratio = temp_score * inv_sum_score;
        tx_count++;
      }
    }
    // Set remaining tx types as pruned
    for (; tx_idx < TX_TYPES; tx_idx++)
      allow_bitmask &= ~(1 << tx_type_table_2D[tx_idx]);
  }
  memcpy(txk_map, tx_type_table_2D, sizeof(tx_type_table_2D));
  *allowed_tx_mask = allow_bitmask;
}

static int64_t search_txk_type(const AV1_COMP *cpi, MACROBLOCK *x, int plane,
                               int block, int blk_row, int blk_col,
                               BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                               const TXB_CTX *const txb_ctx,
                               FAST_TX_SEARCH_MODE ftxs_mode,
                               int use_fast_coef_costing, int skip_trellis,
                               int64_t ref_best_rd, RD_STATS *best_rd_stats) {
  const AV1_COMMON *cm = &cpi->common;
  MACROBLOCKD *xd = &x->e_mbd;
  struct macroblockd_plane *const pd = &xd->plane[plane];
  MB_MODE_INFO *mbmi = xd->mi[0];
  const int is_inter = is_inter_block(mbmi);
  int64_t best_rd = INT64_MAX;
  uint16_t best_eob = 0;
  TX_TYPE best_tx_type = DCT_DCT;
  TX_TYPE last_tx_type = TX_TYPES;
  const int fast_tx_search = ftxs_mode & FTXS_DCT_AND_1D_DCT_ONLY;
  // The buffer used to swap dqcoeff in macroblockd_plane so we can keep dqcoeff
  // of the best tx_type
  DECLARE_ALIGNED(32, tran_low_t, this_dqcoeff[MAX_SB_SQUARE]);
  tran_low_t *orig_dqcoeff = pd->dqcoeff;
  tran_low_t *best_dqcoeff = this_dqcoeff;
  const int tx_type_map_idx =
      plane ? 0 : blk_row * xd->tx_type_map_stride + blk_col;
  int perform_block_coeff_opt = 0;
  av1_invalid_rd_stats(best_rd_stats);

  TXB_RD_INFO *intra_txb_rd_info = NULL;
  uint16_t cur_joint_ctx = 0;
  const int mi_row = xd->mi_row;
  const int mi_col = xd->mi_col;
  const int within_border =
      mi_row >= xd->tile.mi_row_start &&
      (mi_row + mi_size_high[plane_bsize] < xd->tile.mi_row_end) &&
      mi_col >= xd->tile.mi_col_start &&
      (mi_col + mi_size_wide[plane_bsize] < xd->tile.mi_col_end);
  skip_trellis |=
      cpi->optimize_seg_arr[mbmi->segment_id] == NO_TRELLIS_OPT ||
      cpi->optimize_seg_arr[mbmi->segment_id] == FINAL_PASS_TRELLIS_OPT;
  if (within_border && cpi->sf.tx_sf.use_intra_txb_hash &&
      frame_is_intra_only(cm) && !is_inter && plane == 0 &&
      tx_size_wide[tx_size] == tx_size_high[tx_size]) {
    const uint32_t intra_hash =
        get_intra_txb_hash(x, plane, blk_row, blk_col, plane_bsize, tx_size);
    const int intra_hash_idx =
        find_tx_size_rd_info(&x->txb_rd_record_intra, intra_hash);
    intra_txb_rd_info = &x->txb_rd_record_intra.tx_rd_info[intra_hash_idx];

    cur_joint_ctx = (txb_ctx->dc_sign_ctx << 8) + txb_ctx->txb_skip_ctx;
    if (intra_txb_rd_info->entropy_context == cur_joint_ctx &&
        x->txb_rd_record_intra.tx_rd_info[intra_hash_idx].valid) {
      xd->tx_type_map[tx_type_map_idx] = intra_txb_rd_info->tx_type;
      const TX_TYPE ref_tx_type =
          av1_get_tx_type(xd, get_plane_type(plane), blk_row, blk_col, tx_size,
                          cpi->common.reduced_tx_set_used);
      if (ref_tx_type == intra_txb_rd_info->tx_type) {
        best_rd_stats->rate = intra_txb_rd_info->rate;
        best_rd_stats->dist = intra_txb_rd_info->dist;
        best_rd_stats->sse = intra_txb_rd_info->sse;
        best_rd_stats->skip = intra_txb_rd_info->eob == 0;
        x->plane[plane].eobs[block] = intra_txb_rd_info->eob;
        x->plane[plane].txb_entropy_ctx[block] =
            intra_txb_rd_info->txb_entropy_ctx;
        best_rd = RDCOST(x->rdmult, best_rd_stats->rate, best_rd_stats->dist);
        best_eob = intra_txb_rd_info->eob;
        best_tx_type = intra_txb_rd_info->tx_type;
        perform_block_coeff_opt = intra_txb_rd_info->perform_block_coeff_opt;
        skip_trellis |= !perform_block_coeff_opt;
        update_txk_array(xd, blk_row, blk_col, tx_size, best_tx_type);
        goto RECON_INTRA;
      }
    }
  }

  int rate_cost = 0;
  // if txk_allowed = TX_TYPES, >1 tx types are allowed, else, if txk_allowed <
  // TX_TYPES, only that specific tx type is allowed.
  TX_TYPE txk_allowed = TX_TYPES;
  int txk_map[TX_TYPES] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
  };

  if ((!is_inter && x->use_default_intra_tx_type) ||
      (is_inter && x->use_default_inter_tx_type)) {
    txk_allowed =
        get_default_tx_type(0, xd, tx_size, cpi->is_screen_content_type);
  } else if (x->rd_model == LOW_TXFM_RD) {
    if (plane == 0) txk_allowed = DCT_DCT;
  }

  uint8_t best_txb_ctx = 0;
  const TxSetType tx_set_type =
      av1_get_ext_tx_set_type(tx_size, is_inter, cm->reduced_tx_set_used);

  TX_TYPE uv_tx_type = DCT_DCT;
  if (plane) {
    // tx_type of PLANE_TYPE_UV should be the same as PLANE_TYPE_Y
    uv_tx_type = txk_allowed =
        av1_get_tx_type(xd, get_plane_type(plane), blk_row, blk_col, tx_size,
                        cm->reduced_tx_set_used);
  }
  PREDICTION_MODE intra_dir =
      mbmi->filter_intra_mode_info.use_filter_intra
          ? fimode_to_intradir[mbmi->filter_intra_mode_info.filter_intra_mode]
          : mbmi->mode;
  const uint16_t ext_tx_used_flag =
      cpi->sf.tx_sf.tx_type_search.use_reduced_intra_txset &&
              tx_set_type == EXT_TX_SET_DTT4_IDTX_1DDCT
          ? av1_reduced_intra_tx_used_flag[intra_dir]
          : av1_ext_tx_used_flag[tx_set_type];
  if (xd->lossless[mbmi->segment_id] || txsize_sqr_up_map[tx_size] > TX_32X32 ||
      ext_tx_used_flag == 0x0001 ||
      (is_inter && cpi->oxcf.use_inter_dct_only) ||
      (!is_inter && cpi->oxcf.use_intra_dct_only)) {
    txk_allowed = DCT_DCT;
  }
  uint16_t allowed_tx_mask = 0;  // 1: allow; 0: skip.
  if (txk_allowed < TX_TYPES) {
    allowed_tx_mask = 1 << txk_allowed;
    allowed_tx_mask &= ext_tx_used_flag;
  } else if (fast_tx_search) {
    allowed_tx_mask = 0x0c01;  // V_DCT, H_DCT, DCT_DCT
    allowed_tx_mask &= ext_tx_used_flag;
  } else {
    assert(plane == 0);
    allowed_tx_mask = ext_tx_used_flag;
    int num_allowed = 0;
    const FRAME_UPDATE_TYPE update_type = get_frame_update_type(&cpi->gf_group);
    const int *tx_type_probs = cpi->tx_type_probs[update_type][tx_size];
    int i;

    if (cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats) {
      const int thresh = cpi->tx_type_probs_thresh[update_type];
      uint16_t prune = 0;
      int max_prob = -1;
      int max_idx = 0;
      for (i = 0; i < TX_TYPES; i++) {
        if (tx_type_probs[i] > max_prob && (allowed_tx_mask & (1 << i))) {
          max_prob = tx_type_probs[i];
          max_idx = i;
        }
      }

      for (i = 0; i < TX_TYPES; i++) {
        if (tx_type_probs[i] < thresh && i != max_idx) prune |= (1 << i);
      }
      allowed_tx_mask &= (~prune);
    }

    for (i = 0; i < TX_TYPES; i++) {
      if (allowed_tx_mask & (1 << i)) num_allowed++;
    }
    assert(num_allowed > 0);

    int allowed_tx_count = (x->prune_mode == PRUNE_2D_AGGRESSIVE) ? 1 : 5;
    // !fast_tx_search && txk_end != txk_start && plane == 0
    if (x->prune_mode >= PRUNE_2D_ACCURATE && is_inter &&
        num_allowed > allowed_tx_count) {
      prune_tx_2D(x, plane_bsize, tx_size, blk_row, blk_col, tx_set_type,
                  x->prune_mode, txk_map, &allowed_tx_mask);
    }
  }

  if (cpi->oxcf.enable_flip_idtx == 0) {
    for (TX_TYPE tx_type = FLIPADST_DCT; tx_type <= H_FLIPADST; ++tx_type) {
      allowed_tx_mask &= ~(1 << tx_type);
    }
  }

  // Need to have at least one transform type allowed.
  if (allowed_tx_mask == 0) {
    txk_allowed = (plane ? uv_tx_type : DCT_DCT);
    allowed_tx_mask = (1 << txk_allowed);
  }

  const BLOCK_SIZE tx_bsize = txsize_to_bsize[tx_size];
  int64_t block_sse = 0;
  unsigned int block_mse_q8 = UINT_MAX;
  block_sse = pixel_diff_dist(x, plane, blk_row, blk_col, plane_bsize, tx_bsize,
                              &block_mse_q8);
  assert(block_mse_q8 != UINT_MAX);
  if (is_cur_buf_hbd(xd)) {
    block_sse = ROUND_POWER_OF_TWO(block_sse, (xd->bd - 8) * 2);
    block_mse_q8 = ROUND_POWER_OF_TWO(block_mse_q8, (xd->bd - 8) * 2);
  }
  block_sse *= 16;
  // Tranform domain distortion is accurate for higher residuals.
  // TODO(any): Experiment with variance and mean based thresholds
  int use_transform_domain_distortion =
      (x->use_transform_domain_distortion > 0) &&
      (block_mse_q8 >= x->tx_domain_dist_threshold) &&
      // Any 64-pt transforms only preserves half the coefficients.
      // Therefore transform domain distortion is not valid for these
      // transform sizes.
      txsize_sqr_up_map[tx_size] != TX_64X64;
#if CONFIG_DIST_8X8
  if (x->using_dist_8x8) use_transform_domain_distortion = 0;
#endif
  int calc_pixel_domain_distortion_final =
      x->use_transform_domain_distortion == 1 &&
      use_transform_domain_distortion && x->rd_model != LOW_TXFM_RD;
  if (calc_pixel_domain_distortion_final &&
      (txk_allowed < TX_TYPES || allowed_tx_mask == 0x0001))
    calc_pixel_domain_distortion_final = use_transform_domain_distortion = 0;

  const uint16_t *eobs_ptr = x->plane[plane].eobs;

  // Used mse based threshold logic to take decision of R-D of optimization of
  // coeffs. For smaller residuals, coeff optimization would be helpful. For
  // larger residuals, R-D optimization may not be effective.
  // TODO(any): Experiment with variance and mean based thresholds
  perform_block_coeff_opt = (block_mse_q8 <= x->coeff_opt_dist_threshold);
  skip_trellis |= !perform_block_coeff_opt;

  assert(IMPLIES(txk_allowed < TX_TYPES, allowed_tx_mask == 1 << txk_allowed));

  TxfmParam txfm_param;
  QUANT_PARAM quant_param;
  av1_setup_xform(cm, x, tx_size, DCT_DCT, &txfm_param);
  av1_setup_quant(cm, tx_size, !skip_trellis,
                  skip_trellis ? (USE_B_QUANT_NO_TRELLIS ? AV1_XFORM_QUANT_B
                                                         : AV1_XFORM_QUANT_FP)
                               : AV1_XFORM_QUANT_FP,
                  &quant_param);
  int use_qm = !(xd->lossless[mbmi->segment_id] || cm->using_qmatrix == 0);

  for (int idx = 0; idx < TX_TYPES; ++idx) {
    const TX_TYPE tx_type = (TX_TYPE)txk_map[idx];
    if (!(allowed_tx_mask & (1 << tx_type))) continue;
    txfm_param.tx_type = tx_type;
    if (use_qm) {
      av1_setup_qmatrix(cm, x, plane, tx_size, tx_type, &quant_param);
    }
    if (plane == 0) xd->tx_type_map[tx_type_map_idx] = tx_type;
    RD_STATS this_rd_stats;
    av1_invalid_rd_stats(&this_rd_stats);

    av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, &txfm_param,
                    &quant_param);

    if (quant_param.use_optimize_b) {
      if (cpi->sf.rd_sf.optimize_b_precheck && best_rd < INT64_MAX &&
          eobs_ptr[block] >= 4) {
        // Calculate distortion quickly in transform domain.
        dist_block_tx_domain(x, plane, block, tx_size, &this_rd_stats.dist,
                             &this_rd_stats.sse);

        const int64_t best_rd_ = AOMMIN(best_rd, ref_best_rd);
        const int64_t dist_cost_estimate =
            RDCOST(x->rdmult, 0, AOMMIN(this_rd_stats.dist, this_rd_stats.sse));
        if (dist_cost_estimate - (dist_cost_estimate >> 3) > best_rd_) continue;
      }
      av1_optimize_b(cpi, x, plane, block, tx_size, tx_type, txb_ctx,
                     cpi->sf.rd_sf.trellis_eob_fast, &rate_cost);
    } else {
      rate_cost =
          av1_cost_coeffs(x, plane, block, tx_size, tx_type, txb_ctx,
                          use_fast_coef_costing, cm->reduced_tx_set_used);
    }

    // If rd cost based on coeff rate is more than best_rd, skip the calculation
    // of distortion
    int64_t tmp_rd = RDCOST(x->rdmult, rate_cost, 0);
    if (tmp_rd > best_rd) continue;
    if (eobs_ptr[block] == 0) {
      // When eob is 0, pixel domain distortion is more efficient and accurate.
      this_rd_stats.dist = this_rd_stats.sse = block_sse;
    } else if (use_transform_domain_distortion) {
      dist_block_tx_domain(x, plane, block, tx_size, &this_rd_stats.dist,
                           &this_rd_stats.sse);
    } else {
      int64_t sse_diff = INT64_MAX;
      // high_energy threshold assumes that every pixel within a txfm block
      // has a residue energy of at least 25% of the maximum, i.e. 128 * 128
      // for 8 bit, then the threshold is scaled based on input bit depth.
      const int64_t high_energy_thresh =
          ((int64_t)128 * 128 * tx_size_2d[tx_size]) << ((xd->bd - 8) * 2);
      const int is_high_energy = (block_sse >= high_energy_thresh);
      if (tx_size == TX_64X64 || is_high_energy) {
        // Because 3 out 4 quadrants of transform coefficients are forced to
        // zero, the inverse transform has a tendency to overflow. sse_diff
        // is effectively the energy of those 3 quadrants, here we use it
        // to decide if we should do pixel domain distortion. If the energy
        // is mostly in first quadrant, then it is unlikely that we have
        // overflow issue in inverse transform.
        dist_block_tx_domain(x, plane, block, tx_size, &this_rd_stats.dist,
                             &this_rd_stats.sse);
        sse_diff = block_sse - this_rd_stats.sse;
      }
      if (tx_size != TX_64X64 || !is_high_energy ||
          (sse_diff * 2) < this_rd_stats.sse) {
        const int64_t tx_domain_dist = this_rd_stats.dist;
        this_rd_stats.dist = dist_block_px_domain(
            cpi, x, plane, plane_bsize, block, blk_row, blk_col, tx_size);
        // For high energy blocks, occasionally, the pixel domain distortion
        // can be artificially low due to clamping at reconstruction stage
        // even when inverse transform output is hugely different from the
        // actual residue.
        if (is_high_energy && this_rd_stats.dist < tx_domain_dist)
          this_rd_stats.dist = tx_domain_dist;
      } else {
        this_rd_stats.dist += sse_diff;
      }
      this_rd_stats.sse = block_sse;
    }

    this_rd_stats.rate = rate_cost;

    const int64_t rd =
        RDCOST(x->rdmult, this_rd_stats.rate, this_rd_stats.dist);

    if (rd < best_rd) {
      best_rd = rd;
      *best_rd_stats = this_rd_stats;
      best_tx_type = tx_type;
      best_txb_ctx = x->plane[plane].txb_entropy_ctx[block];
      best_eob = x->plane[plane].eobs[block];
      last_tx_type = best_tx_type;

      // Swap qcoeff and dqcoeff buffers
      tran_low_t *const tmp_dqcoeff = best_dqcoeff;
      best_dqcoeff = pd->dqcoeff;
      pd->dqcoeff = tmp_dqcoeff;
    }

#if CONFIG_COLLECT_RD_STATS == 1
    if (plane == 0) {
      PrintTransformUnitStats(cpi, x, &this_rd_stats, blk_row, blk_col,
                              plane_bsize, tx_size, tx_type, rd);
    }
#endif  // CONFIG_COLLECT_RD_STATS == 1

#if COLLECT_TX_SIZE_DATA
    // Generate small sample to restrict output size.
    static unsigned int seed = 21743;
    if (lcg_rand16(&seed) % 200 == 0) {
      FILE *fp = NULL;

      if (within_border) {
        fp = fopen(av1_tx_size_data_output_file, "a");
      }

      if (fp) {
        // Transform info and RD
        const int txb_w = tx_size_wide[tx_size];
        const int txb_h = tx_size_high[tx_size];

        // Residue signal.
        const int diff_stride = block_size_wide[plane_bsize];
        struct macroblock_plane *const p = &x->plane[plane];
        const int16_t *src_diff =
            &p->src_diff[(blk_row * diff_stride + blk_col) * 4];

        for (int r = 0; r < txb_h; ++r) {
          for (int c = 0; c < txb_w; ++c) {
            fprintf(fp, "%d,", src_diff[c]);
          }
          src_diff += diff_stride;
        }

        fprintf(fp, "%d,%d,%d,%" PRId64, txb_w, txb_h, tx_type, rd);
        fprintf(fp, "\n");
        fclose(fp);
      }
    }
#endif  // COLLECT_TX_SIZE_DATA

    if (cpi->sf.tx_sf.adaptive_txb_search_level) {
      if ((best_rd - (best_rd >> cpi->sf.tx_sf.adaptive_txb_search_level)) >
          ref_best_rd) {
        break;
      }
    }

    // Skip transform type search when we found the block has been quantized to
    // all zero and at the same time, it has better rdcost than doing transform.
    if (cpi->sf.tx_sf.tx_type_search.skip_tx_search && !best_eob) break;
  }

  assert(best_rd != INT64_MAX);

  best_rd_stats->skip = best_eob == 0;
  if (plane == 0) update_txk_array(xd, blk_row, blk_col, tx_size, best_tx_type);
  x->plane[plane].txb_entropy_ctx[block] = best_txb_ctx;
  x->plane[plane].eobs[block] = best_eob;

  pd->dqcoeff = best_dqcoeff;

  if (calc_pixel_domain_distortion_final && best_eob) {
    best_rd_stats->dist = dist_block_px_domain(
        cpi, x, plane, plane_bsize, block, blk_row, blk_col, tx_size);
    best_rd_stats->sse = block_sse;
  }

  if (intra_txb_rd_info != NULL) {
    intra_txb_rd_info->valid = 1;
    intra_txb_rd_info->entropy_context = cur_joint_ctx;
    intra_txb_rd_info->rate = best_rd_stats->rate;
    intra_txb_rd_info->dist = best_rd_stats->dist;
    intra_txb_rd_info->sse = best_rd_stats->sse;
    intra_txb_rd_info->eob = best_eob;
    intra_txb_rd_info->txb_entropy_ctx = best_txb_ctx;
    intra_txb_rd_info->perform_block_coeff_opt = perform_block_coeff_opt;
    if (plane == 0) intra_txb_rd_info->tx_type = best_tx_type;
  }

RECON_INTRA:
  if (!is_inter && best_eob &&
      (blk_row + tx_size_high_unit[tx_size] < mi_size_high[plane_bsize] ||
       blk_col + tx_size_wide_unit[tx_size] < mi_size_wide[plane_bsize])) {
    // intra mode needs decoded result such that the next transform block
    // can use it for prediction.
    // if the last search tx_type is the best tx_type, we don't need to
    // do this again
    if (best_tx_type != last_tx_type) {
      TxfmParam txfm_param_intra;
      QUANT_PARAM quant_param_intra;
      av1_setup_xform(cm, x, tx_size, best_tx_type, &txfm_param_intra);
      av1_setup_quant(cm, tx_size, !skip_trellis,
                      skip_trellis
                          ? (USE_B_QUANT_NO_TRELLIS ? AV1_XFORM_QUANT_B
                                                    : AV1_XFORM_QUANT_FP)
                          : AV1_XFORM_QUANT_FP,
                      &quant_param_intra);
      av1_setup_qmatrix(cm, x, plane, tx_size, best_tx_type,
                        &quant_param_intra);
      av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                      &txfm_param_intra, &quant_param_intra);
      if (quant_param_intra.use_optimize_b) {
        av1_optimize_b(cpi, x, plane, block, tx_size, best_tx_type, txb_ctx,
                       cpi->sf.rd_sf.trellis_eob_fast, &rate_cost);
      }
    }

    inverse_transform_block_facade(xd, plane, block, blk_row, blk_col,
                                   x->plane[plane].eobs[block],
                                   cm->reduced_tx_set_used);

    // This may happen because of hash collision. The eob stored in the hash
    // table is non-zero, but the real eob is zero. We need to make sure tx_type
    // is DCT_DCT in this case.
    if (plane == 0 && x->plane[plane].eobs[block] == 0 &&
        best_tx_type != DCT_DCT) {
      update_txk_array(xd, blk_row, blk_col, tx_size, DCT_DCT);
    }
  }
  pd->dqcoeff = orig_dqcoeff;

  return best_rd;
}

// Pick transform type for a transform block of tx_size.
static AOM_INLINE void tx_type_rd(const AV1_COMP *cpi, MACROBLOCK *x,
                                  TX_SIZE tx_size, int blk_row, int blk_col,
                                  int plane, int block, int plane_bsize,
                                  TXB_CTX *txb_ctx, RD_STATS *rd_stats,
                                  FAST_TX_SEARCH_MODE ftxs_mode,
                                  int64_t ref_rdcost,
                                  TXB_RD_INFO *rd_info_array) {
  const struct macroblock_plane *const p = &x->plane[plane];
  const uint16_t cur_joint_ctx =
      (txb_ctx->dc_sign_ctx << 8) + txb_ctx->txb_skip_ctx;
  MACROBLOCKD *xd = &x->e_mbd;
  const int tx_type_map_idx =
      plane ? 0 : blk_row * xd->tx_type_map_stride + blk_col;
  // Look up RD and terminate early in case when we've already processed exactly
  // the same residual with exactly the same entropy context.
  if (rd_info_array != NULL && rd_info_array->valid &&
      rd_info_array->entropy_context == cur_joint_ctx) {
    if (plane == 0) xd->tx_type_map[tx_type_map_idx] = rd_info_array->tx_type;
    const TX_TYPE ref_tx_type =
        av1_get_tx_type(&x->e_mbd, get_plane_type(plane), blk_row, blk_col,
                        tx_size, cpi->common.reduced_tx_set_used);
    if (ref_tx_type == rd_info_array->tx_type) {
      rd_stats->rate += rd_info_array->rate;
      rd_stats->dist += rd_info_array->dist;
      rd_stats->sse += rd_info_array->sse;
      rd_stats->skip &= rd_info_array->eob == 0;
      p->eobs[block] = rd_info_array->eob;
      p->txb_entropy_ctx[block] = rd_info_array->txb_entropy_ctx;
      return;
    }
  }

  RD_STATS this_rd_stats;
  search_txk_type(cpi, x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                  txb_ctx, ftxs_mode, 0, 0, ref_rdcost, &this_rd_stats);

  av1_merge_rd_stats(rd_stats, &this_rd_stats);

  // Save RD results for possible reuse in future.
  if (rd_info_array != NULL) {
    rd_info_array->valid = 1;
    rd_info_array->entropy_context = cur_joint_ctx;
    rd_info_array->rate = this_rd_stats.rate;
    rd_info_array->dist = this_rd_stats.dist;
    rd_info_array->sse = this_rd_stats.sse;
    rd_info_array->eob = p->eobs[block];
    rd_info_array->txb_entropy_ctx = p->txb_entropy_ctx[block];
    if (plane == 0) rd_info_array->tx_type = xd->tx_type_map[tx_type_map_idx];
  }
}

// Finds rd cost for a y block, given the transform size partitions
static void tx_block_yrd(const AV1_COMP *cpi, MACROBLOCK *x, int blk_row,
                         int blk_col, int block, TX_SIZE tx_size,
                         BLOCK_SIZE plane_bsize, int depth,
                         ENTROPY_CONTEXT *above_ctx, ENTROPY_CONTEXT *left_ctx,
                         TXFM_CONTEXT *tx_above, TXFM_CONTEXT *tx_left,
                         int64_t ref_best_rd, RD_STATS *rd_stats,
                         FAST_TX_SEARCH_MODE ftxs_mode) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int max_blocks_high = max_block_high(xd, plane_bsize, 0);
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, 0);

  assert(tx_size < TX_SIZES_ALL);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  const TX_SIZE plane_tx_size = mbmi->inter_tx_size[av1_get_txb_size_index(
      plane_bsize, blk_row, blk_col)];

  int ctx = txfm_partition_context(tx_above + blk_col, tx_left + blk_row,
                                   mbmi->sb_type, tx_size);

  av1_init_rd_stats(rd_stats);
  if (tx_size == plane_tx_size) {
    ENTROPY_CONTEXT *ta = above_ctx + blk_col;
    ENTROPY_CONTEXT *tl = left_ctx + blk_row;
    const TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);
    TXB_CTX txb_ctx;
    get_txb_ctx(plane_bsize, tx_size, 0, ta, tl, &txb_ctx);

    const int zero_blk_rate = x->coeff_costs[txs_ctx][get_plane_type(0)]
                                  .txb_skip_cost[txb_ctx.txb_skip_ctx][1];
    rd_stats->zero_rate = zero_blk_rate;
    tx_type_rd(cpi, x, tx_size, blk_row, blk_col, 0, block, plane_bsize,
               &txb_ctx, rd_stats, ftxs_mode, ref_best_rd, NULL);
    const int mi_width = mi_size_wide[plane_bsize];
    if (RDCOST(x->rdmult, rd_stats->rate, rd_stats->dist) >=
            RDCOST(x->rdmult, zero_blk_rate, rd_stats->sse) ||
        rd_stats->skip == 1) {
      rd_stats->rate = zero_blk_rate;
      rd_stats->dist = rd_stats->sse;
      rd_stats->skip = 1;
      set_blk_skip(x, 0, blk_row * mi_width + blk_col, 1);
      x->plane[0].eobs[block] = 0;
      x->plane[0].txb_entropy_ctx[block] = 0;
      update_txk_array(xd, blk_row, blk_col, tx_size, DCT_DCT);
    } else {
      rd_stats->skip = 0;
      set_blk_skip(x, 0, blk_row * mi_width + blk_col, 0);
    }
    if (tx_size > TX_4X4 && depth < MAX_VARTX_DEPTH)
      rd_stats->rate += x->txfm_partition_cost[ctx][0];
    av1_set_txb_context(x, 0, block, tx_size, ta, tl);
    txfm_partition_update(tx_above + blk_col, tx_left + blk_row, tx_size,
                          tx_size);
  } else {
    const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
    const int bsw = tx_size_wide_unit[sub_txs];
    const int bsh = tx_size_high_unit[sub_txs];
    const int step = bsh * bsw;
    RD_STATS pn_rd_stats;
    int64_t this_rd = 0;
    assert(bsw > 0 && bsh > 0);

    for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh) {
      for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
        const int offsetr = blk_row + row;
        const int offsetc = blk_col + col;

        if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide) continue;

        av1_init_rd_stats(&pn_rd_stats);
        tx_block_yrd(cpi, x, offsetr, offsetc, block, sub_txs, plane_bsize,
                     depth + 1, above_ctx, left_ctx, tx_above, tx_left,
                     ref_best_rd - this_rd, &pn_rd_stats, ftxs_mode);
        if (pn_rd_stats.rate == INT_MAX) {
          av1_invalid_rd_stats(rd_stats);
          return;
        }
        av1_merge_rd_stats(rd_stats, &pn_rd_stats);
        this_rd += RDCOST(x->rdmult, pn_rd_stats.rate, pn_rd_stats.dist);
        block += step;
      }
    }

    if (tx_size > TX_4X4 && depth < MAX_VARTX_DEPTH)
      rd_stats->rate += x->txfm_partition_cost[ctx][1];
  }
}

// Return value 0: early termination triggered, no valid rd cost available;
//              1: rd cost values are valid.
static int inter_block_yrd(const AV1_COMP *cpi, MACROBLOCK *x,
                           RD_STATS *rd_stats, BLOCK_SIZE bsize,
                           int64_t ref_best_rd, FAST_TX_SEARCH_MODE ftxs_mode) {
  MACROBLOCKD *const xd = &x->e_mbd;
  int is_cost_valid = 1;
  int64_t this_rd = 0;

  if (ref_best_rd < 0) is_cost_valid = 0;

  av1_init_rd_stats(rd_stats);

  if (is_cost_valid) {
    const struct macroblockd_plane *const pd = &xd->plane[0];
    const BLOCK_SIZE plane_bsize =
        get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
    const int mi_width = mi_size_wide[plane_bsize];
    const int mi_height = mi_size_high[plane_bsize];
    const TX_SIZE max_tx_size = get_vartx_max_txsize(xd, plane_bsize, 0);
    const int bh = tx_size_high_unit[max_tx_size];
    const int bw = tx_size_wide_unit[max_tx_size];
    const int init_depth = get_search_init_depth(
        mi_width, mi_height, 1, &cpi->sf, x->tx_size_search_method);
    int idx, idy;
    int block = 0;
    int step = tx_size_wide_unit[max_tx_size] * tx_size_high_unit[max_tx_size];
    ENTROPY_CONTEXT ctxa[MAX_MIB_SIZE];
    ENTROPY_CONTEXT ctxl[MAX_MIB_SIZE];
    TXFM_CONTEXT tx_above[MAX_MIB_SIZE];
    TXFM_CONTEXT tx_left[MAX_MIB_SIZE];
    RD_STATS pn_rd_stats;

    av1_get_entropy_contexts(bsize, pd, ctxa, ctxl);
    memcpy(tx_above, xd->above_txfm_context, sizeof(TXFM_CONTEXT) * mi_width);
    memcpy(tx_left, xd->left_txfm_context, sizeof(TXFM_CONTEXT) * mi_height);

    for (idy = 0; idy < mi_height; idy += bh) {
      for (idx = 0; idx < mi_width; idx += bw) {
        av1_init_rd_stats(&pn_rd_stats);
        tx_block_yrd(cpi, x, idy, idx, block, max_tx_size, plane_bsize,
                     init_depth, ctxa, ctxl, tx_above, tx_left,
                     ref_best_rd - this_rd, &pn_rd_stats, ftxs_mode);
        if (pn_rd_stats.rate == INT_MAX) {
          av1_invalid_rd_stats(rd_stats);
          return 0;
        }
        av1_merge_rd_stats(rd_stats, &pn_rd_stats);
        this_rd +=
            AOMMIN(RDCOST(x->rdmult, pn_rd_stats.rate, pn_rd_stats.dist),
                   RDCOST(x->rdmult, pn_rd_stats.zero_rate, pn_rd_stats.sse));
        block += step;
      }
    }
  }

  const int skip_ctx = av1_get_skip_context(xd);
  const int s0 = x->skip_cost[skip_ctx][0];
  const int s1 = x->skip_cost[skip_ctx][1];
  int64_t skip_rd = RDCOST(x->rdmult, s1, rd_stats->sse);
  this_rd = RDCOST(x->rdmult, rd_stats->rate + s0, rd_stats->dist);
  if (skip_rd < this_rd) {
    this_rd = skip_rd;
    rd_stats->rate = 0;
    rd_stats->dist = rd_stats->sse;
    rd_stats->skip = 1;
  }
  if (this_rd > ref_best_rd) is_cost_valid = 0;

  if (!is_cost_valid) {
    // reset cost value
    av1_invalid_rd_stats(rd_stats);
  }
  return is_cost_valid;
}

static AOM_INLINE void select_tx_block(
    const AV1_COMP *cpi, MACROBLOCK *x, int blk_row, int blk_col, int block,
    TX_SIZE tx_size, int depth, BLOCK_SIZE plane_bsize, ENTROPY_CONTEXT *ta,
    ENTROPY_CONTEXT *tl, TXFM_CONTEXT *tx_above, TXFM_CONTEXT *tx_left,
    RD_STATS *rd_stats, int64_t prev_level_rd, int64_t ref_best_rd,
    int *is_cost_valid, FAST_TX_SEARCH_MODE ftxs_mode,
    TXB_RD_INFO_NODE *rd_info_node);

static int64_t select_tx_size_and_type(const AV1_COMP *cpi, MACROBLOCK *x,
                                       RD_STATS *rd_stats, BLOCK_SIZE bsize,
                                       int64_t ref_best_rd,
                                       TXB_RD_INFO_NODE *rd_info_tree) {
  MACROBLOCKD *const xd = &x->e_mbd;
  assert(is_inter_block(xd->mi[0]));
  assert(bsize < BLOCK_SIZES_ALL);

  // TODO(debargha): enable this as a speed feature where the
  // select_inter_block_yrd() function above will use a simplified search
  // such as not using full optimize, but the inter_block_yrd() function
  // will use more complex search given that the transform partitions have
  // already been decided.

  const int fast_tx_search = x->tx_size_search_method > USE_FULL_RD;
  int64_t rd_thresh = ref_best_rd;
  if (fast_tx_search && rd_thresh < INT64_MAX) {
    if (INT64_MAX - rd_thresh > (rd_thresh >> 3)) rd_thresh += (rd_thresh >> 3);
  }
  assert(rd_thresh > 0);

  const FAST_TX_SEARCH_MODE ftxs_mode =
      fast_tx_search ? FTXS_DCT_AND_1D_DCT_ONLY : FTXS_NONE;
  const struct macroblockd_plane *const pd = &xd->plane[0];
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
  assert(plane_bsize < BLOCK_SIZES_ALL);
  const int mi_width = mi_size_wide[plane_bsize];
  const int mi_height = mi_size_high[plane_bsize];
  ENTROPY_CONTEXT ctxa[MAX_MIB_SIZE];
  ENTROPY_CONTEXT ctxl[MAX_MIB_SIZE];
  TXFM_CONTEXT tx_above[MAX_MIB_SIZE];
  TXFM_CONTEXT tx_left[MAX_MIB_SIZE];
  av1_get_entropy_contexts(bsize, pd, ctxa, ctxl);
  memcpy(tx_above, xd->above_txfm_context, sizeof(TXFM_CONTEXT) * mi_width);
  memcpy(tx_left, xd->left_txfm_context, sizeof(TXFM_CONTEXT) * mi_height);

  const int skip_ctx = av1_get_skip_context(xd);
  const int s0 = x->skip_cost[skip_ctx][0];
  const int s1 = x->skip_cost[skip_ctx][1];
  const int init_depth = get_search_init_depth(mi_width, mi_height, 1, &cpi->sf,
                                               x->tx_size_search_method);
  const TX_SIZE max_tx_size = max_txsize_rect_lookup[plane_bsize];
  const int bh = tx_size_high_unit[max_tx_size];
  const int bw = tx_size_wide_unit[max_tx_size];
  const int step = bw * bh;
  int64_t skip_rd = RDCOST(x->rdmult, s1, 0);
  int64_t this_rd = RDCOST(x->rdmult, s0, 0);
  int block = 0;

  av1_init_rd_stats(rd_stats);
  for (int idy = 0; idy < mi_height; idy += bh) {
    for (int idx = 0; idx < mi_width; idx += bw) {
      const int64_t best_rd_sofar =
          (rd_thresh == INT64_MAX) ? INT64_MAX
                                   : (rd_thresh - (AOMMIN(skip_rd, this_rd)));
      int is_cost_valid = 1;
      RD_STATS pn_rd_stats;
      select_tx_block(cpi, x, idy, idx, block, max_tx_size, init_depth,
                      plane_bsize, ctxa, ctxl, tx_above, tx_left, &pn_rd_stats,
                      INT64_MAX, best_rd_sofar, &is_cost_valid, ftxs_mode,
                      rd_info_tree);
      if (!is_cost_valid || pn_rd_stats.rate == INT_MAX) {
        av1_invalid_rd_stats(rd_stats);
        return INT64_MAX;
      }
      av1_merge_rd_stats(rd_stats, &pn_rd_stats);
      skip_rd = RDCOST(x->rdmult, s1, rd_stats->sse);
      this_rd = RDCOST(x->rdmult, rd_stats->rate + s0, rd_stats->dist);
      block += step;
      if (rd_info_tree != NULL) rd_info_tree += 1;
    }
  }

  if (skip_rd <= this_rd) {
    rd_stats->skip = 1;
  } else {
    rd_stats->skip = 0;
  }

  if (rd_stats->rate == INT_MAX) return INT64_MAX;

  // If fast_tx_search is true, only DCT and 1D DCT were tested in
  // select_inter_block_yrd() above. Do a better search for tx type with
  // tx sizes already decided.
  if (fast_tx_search) {
    if (!inter_block_yrd(cpi, x, rd_stats, bsize, ref_best_rd, FTXS_NONE))
      return INT64_MAX;
  }

  int64_t rd;
  if (rd_stats->skip) {
    rd = RDCOST(x->rdmult, s1, rd_stats->sse);
  } else {
    rd = RDCOST(x->rdmult, rd_stats->rate + s0, rd_stats->dist);
    if (!xd->lossless[xd->mi[0]->segment_id])
      rd = AOMMIN(rd, RDCOST(x->rdmult, s1, rd_stats->sse));
  }

  return rd;
}

// Search for best transform size and type for luma inter blocks.
void pick_tx_size_type_yrd(const AV1_COMP *cpi, MACROBLOCK *x,
                           RD_STATS *rd_stats, BLOCK_SIZE bsize,
                           int64_t ref_best_rd) {
  const AV1_COMMON *cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  assert(is_inter_block(xd->mi[0]));

  av1_invalid_rd_stats(rd_stats);

  if (cpi->sf.tx_sf.model_based_prune_tx_search_level &&
      ref_best_rd != INT64_MAX) {
    int model_rate;
    int64_t model_dist;
    int model_skip;
    model_rd_sb_fn[MODELRD_TYPE_TX_SEARCH_PRUNE](
        cpi, bsize, x, xd, 0, 0, &model_rate, &model_dist, &model_skip, NULL,
        NULL, NULL, NULL);
    const int64_t model_rd = RDCOST(x->rdmult, model_rate, model_dist);
    // If the modeled rd is a lot worse than the best so far, breakout.
    // TODO(debargha, urvang): Improve the model and make the check below
    // tighter.
    assert(cpi->sf.tx_sf.model_based_prune_tx_search_level >= 0 &&
           cpi->sf.tx_sf.model_based_prune_tx_search_level <= 2);
    static const int prune_factor_by8[] = { 3, 5 };
    if (!model_skip &&
        ((model_rd *
          prune_factor_by8[cpi->sf.tx_sf.model_based_prune_tx_search_level -
                           1]) >>
         3) > ref_best_rd)
      return;
  }

  uint32_t hash = 0;
  int32_t match_index = -1;
  MB_RD_RECORD *mb_rd_record = NULL;
  const int mi_row = x->e_mbd.mi_row;
  const int mi_col = x->e_mbd.mi_col;
  const int within_border =
      mi_row >= xd->tile.mi_row_start &&
      (mi_row + mi_size_high[bsize] < xd->tile.mi_row_end) &&
      mi_col >= xd->tile.mi_col_start &&
      (mi_col + mi_size_wide[bsize] < xd->tile.mi_col_end);
  const int is_mb_rd_hash_enabled =
      (within_border && cpi->sf.rd_sf.use_mb_rd_hash);
  const int n4 = bsize_to_num_blk(bsize);
  if (is_mb_rd_hash_enabled) {
    hash = get_block_residue_hash(x, bsize);
    mb_rd_record = &x->mb_rd_record;
    match_index = find_mb_rd_info(mb_rd_record, ref_best_rd, hash);
    if (match_index != -1) {
      MB_RD_INFO *tx_rd_info = &mb_rd_record->tx_rd_info[match_index];
      fetch_tx_rd_info(n4, tx_rd_info, rd_stats, x);
      return;
    }
  }

  // If we predict that skip is the optimal RD decision - set the respective
  // context and terminate early.
  int64_t dist;
  if (x->predict_skip_level &&
      predict_skip_flag(x, bsize, &dist, cm->reduced_tx_set_used)) {
    set_skip_flag(x, rd_stats, bsize, dist);
    // Save the RD search results into tx_rd_record.
    if (is_mb_rd_hash_enabled)
      save_tx_rd_info(n4, hash, x, rd_stats, mb_rd_record);
    return;
  }
#if CONFIG_SPEED_STATS
  ++x->tx_search_count;
#endif  // CONFIG_SPEED_STATS

  // Precompute residual hashes and find existing or add new RD records to
  // store and reuse rate and distortion values to speed up TX size search.
  TXB_RD_INFO_NODE matched_rd_info[4 + 16 + 64];
  int found_rd_info = 0;
  if (ref_best_rd != INT64_MAX && within_border &&
      cpi->sf.tx_sf.use_inter_txb_hash) {
    found_rd_info = find_tx_size_rd_records(x, bsize, matched_rd_info);
  }

  int found = 0;
  RD_STATS this_rd_stats;
  av1_init_rd_stats(&this_rd_stats);
  const int64_t rd =
      select_tx_size_and_type(cpi, x, &this_rd_stats, bsize, ref_best_rd,
                              found_rd_info ? matched_rd_info : NULL);

  if (rd < INT64_MAX) {
    *rd_stats = this_rd_stats;
    found = 1;
  }

  // We should always find at least one candidate unless ref_best_rd is less
  // than INT64_MAX (in which case, all the calls to select_tx_size_fix_type
  // might have failed to find something better)
  assert(IMPLIES(!found, ref_best_rd != INT64_MAX));
  if (!found) return;

  // Save the RD search results into tx_rd_record.
  if (is_mb_rd_hash_enabled) {
    assert(mb_rd_record != NULL);
    save_tx_rd_info(n4, hash, x, rd_stats, mb_rd_record);
  }
}

static AOM_INLINE void choose_largest_tx_size(const AV1_COMP *const cpi,
                                              MACROBLOCK *x, RD_STATS *rd_stats,
                                              int64_t ref_best_rd,
                                              BLOCK_SIZE bs) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  mbmi->tx_size = tx_size_from_tx_mode(bs, x->tx_mode_search_type);

  // If tx64 is not enabled, we need to go down to the next available size
  if (!cpi->oxcf.enable_tx64) {
    static const TX_SIZE tx_size_max_32[TX_SIZES_ALL] = {
      TX_4X4,    // 4x4 transform
      TX_8X8,    // 8x8 transform
      TX_16X16,  // 16x16 transform
      TX_32X32,  // 32x32 transform
      TX_32X32,  // 64x64 transform
      TX_4X8,    // 4x8 transform
      TX_8X4,    // 8x4 transform
      TX_8X16,   // 8x16 transform
      TX_16X8,   // 16x8 transform
      TX_16X32,  // 16x32 transform
      TX_32X16,  // 32x16 transform
      TX_32X32,  // 32x64 transform
      TX_32X32,  // 64x32 transform
      TX_4X16,   // 4x16 transform
      TX_16X4,   // 16x4 transform
      TX_8X32,   // 8x32 transform
      TX_32X8,   // 32x8 transform
      TX_16X32,  // 16x64 transform
      TX_32X16,  // 64x16 transform
    };

    mbmi->tx_size = tx_size_max_32[mbmi->tx_size];
  }

  const int skip_ctx = av1_get_skip_context(xd);
  int s0, s1;

  s0 = x->skip_cost[skip_ctx][0];
  s1 = x->skip_cost[skip_ctx][1];

  int64_t skip_rd = INT64_MAX;
  int64_t this_rd = RDCOST(x->rdmult, s0, 0);

  // Skip RDcost is used only for Inter blocks
  if (is_inter_block(xd->mi[0])) skip_rd = RDCOST(x->rdmult, s1, 0);

  txfm_rd_in_plane(x, cpi, rd_stats, ref_best_rd, AOMMIN(this_rd, skip_rd),
                   AOM_PLANE_Y, bs, mbmi->tx_size,
                   cpi->sf.rd_sf.use_fast_coef_costing, FTXS_NONE, 0);
}

static AOM_INLINE void choose_smallest_tx_size(const AV1_COMP *const cpi,
                                               MACROBLOCK *x,
                                               RD_STATS *rd_stats,
                                               int64_t ref_best_rd,
                                               BLOCK_SIZE bs) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];

  mbmi->tx_size = TX_4X4;
  // TODO(any) : Pass this_rd based on skip/non-skip cost
  txfm_rd_in_plane(x, cpi, rd_stats, ref_best_rd, 0, 0, bs, mbmi->tx_size,
                   cpi->sf.rd_sf.use_fast_coef_costing, FTXS_NONE, 0);
}

static AOM_INLINE void try_tx_block_no_split(
    const AV1_COMP *cpi, MACROBLOCK *x, int blk_row, int blk_col, int block,
    TX_SIZE tx_size, int depth, BLOCK_SIZE plane_bsize,
    const ENTROPY_CONTEXT *ta, const ENTROPY_CONTEXT *tl,
    int txfm_partition_ctx, RD_STATS *rd_stats, int64_t ref_best_rd,
    FAST_TX_SEARCH_MODE ftxs_mode, TXB_RD_INFO_NODE *rd_info_node,
    TxCandidateInfo *no_split) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  struct macroblock_plane *const p = &x->plane[0];
  const int bw = mi_size_wide[plane_bsize];

  no_split->rd = INT64_MAX;
  no_split->txb_entropy_ctx = 0;
  no_split->tx_type = TX_TYPES;

  const ENTROPY_CONTEXT *const pta = ta + blk_col;
  const ENTROPY_CONTEXT *const ptl = tl + blk_row;

  const TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);
  TXB_CTX txb_ctx;
  get_txb_ctx(plane_bsize, tx_size, 0, pta, ptl, &txb_ctx);
  const int zero_blk_rate = x->coeff_costs[txs_ctx][PLANE_TYPE_Y]
                                .txb_skip_cost[txb_ctx.txb_skip_ctx][1];
  rd_stats->zero_rate = zero_blk_rate;
  const int index = av1_get_txb_size_index(plane_bsize, blk_row, blk_col);
  mbmi->inter_tx_size[index] = tx_size;
  tx_type_rd(cpi, x, tx_size, blk_row, blk_col, 0, block, plane_bsize, &txb_ctx,
             rd_stats, ftxs_mode, ref_best_rd,
             rd_info_node != NULL ? rd_info_node->rd_info_array : NULL);
  assert(rd_stats->rate < INT_MAX);

  if ((RDCOST(x->rdmult, rd_stats->rate, rd_stats->dist) >=
           RDCOST(x->rdmult, zero_blk_rate, rd_stats->sse) ||
       rd_stats->skip == 1) &&
      !xd->lossless[mbmi->segment_id]) {
#if CONFIG_RD_DEBUG
    av1_update_txb_coeff_cost(rd_stats, 0, tx_size, blk_row, blk_col,
                              zero_blk_rate - rd_stats->rate);
#endif  // CONFIG_RD_DEBUG
    rd_stats->rate = zero_blk_rate;
    rd_stats->dist = rd_stats->sse;
    rd_stats->skip = 1;
    set_blk_skip(x, 0, blk_row * bw + blk_col, 1);
    p->eobs[block] = 0;
    update_txk_array(xd, blk_row, blk_col, tx_size, DCT_DCT);
  } else {
    set_blk_skip(x, 0, blk_row * bw + blk_col, 0);
    rd_stats->skip = 0;
  }

  if (tx_size > TX_4X4 && depth < MAX_VARTX_DEPTH)
    rd_stats->rate += x->txfm_partition_cost[txfm_partition_ctx][0];

  no_split->rd = RDCOST(x->rdmult, rd_stats->rate, rd_stats->dist);
  no_split->txb_entropy_ctx = p->txb_entropy_ctx[block];
  no_split->tx_type =
      xd->tx_type_map[blk_row * xd->tx_type_map_stride + blk_col];
}

static AOM_INLINE void try_tx_block_split(
    const AV1_COMP *cpi, MACROBLOCK *x, int blk_row, int blk_col, int block,
    TX_SIZE tx_size, int depth, BLOCK_SIZE plane_bsize, ENTROPY_CONTEXT *ta,
    ENTROPY_CONTEXT *tl, TXFM_CONTEXT *tx_above, TXFM_CONTEXT *tx_left,
    int txfm_partition_ctx, int64_t no_split_rd, int64_t ref_best_rd,
    FAST_TX_SEARCH_MODE ftxs_mode, TXB_RD_INFO_NODE *rd_info_node,
    RD_STATS *split_rd_stats, int64_t *split_rd) {
  assert(tx_size < TX_SIZES_ALL);
  MACROBLOCKD *const xd = &x->e_mbd;
  const int max_blocks_high = max_block_high(xd, plane_bsize, 0);
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, 0);
  const TX_SIZE sub_txs = sub_tx_size_map[tx_size];
  const int bsw = tx_size_wide_unit[sub_txs];
  const int bsh = tx_size_high_unit[sub_txs];
  const int sub_step = bsw * bsh;
  const int nblks =
      (tx_size_high_unit[tx_size] / bsh) * (tx_size_wide_unit[tx_size] / bsw);
  assert(nblks > 0);
  int blk_idx = 0;
  int64_t tmp_rd = 0;
  *split_rd = INT64_MAX;
  split_rd_stats->rate = x->txfm_partition_cost[txfm_partition_ctx][1];

  for (int r = 0; r < tx_size_high_unit[tx_size]; r += bsh) {
    for (int c = 0; c < tx_size_wide_unit[tx_size]; c += bsw, ++blk_idx) {
      assert(blk_idx < 4);
      const int offsetr = blk_row + r;
      const int offsetc = blk_col + c;
      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide) continue;

      RD_STATS this_rd_stats;
      int this_cost_valid = 1;
      select_tx_block(
          cpi, x, offsetr, offsetc, block, sub_txs, depth + 1, plane_bsize, ta,
          tl, tx_above, tx_left, &this_rd_stats, no_split_rd / nblks,
          ref_best_rd - tmp_rd, &this_cost_valid, ftxs_mode,
          (rd_info_node != NULL) ? rd_info_node->children[blk_idx] : NULL);
      if (!this_cost_valid) return;
      av1_merge_rd_stats(split_rd_stats, &this_rd_stats);
      tmp_rd = RDCOST(x->rdmult, split_rd_stats->rate, split_rd_stats->dist);
      if (no_split_rd < tmp_rd) return;
      block += sub_step;
    }
  }

  *split_rd = tmp_rd;
}

// Search for the best tx partition/type for a given luma block.
static AOM_INLINE void select_tx_block(
    const AV1_COMP *cpi, MACROBLOCK *x, int blk_row, int blk_col, int block,
    TX_SIZE tx_size, int depth, BLOCK_SIZE plane_bsize, ENTROPY_CONTEXT *ta,
    ENTROPY_CONTEXT *tl, TXFM_CONTEXT *tx_above, TXFM_CONTEXT *tx_left,
    RD_STATS *rd_stats, int64_t prev_level_rd, int64_t ref_best_rd,
    int *is_cost_valid, FAST_TX_SEARCH_MODE ftxs_mode,
    TXB_RD_INFO_NODE *rd_info_node) {
  assert(tx_size < TX_SIZES_ALL);
  av1_init_rd_stats(rd_stats);
  if (ref_best_rd < 0) {
    *is_cost_valid = 0;
    return;
  }

  MACROBLOCKD *const xd = &x->e_mbd;
  const int max_blocks_high = max_block_high(xd, plane_bsize, 0);
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, 0);
  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

  const int bw = mi_size_wide[plane_bsize];
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int ctx = txfm_partition_context(tx_above + blk_col, tx_left + blk_row,
                                         mbmi->sb_type, tx_size);
  struct macroblock_plane *const p = &x->plane[0];

  const int try_no_split =
      cpi->oxcf.enable_tx64 || txsize_sqr_up_map[tx_size] != TX_64X64;
  int try_split = tx_size > TX_4X4 && depth < MAX_VARTX_DEPTH;
#if CONFIG_DIST_8X8
  if (x->using_dist_8x8)
    try_split &= tx_size_wide[tx_size] >= 16 && tx_size_high[tx_size] >= 16;
#endif
  TxCandidateInfo no_split = { INT64_MAX, 0, TX_TYPES };

  // TX no split
  if (try_no_split) {
    try_tx_block_no_split(cpi, x, blk_row, blk_col, block, tx_size, depth,
                          plane_bsize, ta, tl, ctx, rd_stats, ref_best_rd,
                          ftxs_mode, rd_info_node, &no_split);

    if (cpi->sf.tx_sf.adaptive_txb_search_level &&
        (no_split.rd -
         (no_split.rd >> (1 + cpi->sf.tx_sf.adaptive_txb_search_level))) >
            ref_best_rd) {
      *is_cost_valid = 0;
      return;
    }

    if (cpi->sf.tx_sf.txb_split_cap) {
      if (p->eobs[block] == 0) try_split = 0;
    }

    if (cpi->sf.tx_sf.adaptive_txb_search_level &&
        (no_split.rd -
         (no_split.rd >> (2 + cpi->sf.tx_sf.adaptive_txb_search_level))) >
            prev_level_rd) {
      try_split = 0;
    }
  }

  if (x->e_mbd.bd == 8 && try_split &&
      !(ref_best_rd == INT64_MAX && no_split.rd == INT64_MAX)) {
    const int threshold = cpi->sf.tx_sf.tx_type_search.ml_tx_split_thresh;
    if (threshold >= 0) {
      const int split_score =
          ml_predict_tx_split(x, plane_bsize, blk_row, blk_col, tx_size);
      if (split_score < -threshold) try_split = 0;
    }
  }

  // TX split
  int64_t split_rd = INT64_MAX;
  RD_STATS split_rd_stats;
  av1_init_rd_stats(&split_rd_stats);
  if (try_split) {
    try_tx_block_split(cpi, x, blk_row, blk_col, block, tx_size, depth,
                       plane_bsize, ta, tl, tx_above, tx_left, ctx, no_split.rd,
                       AOMMIN(no_split.rd, ref_best_rd), ftxs_mode,
                       rd_info_node, &split_rd_stats, &split_rd);
  }

  if (no_split.rd < split_rd) {
    ENTROPY_CONTEXT *pta = ta + blk_col;
    ENTROPY_CONTEXT *ptl = tl + blk_row;
    const TX_SIZE tx_size_selected = tx_size;
    p->txb_entropy_ctx[block] = no_split.txb_entropy_ctx;
    av1_set_txb_context(x, 0, block, tx_size_selected, pta, ptl);
    txfm_partition_update(tx_above + blk_col, tx_left + blk_row, tx_size,
                          tx_size);
    for (int idy = 0; idy < tx_size_high_unit[tx_size]; ++idy) {
      for (int idx = 0; idx < tx_size_wide_unit[tx_size]; ++idx) {
        const int index =
            av1_get_txb_size_index(plane_bsize, blk_row + idy, blk_col + idx);
        mbmi->inter_tx_size[index] = tx_size_selected;
      }
    }
    mbmi->tx_size = tx_size_selected;
    update_txk_array(xd, blk_row, blk_col, tx_size, no_split.tx_type);
    set_blk_skip(x, 0, blk_row * bw + blk_col, rd_stats->skip);
  } else {
    *rd_stats = split_rd_stats;
    if (split_rd == INT64_MAX) *is_cost_valid = 0;
  }
}

void super_block_yrd(const AV1_COMP *const cpi, MACROBLOCK *x,
                     RD_STATS *rd_stats, BLOCK_SIZE bs, int64_t ref_best_rd) {
  MACROBLOCKD *xd = &x->e_mbd;
  av1_init_rd_stats(rd_stats);
  int is_inter = is_inter_block(xd->mi[0]);
  assert(bs == xd->mi[0]->sb_type);

  const int mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
  const int mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);

  uint32_t hash = 0;
  int32_t match_index = -1;
  MB_RD_RECORD *mb_rd_record = NULL;
  const int within_border = mi_row >= xd->tile.mi_row_start &&
                            (mi_row + mi_size_high[bs] < xd->tile.mi_row_end) &&
                            mi_col >= xd->tile.mi_col_start &&
                            (mi_col + mi_size_wide[bs] < xd->tile.mi_col_end);
  const int is_mb_rd_hash_enabled =
      (within_border && cpi->sf.rd_sf.use_mb_rd_hash && is_inter);
  const int n4 = bsize_to_num_blk(bs);
  if (is_mb_rd_hash_enabled) {
    hash = get_block_residue_hash(x, bs);
    mb_rd_record = &x->mb_rd_record;
    match_index = find_mb_rd_info(mb_rd_record, ref_best_rd, hash);
    if (match_index != -1) {
      MB_RD_INFO *tx_rd_info = &mb_rd_record->tx_rd_info[match_index];
      fetch_tx_rd_info(n4, tx_rd_info, rd_stats, x);
      return;
    }
  }

  // If we predict that skip is the optimal RD decision - set the respective
  // context and terminate early.
  int64_t dist;

  if (x->predict_skip_level && is_inter &&
      (!xd->lossless[xd->mi[0]->segment_id]) &&
      predict_skip_flag(x, bs, &dist, cpi->common.reduced_tx_set_used)) {
    // Populate rdstats as per skip decision
    set_skip_flag(x, rd_stats, bs, dist);
    // Save the RD search results into tx_rd_record.
    if (is_mb_rd_hash_enabled)
      save_tx_rd_info(n4, hash, x, rd_stats, mb_rd_record);
    return;
  }

  if (xd->lossless[xd->mi[0]->segment_id]) {
    choose_smallest_tx_size(cpi, x, rd_stats, ref_best_rd, bs);
  } else if (x->tx_size_search_method == USE_LARGESTALL) {
    choose_largest_tx_size(cpi, x, rd_stats, ref_best_rd, bs);
  } else {
    choose_tx_size_type_from_rd(cpi, x, rd_stats, ref_best_rd, bs);
  }

  // Save the RD search results into tx_rd_record.
  if (is_mb_rd_hash_enabled) {
    assert(mb_rd_record != NULL);
    save_tx_rd_info(n4, hash, x, rd_stats, mb_rd_record);
  }
}

// Return value 0: early termination triggered, no valid rd cost available;
//              1: rd cost values are valid.
int super_block_uvrd(const AV1_COMP *const cpi, MACROBLOCK *x,
                     RD_STATS *rd_stats, BLOCK_SIZE bsize,
                     int64_t ref_best_rd) {
  av1_init_rd_stats(rd_stats);
  int is_cost_valid = 1;
  if (ref_best_rd < 0) is_cost_valid = 0;
  if (x->skip_chroma_rd || !is_cost_valid) return is_cost_valid;

  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_U];
  int plane;
  const int is_inter = is_inter_block(mbmi);
  int64_t this_rd = 0, skip_rd = 0;
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);

  if (is_inter && is_cost_valid) {
    for (plane = 1; plane < MAX_MB_PLANE; ++plane)
      av1_subtract_plane(x, plane_bsize, plane);
  }

  if (is_cost_valid) {
    const TX_SIZE uv_tx_size = av1_get_tx_size(AOM_PLANE_U, xd);
    for (plane = 1; plane < MAX_MB_PLANE; ++plane) {
      RD_STATS pn_rd_stats;
      int64_t chroma_ref_best_rd = ref_best_rd;
      // For inter blocks, refined ref_best_rd is used for early exit
      // For intra blocks, even though current rd crosses ref_best_rd, early
      // exit is not recommended as current rd is used for gating subsequent
      // modes as well (say, for angular modes)
      // TODO(any): Extend the early exit mechanism for intra modes as well
      if (cpi->sf.inter_sf.perform_best_rd_based_gating_for_chroma &&
          is_inter && chroma_ref_best_rd != INT64_MAX)
        chroma_ref_best_rd = ref_best_rd - AOMMIN(this_rd, skip_rd);
      txfm_rd_in_plane(x, cpi, &pn_rd_stats, chroma_ref_best_rd, 0, plane,
                       plane_bsize, uv_tx_size,
                       cpi->sf.rd_sf.use_fast_coef_costing, FTXS_NONE, 0);
      if (pn_rd_stats.rate == INT_MAX) {
        is_cost_valid = 0;
        break;
      }
      av1_merge_rd_stats(rd_stats, &pn_rd_stats);
      this_rd = RDCOST(x->rdmult, rd_stats->rate, rd_stats->dist);
      skip_rd = RDCOST(x->rdmult, 0, rd_stats->sse);
      if (AOMMIN(this_rd, skip_rd) > ref_best_rd) {
        is_cost_valid = 0;
        break;
      }
    }
  }

  if (!is_cost_valid) {
    // reset cost value
    av1_invalid_rd_stats(rd_stats);
  }

  return is_cost_valid;
}

int txfm_search(const AV1_COMP *cpi, const TileDataEnc *tile_data,
                MACROBLOCK *x, BLOCK_SIZE bsize, RD_STATS *rd_stats,
                RD_STATS *rd_stats_y, RD_STATS *rd_stats_uv, int mode_rate,
                int64_t ref_best_rd) {
  /*
   * This function combines y and uv planes' transform search processes
   * together, when the prediction is generated. It first does subtraction to
   * obtain the prediction error. Then it calls
   * pick_tx_size_type_yrd/super_block_yrd and super_block_uvrd sequentially and
   * handles the early terminations happening in those functions. At the end, it
   * computes the rd_stats/_y/_uv accordingly.
   */
  const AV1_COMMON *cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int ref_frame_1 = mbmi->ref_frame[1];
  const int64_t mode_rd = RDCOST(x->rdmult, mode_rate, 0);
  const int64_t rd_thresh =
      ref_best_rd == INT64_MAX ? INT64_MAX : ref_best_rd - mode_rd;
  const int skip_ctx = av1_get_skip_context(xd);
  const int skip_flag_cost[2] = { x->skip_cost[skip_ctx][0],
                                  x->skip_cost[skip_ctx][1] };
  const int64_t min_header_rate =
      mode_rate + AOMMIN(skip_flag_cost[0], skip_flag_cost[1]);
  // Account for minimum skip and non_skip rd.
  // Eventually either one of them will be added to mode_rate
  const int64_t min_header_rd_possible = RDCOST(x->rdmult, min_header_rate, 0);
  (void)tile_data;

  if (min_header_rd_possible > ref_best_rd) {
    av1_invalid_rd_stats(rd_stats_y);
    return 0;
  }

  av1_init_rd_stats(rd_stats);
  av1_init_rd_stats(rd_stats_y);
  rd_stats->rate = mode_rate;

  // cost and distortion
  av1_subtract_plane(x, bsize, 0);
  if (x->tx_mode_search_type == TX_MODE_SELECT &&
      !xd->lossless[mbmi->segment_id]) {
    pick_tx_size_type_yrd(cpi, x, rd_stats_y, bsize, rd_thresh);
#if CONFIG_COLLECT_RD_STATS == 2
    PrintPredictionUnitStats(cpi, tile_data, x, rd_stats_y, bsize);
#endif  // CONFIG_COLLECT_RD_STATS == 2
  } else {
    super_block_yrd(cpi, x, rd_stats_y, bsize, rd_thresh);
    memset(mbmi->inter_tx_size, mbmi->tx_size, sizeof(mbmi->inter_tx_size));
    for (int i = 0; i < xd->n4_h * xd->n4_w; ++i)
      set_blk_skip(x, 0, i, rd_stats_y->skip);
  }

  if (rd_stats_y->rate == INT_MAX) {
    // TODO(angiebird): check if we need this
    // restore_dst_buf(xd, *orig_dst, num_planes);
    mbmi->ref_frame[1] = ref_frame_1;
    return 0;
  }

  av1_merge_rd_stats(rd_stats, rd_stats_y);

  const int64_t non_skip_rdcosty =
      RDCOST(x->rdmult, rd_stats->rate + skip_flag_cost[0], rd_stats->dist);
  const int64_t skip_rdcosty =
      RDCOST(x->rdmult, mode_rate + skip_flag_cost[1], rd_stats->sse);
  const int64_t min_rdcosty = AOMMIN(non_skip_rdcosty, skip_rdcosty);
  if (min_rdcosty > ref_best_rd) {
    const int64_t tokenonly_rdy =
        AOMMIN(RDCOST(x->rdmult, rd_stats_y->rate, rd_stats_y->dist),
               RDCOST(x->rdmult, 0, rd_stats_y->sse));
    // Invalidate rd_stats_y to skip the rest of the motion modes search
    if (tokenonly_rdy -
            (tokenonly_rdy >> cpi->sf.inter_sf.prune_motion_mode_level) >
        rd_thresh)
      av1_invalid_rd_stats(rd_stats_y);
    mbmi->ref_frame[1] = ref_frame_1;
    return 0;
  }

  av1_init_rd_stats(rd_stats_uv);
  const int num_planes = av1_num_planes(cm);
  if (num_planes > 1) {
    int64_t ref_best_chroma_rd = ref_best_rd;
    // Calculate best rd cost possible for chroma
    if (cpi->sf.inter_sf.perform_best_rd_based_gating_for_chroma &&
        (ref_best_chroma_rd != INT64_MAX)) {
      ref_best_chroma_rd =
          (ref_best_chroma_rd - AOMMIN(non_skip_rdcosty, skip_rdcosty));
    }
    const int is_cost_valid_uv =
        super_block_uvrd(cpi, x, rd_stats_uv, bsize, ref_best_chroma_rd);
    if (!is_cost_valid_uv) {
      mbmi->ref_frame[1] = ref_frame_1;
      return 0;
    }
    av1_merge_rd_stats(rd_stats, rd_stats_uv);
  }

  if (rd_stats->skip) {
    rd_stats->rate -= rd_stats_uv->rate + rd_stats_y->rate;
    rd_stats_y->rate = 0;
    rd_stats_uv->rate = 0;
    rd_stats->dist = rd_stats->sse;
    rd_stats_y->dist = rd_stats_y->sse;
    rd_stats_uv->dist = rd_stats_uv->sse;
    rd_stats->rate += skip_flag_cost[1];
    mbmi->skip = 1;
    // here mbmi->skip temporarily plays a role as what this_skip2 does

    const int64_t tmprd = RDCOST(x->rdmult, rd_stats->rate, rd_stats->dist);
    if (tmprd > ref_best_rd) {
      mbmi->ref_frame[1] = ref_frame_1;
      return 0;
    }
  } else if (!xd->lossless[mbmi->segment_id] &&
             (RDCOST(x->rdmult,
                     rd_stats_y->rate + rd_stats_uv->rate + skip_flag_cost[0],
                     rd_stats->dist) >=
              RDCOST(x->rdmult, skip_flag_cost[1], rd_stats->sse))) {
    rd_stats->rate -= rd_stats_uv->rate + rd_stats_y->rate;
    rd_stats->rate += skip_flag_cost[1];
    rd_stats->dist = rd_stats->sse;
    rd_stats_y->dist = rd_stats_y->sse;
    rd_stats_uv->dist = rd_stats_uv->sse;
    rd_stats_y->rate = 0;
    rd_stats_uv->rate = 0;
    mbmi->skip = 1;
  } else {
    rd_stats->rate += skip_flag_cost[0];
    mbmi->skip = 0;
  }

  return 1;
}

static AOM_INLINE void block_rd_txfm(int plane, int block, int blk_row,
                                     int blk_col, BLOCK_SIZE plane_bsize,
                                     TX_SIZE tx_size, void *arg) {
  struct rdcost_block_args *args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int is_inter = is_inter_block(xd->mi[0]);
  const AV1_COMP *cpi = args->cpi;
  ENTROPY_CONTEXT *a = args->t_above + blk_col;
  ENTROPY_CONTEXT *l = args->t_left + blk_row;
  const AV1_COMMON *cm = &cpi->common;
  RD_STATS this_rd_stats;

  av1_init_rd_stats(&this_rd_stats);

  if (args->exit_early) {
    args->incomplete_exit = 1;
    return;
  }

  if (!is_inter) {
    av1_predict_intra_block_facade(cm, xd, plane, blk_col, blk_row, tx_size);
    av1_subtract_txb(x, plane, plane_bsize, blk_col, blk_row, tx_size);
  }
  TXB_CTX txb_ctx;
  get_txb_ctx(plane_bsize, tx_size, plane, a, l, &txb_ctx);
  search_txk_type(cpi, x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                  &txb_ctx, args->ftxs_mode, args->use_fast_coef_costing,
                  args->skip_trellis, args->best_rd - args->this_rd,
                  &this_rd_stats);

  if (plane == AOM_PLANE_Y && xd->cfl.store_y) {
    assert(!is_inter || plane_bsize < BLOCK_8X8);
    cfl_store_tx(xd, blk_row, blk_col, tx_size, plane_bsize);
  }

#if CONFIG_RD_DEBUG
  av1_update_txb_coeff_cost(&this_rd_stats, plane, tx_size, blk_row, blk_col,
                            this_rd_stats.rate);
#endif  // CONFIG_RD_DEBUG
  av1_set_txb_context(x, plane, block, tx_size, a, l);

  const int blk_idx =
      blk_row * (block_size_wide[plane_bsize] >> MI_SIZE_LOG2) + blk_col;

  if (plane == 0)
    set_blk_skip(x, plane, blk_idx, x->plane[plane].eobs[block] == 0);
  else
    set_blk_skip(x, plane, blk_idx, 0);

  int64_t rd;
  if (is_inter) {
    const int64_t rd1 =
        RDCOST(x->rdmult, this_rd_stats.rate, this_rd_stats.dist);
    const int64_t rd2 = RDCOST(x->rdmult, 0, this_rd_stats.sse);

    // TODO(jingning): temporarily enabled only for luma component
    rd = AOMMIN(rd1, rd2);
    this_rd_stats.skip &= !x->plane[plane].eobs[block];
  } else {
    // Signal non-skip for Intra blocks
    rd = RDCOST(x->rdmult, this_rd_stats.rate, this_rd_stats.dist);
    this_rd_stats.skip = 0;
  }

  av1_merge_rd_stats(&args->rd_stats, &this_rd_stats);

  args->this_rd += rd;

  if (args->this_rd > args->best_rd) args->exit_early = 1;
}

void txfm_rd_in_plane(MACROBLOCK *x, const AV1_COMP *cpi, RD_STATS *rd_stats,
                      int64_t ref_best_rd, int64_t this_rd, int plane,
                      BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                      int use_fast_coef_casting, FAST_TX_SEARCH_MODE ftxs_mode,
                      int skip_trellis) {
  if (!cpi->oxcf.enable_tx64 && txsize_sqr_up_map[tx_size] == TX_64X64) {
    av1_invalid_rd_stats(rd_stats);
    return;
  }

  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  struct rdcost_block_args args;
  av1_zero(args);
  args.x = x;
  args.cpi = cpi;
  args.best_rd = ref_best_rd;
  args.use_fast_coef_costing = use_fast_coef_casting;
  args.ftxs_mode = ftxs_mode;
  args.this_rd = this_rd;
  args.skip_trellis = skip_trellis;
  av1_init_rd_stats(&args.rd_stats);

  if (plane == 0) xd->mi[0]->tx_size = tx_size;

  av1_get_entropy_contexts(plane_bsize, pd, args.t_above, args.t_left);

  if (args.this_rd > args.best_rd) {
    args.exit_early = 1;
  }

  av1_foreach_transformed_block_in_plane(xd, plane_bsize, plane, block_rd_txfm,
                                         &args);

  MB_MODE_INFO *const mbmi = xd->mi[0];
  const int is_inter = is_inter_block(mbmi);
  const int invalid_rd = is_inter ? args.incomplete_exit : args.exit_early;

  if (invalid_rd) {
    av1_invalid_rd_stats(rd_stats);
  } else {
    *rd_stats = args.rd_stats;
  }
}
