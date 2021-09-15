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

#ifndef AOM_AV1_ENCODER_PARTITION_SEARCH_H_
#define AOM_AV1_ENCODER_PARTITION_SEARCH_H_

#include "av1/encoder/block.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodeframe.h"
#include "av1/encoder/encodeframe_utils.h"
#include "av1/encoder/tokenize.h"

void av1_set_offsets_without_segment_id(const AV1_COMP *const cpi,
                                        const TileInfo *const tile,
                                        MACROBLOCK *const x, int mi_row,
                                        int mi_col, BLOCK_SIZE bsize,
                                        const CHROMA_REF_INFO *chr_ref_info);
void av1_set_offsets(const AV1_COMP *const cpi, const TileInfo *const tile,
                     MACROBLOCK *const x, int mi_row, int mi_col,
                     BLOCK_SIZE bsize, const CHROMA_REF_INFO *chr_ref_info);
void av1_rd_use_partition(AV1_COMP *cpi, ThreadData *td, TileDataEnc *tile_data,
                          MB_MODE_INFO **mib, TokenExtra **tp, int mi_row,
                          int mi_col, BLOCK_SIZE bsize, int *rate,
                          int64_t *dist, int do_recon, PARTITION_TREE *ptree,
                          PC_TREE *pc_tree);
void av1_nonrd_use_partition(AV1_COMP *cpi, ThreadData *td,
                             TileDataEnc *tile_data, MB_MODE_INFO **mib,
                             TokenExtra **tp, int mi_row, int mi_col,
                             BLOCK_SIZE bsize, PC_TREE *pc_tree,
                             PARTITION_TREE *ptree);
bool av1_rd_pick_partition(AV1_COMP *const cpi, ThreadData *td,
                           TileDataEnc *tile_data, TokenExtra **tp, int mi_row,
                           int mi_col, BLOCK_SIZE bsize, RD_STATS *rd_cost,
                           RD_STATS best_rdc, PC_TREE *pc_tree,
#if CONFIG_SDP && CONFIG_EXT_RECUR_PARTITIONS
                           PARTITION_TREE *ptree_luma,
#endif  // CONFIG_SDP && CONFIG_EXT_RECUR_PARTITIONS
                           SIMPLE_MOTION_DATA_TREE *sms_tree, int64_t *none_rd,
                           SB_MULTI_PASS_MODE multi_pass_mode,
                           RD_RECT_PART_WIN_INFO *rect_part_win_info,
                           bool enable_full_recursion);
#if CONFIG_EXT_RECUR_PARTITIONS
void av1_build_partition_tree_fixed_partitioning(AV1_COMMON *const cm,
                                                 int mi_row, int mi_col,
                                                 BLOCK_SIZE bsize,
                                                 PARTITION_TREE *ptree);
#endif  // CONFIG_EXT_RECUR_PARTITIONS
void setup_block_rdmult(const AV1_COMP *const cpi, MACROBLOCK *const x,
                        int mi_row, int mi_col, BLOCK_SIZE bsize,
                        AQ_MODE aq_mode, MB_MODE_INFO *mbmi);

#endif  // AOM_AV1_ENCODER_PARTITION_SEARCH_H_
