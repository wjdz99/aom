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

#include <limits.h>

#include "aom_mem/aom_mem.h"

#include "av1/common/pred_common.h"
#include "av1/common/tile_common.h"

#include "av1/encoder/cost.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/subexp.h"

void av1_enable_segmentation(struct segmentation *seg) {
  seg->enabled = 1;
  seg->update_map = 1;
  seg->update_data = 1;
}

void av1_disable_segmentation(struct segmentation *seg) {
  seg->enabled = 0;
  seg->update_map = 0;
  seg->update_data = 0;
}

void av1_set_segment_data(struct segmentation *seg, signed char *feature_data,
                          unsigned char abs_delta) {
  seg->abs_delta = abs_delta;

  memcpy(seg->feature_data, feature_data, sizeof(seg->feature_data));
}
void av1_disable_segfeature(struct segmentation *seg, int segment_id,
                            SEG_LVL_FEATURES feature_id) {
  seg->feature_mask[segment_id] &= ~(1 << feature_id);
}

void av1_clear_segdata(struct segmentation *seg, int segment_id,
                       SEG_LVL_FEATURES feature_id) {
  seg->feature_data[segment_id][feature_id] = 0;
}
#if !CONFIG_EXT_SEGMENT
// Based on set of segment counts calculate a probability tree
static void calc_segtree_probs(unsigned *segcounts,
                               aom_prob *segment_tree_probs,
                               const aom_prob *cur_tree_probs,
                               const int probwt) {
  // Work out probabilities of each segment
  const unsigned cc[4] = { segcounts[0] + segcounts[1],
                           segcounts[2] + segcounts[3],
                           segcounts[4] + segcounts[5],
                           segcounts[6] + segcounts[7] };
  const unsigned ccc[2] = { cc[0] + cc[1], cc[2] + cc[3] };
  int i;

  segment_tree_probs[0] = get_binary_prob(ccc[0], ccc[1]);
  segment_tree_probs[1] = get_binary_prob(cc[0], cc[1]);
  segment_tree_probs[2] = get_binary_prob(cc[2], cc[3]);
  segment_tree_probs[3] = get_binary_prob(segcounts[0], segcounts[1]);
  segment_tree_probs[4] = get_binary_prob(segcounts[2], segcounts[3]);
  segment_tree_probs[5] = get_binary_prob(segcounts[4], segcounts[5]);
  segment_tree_probs[6] = get_binary_prob(segcounts[6], segcounts[7]);

  for (i = 0; i < 7; i++) {
    const unsigned *ct =
        i == 0 ? ccc : i < 3 ? cc + (i & 2) : segcounts + (i - 3) * 2;
    av1_prob_diff_update_savings_search(ct, cur_tree_probs[i],
                                        &segment_tree_probs[i],
                                        DIFF_UPDATE_PROB, probwt);
  }
}

// Based on set of segment counts and probabilities calculate a cost estimate
static int cost_segmap(unsigned *segcounts, aom_prob *probs) {
  const int c01 = segcounts[0] + segcounts[1];
  const int c23 = segcounts[2] + segcounts[3];
  const int c45 = segcounts[4] + segcounts[5];
  const int c67 = segcounts[6] + segcounts[7];
  const int c0123 = c01 + c23;
  const int c4567 = c45 + c67;

  // Cost the top node of the tree
  int cost = c0123 * av1_cost_zero(probs[0]) + c4567 * av1_cost_one(probs[0]);

  // Cost subsequent levels
  if (c0123 > 0) {
    cost += c01 * av1_cost_zero(probs[1]) + c23 * av1_cost_one(probs[1]);

    if (c01 > 0)
      cost += segcounts[0] * av1_cost_zero(probs[3]) +
              segcounts[1] * av1_cost_one(probs[3]);
    if (c23 > 0)
      cost += segcounts[2] * av1_cost_zero(probs[4]) +
              segcounts[3] * av1_cost_one(probs[4]);
  }

  if (c4567 > 0) {
    cost += c45 * av1_cost_zero(probs[2]) + c67 * av1_cost_one(probs[2]);

    if (c45 > 0)
      cost += segcounts[4] * av1_cost_zero(probs[5]) +
              segcounts[5] * av1_cost_one(probs[5]);
    if (c67 > 0)
      cost += segcounts[6] * av1_cost_zero(probs[6]) +
              segcounts[7] * av1_cost_one(probs[6]);
  }

  return cost;
}
#endif
#if CONFIG_EXT_SEGMENT
static void count_segs(
    const AV1_COMMON *cm, MACROBLOCKD *xd, const TileInfo *tile, MODE_INFO **mi,
    unsigned (*harmonic_pred_seg_counts)[HARMONIC_SPATIAL_PRED_PROBS + 1],
    unsigned (
        *heterogeneous_pred_seg_counts)[HETEROGENEOUS_SPATIAL_PRED_PROBS + 1],
    unsigned (*temporal_predictor_count)[2],
    unsigned (*t_unpred_seg_counts)[MAX_SEGMENTS], SEG_CATEGORIES seg_cat_idx,
    int bw, int bh, int mi_row, int mi_col) {
  int segment_id;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  xd->mi = mi;
  segment_id = xd->mi[0]->mbmi.segment_id[seg_cat_idx];

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  {
    int i, j;
    int seg_map_minb_size =
        1 << (cm->seg[seg_cat_idx].seg_map_minb_size_log2_minus3);
    const int seg_id_written = !(mi_row & (seg_map_minb_size - 1)) &&
                               !(mi_col & (seg_map_minb_size - 1));
    const int last_block_seg_minb =
        (mi_row + mi_size_high[xd->mi[0]->mbmi.sb_type] == cm->mi_rows ||
         !((mi_row + mi_size_high[xd->mi[0]->mbmi.sb_type]) &
           (seg_map_minb_size - 1))) &&
        (mi_col + mi_size_wide[xd->mi[0]->mbmi.sb_type] == cm->mi_cols ||
         !((mi_col + mi_size_wide[xd->mi[0]->mbmi.sb_type]) &
           (seg_map_minb_size - 1)));
    const int seg_id_skip =
        ((xd->mi[0]->mbmi.ref_frame > INTRA_FRAME) ||
         !av1_is_segfeature_enabled(&cm->seg[QUALITY_SEG_IDX],
                                    QUALITY_SEG_LVL_ALT_LF)) &&
                xd->mi[0]->mbmi.skip &&
                !cm->seg[seg_cat_idx].skip_seg_id_disabled
            ? 1
            : 0;
    MODE_INFO *top_mi =
        xd->mi[-((mi_row & (seg_map_minb_size - 1)) + 1) * cm->mi_stride -
               (mi_col & (seg_map_minb_size - 1))];
    MODE_INFO *left_mi =
        xd->mi[-(mi_row & (seg_map_minb_size - 1)) * cm->mi_stride -
               ((mi_col & (seg_map_minb_size - 1)) + 1)];
    const int pred_segment_id = get_segment_id(
        cm, cm->last_frame_seg_map,
        seg_map_minb_size == 1 ? BLOCK_8X8 : seg_map_minb_size == 2
                                                 ? BLOCK_16X16
                                                 : seg_map_minb_size == 4
                                                       ? BLOCK_32X32
                                                       : BLOCK_64X64,
        seg_cat_idx, (mi_row & ~(seg_map_minb_size - 1)),
        (mi_col & ~(seg_map_minb_size - 1)));
    if (seg_id_written) xd->seg_id_is_done[seg_cat_idx] = 0;
    if (xd->seg_id_is_done[seg_cat_idx] ||
        (!xd->seg_id_is_done[seg_cat_idx] && seg_id_skip &&
         !last_block_seg_minb))
      return;
    if (!xd->seg_id_is_done[seg_cat_idx] && seg_id_skip &&
        last_block_seg_minb) {
      for (i = 0; i < AOMMIN(AOMMAX(mi_size_high[xd->mi[0]->mbmi.sb_type],
                                    seg_map_minb_size),
                             cm->mi_rows - (mi_row & ~(seg_map_minb_size - 1)));
           i++) {
        for (j = 0;
             j < AOMMIN(AOMMAX(mi_size_wide[xd->mi[0]->mbmi.sb_type],
                               seg_map_minb_size),
                        cm->mi_cols - (mi_col & ~(seg_map_minb_size - 1)));
             j++) {
          MODE_INFO *l_mi =
              xd->mi[-(mi_row & (seg_map_minb_size - 1)) * cm->mi_stride -
                     ((mi_col & (seg_map_minb_size - 1))) + i * cm->mi_stride +
                     j];
          if (cm->seg[seg_cat_idx].temporal_update &&
              cm->seg[seg_cat_idx].skip_seg_id_replacement ==
                  SKIP_SEG_ID_TEMPORAL_REPLACE) {
            l_mi->mbmi.segment_id[seg_cat_idx] = pred_segment_id;
            l_mi->mbmi.seg_id_predicted[seg_cat_idx] = 1;
            l_mi->mbmi.seg_id_spatial_predicted[seg_cat_idx] = 0;
          } else {
            l_mi->mbmi.segment_id[seg_cat_idx] =
                left_mi ? left_mi->mbmi.segment_id[seg_cat_idx]
                        : top_mi ? top_mi->mbmi.segment_id[seg_cat_idx] : 0;
            l_mi->mbmi.seg_id_spatial_predicted[seg_cat_idx] =
                left_mi ? SEG_ID_COPY_LEFT : top_mi
                                                 ? SEG_ID_COPY_ABOVE
                                                 : SEG_ID_SPATIAL_UNPREDICTED;
            l_mi->mbmi.seg_id_predicted[seg_cat_idx] = 0;
          }
        }
      }
      return;
    }
    if (!seg_id_skip) {
      int pred_flag = 0;
      int spatial_pred_idc = SEG_ID_SPATIAL_UNPREDICTED;
      xd->seg_id_is_done[seg_cat_idx] = 1;
      if (cm->seg[seg_cat_idx].temporal_update) {
        const int pred_context =
            (top_mi ? top_mi->mbmi.seg_id_predicted[seg_cat_idx] : 0) +
            (left_mi ? left_mi->mbmi.seg_id_predicted[seg_cat_idx] : 0);
        pred_flag = pred_segment_id == segment_id;
        xd->mi[0]->mbmi.seg_id_predicted[seg_cat_idx] = pred_flag;
        temporal_predictor_count[pred_context][pred_flag]++;
      }
      if (!pred_flag) {
        const int spatial_seg_pred_ctx =
            top_mi && left_mi &&
                    top_mi->mbmi.segment_id[seg_cat_idx] ==
                        left_mi->mbmi.segment_id[seg_cat_idx]
                ? 1
                : 0;
        if (spatial_seg_pred_ctx) {
          spatial_pred_idc = segment_id == left_mi->mbmi.segment_id[seg_cat_idx]
                                 ? SEG_ID_COPY_LEFT
                                 : SEG_ID_SPATIAL_UNPREDICTED;
          harmonic_pred_seg_counts[0][spatial_pred_idc]++;
        } else {
          if (left_mi)
            spatial_pred_idc =
                segment_id == left_mi->mbmi.segment_id[seg_cat_idx]
                    ? SEG_ID_COPY_LEFT
                    : SEG_ID_SPATIAL_UNPREDICTED;
          if (top_mi && spatial_pred_idc == SEG_ID_SPATIAL_UNPREDICTED)
            spatial_pred_idc =
                segment_id == top_mi->mbmi.segment_id[seg_cat_idx]
                    ? SEG_ID_COPY_ABOVE
                    : SEG_ID_SPATIAL_UNPREDICTED;
          heterogeneous_pred_seg_counts[0][spatial_pred_idc]++;
        }
        if (spatial_pred_idc == SEG_ID_SPATIAL_UNPREDICTED) {
          t_unpred_seg_counts[0][segment_id]++;
        }
      }
      for (i = 0; i < AOMMIN(AOMMAX(mi_size_high[xd->mi[0]->mbmi.sb_type],
                                    seg_map_minb_size),
                             cm->mi_rows - (mi_row & ~(seg_map_minb_size - 1)));
           i++) {
        for (j = 0;
             j < AOMMIN(AOMMAX(mi_size_wide[xd->mi[0]->mbmi.sb_type],
                               seg_map_minb_size),
                        cm->mi_cols - (mi_col & ~(seg_map_minb_size - 1)));
             j++) {
          MODE_INFO *l_mi =
              xd->mi[-(mi_row & (seg_map_minb_size - 1)) * cm->mi_stride -
                     ((mi_col & (seg_map_minb_size - 1))) + i * cm->mi_stride +
                     j];
          l_mi->mbmi.segment_id[seg_cat_idx] = segment_id;
          l_mi->mbmi.seg_id_predicted[seg_cat_idx] = pred_flag;
          l_mi->mbmi.seg_id_spatial_predicted[seg_cat_idx] = spatial_pred_idc;
        }
      }
      return;
    }
  }
}
#else
static void count_segs(const AV1_COMMON *cm, MACROBLOCKD *xd,
                       const TileInfo *tile, MODE_INFO **mi,
                       unsigned *no_pred_segcounts,
                       unsigned (*temporal_predictor_count)[2],
                       unsigned *t_unpred_seg_counts, int bw, int bh,
                       int mi_row, int mi_col) {
  int segment_id;

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

  xd->mi = mi;
  segment_id = xd->mi[0]->mbmi.segment_id;

  set_mi_row_col(xd, tile, mi_row, bh, mi_col, bw, cm->mi_rows, cm->mi_cols);

  // Count the number of hits on each segment with no prediction
  no_pred_segcounts[segment_id]++;

  // Temporal prediction not allowed on key frames
  if (cm->frame_type != KEY_FRAME) {
    const BLOCK_SIZE bsize = xd->mi[0]->mbmi.sb_type;
    // Test to see if the segment id matches the predicted value.

    const int pred_segment_id =
        get_segment_id(cm, cm->last_frame_seg_map, bsize, mi_row, mi_col);
    const int pred_flag = pred_segment_id == segment_id;
    const int pred_context = av1_get_pred_context_seg_id(xd);

    // Store the prediction status for this mb and update counts
    // as appropriate
    xd->mi[0]->mbmi.seg_id_predicted = pred_flag;
    temporal_predictor_count[pred_context][pred_flag]++;

    // Update the "unpredicted" segment count
    if (!pred_flag) t_unpred_seg_counts[segment_id]++;
  }
}
#endif
static void count_segs_sb(
    const AV1_COMMON *cm, MACROBLOCKD *xd, const TileInfo *tile, MODE_INFO **mi,
#if CONFIG_EXT_SEGMENT
    unsigned (*harmonic_pred_seg_counts)[HARMONIC_SPATIAL_PRED_PROBS + 1],
    unsigned (
        *heterogeneous_pred_seg_counts)[HETEROGENEOUS_SPATIAL_PRED_PROBS + 1],
    unsigned (*temporal_predictor_count)[2],
    unsigned (*t_unpred_seg_counts)[MAX_SEGMENTS], SEG_CATEGORIES seg_cat_idx,
#else
    unsigned *no_pred_segcounts, unsigned (*temporal_predictor_count)[2],
    unsigned *t_unpred_seg_counts,
#endif
    int mi_row, int mi_col, BLOCK_SIZE bsize) {
  const int mis = cm->mi_stride;
  const int bs = mi_size_wide[bsize], hbs = bs / 2;
#if CONFIG_EXT_PARTITION_TYPES
  PARTITION_TYPE partition;
#else
  int bw, bh;
#endif  // CONFIG_EXT_PARTITION_TYPES

  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols) return;

#if CONFIG_EXT_PARTITION_TYPES
  if (bsize == BLOCK_8X8)
    partition = PARTITION_NONE;
  else
    partition = get_partition(cm, mi_row, mi_col, bsize);
  switch (partition) {
    case PARTITION_NONE:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, bs, bs, mi_row, mi_col);
      break;
    case PARTITION_HORZ:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, bs, hbs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, bs, hbs,
                 mi_row + hbs, mi_col);
      break;
    case PARTITION_VERT:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, hbs, bs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, bs, mi_row,
                 mi_col + hbs);
      break;
    case PARTITION_HORZ_A:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, hbs, hbs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row, mi_col + hbs);
      count_segs(cm, xd, tile, mi + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, bs, hbs,
                 mi_row + hbs, mi_col);
      break;
    case PARTITION_HORZ_B:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, bs, hbs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row + hbs, mi_col);
      count_segs(cm, xd, tile, mi + hbs + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row + hbs, mi_col + hbs);
      break;
    case PARTITION_VERT_A:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, hbs, hbs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row + hbs, mi_col);
      count_segs(cm, xd, tile, mi + hbs, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, bs, mi_row,
                 mi_col + hbs);
      break;
    case PARTITION_VERT_B:
      count_segs(cm, xd, tile, mi, no_pred_segcounts, temporal_predictor_count,
                 t_unpred_seg_counts, hbs, bs, mi_row, mi_col);
      count_segs(cm, xd, tile, mi + hbs, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row, mi_col + hbs);
      count_segs(cm, xd, tile, mi + hbs + hbs * mis, no_pred_segcounts,
                 temporal_predictor_count, t_unpred_seg_counts, hbs, hbs,
                 mi_row + hbs, mi_col + hbs);
      break;
    case PARTITION_SPLIT: {
      const BLOCK_SIZE subsize = subsize_lookup[PARTITION_SPLIT][bsize];
      int n;

      assert(num_8x8_blocks_wide_lookup[mi[0]->mbmi.sb_type] < bs &&
             num_8x8_blocks_high_lookup[mi[0]->mbmi.sb_type] < bs);

      for (n = 0; n < 4; n++) {
        const int mi_dc = hbs * (n & 1);
        const int mi_dr = hbs * (n >> 1);

        count_segs_sb(cm, xd, tile, &mi[mi_dr * mis + mi_dc], no_pred_segcounts,
                      temporal_predictor_count, t_unpred_seg_counts,
                      mi_row + mi_dr, mi_col + mi_dc, subsize);
      }
    } break;
    default: assert(0);
  }
#else
  bw = mi_size_wide[mi[0]->mbmi.sb_type];
  bh = mi_size_high[mi[0]->mbmi.sb_type];

  if (bw == bs && bh == bs) {
    count_segs(cm, xd, tile, mi,
#if CONFIG_EXT_SEGMENT
               harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
               temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
               no_pred_segcounts, temporal_predictor_count, t_unpred_seg_counts,
#endif
               bs, bs, mi_row, mi_col);
  } else if (bw == bs && bh < bs) {
    count_segs(cm, xd, tile, mi,
#if CONFIG_EXT_SEGMENT
               harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
               temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
               no_pred_segcounts, temporal_predictor_count, t_unpred_seg_counts,
#endif
               bs, hbs, mi_row, mi_col);
    count_segs(cm, xd, tile, mi + hbs * mis,
#if CONFIG_EXT_SEGMENT
               harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
               temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
               no_pred_segcounts, temporal_predictor_count, t_unpred_seg_counts,
#endif
               bs, hbs, mi_row + hbs, mi_col);
  } else if (bw < bs && bh == bs) {
    count_segs(cm, xd, tile, mi,
#if CONFIG_EXT_SEGMENT
               harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
               temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
               no_pred_segcounts, temporal_predictor_count, t_unpred_seg_counts,
#endif
               hbs, bs, mi_row, mi_col);
    count_segs(cm, xd, tile, mi + hbs,
#if CONFIG_EXT_SEGMENT
               harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
               temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
               no_pred_segcounts, temporal_predictor_count, t_unpred_seg_counts,
#endif
               hbs, bs, mi_row, mi_col + hbs);
  } else {
    const BLOCK_SIZE subsize = subsize_lookup[PARTITION_SPLIT][bsize];
    int n;

    assert(bw < bs && bh < bs);

    for (n = 0; n < 4; n++) {
      const int mi_dc = hbs * (n & 1);
      const int mi_dr = hbs * (n >> 1);

      count_segs_sb(cm, xd, tile, &mi[mi_dr * mis + mi_dc],
#if CONFIG_EXT_SEGMENT
                    harmonic_pred_seg_counts, heterogeneous_pred_seg_counts,
                    temporal_predictor_count, t_unpred_seg_counts, seg_cat_idx,
#else
                    no_pred_segcounts, temporal_predictor_count,
                    t_unpred_seg_counts,
#endif
                    mi_row + mi_dr, mi_col + mi_dc, subsize);
    }
  }
#endif  // CONFIG_EXT_PARTITION_TYPES
}

#if CONFIG_EXT_SEGMENT
void av1_choose_segmap_coding_method(AV1_COMMON *cm, MACROBLOCKD *xd,
                                     SEG_CATEGORIES seg_cat_idx) {
  struct segmentation *seg = &cm->seg[seg_cat_idx];
  struct segmentation_probs *segp = &cm->fc->seg[seg_cat_idx];
#else
void av1_choose_segmap_coding_method(AV1_COMMON *cm, MACROBLOCKD *xd) {
  struct segmentation *seg = &cm->seg;
  struct segmentation_probs *segp = &cm->fc->seg;
  int no_pred_cost = 0;
#endif
  int t_pred_cost = INT_MAX;
  int i, tile_col, tile_row, mi_row, mi_col;
#if CONFIG_TILE_GROUPS
  const int probwt = cm->num_tg;
#else
  const int probwt = 1;
#endif
#if CONFIG_EXT_SEGMENT
  unsigned temporal_predictor_count[PREDICTION_PROBS][2];
  unsigned harmonic_pred_seg_counts[2][HARMONIC_SPATIAL_PRED_PROBS + 1];
  unsigned heterogeneous_pred_seg_counts[2]
                                        [HETEROGENEOUS_SPATIAL_PRED_PROBS + 1];
  unsigned unpred_seg_counts[2][MAX_SEGMENTS];
  aom_prob no_pred_tree[2][SEG_TREE_PROBS];
  aom_prob harmonic_pred_prob[2][HARMONIC_SPATIAL_PRED_PROBS];
  aom_prob heterogeneous_pred_prob[2][HETEROGENEOUS_SPATIAL_PRED_PROBS];
  aom_prob t_nopred_prob[PREDICTION_PROBS];

  av1_zero(temporal_predictor_count);
  av1_zero(harmonic_pred_seg_counts);
  av1_zero(heterogeneous_pred_seg_counts);
  av1_zero(unpred_seg_counts);
#else
  unsigned(*temporal_predictor_count)[2] = cm->counts.seg.pred;
  unsigned *no_pred_segcounts = cm->counts.seg.tree_total;
  unsigned *t_unpred_seg_counts = cm->counts.seg.tree_mispred;

  aom_prob no_pred_tree[SEG_TREE_PROBS];
  aom_prob t_pred_tree[SEG_TREE_PROBS];
  aom_prob t_nopred_prob[PREDICTION_PROBS];
#endif

  (void)xd;

// We are about to recompute all the segment counts, so zero the accumulators.
#if CONFIG_EXT_SEGMENT
  av1_zero(cm->counts.seg[seg_cat_idx]);
  if (!frame_is_intra_only(cm) && !cm->error_resilient_mode)
    cm->seg[seg_cat_idx].temporal_update = 1;
  else
    cm->seg[seg_cat_idx].temporal_update = 0;
#else
  av1_zero(cm->counts.seg);
#endif

  // First of all generate stats regarding how well the last segment map
  // predicts this one
  for (tile_row = 0; tile_row < cm->tile_rows; tile_row++) {
    TileInfo tile_info;
    av1_tile_set_row(&tile_info, cm, tile_row);
    for (tile_col = 0; tile_col < cm->tile_cols; tile_col++) {
      MODE_INFO **mi_ptr;
      av1_tile_set_col(&tile_info, cm, tile_col);
      mi_ptr = cm->mi_grid_visible + tile_info.mi_row_start * cm->mi_stride +
               tile_info.mi_col_start;
      for (mi_row = tile_info.mi_row_start; mi_row < tile_info.mi_row_end;
           mi_row += cm->mib_size, mi_ptr += cm->mib_size * cm->mi_stride) {
        MODE_INFO **mi = mi_ptr;
        for (mi_col = tile_info.mi_col_start; mi_col < tile_info.mi_col_end;
             mi_col += cm->mib_size, mi += cm->mib_size) {
#if CONFIG_EXT_SEGMENT
          count_segs_sb(cm, xd, &tile_info, mi, harmonic_pred_seg_counts,
                        heterogeneous_pred_seg_counts, temporal_predictor_count,
                        unpred_seg_counts, seg_cat_idx, mi_row, mi_col,
                        cm->sb_size);
#else
          count_segs_sb(cm, xd, &tile_info, mi, no_pred_segcounts,
                        temporal_predictor_count, t_unpred_seg_counts, mi_row,
                        mi_col, cm->sb_size);
#endif
        }
      }
    }
  }

#if CONFIG_EXT_SEGMENT
  {
    // Work out probability tree for coding segments without temporal prediction
    // and the cost.
    unsigned int branch_ct[32][2];
    int count0;
    int count1;

    // calc_segtree_probs(unpred_seg_counts[0], no_pred_tree[0],
    // segp->tree_probs, probwt);
    // no_pred_cost += cost_segmap(unpred_seg_counts[0], no_pred_tree[0]);

    // Key frames cannot use temporal prediction

    t_pred_cost = 0;
    // Work out probability tree for coding those segments not
    // predicted using the temporal method and the cost.
    count0 = harmonic_pred_seg_counts[0][0];
    count1 = harmonic_pred_seg_counts[0][1];

    harmonic_pred_prob[0][0] = get_binary_prob(count0, count1);
    av1_prob_diff_update_savings_search(
        harmonic_pred_seg_counts[0], segp->harmonic_spatial_pred_probs[0],
        &harmonic_pred_prob[0][0], DIFF_UPDATE_PROB, probwt);

    // Add in the predictor signaling cost
    t_pred_cost += count0 * av1_cost_zero(harmonic_pred_prob[0][0]) +
                   count1 * av1_cost_one(harmonic_pred_prob[0][0]);

    av1_tree_probs_from_distribution(av1_seg_id_spatial_prediction_tree,
                                     branch_ct,
                                     heterogeneous_pred_seg_counts[0]);

    for (i = 0; i < HETEROGENEOUS_SPATIAL_PRED_PROBS; i++) {
      count0 = branch_ct[i][0];
      count1 = branch_ct[i][1];

      heterogeneous_pred_prob[0][i] = get_binary_prob(count0, count1);
      av1_prob_diff_update_savings_search(
          branch_ct[i], segp->heterogeneous_spatial_pred_probs[i],
          &heterogeneous_pred_prob[0][i], DIFF_UPDATE_PROB, probwt);

      // Add in the predictor signaling cost
      t_pred_cost += count0 * av1_cost_zero(heterogeneous_pred_prob[0][i]) +
                     count1 * av1_cost_one(heterogeneous_pred_prob[0][i]);
    }

    av1_tree_probs_from_distribution(av1_segment_tree[seg->num_seg - 2],
                                     branch_ct, unpred_seg_counts[0]);

    for (i = 0; i < seg->num_seg - 1; i++) {
      count0 = branch_ct[i][0];
      count1 = branch_ct[i][1];

      no_pred_tree[0][i] = get_binary_prob(count0, count1);
      av1_prob_diff_update_savings_search(branch_ct[i], segp->tree_probs[i],
                                          &no_pred_tree[0][i], DIFF_UPDATE_PROB,
                                          probwt);

      // Add in the predictor signaling cost
      t_pred_cost += count0 * av1_cost_zero(no_pred_tree[0][i]) +
                     count1 * av1_cost_one(no_pred_tree[0][i]);
    }

    // Add in the cost of the signaling for each prediction context.
    if (cm->seg[seg_cat_idx].temporal_update) {
      for (i = 0; i < PREDICTION_PROBS; i++) {
        count0 = temporal_predictor_count[i][0];
        count1 = temporal_predictor_count[i][1];

        t_nopred_prob[i] = get_binary_prob(count0, count1);
        av1_prob_diff_update_savings_search(
            temporal_predictor_count[i], segp->pred_probs[i], &t_nopred_prob[i],
            DIFF_UPDATE_PROB, probwt);

        // Add in the predictor signaling cost
        t_pred_cost += count0 * av1_cost_zero(t_nopred_prob[i]) +
                       count1 * av1_cost_one(t_nopred_prob[i]);
      }
    }
    if (cm->seg[seg_cat_idx].temporal_update)
      av1_copy(cm->counts.seg[seg_cat_idx].temp_pred, temporal_predictor_count);
    av1_copy(cm->counts.seg[seg_cat_idx].tree_mispred, unpred_seg_counts[0]);
    av1_copy(cm->counts.seg[seg_cat_idx].harmonic_spatial_pred,
             harmonic_pred_seg_counts[0]);
    av1_copy(cm->counts.seg[seg_cat_idx].heterogeneous_spatial_pred,
             heterogeneous_pred_seg_counts[0]);
  }
#else
  // Work out probability tree for coding segments without prediction
  // and the cost.
  calc_segtree_probs(no_pred_segcounts, no_pred_tree, segp->tree_probs, probwt);
  no_pred_cost = cost_segmap(no_pred_segcounts, no_pred_tree);

  // Key frames cannot use temporal prediction
  if (!frame_is_intra_only(cm) && !cm->error_resilient_mode) {
    // Work out probability tree for coding those segments not
    // predicted using the temporal method and the cost.
    calc_segtree_probs(t_unpred_seg_counts, t_pred_tree, segp->tree_probs,
                       probwt);
    t_pred_cost = cost_segmap(t_unpred_seg_counts, t_pred_tree);

    // Add in the cost of the signaling for each prediction context.
    for (i = 0; i < PREDICTION_PROBS; i++) {
      const int count0 = temporal_predictor_count[i][0];
      const int count1 = temporal_predictor_count[i][1];

      t_nopred_prob[i] = get_binary_prob(count0, count1);
      av1_prob_diff_update_savings_search(
          temporal_predictor_count[i], segp->pred_probs[i], &t_nopred_prob[i],
          DIFF_UPDATE_PROB, probwt);

      // Add in the predictor signaling cost
      t_pred_cost += count0 * av1_cost_zero(t_nopred_prob[i]) +
                     count1 * av1_cost_one(t_nopred_prob[i]);
    }
  }

  // Now choose which coding method to use.
  if (t_pred_cost < no_pred_cost) {
    assert(!cm->error_resilient_mode);
    seg->temporal_update = 1;
  } else {
    seg->temporal_update = 0;
  }
#endif
}

#if CONFIG_EXT_SEGMENT
void av1_reset_segment_features(AV1_COMMON *cm, SEG_CATEGORIES seg_cat_idx) {
  struct segmentation *seg = &cm->seg[seg_cat_idx];

  // Set up default state for MB feature flags
  seg->enabled = 0;
  seg->update_map = 0;
  seg->update_data = 0;
  av1_clearall_segfeatures(seg);
}
#else
void av1_reset_segment_features(AV1_COMMON *cm) {
  struct segmentation *seg = &cm->seg;

  // Set up default state for MB feature flags
  seg->enabled = 0;
  seg->update_map = 0;
  seg->update_data = 0;
  av1_clearall_segfeatures(seg);
}
#endif