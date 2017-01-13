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

#include <assert.h>

#include "av1/common/blockd.h"
#include "av1/common/loopfilter.h"
#include "av1/common/seg_common.h"
#include "av1/common/quant_common.h"

#if CONFIG_EXT_SEGMENT
static const int seg_feature_data_signed[NUM_SEG_CATEGORIES]
                                        [SEG_LVL_MAX] = { {0, 0}, {1, 1} };

static const int seg_feature_data_max[NUM_SEG_CATEGORIES][SEG_LVL_MAX] = {
  { 3, 0 }, {MAXQ, MAX_LOOP_FILTER}
};

// These functions provide access to new segment level features.
// Eventually these function may be "optimized out" but for the moment,
// the coding mechanism is still subject to change so these provide a
// convenient single point of change.

void av1_clearall_segfeatures(struct segmentation *seg) {
  av1_zero(seg->feature_data);
  av1_zero(seg->feature_mask);
}

void av1_enable_segfeature(struct segmentation *seg, int segment_id,
                           SEG_LVL_FEATURES feature_id) {
  seg->feature_mask[segment_id] |= 1 << feature_id;
}

int av1_seg_feature_data_max(SEG_LVL_FEATURES feature_id,
                             SEG_CATEGORIES seg_cat_idx) {
  return seg_feature_data_max[seg_cat_idx][feature_id];
}

int av1_is_segfeature_signed(SEG_LVL_FEATURES feature_id,
                             SEG_CATEGORIES seg_cat_idx) {
  return seg_feature_data_signed[seg_cat_idx][feature_id];
}

void av1_set_segdata(struct segmentation *seg, int segment_id,
                     SEG_LVL_FEATURES feature_id, int seg_data,
                     SEG_CATEGORIES seg_cat_idx) {
  assert(seg_data <= seg_feature_data_max[seg_cat_idx][feature_id]);
  if (seg_data < 0) {
    assert(seg_feature_data_signed[seg_cat_idx][feature_id]);
    assert(-seg_data <= seg_feature_data_max[seg_cat_idx][feature_id]);
  }

  seg->feature_data[segment_id][feature_id] = seg_data;
}
int av1_extract_segment_id_from_map(int segment_id_in_map,
                                    SEG_CATEGORIES seg_cat_idx) {
  return seg_cat_idx == ACTIVE_SEG_IDX
             ? (segment_id_in_map & ((1 << MAX_LOG2_ACTIVE_SEGMENTS) - 1))
             : ((segment_id_in_map >> MAX_LOG2_ACTIVE_SEGMENTS) &
                ((1 << MAX_LOG2_QUALITY_SEGMENTS) - 1));
}
void av1_store_segment_id_into_map(int segment_id, uint8_t* segment_id_in_map,
                                   SEG_CATEGORIES seg_cat_idx) {
  if (seg_cat_idx == ACTIVE_SEG_IDX) {
    *segment_id_in_map =
        (*segment_id_in_map & ~((1 << MAX_LOG2_ACTIVE_SEGMENTS) - 1)) |
        (segment_id & ((1 << MAX_LOG2_ACTIVE_SEGMENTS) - 1));
  } else {
    *segment_id_in_map =
        (*segment_id_in_map & ((~((1 << MAX_LOG2_QUALITY_SEGMENTS) - 1))
                               << MAX_LOG2_ACTIVE_SEGMENTS)) |
        ((segment_id & ((1 << MAX_LOG2_QUALITY_SEGMENTS) - 1))
         << MAX_LOG2_ACTIVE_SEGMENTS);
  }
}
int av1_is_segfeature_enabled(const struct segmentation *seg,
                              SEG_LVL_FEATURES feature_id) {
  int is_segfeature_enabled = seg->enabled;
  int i;
  for (i = 0; i < seg->num_seg; ++i) {
    if (!(seg->feature_mask[i] & (1 << feature_id))) {
      is_segfeature_enabled = 0;
      break;
    }
  }
  return is_segfeature_enabled;
}
const aom_tree_index av1_segment_tree[MAX_SEGMENTS - 1][TREE_SIZE(
    MAX_SEGMENTS)] = { 
                       { 0, -1 },
                       { 0, 2, -1, -2 },
                       { 2, 4, 0, -1, -2, -3 },
                       { 2, 4, 0, -1, -2, 6, -3, -4 },
                       { 2, 4, 0, -1, 6, 8, -2, -3, -4, -5 },
                       { 2, 4, 0, 6, 8, 10, -1, -2, -3, -4, -5, -6 },
                       { 2, 4, 6, 8, 10, 12, 0, -1, -2, -3, -4, -5, -6, -7 }};
#else
static const int seg_feature_data_signed[SEG_LVL_MAX] = { 1, 1, 0, 0 };

static const int seg_feature_data_max[SEG_LVL_MAX] = { MAXQ, MAX_LOOP_FILTER, 3,
                                                       0 };

// These functions provide access to new segment level features.
// Eventually these function may be "optimized out" but for the moment,
// the coding mechanism is still subject to change so these provide a
// convenient single point of change.

void av1_clearall_segfeatures(struct segmentation *seg) {
  av1_zero(seg->feature_data);
  av1_zero(seg->feature_mask);
}

void av1_enable_segfeature(struct segmentation *seg, int segment_id,
                           SEG_LVL_FEATURES feature_id) {
  seg->feature_mask[segment_id] |= 1 << feature_id;
}

int av1_seg_feature_data_max(SEG_LVL_FEATURES feature_id) {
  return seg_feature_data_max[feature_id];
}

int av1_is_segfeature_signed(SEG_LVL_FEATURES feature_id) {
  return seg_feature_data_signed[feature_id];
}

void av1_set_segdata(struct segmentation *seg, int segment_id,
                     SEG_LVL_FEATURES feature_id, int seg_data) {
  assert(seg_data <= seg_feature_data_max[feature_id]);
  if (seg_data < 0) {
    assert(seg_feature_data_signed[feature_id]);
    assert(-seg_data <= seg_feature_data_max[feature_id]);
  }

  seg->feature_data[segment_id][feature_id] = seg_data;
}

const aom_tree_index av1_segment_tree[TREE_SIZE(MAX_SEGMENTS)] = {
  2, 4, 6, 8, 10, 12, 0, -1, -2, -3, -4, -5, -6, -7
};
#endif

// TBD? Functions to read and write segment data with range / validity checking
