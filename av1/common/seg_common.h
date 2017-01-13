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

#ifndef AV1_COMMON_SEG_COMMON_H_
#define AV1_COMMON_SEG_COMMON_H_

#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif
#if CONFIG_EXT_SEGMENT
/****************************************
*seg common is defined for both active segments 
*and quality segments. So there is redundant in 
*structure definition, e.g. segment feature data defined
*for active segments.
****************************************/
#endif
#define SEGMENT_DELTADATA 0
#define SEGMENT_ABSDATA 1

#if CONFIG_EXT_SEGMENT
#define MAX_LOG2_ACTIVE_SEGMENTS 2
#define MAX_ACTIVE_SEGMENTS  (1<<MAX_LOG2_ACTIVE_SEGMENTS)
#define MAX_LOG2_QUALITY_SEGMENTS 4
#define MAX_QUALITY_SEGMENTS (1<<MAX_LOG2_QUALITY_SEGMENTS)
#define MAX_SEGMENTS  MAX_QUALITY_SEGMENTS
#else
#define MAX_SEGMENTS 8
#endif

#define SEG_TREE_PROBS (MAX_SEGMENTS - 1)

#define PREDICTION_PROBS 3
#if CONFIG_EXT_SEGMENT
#define HARMONIC_SPATIAL_PRED_PROBS 1
#define HETEROGENEOUS_SPATIAL_PRED_PROBS 2
#endif

// Segment level features.
#if CONFIG_EXT_SEGMENT
typedef enum {
ACTIVE_SEG_IDX = 0,
QUALITY_SEG_IDX = 1,
NUM_SEG_CATEGORIES
}SEG_CATEGORIES;
typedef enum {
  QUALITY_SEG_LVL_ALT_Q = 0,      // Use alternate Quantizer ....
  QUALITY_SEG_LVL_ALT_LF = 1,     // Use alternate loop filter value...    
  QUALITY_SEG_LVL_MAX
} QUALITY_SEG_LVL_FEATURES;
typedef enum {
  ACTIVE_SEG_LVL_REF_FRAME = 0,  // Optional Segment reference frame
  ACTIVE_SEG_LVL_SKIP = 1,       // Optional Segment (0,0) + skip mode 
  ACTIVE_SEG_LVL_MAX
} ACTIVE_SEG_LVL_FEATURES;
#define SEG_LVL_MAX  2    // Number of features (max of active segments and quality segments) supported
#define SEG_LVL_FEATURES uint8_t
#else
typedef enum {
  SEG_LVL_ALT_Q = 0,      // Use alternate Quantizer ....
  SEG_LVL_ALT_LF = 1,     // Use alternate loop filter value...
  SEG_LVL_REF_FRAME = 2,  // Optional Segment reference frame
  SEG_LVL_SKIP = 3,       // Optional Segment (0,0) + skip mode
  SEG_LVL_MAX = 4         // Number of features supported
} SEG_LVL_FEATURES;
#endif
struct segmentation {
  uint8_t enabled;
  uint8_t update_map;
  uint8_t update_data;
  uint8_t abs_delta;
  uint8_t temporal_update;
#if CONFIG_EXT_SEGMENT
  uint8_t seg_map_minb_size_log2_minus3;
  uint8_t num_seg;
#endif

  int16_t feature_data[MAX_SEGMENTS][SEG_LVL_MAX];
  unsigned int feature_mask[MAX_SEGMENTS];
};

struct segmentation_probs {
  aom_prob tree_probs[SEG_TREE_PROBS];
#if CONFIG_EC_MULTISYMBOL
  aom_cdf_prob tree_cdf[MAX_SEGMENTS];
#endif
#if CONFIG_EXT_SEGMENT
  aom_prob harmonic_spatial_pred_probs[HARMONIC_SPATIAL_PRED_PROBS];  
#if CONFIG_EC_MULTISYMBOL
  aom_cdf_prob heterogeneous_spatial_pred_cdf[HETEROGENEOUS_SPATIAL_PRED_PROBS + 1];
#endif
  aom_prob heterogeneous_spatial_pred_probs[HETEROGENEOUS_SPATIAL_PRED_PROBS];
#endif
  aom_prob pred_probs[PREDICTION_PROBS]; //temporal prediction

};

static INLINE int segfeature_active(const struct segmentation *seg,
                                    int segment_id,
                                    SEG_LVL_FEATURES feature_id) {
  return seg->enabled && (seg->feature_mask[segment_id] & (1 << feature_id));
}

void av1_clearall_segfeatures(struct segmentation *seg);

void av1_enable_segfeature(struct segmentation *seg, int segment_id,
                           SEG_LVL_FEATURES feature_id);
#if CONFIG_EXT_SEGMENT
int av1_seg_feature_data_max(SEG_LVL_FEATURES feature_id, SEG_CATEGORIES seg_cat_idx);

int av1_is_segfeature_signed(SEG_LVL_FEATURES feature_id, SEG_CATEGORIES seg_cat_idx);

void av1_set_segdata(struct segmentation *seg, int segment_id,
                     SEG_LVL_FEATURES feature_id, int seg_data, SEG_CATEGORIES seg_cat_idx);
#else
int av1_seg_feature_data_max(SEG_LVL_FEATURES feature_id);

int av1_is_segfeature_signed(SEG_LVL_FEATURES feature_id);

void av1_set_segdata(struct segmentation *seg, int segment_id,
                     SEG_LVL_FEATURES feature_id, int seg_data);
#endif

static INLINE int get_segdata(const struct segmentation *seg, int segment_id,
                              SEG_LVL_FEATURES feature_id) {
  return seg->feature_data[segment_id][feature_id];
}
#if CONFIG_EXT_SEGMENT
int av1_extract_segment_id_from_map(int segment_id_in_map, SEG_CATEGORIES seg_cat_idx);
void av1_store_segment_id_into_map(int segment_id, uint8_t* segment_id_in_map, SEG_CATEGORIES seg_cat_idx);
int av1_is_segfeature_enabled(const struct segmentation *seg, SEG_LVL_FEATURES feature_id);
extern const aom_tree_index av1_segment_tree[MAX_SEGMENTS - 1][TREE_SIZE(MAX_SEGMENTS)];
#else
extern const aom_tree_index av1_segment_tree[TREE_SIZE(MAX_SEGMENTS)];
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_SEG_COMMON_H_
