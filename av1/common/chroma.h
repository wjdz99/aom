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

#ifndef AOM_AV1_COMMON_CHROMA_H_
#define AOM_AV1_COMMON_CHROMA_H_

#include <assert.h>
#include <stdbool.h>

#include "config/aom_config.h"

#include "av1/common/common_data.h"
#include "av1/common/enums.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CHROMA_REF_INFO {
  int is_chroma_ref;
  int offset_started;
  int mi_row_chroma_base;
  int mi_col_chroma_base;
  BLOCK_SIZE bsize;
  BLOCK_SIZE bsize_base;
} CHROMA_REF_INFO;

static INLINE void initialize_chr_ref_info(int mi_row, int mi_col,
                                           BLOCK_SIZE bsize,
                                           CHROMA_REF_INFO *info) {
  info->is_chroma_ref = 1;
  info->offset_started = 0;
  info->mi_row_chroma_base = mi_row;
  info->mi_col_chroma_base = mi_col;
  info->bsize = bsize;
  info->bsize_base = bsize;
}

// Decide whether a block needs coding multiple chroma coding blocks in it at
// once to get around sub-4x4 coding.
static INLINE int have_nz_chroma_ref_offset(BLOCK_SIZE bsize,
                                            PARTITION_TYPE partition,
                                            int subsampling_x,
                                            int subsampling_y) {
  const int bw = block_size_wide[bsize] >> subsampling_x;
  const int bh = block_size_high[bsize] >> subsampling_y;
  const int bw_less_than_4 = bw < 4;
  const int bh_less_than_4 = bh < 4;
  const int hbw_less_than_4 = bw < 8;
  const int hbh_less_than_4 = bh < 8;
  const int qbw_less_than_4 = bw < 16;
  const int qbh_less_than_4 = bh < 16;

  switch (partition) {
    case PARTITION_NONE: return bw_less_than_4 || bh_less_than_4;
    case PARTITION_HORZ: return bw_less_than_4 || hbh_less_than_4;
    case PARTITION_VERT: return hbw_less_than_4 || bh_less_than_4;
    case PARTITION_SPLIT: return hbw_less_than_4 || hbh_less_than_4;
#if !CONFIG_EXT_RECUR_PARTITIONS
    case PARTITION_HORZ_A:
    case PARTITION_HORZ_B:
    case PARTITION_VERT_A:
    case PARTITION_VERT_B: return hbw_less_than_4 || hbh_less_than_4;
#endif  // !CONFIG_EXT_RECUR_PARTITIONS
#if CONFIG_EXT_PARTITIONS
    case PARTITION_HORZ_3: return bw_less_than_4 || qbh_less_than_4;
    case PARTITION_VERT_3: return qbw_less_than_4 || bh_less_than_4;
#else
    case PARTITION_HORZ_4: return bw_less_than_4 || qbh_less_than_4;
    case PARTITION_VERT_4: return qbw_less_than_4 || bh_less_than_4;
#endif  // CONFIG_EXT_PARTITIONS
    default:
      assert(0 && "Invalid partition type!");
      return 0;
      break;
  }
}

// Decide whether a subblock is the main chroma reference when its parent block
// needs coding multiple chroma coding blocks at once. The function returns a
// flag indicating whether the mode info used for the combined chroma block is
// located in the subblock.
static INLINE int is_sub_partition_chroma_ref(PARTITION_TYPE partition,
                                              int index, BLOCK_SIZE bsize,
                                              BLOCK_SIZE parent_bsize, int ss_x,
                                              int ss_y, int is_offset_started) {
  (void)is_offset_started;
  (void)parent_bsize;
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  const int pw = bw >> ss_x;
  const int ph = bh >> ss_y;
  const int pw_less_than_4 = pw < 4;
  const int ph_less_than_4 = ph < 4;

  switch (partition) {
    case PARTITION_NONE: return 1;
    case PARTITION_HORZ:
    case PARTITION_VERT: return index == 1;
    case PARTITION_SPLIT:
      if (is_offset_started) {
        return index == 3;
      } else {
        if (pw_less_than_4 && ph_less_than_4)
          return index == 3;
        else if (pw_less_than_4)
          return index == 1 || index == 3;
        else if (ph_less_than_4)
          return index == 2 || index == 3;
        else
          return 1;
      }
#if !CONFIG_EXT_RECUR_PARTITIONS
    case PARTITION_HORZ_A:
    case PARTITION_HORZ_B:
    case PARTITION_VERT_A:
    case PARTITION_VERT_B:
      if (is_offset_started) {
        return index == 2;
      } else {
        const int smallest_w = block_size_wide[parent_bsize] >> (ss_x + 1);
        const int smallest_h = block_size_high[parent_bsize] >> (ss_y + 1);
        const int smallest_w_less_than_4 = smallest_w < 4;
        const int smallest_h_less_than_4 = smallest_h < 4;

        if (smallest_w_less_than_4 && smallest_h_less_than_4) {
          return index == 2;
        } else if (smallest_w_less_than_4) {
          if (partition == PARTITION_VERT_A || partition == PARTITION_VERT_B) {
            return index == 2;
          } else if (partition == PARTITION_HORZ_A) {
            return index == 1 || index == 2;
          } else {
            return index == 0 || index == 2;
          }
        } else if (smallest_h_less_than_4) {
          if (partition == PARTITION_HORZ_A || partition == PARTITION_HORZ_B) {
            return index == 2;
          } else if (partition == PARTITION_VERT_A) {
            return index == 1 || index == 2;
          } else {
            return index == 0 || index == 2;
          }
        } else {
          return 1;
        }
      }
#endif  // !CONFIG_EXT_RECUR_PARTITIONS
#if CONFIG_EXT_PARTITIONS
    case PARTITION_VERT_3:
    case PARTITION_HORZ_3: return index == 2;
#else
    case PARTITION_HORZ_4:
    case PARTITION_VERT_4:
      if (is_offset_started) {
        return index == 3;
      } else {
        if ((partition == PARTITION_HORZ_4 && ph_less_than_4) ||
            (partition == PARTITION_VERT_4 && pw_less_than_4)) {
          return index == 1 || index == 3;
        } else {
          return 1;
        }
      }
#endif  // CONFIG_EXT_PARTITIONS
    default:
      assert(0 && "Invalid partition type!");
      return 0;
      break;
  }
}

static INLINE void set_chroma_ref_offset_size(int mi_row, int mi_col,
                                              PARTITION_TYPE partition,
                                              BLOCK_SIZE bsize,
                                              BLOCK_SIZE parent_bsize, int ss_x,
                                              int ss_y, CHROMA_REF_INFO *info,
                                              CHROMA_REF_INFO *parent_info) {
  const int pw = block_size_wide[bsize] >> ss_x;
  const int ph = block_size_high[bsize] >> ss_y;
  const int pw_less_than_4 = pw < 4;
  const int ph_less_than_4 = ph < 4;
#if !CONFIG_EXT_RECUR_PARTITIONS
  const int hppw = block_size_wide[parent_bsize] >> (ss_x + 1);
  const int hpph = block_size_high[parent_bsize] >> (ss_y + 1);
  const int hppw_less_than_4 = hppw < 4;
  const int hpph_less_than_4 = hpph < 4;
#endif  // !CONFIG_EXT_RECUR_PARTITIONS
#if !CONFIG_EXT_PARTITIONS
  const int mi_row_mid_point =
      parent_info->mi_row_chroma_base + (mi_size_high[parent_bsize] >> 1);
  const int mi_col_mid_point =
      parent_info->mi_col_chroma_base + (mi_size_wide[parent_bsize] >> 1);
#endif  // !CONFIG_EXT_PARTITIONS

  assert(parent_info->offset_started == 0);

  switch (partition) {
    case PARTITION_NONE:
    case PARTITION_HORZ:
    case PARTITION_VERT:
#if CONFIG_EXT_PARTITIONS
    case PARTITION_VERT_3:
    case PARTITION_HORZ_3:
#endif  // CONFIG_EXT_PARTITIONS
      info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
      info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
      info->bsize_base = parent_bsize;
      break;
    case PARTITION_SPLIT:
      if (pw_less_than_4 && ph_less_than_4) {
        info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        info->bsize_base = parent_bsize;
      } else if (pw_less_than_4) {
        info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_HORZ);
        info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        if (mi_row == parent_info->mi_row_chroma_base) {
          info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        } else {
          info->mi_row_chroma_base =
              parent_info->mi_row_chroma_base + mi_size_high[bsize];
        }
      } else {
        assert(ph_less_than_4);
        info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_VERT);
        info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        if (mi_col == parent_info->mi_col_chroma_base) {
          info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        } else {
          info->mi_col_chroma_base =
              parent_info->mi_col_chroma_base + mi_size_wide[bsize];
        }
      }
      break;
#if !CONFIG_EXT_RECUR_PARTITIONS
    case PARTITION_HORZ_A:
    case PARTITION_HORZ_B:
    case PARTITION_VERT_A:
    case PARTITION_VERT_B:
      if ((hppw_less_than_4 && hpph_less_than_4) ||
          (hppw_less_than_4 &&
           (partition == PARTITION_VERT_A || partition == PARTITION_VERT_B)) ||
          (hpph_less_than_4 &&
           (partition == PARTITION_HORZ_A || partition == PARTITION_HORZ_B))) {
        info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        info->bsize_base = parent_bsize;
      } else if (hppw_less_than_4) {
        info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_HORZ);
        info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        if (mi_row == parent_info->mi_row_chroma_base) {
          info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        } else {
          info->mi_row_chroma_base = parent_info->mi_row_chroma_base +
                                     (mi_size_high[parent_bsize] >> 1);
        }
      } else {
        assert(hpph_less_than_4);

        info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_VERT);
        info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
        if (mi_col == parent_info->mi_col_chroma_base) {
          info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
        } else {
          info->mi_col_chroma_base = parent_info->mi_col_chroma_base +
                                     (mi_size_wide[parent_bsize] >> 1);
        }
      }
      break;
#endif  // !CONFIG_EXT_RECUR_PARTITIONS
#if !CONFIG_EXT_PARTITIONS
    case PARTITION_HORZ_4:
      info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_HORZ);
      info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
      if (mi_row < mi_row_mid_point) {
        info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
      } else {
        info->mi_row_chroma_base = mi_row_mid_point;
      }
      break;
    case PARTITION_VERT_4:
      info->bsize_base = get_partition_subsize(parent_bsize, PARTITION_VERT);
      info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
      if (mi_col < mi_col_mid_point) {
        info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
      } else {
        info->mi_col_chroma_base = mi_col_mid_point;
      }
      break;
#endif  // !CONFIG_EXT_PARTITIONS
    default: assert(0 && "Invalid partition type!"); break;
  }
}

static INLINE void set_chroma_ref_info(int mi_row, int mi_col, int index,
                                       BLOCK_SIZE bsize, CHROMA_REF_INFO *info,
                                       CHROMA_REF_INFO *parent_info,
                                       BLOCK_SIZE parent_bsize,
                                       PARTITION_TYPE parent_partition,
                                       int ss_x, int ss_y) {
  assert(bsize < BLOCK_SIZES_ALL);
  initialize_chr_ref_info(mi_row, mi_col, bsize, info);

  if (parent_info == NULL) return;

  if (parent_info->is_chroma_ref) {
    if (parent_info->offset_started) {
      if (is_sub_partition_chroma_ref(parent_partition, index, bsize,
                                      parent_bsize, ss_x, ss_y, 1)) {
        info->is_chroma_ref = 1;
      } else {
        info->is_chroma_ref = 0;
      }
      info->offset_started = 1;
      info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
      info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
      info->bsize_base = parent_info->bsize_base;
    } else if (have_nz_chroma_ref_offset(parent_bsize, parent_partition, ss_x,
                                         ss_y)) {
      info->offset_started = 1;
      info->is_chroma_ref = is_sub_partition_chroma_ref(
          parent_partition, index, bsize, parent_bsize, ss_x, ss_y, 0);
      set_chroma_ref_offset_size(mi_row, mi_col, parent_partition, bsize,
                                 parent_bsize, ss_x, ss_y, info, parent_info);
    }
  } else {
    info->is_chroma_ref = 0;
    info->offset_started = 1;
    info->mi_row_chroma_base = parent_info->mi_row_chroma_base;
    info->mi_col_chroma_base = parent_info->mi_col_chroma_base;
    info->bsize_base = parent_info->bsize_base;
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_CHROMA_H_
