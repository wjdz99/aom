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

#include "aom_ports/system_state.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/level.h"

#define UNDEFINED_LEVEL                                                 \
  {                                                                     \
    .level = SEQ_LEVEL_MAX, .max_picture_size = 0, .max_h_size = 0,     \
    .max_v_size = 0, .max_display_rate = 0, .max_decode_rate = 0,       \
    .max_header_rate = 0, .main_mbps = 0, .high_mbps = 0, .main_cr = 0, \
    .high_cr = 0, .max_tiles = 0, .max_tile_cols = 0                    \
  }

static const AV1LevelSpec av1_level_defs[SEQ_LEVELS] = {
  { .level = SEQ_LEVEL_2_0,
    .max_picture_size = 147456,
    .max_h_size = 2048,
    .max_v_size = 1152,
    .max_display_rate = 4423680L,
    .max_decode_rate = 5529600L,
    .max_header_rate = 150,
    .main_mbps = 1.5,
    .high_mbps = 0,
    .main_cr = 2.0,
    .high_cr = 0,
    .max_tiles = 8,
    .max_tile_cols = 4 },
  { .level = SEQ_LEVEL_2_1,
    .max_picture_size = 278784,
    .max_h_size = 2816,
    .max_v_size = 1584,
    .max_display_rate = 8363520L,
    .max_decode_rate = 10454400L,
    .max_header_rate = 150,
    .main_mbps = 3.0,
    .high_mbps = 0,
    .main_cr = 2.0,
    .high_cr = 0,
    .max_tiles = 8,
    .max_tile_cols = 4 },
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
  { .level = SEQ_LEVEL_3_0,
    .max_picture_size = 665856,
    .max_h_size = 4352,
    .max_v_size = 2448,
    .max_display_rate = 19975680L,
    .max_decode_rate = 24969600L,
    .max_header_rate = 150,
    .main_mbps = 6.0,
    .high_mbps = 0,
    .main_cr = 2.0,
    .high_cr = 0,
    .max_tiles = 16,
    .max_tile_cols = 6 },
  { .level = SEQ_LEVEL_3_1,
    .max_picture_size = 1065024,
    .max_h_size = 5504,
    .max_v_size = 3096,
    .max_display_rate = 31950720L,
    .max_decode_rate = 39938400L,
    .max_header_rate = 150,
    .main_mbps = 10.0,
    .high_mbps = 0,
    .main_cr = 2.0,
    .high_cr = 0,
    .max_tiles = 16,
    .max_tile_cols = 6 },
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
  { .level = SEQ_LEVEL_4_0,
    .max_picture_size = 2359296,
    .max_h_size = 6144,
    .max_v_size = 3456,
    .max_display_rate = 70778880L,
    .max_decode_rate = 77856768L,
    .max_header_rate = 300,
    .main_mbps = 12.0,
    .high_mbps = 30.0,
    .main_cr = 4.0,
    .high_cr = 4.0,
    .max_tiles = 32,
    .max_tile_cols = 8 },
  { .level = SEQ_LEVEL_4_1,
    .max_picture_size = 2359296,
    .max_h_size = 6144,
    .max_v_size = 3456,
    .max_display_rate = 141557760L,
    .max_decode_rate = 155713536L,
    .max_header_rate = 300,
    .main_mbps = 20.0,
    .high_mbps = 50.0,
    .main_cr = 4.0,
    .high_cr = 4.0,
    .max_tiles = 32,
    .max_tile_cols = 8 },
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
  { .level = SEQ_LEVEL_5_0,
    .max_picture_size = 8912896,
    .max_h_size = 8192,
    .max_v_size = 4352,
    .max_display_rate = 267386880L,
    .max_decode_rate = 273715200L,
    .max_header_rate = 300,
    .main_mbps = 30.0,
    .high_mbps = 100.0,
    .main_cr = 6.0,
    .high_cr = 4.0,
    .max_tiles = 64,
    .max_tile_cols = 8 },
  { .level = SEQ_LEVEL_5_1,
    .max_picture_size = 8912896,
    .max_h_size = 8192,
    .max_v_size = 4352,
    .max_display_rate = 534773760L,
    .max_decode_rate = 547430400L,
    .max_header_rate = 300,
    .main_mbps = 40.0,
    .high_mbps = 160.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 64,
    .max_tile_cols = 8 },
  { .level = SEQ_LEVEL_5_2,
    .max_picture_size = 8912896,
    .max_h_size = 8192,
    .max_v_size = 4352,
    .max_display_rate = 1069547520L,
    .max_decode_rate = 1094860800L,
    .max_header_rate = 300,
    .main_mbps = 60.0,
    .high_mbps = 240.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 64,
    .max_tile_cols = 8 },
  { .level = SEQ_LEVEL_5_3,
    .max_picture_size = 8912896,
    .max_h_size = 8192,
    .max_v_size = 4352,
    .max_display_rate = 1069547520L,
    .max_decode_rate = 1176502272L,
    .max_header_rate = 300,
    .main_mbps = 60.0,
    .high_mbps = 240.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 64,
    .max_tile_cols = 8 },
  { .level = SEQ_LEVEL_6_0,
    .max_picture_size = 35651584,
    .max_h_size = 16384,
    .max_v_size = 8704,
    .max_display_rate = 1069547520L,
    .max_decode_rate = 1176502272L,
    .max_header_rate = 300,
    .main_mbps = 60.0,
    .high_mbps = 240.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 128,
    .max_tile_cols = 16 },
  { .level = SEQ_LEVEL_6_1,
    .max_picture_size = 35651584,
    .max_h_size = 16384,
    .max_v_size = 8704,
    .max_display_rate = 2139095040L,
    .max_decode_rate = 2189721600L,
    .max_header_rate = 300,
    .main_mbps = 100.0,
    .high_mbps = 480.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 128,
    .max_tile_cols = 16 },
  { .level = SEQ_LEVEL_6_2,
    .max_picture_size = 35651584,
    .max_h_size = 16384,
    .max_v_size = 8704,
    .max_display_rate = 4278190080L,
    .max_decode_rate = 4379443200L,
    .max_header_rate = 300,
    .main_mbps = 160.0,
    .high_mbps = 800.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 128,
    .max_tile_cols = 16 },
  { .level = SEQ_LEVEL_6_3,
    .max_picture_size = 35651584,
    .max_h_size = 16384,
    .max_v_size = 8704,
    .max_display_rate = 4278190080L,
    .max_decode_rate = 4706009088L,
    .max_header_rate = 300,
    .main_mbps = 160.0,
    .high_mbps = 800.0,
    .main_cr = 8.0,
    .high_cr = 4.0,
    .max_tiles = 128,
    .max_tile_cols = 16 },
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
  UNDEFINED_LEVEL,
};

typedef enum {
  LUMA_PIC_SIZE_TOO_LARGE,
  LUMA_PIC_H_SIZE_TOO_LARGE,
  LUMA_PIC_V_SIZE_TOO_LARGE,
  TOO_MANY_TILE_COLUMNS,
  TOO_MANY_TILES,
  TILE_TOO_LARGE,
  CROPPED_TILE_WIDTH_TOO_SMALL,
  CROPPED_TILE_HEIGHT_TOO_SMALL,
  TILE_WIDTH_INVALID,

  TARGET_LEVEL_FAIL_IDS
} TARGET_LEVEL_FAIL_ID;

static const char *level_fail_messages[TARGET_LEVEL_FAIL_IDS] = {
  "The picture size is too large.",
  "The picture width is too large.",
  "The picture height is too large.",
  "Too many tile columns are used.",
  "Too many tiles are used.",
  "The tile size is too large.",
  "The cropped tile width is less than 8",
  "The cropped tile height is less than 8",
  "The tile width is invalid",
};

static void check_level_constraints(AV1_COMP *cpi, int operating_point_idx,
                                    const AV1LevelSpec *level_spec) {
  const AV1_LEVEL target_seq_level_idx =
      cpi->target_seq_level_idx[operating_point_idx];
  if (target_seq_level_idx >= SEQ_LEVELS) return;
  TARGET_LEVEL_FAIL_ID fail_id = TARGET_LEVEL_FAIL_IDS;
  // Check level conformance
  // TODO(kyslov@) implement all constraints
  do {
    const AV1LevelSpec *const target_level_spec =
        av1_level_defs + target_seq_level_idx;
    if (level_spec->max_picture_size > target_level_spec->max_picture_size) {
      fail_id = LUMA_PIC_SIZE_TOO_LARGE;
      break;
    }

    if (level_spec->max_h_size > target_level_spec->max_h_size) {
      fail_id = LUMA_PIC_H_SIZE_TOO_LARGE;
      break;
    }

    if (level_spec->max_v_size > target_level_spec->max_v_size) {
      fail_id = LUMA_PIC_V_SIZE_TOO_LARGE;
      break;
    }

    if (level_spec->max_tile_cols > target_level_spec->max_tile_cols) {
      fail_id = TOO_MANY_TILE_COLUMNS;
      break;
    }

    if (level_spec->max_tiles > target_level_spec->max_tiles) {
      fail_id = TOO_MANY_TILES;
      break;
    }

    if (level_spec->max_tile_size > 4096 * 2304) {
      fail_id = TILE_TOO_LARGE;
      break;
    }

    if (level_spec->min_cropped_tile_width < 8) {
      fail_id = CROPPED_TILE_WIDTH_TOO_SMALL;
      break;
    }

    if (level_spec->min_cropped_tile_height < 8) {
      fail_id = CROPPED_TILE_HEIGHT_TOO_SMALL;
      break;
    }

    if (!level_spec->tile_width_is_valid) {
      fail_id = TILE_WIDTH_INVALID;
      break;
    }
  } while (0);

  if (fail_id != TARGET_LEVEL_FAIL_IDS) {
    AV1_COMMON *const cm = &cpi->common;
    const int target_level_major = 2 + (target_seq_level_idx >> 2);
    const int target_level_minor = target_seq_level_idx & 3;
    aom_internal_error(&cm->error, AOM_CODEC_ERROR,
                       "Failed to encode to the target level %d_%d. %s",
                       target_level_major, target_level_minor,
                       level_fail_messages[fail_id]);
  }
}

static INLINE int is_in_operating_point(int operating_point,
                                        int temporal_layer_id,
                                        int spatial_layer_id) {
  if (!operating_point) return 1;

  return ((operating_point >> temporal_layer_id) & 1) &&
         ((operating_point >> (spatial_layer_id + 8)) & 1);
}

static void get_tile_stats(const AV1_COMP *const cpi, int *max_tile_size,
                           int *min_cropped_tile_width,
                           int *min_cropped_tile_height,
                           int *tile_width_valid) {
  const AV1_COMMON *const cm = &cpi->common;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;

  *max_tile_size = 0;
  *min_cropped_tile_width = INT_MAX;
  *min_cropped_tile_height = INT_MAX;
  *tile_width_valid = 1;

  for (int tile_row = 0; tile_row < tile_rows; ++tile_row) {
    for (int tile_col = 0; tile_col < tile_cols; ++tile_col) {
      const TileInfo *const tile_info =
          &cpi->tile_data[tile_row * cm->tile_cols + tile_col].tile_info;
      const int tile_width =
          (tile_info->mi_col_end - tile_info->mi_col_start) * MI_SIZE;
      const int tile_height =
          (tile_info->mi_row_end - tile_info->mi_row_start) * MI_SIZE;
      const int tile_size = tile_width * tile_height;
      *max_tile_size = AOMMAX(*max_tile_size, tile_size);

      const int cropped_tile_width =
          cm->width - tile_info->mi_col_start * MI_SIZE;
      const int cropped_tile_height =
          cm->height - tile_info->mi_row_start * MI_SIZE;
      *min_cropped_tile_width =
          AOMMIN(*min_cropped_tile_width, cropped_tile_width);
      *min_cropped_tile_height =
          AOMMIN(*min_cropped_tile_height, cropped_tile_height);

      const int is_right_most_tile = tile_info->mi_col_end == cm->mi_cols;
      if (!is_right_most_tile) {
        if (av1_superres_scaled(cm))
          *tile_width_valid &= tile_width >= 128;
        else
          *tile_width_valid &= tile_width >= 64;
      }
    }
  }
}

void av1_update_level_info(AV1_COMP *cpi, size_t size, int64_t ts_start,
                           int64_t ts_end) {
  const AV1_COMMON *const cm = &cpi->common;
  const int upscaled_width = cm->superres_upscaled_width;
  const int height = cm->height;
  const int tile_cols = cm->tile_cols;
  const int tile_rows = cm->tile_rows;
  const int tiles = tile_cols * tile_rows;
  const int luma_pic_size = upscaled_width * height;

  // Store info. of current frame into FrameWindowBuffer.
  int new_idx = 0;
  FrameWindowBuffer *const buffer = &cpi->frame_window_buffer;
  if (cpi->frame_window_buffer.num < FRAME_WINDOW_SIZE) {
    new_idx = (buffer->start + buffer->num++) % FRAME_WINDOW_SIZE;
  } else {
    new_idx = buffer->start;
    buffer->start = (new_idx + 1) % FRAME_WINDOW_SIZE;
  }
  FrameRecord *const record = &buffer->buf[new_idx];
  record->ts_start = ts_start;
  record->ts_end = ts_end;
  record->pic_size = luma_pic_size;
  record->frame_header_count = cpi->frame_header_count;
  record->show_frame = cm->show_frame;
  record->show_existing_frame = cm->show_existing_frame;

  int max_tile_size;
  int min_cropped_tile_width;
  int min_cropped_tile_height;
  int tile_width_is_valid;
  get_tile_stats(cpi, &max_tile_size, &min_cropped_tile_width,
                 &min_cropped_tile_height, &tile_width_is_valid);

  const int pic_size_profile_factor =
      cm->seq_params.profile == PROFILE_0
          ? 15
          : (cm->seq_params.profile == PROFILE_1 ? 30 : 36);
  const size_t frame_compressed_size = (size > 129 ? size - 128 : 1);
  const size_t frame_uncompressed_size =
      (luma_pic_size * pic_size_profile_factor) >> 3;

  aom_clear_system_state();
  const double compression_ratio =
      frame_uncompressed_size / (double)frame_compressed_size;
  const double total_time_encoded =
      (cpi->last_end_time_stamp_seen - cpi->first_time_stamp_ever) /
      (double)TICKS_PER_SEC;

  const SequenceHeader *const seq = &cm->seq_params;
  const int temporal_layer_id = cm->temporal_layer_id;
  const int spatial_layer_id = cm->spatial_layer_id;
  // update level_stats
  // TODO(kyslov@) fix the implementation according to buffer model
  for (int i = 0; i < seq->operating_points_cnt_minus_1 + 1; ++i) {
    if (!is_in_operating_point(seq->operating_point_idc[i], temporal_layer_id,
                               spatial_layer_id)) {
      continue;
    }

    AV1LevelInfo *const level_info = &cpi->level_info[i];
    AV1LevelStats *const level_stats = &level_info->level_stats;

    level_stats->total_compressed_size += frame_compressed_size;
    if (cm->show_frame) level_stats->total_time_encoded = total_time_encoded;

    // update level_spec
    // TODO(kyslov@) update all spec fields
    AV1LevelSpec *const level_spec = &level_info->level_spec;
    level_spec->max_picture_size =
        AOMMAX(level_spec->max_picture_size, luma_pic_size);
    level_spec->max_h_size =
        AOMMAX(level_spec->max_h_size, cm->superres_upscaled_width);
    level_spec->max_v_size = AOMMAX(level_spec->max_v_size, height);
    level_spec->max_tile_cols = AOMMAX(level_spec->max_tile_cols, tile_cols);
    level_spec->max_tiles = AOMMAX(level_spec->max_tiles, tiles);
    level_spec->max_tile_size =
        AOMMAX(level_spec->max_tile_size, max_tile_size);
    level_spec->min_cropped_tile_width =
        AOMMIN(level_spec->min_cropped_tile_width, min_cropped_tile_width);
    level_spec->min_cropped_tile_height =
        AOMMIN(level_spec->min_cropped_tile_height, min_cropped_tile_height);
    level_spec->tile_width_is_valid &= tile_width_is_valid;

    // TODO(kyslov@) These are needed for further level stat calculations
    (void)compression_ratio;
    (void)ts_start;
    (void)ts_end;

    check_level_constraints(cpi, i, level_spec);
  }
}

aom_codec_err_t av1_get_seq_level_idx(const AV1_COMP *cpi, int *seq_level_idx) {
  // TODO(chiyotsai@, kyslov@, huisu@): put in real implementations.
  (void)cpi;
  for (int i = 0; i < MAX_NUM_OPERATING_POINTS; ++i) {
    seq_level_idx[i] = (int)SEQ_LEVEL_MAX;
  }

  return AOM_CODEC_OK;
}
