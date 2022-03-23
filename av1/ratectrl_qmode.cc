/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "av1/ratectrl_qmode.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#include "av1/encoder/pass2_strategy.h"

namespace aom {

static GopFrame gop_frame_basic(int coding_idx, int order_idx,
                                bool is_key_frame, bool is_arf_frame,
                                bool is_golden_frame, bool is_show_frame) {
  GopFrame gop_frame;
  gop_frame.coding_idx = coding_idx;
  gop_frame.order_idx = order_idx;
  gop_frame.is_key_frame = is_key_frame;
  gop_frame.is_arf_frame = is_arf_frame;
  gop_frame.is_golden_frame = is_golden_frame;
  gop_frame.is_show_frame = is_show_frame;
  gop_frame.encode_ref_mode = EncodeRefMode::kRegular;
  gop_frame.colocated_ref_idx = -1;
  gop_frame.update_ref_idx = -1;
  return gop_frame;
}

// This function create gop frames with indices of display order from
// order_start to order_end - 1. The function will recursively introduce
// intermediate ARF untill maximum depth is met or the number of regular frames
// in between two ARFs are less than 3. Than the regular frames will be added
// into the gop_struct.
void construct_gop_multi_layer(GopStruct *gop_struct,
                               RefFrameManager *ref_frame_manager,
                               int max_depth, int depth, int order_start,
                               int order_end) {
  int coding_idx = static_cast<int>(gop_struct->gop_frame_list.size());
  GopFrame gop_frame;
  int num_frames = order_end - order_start;
  // If there are less than 3 frames, stop introducing ARF
  if (depth < max_depth && num_frames < 3) {
    int order_mid = (order_start + order_end) / 2;
    // intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 1, 0, 0);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kForward,
                                   EncodeRefMode::kRegular);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_start, order_mid);
    // show existing intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 0, 0, 1);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kNone,
                                   EncodeRefMode::kShowExisting);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_mid + 1, order_end);
  } else {
    // regular frame
    for (int i = order_start; i < order_end; ++i) {
      coding_idx = static_cast<int>(gop_struct->gop_frame_list.size());
      gop_frame = gop_frame_basic(coding_idx, i, 0, 0, 0, 1);
      ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kLast,
                                     EncodeRefMode::kRegular);
      gop_struct->gop_frame_list.push_back(gop_frame);
    }
  }
}

GopStruct construct_gop(RefFrameManager *ref_frame_manager,
                        int show_frame_count, bool has_key_frame) {
  GopStruct gop_struct;
  gop_struct.show_frame_count = show_frame_count;
  int order_start = 0;
  int order_arf = show_frame_count - 1;
  int coding_idx;
  GopFrame gop_frame;
  if (has_key_frame) {
    ref_frame_manager->Reset();
    coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
    gop_frame = gop_frame_basic(coding_idx, order_start, 1, 0, 1, 1);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kBackward,
                                   EncodeRefMode::kRegular);
    gop_struct.gop_frame_list.push_back(gop_frame);
    order_start++;
  }
  // ARF
  coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 1, 1, 0);
  ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kForward,
                                 EncodeRefMode::kRegular);
  gop_struct.gop_frame_list.push_back(gop_frame);
  construct_gop_multi_layer(&gop_struct, ref_frame_manager,
                            ref_frame_manager->ForwardMaxSize(), 1, order_start,
                            order_arf);
  // Overlay
  coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 0, 0, 1);
  ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kNone,
                                 EncodeRefMode::kOverlay);
  gop_struct.gop_frame_list.push_back(gop_frame);
  return gop_struct;
}

void AV1RateControlQMode::SetRcParam(const RateControlParam &rc_param) {
  rc_param_ = rc_param;
}

GopStructList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  std::vector<REGIONS> regions_list(MAX_FIRSTPASS_ANALYSIS_FRAMES);
  int total_regions = 0;
  // TODO(jianj): firstpass_info.size() should eventually be replaced
  // by the number of frames to the next KF.
  assert(firstpass_info.size() <= std::numeric_limits<int>::max());
  av1_identify_regions(firstpass_info.data(),
                       std::min(static_cast<int>(firstpass_info.size()),
                                MAX_FIRSTPASS_ANALYSIS_FRAMES),
                       0, regions_list.data(), &total_regions);

  // A temporary simple implementation
  const int max_gop_show_frame_count = 16;
  int remaining_show_frame_count = static_cast<int>(firstpass_info.size());
  GopStructList gop_list;

  RefFrameManager ref_frame_manager(rc_param_.max_ref_frames);

  while (remaining_show_frame_count > 0) {
    int show_frame_count =
        std::min(remaining_show_frame_count, max_gop_show_frame_count);
    // TODO(angiebird): determine gop show frame count based on first pass stats
    // here.
    bool has_key_frame = gop_list.size() == 0;
    GopStruct gop =
        construct_gop(&ref_frame_manager, show_frame_count, has_key_frame);
    gop_list.push_back(gop);
    remaining_show_frame_count -= show_frame_count;
  }
  return gop_list;
}

TplFrameDepStats create_tpl_frame_dep_stats_wo_propagation(
    const TplFrameStats &frame_stats) {
  const int frame_height = frame_stats.frame_height;
  const int frame_width = frame_stats.frame_width;
  const int min_block_size = frame_stats.min_block_size;
  const int unit_rows =
      frame_height / min_block_size + !!(frame_height % min_block_size);
  const int unit_cols =
      frame_width / min_block_size + !!(frame_width % min_block_size);
  TplFrameDepStats frame_dep_stats;
  frame_dep_stats.unit_size = min_block_size;
  frame_dep_stats.unit_stats = std::vector<std::vector<double>>(
      unit_rows, std::vector<double>(unit_cols, 0));
  for (const TplBlockStats &block_stats : frame_stats.block_stats_list) {
    int block_unit_rows = block_stats.height / min_block_size;
    int block_unit_cols = block_stats.width / min_block_size;
    int unit_count = block_unit_rows * block_unit_cols;
    int block_unit_row = block_stats.row / min_block_size;
    int block_unit_col = block_stats.col / min_block_size;
    double cost_diff =
        (block_stats.inter_cost - block_stats.intra_cost) * 1.0 / unit_count;
    for (int r = 0; r < block_unit_rows; r++) {
      for (int c = 0; c < block_unit_cols; c++) {
        frame_dep_stats.unit_stats[block_unit_row + r][block_unit_col + c] =
            cost_diff;
      }
    }
  }
  return frame_dep_stats;
}

int get_ref_coding_idx_list(const TplBlockStats &block_stats,
                            const RefFrameTable &ref_frame_table,
                            int *ref_coding_idx_list) {
  int ref_frame_count = 0;
  for (int i = 0; i < 2; ++i) {
    ref_coding_idx_list[i] = -1;
    int ref_frame_index = block_stats.ref_frame_index[i];
    if (ref_frame_index != -1) {
      ref_coding_idx_list[i] = ref_frame_table[ref_frame_index].coding_idx;
      ref_frame_count++;
    }
  }
  return ref_frame_count;
}

int get_block_overlap_area(int r0, int c0, int r1, int c1, int size) {
  int r_low = std::max(r0, r1);
  int r_high = std::min(r0 + size, r1 + size);
  int c_low = std::max(c0, c1);
  int c_high = std::min(c0 + size, c1 + size);
  return (r_high - r_low) * (c_high - c_low);
}

void tpl_frame_stats_propagate(const TplFrameStats &frame_stats,
                               const RefFrameTable &ref_frame_table,
                               TplGopDepStats *tpl_gop_dep_stats) {
  int min_block_size = frame_stats.min_block_size;
  for (const TplBlockStats &block_stats : frame_stats.block_stats_list) {
    int ref_coding_idx_list[2] = { -1, -1 };
    int ref_frame_count = get_ref_coding_idx_list(block_stats, ref_frame_table,
                                                  ref_coding_idx_list);
    if (ref_frame_count > 0) {
      double propagation_ratio = 1.0 / ref_frame_count;
      for (int i = 0; i < 2; ++i) {
        if (ref_coding_idx_list[i] != -1) {
          auto &ref_frame_dep_stats =
              tpl_gop_dep_stats->frame_dep_stats_list[ref_coding_idx_list[i]];
          int mv_row = block_stats.mv[i].row >> 3;
          int mv_col = block_stats.mv[i].col >> 3;
          int block_unit_rows = block_stats.height / min_block_size;
          int block_unit_cols = block_stats.width / min_block_size;
          int unit_count = block_unit_rows * block_unit_cols;
          double cost_diff =
              (block_stats.intra_cost - block_stats.inter_cost) / unit_count;
          for (int r = 0; r < block_unit_rows; r++) {
            for (int c = 0; c < block_unit_cols; c++) {
              int ref_block_row = block_stats.row + r * min_block_size + mv_row;
              int ref_block_col = block_stats.col + c * min_block_size + mv_col;
              int ref_unit_row_low = ref_block_row / min_block_size;
              int ref_unit_col_low = ref_block_col / min_block_size;
              for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                  int unit_row = ref_unit_row_low + j;
                  int unit_col = ref_unit_col_low + k;
                  int overlap_area = get_block_overlap_area(
                      unit_row * min_block_size, unit_col * min_block_size,
                      ref_block_row, ref_block_col, min_block_size);
                  double overlap_ratio =
                      overlap_area * 1.0 / (min_block_size * min_block_size);
                  ref_frame_dep_stats.unit_stats[unit_row][unit_col] +=
                      cost_diff * overlap_ratio * propagation_ratio;
                }
              }
            }
          }
        }
      }
    }
  }
}

void TplPropagation(const GopStruct &gop_struct,
                    const RefFrameTable &ref_frame_table,
                    const TplGopStats &tpl_gop_stats) {
  size_t frame_count = gop_struct.gop_frame_list.size();
  assert(tpl_gop_stats.frame_stats_list.size() == frame_count);

  // Generate ref_frame_table_list
  std::vector<RefFrameTable> ref_frame_table_list;
  ref_frame_table_list.push_back(ref_frame_table);
  for (size_t coding_idx = 0; coding_idx < frame_count; coding_idx++) {
    auto next_ref_frame_table = ref_frame_table;
    const auto &gop_frame = gop_struct.gop_frame_list[coding_idx];
    if (gop_frame.update_ref_idx != -1) {
      next_ref_frame_table[gop_frame.update_ref_idx] = gop_frame;
    }
    ref_frame_table_list.push_back(next_ref_frame_table);
  }

  // Create the struct to store TPL dependency stats
  TplGopDepStats tpl_gop_dep_stats;
  for (size_t coding_idx = 0; coding_idx < frame_count; coding_idx++) {
    TplFrameDepStats frame_dep_stats =
        create_tpl_frame_dep_stats_wo_propagation(
            tpl_gop_stats.frame_stats_list[coding_idx]);
    tpl_gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats);
  }

  // Back propagation
  for (int coding_idx = frame_count - 1; coding_idx >= 0; coding_idx--) {
    tpl_frame_stats_propagate(tpl_gop_stats.frame_stats_list[coding_idx],
                              ref_frame_table, &tpl_gop_dep_stats);
  }
}

std::vector<FrameEncodeParameters> AV1RateControlQMode::GetGopEncodeInfo(
    const GopStruct &gop_struct, const TplGopStats &tpl_stats_list) {
  std::vector<FrameEncodeParameters> frame_encoder_param;
  (void)tpl_stats_list;
  (void)gop_struct;
  return frame_encoder_param;
}

}  // namespace aom
