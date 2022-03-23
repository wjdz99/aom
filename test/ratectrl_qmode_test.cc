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

#include <memory>

#include "av1/ratectrl_qmode.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
namespace aom {

void test_gop_display_order(const GopStruct &gop_struct) {
  // Test whether show frames' order indices are sequential
  int ref_order_idx = 0;
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    if (gop_frame.is_show_frame) {
      EXPECT_EQ(gop_frame.order_idx, ref_order_idx);
      ref_order_idx++;
    }
  }
}

void test_colocated_show_frame(const GopStruct &gop_struct) {
  // Test whether each non show frame has a colocated show frame
  size_t gop_size = gop_struct.gop_frame_list.size();
  for (size_t gop_idx = 0; gop_idx < gop_size; ++gop_idx) {
    auto &gop_frame = gop_struct.gop_frame_list[gop_idx];
    if (gop_frame.is_show_frame == 0) {
      bool found_colocated_ref_frame = false;
      for (size_t i = gop_idx + 1; i < gop_size; ++i) {
        auto &next_gop_frame = gop_struct.gop_frame_list[i];
        if (gop_frame.order_idx == next_gop_frame.order_idx) {
          found_colocated_ref_frame = true;
          EXPECT_EQ(gop_frame.update_ref_idx, next_gop_frame.colocated_ref_idx);
          EXPECT_TRUE(next_gop_frame.is_show_frame);
        }
        if (gop_frame.update_ref_idx == next_gop_frame.update_ref_idx) {
          break;
        }
      }
      EXPECT_TRUE(found_colocated_ref_frame);
    }
  }
}

TEST(RateControlQModeTest, ConstructGopARF) {
  int show_frame_count = 16;
  const int max_ref_frames = 7;
  const bool has_key_frame = false;
  RefFrameManager ref_frame_manager(max_ref_frames);
  GopStruct gop_struct =
      construct_gop(&ref_frame_manager, show_frame_count, has_key_frame);
  test_gop_display_order(gop_struct);
  test_colocated_show_frame(gop_struct);
}

TEST(RateControlQModeTest, ConstructGopKey) {
  int show_frame_count = 16;
  int max_ref_frames = 7;
  int has_key_frame = 1;
  RefFrameManager ref_frame_manager(max_ref_frames);
  GopStruct gop_struct =
      construct_gop(&ref_frame_manager, show_frame_count, has_key_frame);
  test_gop_display_order(gop_struct);
  test_colocated_show_frame(gop_struct);
}

static TplBlockStats create_toy_tpl_block_stats(int h, int w, int r, int c,
                                                int cost_diff) {
  TplBlockStats tpl_block_stats;
  tpl_block_stats.height = h;
  tpl_block_stats.width = w;
  tpl_block_stats.row = r;
  tpl_block_stats.col = c;
  // A random trick that makes inter_cost - intra_cost = cost_diff;
  tpl_block_stats.intra_cost = cost_diff / 2;
  tpl_block_stats.inter_cost = cost_diff + cost_diff / 2;
  tpl_block_stats.mv[0] = { 0, 0 };
  tpl_block_stats.mv[1] = { 0, 0 };
  tpl_block_stats.ref_frame_index[0] = -1;
  tpl_block_stats.ref_frame_index[1] = -1;
  return tpl_block_stats;
}

static TplFrameStats create_toy_tpl_frame_stats() {
  TplFrameStats frame_stats;
  int max_h = 16;
  int max_w = max_h;
  int count = 2;
  frame_stats.min_block_size = max_h >> (count - 1);
  frame_stats.frame_height = max_h * 2;
  frame_stats.frame_width = max_w * 2;
  std::srand(0);
  for (int i = 0; i < count; ++i) {
    for (int j = 0; j < count; ++j) {
      int h = max_h >> i;
      int w = max_w >> j;
      for (int u = 0; u * h < max_h; ++u) {
        for (int v = 0; v * w < max_w; ++v) {
          int r = max_h * i + h * u;
          int c = max_w * j + w * v;
          int cost_diff = std::rand() % 16;
          TplBlockStats block_stats =
              create_toy_tpl_block_stats(h, w, r, c, cost_diff);
          frame_stats.block_stats_list.push_back(block_stats);
        }
      }
    }
  }
  return frame_stats;
}

TEST(RateControlQModeTest, CreateTplFrameDepStats) {
  TplFrameStats frame_stats = create_toy_tpl_frame_stats();
  TplFrameDepStats frame_dep_stats =
      create_tpl_frame_dep_stats_wo_propagation(frame_stats);
  EXPECT_EQ(frame_stats.min_block_size, frame_dep_stats.unit_size);
  int unit_rows = frame_dep_stats.unit_stats.size();
  int unit_cols = frame_dep_stats.unit_stats[0].size();
  EXPECT_EQ(frame_stats.frame_height, unit_rows * frame_dep_stats.unit_size);
  EXPECT_EQ(frame_stats.frame_width, unit_cols * frame_dep_stats.unit_size);
  double sum_cost_diff = 0;
  for (int r = 0; r < unit_rows; ++r) {
    for (int c = 0; c < unit_cols; ++c) {
      sum_cost_diff += frame_dep_stats.unit_stats[r][c];
    }
  }

  double ref_sum_cost_diff = 0;
  for (auto &block_stats : frame_stats.block_stats_list) {
    ref_sum_cost_diff += block_stats.inter_cost - block_stats.intra_cost;
  }
  EXPECT_NEAR(sum_cost_diff, ref_sum_cost_diff, 0.0000001);
}

}  // namespace aom

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
