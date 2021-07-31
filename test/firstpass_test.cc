/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdio.h>
#include <stdlib.h>

#include "av1/common/common.h"
#include "av1/encoder/firstpass.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {
constexpr int ref_static_buf_size = FIRSTPASS_INFO_STATIC_BUF_SIZE;

TEST(FirstpassTest, FirstpassInfoInitWithExtBuf) {
  FIRSTPASS_INFO firstpass_info;
  FIRSTPASS_STATS ext_stats_buf[10];
  const int ref_stats_size = 10;
  for (int i = 0; i < ref_stats_size; ++i) {
    av1_zero(ext_stats_buf[i]);
    ext_stats_buf[i].frame = i;
  }
  aom_codec_err_t ret =
      av1_firstpass_info_init(&firstpass_info, ext_stats_buf, 10);
  EXPECT_EQ(firstpass_info.stats_count, ref_stats_size);
  EXPECT_EQ(firstpass_info.curr_stats_count_future +
                firstpass_info.curr_stats_count_past,
            firstpass_info.stats_count);
  EXPECT_EQ(firstpass_info.curr_index, 0);
  EXPECT_EQ(ret, AOM_CODEC_OK);
}

TEST(FirstpassTest, FirstpassInfoInitWithStaticBuf) {
  FIRSTPASS_INFO firstpass_info;
  aom_codec_err_t ret = av1_firstpass_info_init(&firstpass_info, NULL, 0);
  EXPECT_EQ(firstpass_info.stats_count, 0);
  EXPECT_EQ(firstpass_info.curr_index, 0);
  EXPECT_EQ(ret, AOM_CODEC_OK);
}

TEST(FirstpassTest, FirstpassInfoPushPop) {
  FIRSTPASS_INFO firstpass_info;
  av1_firstpass_info_init(&firstpass_info, NULL, 0);
  EXPECT_EQ(firstpass_info.stats_buf_size, ref_static_buf_size);
  for (int i = 0; i < ref_static_buf_size; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    stats.frame = i;
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, ref_static_buf_size);
  const int pop_count = ref_static_buf_size / 2;
  for (int i = 0; i < pop_count; ++i) {
    const FIRSTPASS_STATS *stats = av1_firstpass_info_peek(&firstpass_info, 0);
    aom_codec_err_t ret = av1_firstpass_info_move_curr_and_pop(&firstpass_info);
    EXPECT_NE(stats, nullptr);
    EXPECT_EQ(stats->frame, i);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, ref_static_buf_size - pop_count);

  const int push_count = ref_static_buf_size / 2;
  for (int i = 0; i < push_count; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, ref_static_buf_size);

  EXPECT_EQ(firstpass_info.stats_count, firstpass_info.stats_buf_size);
  // Pusht a stats when the queue is full.
  FIRSTPASS_STATS stats;
  av1_zero(stats);
  aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
  EXPECT_EQ(ret, AOM_CODEC_ERROR);
}

TEST(FirstpassTest, FirstpassInfoTotalStats) {
  FIRSTPASS_INFO firstpass_info;
  av1_firstpass_info_init(&firstpass_info, NULL, 0);
  EXPECT_EQ(firstpass_info.total_stats.frame, 0);
  for (int i = 0; i < 10; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    stats.count = 1;
    av1_firstpass_info_push(&firstpass_info, &stats);
  }
  EXPECT_EQ(firstpass_info.total_stats.count, 10);
}

TEST(FirstpassTest, FirstpassInfoMoveCurr) {
  FIRSTPASS_INFO firstpass_info;
  av1_firstpass_info_init(&firstpass_info, NULL, 0);
  int frame_cnt = 0;
  EXPECT_EQ(firstpass_info.stats_buf_size, ref_static_buf_size);
  for (int i = 0; i < ref_static_buf_size; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    stats.frame = frame_cnt;
    ++frame_cnt;
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.curr_index, firstpass_info.start_index);
  aom_codec_err_t ret = av1_firstpass_info_pop(&firstpass_info);
  // We cannot pop when curr_index == start_index
  EXPECT_EQ(ret, AOM_CODEC_ERROR);
  int ref_frame_cnt = 0;
  const int move_count = ref_static_buf_size * 2 / 3;
  for (int i = 0; i < move_count; ++i) {
    const FIRSTPASS_STATS *this_stats =
        av1_firstpass_info_peek(&firstpass_info, 0);
    EXPECT_EQ(this_stats->frame, ref_frame_cnt);
    ++ref_frame_cnt;
    av1_firstpass_info_move_curr(&firstpass_info);
  }
  EXPECT_EQ(firstpass_info.curr_stats_count_future,
            ref_static_buf_size - move_count);
  EXPECT_EQ(firstpass_info.curr_stats_count_past, move_count);
  EXPECT_EQ(firstpass_info.stats_count, ref_static_buf_size);

  const int test_count = ref_static_buf_size / 2;
  for (int i = 0; i < test_count; ++i) {
    aom_codec_err_t ret = av1_firstpass_info_pop(&firstpass_info);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }

  // Pop #test_count stats
  for (int i = 0; i < test_count; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    stats.frame = frame_cnt;
    ++frame_cnt;
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }

  // peek and move #test_count stats
  for (int i = 0; i < test_count; ++i) {
    const FIRSTPASS_STATS *this_stats =
        av1_firstpass_info_peek(&firstpass_info, 0);
    EXPECT_EQ(this_stats->frame, ref_frame_cnt);
    ++ref_frame_cnt;
    av1_firstpass_info_move_curr(&firstpass_info);
  }

  // pop #test_count stats
  for (int i = 0; i < test_count; ++i) {
    aom_codec_err_t ret = av1_firstpass_info_pop(&firstpass_info);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
}

}  // namespace
