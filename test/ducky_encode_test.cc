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

#include <array>
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "common/ivfenc.h"
#include "av1/ducky_encode.h"
#include "av1/encoder/encoder.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace aom {
TEST(DuckyEncodeTest, ComputeFirstPassStats) {
  aom_rational_t frame_rate = { 30, 1 };
  VideoInfo video_info = { 352,        288,
                           frame_rate, AOM_IMG_FMT_I420,
                           1,          "bus_352x288_420_f20_b8.yuv" };
  DuckyEncode ducky_encode(video_info);
  std::vector<FIRSTPASS_STATS> frame_stats =
      ducky_encode.ComputeFirstPassStats();
  EXPECT_EQ(frame_stats.size(), static_cast<size_t>(video_info.frame_count));
  for (size_t i = 0; i < frame_stats.size(); ++i) {
    // FIRSTPASS_STATS's first element is frame
    EXPECT_EQ(frame_stats[i].frame, i);
  }
}

TEST(DuckyEncodeTest, EncodeFrame) {
  aom_rational_t frame_rate = { 30, 1 };
  VideoInfo video_info = { 352,        288,
                           frame_rate, AOM_IMG_FMT_I420,
                           17,         "bus_352x288_420_f20_b8.yuv" };
  FILE *outfile = fopen("output.ivf", "w");
  DuckyEncode ducky_encode(video_info);
  std::vector<FIRSTPASS_STATS> frame_stats =
      ducky_encode.ComputeFirstPassStats();
  int coding_idx = 0;
  aom_rational_t time_base = { 1, TICKS_PER_SEC };
  uint32_t fourcc = AV1_FOURCC;
  ivf_write_file_header_with_video_info(outfile, fourcc, video_info.frame_count,
                                        video_info.frame_width,
                                        video_info.frame_height, time_base);
  ducky_encode.StartEncode(frame_stats);
  for (int i = 0; i < 1; ++i) {
    EncodeFrameResult encode_frame_result = ducky_encode.EncodeFrame();
    ivf_write_frame_header(outfile, coding_idx,
                           encode_frame_result.bitstream_buf.size());
    fwrite(encode_frame_result.bitstream_buf.data(), 1,
           encode_frame_result.bitstream_buf.size(), outfile);
  }
  fclose(outfile);
  // encode_frame_result.bitstream_buf
  ducky_encode.EndEncode();
}

}  // namespace aom
