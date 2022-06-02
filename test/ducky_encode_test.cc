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

#include "av1/ducky_encode.h"
#include "av1/encoder/encoder.h"
#include "test/video_source.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#define AV1_FOURCC 0x31305641

namespace aom {

void Lay16Bits(uint8_t *array, int index, int value) {
  array[index] = (uint8_t)(value);
  array[index + 1] = (uint8_t)(value >> 8);
}

void Lay32Bits(uint8_t *array, int index, int value) {
  for (int i = 0; i < 4; i++) {
    array[index + i] = (uint8_t)(value >> (i * 8));
  }
}

void Lay64Bits(uint8_t *array, int index, int64_t value) {
  for (int i = 0; i < 8; i++) {
    array[index + i] = (uint8_t)(value >> (i * 8));
  }
}

TEST(DuckyEncodeTest, ComputeFirstPassStats) {
  aom_rational_t frame_rate = { 30, 1 };
  VideoInfo video_info = { 352,        288,
                           frame_rate, AOM_IMG_FMT_I420,
                           1,          "bus_352x288_420_f20_b8.yuv" };
  video_info.file_path =
      libaom_test::GetDataPath() + "/" + video_info.file_path;
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
  video_info.file_path =
      libaom_test::GetDataPath() + "/" + video_info.file_path;
  DuckyEncode ducky_encode(video_info);
  std::vector<FIRSTPASS_STATS> frame_stats =
      ducky_encode.ComputeFirstPassStats();
  ducky_encode.StartEncode(frame_stats);
  // We set coding_frame_count to a arbitrary number that smaller than
  // 17 here.
  // TODO(angiebird): Set coding_frame_count properly, once the DuckyEncode can
  // provide proper information.
  int coding_frame_count = 5;
  EncodeFrameDecision decision = { EncodeFrameMode::kNone, {} };

  FILE *bitstream_file = fopen(
      "/usr/local/google/home/jianj/Codes/aom/linux_build/bitstream.ivf", "wb");
  assert(bitstream_file != nullptr);
  uint8_t ivf_header_[32];
  memset(ivf_header_, 0, 32);

  ivf_header_[0] = 'D';
  ivf_header_[1] = 'K';
  ivf_header_[2] = 'I';
  ivf_header_[3] = 'F';
  Lay16Bits(ivf_header_, 4, 0);   // version
  Lay16Bits(ivf_header_, 6, 32);  // header size
  ivf_header_[8] = 'A';           // fourcc
  ivf_header_[9] = 'V';
  ivf_header_[10] = '0';
  ivf_header_[11] = '1';

  Lay16Bits(ivf_header_, 12, 352);
  Lay16Bits(ivf_header_, 14, 288);
  Lay32Bits(ivf_header_, 16, 0);
  Lay32Bits(ivf_header_, 20, 0);
  Lay32Bits(ivf_header_, 24, 10);
  Lay32Bits(ivf_header_, 28, 0);

  for (int i = 0; i < coding_frame_count; ++i) {
    EncodeFrameResult encode_frame_result = ducky_encode.EncodeFrame(decision);
    uint8_t frame_header[12];
    size_t frame_size = encode_frame_result.bitstream_buf.size();
    Lay32Bits(frame_header, 0, (int)frame_size);
    Lay64Bits(frame_header, 4, 0);

    if (i == 0) fwrite(ivf_header_, 32, 1, bitstream_file);

    fwrite(frame_header, 12, 1, bitstream_file);
    fwrite(encode_frame_result.bitstream_buf.data(), 1, frame_size,
           bitstream_file);
  }
  ducky_encode.EndEncode();
  fclose(bitstream_file);
}

// TODO(b/231517281): Re-enable after fix.
TEST(DuckyEncodeTest, DISABLED_EncodeFrameWithQindex) {
  aom_rational_t frame_rate = { 30, 1 };
  VideoInfo video_info = { 352,        288,
                           frame_rate, AOM_IMG_FMT_I420,
                           17,         "bus_352x288_420_f20_b8.yuv" };
  video_info.file_path =
      libaom_test::GetDataPath() + "/" + video_info.file_path;
  DuckyEncode ducky_encode(video_info);
  std::vector<FIRSTPASS_STATS> frame_stats =
      ducky_encode.ComputeFirstPassStats();
  ducky_encode.StartEncode(frame_stats);
  // We set coding_frame_count to a arbitrary number that smaller than
  // 17 here.
  // TODO(angiebird): Set coding_frame_count properly, once the DuckyEncode can
  // provide proper information.
  int coding_frame_count = 5;
  int q_index = 0;
  EncodeFrameDecision decision = { EncodeFrameMode::kQindex, { q_index, -1 } };
  for (int i = 0; i < coding_frame_count; ++i) {
    EncodeFrameResult encode_frame_result = ducky_encode.EncodeFrame(decision);
    // TODO(angiebird): Check why distortion is not zero when q_index = 0
    EXPECT_EQ(encode_frame_result.dist, 0);
  }
  ducky_encode.EndEncode();
}

TEST(DuckyEncodeTest, EncodeFrameMode) {
  EXPECT_EQ(DUCKY_ENCODE_FRAME_MODE_NONE,
            static_cast<DUCKY_ENCODE_FRAME_MODE>(EncodeFrameMode::kNone));
  EXPECT_EQ(DUCKY_ENCODE_FRAME_MODE_QINDEX,
            static_cast<DUCKY_ENCODE_FRAME_MODE>(EncodeFrameMode::kQindex));
  EXPECT_EQ(
      DUCKY_ENCODE_FRAME_MODE_QINDEX_RDMULT,
      static_cast<DUCKY_ENCODE_FRAME_MODE>(EncodeFrameMode::kQindexRdmult));
}

}  // namespace aom
