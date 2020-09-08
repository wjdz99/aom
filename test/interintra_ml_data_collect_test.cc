/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <cstring>
#include <stdbool.h>
#include "test/util.h"
#include "av1/encoder/interintra_ml_data_collect.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#if CONFIG_INTERINTRA_ML_DATA_COLLECT

namespace {

// Testing infrastructure for the L-shape copy function. In all test
// cases, the source buffer has 256 elements with src_[i] == i.
// If the high-bitdepth version of it is requested, then src16_[i] == (i << 4).
// A high-bitdepth but 8-bit version can be requested with LowbdSrc16,
// which encodes lowbd_src16_[i] == i.
class InterIntraMLDataCollectTest : public ::testing::Test {
 protected:
  void SetUp() override {
    for (int i = 0; i < 256; ++i) {
      src_[i] = i;
      lowbd_src16_[i] = i;
      src16_[i] = i << 4;
    }
  }

  // Get a reference to the buffer, offset so the border is grabbed at a
  // negative offset.
  const uint8_t *Src(int border) const {
    return src_ + border + SrcStride() * border;
  }

  const uint8_t *Src16(int border) const {
    return CONVERT_TO_BYTEPTR(src16_ + border + SrcStride() * border);
  }

  const uint8_t *LowbdSrc16(int border) const {
    return CONVERT_TO_BYTEPTR(lowbd_src16_ + border + SrcStride() * border);
  }

  int SrcStride() const { return 16; }

 private:
  uint8_t src_[256];
  uint16_t src16_[256];
  uint16_t lowbd_src16_[256];
};

// 4x8 block with a 1 pixel border.
TEST_F(InterIntraMLDataCollectTest, ComputeLShape4x8) {
  const int border = 1;
  const int width = 4;
  const int height = 8;
  uint8_t dst[border + width + height];

  // Low-bitdepth test.
  int bitdepth = 8;
  bool is_src_hbd = false;
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst, Src(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(SrcStride(), 16);
  /* clang-format off */
  uint8_t expected[13] = { 0,  1,  2,  3,  4,
                           16,
                           32,
                           48,
                           64,
                           80,
                           96,
                           112,
                           128 };
  /* clang-format on */
  ASSERT_EQ(0, memcmp(expected, dst, 13));

  // Same results expected if the pipeline is high bit-depth but the
  // actual data is 8-bit.
  memset(dst, 0, 13);
  ASSERT_FALSE(0 == memcmp(expected, dst, 13));
  is_src_hbd = true;
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst, LowbdSrc16(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(0, memcmp(expected, dst, 13));

  // High-bitdepth test.
  uint16_t expected16[13];
  for (int i = 0; i < 13; ++i) {
    expected16[i] = expected[i] << 4;
  }
  bitdepth = 12;
  is_src_hbd = true;
  uint8_t dst16[26];
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst16, Src16(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(0, memcmp(expected16, dst16, 26));
}

// 4x4 block with a 4 pixel border, low bitdepth.
TEST_F(InterIntraMLDataCollectTest, ComputeLShape4x4) {
  const int border = 4;
  const int width = 4;
  const int height = 4;
  uint8_t dst[(border + width) * border + border * height];
  int bitdepth = 8;
  bool is_src_hbd = false;
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst, Src(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(SrcStride(), 16);
  /* clang-format off */
  uint8_t expected[48] = { 0,  1,  2,  3,  4,  5,  6,  7,
                           16, 17, 18, 19, 20, 21, 22, 23,
                           32, 33, 34, 35, 36, 37, 38, 39,
                           48, 49, 50, 51, 52, 53, 54, 55,
                           64, 65, 66, 67,
                           80, 81, 82, 83,
                           96, 97, 98, 99,
                           112, 113, 114, 115 };
  /* clang-format on */
  ASSERT_EQ(0, memcmp(expected, dst, 48));

  // Same results expected if the pipeline is high bit-depth but the
  // actual data is 8-bit.
  memset(dst, 0, 48);
  ASSERT_FALSE(0 == memcmp(expected, dst, 48));
  is_src_hbd = true;
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst, LowbdSrc16(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(0, memcmp(expected, dst, 48));

  // High-bitdepth test.
  uint16_t expected16[48];
  for (int i = 0; i < 48; ++i) {
    expected16[i] = expected[i] << 4;
  }
  bitdepth = 12;
  is_src_hbd = true;
  uint8_t dst16[96];
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst16, Src16(border), SrcStride(), width, height, border, bitdepth,
      is_src_hbd);
  ASSERT_EQ(0, memcmp(expected16, dst16, 96));
}

TEST_F(InterIntraMLDataCollectTest, SerializeLowbd) {
  /* clang-format off */
  uint8_t source_image[] = { 16, 17, 18, 19,
                             20, 21, 22, 23,
                             24, 25, 26, 27,
                             28, 29, 30, 31 };
  uint8_t intrapred_lshape[] = { 91, 92, 93, 94, 95,
                                 96,
                                 97,
                                 98,
                                 99 };
  uint8_t interpred[] = { 35, 36, 37, 38, 39,
                          40, 41, 42, 43, 44,
                          45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54,
                          55, 56, 57, 58, 59 };
  uint8_t expected[73] = { 4 /* width */,
                           4 /* height */,
                           2 /* plane */,
                           8 /* bitdepth */,
                           1 /* border */,
                           0, 0, 2, 43, /* x */
                           0, 0, 0, 26, /* y */
                           0, 0, 0, 7, /* frame_order_hint */
                           0, 0, 48, 57, /* lambda */
                           25 /* ref_q */,
                           13 /* base_q */,
                           /* source_image */
                           16, 17, 18, 19,
                           20, 21, 22, 23,
                           24, 25, 26, 27,
                           28, 29, 30, 31,
                           /* intrapred_lshape */
                           91, 92, 93, 94, 95,
                           96,
                           97,
                           98,
                           99,
                           /* interpred */
                           35, 36, 37, 38, 39,
                           40, 41, 42, 43, 44,
                           45, 46, 47, 48, 49,
                           50, 51, 52, 53, 54,
                           55, 56, 57, 58, 59 };
  /* clang-format on */
  IIMLPlaneInfo info = { .width = 4,
                         .height = 4,
                         .plane = 2,
                         .bitdepth = 8,
                         .border = 1,
                         .x = 555,
                         .y = 26,
                         .frame_order_hint = 7,
                         .lambda = 12345,
                         .ref_q = 25,
                         .base_q = 13,
                         .source_image = source_image,
                         .intrapred_lshape = intrapred_lshape,
                         .interpred = interpred };
  uint8_t *buf;
  size_t buf_size;
  av1_interintra_ml_data_collect_serialize(&info, &buf, &buf_size);
  ASSERT_EQ(73U, buf_size);
  ASSERT_EQ(0, memcmp(expected, buf, buf_size));
  free(buf);
}

TEST_F(InterIntraMLDataCollectTest, SerializeHighbd) {
  /* clang-format off */
  uint16_t source_image[] = { 160, 170,
                              200, 210,
                              240, 250,
                              280, 290 };
  uint16_t intrapred_lshape[] = { 91, 92, 930,
                                  960,
                                  970,
                                  980,
                                  99 };
  uint16_t interpred[] = { 35, 36, 37,
                          40, 41, 42,
                          45, 46, 470,
                          50, 51, 520,
                           55, 56, 57 };
  uint8_t expected[83] = { 2 /* width */,
                           4 /* height */,
                           0 /* plane */,
                           10 /* bitdepth */,
                           1 /* border */,
                           0, 0, 2, 43, /* x */
                           0, 0, 0, 26, /* y */
                           0, 0, 0, 7, /* frame_order_hint */
                           0, 0, 48, 57, /* lambda */
                           25 /* ref_q */,
                           13 /* base_q */,
                           /* source_image */
                           0, 160, 0, 170,
                           0, 200, 0, 210,
                           0, 240, 0, 250,
                           1, 24,  1, 34,
                           /* intrapred_lshape */
                           0, 91, 0, 92, 3, 162,
                           3, 192,
                           3, 202,
                           3, 212,
                           0, 99,
                           /* interpred */
                           0, 35, 0, 36, 0, 37,
                           0, 40, 0, 41, 0, 42,
                           0, 45, 0, 46, 1, 214,
                           0, 50, 0, 51, 2, 8,
                           0, 55, 0, 56, 0, 57 };
  /* clang-format on */
  IIMLPlaneInfo info = { .width = 2,
                         .height = 4,
                         .plane = 0,
                         .bitdepth = 10,
                         .border = 1,
                         .x = 555,
                         .y = 26,
                         .frame_order_hint = 7,
                         .lambda = 12345,
                         .ref_q = 25,
                         .base_q = 13,
                         .source_image =
                             reinterpret_cast<uint8_t *>(source_image),
                         .intrapred_lshape =
                             reinterpret_cast<uint8_t *>(intrapred_lshape),
                         .interpred = reinterpret_cast<uint8_t *>(interpred) };
  uint8_t *buf;
  size_t buf_size;
  av1_interintra_ml_data_collect_serialize(&info, &buf, &buf_size);
  ASSERT_EQ(83U, buf_size);
  ASSERT_EQ(0, memcmp(expected, buf, buf_size));
  free(buf);
}

}  // namespace

#endif  // CONFIG_INTERINTRA_ML_DATA_COLLECT
