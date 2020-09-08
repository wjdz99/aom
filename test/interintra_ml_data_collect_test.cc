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

TEST(InterIntraMLDataCollectTest, ComputeLShape) {
  const int border = 4;
  const int width = 4;
  const int height = 4;
  const int src_stride = 8;
  uint8_t dst[8 * 4 + 4 * 4];
  uint8_t src[256];
  for (int i = 0; i < 256; ++i) {
    src[i] = static_cast<uint8_t>(i);
  }

  const int bitdepth = 8;
  const bool is_src_hbd = false;
  av1_interintra_ml_data_collect_copy_intrapred_lshape(
      dst, src + 4 + 20 * 4, src_stride, width, height, border, bitdepth,
      is_src_hbd);
  /* clang-format off */
  uint8_t expected[48] = { 0,  1,  2,  3,  4,  5,  6,  7,
                           8,  9,  10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23,
                           24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35,
                           40, 41, 42, 43,
                           48, 49, 50, 51,
                           56, 57, 58, 59 };
  /* clang-format on */
  ASSERT_EQ(0, memcmp(expected, dst, 48));
}

}  // namespace

#endif  // CONFIG_INTERINTRA_ML_DATA_COLLECT
