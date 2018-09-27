/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/encoder/rdopt.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

class GaussianBlurTest : public ::testing::Test {};

TEST(GaussianBlurTest, UniformBrightness) {
  // Generate images ranging in size from 3x3 to 6x6, with
  // varying levels of brightness. In all cases, the algorithm should
  // produce the same output.
  for (int width = 3; width < 7; ++width) {
    for (int height = 3; height < 7; ++height) {
      for (int brightness = 0; brightness < 255; ++brightness) {
        // Allocate an extra 4 pixel border.
        int stride = width + 8;
        uint8_t *luma = (uint8_t *)malloc(stride * (8 + height));
        for (int i = 0; i < stride * (8 + height); ++i) {
          luma[i] = brightness;
        }
        uint8_t *output = (uint8_t *)malloc(width * height);
        uint8_t *offset = luma + stride * 4 + 4;
        gaussian_blur(offset, stride, width, height, output);
        for (int i = 0; i < width * height; ++i) {
          ASSERT_EQ(brightness, output[i]);
        }
        free(luma);
        free(output);
      }
    }
  }
}

TEST(GaussianBlurTest, SimpleExample) {
  // Randomly generated 5x5 with repeated 4-pixel border.
  uint8_t luma[] = {
    5,   5,   5,   5,   5,  3,   2,   100, 9,   9,   9,   9,   9,   5,   5,
    5,   5,   5,   3,   2,  100, 9,   9,   9,   9,   9,   5,   5,   5,   5,
    5,   3,   2,   100, 9,  9,   9,   9,   9,   5,   5,   5,   5,   5,   3,
    2,   100, 9,   9,   9,  9,   9,   5,   5,   5,   5,   5,   3,   2,   100,
    9,   9,   9,   9,   9,  15,  15,  15,  15,  15,  253, 16,  207, 23,  23,
    23,  23,  23,  44,  44, 44,  44,  44,  129, 210, 248, 71,  71,  71,  71,
    71,  61,  61,  61,  61, 61,  209, 240, 147, 192, 192, 192, 192, 192, 174,
    174, 174, 174, 174, 14, 242, 0,   51,  51,  51,  51,  51,  174, 174, 174,
    174, 174, 14,  242, 0,  51,  51,  51,  51,  51,  174, 174, 174, 174, 174,
    14,  242, 0,   51,  51, 51,  51,  51,  174, 174, 174, 174, 174, 14,  242,
    0,   51,  51,  51,  51, 51,  174, 174, 174, 174, 174, 14,  242, 0,   51,
    51,  51,  51,  51
  };
  uint8_t output[25];
  uint8_t expected[] = { 31,  48,  62,  65,  51,  58,  85,  101, 100,
                         79,  86,  116, 134, 129, 107, 110, 130, 143,
                         132, 113, 127, 129, 129, 111, 93 };
  const int width = 5;
  const int height = 5;
  const int stride = width + 8;
  gaussian_blur(luma + stride * 4 + 4, stride, width, height, output);
  for (int i = 0; i < 25; ++i) {
    ASSERT_EQ(expected[i], output[i]);
  }
}

class EdgeDetectTest : public ::testing::Test {};

TEST(EdgeDetectTest, UniformBrightness) {
  // Generate images ranging in size from 0x0 to 9x9, with
  // varying levels of brightness. In all cases, the algorithm should
  // produce the same output.
  for (int width = 0; width < 10; ++width) {
    for (int height = 0; height < 10; ++height) {
      for (int brightness = 0; brightness < 255; ++brightness) {
        // Allocate an extra 4 pixel border.
        int stride = width + 8;
        uint8_t *luma = (uint8_t *)malloc(stride * (8 + height));
        for (int i = 0; i < stride * (8 + height); ++i) {
          luma[i] = brightness;
        }
        uint8_t *offset = luma + stride * 4 + 4;
        ASSERT_EQ(0, av1_edge_exists(offset, stride, width, height));
        free(luma);
      }
    }
  }
}

// Generate images ranging in size from 0x0 to 9x9, black on one side
// and white on the other. All of them have a 4 pixel border.
TEST(EdgeDetectTest, BlackWhite) {
  for (int width = 0; width < 10; ++width) {
    for (int height = 0; height < 10; ++height) {
      int stride = width + 8;
      int padded_h = height + 8;
      uint8_t *luma = (uint8_t *)malloc(sizeof(uint8_t) * stride * padded_h);
      for (int i = 0; i < stride; ++i) {
        for (int j = 0; j < padded_h; ++j) {
          if (i < stride / 2) {
            luma[i + j * stride] = 0;
          } else {
            luma[i + j * stride] = 255;
          }
        }
      }
      uint8_t *src = luma + stride * 4 + 4;
      if (height < 3 || width < 3) {
        ASSERT_LE(0, av1_edge_exists(src, stride, width, height));
      } else {
        ASSERT_LE(548, av1_edge_exists(src, stride, width, height));
      }
      free(luma);
    }
  }
}

}  // namespace
