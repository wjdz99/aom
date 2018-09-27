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
  // Generate images ranging in size from 0x0 to 6x6, with
  // varying levels of brightness. In all cases, the algorithm should
  // produce the same output.
  for (int width = 0; width < 7; ++width) {
    for (int height = 0; height < 7; ++height) {
      for (int brightness = 0; brightness < 255; ++brightness) {
        uint8_t *luma = (uint8_t *)malloc(width * height);
        for (int i = 0; i < width * height; ++i) {
          luma[i] = brightness;
        }
        uint8_t *output = (uint8_t *)malloc(width * height);
        gaussian_blur(luma, height, width, output);

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
  // Randomly generated 5x5.
  uint8_t luma[] = { 5,   3,   2,   100, 9,   15, 253, 16,  207,
                     23,  44,  129, 210, 248, 71, 61,  209, 240,
                     147, 192, 174, 14,  242, 0,  51 };
  uint8_t output[25];
  uint8_t expected[] = { 29,  51,  64,  66,  50,  56,  86,  104, 104,
                         83,  87,  114, 132, 129, 107, 112, 130, 141,
                         133, 115, 129, 130, 126, 112, 98 };
  gaussian_blur(luma, 5, 5, output);
  for (int i = 0; i < 25; ++i) {
    ASSERT_EQ(expected[i], output[i]);
  }
}

class SoebelTest : public ::testing::Test {};

TEST(SoebelTest, UniformBrightness) {
  // Generate images ranging in size from 0x0 to 4x4, with
  // varying levels of brightness. In all cases, the algorithm should
  // produce no magnitude for Soebel.
  for (int width = 0; width < 5; ++width) {
    for (int height = 0; height < 5; ++height) {
      for (int brightness = 0; brightness < 255; ++brightness) {
        uint8_t *luma = (uint8_t *)malloc(sizeof(uint8_t) * width * height);
        for (int i = 0; i < width * height; ++i) {
          luma[i] = brightness;
        }
        int *output = NULL;
        output = (int *)malloc(sizeof(*output) * width * height);
        soebel(luma, height, width, output);
        for (int i = 0; i < width * height; ++i) {
          ASSERT_EQ(0, output[i]);
        }
        free(luma);
        free(output);
      }
    }
  }
}

TEST(SoebelTest, SimpleExample) {
  // Randomly generated 3x3.
  uint8_t luma[] = { 5, 3, 2, 15, 253, 16, 44, 129, 210 };
  int output[9];
  int expected[] = { 363, 524, 377, 609, 525, 847, 494, 499, 458 };
  soebel(luma, 3, 3, output);
  for (int i = 0; i < 9; ++i) {
    ASSERT_EQ(expected[i], output[i]);
  }
}

class EdgeDetectTest : public ::testing::Test {};

TEST(EdgeDetectTest, UniformBrightness) {
  // Generate images ranging in size from 0x0 to 9x9, with
  // varying levels of brightness. In all cases, the algorithm should
  // detect absolutely no edge.
  for (int width = 0; width < 10; ++width) {
    for (int height = 0; height < 10; ++height) {
      for (int brightness = 0; brightness < 255; ++brightness) {
        uint8_t *luma = (uint8_t *)malloc(sizeof(uint8_t) * width * height);
        for (int i = 0; i < width * height; ++i) {
          luma[i] = brightness;
        }
        ASSERT_FLOAT_EQ(0, av1_edge_exists(luma, height, width));
        free(luma);
      }
    }
  }
}

// Generate images ranging in size from 2x2 to 9x9, black on one side
// and white on the other.
TEST(EdgeDetectTest, BlackWhite) {
  for (int width = 2; width < 10; ++width) {
    for (int height = 2; height < 10; ++height) {
      uint8_t *luma = (uint8_t *)malloc(sizeof(uint8_t) * width * height);

      for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
          if (i < width / 2) {
            luma[i + j * width] = 0;
          } else {
            luma[i + j * width] = 255;
          }
        }
      }
      ASSERT_FLOAT_EQ(1, av1_edge_exists(luma, height, width));
      free(luma);
    }
  }
}

}  // namespace
