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

#include "av1/common/restoration.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

TEST(AV1RstMergeCoeffTest, CheckWienerEq) {
  WienerInfo w1 = { .vfilter = { 1, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 1, 2, 3, 4, 5, 6, 7, 8 } };
  WienerInfo w2 = { .vfilter = { 1, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 1, 2, 3, 4, 5, 6, 7, 8 } };
  EXPECT_TRUE(check_wiener_eq(&w1, &w2));
  EXPECT_TRUE(check_wiener_eq(&w2, &w1));
}

TEST(AV1RstMergeCoeffTest, CheckWienerNeqHFilter) {
  WienerInfo w1 = { .vfilter = { 1, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 1, 2, 3, 4, 5, 6, 7, 8 } };
  WienerInfo w2 = { .vfilter = { 1, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 2, 2, 3, 4, 5, 6, 7, 8 } };
  EXPECT_FALSE(check_wiener_eq(&w1, &w2));
  EXPECT_FALSE(check_wiener_eq(&w2, &w1));
}

TEST(AV1RstMergeCoeffTest, CheckWienerNeqVFilter) {
  WienerInfo w1 = { .vfilter = { 1, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 1, 2, 3, 4, 5, 6, 7, 8 } };
  WienerInfo w2 = { .vfilter = { 2, 2, 3, 4, 5, 6, 7, 8 },
                    .hfilter = { 1, 2, 3, 4, 5, 6, 7, 8 } };
  EXPECT_FALSE(check_wiener_eq(&w1, &w2));
  EXPECT_FALSE(check_wiener_eq(&w2, &w1));
}

}  // namespace
