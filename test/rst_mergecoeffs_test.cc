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

#include <math.h>
#include "av1/encoder/pickrst.h"
#include "third_party/vector/vector.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

TEST(AV1RstMergeCoeffsTest, TestUntypedGraphSearch) {
  // Initialize graph.
  double graph[25];
  for (int i = 0; i < 25; ++i) {
    graph[i] = INFINITY;
  }
  graph[1] = -1;
  graph[7] = -2;
  graph[10] = -3;
  graph[15] = -1;
  graph[17] = -1;
  graph[23] = 2;
  // Initialize best_path vector.
  Vector best_path;
  aom_vector_setup(&best_path, 1, sizeof(int));

  double cost = min_cost_type_path(4, 1, 5, graph, &best_path, false);
  int correct_path[] = { 4, 3, 2, 0, 1 };
  int i = 0;
  VECTOR_FOR_EACH(&best_path, listed_unit) {
    int node = *(int *)(listed_unit.pointer);
    EXPECT_TRUE(correct_path[i] == node);
    ++i;
  }
  EXPECT_TRUE(cost == -3);
  aom_vector_destroy(&best_path);
}

TEST(AV1RstMergeCoeffsTest, TestTypedGraphSearch) {
  // Initialize graph.
  double graph[] = {
    8, 8,        8,        8,        6,        7,        8, 6,        7,
    8, 6,        7,        1,        6,        8,        8, 6,        8,
    8, 6,        8,        1,        8,        8,        8, 8,        8,
    8, 8,        1,        8,        7,        8,        8, 7,        8,
    8, 7,        1,        8,        6,        8,        8, 6,        8,
    8, 6,        1,        0,        INFINITY, INFINITY, 0, INFINITY, INFINITY,
    0, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY
  };
  // Initialize best_path vector.
  Vector best_path;
  aom_vector_setup(&best_path, 1, sizeof(int));

  double cost = min_cost_type_path(0, 19, 3, graph, &best_path, true);
  int correct_path[] = { 0, 1, 5, 9, 12, 15, 18, 19 };
  int i = 0;
  VECTOR_FOR_EACH(&best_path, listed_unit) {
    int node = *(int *)(listed_unit.pointer);
    EXPECT_TRUE(correct_path[i] == node);
    ++i;
  }
  EXPECT_TRUE(cost == 25);
  aom_vector_destroy(&best_path);
}

}  // namespace
