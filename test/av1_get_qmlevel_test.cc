/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <stdlib.h>

#include "gtest/gtest.h"

#include "av1/av1_cx_iface.h"

namespace {

TEST(GetQmLevelTest, Regular) {
  // Get some QM levels from representative qindex values.
  const int qindex_1_qmlevel =
      aom_get_qmlevel(1, DEFAULT_QM_FIRST, DEFAULT_QM_LAST);
  const int qindex_60_qmlevel =
      aom_get_qmlevel(60, DEFAULT_QM_FIRST, DEFAULT_QM_LAST);
  const int qindex_120_qmlevel =
      aom_get_qmlevel(120, DEFAULT_QM_FIRST, DEFAULT_QM_LAST);
  const int qindex_180_qmlevel =
      aom_get_qmlevel(180, DEFAULT_QM_FIRST, DEFAULT_QM_LAST);
  const int qindex_255_qmlevel =
      aom_get_qmlevel(255, DEFAULT_QM_FIRST, DEFAULT_QM_LAST);

  // Extreme qindex values should also result in extreme QM levels.
  ASSERT_EQ(qindex_1_qmlevel, DEFAULT_QM_FIRST);
  ASSERT_EQ(qindex_255_qmlevel, DEFAULT_QM_LAST);

  // aom_get_qmlevel() QMs become steeper (i.e. levels become lower) the lower
  // the qindex.
  ASSERT_GE(qindex_255_qmlevel, qindex_180_qmlevel);
  ASSERT_GE(qindex_180_qmlevel, qindex_120_qmlevel);
  ASSERT_GE(qindex_120_qmlevel, qindex_60_qmlevel);
  ASSERT_GE(qindex_60_qmlevel, qindex_1_qmlevel);
  ASSERT_GT(qindex_255_qmlevel, qindex_1_qmlevel);

  // Set min and max QM levels to be lower than DEFAULT_QM_FIRST.
  const int qindex_1_qmlevel_belowfirst =
      aom_get_qmlevel_allintra(1, 1, DEFAULT_QM_FIRST - 1);
  const int qindex_255_qmlevel_belowfirst =
      aom_get_qmlevel_allintra(255, 1, DEFAULT_QM_FIRST - 1);

  // Formula should always respect QM level boundaries, even when they're below
  // DEFAULT_QM_FIRST.
  ASSERT_LT(qindex_1_qmlevel_belowfirst, DEFAULT_QM_FIRST);
  ASSERT_LT(qindex_255_qmlevel_belowfirst, DEFAULT_QM_FIRST);

  // Set min and max QM levels to be higher than DEFAULT_QM_LAST.
  const int qindex_1_qmlevel_abovelast =
      aom_get_qmlevel_allintra(1, DEFAULT_QM_LAST + 1, 15);
  const int qindex_255_qmlevel_abovelast =
      aom_get_qmlevel_allintra(255, DEFAULT_QM_LAST + 1, 15);

  // Formula should always respect QM level boundaries, even when they're above
  // DEFAULT_QM_LAST.
  ASSERT_GT(qindex_1_qmlevel_abovelast, DEFAULT_QM_LAST);
  ASSERT_GT(qindex_255_qmlevel_abovelast, DEFAULT_QM_LAST);
}

TEST(GetQmLevelTest, AllIntra) {
  // Get some QM levels from representative qindex values.
  const int qindex_1_qmlevel = aom_get_qmlevel_allintra(
      1, DEFAULT_QM_FIRST_ALLINTRA, DEFAULT_QM_LAST_ALLINTRA);
  const int qindex_60_qmlevel = aom_get_qmlevel_allintra(
      60, DEFAULT_QM_FIRST_ALLINTRA, DEFAULT_QM_LAST_ALLINTRA);
  const int qindex_120_qmlevel = aom_get_qmlevel_allintra(
      120, DEFAULT_QM_FIRST_ALLINTRA, DEFAULT_QM_LAST_ALLINTRA);
  const int qindex_180_qmlevel = aom_get_qmlevel_allintra(
      180, DEFAULT_QM_FIRST_ALLINTRA, DEFAULT_QM_LAST_ALLINTRA);
  const int qindex_255_qmlevel = aom_get_qmlevel_allintra(
      255, DEFAULT_QM_FIRST_ALLINTRA, DEFAULT_QM_LAST_ALLINTRA);

  // Extreme qindex values should also result in extreme QM levels.
  ASSERT_EQ(qindex_1_qmlevel, DEFAULT_QM_LAST_ALLINTRA);
  ASSERT_EQ(qindex_255_qmlevel, DEFAULT_QM_FIRST_ALLINTRA);

  // Unlike with aom_get_qmlevel(), aom_get_qmlevel_allintra() QMs become
  // flatter (i.e. levels become higher) the lower the qindex.
  ASSERT_LE(qindex_255_qmlevel, qindex_180_qmlevel);
  ASSERT_LE(qindex_180_qmlevel, qindex_120_qmlevel);
  ASSERT_LE(qindex_120_qmlevel, qindex_60_qmlevel);
  ASSERT_LE(qindex_60_qmlevel, qindex_1_qmlevel);
  ASSERT_LT(qindex_255_qmlevel, qindex_1_qmlevel);

  // Set min and max QM levels to be lower than DEFAULT_QM_FIRST_ALLINTRA.
  const int qindex_1_qmlevel_belowfirst =
      aom_get_qmlevel_allintra(1, 1, DEFAULT_QM_FIRST_ALLINTRA - 1);
  const int qindex_255_qmlevel_belowfirst =
      aom_get_qmlevel_allintra(255, 1, DEFAULT_QM_FIRST_ALLINTRA - 1);

  // Formula should always respect QM level boundaries, even when they're below
  // DEFAULT_QM_FIRST_ALLINTRA.
  ASSERT_EQ(qindex_1_qmlevel_belowfirst, DEFAULT_QM_FIRST_ALLINTRA - 1);
  ASSERT_EQ(qindex_255_qmlevel_belowfirst, DEFAULT_QM_FIRST_ALLINTRA - 1);

  // Set min and max QM levels to be higher than DEFAULT_QM_LAST_ALLINTRA.
  const int qindex_1_qmlevel_abovelast =
      aom_get_qmlevel_allintra(1, DEFAULT_QM_LAST_ALLINTRA + 1, 15);
  const int qindex_255_qmlevel_abovelast =
      aom_get_qmlevel_allintra(255, DEFAULT_QM_LAST_ALLINTRA + 1, 15);

  // Formula should always respect QM level boundaries, even when they're above
  // DEFAULT_QM_LAST_ALLINTRA.
  ASSERT_EQ(qindex_1_qmlevel_abovelast, DEFAULT_QM_LAST_ALLINTRA + 1);
  ASSERT_EQ(qindex_255_qmlevel_abovelast, DEFAULT_QM_LAST_ALLINTRA + 1);
}

}  // namespace
