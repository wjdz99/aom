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

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "av1/encoder/tpl_model.h"

TEST(TplModelTest, TransformCoeffEntropyTest1) {
  // Check the consistency between av1_est_coeff_entropy() and
  // av1_laplace_prob()
  double b = 1;
  double q_step = 1;
  double zero_bin_ratio = 2;
  for (int qcoeff = -10; qcoeff < 10; ++qcoeff) {
    double rate = av1_est_coeff_entropy(q_step, b, zero_bin_ratio, qcoeff);
    double prob = av1_laplace_prob(q_step, b, zero_bin_ratio, qcoeff);
    double ref_rate = -log2(prob);
    EXPECT_DOUBLE_EQ(rate, ref_rate);
  }
}

TEST(TplModelTest, TransformCoeffEntropyTest2) {
  // Check the consistency between av1_est_coeff_entropy(), av1_laplace_prob()
  // and
  double b = 1;
  double q_step = 1;
  double zero_bin_ratio = 2;
  double est_expected_rate = 0;
  for (int qcoeff = -20; qcoeff < 20; ++qcoeff) {
    double rate = av1_est_coeff_entropy(q_step, b, zero_bin_ratio, qcoeff);
    double prob = av1_laplace_prob(q_step, b, zero_bin_ratio, qcoeff);
    est_expected_rate += prob * rate;
  }
  double expected_rate = av1_laplace_entropy(q_step, b, zero_bin_ratio);
  EXPECT_NEAR(expected_rate, est_expected_rate, 0.001);
}

TEST(TplModelTest, av1_delta_rate_cost1) {
  // When srcrf_dist equal to recrf_dist, av1_delta_rate_cost should return 0
  int64_t srcrf_dist = 256;
  int64_t recrf_dist = 256;
  int64_t delta_rate = 512;
  int pixel_num = 256;
  int64_t rate_cost =
      av1_delta_rate_cost(delta_rate, recrf_dist, srcrf_dist, pixel_num);
  EXPECT_EQ(rate_cost, 0);
}
