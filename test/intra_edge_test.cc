/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/register_state_check.h"
#include "test/function_equivalence_test.h"

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "aom/aom_integer.h"

#include "./av1_rtcd.h"

#include "av1/common/enums.h"

using libaom_test::FunctionEquivalenceTest;

namespace {

template <typename F, typename T>
class UpsampleTest : public FunctionEquivalenceTest<F> {
 protected:
  static const int kIterations = 1e6;
  static const int kMinEdge = 4;
  static const int kMaxEdge = 24;
  static const int kBufSize = 2 * 64 + 32;
  static const int kOffset = 16;

  virtual ~UpsampleTest() {}

  virtual void Execute(T *edge_tst) = 0;

  void Common() {
    edge_ref_ = &edge_ref_data_[kOffset];
    edge_tst_ = &edge_tst_data_[kOffset];

    Execute(edge_tst_);

    int max_idx = (sz_ - 1) * 2;
    for (int r = -2; r <= max_idx; ++r) {
      ASSERT_EQ(edge_ref_[r], edge_tst_[r]);
    }
  }

  T edge_ref_data_[kBufSize];
  T edge_tst_data_[kBufSize];

  T *edge_ref_;
  T *edge_tst_;

  int sz_;
};

//////////////////////////////////////////////////////////////////////////////
// 8 bit version
//////////////////////////////////////////////////////////////////////////////

typedef void (*UP8B)(uint8_t *p, int sz);
typedef libaom_test::FuncParam<UP8B> TestFuncs;

class UpsampleTest8B : public UpsampleTest<UP8B, uint8_t> {
 protected:
  void Execute(uint8_t *edge_tst) {
    params_.ref_func(edge_ref_, sz_);
    ASM_REGISTER_STATE_CHECK(params_.tst_func(edge_tst, sz_));
  }
};

TEST_P(UpsampleTest8B, RandomValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    sz_ = 4 * (this->rng_(4) + 1);

    int i, pix = 0;
    for (i = 0; i < kOffset + sz_; ++i) {
      pix = rng_.Rand8();
      edge_ref_data_[i] = pix;
      edge_tst_data_[i] = edge_ref_data_[i];
    }

    // Extend final sample
    while (i < kBufSize) {
      edge_ref_data_[i] = pix;
      edge_tst_data_[i] = pix;
      i++;
    }

    Common();
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, UpsampleTest8B,
    ::testing::Values(TestFuncs(av1_upsample_intra_edge_c,
                                av1_upsample_intra_edge_sse4_1)));
#endif  // HAVE_SSE4_1

//////////////////////////////////////////////////////////////////////////////
// High bit-depth version
//////////////////////////////////////////////////////////////////////////////

typedef void (*UPHB)(uint16_t *p, int sz, int bd);
typedef libaom_test::FuncParam<UPHB> TestFuncsHBD;

class UpsampleTestHB : public UpsampleTest<UPHB, uint16_t> {
 protected:
  void Execute(uint16_t *edge_tst) {
    params_.ref_func(edge_ref_, sz_, bit_depth_);
    ASM_REGISTER_STATE_CHECK(params_.tst_func(edge_tst, sz_, bit_depth_));
  }
  int bit_depth_;
};

TEST_P(UpsampleTestHB, RandomValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    switch (rng_(3)) {
      case 0: bit_depth_ = 8; break;
      case 1: bit_depth_ = 10; break;
      default: bit_depth_ = 12; break;
    }
    const int hi = 1 << bit_depth_;

    sz_ = 4 * (this->rng_(4) + 1);

    int i, pix = 0;
    for (i = 0; i < kOffset + sz_; ++i) {
      pix = rng_(hi);
      edge_ref_data_[i] = pix;
      edge_tst_data_[i] = pix;
    }

    // Extend final sample
    while (i < kBufSize) {
      edge_ref_data_[i] = pix;
      edge_tst_data_[i] = pix;
      i++;
    }

    Common();
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, UpsampleTestHB,
    ::testing::Values(TestFuncsHBD(av1_upsample_intra_edge_high_c,
                                   av1_upsample_intra_edge_high_sse4_1)));
#endif  // HAVE_SSE4_1

}  // namespace
