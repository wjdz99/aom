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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "test/function_equivalence_test.h"
#include "test/register_state_check.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "aom/aom_integer.h"
#include "aom_ports/aom_timer.h"

#include "./av1_rtcd.h"

#include "av1/common/enums.h"

using libaom_test::FunctionEquivalenceTest;

namespace {
static const TX_SIZE kTxSize[] = {
    TX_4X4,  TX_8X8,   TX_16X16, TX_32X32, TX_4X8,  TX_8X4,  TX_8X16,
    TX_16X8, TX_16X32, TX_32X16, TX_4X16,  TX_16X4, TX_8X32, TX_32X8};

static const char *const kTxSizeStrings[] = {
    "TX_4X4",  "TX_8X8",  "TX_16X16", "TX_32X32", "TX_4X8",
    "TX_8X4",  "TX_8X16", "TX_16X8",  "TX_16X32", "TX_32X16",
    "TX_4X16", "TX_16X4", "TX_8X32",  "TX_32X8"};

template <typename F, typename T>
class DrPredictionTest : public FunctionEquivalenceTest<F> {
 protected:
  static const int kIterations = 10;
  static const int kBufSize = 2 * 64 + 32 + 0;
  static const int kDstStride = 64;
  static const int kDstSize = kDstStride * kDstStride;
  static const int kOffset = 16;

  virtual ~DrPredictionTest() {}

  virtual void Execute() = 0;

  void Common() {
    dst_ref_ = &dst_ref_data_[0];
    dst_tst_ = &dst_tst_data_[0];
    dst_stride_ = kDstStride;
    above_ = &above_data_[kOffset];
    left_ = &left_data_[kOffset];

    upsample_above_ = 0;
    upsample_left_ = 0;

    for (int i = 0; i < kBufSize; ++i) {
      above_data_[i] = this->rng_.Rand8();
      left_data_[i] = this->rng_.Rand8();
    }
  }

  void RunTest() {
    for (int tx = 0; tx < 14; ++tx) {
      for (int i = 0; i < kDstSize; ++i) {
        dst_ref_[i] = dst_tst_[i] = 0x0;
      }

      bw_ = tx_size_wide[kTxSize[tx]];
      bh_ = tx_size_high[kTxSize[tx]];

      Execute();

      for (int r = 0; r < bh_; ++r) {
        for (int c = 0; c < bw_; ++c) {
          ASSERT_EQ(dst_ref_[r * dst_stride_ + c],
                    dst_tst_[r * dst_stride_ + c])
              << bw_ << "x" << bh_ << " r: " << r << " c: " << c
              << " dx: " << dx_ << " dy: " << dy_
              << " upsample_above: " << upsample_above_
              << " upsample_left: " << upsample_left_;
        }
      }
    }
  }

  T dst_ref_data_[kDstSize];
  T dst_tst_data_[kDstSize];

  T above_data_[kBufSize];
  T left_data_[kBufSize];

  T *dst_ref_;
  T *dst_tst_;
  T *above_;
  T *left_;

  int size_;
  int bw_;
  int bh_;
  int upsample_above_;
  int upsample_left_;
  int dst_stride_;
  int dx_;
  int dy_;
};

//////////////////////////////////////////////////////////////////////////////
// av1_dr_prediction_z1
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ1)(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                     const uint8_t *above, const uint8_t *left,
                     int upsample_above, int dx, int dy);
typedef libaom_test::FuncParam<DRZ1> TestFuncsZ1;

class DrPredZ1 : public DrPredictionTest<DRZ1, uint8_t> {
 protected:
  void Execute() {
    params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                     upsample_above_, dx_, dy_);
    ASM_REGISTER_STATE_CHECK(params_.tst_func(dst_tst_, dst_stride_, bw_, bh_,
                                              above_, left_, upsample_above_,
                                              dx_, dy_));
  }
};

TEST_P(DrPredZ1, RandomValues) {
  Common();

  dy_ = 1;

  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    for (int angle = 0; angle < 90; ++angle) {
      dx_ = dr_intra_derivative[angle];
      // pass in all dx values even though zero values are never used
      RunTest();
    }
  }
}

// #define CHECK_ALL_ANGLES

TEST_P(DrPredZ1, DISABLED_Speed) {
  const int test_count = 100000;
#ifndef CHECK_ALL_ANGLES
  const int angles[] = {3, 45, 87};
#endif
  Common();

  dy_ = 1;
  dx_ = 64;
#ifndef CHECK_ALL_ANGLES
  for (int i = 0; i < 3; ++i) {
    dx_= dr_intra_derivative[angles[i]];
    printf("angle %2d\n", angles[i]);
#else
  for (int i = 0; i < 90; i += 1) {
    dx_= dr_intra_derivative[i];
#endif
    printf("dx_ %2d\n", dx_);
    if (dx_)
    for (int tx = 0; tx < 14; ++tx) {
      bw_ = tx_size_wide[kTxSize[tx]];
      bh_ = tx_size_high[kTxSize[tx]];

      aom_usec_timer timer;

      aom_usec_timer_start(&timer);
      for (int n = 0; n < test_count; ++n) {
        params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                         upsample_above_, dx_, dy_);
      }
      aom_usec_timer_mark(&timer);
      const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

      aom_usec_timer_start(&timer);
      for (int n = 0; n < test_count; ++n) {
        params_.tst_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                         upsample_above_, dx_, dy_);
      }
      aom_usec_timer_mark(&timer);
      const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));
      float x = (float)ref_time / (float)tst_time;
      if (ref_time < tst_time)
        x *= -1;

      printf("\t[%8s] :: ref time %6d, tst time %6d     %3.2f\n",
          kTxSizeStrings[tx], ref_time, tst_time, x);
    }
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, DrPredZ1,
    ::testing::Values(TestFuncsZ1(av1_dr_prediction_z1_c,
                                  av1_dr_prediction_z1_sse4_1)));
#endif  // HAVE_SSE4_1

//////////////////////////////////////////////////////////////////////////////
// av1_dr_prediction_z2
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ2)(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                     const uint8_t *above, const uint8_t *left,
                     int upsample_above, int upsample_left, int dx, int dy);
typedef libaom_test::FuncParam<DRZ2> TestFuncsZ2;

class DrPredZ2 : public DrPredictionTest<DRZ2, uint8_t> {
 protected:
  void Execute() {
    params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                     upsample_above_, upsample_left_, dx_, dy_);
    ASM_REGISTER_STATE_CHECK(params_.tst_func(dst_tst_, dst_stride_, bw_, bh_,
                                              above_, left_, upsample_above_,
                                              upsample_left_, dx_, dy_));
  }
};

TEST_P(DrPredZ2, RandomValues) {
  Common();

  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    upsample_left_ = iter & 1;
    for (int angle = 90; angle < 180; ++angle) {
      dx_ = dr_intra_derivative[180 - angle];
      dy_ = dr_intra_derivative[angle - 90];
      // pass in all dx, dy values even though zero values are never used
      RunTest();
    }
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, DrPredZ2,
    ::testing::Values(TestFuncsZ2(av1_dr_prediction_z2_c,
                                  av1_dr_prediction_z2_sse4_1)));
#endif  // HAVE_SSE4_1

//////////////////////////////////////////////////////////////////////////////
// av1_dr_prediction_z3
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ3)(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                     const uint8_t *above, const uint8_t *left,
                     int upsample_left, int dx, int dy);
typedef libaom_test::FuncParam<DRZ3> TestFuncsZ3;

class DrPredZ3 : public DrPredictionTest<DRZ3, uint8_t> {
 protected:
  void Execute() {
    params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                     upsample_left_, dx_, dy_);
    ASM_REGISTER_STATE_CHECK(params_.tst_func(dst_tst_, dst_stride_, bw_, bh_,
                                              above_, left_, upsample_left_,
                                              dx_, dy_));
  }
};

TEST_P(DrPredZ3, RandomValues) {
  Common();

  dx_ = 1;

  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_left_ = iter & 1;
    for (int angle = 180; angle < 270; ++angle) {
      dy_ = dr_intra_derivative[270 - angle];
      // pass in all dy values even though zero values are never used
      RunTest();
    }
  }
}

#if HAVE_SSE4_1
INSTANTIATE_TEST_CASE_P(
    SSE4_1, DrPredZ3,
    ::testing::Values(TestFuncsZ3(av1_dr_prediction_z3_c,
                                  av1_dr_prediction_z3_sse4_1)));
#endif  // HAVE_SSE4_1

}  // namespace
