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

#include "aom/aom_integer.h"
#include "aom_ports/aom_timer.h"
#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "config/av1_rtcd.h"

#include "av1/common/enums.h"

using libaom_test::FunctionEquivalenceTest;

namespace {

static const TX_SIZE kTxSize[] = { TX_4X4,   TX_8X8,   TX_16X16, TX_32X32,
                                   TX_64X64, TX_4X8,   TX_8X4,   TX_8X16,
                                   TX_16X8,  TX_16X32, TX_32X16, TX_32X64,
                                   TX_64X32, TX_4X16,  TX_16X4,  TX_8X32,
                                   TX_32X8,  TX_16X64, TX_64X16 };

static const char *const kTxSizeStrings[] = {
  "TX_4X4",   "TX_8X8",   "TX_16X16", "TX_32X32", "TX_64X64",
  "TX_4X8",   "TX_8X4",   "TX_8X16",  "TX_16X8",  "TX_16X32",
  "TX_32X16", "TX_32X64", "TX_64X32", "TX_4X16",  "TX_16X4",
  "TX_8X32",  "TX_32X8",  "TX_16X64", "TX_64X16"
};

template <typename PredictorFunc, typename Pixel>
class DrPredictionTest : public FunctionEquivalenceTest<PredictorFunc> {
 protected:
  static const int kMaxNumTests = 100000;
  static const int kIterations = 10;
  static const int kDstStride = 64;
  static const int kDstSize = kDstStride * kDstStride;
  static const int kOffset = 16;
  static const int kBufSize = (2 * MAX_TX_SIZE) << 1 + 16;

  DrPredictionTest() {
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

    for (int i = 0; i < kDstSize; ++i) {
      dst_ref_[i] = 0;
    }

    bw_ = 0;
    bh_ = 0;
    dx_ = 1;
    dy_ = 1;

    bd_ = 8;
    txsize_ = TX_4X4;
  }

  virtual ~DrPredictionTest() {}

  virtual void Execute(int speedtest, int tx) = 0;

  void RunTest(int speedtest) {
    for (int i = 0; i < kBufSize; ++i) {
      above_data_[i] = left_data_[i] = (1 << bd_) - 1;
    }

    for (int tx = 0; tx < TX_SIZES_ALL; ++tx) {
      if (this->params_.tst_func == NULL) {
        for (int i = 0; i < kDstSize; ++i) {
          dst_tst_[i] = (1 << bd_) - 1;
        }
      } else {
        for (int i = 0; i < kDstSize; ++i) {
          dst_tst_[i] = 0;
        }
      }

      bw_ = tx_size_wide[kTxSize[tx]];
      bh_ = tx_size_high[kTxSize[tx]];

      Execute(speedtest, tx);

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

  void RunTest_HBD(int speedtest) {
    bd_ = this->params_.bit_depth;
    RunTest(speedtest);
  }

  void OutputTimes(int kNumTests, int ref_time, int tst_time, int tx) {
    if (kNumTests > 1) {
      if (this->params_.tst_func) {
        float x = (float)ref_time / (float)tst_time;
        printf("\t[%8s] :: ref time %6d, tst time %6d     %3.2f\n",
               kTxSizeStrings[tx], ref_time, tst_time, x);
      } else {
        printf("\t[%8s] :: ref time %6d\n", kTxSizeStrings[tx], ref_time);
      }
    }
  }

  Pixel dst_ref_data_[kDstSize];
  Pixel dst_tst_data_[kDstSize];

  Pixel left_data_[kBufSize];
  Pixel dummy_data_[kBufSize];
  Pixel above_data_[kBufSize];

  Pixel *dst_ref_;
  Pixel *dst_tst_;
  Pixel *above_;
  Pixel *left_;

  int bw_;
  int bh_;
  int upsample_above_;
  int upsample_left_;
  int dst_stride_;
  int dx_;
  int dy_;
  int bd_;
  TX_SIZE txsize_;
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
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_above_, dx_, dy_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(params_.tst_func(dst_tst_, dst_stride_, bw_,
                                                  bh_, above_, left_,
                                                  upsample_above_, dx_, dy_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ1, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    for (int angle = 0; angle < 90; ++angle) {
      dx_ = dr_intra_derivative[angle];
      if (dx_) RunTest(0);
    }
  }
}

TEST_P(DrPredZ1, DISABLED_Speed) {
  const int angles[] = { 3, 45, 87 };
  for (upsample_above_ = 0; upsample_above_ < 2; ++upsample_above_) {
    for (int i = 0; i < 3; ++i) {
      dx_ = dr_intra_derivative[angles[i]];
      if (dx_) RunTest(1);
      printf("upsample_above: %d angle: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_above_, angles[i]);
    }
  }
}

INSTANTIATE_TEST_CASE_P(C, DrPredZ1,
                        ::testing::Values(TestFuncsZ1(av1_dr_prediction_z1_c,
                                                      NULL)));

//////////////////////////////////////////////////////////////////////////////
// av1_dr_prediction_z2
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ2)(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                     const uint8_t *above, const uint8_t *left,
                     int upsample_above, int upsample_left, int dx, int dy);
typedef libaom_test::FuncParam<DRZ2> TestFuncsZ2;

class DrPredZ2 : public DrPredictionTest<DRZ2, uint8_t> {
 protected:
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_above_, upsample_left_, dx_, dy_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(
            params_.tst_func(dst_tst_, dst_stride_, bw_, bh_, above_, left_,
                             upsample_above_, upsample_left_, dx_, dy_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ2, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    upsample_left_ = iter & 1;
    for (int angle = 90; angle < 180; ++angle) {
      dx_ = dr_intra_derivative[180 - angle];
      dy_ = dr_intra_derivative[angle - 90];
      if (dx_ && dy_) RunTest(0);
    }
  }
}

TEST_P(DrPredZ2, DISABLED_Speed) {
  const int angles[] = { 3 + 90, 45 + 90, 87 + 90 };
  for (upsample_above_ = 0; upsample_above_ < 2; ++upsample_above_) {
    upsample_left_ = upsample_above_;
    for (int i = 0; i < 3; ++i) {
      dx_ = dr_intra_derivative[180 - angles[i]];
      dy_ = dr_intra_derivative[angles[i] - 90];
      if (dx_ && dy_) RunTest(1);
      printf("upsample: %d angle: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_above_, angles[i]);
    }
  }
}

INSTANTIATE_TEST_CASE_P(C, DrPredZ2,
                        ::testing::Values(TestFuncsZ2(av1_dr_prediction_z2_c,
                                                      NULL)));

//////////////////////////////////////////////////////////////////////////////
// av1_dr_prediction_z3
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ3)(uint8_t *dst, ptrdiff_t stride, int bw, int bh,
                     const uint8_t *above, const uint8_t *left,
                     int upsample_left, int dx, int dy);
typedef libaom_test::FuncParam<DRZ3> TestFuncsZ3;

class DrPredZ3 : public DrPredictionTest<DRZ3, uint8_t> {
 protected:
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_left_, dx_, dy_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(params_.tst_func(dst_tst_, dst_stride_, bw_,
                                                  bh_, above_, left_,
                                                  upsample_left_, dx_, dy_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ3, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_left_ = iter & 1;
    for (int angle = 180; angle < 270; ++angle) {
      dy_ = dr_intra_derivative[270 - angle];
      if (dy_) RunTest(0);
    }
  }
}

TEST_P(DrPredZ3, DISABLED_Speed) {
  const int angles[] = { 3 + 180, 45 + 180, 87 + 180 };
  for (upsample_left_ = 0; upsample_left_ < 2; ++upsample_left_) {
    for (int i = 0; i < 3; ++i) {
      dy_ = dr_intra_derivative[270 - angles[i]];
      if (dy_) RunTest(1);
      printf("upsample_left_: %d angle: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_left_, angles[i]);
    }
  }
}

INSTANTIATE_TEST_CASE_P(C, DrPredZ3,
                        ::testing::Values(TestFuncsZ3(av1_dr_prediction_z3_c,
                                                      NULL)));

//////////////////////////////////////////////////////////////////////////////
// av1_highbd_dr_prediction_z1
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ1_HBD)(uint16_t *dst, ptrdiff_t stride, int bw, int bh,
                         const uint16_t *above, const uint16_t *left,
                         int upsample_above, int dx, int dy, int bd);
typedef libaom_test::FuncParam<DRZ1_HBD> TestFuncsZ1_HBD;

class DrPredZ1_HBD : public DrPredictionTest<DRZ1_HBD, uint16_t> {
 protected:
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_above_, dx_, dy_, bd_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(
            params_.tst_func(dst_tst_, dst_stride_, bw_, bh_, above_, left_,
                             upsample_above_, dx_, dy_, bd_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ1_HBD, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    for (int angle = 0; angle < 90; ++angle) {
      dx_ = dr_intra_derivative[angle];
      if (dx_) RunTest_HBD(0);
    }
  }
}

TEST_P(DrPredZ1_HBD, DISABLED_Speed) {
  const int angles[] = { 3, 45, 87 };
  for (upsample_above_ = 0; upsample_above_ < 2; ++upsample_above_) {
    for (int i = 0; i < 3; ++i) {
      dx_ = dr_intra_derivative[angles[i]];
      if (dx_) RunTest_HBD(1);
      printf("upsample_above: %d angle: %d bd: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_above_, angles[i], bd_);
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    C, DrPredZ1_HBD,
    ::testing::Values(
        TestFuncsZ1_HBD(av1_highbd_dr_prediction_z1_c, NULL, AOM_BITS_8),
        TestFuncsZ1_HBD(av1_highbd_dr_prediction_z1_c, NULL, AOM_BITS_10),
        TestFuncsZ1_HBD(av1_highbd_dr_prediction_z1_c, NULL, AOM_BITS_12)));

//////////////////////////////////////////////////////////////////////////////
// av1_highbd_dr_prediction_z2
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ2_HBD)(uint16_t *dst, ptrdiff_t stride, int bw, int bh,
                         const uint16_t *above, const uint16_t *left,
                         int upsample_above, int upsample_left, int dx, int dy,
                         int bd);
typedef libaom_test::FuncParam<DRZ2_HBD> TestFuncsZ2_HBD;

class DrPredZ2_HBD : public DrPredictionTest<DRZ2_HBD, uint16_t> {
 protected:
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_above_, upsample_left_, dx_, dy_, bd_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(
            params_.tst_func(dst_tst_, dst_stride_, bw_, bh_, above_, left_,
                             upsample_above_, upsample_left_, dx_, dy_, bd_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ2_HBD, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_above_ = iter & 1;
    upsample_left_ = iter & 1;
    for (int angle = 90; angle < 180; ++angle) {
      dx_ = dr_intra_derivative[180 - angle];
      dy_ = dr_intra_derivative[angle - 90];
      if (dx_ && dy_) RunTest_HBD(0);
    }
  }
}

TEST_P(DrPredZ2_HBD, DISABLED_Speed) {
  const int angles[] = { 3 + 90, 45 + 90, 87 + 90 };
  for (upsample_above_ = 0; upsample_above_ < 2; ++upsample_above_) {
    upsample_left_ = upsample_above_;
    for (int i = 0; i < 3; ++i) {
      dx_ = dr_intra_derivative[180 - angles[i]];
      dy_ = dr_intra_derivative[angles[i] - 90];
      if (dx_ && dy_) RunTest_HBD(1);
      printf("upsample: %d angle: %d bd: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_above_, angles[i], bd_);
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    C, DrPredZ2_HBD,
    ::testing::Values(
        TestFuncsZ2_HBD(av1_highbd_dr_prediction_z2_c, NULL, AOM_BITS_8),
        TestFuncsZ2_HBD(av1_highbd_dr_prediction_z2_c, NULL, AOM_BITS_10),
        TestFuncsZ2_HBD(av1_highbd_dr_prediction_z2_c, NULL, AOM_BITS_12)));

//////////////////////////////////////////////////////////////////////////////
// av1_highbd_dr_prediction_z3
//////////////////////////////////////////////////////////////////////////////

typedef void (*DRZ3_HBD)(uint16_t *dst, ptrdiff_t stride, int bw, int bh,
                         const uint16_t *above, const uint16_t *left,
                         int upsample_left, int dx, int dy, int bd);
typedef libaom_test::FuncParam<DRZ3_HBD> TestFuncsZ3_HBD;

class DrPredZ3_HBD : public DrPredictionTest<DRZ3_HBD, uint16_t> {
 protected:
  void Execute(int speedtest, int tx) {
    const int kNumTests = speedtest ? kMaxNumTests : 1;
    aom_usec_timer timer;

    aom_usec_timer_start(&timer);

    for (int k = 0; k < kNumTests; ++k) {
      params_.ref_func(dst_ref_, dst_stride_, bw_, bh_, above_, left_,
                       upsample_left_, dx_, dy_, bd_);
    }
    aom_usec_timer_mark(&timer);
    const int ref_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    aom_usec_timer_start(&timer);

    if (params_.tst_func) {
      for (int k = 0; k < kNumTests; ++k) {
        ASM_REGISTER_STATE_CHECK(
            params_.tst_func(dst_tst_, dst_stride_, bw_, bh_, above_, left_,
                             upsample_left_, dx_, dy_, bd_));
      }
    }
    aom_usec_timer_mark(&timer);
    const int tst_time = static_cast<int>(aom_usec_timer_elapsed(&timer));

    OutputTimes(kNumTests, ref_time, tst_time, tx);
  }
};

TEST_P(DrPredZ3_HBD, SaturatedValues) {
  for (int iter = 0; iter < kIterations && !HasFatalFailure(); ++iter) {
    upsample_left_ = iter & 1;
    for (int angle = 180; angle < 270; ++angle) {
      dy_ = dr_intra_derivative[270 - angle];
      if (dy_) RunTest_HBD(0);
    }
  }
}

TEST_P(DrPredZ3_HBD, DISABLED_Speed) {
  const int angles[] = { 3 + 180, 45 + 180, 87 + 180 };
  for (upsample_left_ = 0; upsample_left_ < 2; ++upsample_left_) {
    for (int i = 0; i < 3; ++i) {
      dy_ = dr_intra_derivative[270 - angles[i]];
      if (dy_) RunTest_HBD(1);
      printf("upsample_left_: %d angle: %d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
             upsample_left_, angles[i]);
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    C, DrPredZ3_HBD,
    ::testing::Values(
        TestFuncsZ3_HBD(av1_highbd_dr_prediction_z3_c, NULL, AOM_BITS_8),
        TestFuncsZ3_HBD(av1_highbd_dr_prediction_z3_c, NULL, AOM_BITS_10),
        TestFuncsZ3_HBD(av1_highbd_dr_prediction_z3_c, NULL, AOM_BITS_12)));
}  // namespace
