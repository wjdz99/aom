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
#include <stdio.h>
#include <stdlib.h>

#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom_mem/aom_mem.h"
#include "aom_ports/aom_timer.h"
#include "aom_ports/mem.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/util.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace AV1CompRoundShift {

typedef void (*comp_round_shift_array_func)(int32_t *arr, int size, int bit);

#if HAVE_SSSE3 || HAVE_NEON
const BLOCK_SIZE kValidBlockSize[] = {
  BLOCK_4X4,   BLOCK_8X8,   BLOCK_16X16, BLOCK_32X32, BLOCK_64X64,
  BLOCK_4X8,   BLOCK_8X4,   BLOCK_8X16,  BLOCK_16X8,  BLOCK_16X32,
  BLOCK_32X16, BLOCK_32X64, BLOCK_64X32, BLOCK_4X16,  BLOCK_16X4,
  BLOCK_8X32,  BLOCK_32X8,  BLOCK_16X64, BLOCK_64X16,
};
const int kValidBitCheck[] = {
  -4, -3, -2, -1, 0, 1, 2, 3, 4,
};
#endif

typedef ::testing::tuple<comp_round_shift_array_func, BLOCK_SIZE, int>
    CompRoundShiftParam;

class AV1CompRoundShiftTest
    : public ::testing::TestWithParam<CompRoundShiftParam> {
 public:
  ~AV1CompRoundShiftTest();

  void SetUp();

  void TearDown();

 protected:
  void RunCheckOutput(comp_round_shift_array_func test_impl, BLOCK_SIZE bsize,
                      int bit);
  void RunSpeedTest(comp_round_shift_array_func test_impl, BLOCK_SIZE bsize,
                    int bit);

  libaom_test::ACMRandom rnd_;
};

AV1CompRoundShiftTest::~AV1CompRoundShiftTest() { ; }

void AV1CompRoundShiftTest::SetUp() {
  rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
}

void AV1CompRoundShiftTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1CompRoundShiftTest::RunCheckOutput(
    comp_round_shift_array_func test_impl, BLOCK_SIZE bsize, int bit) {
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  const int blk_size = 64;
  DECLARE_ALIGNED(32, int32_t, pred_[blk_size]) = { 0 };
  DECLARE_ALIGNED(32, int32_t, ref_buffer_[blk_size]) = { 0 };
  for (int i = 0; i < (blk_size); ++i) {
    pred_[i] = rnd_.Rand31();
    ref_buffer_[i] = pred_[i];
  }
  av1_round_shift_array_c(ref_buffer_, w, bit);
  test_impl(pred_, w, bit);
  for (int x = 0; x < w; ++x) {
    ASSERT_EQ(ref_buffer_[x], pred_[x]) << w << "x" << h << "mismatch @"
                                        << "(" << x << ")";
  }
}

void AV1CompRoundShiftTest::RunSpeedTest(comp_round_shift_array_func test_impl,
                                         BLOCK_SIZE bsize, int bit) {
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  const int blk_size = 64;
  DECLARE_ALIGNED(32, int32_t, ref_buffer_[blk_size]) = { 0 };

  for (int i = 0; i < (blk_size); ++i) {
    ref_buffer_[i] = rnd_.Rand31();
  }

  const int num_loops = 1000000000 / (w + h);
  comp_round_shift_array_func funcs[2] = { av1_round_shift_array_c, test_impl };
  double elapsed_time[2] = { 0 };
  for (int i = 0; i < 2; ++i) {
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);
    comp_round_shift_array_func func = funcs[i];
    for (int j = 0; j < num_loops; ++j) {
      func(ref_buffer_, w, bit);
    }
    aom_usec_timer_mark(&timer);
    double time = static_cast<double>(aom_usec_timer_elapsed(&timer));
    elapsed_time[i] = 1000.0 * time / num_loops;
  }
  printf("compRoundShift %3dx%-3d: %7.2f/%7.2fns", w, h, elapsed_time[0],
         elapsed_time[1]);
  printf("(%3.2f)\n", elapsed_time[0] / elapsed_time[1]);
}

TEST_P(AV1CompRoundShiftTest, CheckOutput) {
  RunCheckOutput(GET_PARAM(0), GET_PARAM(1), GET_PARAM(2));
}

TEST_P(AV1CompRoundShiftTest, Speed) {
  RunSpeedTest(GET_PARAM(0), GET_PARAM(1), GET_PARAM(2));
}

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
    SSSE3, AV1CompRoundShiftTest,
    ::testing::Combine(::testing::Values(&av1_round_shift_array_sse4_1),
                       ::testing::ValuesIn(kValidBlockSize),
                       ::testing::ValuesIn(kValidBitCheck)));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(
    NEON, AV1CompRoundShiftTest,
    ::testing::Combine(::testing::Values(&av1_round_shift_array_neon),
                       ::testing::ValuesIn(kValidBlockSize),
                       ::testing::ValuesIn(kValidBitCheck)));
#endif

};  // namespace AV1CompRoundShift
