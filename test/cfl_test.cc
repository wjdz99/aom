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
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom_ports/aom_timer.h"
#include "./av1_rtcd.h"
#include "test/util.h"
#include "test/acm_random.h"

using std::tr1::make_tuple;

using libaom_test::ACMRandom;

#define AVERAGE_ITERATIONS (100)

namespace {
typedef void (*subtract_fn)(int16_t *pred_buf_q3, int width, int height,
                            int16_t avg_q3);

typedef std::tr1::tuple<int, int, subtract_fn> CFL_Param;

#define ALL_SIZES_CFL(subtract_fn)                                        \
  make_tuple(4, 4, &subtract_fn), make_tuple(8, 4, &subtract_fn),         \
      make_tuple(4, 8, &subtract_fn), make_tuple(8, 8, &subtract_fn),     \
      make_tuple(16, 8, &subtract_fn), make_tuple(8, 16, &subtract_fn),   \
      make_tuple(16, 16, &subtract_fn), make_tuple(32, 16, &subtract_fn), \
      make_tuple(16, 32, &subtract_fn), make_tuple(32, 32, &subtract_fn)

class CFLAverageTest : public ::testing::TestWithParam<CFL_Param> {
 public:
  virtual ~CFLAverageTest() {}
  virtual void SetUp() { subtract = GET_PARAM(2); }

 protected:
  subtract_fn subtract;
  int Width() const { return GET_PARAM(0); }
  int Height() const { return GET_PARAM(1); }
};

class CFLAverageSpeedTest : public ::testing::TestWithParam<CFL_Param> {
 public:
  virtual ~CFLAverageSpeedTest() {}
  virtual void SetUp() {
    subtract = GET_PARAM(2);
    int16_t *pred_buf_q3 = pred_buf_q3_data;
    int k = 0;
    for (int j = 0; j < CFL_BUF_LINE; j++) {
      for (int i = 0; i < CFL_BUF_LINE; i++) {
        pred_buf_q3[i] = k++;
      }
      pred_buf_q3 += CFL_BUF_LINE;
    }
  }

 protected:
  subtract_fn subtract;
  int Width() const { return GET_PARAM(0); }
  int Height() const { return GET_PARAM(1); }
  int16_t pred_buf_q3_data[CFL_BUF_SQUARE];
};

static INLINE void setup(int16_t *pred_buf_q3) {
  int k = 0;
  for (int j = 0; j < CFL_BUF_LINE; j++) {
    for (int i = 0; i < CFL_BUF_LINE; i++) {
      pred_buf_q3[i] = k++;
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}

TEST_P(CFLAverageSpeedTest, DISABLED_AverageSpeedTest) {
  const int width = Width();
  const int height = Height();
  const int iterations = INT16_MAX;
  int16_t *pred_buf_q3 = pred_buf_q3_data;

  aom_usec_timer ref_timer;
  aom_usec_timer timer;

  setup(pred_buf_q3);
  aom_usec_timer_start(&ref_timer);
  for (int k = 0; k < iterations; k++) {
    av1_subtract_average_c(pred_buf_q3, width, height, k);
  }
  aom_usec_timer_mark(&ref_timer);
  int ref_elapsed_time = (int)aom_usec_timer_elapsed(&ref_timer);

  setup(pred_buf_q3);
  aom_usec_timer_start(&timer);
  for (int k = 0; k < iterations; k++) {
    subtract(pred_buf_q3, width, height, k);
  }
  aom_usec_timer_mark(&timer);
  int elapsed_time = (int)aom_usec_timer_elapsed(&timer);

#if 1
  std::cout.precision(2);
  std::cout << "[          ] " << width << "x" << height
            << ": C time = " << ref_elapsed_time
            << " us, SIMD time = " << elapsed_time << " us"
            << " (~" << ref_elapsed_time / (double)elapsed_time << "x) "
            << std::endl;
#endif

  EXPECT_GT(ref_elapsed_time, elapsed_time)
      << "Error: CDEFSpeedTest, SIMD slower than C." << std::endl
      << "C time: " << ref_elapsed_time << " us" << std::endl
      << "SIMD time: " << elapsed_time << " us" << std::endl;
}

TEST_P(CFLAverageTest, AverageTest) {
  const int width = Width();
  const int height = Height();

  ACMRandom rnd(ACMRandom::DeterministicSeed());

  // 32x32 CfL Prediction Buffer
  int16_t pred_buf_q3_data[CFL_BUF_SQUARE];

  for (int it = 0; it < AVERAGE_ITERATIONS; it++) {
    int16_t *pred_buf_q3 = pred_buf_q3_data;
    int16_t k = 0;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        pred_buf_q3[i] = k++;
      }
      pred_buf_q3 += CFL_BUF_LINE;
    }
    k = rnd.Rand15Signed();
    subtract(pred_buf_q3_data, width, height, k);
    pred_buf_q3 = pred_buf_q3_data;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        ASSERT_EQ(pred_buf_q3[i], -k);
        k--;
      }
      pred_buf_q3 += CFL_BUF_LINE;
    }
  }
}

const CFL_Param kArrayCFL_c[] = { ALL_SIZES_CFL(av1_subtract_average_c) };
INSTANTIATE_TEST_CASE_P(C, CFLAverageTest, ::testing::ValuesIn(kArrayCFL_c));

#if HAVE_SSE2
const CFL_Param kArrayCFL_sse2[] = { ALL_SIZES_CFL(av1_subtract_average_sse2) };
INSTANTIATE_TEST_CASE_P(SSE2, CFLAverageTest,
                        ::testing::ValuesIn(kArrayCFL_sse2));

INSTANTIATE_TEST_CASE_P(SSE2, CFLAverageSpeedTest,
                        ::testing::ValuesIn(kArrayCFL_sse2));
#endif

#if HAVE_AVX2
const CFL_Param kArrayCFL_avx2[] = { ALL_SIZES_CFL(av1_subtract_average_avx2) };
INSTANTIATE_TEST_CASE_P(AVX2, CFLAverageTest,
                        ::testing::ValuesIn(kArrayCFL_avx2));

INSTANTIATE_TEST_CASE_P(AVX2, CFLAverageSpeedTest,
                        ::testing::ValuesIn(kArrayCFL_avx2));
#endif
}  // namespace
