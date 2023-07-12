/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "config/aom_dsp_rtcd.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/util.h"
#include "test/register_state_check.h"
#include "test/yuv_video_source.h"

#include "aom_dsp/flow_estimation/disflow.h"

namespace {

typedef void (*compute_flow_at_point_func)(const uint8_t *src,
                                           const uint8_t *ref, int x, int y,
                                           int width, int height, int stride,
                                           double *u, double *v);

class ComputeFlowTest
    : public ::testing::TestWithParam<compute_flow_at_point_func> {
 public:
  ~ComputeFlowTest();
  void SetUp();

  void TearDown();

 protected:
  void RunCheckOutput(int run_times);
  compute_flow_at_point_func target_func;
};
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(ComputeFlowTest);

ComputeFlowTest::~ComputeFlowTest() {}

void ComputeFlowTest::SetUp() { target_func = GetParam(); }

void ComputeFlowTest::TearDown() {}

void ComputeFlowTest::RunCheckOutput(int run_times) {
  ::libaom_test::YUVVideoSource video("bus_352x288_420_f20_b8.yuv",
                                      AOM_IMG_FMT_I420, 352, 288, 30, 1, 0, 2);
  // Use Y (Luminance) plane.
  video.Begin();
  uint8_t *src = video.img()->planes[0];
  assert(src != nullptr);
  video.Next();
  uint8_t *ref = video.img()->planes[0];
  assert(ref != nullptr);

  int width = 352;
  int height = 288;
  int stride = 288;

  // Pick a random value between -5 and 5. The range was chosen arbitrarily as
  // u and v can take any kind of value in practise, but it shouldn't change the
  // outcome of the tests.
  double u_rand = ((double)rand() / RAND_MAX) * 10 - 5;
  double u_ref = u_rand;
  double u_test = u_rand;

  double v_rand = ((double)rand() / RAND_MAX) * 10 - 5;
  double v_ref = v_rand;
  double v_test = v_rand;

  // Pick a random point in the frame. If the frame is 352x288, that means we
  // can call the function on all values of x comprised between 8 and 344, and
  // all values of y comprised between 8 and 280.
  int x = (rand() % (344 - 8 + 1)) + 8;
  int y = (rand() % (280 - 8 + 1)) + 8;

  aom_usec_timer ref_timer, test_timer;

  aom_compute_flow_at_point_c(src, ref, x, y, width, height, stride, &u_ref,
                              &v_ref);

  target_func(src, ref, x, y, width, height, stride, &u_test, &v_test);

  if (run_times > 1) {
    aom_usec_timer_start(&ref_timer);
    for (int i = 0; i < run_times; i++) {
      aom_compute_flow_at_point_c(src, ref, x, y, width, height, stride, &u_ref,
                                  &v_ref);
    }
    aom_usec_timer_mark(&ref_timer);
    const double elapsed_time_c =
        static_cast<double>(aom_usec_timer_elapsed(&ref_timer));

    aom_usec_timer_start(&test_timer);
    for (int i = 0; i < run_times; i++) {
      target_func(src, ref, x, y, width, height, stride, &u_test, &v_test);
    }
    aom_usec_timer_mark(&test_timer);
    const double elapsed_time_simd =
        static_cast<double>(aom_usec_timer_elapsed(&test_timer));

    printf("c_time=%fns \t simd_time=%fns \t speedup=%.2f\n", elapsed_time_c,
           elapsed_time_simd, (elapsed_time_c / elapsed_time_simd));
  } else {
    ASSERT_EQ(u_ref, u_test);
    ASSERT_EQ(v_ref, v_test);
  }
}

TEST_P(ComputeFlowTest, CheckOutput) { RunCheckOutput(1); }

TEST_P(ComputeFlowTest, DISABLED_Speed) { RunCheckOutput(10000000); }

#if HAVE_SSE4_1
INSTANTIATE_TEST_SUITE_P(SSE4_1, ComputeFlowTest,
                         ::testing::Values(aom_compute_flow_at_point_sse4_1));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_SUITE_P(NEON, ComputeFlowTest,
                         ::testing::Values(aom_compute_flow_at_point_neon));
#endif

}  // namespace
