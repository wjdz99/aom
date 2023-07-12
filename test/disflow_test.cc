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
#include "test/y4m_video_source.h"

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
  ::libaom_test::Y4mVideoSource video("Vertical_Bayshore_270x480_2997.y4m", 0,
                                      2);
  // Use Y (Luminance) plane.
  video.Begin();
  uint8_t *src = video.img()->planes[0];
  assert(src != nullptr);
  video.Next();
  uint8_t *ref = video.img()->planes[0];
  assert(ref != nullptr);

  int width = 270;
  int height = 480;
  int stride = 270;

  double u_init = 10.678;
  double u_ref = u_init;
  double u_test = u_init;

  double v_init = 0.5674;
  double v_ref = v_init;
  double v_test = v_init;

  // Pick a random point in the frame. If the frame is 270x480, that means we
  // can call the function on all values of x comprised between 8 and 262, and
  // all values of y comprised between 8 and 472.
  int x = (rand() % (262 - 8 + 1)) + 8;
  int y = (rand() % (472 - 8 + 1)) + 8;

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

}  // namespace
