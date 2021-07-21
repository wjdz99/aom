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

#include <tuple>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom/aom_integer.h"
#include "aom_ports/aom_timer.h"
#include "av1/encoder/ml.h"
#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"
#include "test/util.h"
#include "test/register_state_check.h"
#include "test/acm_random.h"

namespace {
typedef void (*FastSoftmaxFn)(const float *const input, float *output);

typedef std::tuple<const FastSoftmaxFn, int> FastSoftmaxTestParams;

const float rel_epsilon = 5e-2f;  // Error threshold for functional equivalence
const float abs_epsilon = 5e-3f;  // Error threshold for functional equivalence

class FastSoftmaxTest : public ::testing::TestWithParam<FastSoftmaxTestParams> {
 public:
  FastSoftmaxTest() : _target_fn{ GET_PARAM(0) }, _num_classes(GET_PARAM(1)) {
    _ref_buf = (float *)calloc(_num_classes, sizeof(*_ref_buf));
    _dst_buf = (float *)calloc(_num_classes, sizeof(*_dst_buf));
    _input = (float *)calloc(_num_classes, sizeof(*_input));
  }
  ~FastSoftmaxTest() {
    free(_ref_buf);
    free(_dst_buf);
    free(_input);
  }
  void RunSoftmaxTest();
  void RunSoftmaxSpeedTest(const int run_times);
  void FillInputBuf();

 private:
  const FastSoftmaxFn _target_fn;
  const int _num_classes;
  float *_ref_buf, *_dst_buf, *_input;
  libaom_test::ACMRandom rng_;
};

void FastSoftmaxTest::FillInputBuf() {
  for (int idx = 0; idx < _num_classes; idx++) {
    _input[idx] = ((float)rng_.Rand31() - (1 << 30)) / (1u << 30);
  }
}

void FastSoftmaxTest::RunSoftmaxTest() {
  av1_nn_softmax(_input, _ref_buf, _num_classes);
  _target_fn(_input, _dst_buf);

  for (int idx = 0; idx < _num_classes; idx++) {
    if (_ref_buf[idx] < abs_epsilon) {
      ASSERT_LE(_dst_buf[idx], abs_epsilon)
          << "Reference output was near-zero, test output was not" << std::endl;
    } else {
      const float error = _dst_buf[idx] - _ref_buf[idx];
      const float relative_error = fabsf(error / _ref_buf[idx]);
      ASSERT_LE(relative_error, rel_epsilon)
          << "Excessive relative error between reference and test output"
          << std::endl;
      ASSERT_LE(error, abs_epsilon)
          << "Excessive absolute error between reference and test output"
          << std::endl;
    }
  }
}

void FastSoftmaxTest::RunSoftmaxSpeedTest(const int run_times) {
  aom_usec_timer timer;
  aom_usec_timer_start(&timer);
  for (int idx = 0; idx < run_times; idx++) {
    _target_fn(_input, _dst_buf);
  }
  aom_usec_timer_mark(&timer);
  const int64_t time = aom_usec_timer_elapsed(&timer);
  std::cout << "Test with " << _num_classes << " classes took " << time
            << " us." << std::endl;
}

TEST_P(FastSoftmaxTest, RandomValues) {
  FillInputBuf();
  RunSoftmaxTest();
}

TEST_P(FastSoftmaxTest, DISABLED_Speed) {
  constexpr int kNumTimes = 1000000;
  RunSoftmaxSpeedTest(kNumTimes);
}

void AnchorSoftmax16Fn(const float *input, float *output) {
  av1_nn_softmax(input, output, 16);
}

const FastSoftmaxTestParams kArrayParams_c[] = {
  { AnchorSoftmax16Fn, 16 }, { av1_nn_fast_softmax_16_c, 16 }
};
INSTANTIATE_TEST_SUITE_P(C, FastSoftmaxTest,
                         ::testing::ValuesIn(kArrayParams_c));

#if HAVE_SSE3 && !CONFIG_EXCLUDE_SIMD_MISMATCH
INSTANTIATE_TEST_SUITE_P(
    SSE3, FastSoftmaxTest,
    ::testing::Values(FastSoftmaxTestParams(av1_nn_fast_softmax_16_sse3, 16)));
#endif
}  // namespace
