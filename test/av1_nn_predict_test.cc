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

#include <vector>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "test/function_equivalence_test.h"
#include "test/register_state_check.h"

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom/aom_integer.h"
#include "av1/encoder/ml.h"

using libaom_test::FunctionEquivalenceTest;

namespace {
typedef void (*NnPredict_Func)(const float *const input_nodes,
                               const NN_CONFIG *const nn_config,
                               float *const output);

typedef libaom_test::FuncParam<NnPredict_Func> TestFuncs;

typedef ::testing::tuple<const NnPredict_Func> NnPredictTestParam;

const float epsilon = 1e-3f;  // Error threshold for functional equivalence

class NnPredictTest : public ::testing::TestWithParam<NnPredictTestParam> {
 public:
  virtual void SetUp() {
    const int MAX_NODES2 = NN_MAX_NODES_PER_LAYER * NN_MAX_NODES_PER_LAYER;
    // Allocate two massive buffers on the heap for edge weights and node bias
    // Then set-up the double-dimension arrays pointing into the big buffers
    weights_buf = (float *)aom_malloc(MAX_NODES2 * (NN_MAX_HIDDEN_LAYERS + 1) *
                                      sizeof(*weights_buf));
    bias_buf =
        (float *)aom_malloc(NN_MAX_NODES_PER_LAYER *
                            (NN_MAX_HIDDEN_LAYERS + 1) * sizeof(*bias_buf));
    ASSERT_NE(weights_buf, nullptr);
    ASSERT_NE(bias_buf, nullptr);
    for (int i = 0; i < NN_MAX_HIDDEN_LAYERS + 1; i++) {
      weights[i] = &weights_buf[i * MAX_NODES2];
      bias[i] = &bias_buf[i * NN_MAX_NODES_PER_LAYER];
    }
    target_func_ = GET_PARAM(0);
  }
  virtual void TearDown() {
    aom_free(weights_buf);
    aom_free(bias_buf);
  }
  void runNnPredictTest(const NN_CONFIG *const shape);
  void runNnPredictSpeedTest(const NN_CONFIG *const shape, const int run_times);
  void runNnPredictTest_all(const NN_CONFIG *const shapes,
                            const int num_shapes);
  void runNnPredictSpeedTest_all(const NN_CONFIG *const shapes,
                                 const int num_shapes, const int run_times);
  void predict_original(const float *features, const NN_CONFIG *nn_config,
                        float *output);

 private:
  NnPredict_Func target_func_;
  ACMRandom rng_;
  float *weights[NN_MAX_HIDDEN_LAYERS + 1] = { 0 };
  float *bias[NN_MAX_HIDDEN_LAYERS + 1] = { 0 };
  float *weights_buf = nullptr, *bias_buf = nullptr;
};

void NnPredictTest::runNnPredictTest(const NN_CONFIG *const shape) {
  float inputs[NN_MAX_NODES_PER_LAYER] = { 0 };
  float outputs_test[NN_MAX_NODES_PER_LAYER] = { 0 };
  float outputs_ref[NN_MAX_NODES_PER_LAYER] = { 0 };

  NN_CONFIG nn_config;
  memcpy(&nn_config, shape, sizeof(NN_CONFIG));

  for (int i = 0; i < NN_MAX_HIDDEN_LAYERS + 1; i++) {
    nn_config.weights[i] = weights[i];
    nn_config.bias[i] = bias[i];
  }

  for (int iter = 0; iter < 10000 && !HasFatalFailure(); ++iter) {
    for (int node = 0; node < shape->num_inputs; node++) {
      inputs[node] = (int32_t)rng_.Rand31() - (1 << 30);
      inputs[node] /= (1 << 31);
    }
    for (int layer = 0; layer < shape->num_hidden_layers; layer++) {
      for (int node = 0; node < NN_MAX_NODES_PER_LAYER; node++) {
        bias[layer][node] = (int32_t)rng_.Rand31() - (1 << 30);
        bias[layer][node] /= (1 << 31);
      }
      for (int node = 0; node < NN_MAX_NODES_PER_LAYER * NN_MAX_NODES_PER_LAYER;
           node++) {
        weights[layer][node] = (int32_t)rng_.Rand31() - (1 << 30);
        weights[layer][node] /= (1 << 31);
      }
    }
    // Now the outputs:
    int layer = shape->num_hidden_layers;
    for (int node = 0; node < NN_MAX_NODES_PER_LAYER; node++) {
      bias[layer][node] = (int32_t)rng_.Rand31() - (1 << 30);
      bias[layer][node] /= (1 << 31);
    }
    for (int node = 0; node < NN_MAX_NODES_PER_LAYER * NN_MAX_NODES_PER_LAYER;
         node++) {
      weights[layer][node] = (int32_t)rng_.Rand31() - (1 << 30);
      weights[layer][node] /= (1 << 31);
    }

    av1_nn_predict_c(inputs, &nn_config, outputs_ref);
    target_func_(inputs, &nn_config, outputs_test);

    for (int node = 0; node < shape->num_outputs; node++) {
      if (outputs_ref[node] < epsilon) {
        ASSERT_LE(outputs_test[node], epsilon)
            << "Reference output was near-zero, test output was not";
      } else {
        const float error = outputs_ref[node] - outputs_test[node];
        const float relative_error = fabs(error / outputs_ref[node]);
        ASSERT_LE(relative_error, epsilon)
            << "Excessive relative error between reference and test";
      }
    }
  }
  printf("Passed shape %d", shape->num_inputs);
  for (int layer = 0; layer < shape->num_hidden_layers; layer++)
    printf("x%d", shape->num_hidden_nodes[layer]);
  printf("x%d\n", shape->num_outputs);
}

void NnPredictTest::runNnPredictSpeedTest(const NN_CONFIG *const shape,
                                          const int run_times) {
  float inputs[NN_MAX_NODES_PER_LAYER] = { 0 };
  float outputs_test[NN_MAX_NODES_PER_LAYER] = { 0 };
  float outputs_ref[NN_MAX_NODES_PER_LAYER] = { 0 };

  NN_CONFIG nn_config;
  memcpy(&nn_config, shape, sizeof(NN_CONFIG));

  for (int i = 0; i < NN_MAX_HIDDEN_LAYERS; i++) {
    nn_config.weights[i] = weights[i];
    nn_config.bias[i] = bias[i];
  }
  // Don't bother actually changing the values for inputs/weights/bias: it
  // shouldn't make any difference for a speed test.

  aom_usec_timer timer;
  aom_usec_timer_start(&timer);
  for (int i = 0; i < run_times; ++i) {
    av1_nn_predict_c(inputs, &nn_config, outputs_ref);
  }
  aom_usec_timer_mark(&timer);
  const double time1 = static_cast<double>(aom_usec_timer_elapsed(&timer));
  aom_usec_timer_start(&timer);
  for (int i = 0; i < run_times; ++i) {
    target_func_(inputs, &nn_config, outputs_test);
  }
  aom_usec_timer_mark(&timer);
  const double time2 = static_cast<double>(aom_usec_timer_elapsed(&timer));

  printf("%d", shape->num_inputs);
  for (int layer = 0; layer < shape->num_hidden_layers; layer++)
    printf("x%d", shape->num_hidden_nodes[layer]);
  printf("x%d: ", shape->num_outputs);
  printf("%7.2f/%7.2fns (%3.2f)\n", time1, time2, time1 / time2);
}

static const NN_CONFIG shapes[] = {
  { .num_inputs = 10,
    .num_outputs = 16,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 64 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 12,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 12 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 12,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 24 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 12,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 32 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 18,
    .num_outputs = 4,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 24 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 18,
    .num_outputs = 4,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 32 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 4,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 16 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 8,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 16 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 8,
    .num_outputs = 4,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 16 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 8,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 24 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 8,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 32 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 8,
    .num_outputs = 1,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 64 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 9,
    .num_outputs = 3,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 32 },
    .weights = { 0 },
    .bias = { 0 } },
  { .num_inputs = 4,
    .num_outputs = 4,
    .num_hidden_layers = 1,
    .num_hidden_nodes = { 8 },
    .weights = { 0 },
    .bias = { 0 } },
};

void NnPredictTest::runNnPredictTest_all(const NN_CONFIG *const shapes,
                                         const int num_shapes) {
  for (int i = 0; i < num_shapes; i++) runNnPredictTest(&shapes[i]);
}

void NnPredictTest::runNnPredictSpeedTest_all(const NN_CONFIG *const shapes,
                                              const int num_shapes,
                                              const int run_times) {
  for (int i = 0; i < num_shapes; i++)
    NnPredictTest::runNnPredictSpeedTest(&shapes[i], run_times);
}

void predict_original(const float *features, const NN_CONFIG *nn_config,
                      float *output) {
  int num_input_nodes = nn_config->num_inputs;
  int buf_index = 0;
  float buf[2][NN_MAX_NODES_PER_LAYER];
  const float *input_nodes = features;

  // Propagate hidden layers.
  const int num_layers = nn_config->num_hidden_layers;
  assert(num_layers <= NN_MAX_HIDDEN_LAYERS);
  for (int layer = 0; layer < num_layers; ++layer) {
    const float *weights = nn_config->weights[layer];
    const float *bias = nn_config->bias[layer];
    float *output_nodes = buf[buf_index];
    const int num_output_nodes = nn_config->num_hidden_nodes[layer];
    assert(num_output_nodes < NN_MAX_NODES_PER_LAYER);
    for (int node = 0; node < num_output_nodes; ++node) {
      float val = 0.0f;
      for (int i = 0; i < num_input_nodes; ++i)
        val += weights[i] * input_nodes[i];
      val += bias[node];
      // ReLU as activation function.
      val = val > 0.0f ? val : 0.0f;  // Could use AOMMAX().
      output_nodes[node] = val;
      weights += num_input_nodes;
    }
    num_input_nodes = num_output_nodes;
    input_nodes = output_nodes;
    buf_index = 1 - buf_index;
  }

  // Final output layer.
  const float *weights = nn_config->weights[num_layers];
  for (int node = 0; node < nn_config->num_outputs; ++node) {
    const float *bias = nn_config->bias[num_layers];
    float val = 0.0f;
    for (int i = 0; i < num_input_nodes; ++i)
      val += weights[i] * input_nodes[i];
    output[node] = val + bias[node];
    weights += num_input_nodes;
  }
}

TEST_P(NnPredictTest, RandomValues) {
  runNnPredictTest_all(shapes, sizeof(shapes) / sizeof(NN_CONFIG));
}

TEST_P(NnPredictTest, DISABLED_Speed) {
  runNnPredictSpeedTest_all(shapes, sizeof(shapes) / sizeof(NN_CONFIG),
                            10000000);
}

INSTANTIATE_TEST_CASE_P(C, NnPredictTest, ::testing::Values(predict_original));

#if HAVE_AVX
INSTANTIATE_TEST_CASE_P(AVX2, NnPredictTest,
                        ::testing::Values(av1_nn_predict_avx2));
#endif  // HAVE_AVX

}  // namespace
