/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "av1/common/cnn.h"
#include "config/av1_rtcd.h"

#define SQR(x) ((x) * (x))
#define FLOAT_TOL 1E-6
#define INT_TOL 0

namespace {

class CNNTest : public ::testing::Test {
 protected:
  static void RunCNNTest(int image_width, int image_height, float *input,
                         float *expected, CNN_CONFIG cnn_config, int in_stride,
                         double tolerance, int use_rounding) {
    int out_width, out_height;
    av1_find_cnn_output_size(image_width, image_height, &cnn_config, &out_width,
                             &out_height);
    const int out_size = out_width * out_height;
    float *output = (float *)aom_malloc(sizeof(*output) * out_size);
    const int out_stride = out_width;

    av1_cnn_predict((const float **)&input, image_width, image_height,
                    in_stride, &cnn_config, &output, out_stride);

    if (use_rounding) {
      for (int i = 0; i < out_size; ++i) {
        output[i] = roundf(output[i]);
      }
    }

    double mse = 0;
    for (int i = 0; i < out_size; ++i) {
      EXPECT_LE(fabsf(expected[i] - output[i]), 1)
          << i << ": " << expected[i] << "/" << output[i] << std::endl;
      mse += SQR(expected[i] - output[i]);
    }
    mse /= out_size;
    EXPECT_LE(mse, tolerance);

    aom_free(output);
  }

  static void AssignLayerWeightsBiases(CNN_CONFIG *cnn_config, float *weights,
                                       float *bias) {
    size_t weight_offset = 0;
    size_t bias_offset = 0;
    for (int layer = 0; layer < cnn_config->num_layers; ++layer) {
      CNN_LAYER_CONFIG *layer_config = &cnn_config->layer_config[layer];
      layer_config->weights = weights + weight_offset;
      layer_config->bias = bias + bias_offset;
      weight_offset += layer_config->filter_width *
                       layer_config->filter_height * layer_config->in_channels *
                       layer_config->out_channels;
      bias_offset += layer_config->out_channels;

      ASSERT_NE(layer_config->weights, nullptr);
      ASSERT_NE(layer_config->bias, nullptr);
    }
  }
};

}  // namespace

TEST_F(CNNTest, TestMultilayerConvolution) {
  int image_height = 16;
  int image_width = 16;
  int filter_height = 5;
  int filter_width = 4;

  float input[] = {
    2,  -1, -4, 0,  0,  -4, -4, -5, 3,  -1, 3,  1,  2,  2,  0,  0,  0,  -2, 3,
    4,  -3, -1, -1, 3,  3,  1,  -3, -3, 3,  4,  3,  0,  -4, 0,  4,  -2, -1, -2,
    0,  -5, 3,  -1, -2, -2, 0,  -3, -2, -4, 2,  -2, -1, 3,  -2, -3, 3,  2,  -4,
    -3, 4,  -3, -4, -4, -3, -1, -1, 3,  -1, 0,  -2, -4, 2,  3,  0,  1,  4,  3,
    1,  1,  2,  -1, -2, -5, 1,  -5, 1,  0,  -5, -5, -2, -3, -1, -1, -3, 1,  -4,
    -5, -1, 4,  3,  -1, -1, 2,  -1, -3, 3,  -5, 3,  -4, -5, -5, -5, -3, 2,  0,
    -2, -2, 2,  3,  2,  1,  2,  -4, -1, -1, 2,  -5, 3,  -5, 0,  -4, 1,  3,  -3,
    3,  -4, -4, -4, 2,  0,  4,  4,  4,  -1, 1,  -4, -2, 0,  -2, -4, -2, 2,  0,
    -3, 3,  -1, -4, -4, 1,  3,  -3, -3, -2, -3, -5, -4, -1, -3, -4, -3, -4, 4,
    2,  3,  2,  -4, 3,  0,  -1, -4, 0,  -1, 2,  -1, 3,  1,  -1, -5, -3, -2, 3,
    3,  4,  -2, -3, 3,  -1, -4, -5, 4,  -4, 1,  1,  4,  2,  -2, -3, -4, -1, 1,
    -4, 2,  -5, 0,  -5, 1,  4,  -3, -5, -4, 1,  2,  4,  -5, 4,  1,  3,  2,  3,
    -3, 1,  -5, 2,  -3, 2,  -1, 4,  0,  -2, 1,  -2, 4,  -2, -3, 0,  2,  -1, -1,
    3,  -5, 0,  3,  0,  -4, 3,  -5, 3,
  };

  float weights[] = {
    3,  -5, 2,  4,  2,  4,  3,  -3, 3,  1,  0,  1,  3,  2,  1,  4,  -1, 3,  2,
    3,  1,  4,  1,  2,  -1, 3,  0,  4,  3,  4,  -2, -3, 2,  -3, -5, 4,  1,  2,
    -1, 4,  -4, -3, -4, -5, 2,  3,  4,  2,  -4, -2, 2,  4,  -1, 0,  -4, 3,  1,
    -3, -5, 1,  2,  -1, -4, -1, -2, -1, -3, 2,  -3, -1, -1, -5, 4,  1,  3,  1,
    4,  4,  4,  3,  3,  0,  3,  -1, -4, -1, 4,  1,  -2, 1,  -5, 1,  -5, 3,  1,
    2,  -2, 3,  -4, -4, -2, 4,  -4, 2,  4,  -3, -5, -5, -3, -2, -5, 2,  0,  3,
    -5, -3, 1,  -4, -2, 0,  -5, 2,  4,  0,  -2, -1, -4, 3,  0,  3,  2,  1,  1,
    3,  -4, -1, -2, -5, -4, -4, -1, -4, -1, 3,  -4, 3,  -3, -5, -4, -5, 3,  -1,
    -1, 0,  2,  0,  -1, -1, -3, 4,  -2, -4, -3, -4, -3, -5, 1,  0,  -3, -3, -3,
    -5, 3,  3,  -2, 2,  -3, 4,  1,  -2, -2, -2, -2, -4, 4,  -2, 1,  3,  -5, -5,
    4,  -3, -5, -1, -4, 3,  3,  -5, 3,  -3, -2, -1, -5, 2,  -5, 1,  -5, 1,  2,
    4,  -3, 4,  -5, -5, -5, -4, 0,  1,  4,  3,  3,  1,  4,  -2, 1,  2,  0,  0,
    2,  2,  3,  -2, -2, -3, 4,  -1, 4,  -4, -5, 2,  -5, 2,  -4, 4,  -1, -5, -4,
    -1, 0,  3,  2,  -4, 4,  -2, 1,  2,  2,  -2, 3,  -5, 2,  -2, -2, -4, -1, -1,
    -4, -3, 1,  -3, -5, 4,  -5, -4, 1,  -2, -3, 1,  3,  1,  4,  -4, 2,  1,  3,
    0,  2,  3,  -3, 4,  2,  4,  -5, 3,  3,  -4, 1,  -1, 3,  1,
  };

  float bias[] = {
    2, -1, 4, 2, -3, -1, 4,
  };

  float expected_same[] = {
    2412,   -1353,  -36855, -13267, -1353,  -7223,  -28091, -10749, -1390,
    14674,  8070,   -7313,  5838,   18967,  19163,  12784,  -5344,  -7332,
    15894,  6061,   7721,   8502,   -4010,  -312,   31069,  45410,  15857,
    34720,  35293,  38747,  22561,  7262,   -11152, -11603, 4140,   -27075,
    -37735, -26126, -30887, 3408,   776,    4609,   7914,   20675,  15715,
    2733,   -5979,  2130,   31779,  30184,  -13805, 6444,   -1153,  -12613,
    14648,  -4999,  -4238,  -19152, -18056, 15026,  9023,   -7714,  15410,
    8405,   -1379,  -12796, 8300,   25404,  8627,   12156,  3264,   -9059,
    -13648, 34210,  24965,  7790,   -992,   11812,  18209,  4691,   -9101,
    -3233,  -15160, -23722, -38272, -32806, -41914, -22635, -6012,  19254,
    -20767, -22402, -10091, -17617, -28073, -16045, 2587,   7860,   -21220,
    -14177, 10389,  -4491,  248,    19575,  15685,  -1359,  -17640, -3632,
    -13841, -49222, -32822, -26854, 5381,   4693,   -8310,  9519,   3587,
    -16949, -24650, -27537, -26269, -29909, -30116, -45298, -36634, -45915,
    -18383, -6127,  24443,  24706,  4035,   23397,  22435,  -15255, 11452,
    -5270,  -27475, -260,   -42195, -7682,  -37525, -22571, -12201, -33152,
    5887,   17927,  40896,  13472,  26607,  10167,  -46,    -25187, -29185,
    -18555, 16337,  -16873, -18328, -19878, -32956, -10003, -17438, 5109,
    -29156, -63238, -20460, -27501, -20033, -22191, -336,   21400,  4625,
    14473,  27766,  13775,  30088,  -4718,  206,    -15532, -46812, -7503,
    -22897, 9217,   -24039, -28894, -2080,  -13868, 20672,  17776,  9825,
    -13053, 13383,  4378,   -58379, -49322, -41920, -38439, -31943, -46785,
    -18418, -8618,  -11034, -18459, -12504, -16734, 13054,  14207,  12332,
    -4743,  -25130, -15763, -10218, -29675, -38421, -16386, -12881, -32016,
    -12225, -24790, -30425, 8730,   24811,  173,    17422,  -450,   -3268,
    11600,  -36617, -33881, -17016, -14122, -24865, -731,   -1673,  -22517,
    10455,  -12649, -7425,  -16148, -7160,  -3403,  -943,   663,    -28774,
    -10435, -35024, -4955,  -14204, -42191, -12220, -5036,  -14723, 3055,
    -6123,  -11845, 20722,  -13351,
  };

  float expected_valid[] = {
    -14177, 10389, -4491,  248,    19575,  15685,  -1359,
    9519,   3587,  -16949, -24650, -27537, -26269, -29909,
    23397,  22435, -15255, 11452,  -5270,  -27475, -260,
    13472,  26607, 10167,  -46,    -25187, -29185, -18555,
  };

  CNN_CONFIG cnn_config = { .num_layers = 3,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = filter_width,
                                    .filter_height = filter_height,
                                    .out_channels = 3,
                                    .skip_width = 1,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = nullptr,
                                    .bias = nullptr,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                                {
                                    .deconvolve = 0,
                                    .in_channels = 3,
                                    .filter_width = filter_width,
                                    .filter_height = filter_height,
                                    .out_channels = 3,
                                    .skip_width = 1,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = nullptr,
                                    .bias = nullptr,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                                {
                                    .deconvolve = 0,
                                    .in_channels = 3,
                                    .filter_width = filter_width,
                                    .filter_height = filter_height,
                                    .out_channels = 1,
                                    .skip_width = 1,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = nullptr,
                                    .bias = nullptr,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  // Weights and biases need to be specified separately because
  // of the offset.
  AssignLayerWeightsBiases(&cnn_config, weights, bias);

  RunCNNTest(image_width, image_height, input, expected_same, cnn_config,
             image_width, INT_TOL, 1);

  for (int i = 0; i < cnn_config.num_layers; ++i) {
    cnn_config.layer_config[i].pad = PADDING_VALID;
  }

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, INT_TOL, 1);
}

TEST_F(CNNTest, TestRELUSingleLayer) {
  int image_width = 8;
  int image_height = 8;
  int filter_height = 5;
  int filter_width = 4;
  float input[] = {
    2,  -5, 4,  -2, 1, 1,  -3, 0,  2,  -4, -4, -5, 4,  0, 3,  4,
    0,  -4, -1, -5, 4, -2, -3, 3,  -5, 1,  -2, -4, 1,  2, 0,  4,
    1,  0,  4,  -4, 2, 4,  -3, 4,  -2, 3,  4,  -3, -4, 4, -2, 4,
    -3, 3,  -3, 4,  2, -3, 2,  -2, -5, -3, 2,  3,  0,  2, 0,  0,
  };
  float expected_same[] = {
    67, 0,  63, 0,  54, 0,  0, 0,  33, 36, 77, 0, 50, 0,  0,  52,
    11, 63, 71, 0,  43, 0,  0, 18, 0,  0,  89, 0, 22, 5,  0,  31,
    0,  0,  0,  40, 0,  17, 0, 14, 0,  0,  0,  3, 0,  0,  0,  6,
    0,  58, 0,  0,  39, 0,  0, 13, 24, 0,  0,  6, 16, 11, 22, 0,
  };
  float expected_valid[] = {
    63, 71, 0, 43, 0, 0, 89, 0, 22, 5, 0, 0, 40, 0, 17, 0, 0, 3, 0, 0,
  };
  float weights[] = {
    -4, -1, 4, 1, -3, 0, -3, -4, 2, 3, -5, 1, -5, -1, -5, 1, -1, 3, -3, -5,
  };
  float bias[] = { 1 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                .deconvolve = 0,
                                .in_channels = 1,
                                .filter_width = filter_width,
                                .filter_height = filter_height,
                                .out_channels = 1,
                                .skip_width = 1,
                                .skip_height = 1,
                                .maxpool = 0,
                                .weights = weights,
                                .bias = bias,
                                .pad = PADDING_SAME_ZERO,
                                .activation = RELU,
                                .input_copy = 0,
                                .skip_combine = SKIP_NONE,
                            } } };

  // Run multilayer same_zero test
  RunCNNTest(image_width, image_height, input, expected_same, cnn_config,
             image_width, INT_TOL, 1);

  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, INT_TOL, 1);
}

TEST_F(CNNTest, TestVaryingStridesVaryingDimImages) {
  float weights[] = {
    1,  -5, -3, -4, -1, 1,  2,  -3, 2,  2,  -1, 1,  -5, 1,  1,
    -3, -5, 3,  1,  4,  -2, -5, -2, -3, -5, 0,  -1, -5, 2,  -2,
    -2, 1,  -2, -4, 1,  3,  -2, 2,  0,  -3, 2,  -3, -2, -3,
  };
  float bias[] = { 2 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 4,
                                    .filter_height = 11,
                                    .out_channels = 1,
                                    .skip_width = 7,
                                    .skip_height = 6,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  int image_height = 24;
  int image_width = 17;
  float input[] = {
    -1, -3, 4,  4,  -5, 4,  3,  -5, -1, -3, 4,  -4, 2,  -3, 3,  -5, 2,  -1, -5,
    1,  -1, 3,  1,  -3, -3, 4,  0,  2,  -3, -5, -5, -4, 0,  -5, -2, -3, -1, -2,
    2,  -5, 4,  4,  0,  -4, -3, 1,  -3, -5, -4, -4, 1,  -2, -3, 3,  -3, -3, -1,
    -5, -5, -2, 3,  1,  -1, -5, -5, 1,  -4, -2, -1, -2, -4, -4, 2,  -2, 2,  1,
    -2, -4, -1, 1,  -2, -5, 3,  -2, -1, -1, -5, -3, 1,  -2, -2, -3, -1, -2, -4,
    -2, 1,  -4, -1, 4,  3,  -4, 0,  4,  2,  2,  4,  -3, -5, 2,  2,  1,  -1, -4,
    -2, 1,  3,  2,  0,  4,  -1, -3, 2,  1,  -4, 2,  2,  -4, -2, 0,  -2, -1, 4,
    4,  2,  3,  -4, 2,  -4, -5, 4,  -1, -3, -1, 0,  -4, 1,  3,  -1, -3, -5, 3,
    -2, -4, 1,  2,  -2, -3, -3, -5, 1,  -3, -1, 0,  -1, 3,  -4, -1, -5, -5, 1,
    0,  0,  -2, -2, 2,  -2, 0,  0,  2,  0,  -3, 0,  -1, -4, -4, -1, 3,  -4, -4,
    -1, 0,  -5, -3, -2, 4,  -3, -4, -4, 0,  -5, 1,  -2, -3, -3, -4, 4,  3,  4,
    3,  3,  -1, 3,  1,  -3, -2, 3,  3,  0,  2,  -4, -3, 2,  2,  0,  -2, 4,  -2,
    2,  -2, -1, -4, -2, 2,  -4, 3,  -1, 4,  1,  1,  4,  -1, -4, -4, 1,  1,  -2,
    4,  -1, 3,  2,  -3, 4,  3,  1,  4,  0,  -4, 2,  0,  2,  4,  -2, -2, 4,  2,
    -1, -2, 1,  -3, 2,  3,  -5, -3, 4,  4,  2,  -5, -4, -5, -2, -4, 2,  0,  2,
    -5, 4,  -4, -2, -5, 2,  1,  0,  4,  1,  -2, -3, -4, -3, -4, 3,  3,  2,  0,
    -3, 1,  -5, 4,  0,  4,  -1, 3,  -5, -5, -2, -1, -1, 4,  3,  3,  4,  3,  -4,
    4,  -3, -3, -1, -4, -1, -4, -1, -2, 4,  -2, -4, 4,  4,  -3, -4, -1, 1,  2,
    -1, -2, -2, 3,  2,  2,  -3, 0,  -1, 0,  3,  2,  -5, 0,  -4, 0,  0,  2,  -4,
    -1, -1, 0,  -2, 0,  1,  0,  0,  4,  -5, -1, -5, 2,  -1, 0,  2,  -1, 1,  3,
    -3, -5, -2, -3, 4,  -2, -2, -1, -3, -4, -1, -2, -4, 1,  4,  -3, -2, -1, 3,
    -3, -2, 3,  2,  1,  -4, -3, -5, 1,
  };
  float expected_1[] = {
    41, -26, 5, 76, 13, 83, -21, 53, -54, -14, 21, 121,
  };

  RunCNNTest(image_width, image_height, input, expected_1, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 6;
  cnn_config.layer_config[0].skip_height = 7;

  float expected_2[] = {
    21, -50, 41, 20, 72, 127, -21, 103, 62, -37, 83, -3,
  };
  RunCNNTest(image_width, image_height, input, expected_2, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 3;
  cnn_config.layer_config[0].skip_height = 10;

  float expected_3[] = {
    -26, -21, -35, 69, 49,  4,  -51, -43, -56,
    -41, 15,  -44, 40, -62, 63, 38,  27,  47,
  };
  RunCNNTest(image_width, image_height, input, expected_3, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 10;
  cnn_config.layer_config[0].skip_height = 3;

  float expected_4[] = {
    21, 49, 28, 87, 50, 40, 102, 81, 58, 85, 51, 66, 36, 19, -37, -45,
  };

  RunCNNTest(image_width, image_height, input, expected_4, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestMaxPool) {
  int image_width = 8;
  int image_height = 8;
  int stride = 3;
  float input[] = {
    1,  -4, -4, 8, 0, 7, -5, -2, 8, 2, 2, 8,  5,  -1, -1, 9,
    -3, 0,  -2, 0, 6, 3, -4, 8,  7, 8, 7, -1, 4,  -1, 0,  2,
    -5, -2, 8,  5, 5, 4, 2,  7,  4, 6, 2, 8,  8,  -4, -3, -4,
    -3, -1, 2,  3, 3, 6, -5, 8,  9, 5, 0, -2, -1, 6,  5,  7,
  };

  float expected[] = {
    49, 58, 70, 68, 68, 70, 48, 57, 88,
  };

  float weights[] = {
    3, 1, 3, 4, -1, 5, -2, 1, -4,
  };

  float bias[] = {
    -3,
  };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = stride,
                                    .maxpool = 1,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestDeconvolveNonActivationSingleLayerSingleKernel) {
  int image_width = 4;
  int image_height = 7;
  float input[] = {
    9,  6,   181, 9,  218, 30, 80,  108, 68,  216, 70, 128, 179, 228,
    33, 212, 34,  14, 48,  27, 230, 23,  202, 113, 80, 56,  122, 112,
  };

  float expected_1_same[] = {
    15,   -30,  36,   -525,  377, -193, 558, 531,  6,   -24,  -15,  124,
    166,  -561, -356, -754,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    433,  -311, 711,  381,   247, -317, 453, 129,  215, -627, -409, -885,
    17,   -255, -55,  -647,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    133,  -719, 633,  -225,  785, 191,  463, 79,   65,  9,    77,   -853,
    -365, -949, -15,  -667,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    355,  -866, 990,  207,   747, 12,   520, -116, 176, -312, -133, -1370,
    -426, -802, 143,  -771,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    65,   -79,  127,  -59,   135, -90,  195, 114,  31,  -91,  -57,  -133,
    17,   -176, -72,  -276,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    457,  -302, 733,  58,    470, -475, 829, 490,  227, -670, -440, -790,
    153,  -588, -294, -1150, -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
    157,  -251, 349,  -185,  409, -293, 587, 251,  77,  -187, -107, -369,
    7,    -481, -135, -827,  -3,  -3,   -3,  -3,   -3,  -3,   -3,   -3,
  };
  float expected_1_valid[] = {
    -30,  15,   -30,  36,   -525,  377,  -193,  558,  531,  24,   24,   6,
    6,    -24,  -15,  124,  166,   -561, -356,  -754, -21,  -39,  -3,   -3,
    -3,   -3,   -3,   -3,   -3,    -3,   -3,    -3,   -3,   -657, 433,  -311,
    711,  381,  247,  -317, 453,   129,  321,   321,  215,  215,  -627, -409,
    -885, 17,   -255, -55,  -647,  -219, -435,  -3,   -3,   -3,   -3,   -3,
    -3,   -3,   -3,   -3,   -3,    -3,   -207,  133,  -719, 633,  -225, 785,
    191,  463,  79,   381,  381,   65,   65,    9,    77,   -853, -365, -949,
    -15,  -667, -259, -515, -3,    -3,   -3,    -3,   -3,   -3,   -3,   -3,
    -3,   -3,   -3,   -540, 355,   -866, 990,   207,  747,  12,   520,  -116,
    633,  633,  176,  176,  -312,  -133, -1370, -426, -802, 143,  -771, -427,
    -851, -3,   -3,   -3,   -3,    -3,   -3,    -3,   -3,   -3,   -3,   -3,
    -105, 65,   -79,  127,  -59,   135,  -90,   195,  114,  78,   78,   31,
    31,   -91,  -57,  -133, 17,    -176, -72,   -276, -57,  -111, -3,   -3,
    -3,   -3,   -3,   -3,   -3,    -3,   -3,    -3,   -3,   -693, 457,  -302,
    733,  58,   470,  -475, 829,   490,  336,   336,  227,  227,  -670, -440,
    -790, 153,  -588, -294, -1150, -229, -455,  -3,   -3,   -3,   -3,   -3,
    -3,   -3,   -3,   -3,   -3,    -3,   -243,  157,  -251, 349,  -185, 409,
    -293, 587,  251,  333,  333,   77,   77,    -187, -107, -369, 7,    -481,
    -135, -827, -227, -451,
  };
  float weights_1[] = { -3, 2, -1, 3, 3, 1, 1, -3, -2, -4 };
  float bias_1[] = { -3 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                .deconvolve = 1,
                                .in_channels = 1,
                                .filter_width = 5,
                                .filter_height = 2,
                                .out_channels = 1,
                                .skip_width = 2,
                                .skip_height = 3,
                                .maxpool = 0,
                                .weights = weights_1,
                                .bias = bias_1,
                                .pad = PADDING_SAME_ZERO,
                                .activation = NONE,
                                .input_copy = 0,
                                .skip_combine = SKIP_NONE,
                            } } };

  RunCNNTest(image_width, image_height, input, expected_1_same, cnn_config,
             image_width, INT_TOL, 1);

  // Change padding to valid
  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_1_valid, cnn_config,
             image_width, INT_TOL, 1);

  float expected_12_same[] = {
    15,  -12,  6,    36,   -9,   -528, 377,  -184, 513,  558,  -12,  24,
    6,   -30,  -15,  -33,  -21,  166,  154,  -546, -356, -718, -30,  -21,
    433, -221, 561,  711,  -33,  -153, 247,  -83,  -87,  453,  -111, 321,
    215, -657, -409, -845, -93,  17,   -43,  -243, -55,  -215, -327, -219,
    133, -71,  -447, 633,  -219, 435,  785,  -73,  -177, 463,  -131, 381,
    65,  -207, 77,   -59,  -651, -365, -797, -213, -15,  -155, -387, -259,
    355, -182, -150, 990,  -231, 582,  747,  -36,  -540, 520,  -215, 633,
    176, -540, -133, -491, -687, -426, -882, -102, 143,  77,   -639, -427,
    65,  -37,  57,   127,  -17,  -105, 135,  -51,  60,   195,  -30,  78,
    31,  -105, -57,  -125, -45,  17,   -11,  -147, -72,  -168, -84,  -57,
    457, -233, 618,  733,  -26,  -540, 470,  -205, 264,  829,  -116, 336,
    227, -693, -440, -900, -72,  153,  107,  -609, -294, -698, -342, -229,
    157, -83,  69,   349,  -59,  -201, 409,  -125, 27,   587,  -115, 333,
    77,  -243, -107, -267, -171, 7,    -105, -369, -135, -379, -339, -227,
  };
  float expected_12_valid[] = {
    -30,  15,   -12,  6,    36,   -9,   -528, 377,  -184, 513,  558,  -12,
    24,   24,   6,    6,    -30,  -15,  -33,  -21,  166,  154,  -546, -356,
    -718, -30,  -21,  -39,  -657, 433,  -221, 561,  711,  -33,  -153, 247,
    -83,  -87,  453,  -111, 321,  321,  215,  215,  -657, -409, -845, -93,
    17,   -43,  -243, -55,  -215, -327, -219, -435, -207, 133,  -71,  -447,
    633,  -219, 435,  785,  -73,  -177, 463,  -131, 381,  381,  65,   65,
    -207, 77,   -59,  -651, -365, -797, -213, -15,  -155, -387, -259, -515,
    -540, 355,  -182, -150, 990,  -231, 582,  747,  -36,  -540, 520,  -215,
    633,  633,  176,  176,  -540, -133, -491, -687, -426, -882, -102, 143,
    77,   -639, -427, -851, -105, 65,   -37,  57,   127,  -17,  -105, 135,
    -51,  60,   195,  -30,  78,   78,   31,   31,   -105, -57,  -125, -45,
    17,   -11,  -147, -72,  -168, -84,  -57,  -111, -693, 457,  -233, 618,
    733,  -26,  -540, 470,  -205, 264,  829,  -116, 336,  336,  227,  227,
    -693, -440, -900, -72,  153,  107,  -609, -294, -698, -342, -229, -455,
    -243, 157,  -83,  69,   349,  -59,  -201, 409,  -125, 27,   587,  -115,
    333,  333,  77,   77,   -243, -107, -267, -171, 7,    -105, -369, -135,
    -379, -339, -227, -451,
  };

  // Change skip_width, skip_height to {2, 3}
  cnn_config.layer_config[0].skip_width = 3;
  cnn_config.layer_config[0].skip_height = 2;
  // Set padding to same
  cnn_config.layer_config[0].pad = PADDING_SAME_ZERO;

  RunCNNTest(image_width, image_height, input, expected_12_same, cnn_config,
             image_width, INT_TOL, 1);

  // Change padding to valid
  cnn_config.layer_config[0].pad = PADDING_VALID;
  RunCNNTest(image_width, image_height, input, expected_12_valid, cnn_config,
             image_width, INT_TOL, 1);

  cnn_config.layer_config[0].filter_width = 4;
  cnn_config.layer_config[0].filter_height = 3;
  float weights_2[] = { -1, -3, -1, -3, 0, 2, -2, 4, 3, 0, 1, 4 };
  float bias_2[] = { -4 };
  cnn_config.layer_config[0].weights = weights_2;
  cnn_config.layer_config[0].bias = bias_2;

  cnn_config.layer_config[0].skip_width = 5;
  cnn_config.layer_config[0].skip_height = 2;
  float expected_2_same[] = {
    -13,  -31,  -13,  -31,  -4,   -10,  -22,  -10,  -22,  -4,   -185, -547,
    -185, -547, -4,   -13,  -31,  -13,  -31,  -4,   -4,   14,   -22,  32,
    -4,   -4,   8,    -16,  20,   -4,   -4,   358,  -366, 720,  -4,   -4,
    14,   -22,  32,   -4,   -195, -658, -213, -622, -4,   -16,  -94,  -28,
    -70,  -4,   459,  -244, 97,   480,  -4,   -85,  -328, -103, -292, -4,
    -4,   432,  -440, 868,  -4,   -4,   56,   -64,  116,  -4,   -4,   156,
    -164, 316,  -4,   -4,   212,  -220, 428,  -4,   582,  -208, 146,  664,
    -4,   -130, -652, -190, -532, -4,   166,  -214, 6,    106,  -4,   192,
    -388, -24,  44,   -4,   -4,   132,  -140, 268,  -4,   -4,   428,  -436,
    860,  -4,   -4,   136,  -144, 276,  -4,   -4,   252,  -260, 508,  -4,
    21,   -541, -115, -269, -4,   416,  -688, -16,  176,  -4,   173,  -103,
    33,   177,  -4,   168,  -640, -88,  -128, -4,   -4,   354,  -362, 712,
    -4,   -4,   452,  -460, 908,  -4,   -4,   62,   -70,  128,  -4,   -4,
    420,  -428, 844,  -4,   499,  -106, 141,  610,  -4,   666,  -46,  210,
    866,  -4,   47,   -148, -19,  -16,  -4,   605,  -85,  181,  763,  -4,
    -4,   64,   -72,  132,  -4,   -4,   24,   -32,  52,   -4,   -4,   92,
    -100, 188,  -4,   -4,   50,   -58,  104,  -4,   -132, -694, -200, -558,
    -4,   15,   -73,  -13,  -17,  -4,   -62,  -610, -158, -418, -4,   -36,
    -343, -90,  -235, -4,   -4,   456,  -464, 916,  -4,   -4,   42,   -50,
    88,   -4,   -4,   400,  -408, 804,  -4,   -4,   222,  -230, 448,  -4,
    606,  -244, 146,  676,  -4,   9,    -172, -37,  -80,  -4,   480,  -370,
    76,   438,  -4,   223,  -340, -3,   112,  -4,   -4,   156,  -164, 316,
    -4,   -4,   108,  -116, 220,  -4,   -4,   240,  -248, 484,  -4,   -4,
    220,  -228, 444,  -4,
  };
  float expected_2_valid[] = {
    -13,  -31,  -13,  -31,  -4,   -10,  -22,  -10,  -22,  -4,   -185, -547,
    -185, -547, -4,   -13,  -31,  -13,  -31,  -4,   14,   -22,  32,   -4,
    -4,   8,    -16,  20,   -4,   -4,   358,  -366, 720,  -4,   -4,   14,
    -22,  32,   -195, -658, -213, -622, -4,   -16,  -94,  -28,  -70,  -4,
    459,  -244, 97,   480,  -4,   -85,  -328, -103, -292, -4,   432,  -440,
    868,  -4,   -4,   56,   -64,  116,  -4,   -4,   156,  -164, 316,  -4,
    -4,   212,  -220, 428,  582,  -208, 146,  664,  -4,   -130, -652, -190,
    -532, -4,   166,  -214, 6,    106,  -4,   192,  -388, -24,  44,   -4,
    132,  -140, 268,  -4,   -4,   428,  -436, 860,  -4,   -4,   136,  -144,
    276,  -4,   -4,   252,  -260, 508,  21,   -541, -115, -269, -4,   416,
    -688, -16,  176,  -4,   173,  -103, 33,   177,  -4,   168,  -640, -88,
    -128, -4,   354,  -362, 712,  -4,   -4,   452,  -460, 908,  -4,   -4,
    62,   -70,  128,  -4,   -4,   420,  -428, 844,  499,  -106, 141,  610,
    -4,   666,  -46,  210,  866,  -4,   47,   -148, -19,  -16,  -4,   605,
    -85,  181,  763,  -4,   64,   -72,  132,  -4,   -4,   24,   -32,  52,
    -4,   -4,   92,   -100, 188,  -4,   -4,   50,   -58,  104,  -132, -694,
    -200, -558, -4,   15,   -73,  -13,  -17,  -4,   -62,  -610, -158, -418,
    -4,   -36,  -343, -90,  -235, -4,   456,  -464, 916,  -4,   -4,   42,
    -50,  88,   -4,   -4,   400,  -408, 804,  -4,   -4,   222,  -230, 448,
    606,  -244, 146,  676,  -4,   9,    -172, -37,  -80,  -4,   480,  -370,
    76,   438,  -4,   223,  -340, -3,   112,  -4,   156,  -164, 316,  -4,
    -4,   108,  -116, 220,  -4,   -4,   240,  -248, 484,  -4,   -4,   220,
    -228, 444,  236,  -4,   76,   316,  -4,   164,  -4,   52,   220,  -4,
    362,  -4,   118,  484,  -4,   332,  -4,   108,  444,
  };
  // Set padding to same
  cnn_config.layer_config[0].pad = PADDING_SAME_ZERO;

  RunCNNTest(image_width, image_height, input, expected_2_same, cnn_config,
             image_width, INT_TOL, 1);

  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_2_valid, cnn_config,
             image_width, INT_TOL, 1);

  cnn_config.layer_config[0].skip_width = 2;
  cnn_config.layer_config[0].skip_height = 5;
  float expected_21_same[] = {
    -31,  -19,  -49,   -191, -565, -194, -574, -13,  14,   -22,  44,   -16,
    382,  -366, 738,   -22,  -4,   23,   32,   545,  20,   204,  720,  5,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -658, -252, -748, -114, -334, -192, -568, -112,
    432,  -440, 928,   -64,  276,  -164, 532,  -220, -4,   304,  868,  266,
    116,  400,  316,   104,  -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -208, -288, -856, -290,
    -862, -202, -598,  -132, 132,  -140, 700,  -436, 1000, -144, 532,  -260,
    -4,   712,  268,   422,  860,  450,  276,  124,  -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -541, -411, -1225, -265, -787, -249, -739, -216, 354,  -362, 1168, -460,
    974,  -70,  552,   -428, -4,   859,  712,  323,  908,  665,  128,  208,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -106, -52,  -148, -66,  -190, -79,  -229, -31,
    64,   -72,  160,   -32,  148,  -100, 242,  -58,  -4,   72,   132,  154,
    52,   125,  188,   23,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -694, -257, -763, -229,
    -679, -319, -949,  -117, 456,  -464, 962,  -50,  492,  -408, 1030, -230,
    -4,   295,  916,   625,  88,   537,  804,  109,  -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -244, -140, -412,  -182, -538, -238, -706, -116, 156,  -164, 428,  -116,
    464,  -248, 708,   -228, -4,   244,  316,  418,  220,  454,  484,  108,
    -4,   -4,   -4,    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,
    -4,   -4,   -4,    -4,
  };
  float expected_21_valid[] = {
    -13,  -31,  -19,  -49,  -191, -565, -194, -574, -13,  -31,   -4,   14,
    -22,  44,   -16,  382,  -366, 738,  -22,  32,   23,   -4,    23,   32,
    545,  20,   204,  720,  5,    32,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -222, -658, -252, -748, -114, -334, -192, -568,  -112, -328,
    -4,   432,  -440, 928,  -64,  276,  -164, 532,  -220, 428,   650,  -4,
    304,  868,  266,  116,  400,  316,  104,  428,  -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -72,  -208, -288, -856, -290, -862,  -202, -598,
    -132, -388, -4,   132,  -140, 700,  -436, 1000, -144, 532,   -260, 508,
    200,  -4,   712,  268,  422,  860,  450,  276,  124,  508,   -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -183, -541, -411, -1225, -265, -787,
    -249, -739, -216, -640, -4,   354,  -362, 1168, -460, 974,   -70,  552,
    -428, 844,  533,  -4,   859,  712,  323,  908,  665,  128,   208,  844,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -38,  -106,  -52,  -148,
    -66,  -190, -79,  -229, -31,  -85,  -4,   64,   -72,  160,   -32,  148,
    -100, 242,  -58,  104,  98,   -4,   72,   132,  154,  52,    125,  188,
    23,   104,  -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -234, -694,
    -257, -763, -229, -679, -319, -949, -117, -343, -4,   456,   -464, 962,
    -50,  492,  -408, 1030, -230, 448,  686,  -4,   295,  916,   625,  88,
    537,  804,  109,  448,  -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,   -4,    -4,   -4,
    -84,  -244, -140, -412, -182, -538, -238, -706, -116, -340,  -4,   156,
    -164, 428,  -116, 464,  -248, 708,  -228, 444,  236,  -4,    244,  316,
    418,  220,  454,  484,  108,  444,
  };

  cnn_config.layer_config[0].pad = PADDING_SAME_ZERO;

  RunCNNTest(image_width, image_height, input, expected_21_same, cnn_config,
             image_width, INT_TOL, 1);

  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_21_valid, cnn_config,
             image_width, INT_TOL, 1);
}

TEST_F(CNNTest, TestSingleLayerEvenStrideEvenDimImage) {
  int image_width = 4;
  int image_height = 4;
  int stride = 2;
  float input[] = {
    -4, 4, -2, 8, -5, 9, 4, -2, 8, 7, 8, 4, -1, -5, 2, 6,
  };
  float expected_horizontal[] = {
    57, -42, 20, 14, -85, -52, -23, -60,
  };
  float expected_vertical[] = {
    13, 20, -81, 14, 7, -23, -21, -60,
  };
  float expected_ver_hor[] = {
    20,
    14,
    -23,
    -60,
  };
  float weights[] = {
    -4, 2, -3, -5, -4, -2, -2, 5, -1,
  };
  float bias[] = { -2 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected_horizontal, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_vertical, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = stride;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_ver_hor, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestSingleLayerEvenStrideOddDimImage) {
  int image_width = 5;
  int image_height = 5;
  int stride = 2;
  float input[] = {
    8, 9, 0,  5, 9,  7,  -4, 7, 0, 1, 2, 9,  6,
    9, 4, -5, 8, -2, -1, 2,  5, 3, 6, 8, -5,
  };
  float expected_horizontal[] = {
    1, -27, -59, -27, 6, 42, -136, -42, -50, 54, -15, -43, 2, -109, -5,
  };
  float expected_vertical[] = {
    1, -99, -27, -16, -59, -136, -39, -42, -142, -50, 2, -45, -109, -31, -5,
  };
  float expected_ver_hor[] = {
    1, -27, -59, -136, -42, -50, 2, -109, -5,
  };
  float weights[] = {
    -5, -1, 4, -4, -5, -3, 6, 7, -5,
  };
  float bias[] = { -1 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected_horizontal, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_vertical, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = stride;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_ver_hor, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestSingleLayerOddStrideEvenDimImage) {
  int image_width = 8;
  int image_height = 8;
  int stride = 3;
  float input[] = {
    1, 3, 1,  -3, 1,  8,  7,  -1, 4,  2,  -1, -2, 0,  -4, -5, -3,
    5, 0, -3, -1, -1, -3, 0,  7,  -4, 1,  8,  2,  1,  0,  1,  6,
    5, 1, 2,  8,  0,  0,  -2, 2,  -1, -3, -1, 8,  6,  6,  -3, 4,
    5, 9, 9,  2,  -3, -4, 2,  7,  1,  6,  0,  -4, -5, 1,  1,  3,
  };
  float expected_horizontal[] = {
    22, -28, -49, 46,  86, -20, -53, -11, 18, -31, 12,  23,
    38, -16, 1,   -15, 0,  25,  -39, 36,  32, 104, -56, 22,
  };
  float expected_vertical[] = {
    19,  46, -15, -20, 86, 67,  -31, -20, 26,  38, -44, 2,
    -16, 40, 44,  1,   44, 104, 30,  -16, -56, 17, 36,  22,
  };
  float expected_ver_hor[] = {
    46, 86, -20, 38, -16, 1, 104, -56, 22,
  };
  float weights[] = {
    0, 2, 7, -3, 5, -5, 5, -2, -3,
  };
  float bias[] = { -4 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected_horizontal, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_vertical, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = stride;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_ver_hor, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestSingleLayerOddStrideOddDimImage) {
  int image_width = 7;
  int image_height = 7;
  int stride = 3;
  float input[] = {
    -5, -3, -1, 1, 4,  4, 6,  4,  0,  0,  3,  2, -5, 3, -5, 2, 7,
    8,  -4, 7,  5, 9,  2, -5, -3, -2, -2, -5, 9, 8,  7, 9,  7, -1,
    -3, 9,  1,  4, -3, 1, -2, 9,  5,  8,  -1, 4, 2,  7, 9,
  };
  float expected_horizontal[] = {
    -2, 4,  -39, -14, -2,  5,  19,  -37, 7,   -102, 14,
    57, -3, 26,  -37, -72, 35, -65, 50,  -12, 2,
  };
  float expected_vertical[] = {
    -2,  9, -4, 4,  28,  43, -39, -102, -53, -7, 14,
    -16, 8, 57, 50, -19, 79, -12, 50,   -3,  2,
  };
  float expected_ver_hor[] = {
    -2, 4, -39, -102, 14, 57, 50, -12, 2,
  };
  float weights[] = {
    2, 3, -3, 3, -5, 6, 3, -3, -3,
  };
  float bias[] = { 3 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_SAME_ZERO,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected_horizontal, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_vertical, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = stride;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_ver_hor, cnn_config,
             image_width, INT_TOL, 0);
}

TEST_F(CNNTest, TestValidPadding) {
  int image_width = 8;
  int image_height = 8;
  int stride = 2;
  float input[] = {
    -4, 9,  8,  0,  4,  5, -1, -4, 6,  9,  -5, 9,  6,  6, 7, 2,
    -1, 0,  -2, -5, 8,  0, 2,  7,  4,  2,  1,  7,  6,  1, 1, 3,
    3,  -1, -1, 7,  -2, 8, 4,  1,  -4, 4,  -1, -3, 5,  3, 5, -2,
    -5, 8,  3,  -1, 3,  0, 2,  5,  -5, -3, 8,  6,  -3, 9, 7, 0,
  };

  float expected_hor_stride[] = {
    11, 50, 68, -1, -88, 65, 29, 37, 16, -2, 2, 48, -33, -29, 6, -56, 50, -29,
  };

  float expected_ver_stride[] = {
    11, 63, 50, 0, 68, 25, 29, -15, 37, -10, 16, 36, -33, 45, -29, 60, 6, 27,
  };

  float expected_hor_ver_stride[] = {
    11, 50, 68, 29, 37, 16, -33, -29, 6,
  };

  float weights[] = {
    5, -3, 2, 2, 3, 1, 4, -3, -4,
  };

  float bias[] = {
    4,
  };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            {
                                {
                                    .deconvolve = 0,
                                    .in_channels = 1,
                                    .filter_width = 3,
                                    .filter_height = 3,
                                    .out_channels = 1,
                                    .skip_width = stride,
                                    .skip_height = 1,
                                    .maxpool = 0,
                                    .weights = weights,
                                    .bias = bias,
                                    .pad = PADDING_VALID,
                                    .activation = NONE,
                                    .input_copy = 0,
                                    .skip_combine = SKIP_NONE,
                                },
                            } };

  RunCNNTest(image_width, image_height, input, expected_hor_stride, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_ver_stride, cnn_config,
             image_width, INT_TOL, 0);

  cnn_config.layer_config[0].skip_width = stride;
  cnn_config.layer_config[0].skip_height = stride;

  RunCNNTest(image_width, image_height, input, expected_hor_ver_stride,
             cnn_config, image_width, INT_TOL, 0);
}
