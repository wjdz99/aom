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
                         double tolerance, int use_rounding, int verbose) {
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

    // Helpful when debugging floating point inputs.
    if (verbose) {
      printf("exp\tout\tdiff\n");
      for (int i = 0; i < out_size; ++i) {
        printf("%.3f,\t%.3f,\t%.3f\n", expected[i], output[i],
               fabsf(expected[i] - output[i]));
      }
      printf("\n");
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
    -3, 1,  -3, 2,  -2, -2, 2,  -2, 1,  -2, -3, 1,  2,  2,  2,  -2, 0,  1,  -1,
    -3, -1, -1, 1,  0,  -3, 1,  0,  -1, 1,  0,  0,  -3, -3, -3, 0,  2,  1,  -1,
    2,  0,  1,  -3, -1, 2,  2,  1,  -2, 0,  -1, 0,  -2, -2, -1, 1,  0,  0,  0,
    -2, -2, -2, 1,  1,  -2, 1,  1,  -2, -2, 1,  -2, -1, -2, -3, 2,  -3, -1, 1,
    0,  -2, -2, -2, 1,  -2, -2, -1, -1, 2,  2,  2,  -1, 1,  -3, -3, 0,  2,  0,
    2,  1,  -3, -3, 1,  2,  2,  1,  -2, -3, 0,  -3, 0,  -3, -2, 0,  1,  1,  0,
    -3, 2,  -1, 2,  1,  0,  1,  -2, 1,  -1, -1, 2,  0,  -2, -3, 1,  1,  -2, -1,
    -3, -3, -1, 0,  -3, -2, 0,  0,  1,  0,  -3, -2, -1, 1,  0,  2,  1,  0,  -3,
    -2, -3, -3, -1, 0,  -2, 2,  -1, -3, 0,  -1, -1, 2,  0,  -3, -2, -1, 0,  0,
    1,  -2, 1,  2,  1,  2,  2,  -3, 2,  -1, 0,  0,  -1, 0,  2,  2,  -1, 2,  -2,
    1,  1,  -3, -3, 1,  -1, -1, -2, 2,  -2, -2, 2,  -1, -3, 2,  -3, 1,  -1, -1,
    -3, 1,  -1, 1,  0,  -3, -3, 1,  -3, -3, 0,  2,  2,  -2, -1, 2,  0,  2,  1,
    -1, -3, 0,  0,  -1, -1, 1,  0,  2,  0,  -3, 2,  1,  0,  1,  -3, 2,  -3, -3,
    -1, -3, -3, 2,  0,  2,  -2, 1,  -1,
  };

  float weights[] = {
    -2, 2,  -2, 2,  -1, -3, 2,  2,  0,  0,  -3, -1, -2, -3, 1,  -1, 0,  0,  0,
    2,  -2, 2,  -2, -3, 1,  1,  1,  -3, -1, 0,  1,  2,  -2, 0,  -1, -3, -1, -2,
    2,  -3, -3, 1,  -2, -3, 0,  2,  1,  -3, -3, -1, -3, -2, -1, -3, -1, -3, -2,
    -1, -3, -1, -2, -2, -3, 2,  0,  -3, 0,  -3, -3, 1,  -3, -1, 0,  -1, 1,  1,
    -1, 1,  -2, 0,  2,  0,  -3, 1,  -1, -1, 2,  0,  1,  -3, -3, 1,  2,  -3, -3,
    1,  -3, 2,  0,  -3, 1,  2,  2,  -2, -1, -2, 1,  1,  0,  -2, -2, 1,  2,  -1,
    -3, 1,  -2, 2,  -3, -2, -3, 2,  1,  0,  -2, 0,  1,  -3, 2,  -2, -2, 0,  2,
    -3, 2,  0,  0,  1,  -2, 1,  1,  -2, -1, -2, 1,  -2, 0,  -2, -2, 0,  -1, -1,
    -3, -3, -3, 1,  -3, -2, 2,  -1, 2,  0,  2,  -2, 2,  -2, 1,  -3, -3, -1, 0,
    2,  2,  1,  -1, -3, -1, -3, 2,  1,  -2, 0,  -3, -1, -3, -1, 2,  1,  0,  2,
    -1, 1,  0,  1,  2,  -1, -2, 2,  1,  -3, -1, -3, 0,  1,  -2, 0,  -2, -3, 0,
    -2, 2,  2,  0,  0,  2,  -3, 2,  -3, -2, 1,  2,  -3, -3, -1, -3, 0,  -3, -3,
    -2, -2, -2, 0,  0,  1,  0,  0,  -1, 0,  0,  -3, 0,  -3, -1, -2, 1,  -2, -1,
    2,  -2, 0,  0,  1,  0,  -2, -1, 0,  -3, 1,  0,  -1, -3, 1,  -1, 1,  -1, -3,
    1,  0,  1,  1,  -1, 2,  2,  0,  0,  1,  -3, 2,  -2, -2, -3, -2, -1, -2, 2,
    0,  2,  -2, -3, -1, -3, 2,  2,  -1, 2,  2,  -1, 0,  -3, 1,
  };

  float bias[] = {
    1, -1, 0, 1, 1, 1, -2,
  };

  float expected_same[] = {
    -1125, 2926,  6406,  631,   -1244, 97,    -1454, 2526,  1065,  3292,  3464,
    2553,  -330,  532,   1038,  1182,  -402,  3758,  3392,  9854,  4365,  1408,
    4736,  3134,  3838,  2409,  3221,  4350,  6750,  4045,  815,   1188,  2959,
    9802,  9590,  4572,  5740,  4253,  1701,  7974,  7012,  6854,  7093,  3907,
    4539,  3886,  4267,  3505,  465,   7824,  9219,  10026, 7968,  957,   2295,
    5594,  10811, 9641,  5950,  10043, 8783,  3132,  1421,  1110,  4108,  13929,
    10660, -84,   -61,   3932,  -180,  6811,  13393, 15147, 15640, 9337,  6961,
    3808,  1604,  1398,  1047,  6739,  10144, 6517,  4698,  2678,  7389,  2595,
    5248,  12075, 11272, 13951, 8820,  1090,  2199,  2206,  2788,  12116, 6683,
    2612,  -291,  3183,  9414,  12316, 14524, 12333, 13208, 7832,  4664,  4657,
    3534,  1298,  -666,  4250,  7707,  9103,  5760,  688,   9571,  15782, 14203,
    14878, 17339, 14684, 8690,  5671,  875,   1429,  1531,  6173,  2984,  5558,
    2996,  7928,  6733,  16117, 15262, 12757, 7980,  3923,  4795,  5973,  2051,
    455,   -1922, 1816,  5906,  3321,  10908, 10910, 7377,  12204, 12809, 11195,
    7451,  6666,  74,    -1645, -35,   -391,  3813,  7324,  892,   1656,  6095,
    12193, 14648, 12156, 14663, 10251, 10325, 7821,  3925,  323,   697,   442,
    1324,  4669,  7002,  5485,  5171,  5086,  10582, 11053, 9709,  11353, 8543,
    5256,  2873,  235,   -628,  1496,  1878,  -867,  3420,  6865,  5937,  10182,
    13277, 10069, 10789, 5998,  624,   -2082, 4417,  1258,  -1080, -819,  -1430,
    1033,  5220,  6335,  8471,  8980,  11908, 14430, 12584, 8404,  1576,  -803,
    985,   1481,  1367,  -193,  873,   3684,  2288,  6676,  9477,  11155, 9602,
    9707,  10507, 4739,  3174,  -575,  -178,  3002,  1710,  423,   -477,  554,
    3088,  2029,  5113,  5000,  3771,  6090,  5365,  1185,  2855,  399,   -312,
    -1577, 176,   955,
  };

  float expected_replicate[] = {
    13768, 13528, 12999, 6906,  4618,  4043,  2611,  9955,  6685,  4776,  2753,
    1036,  3063,  4544,  5183,  7349,  12451, 12501, 9131,  12753, 8908,  4058,
    6299,  7542,  7115,  3307,  3360,  3543,  9754,  7808,  5991,  9019,  14320,
    14919, 12492, 6871,  7373,  3336,  2085,  10604, 9377,  6882,  5009,  3103,
    6220,  6278,  7588,  10196, 11045, 11563, 11842, 11911, 8279,  2030,  1858,
    6368,  12123, 9909,  6347,  10345, 9365,  4038,  1673,  3051,  16492, 16649,
    12276, 408,   -301,  4122,  -654,  7864,  14038, 15279, 15315, 9744,  8243,
    5298,  746,   380,   9824,  9124,  10895, 6640,  4712,  2669,  6980,  2759,
    5385,  12345, 11336, 13129, 8600,  2370,  3682,  5219,  12407, 13123, 6784,
    2612,  -291,  3183,  9414,  12316, 14524, 12333, 13397, 7543,  3916,  4153,
    4477,  4314,  7983,  8418,  9163,  9103,  5760,  688,   9571,  15782, 14203,
    14878, 17718, 14570, 7940,  6642,  5094,  7133,  9964,  10219, 3224,  5558,
    2996,  7928,  6733,  16117, 15262, 12757, 7958,  4401,  5187,  5476,  5529,
    6055,  2206,  3909,  6015,  3321,  10908, 10910, 7377,  12204, 12809, 11195,
    6967,  6840,  481,   -1600, 274,   1,     10373, 8514,  1123,  2117,  6758,
    12736, 16223, 13585, 15988, 11771, 10600, 7918,  4156,  2840,  3111,  3287,
    6359,  7652,  8813,  6530,  6967,  7789,  13671, 13990, 13247, 13241, 9836,
    5251,  3024,  2313,  1834,  4187,  2637,  -1312, 2139,  7378,  7665,  11933,
    15591, 15314, 15678, 9531,  2820,  -1516, 3400,  1314,  22,    363,   -2896,
    -898,  5906,  7308,  10650, 12975, 16978, 20370, 18817, 12381, 4118,  -861,
    -137,  236,   1802,  1632,  -350,  2334,  3400,  8680,  14064, 18216, 18675,
    21765, 22871, 11491, 4937,  -1555, -11,   1669,  2392,  3265,  -5254, -217,
    5001,  8063,  13444, 18884, 19706, 22794, 21064, 9545,  6689,  -7,    289,
    -2021, 504,   2347,
  };

  float expected_valid[] = {
    2612,  -291,  3183,  9414,  12316, 14524, 12333, 9103,  5760,  688,
    9571,  15782, 14203, 14878, 5558,  2996,  7928,  6733,  16117, 15262,
    12757, 3321,  10908, 10910, 7377,  12204, 12809, 11195,
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
             image_width, INT_TOL, 1, 0);

  for (int i = 0; i < cnn_config.num_layers; ++i) {
    cnn_config.layer_config[i].pad = PADDING_SAME_REPLICATE;
  }

  RunCNNTest(image_width, image_height, input, expected_replicate, cnn_config,
             image_width, INT_TOL, 1, 0);

  for (int i = 0; i < cnn_config.num_layers; ++i) {
    cnn_config.layer_config[i].pad = PADDING_VALID;
  }

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, INT_TOL, 1, 0);
}

TEST_F(CNNTest, TestRELUSingleLayer) {
  int image_width = 8;
  int image_height = 8;
  int filter_height = 5;
  int filter_width = 4;
  float input[] = {
    0, -2, -3, 1,  -1, 2,  -2, 1,  -3, -1, 0,  1,  -2, -3, -2, -2,
    1, -3, 2,  -3, -1, -1, 2,  0,  -2, -3, 0,  -2, -3, 1,  -1, -1,
    2, -2, 0,  -2, -3, -3, 1,  1,  -1, 1,  0,  1,  -3, 0,  2,  2,
    0, -3, 1,  -3, 2,  -2, 1,  -1, -1, -2, -3, -2, -1, -3, -2, -1,
  };
  float expected_same[] = {
    9,  0,  1,  1,  0,  3,  0,  19, 0,  12, 10, 0,  0,  0,  5, 0,
    0,  18, 21, 7,  19, 4,  3,  0,  0,  9,  16, 0,  11, 16, 0, 11,
    12, 2,  0,  11, 0,  16, 6,  0,  8,  22, 13, 10, 12, 0,  0, 0,
    0,  1,  2,  12, 29, 6,  10, 0,  13, 0,  0,  5,  8,  10, 0, 0,
  };
  float expected_replicate[] = {
    18, 17, 12, 2,  0,  0,  5,  11, 0,  17, 22, 6,  0,  0,  17, 0,
    0,  18, 21, 7,  19, 4,  3,  5,  3,  9,  16, 0,  11, 16, 0,  3,
    3,  2,  0,  11, 0,  16, 6,  0,  17, 22, 13, 10, 12, 0,  0,  0,
    0,  4,  1,  10, 30, 7,  10, 0,  23, 8,  0,  13, 15, 19, 8,  10,
  };
  float expected_valid[] = {
    18, 21, 7, 19, 4, 9, 16, 0, 11, 16, 2, 0, 11, 0, 16, 22, 13, 10, 12, 0,
  };
  float weights[] = {
    -2, -3, 1, 2, 2, -2, -3, 0, -3, 2, 2, -3, -3, -2, 0, 1, 2, 0, -1, -1,
  };
  float bias[] = { -3 };

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

  RunCNNTest(image_width, image_height, input, expected_same, cnn_config,
             image_width, INT_TOL, 1, 0);

  cnn_config.layer_config[0].pad = PADDING_SAME_REPLICATE;

  RunCNNTest(image_width, image_height, input, expected_replicate, cnn_config,
             image_width, INT_TOL, 1, 0);

  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, INT_TOL, 1, 0);
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
             image_width, INT_TOL, 0, 0);

  cnn_config.layer_config[0].skip_width = 6;
  cnn_config.layer_config[0].skip_height = 7;

  float expected_2[] = {
    21, -50, 41, 20, 72, 127, -21, 103, 62, -37, 83, -3,
  };
  RunCNNTest(image_width, image_height, input, expected_2, cnn_config,
             image_width, INT_TOL, 0, 0);

  cnn_config.layer_config[0].skip_width = 3;
  cnn_config.layer_config[0].skip_height = 10;

  float expected_3[] = {
    -26, -21, -35, 69, 49,  4,  -51, -43, -56,
    -41, 15,  -44, 40, -62, 63, 38,  27,  47,
  };
  RunCNNTest(image_width, image_height, input, expected_3, cnn_config,
             image_width, INT_TOL, 0, 0);

  cnn_config.layer_config[0].skip_width = 10;
  cnn_config.layer_config[0].skip_height = 3;

  float expected_4[] = {
    21, 49, 28, 87, 50, 40, 102, 81, 58, 85, 51, 66, 36, 19, -37, -45,
  };

  RunCNNTest(image_width, image_height, input, expected_4, cnn_config,
             image_width, INT_TOL, 0, 0);
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
             image_width, INT_TOL, 0, 0);
}

TEST_F(CNNTest, TestDeconvolveNonActivationSingleLayerSingleKernel) {
  int image_width = 4;
  int image_height = 4;
  float input[] = {
    160, 80, 118, 226, 117, 171, 134, 77, 120, 117, 76, 137, 175, 187, 56, 137,
  };
  float expected_same[] = {
    -475, -475,  -235,  -547, -349,  -1017, -673,  -221, 1784, 162,  1897,
    930,  1196,  887,   1820, -134,  -666,  -1436, -988, -866, -793, -575,
    -914, -1202, 1547,  602,  1835,  427,   1611,  15,   1716, 733,  -589,
    -826, -922,  -1003, -833, -1135, -828,  -517,  1950, 730,  1741, 1003,
    503,  169,   1886,  553,  -760,  -1284, -1030, -839, -549, -705, -832,
    -817, 1055,  -146,  2002, -444,  1276,  111,   1107, -406,
  };
  float expected_valid[] = {
    -315, 1125, 965,   -235,  329,   431,  379,   997,  1587,  -1125,
    -635, -475, -475,  -235,  -547,  -349, -1017, -673, -221,  5,
    91,   1784, 162,   1897,  930,   1196, 887,   1820, -134,  750,
    -143, -666, -1436, -988,  -866,  -793, -575,  -914, -1202, -447,
    -1,   1547, 602,   1835,  427,   1611, 15,    1716, 733,   -295,
    -241, -589, -826,  -922,  -1003, -833, -1135, -828, -517,  -149,
    -105, 1950, 730,   1741,  1003,  503,  169,   1886, 553,   5,
    -455, -760, -1284, -1030, -839,  -549, -705,  -832, -817,  -269,
    355,  1055, -146,  2002,  -444,  1276, 111,   1107, -406,  690,
    355,  -345, -496,  -719,  -818,  -481, -1,    -381, -680,  -269,
  };
  float weights[] = { -2, 7, 7, -5, -4, -3, -1, 0, 2, 6, -3, 5, 2, -2, -5, -2 };
  float bias[] = { 5 };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                .deconvolve = 1,
                                .in_channels = 1,
                                .filter_width = 4,
                                .filter_height = 4,
                                .out_channels = 1,
                                .skip_width = 2,
                                .skip_height = 2,
                                .maxpool = 0,
                                .weights = weights,
                                .bias = bias,
                                .pad = PADDING_SAME_ZERO,
                                .activation = NONE,
                                .input_copy = 0,
                                .skip_combine = SKIP_NONE,
                            } } };

  RunCNNTest(image_width, image_height, input, expected_same, cnn_config,
             image_width, INT_TOL, 1, 0);

  // Change padding to valid
  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, INT_TOL, 1, 0);
}

TEST_F(CNNTest, TestLargeKernelsAndStrides) {
  float input_10x11[] = {
    4,  4,  2,  4,  2,  -5, -2, 3, -1, 0,  0,  1,  2,  0,  -5, -2, -5, 1,  -3,
    -1, 4,  -3, 2,  -2, 1,  0,  1, -3, -3, -4, -2, -2, 1,  -4, -1, 4,  1,  -4,
    -4, -4, 3,  2,  -5, 3,  -5, 1, 2,  -4, 1,  -1, 3,  4,  -2, 3,  -3, 3,  0,
    2,  -4, -5, -5, -2, -1, -2, 1, 1,  1,  -2, 4,  -5, 4,  -1, -1, 2,  3,  -4,
    2,  2,  3,  0,  0,  1,  0,  3, 2,  3,  1,  -2, 3,  -4, 3,  2,  4,  -2, 0,
    4,  -4, 1,  -3, -3, -3, -5, 1, -3, -5, 0,  4,  -1, -3, 2,
  };

  float weights_10x11[] = {
    -3, 4,  -4, -3, -5, 1,  -2, 3,  1,  -4, -4, 0,  -1, 0,  3,  1,  -3, -2, 0,
    -1, 1,  3,  -4, -4, -3, -3, -2, 4,  3,  -5, 4,  2,  -3, 4,  -2, -1, 2,  -1,
    -5, 0,  -3, 0,  3,  -5, -5, 3,  -4, -1, -5, 3,  4,  0,  4,  -5, 2,  -1, 2,
    -1, -1, -1, -5, 0,  -4, 3,  -1, 1,  1,  -1, 3,  2,  -5, -4, 0,  -4, 4,  -5,
    -3, 4,  -5, 2,  -5, -4, -4, -1, 3,  3,  0,  2,  -4, 1,  -2, 1,  1,  0,  3,
    -2, 0,  1,  2,  4,  -3, -1, -5, -5, 2,  -4, 1,  1,  2,  -4, -2, -2, 2,  1,
    3,  4,  -5, 1,  -1, -3, -3, -1, -2, -5, 1,  -1, 0,  1,  4,  4,  0,  0,  4,
    -3, -1, -5, -3, 0,  1,  1,  1,  -5, 3,  4,  3,  -5, 3,  -2, -2, 0,  -4, 0,
    0,  -2, 1,  -4, -1, 0,  -5, -2, -2, -5, -3, -3, 1,  1,  -3, 2,  4,  2,  4,
    -4, -3, 3,  1,  1,  3,  -4, 4,  -2, -3, -3, -3, -3, -4, -2, 3,  -5, 2,  4,
    -1, -4, -4, 4,  -2, -1, 3,  -3, -4, -4, -2, 4,  1,  0,  2,  -1, 4,  -3, 1,
    4,  -3, 4,  4,  0,  -4, 3,  -2, -3, 2,  3,  -1, -3, 2,  1,  4,  -2, -3, 1,
    4,  -2, 2,  -2, -5, -2, 1,  4,  -1, -4, 4,  -5, 2,  -5, -4, -1, -2, 3,  1,
    2,  1,  -5, 1,  -5, -4, -1, -2, 2,  -2, -4, -3, -2, -2, 4,  -1, 2,  2,  -4,
    2,  -2, 4,  -4, -2, -2, 1,  -1, 1,  1,  1,  -4, -5, -2, 3,  -4, -1, 3,  -2,
    3,  2,  -5, -4, 0,  3,  -2, -4, -5, 3,  -2, -4, 2,  -2, 1,  -4, 0,  2,  -5,
    1,  -4, -1, -1, 4,  -5, -4, 0,  -5, -4, -3, -5, -4, 0,  2,  0,  -4, 2,  -2,
    1,  1,  -3, 2,  0,  -4, 0,  -4, 1,  0,  -5, -1, -1, -1, -5, 4,  2,  2,  -4,
    3,  -2, -2, 2,  -3, -2, -1, 2,  -4, -5, 2,  -2, -4, -5, -5, -1, 2,  -1, 0,
    -5, -2, -2, -5, 0,  1,  -1, -5, 0,  3,  2,  3,  0,  -3, -2, 0,  -5, -1, -2,
    2,  -4, -1, 2,  2,  -5, 2,  -4, 0,  3,  -3, 1,  0,  0,  1,  -5, -3, 1,  -1,
    0,  -4, -3, 2,  -4, -4, 4,  -1, 0,  1,  2,  -4, -5, 4,  -2, 1,  -4, -4, -3,
    -1, -1, 1,  -1, -4, -1, -4, -3, 2,  -1, -2, -4, 1,  1,  0,  -2, 0,  -4, 3,
    -3, 0,  -4, -1, -4, 2,  -1, -2, -5, -1, -2, -3, 3,  -1, 0,  -3, 0,  1,  -5,
    1,  -5, 0,  1,
  };

  float bias_10x11[] = { 3 };

  float expected_10x11[] = {
    118,
  };

  CNN_CONFIG cnn_config = { .num_layers = 1,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                .deconvolve = 0,
                                .in_channels = 1,
                                .filter_width = 23,
                                .filter_height = 20,
                                .out_channels = 1,
                                .skip_width = 15,
                                .skip_height = 20,
                                .maxpool = 0,
                                .weights = weights_10x11,
                                .bias = bias_10x11,
                                .pad = PADDING_SAME_ZERO,
                                .activation = NONE,
                                .input_copy = 0,
                                .skip_combine = SKIP_NONE,
                            } } };

  int image_height = 10;
  int image_width = 11;

  RunCNNTest(image_width, image_height, input_10x11, expected_10x11, cnn_config,
             image_width, INT_TOL, 0, 0);

  float input_11x10[] = {
    -2, -2, 3,  -5, -1, -3, 1,  3,  2,  1,  1,  -5, 4,  1,  3,  -5, 3,  -3, -5,
    0,  -1, -3, -3, 1,  1,  -5, -1, -5, -5, -3, 0,  1,  -3, -1, -3, -3, 0,  3,
    4,  -4, -1, 3,  -3, -1, -3, 1,  -3, -2, -1, -4, -3, 2,  -4, 1,  -4, -1, -3,
    -5, -1, 2,  3,  0,  2,  2,  -5, 4,  1,  2,  -1, -4, 4,  -4, -4, 0,  -1, 1,
    -1, 1,  -3, -3, -2, 1,  2,  4,  4,  4,  -3, -3, 0,  1,  0,  1,  4,  1,  3,
    4,  -3, -2, -4, 4,  2,  0,  3,  4,  -1, 2,  -2, 1,  -3, -2,
  };

  float weights_11x10[] = {
    4,  -1, 1,  -1, 2,  4,  3,  3,  -4, 3,  -5, 1,  -1, -1, -2, -2, 0,  2,  -3,
    -2, 3,  -5, -1, 0,  -1, -2, -2, -1, 2,  4,  3,  1,  0,  0,  -3, 3,  -4, -1,
    -5, 4,  -2, -2, 1,  2,  -1, -3, 1,  2,  -5, 1,  -3, 3,  3,  0,  -4, -4, -5,
    -3, -4, -4, 4,  -2, 4,  4,  -2, 2,  -5, -1, -2, -5, -1, 4,  -3, 3,  -2, 0,
    -4, -3, 0,  -1, -2, 4,  2,  0,  -2, -5, -4, 1,  4,  -4, -2, 2,  -2, 1,  1,
    -4, 1,  -4, -4, -2, 4,  2,  -1, -5, -5, 1,  -3, -3, 3,  -3, -5, -3, 4,  -1,
    -1, -3, 0,  -4, 3,  -1, 0,  -2, 0,  -5, -2, -5, 2,  0,  -5, 2,  3,  -2, 2,
    4,  -1, 1,  -3, 2,  3,  2,  0,  -5, -4, -5, 2,  1,  1,  -1, -2, 3,  4,  2,
    -2, 4,  -2, 3,  1,  -4, -3, -1, 4,  4,  -3, -5, -2, 2,  0,  3,  -2, 3,  -1,
    -4, 0,  -2, 0,  3,  4,  -2, -3, -2, 0,  3,  4,  2,  -4, 0,  1,  2,  2,  -1,
    -1, 4,  1,  4,  -2, -1, -1, -5, 1,  -3, 3,  3,  -1, -4, 3,  -5, 0,  0,  -1,
    -4, -1, -2, 4,  -2, 3,  3,  -3, 1,  -1, 2,  -1, 4,  4,  -2, -2, 4,  -2, 0,
    3,  -3, -5, -1, -2, 4,  -4, 2,  -4, 0,  -2, 3,  -3, 2,  2,  -2, -5, -1, 4,
    3,  -2, -1, 3,  3,  -1, 3,  0,  -3, 0,  4,  2,  0,  -1, 4,  1,  1,  2,  1,
    3,  1,  1,  1,  -3, -5, -4, 4,  -4, 2,  0,  0,  -4, 1,  4,  -5, 4,  4,  0,
    1,  0,  -2, -4, -4, -3, 0,  1,  -5, 4,  0,  -3, -2, -4, 2,  4,  1,  -5, 1,
    -4, 1,  0,  -3, -3, 0,  2,  -5, 4,  3,  -2, -5, 3,  1,  -1, 0,  3,  -2, -2,
    3,  -2, -5, 4,  1,  -2, 2,  -1, 0,  4,  0,  -5, 3,  -2, 1,  2,  1,  -5, -3,
    -2, -5, 4,  -4, 0,  3,  2,  -1, -4, -1, 2,  1,  -2, 3,  -1, -4, 2,  0,  -3,
    1,  -1, 2,  -5, -4, -1, -5, 1,  4,  3,  4,  2,  -3, 1,  -5, -1, 3,  0,  -1,
    -4, 3,  4,  -5, 4,  4,  -3, 2,  -3, -1, -3, -5, -3, 2,  -3, -2, 1,  1,  0,
    -5, 3,  2,  1,  -5, 1,  1,  1,  3,  4,  -4, -1, -2, 0,  -5, -3, -5, -2, -4,
    3,  3,  3,  4,  0,  -4, -1, -5, 0,  -3, 1,  4,  4,  -4, 4,  -5, -5, -1, -2,
    -5, 3,  -4, 4,  3,  0,  -3, 2,  -2, 0,  0,  4,  4,  0,  -2, 1,  -1, -3, 2,
    -1, 1,  -3, -5,
  };

  float bias_11x10[] = {
    -5,
  };

  float expected_11x10[] = {
    36,  -84,  95,   45,  18,   46,   77,  -54, -99,  -149, 66,  49,  161, 11,
    39,  61,   -66,  61,  4,    -3,   34,  -44, -23,  31,   64,  29,  47,  72,
    -27, -27,  121,  -3,  100,  1,    30,  -78, -12,  -89,  -59, 8,   -16, 112,
    91,  -102, -26,  -4,  30,   54,   4,   -84, -24,  -58,  27,  -53, -33, 5,
    53,  -26,  63,   50,  -103, -130, -23, 6,   -104, -207, 73,  23,  77,  132,
    38,  32,   -130, -44, -60,  7,    27,  176, 45,   -32,  -2,  99,  -97, 63,
    69,  126,  47,   63,  136,  -57,  5,   16,  -40,  -157, 8,   38,  -44, -10,
    91,  7,    122,  140, 30,   -105, 4,   -1,  113,  64,   180, 141,
  };

  cnn_config.layer_config[0].weights = weights_11x10;
  cnn_config.layer_config[0].bias = bias_11x10;
  cnn_config.layer_config[0].filter_width = 20;
  cnn_config.layer_config[0].filter_height = 23;
  cnn_config.layer_config[0].skip_width = 1;
  cnn_config.layer_config[0].skip_height = 1;
  image_height = 11;
  image_width = 10;

  RunCNNTest(image_width, image_height, input_11x10, expected_11x10, cnn_config,
             image_width, INT_TOL, 0, 0);
}

TEST_F(CNNTest, TestSoftsignSingleLayer) {
  int image_width = 8;
  int image_height = 8;
  int filter_height = 5;
  int filter_width = 4;
  float input[] = {
    -0.5220f, 0.8410f,  -0.8990f, -0.0090f, 0.6710f,  -0.9470f, -0.8240f,
    -0.0870f, 0.5380f,  0.4750f,  0.570f,   -0.3760f, -0.6960f, -0.5940f,
    -0.3830f, 0.080f,   -0.0980f, -0.4940f, -0.4030f, 0.9460f,  -0.6020f,
    0.4220f,  0.6190f,  0.6640f,  -0.9210f, -0.1470f, -0.2480f, -0.1120f,
    -0.580f,  -0.0650f, 0.3330f,  0.9860f,  -0.7430f, 0.7610f,  0.4840f,
    0.1030f,  0.9570f,  0.6120f,  -0.5240f, -0.1220f, -0.5850f, -0.270f,
    0.7840f,  -0.9790f, 0.7290f,  -0.30f,   -0.6460f, 0.0780f,  0.4750f,
    -0.0510f, 0.4550f,  0.3850f,  -0.7230f, 0.4460f,  -0.6260f, -0.810f,
    0.8720f,  -0.2120f, -0.580f,  -0.9510f, -0.8430f, -0.1340f, -0.0850f,
    0.9190f,
  };
  float expected_same[] = {
    0.430f,   0.660f,  0.5510f,  -0.610f,  0.450f,  -0.1610f, 0.0520f,  0.3240f,
    0.6820f,  0.3820f, 0.6360f,  0.7480f,  0.3080f, 0.090f,   0.3910f,  0.1730f,
    0.340f,   0.6660f, -0.4990f, 0.4280f,  0.1540f, 0.120f,   0.4670f,  0.6150f,
    -0.3880f, 0.7590f, 0.4190f,  0.7350f,  0.5310f, -0.5160f, -0.1760f, 0.6790f,
    -0.6780f, 0.5470f, 0.5750f,  -0.6420f, 0.7210f, -0.4620f, 0.5430f,  0.770f,
    -0.1990f, 0.3950f, 0.7860f,  -0.4380f, 0.7540f, 0.2640f,  -0.6430f, 0.4510f,
    -0.1260f, 0.1590f, -0.2110f, -0.0560f, 0.6570f, 0.680f,   0.5870f,  0.4720f,
    0.4040f,  0.3630f, 0.670f,   0.2360f,  0.410f,  0.6980f,  -0.5350f, 0.3940f,
  };
  float expected_replicate[] = {
    0.540f,   0.7230f,  -0.3530f, -0.2130f, 0.7440f,  -0.4470f, -0.6260f,
    -0.2050f, 0.7230f,  0.4630f,  0.5920f,  0.7440f,  0.6080f,  0.3130f,
    -0.5670f, -0.4720f, 0.5480f,  0.6660f,  -0.4990f, 0.4280f,  0.1540f,
    0.120f,   0.3390f,  0.6090f,  0.4160f,  0.7590f,  0.4190f,  0.7350f,
    0.5310f,  -0.5160f, -0.490f,  0.4450f,  -0.610f,  0.5470f,  0.5750f,
    -0.6420f, 0.7210f,  -0.4620f, 0.3150f,  0.7370f,  -0.5820f, 0.3950f,
    0.7860f,  -0.4380f, 0.7540f,  0.2640f,  -0.7430f, -0.5340f, -0.6270f,
    0.4430f,  0.4730f,  0.4570f,  0.7450f,  0.630f,   0.2620f,  0.3140f,
    -0.1840f, 0.1810f,  0.7210f,  0.2760f,  0.6430f,  0.6720f,  -0.4390f,
    0.2040f,
  };
  float expected_valid[] = {
    0.6660f,  -0.4990f, 0.4280f,  0.1540f,  0.120f,  0.7590f,  0.4190f,
    0.7350f,  0.5310f,  -0.5160f, 0.5470f,  0.5750f, -0.6420f, 0.7210f,
    -0.4620f, 0.3950f,  0.7860f,  -0.4380f, 0.7540f, 0.2640f,
  };
  float weights[] = {
    0.6210f,  0.3710f,  -0.2770f, -0.7230f, -0.2450f, 0.6770f,  0.3080f,
    -0.9880f, -0.080f,  0.7190f,  -0.6760f, -0.0170f, -0.8970f, 0.8260f,
    0.7390f,  -0.4550f, -0.4260f, -0.6330f, 0.0880f,  -0.9390f,
  };
  float bias[] = {
    0.750f,
  };

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
                                .activation = SOFTSIGN,
                                .input_copy = 0,
                                .skip_combine = SKIP_NONE,
                            } } };

  RunCNNTest(image_width, image_height, input, expected_same, cnn_config,
             image_width, FLOAT_TOL, 0, 0);

  cnn_config.layer_config[0].pad = PADDING_SAME_REPLICATE;

  RunCNNTest(image_width, image_height, input, expected_replicate, cnn_config,
             image_width, FLOAT_TOL, 0, 0);

  cnn_config.layer_config[0].pad = PADDING_VALID;

  RunCNNTest(image_width, image_height, input, expected_valid, cnn_config,
             image_width, FLOAT_TOL, 0, 0);
}

TEST_F(CNNTest, TestSkipTensorAdd) {
  int filter_width = 3;
  int filter_height = 3;

  int image_width = 6;
  int image_height = 6;

  float input[] = {
    -2, -2, -3, 2,  0,  -1, -1, 0, 1,  -2, 1, -1, -2, 0, 1, 2,  -2, 2,
    -3, 2,  -1, -3, -3, 0,  0,  2, -2, 0,  2, -1, -1, 0, 2, -1, -2, -2,
  };

  float weights[] = {
    1,  -1, -1, -3, -3, -1, 1,  -3, 1,  -2, 0,  2,  -1, -1, -2, 0,  2,  0,
    -3, -1, -1, -3, -1, 2,  -2, -1, 1,  0,  -3, 1,  -2, 0,  -2, -2, -3, -2,
    -2, -1, -1, 2,  -2, 2,  -1, -2, -1, -1, 2,  0,  -2, 0,  -2, -2, 1,  -2,
    2,  -3, -1, -2, 0,  -2, -1, 1,  2,  0,  0,  -1, 1,  0,  0,  -3, 2,  0,
    2,  0,  0,  -3, 2,  -2, -1, -3, 0,  2,  2,  0,  1,  -3, -2, -2, 1,  -1,
    2,  -2, -2, 1,  0,  -1, 0,  0,  1,  0,  0,  -1, 2,  -3, -3, 0,  0,  -2,
  };

  float bias[] = {
    1, -1, 1, -1, -3, -1, 1,
  };

  float expected[] = {
    -2971, -3879, -4326, -7119, -835,  -521,  -1401, -4918, -9286,
    -6980, -1696, -611,  -157,  -5482, -9842, -7264, -1820, -3158,
    489,   -6367, -7028, -4763, -1629, -1987, 1095,  -3985, -4128,
    -2136, -2574, -4071, 880,   -2678, -1462, -1052, -377,  -1370,
  };

  int channels = 2;

  CNN_CONFIG cnn_config = { .num_layers = 4,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                  .deconvolve = 0,
                                  .in_channels = 1,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = weights,
                                  .bias = bias,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 0,
                                  .skip_combine = SKIP_NONE,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = nullptr,
                                  .bias = nullptr,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 1,
                                  .skip_combine = SKIP_NONE,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = nullptr,
                                  .bias = nullptr,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 0,
                                  .skip_combine = SKIP_ADD,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels,
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
                              } } };

  // Weights and biases need to be specified separately because
  // of the offset.
  AssignLayerWeightsBiases(&cnn_config, weights, bias);

  RunCNNTest(image_width, image_height, input, expected, cnn_config,
             image_width, FLOAT_TOL, 0, 0);
}

TEST_F(CNNTest, TestSkipTensorConcatenation) {
  int filter_width = 3;
  int filter_height = 3;

  int image_width = 6;
  int image_height = 6;

  float input[] = {
    0, 1,  -1, 0,  -2, -3, 0,  2, 1,  1,  -2, -2, 1,  1,  0,  1,  -1, 0,
    2, -3, 2,  -2, 0,  -2, -3, 2, -3, -1, 1,  2,  -2, -3, -3, -2, 1,  -2,
  };

  float weights[] = {
    -2, 0,  2,  -2, 1,  1,  -3, -1, -1, -1, 1,  1,  -1, -3, -1, 0,  0,  0,
    0,  0,  2,  2,  0,  0,  -1, -2, -1, 1,  1,  1,  2,  1,  -1, 2,  -2, 2,
    -1, -3, -2, -2, 0,  -2, 1,  1,  0,  0,  0,  -1, 1,  0,  2,  -3, 1,  -1,
    -1, -1, 1,  -2, -2, -2, -3, 2,  2,  1,  -1, -3, 0,  -2, 1,  -2, -2, -3,
    2,  1,  -3, -2, -2, 1,  1,  2,  2,  -1, -1, -1, 2,  2,  -2, -2, 2,  1,
    1,  -2, 2,  2,  2,  -1, 2,  2,  2,  0,  2,  -1, 0,  -2, -3, 2,  -3, 1,
    1,  -3, 0,  2,  2,  2,  -1, -1, -2, -2, -3, -2, 2,  1,  -3, -3, 1,  2,
  };

  float bias[] = {
    -2, 0, 1, 1, -1, 0, 2,
  };

  float expected[] = {
    4032,  3724,  -21,   -2743, -1811, 684,   5068,  9168,  7067,
    1217,  -4681, -1320, 2373,  6161,  7928,  6649,  -1183, -1995,
    -1976, -2717, 2881,  5825,  2324,  -187,  -5585, -4958, -4390,
    1053,  2734,  987,   881,   -409,  -1173, -1899, 340,   -1015,
  };

  int channels = 2;

  CNN_CONFIG cnn_config = { .num_layers = 4,
                            .is_residue = 0,
                            .ext_width = 0,
                            .ext_height = 0,
                            .strict_bounds = 0,
                            { {
                                  .deconvolve = 0,
                                  .in_channels = 1,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = weights,
                                  .bias = bias,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 0,
                                  .skip_combine = SKIP_NONE,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = nullptr,
                                  .bias = nullptr,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 1,
                                  .skip_combine = SKIP_NONE,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels,
                                  .filter_width = filter_width,
                                  .filter_height = filter_height,
                                  .out_channels = channels,
                                  .skip_width = 1,
                                  .skip_height = 1,
                                  .maxpool = 0,
                                  .weights = nullptr,
                                  .bias = nullptr,
                                  .pad = PADDING_SAME_ZERO,
                                  .activation = NONE,
                                  .input_copy = 0,
                                  .skip_combine = SKIP_CAT,
                              },
                              {
                                  .deconvolve = 0,
                                  .in_channels = channels + channels,
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
                              } } };

  // Weights and biases need to be specified separately because
  // of the offset.
  AssignLayerWeightsBiases(&cnn_config, weights, bias);

  RunCNNTest(image_width, image_height, input, expected, cnn_config,
             image_width, FLOAT_TOL, 0, 0);
}
