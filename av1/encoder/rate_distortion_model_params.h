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

#ifndef AOM_AV1_ENCODER_RATE_DISTORTION_MODEL_PARAMS_H_
#define AOM_AV1_ENCODER_RATE_DISTORTION_MODEL_PARAMS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/ml.h"

// 11 float features +
// 2 categorical features with 4 possible values, converted to one-hot vectors.
// So, total 11 + 2 * 4 = 19 features.
#define NUM_FEATURES 19
#define NUM_HIDDEN_LAYERS 2
#define NUM_HIDDEN_NODES1 8
#define NUM_HIDDEN_NODES2 12
#define NUM_OUTPUTS 1

//------------------------------------------------------------------------------
// RDCost model

static const float av1_rdcost_model_nn_weights_layer0[NUM_FEATURES *
                                                      NUM_HIDDEN_NODES1] = {
  -0.3194f,  0.0788f,   -0.4605f,  0.2537f,   0.2057f,   -0.4373f,  -0.1677f,
  0.3005f,   -0.1586f,  0.4634f,   0.2204f,   -0.1809f,  -0.4144f,  0.0385f,
  -0.1187f,  -0.3390f,  -0.0652f,  -0.3108f,  0.1935f,   1.8068f,   4.8725f,
  20.7078f,  -41.5987f, -2.5585f,  -0.0555f,  101.0366f, 38.0223f,  39.2226f,
  70.6736f,  22.6017f,  38.0128f,  40.3234f,  64.1773f,  19.2249f,  46.6930f,
  2.5733f,   4.9184f,   16.0251f,  -0.2709f,  11.3126f,  29.5035f,  -31.0721f,
  -2.0760f,  -0.0155f,  148.2460f, 34.6513f,  34.0424f,  108.0908f, 28.9256f,
  39.1118f,  41.5133f,  113.1429f, 19.3401f,  54.6402f,  3.4706f,   10.4871f,
  24.2695f,  38.6317f,  54.2707f,  69.8732f,  -58.6257f, -3.1766f,  0.0121f,
  166.2681f, 75.7714f,  76.0826f,  153.4427f, 67.2724f,  73.6636f,  73.3582f,
  138.7363f, 60.3636f,  155.7439f, 38.7893f,  53.5311f,  66.7174f,  -0.1354f,
  0.2755f,   0.3934f,   -0.2216f,  0.4563f,   -0.4256f,  -0.3759f,  -0.0034f,
  0.0993f,   -0.2628f,  0.1446f,   0.3274f,   -0.2631f,  -0.0120f,  0.2464f,
  -0.0764f,  0.4114f,   -0.1336f,  0.2765f,   0.3621f,   -0.2489f,  0.1145f,
  0.4238f,   -0.1868f,  -0.4576f,  -0.3150f,  -0.2988f,  -0.0370f,  0.4271f,
  -0.4276f,  0.4036f,   0.2807f,   -0.0046f,  0.4578f,   0.3065f,   0.1249f,
  -0.1817f,  -0.1324f,  -0.0633f,  0.7054f,   2.5265f,   -12.4316f, 3.8875f,
  0.0580f,   268.6045f, 8.9126f,   6.6741f,   110.9375f, 5.3596f,   9.7548f,
  8.6596f,   105.4006f, 2.9038f,   157.8998f, 0.2947f,   0.2982f,   2.4534f,
  -30.0583f, -34.9159f, -40.2105f, 19.7866f,  -0.2511f,  0.0013f,   340.6193f,
  -29.8182f, -29.8033f, 93.3671f,  -32.0211f, -29.0984f, -29.8956f, 80.6754f,
  -37.1080f, 133.1892f, -28.7440f, -35.2644f, -39.7319f,
};

static const float av1_rdcost_model_nn_biases_layer0[NUM_HIDDEN_NODES1] = {
  0.f, 49.133869f, 57.746731f, 82.337738f, 0.f, 0.f, 5.922832f, -49.669041f,
};

static const float av1_rdcost_model_nn_weights_layer1[NUM_HIDDEN_NODES1 *
                                                      NUM_HIDDEN_NODES2] = {
  0.4616f,  -0.2140f, -0.4262f,  -0.6322f,  0.2519f,   -0.4138f,  -0.2531f,
  -0.4232f, 0.2007f,  -0.5108f,  -1.2897f,  -1.5063f,  0.5018f,   -0.0987f,
  -0.1466f, -1.7344f, 0.4355f,   -3.4744f,  -0.2929f,  -6.0399f,  0.4195f,
  0.0363f,  -4.3033f, -25.5614f, 0.0663f,   -35.1500f, -27.2651f, 0.3588f,
  -0.4110f, 0.4784f,  0.2971f,   22.9029f,  0.4231f,   -0.1536f,  -1.2862f,
  -1.0059f, -0.2625f, 0.1637f,   -0.1976f,  -0.1927f,  0.4812f,   -5.3224f,
  1.1864f,  -4.6968f, -0.3280f,  0.2552f,   -9.8238f,  -39.6380f, -0.3323f,
  -6.1605f, -5.6848f, -36.7243f, 0.5441f,   -0.2701f,  -0.7609f,  24.7130f,
  -0.2551f, -0.5274f, -0.4160f,  -0.0890f,  -0.1762f,  0.2924f,   -0.4012f,
  -0.1150f, 0.0185f,  -0.3240f,  -0.0205f,  -0.1754f,  0.4196f,   0.2582f,
  -0.4914f, -0.5665f, 0.4401f,   -0.2059f,  -0.0437f,  -0.6585f,  0.4588f,
  -0.4420f, -0.2339f, -0.2761f,  -0.5063f,  -2.7243f,  0.8213f,   -8.6805f,
  0.0302f,  -0.3381f, -5.5592f,  -36.4037f, 0.0424f,   -53.5805f, -2.0345f,
  0.3038f,  0.3710f,  -0.3529f,  -35.0580f, 12.2976f,
};

static const float av1_rdcost_model_nn_biases_layer1[NUM_HIDDEN_NODES2] = {
  -0.881873f,  -0.53641f,  -1.430711f, -59.227097f, -0.208408f, -1.803418f,
  -59.931042f, -0.037537f, -0.581957f, -0.235211f,  -1.533446f, -33.803391f,
};

static const float
    av1_rdcost_model_nn_weights_logit_layer[NUM_HIDDEN_NODES2 * NUM_OUTPUTS] = {
      0.28319f,   0.123867f,  -2.284281f, 18.444811f, -0.671599f, -3.94912f,
      48.863983f, -0.395265f, -0.16088f,  -0.253527f, -3.3778f,   14.326732f,
    };

static const float av1_rdcost_model_nn_biases_logit_layer[NUM_OUTPUTS] = {
  49.332661f,
};

static const NN_CONFIG av1_rdcost_model_nnconfig = {
  NUM_FEATURES,
  NUM_OUTPUTS,
  NUM_HIDDEN_LAYERS,
  { NUM_HIDDEN_NODES1, NUM_HIDDEN_NODES2 },
  {
      av1_rdcost_model_nn_weights_layer0,
      av1_rdcost_model_nn_weights_layer1,
      av1_rdcost_model_nn_weights_logit_layer,
  },
  {
      av1_rdcost_model_nn_biases_layer0,
      av1_rdcost_model_nn_biases_layer1,
      av1_rdcost_model_nn_biases_logit_layer,
  },
};

//------------------------------------------------------------------------------

#undef NUM_FEATURES
#undef NUM_HIDDEN_LAYERS
#undef NUM_HIDDEN_NODES
#undef NUM_OUTPUTS

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_RATE_DISTORTION_MODEL_PARAMS_H_
