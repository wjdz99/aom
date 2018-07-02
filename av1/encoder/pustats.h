/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AV1_ENCODER_PUSTATS_H_
#define AV1_ENCODER_PUSTATS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/ml.h"

#define NUM_FEATURES 12
#define NUM_HIDDEN_LAYERS 2
#define HIDDEN_LAYERS_0_NODES 10
#define HIDDEN_LAYERS_1_NODES 10
#define LOGITS_NODES 1

static const float
    av1_pustats_rate_hiddenlayer_0_kernel[NUM_FEATURES *
                                          HIDDEN_LAYERS_0_NODES] = {
      -52.9848f, -3.8058f,  -0.0282f, -0.0054f, 3.5367f,   -0.9193f,
      -0.9668f,  -0.9938f,  -0.9016f, -1.6092f, 5.3714f,   -57.6916f,
      101.6548f, 24.1252f,  0.0247f,  0.0021f,  12.5817f,  0.0826f,
      0.0865f,   0.0863f,   0.0851f,  -0.3443f, -0.0030f,  125.1990f,
      8.0394f,   33.9608f,  -0.1035f, -0.0418f, -13.6075f, 0.0522f,
      0.0570f,   0.0555f,   0.0539f,  -0.1967f, 0.1307f,   48.7033f,
      6.1298f,   -14.1249f, -1.1600f, -0.1019f, -15.8744f, -0.0388f,
      0.0577f,   -0.1058f,  -0.0021f, -0.1957f, 0.0631f,   14.2254f,
      -50.5693f, 14.3733f,  0.4533f,  -0.0918f, -2.1714f,  -0.0448f,
      -0.0633f,  -0.0821f,  -2.7087f, 0.1475f,  0.5090f,   -6.2023f,
      -67.2794f, 7.3298f,   0.0942f,  -0.1726f, -48.4015f, 0.3683f,
      0.3698f,   0.3700f,   0.3686f,  -1.4749f, 0.0488f,   -48.5977f,
      10.8774f,  -2.2353f,  -0.4552f, -0.1260f, -4.4067f,  -0.0928f,
      -0.1082f,  -0.4576f,  -0.2744f, 0.0842f,  0.0779f,   11.0408f,
      -29.7863f, 35.6904f,  0.0205f,  0.0020f,  -0.0657f,  0.1323f,
      0.1330f,   0.1340f,   0.1356f,  -0.5741f, -0.0124f,  -56.8193f,
      -0.1225f,  -0.2490f,  -0.3611f, -0.5065f, 0.0671f,   0.0470f,
      -0.3119f,  0.0031f,   -0.3892f, 0.0209f,  -0.3374f,  0.3125f,
      -84.2314f, 26.4263f,  -0.5345f, -0.0701f, -2.9544f,  -2.8927f,
      0.0248f,   0.0265f,   0.0326f,  -0.2016f, 0.3555f,   -81.5284f,
    };

static const float av1_pustats_rate_hiddenlayer_0_bias[HIDDEN_LAYERS_0_NODES] =
    {
      49.6087f,  -97.7676f, 63.262f,  -60.7973f, 37.9911f,
      114.9542f, -36.5848f, 44.9115f, 0.f,       -19.0169f,
    };

static const float
    av1_pustats_rate_hiddenlayer_1_kernel[HIDDEN_LAYERS_0_NODES *
                                          HIDDEN_LAYERS_1_NODES] = {
      -0.0319f, -0.5428f, -1.0271f, -0.5222f, -0.5177f, -0.1878f, -0.7222f,
      -0.3883f, -0.2759f, -0.3861f, -0.2250f, -0.1830f, 0.0911f,  0.0099f,
      -0.0941f, 0.2999f,  -0.2001f, -0.1672f, -0.4494f, -0.6745f, -0.9414f,
      -0.3108f, -0.4178f, -0.0250f, -1.5476f, -0.7486f, -0.0949f, -1.0776f,
      -0.1162f, -0.0393f, -0.2844f, -1.5018f, 1.0457f,  -0.4873f, 0.0520f,
      0.0325f,  -0.7639f, 0.0901f,  -0.4397f, -0.2217f, -0.3215f, 0.0028f,
      0.1817f,  -0.3343f, -0.0285f, 0.0747f,  -0.4363f, -0.1085f, 0.0262f,
      0.0557f,  -0.0426f, -0.0010f, -0.1310f, -0.0373f, -0.0478f, 0.8759f,
      0.3227f,  -0.3028f, -0.1272f, 0.1649f,  -0.0218f, 0.0020f,  0.4598f,
      -0.5867f, -0.1339f, -0.2585f, -0.5198f, -0.3986f, 0.3847f,  -0.9327f,
      -2.9352f, -1.3549f, 0.0275f,  -0.3291f, 0.0399f,  1.7515f,  -0.4562f,
      -0.5511f, -0.1604f, -5.1808f, -0.5439f, -0.1537f, 0.0408f,  0.0433f,
      -0.1239f, 0.0982f,  -0.1869f, 0.0367f,  0.1234f,  -0.1579f, -6.5508f,
      -4.2091f, 0.4521f,  -0.1081f, 1.1255f,  -0.1830f, 0.1175f,  3.0831f,
      -0.0361f, -7.9087f,
    };

static const float av1_pustats_rate_hiddenlayer_1_bias[HIDDEN_LAYERS_1_NODES] =
    {
      -0.487f,  73.2381f,  -7.2221f,  67.742f,  -94.2112f,
      86.6603f, -78.5546f, -53.7831f, 74.8688f, 62.3852f,
    };

static const float
    av1_pustats_rate_logits_kernel[HIDDEN_LAYERS_1_NODES * LOGITS_NODES] = {
      0.0934f, 5.1931f,  0.0059f,  1.9797f, -7.2527f,
      3.175f,  -3.5956f, -0.7743f, 3.6156f, -0.4359f,
    };

static const float av1_pustats_rate_logits_bias[LOGITS_NODES] = {
  0.0f,
};

static const NN_CONFIG av1_pustats_rate_nnconfig = {
  NUM_FEATURES,                                      // num_inputs
  LOGITS_NODES,                                      // num_outputs
  NUM_HIDDEN_LAYERS,                                 // num_hidden_layers
  { HIDDEN_LAYERS_0_NODES, HIDDEN_LAYERS_1_NODES },  // num_hidden_nodes
  {
      av1_pustats_rate_hiddenlayer_0_kernel,
      av1_pustats_rate_hiddenlayer_1_kernel,
      av1_pustats_rate_logits_kernel,
  },
  {
      av1_pustats_rate_hiddenlayer_0_bias,
      av1_pustats_rate_hiddenlayer_1_bias,
      av1_pustats_rate_logits_bias,
  },
};

static const float
    av1_pustats_dist_hiddenlayer_0_kernel[NUM_FEATURES *
                                          HIDDEN_LAYERS_0_NODES] = {
      -56.6017f, 15.6621f,  3.7010f,  0.0218f,  -0.0018f,  -0.9096f,
      -0.9452f,  -1.1043f,  -1.0458f, -0.4551f, 3.6102f,   -66.0748f,
      4.5685f,   2.2748f,   6.6873f,  0.0012f,  -0.1211f,  0.1272f,
      0.1572f,   0.1432f,   -0.1063f, -0.3305f, -23.2573f, 4.9196f,
      -54.4258f, -31.9244f, 9.2399f,  -0.0801f, -0.0751f,  0.0453f,
      0.0063f,   0.0027f,   0.0427f,  -0.2933f, -0.2374f,  -54.3016f,
      -2.2217f,  -13.5808f, -2.2862f, 0.0104f,  -0.2515f,  -0.8129f,
      -0.5518f,  -0.6280f,  -0.7046f, -0.8152f, -1.0145f,  -1.1242f,
      44.3947f,  10.6982f,  -5.9200f, -0.0189f, 0.1415f,   -0.0944f,
      -0.0716f,  -0.1281f,  -0.0901f, -0.4937f, 0.7867f,   45.5027f,
      -4.5371f,  -2.0040f,  5.6031f,  0.0011f,  -0.5381f,  -0.0243f,
      -0.0368f,  -0.0206f,  0.0020f,  0.0506f,  -0.9719f,  -5.6801f,
      -19.7015f, -5.3514f,  12.4766f, -0.0045f, -0.2037f,  0.0521f,
      0.0282f,   -0.1727f,  -0.2473f, -0.2590f, 0.2071f,   -22.6131f,
      -4.9358f,  -3.9331f,  3.5388f,  0.0008f,  -4.2155f,  0.1522f,
      0.1481f,   0.1467f,   0.1513f,  -0.5977f, 0.0034f,   -10.6549f,
      -30.6004f, -2.1210f,  14.3029f, -0.0014f, 0.2949f,   -0.0915f,
      -0.0850f,  -0.0853f,  -0.0979f, 0.1657f,  0.0993f,   -29.6305f,
      79.3069f,  29.7960f,  -5.0432f, -0.0022f, 0.0707f,   0.0651f,
      0.0940f,   0.1006f,   0.0707f,  -0.3220f, 0.0070f,   90.5679f,
    };

static const float av1_pustats_dist_hiddenlayer_0_bias[HIDDEN_LAYERS_0_NODES] =
    {
      -6.9379f, 5.8412f,  -46.0058f, -8.1814f, 27.7744f,
      21.086f,  22.0872f, 67.3089f,  -2.7895f, 11.0336f,
    };

static const float
    av1_pustats_dist_hiddenlayer_1_kernel[HIDDEN_LAYERS_0_NODES *
                                          HIDDEN_LAYERS_1_NODES] = {
      0.1014f,  0.6381f,  -0.6568f, -0.3215f, -3.4902f, -0.6226f, 0.4069f,
      -0.9190f, -1.1261f, 0.0637f,  0.2580f,  0.2632f,  -0.2919f, -0.1353f,
      -0.6453f, -0.0520f, -0.7564f, 0.2775f,  -0.0796f, 0.5749f,  -0.0057f,
      0.0386f,  -0.0886f, -0.0232f, -0.0206f, -0.5781f, 0.0161f,  -0.0295f,
      -0.0142f, -0.0237f, 0.0169f,  0.2226f,  -0.1543f, -0.4029f, 0.0012f,
      -0.2722f, 0.2627f,  0.4118f,  -0.4138f, -0.0527f, -1.6318f, -0.0647f,
      -1.3425f, -0.8652f, -0.8812f, -1.3089f, -0.1810f, 0.2734f,  -4.9250f,
      0.0292f,  -0.1184f, -0.0521f, -0.2676f, 0.0579f,  -0.1999f, -0.5473f,
      -0.5898f, -0.6735f, -0.3334f, -0.2583f, -0.1262f, -0.1793f, -0.5208f,
      -0.2366f, -0.4227f, 0.6320f,  0.2212f,  -0.1689f, 0.1385f,  -0.5496f,
      -0.2631f, -0.4962f, -0.2349f, -0.7163f, 0.2350f,  -1.2822f, -0.2819f,
      0.8086f,  -0.0526f, 0.0021f,  0.1869f,  -0.1447f, -5.4535f, -0.0396f,
      -0.1659f, -0.6802f, -0.0303f, 0.2985f,  -0.0521f, 0.0192f,  0.1758f,
      0.5731f,  -0.1088f, -0.0085f, -0.6895f, 1.6736f,  -0.8303f, 0.1120f,
      -0.0766f, 0.7975f,
    };

static const float av1_pustats_dist_hiddenlayer_1_bias[HIDDEN_LAYERS_1_NODES] =
    {
      -1.1869f, -8.1852f,  15.5678f, 27.2187f, 7.269f,
      -0.4661f, -17.8186f, 16.6073f, 28.1326f, 44.6699f,
    };

static const float
    av1_pustats_dist_logits_kernel[HIDDEN_LAYERS_1_NODES * LOGITS_NODES] = {
      0.061f,  0.5295f, -0.6127f, -0.1931f, -0.0029f,
      0.0613f, 0.9093f, -0.125f,  -0.386f,  -0.4385f,
    };

static const float av1_pustats_dist_logits_bias[LOGITS_NODES] = {
  67.2089f,
};

static const NN_CONFIG av1_pustats_dist_nnconfig = {
  NUM_FEATURES,                                      // num_inputs
  LOGITS_NODES,                                      // num_outputs
  NUM_HIDDEN_LAYERS,                                 // num_hidden_layers
  { HIDDEN_LAYERS_0_NODES, HIDDEN_LAYERS_1_NODES },  // num_hidden_nodes
  {
      av1_pustats_dist_hiddenlayer_0_kernel,
      av1_pustats_dist_hiddenlayer_1_kernel,
      av1_pustats_dist_logits_kernel,
  },
  {
      av1_pustats_dist_hiddenlayer_0_bias,
      av1_pustats_dist_hiddenlayer_1_bias,
      av1_pustats_dist_logits_bias,
  },
};

#undef NUM_FEATURES
#undef NUM_HIDDEN_LAYERS
#undef HIDDEN_LAYERS_0_NODES
#undef HIDDEN_LAYERS_1_NODES
#undef LOGITS_NODES

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_PUSTATS_H_
