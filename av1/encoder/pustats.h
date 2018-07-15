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

#define NUM_FEATURES 8
#define NUM_HIDDEN_LAYERS 2
#define HIDDEN_LAYERS_0_NODES 12
#define HIDDEN_LAYERS_1_NODES 10
#define LOGITS_NODES 1

static const float
    av1_pustats_rate_hiddenlayer_0_kernel[NUM_FEATURES *
                                          HIDDEN_LAYERS_0_NODES] = {
      -2.2768f,  -0.2548f,  -8.3175f, -3.5508f, 0.5796f,   -0.2358f, 0.4757f,
      -1.0083f,  -0.7177f,  -0.0262f, -0.3647f, -79.6802f, -0.1941f, -0.0910f,
      -0.0927f,  -0.3445f,  -0.4002f, -0.4696f, -0.2362f,  -0.1748f, 0.3036f,
      0.1301f,   -0.1491f,  -0.3097f, 11.2431f, -0.2269f,  3.3231f,  0.7012f,
      0.5869f,   0.7722f,   -0.3415f, 9.2937f,  -2.8292f,  0.0160f,  -1.0942f,
      -35.2705f, 0.4508f,   0.1868f,  0.4050f,  -3.0530f,  -0.6729f, 0.2965f,
      -2.6363f,  -0.6842f,  11.5496f, 11.5213f, -3.9689f,  -2.4595f, -1.1458f,
      0.2378f,   10.6452f,  2.4710f,  -0.2131f, -0.1031f,  -0.3490f, -0.5213f,
      0.3037f,   -1.4847f,  0.1297f,  0.1229f,  -0.5972f,  -0.1977f, -0.3945f,
      0.3176f,   1.5463f,   -0.2800f, 6.8123f,  -1.5447f,  1.3385f,  1.4477f,
      0.3816f,   -0.0203f,  9.6756f,  0.1166f,  0.1628f,   9.2116f,  0.1288f,
      -0.0193f,  -0.0909f,  -3.3687f, -1.6184f, 0.5490f,   -3.8993f, 0.0095f,
      1.2125f,   -45.9498f, 2.8108f,  -3.5735f, -1.2241f,  0.0911f,  0.4030f,
      14.4361f,  0.1872f,   0.0396f,  0.0439f,  7.0700f,
    };

static const float av1_pustats_rate_hiddenlayer_0_bias[HIDDEN_LAYERS_0_NODES] =
    {
      12.2774f, 2.593f, 0.f,      -3.4538f, 8.4591f, -8.1525f,
      -8.892f,  9.377f, 12.4419f, -6.3565f, 1.955f,  -6.1687f,
    };

static const float
    av1_pustats_rate_hiddenlayer_1_kernel[HIDDEN_LAYERS_0_NODES *
                                          HIDDEN_LAYERS_1_NODES] = {
      -0.4061f,  0.1307f,  0.3658f,   -0.4328f, -0.6501f,  -0.4002f,  -0.2277f,
      -0.5164f,  -0.0212f, -0.1756f,  0.2438f,  -0.1566f,  -1.2301f,  -10.4770f,
      0.3744f,   3.8508f,  -11.2213f, 0.5725f,  3.8866f,   -10.8543f, 4.9122f,
      -0.8831f,  4.1796f,  -0.3155f,  2.8838f,  16.1998f,  0.4025f,   -0.6593f,
      2.6392f,   -2.6271f, -20.6320f, -0.6707f, 2.6860f,   -3.0074f,  -2.6359f,
      -2.4866f,  -2.6242f, 16.6459f,  -0.4243f, 2.9406f,   12.8957f,  -2.2209f,
      -13.2353f, 1.3796f,  -2.6115f,  -7.9907f, -1.6771f,  -22.0804f, 2.6535f,
      -1.1387f,  0.1484f,  3.4947f,   0.4213f,  5.6983f,   -13.3756f, -2.1722f,
      -3.1575f,  -5.5585f, 5.4626f,   3.0043f,  -0.2205f,  -8.2210f,  -0.4475f,
      5.5603f,   -1.1792f, -8.6812f,  -4.1572f, -3.5585f,  -2.0427f,  -1.0540f,
      -4.7350f,  -0.0474f, 6.1540f,   10.3184f, -0.4580f,  4.1306f,   -0.1982f,
      -3.3060f,  0.9139f,  -0.9601f,  5.8850f,  -0.3108f,  4.6061f,   -0.2691f,
      7.2848f,   11.5739f, -0.4012f,  -0.7281f, -0.0458f,  -8.0863f,  -7.8078f,
      0.8343f,   3.5845f,  -0.1278f,  -1.9267f, 0.1300f,   2.1685f,   13.9987f,
      0.3020f,   1.6483f,  3.7326f,   -3.0284f, -20.4863f, 0.4174f,   1.2668f,
      -10.3208f, -3.2476f, -11.9617f, 0.0759f,  5.9590f,   -0.2098f,  3.6788f,
      2.0327f,   -1.0612f, -7.4309f,  -3.2946f, -2.8528f,  6.6583f,   0.7427f,
      -5.1442f,
    };

static const float av1_pustats_rate_hiddenlayer_1_bias[HIDDEN_LAYERS_1_NODES] =
    {
      -0.2841f,  -3.3546f, 12.1878f, 0.769f,  -10.0966f,
      -12.0991f, 18.3405f, 12.0135f, 9.8429f, -12.8483f,
    };

static const float
    av1_pustats_rate_logits_kernel[HIDDEN_LAYERS_1_NODES * LOGITS_NODES] = {
      -0.3018f, -4.4424f, 4.2504f, 12.893f, -6.386f,
      -4.7677f, 3.779f,   4.3855f, 5.2613f, -12.6303f,
    };

static const float av1_pustats_rate_logits_bias[LOGITS_NODES] = {
  5.651f,
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
      -11.5748f, 0.7135f,  -11.5758f, -0.0247f, -2.6574f, -1.3213f, -2.2256f,
      -6.7623f,  2.5226f,  -0.5244f,  -8.3804f, -0.0030f, -0.2456f, -0.1897f,
      0.3078f,   3.3137f,  -0.7183f,  -1.2425f, -3.3156f, 0.0413f,  0.5014f,
      0.4676f,   0.1161f,  -0.4057f,  1.1255f,  0.0936f,  2.7856f,  -1.0967f,
      0.3070f,   0.0744f,  0.2993f,   1.4363f,  -0.2683f, -0.2982f, -1.8316f,
      0.1891f,   -0.2623f, -0.1205f,  -0.1713f, -0.3573f, 0.2467f,  0.8551f,
      13.5079f,  -0.0156f, 0.0786f,   0.1090f,  0.0942f,  4.0608f,  -1.6618f,
      0.1554f,   -1.7668f, -0.0166f,  -0.5786f, -0.3844f, -0.3761f, -0.7653f,
      -2.6572f,  -0.2151f, 0.1682f,   -0.0979f, 5.9708f,  5.3669f,  -2.3301f,
      -0.0091f,  4.9852f,  -0.2900f,  -1.2859f, -0.1975f, 1.1933f,  1.1571f,
      0.7575f,   -4.2284f, 5.5757f,   -0.4273f, -5.6879f, -0.0311f, 0.3892f,
      0.3571f,   0.4396f,  4.1854f,   -1.2371f, 0.0280f,  -0.4446f, -0.0492f,
      -10.7118f, -1.9355f, -0.8766f,  -0.1715f, 0.0515f,  -0.3483f, -0.3994f,
      -0.0988f,  0.0468f,  0.1311f,   -0.4209f, -0.5312f,
    };

static const float av1_pustats_dist_hiddenlayer_0_bias[HIDDEN_LAYERS_0_NODES] =
    {
      0.697f,   4.6334f, 7.6687f, 4.2125f,  1.0119f, -5.6733f,
      -1.8012f, 0.1909f, 1.6774f, -0.8309f, 3.1478f, 0.f,
    };

static const float
    av1_pustats_dist_hiddenlayer_1_kernel[HIDDEN_LAYERS_0_NODES *
                                          HIDDEN_LAYERS_1_NODES] = {
      -0.4904f, -0.4618f,  -2.1756f, -2.6557f, 0.3524f,  1.1167f,  0.0595f,
      3.8899f,  1.7736f,   1.5996f,  5.5467f,  0.3871f,  0.0042f,  -0.5585f,
      0.0685f,  -0.1707f,  -0.2692f, -0.1710f, -0.1626f, -0.4055f, -0.0957f,
      0.0881f,  -0.0040f,  0.1965f,  0.3422f,  -0.0724f, 0.7483f,  -1.1085f,
      -0.2173f, -0.1700f,  0.0096f,  -0.6435f, -0.8050f, -0.4250f, -1.3209f,
      -0.3563f, -1.1345f,  1.2263f,  -3.1651f, -2.5455f, -0.4052f, 1.5089f,
      -0.5873f, 0.5438f,   1.6792f,  0.6535f,  0.7872f,  -0.3771f, -0.6016f,
      -0.1101f, -0.0783f,  -0.1199f, -0.5281f, -0.1894f, -0.3710f, -0.7708f,
      -0.5350f, 0.1006f,   -0.3688f, -0.0301f, -0.0932f, 0.7416f,  0.0688f,
      -0.9029f, 0.4990f,   0.1509f,  -0.1793f, 0.4524f,  0.2629f,  -0.5012f,
      1.0667f,  -0.2882f,  0.0782f,  1.3276f,  0.7715f,  1.2118f,  -3.6182f,
      -0.0967f, 0.7434f,   -0.0317f, -0.5569f, -1.1145f, 0.1442f,  -0.4281f,
      -0.8192f, 0.2492f,   0.5602f,  -0.5694f, 0.5241f,  0.1024f,  0.0391f,
      -4.7480f, -0.6323f,  -0.0046f, -3.0745f, -0.1156f, -2.7569f, 1.0767f,
      -0.5614f, -10.9421f, 0.2740f,  0.7001f,  -0.2277f, -6.3701f, -2.1802f,
      0.4267f,  -5.9207f,  0.3395f,  0.4187f,  -0.4590f, -0.5815f, 0.7861f,
      -4.7339f, -0.5587f,  -0.3380f, 1.3454f,  1.6829f,  -0.1177f, 2.8897f,
      0.0797f,
    };

static const float av1_pustats_dist_hiddenlayer_1_bias[HIDDEN_LAYERS_1_NODES] =
    {
      -2.4977f, -0.1337f, 4.5335f, -0.6401f, -0.2269f,
      4.5798f,  0.5122f,  -5.581f, -8.8993f, -0.9924f,
    };

static const float
    av1_pustats_dist_logits_kernel[HIDDEN_LAYERS_1_NODES * LOGITS_NODES] = {
      -0.2442f, -0.3422f, 0.5221f,  -0.3211f, 0.0805f,
      0.6572f,  -0.5178f, -0.4464f, -0.3945f, -0.2623f,
    };

static const float av1_pustats_dist_logits_bias[LOGITS_NODES] = {
  5.3932f,
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
