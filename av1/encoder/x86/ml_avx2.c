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

#include <stdbool.h>
#include <assert.h>
#include <immintrin.h>
#include "aom_ports/system_state.h"

#include "config/av1_rtcd.h"
#include "av1/encoder/ml.h"

static float nn_propagate_8to1(const float *const inputs,
                               const float *const weights) {
  const __m256 inputs256 = _mm256_loadu_ps(inputs);
  // [i7 i6 i5 i4] [i3 i2 i1 i0]

  const __m256 weights256 = _mm256_loadu_ps(weights);
  // [w7 w6 w5 w4] [w3 w2 w1 w0]

  const __m256 mul = _mm256_mul_ps(inputs256, weights256);
  // [7 6 5 4] [3 2 1 0]

  const __m128 mul_h = _mm256_extractf128_ps(mul, 1);
  // [7 6 5 4]
  const __m128 mul_l = _mm256_castps256_ps128(mul);
  // [3 2 1 0]

  const __m128 vadd = _mm_add_ps(mul_l, mul_h);
  // [7+3 6+2 5+1 4+0]
  const __m128 hadd1 = _mm_hadd_ps(vadd, vadd);
  // [7+3+6+2 5+1+4+0 7+3+6+2 5+1+4+0]
  const __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
  // [7+3+6+2+5+1+4+0 7+3+6+2+5+1+4+0 7+3+6+2+5+1+4+0 7+3+6+2+5+1+4+0]
  return _mm_cvtss_f32(hadd2);
}

static float nn_propagate_4to1(const float *const inputs,
                               const float *const weights) {
  const __m128 inputs128 = _mm_loadu_ps(inputs);
  // [i3 i2 i1 i0]

  const __m128 weights128 = _mm_loadu_ps(weights);
  // [w3 w2 w1 w0]

  const __m128 mul = _mm_mul_ps(inputs128, weights128);
  // [3 2 1 0]

  const __m128 hadd1 = _mm_hadd_ps(mul, mul);
  // [3+2 1+0 3+2 1+0]
  const __m128 hadd2 = _mm_hadd_ps(hadd1, hadd1);
  // [3+2+1+0 3+2+1+0 3+2+1+0 3+2+1+0]
  return _mm_cvtss_f32(hadd2);
}

static void nn_propagate_8to8(const float *const inputs,
                              const float *const weights, __m256 *const outputs,
                              const int num_inputs) {
  const __m256 inputs256 = _mm256_loadu_ps(inputs);
  // [i7 i6 i5 i4] [i3 i2 i1 i0]

  __m256 hadd[4];
  for (int i = 0; i < 4; i++) {
    const __m256 weight0 = _mm256_loadu_ps(&weights[2 * i * num_inputs]);
    // [w7 w6 w5 w4] [w3 w2 w1 w0]
    const __m256 weight1 = _mm256_loadu_ps(&weights[(2 * i + 1) * num_inputs]);
    // [w15 w14 w13 w12] [w11 w10 w9 w8]
    const __m256 mul0 = _mm256_mul_ps(weight0, inputs256);
    const __m256 mul1 = _mm256_mul_ps(weight1, inputs256);
    hadd[i] = _mm256_hadd_ps(mul0, mul1);
    // [15+14 13+12 7+6 5+4] [11+10 9+8 3+2 1+0]
  }
  // hadd[0] = [15+14 13+12  7+6   5+4 ] [11+10  9+8   3+2   1+0 ]
  // hadd[1] = [31+30 29+28 23+22 21+20] [27+26 25+24 19+18 17+16]
  // hadd[2] = [47+46 45+44 39+38 37+36] [43+42 41+40 35+34 33+32]
  // hadd[3] = [63+62 61+60 55+54 53+52] [59+58 57+56 51+50 49+48]

  const __m256 hh0 = _mm256_hadd_ps(hadd[0], hadd[1]);
  const __m256 hh1 = _mm256_hadd_ps(hadd[2], hadd[3]);
  // hh0 = [31+30+29+28 23+22+21+20 15+14+13+12 7+6+5+4]
  //       [27+26+25+24 19+18+17+16 11+10+9+8   3+2+1+0]
  // hh1 = [63+62+61+60 55+54+53+52 47+46+45+44 39+38+37+36]
  //       [59+58+57+56 51+50+49+48 43+42+41+40 35+34+33+32]

  // Permute control: [7 6 5 4][3 2 1 0] => [7 3 6 2][5 1 4 0]
  const __m256i ctrl = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
  const __m256 hh0p = _mm256_permutevar8x32_ps(hh0, ctrl);
  const __m256 hh1p = _mm256_permutevar8x32_ps(hh1, ctrl);
  // hh0p = [31+30+29+28 27+26+25+24 23+22+21+20 19+18+17+16]
  //        [15+14+13+12 11+10+9+8 7+6+5+4 3+2+1+0]
  // hh1p = [63+62+61+60 59+58+57+56 55+54+53+52 51+50+49+48]
  //        [47+46+45+44 43+42+41+40 39+38+37+36 35+34+33+32]

  const __m256 hhh = _mm256_hadd_ps(hh0p, hh1p);
  // [63+62+61+60+59+58+57+56 55+54+53+52+51+50+49+48
  //  31+30+29+28+27+26+25+24 23+22+21+20+19+18+17+16]
  // [47+46+45+44+43+42+41+40 39+38+37+36+35+34+33+32
  //  15+14+13+12+11+10+9+8 7+6+5+4+3+2+1+0]

  // Permute control: [7 6 5 4][3 2 1 0] => [7 6 3 2][5 4 1 0]
  const __m256i ctrl2 = _mm256_set_epi32(7, 6, 3, 2, 5, 4, 1, 0);
  const __m256 hhhp = _mm256_permutevar8x32_ps(hhh, ctrl2);
  // [63+62+61+60+59+58+57+56 55+54+53+52+51+50+49+48
  //  47+46+45+44+43+42+41+40 39+38+37+36+35+34+33+32]
  // [31+30+29+28+27+26+25+24 23+22+21+20+19+18+17+16
  //  15+14+13+12+11+10+9+8 7+6+5+4+3+2+1+0]

  *outputs = _mm256_add_ps(*outputs, hhhp);
}

static void nn_propagate_4to8(const float *const inputs,
                              const float *const weights, __m256 *const outputs,
                              const int num_inputs) {
  const __m128 inputs128 = _mm_loadu_ps(inputs);
  const __m256 inputs256 =
      _mm256_insertf128_ps(_mm256_castps128_ps256(inputs128), inputs128, 1);
  // [i3 i2 i1 i0] [i3 i2 i1 i0]

  __m256 weight0, weight1, weight2, weight3;
  if (num_inputs == 4) {
    weight0 = _mm256_loadu_ps(weights);
    weight1 = _mm256_loadu_ps(&weights[8]);
    weight2 = _mm256_loadu_ps(&weights[16]);
    weight3 = _mm256_loadu_ps(&weights[24]);
  } else {
    weight0 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(&weights[0 * num_inputs])),
        _mm_loadu_ps(&weights[1 * num_inputs]), 1);
    weight1 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(&weights[2 * num_inputs])),
        _mm_loadu_ps(&weights[3 * num_inputs]), 1);
    weight2 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(&weights[4 * num_inputs])),
        _mm_loadu_ps(&weights[5 * num_inputs]), 1);
    weight3 = _mm256_insertf128_ps(
        _mm256_castps128_ps256(_mm_loadu_ps(&weights[6 * num_inputs])),
        _mm_loadu_ps(&weights[7 * num_inputs]), 1);
  }
  // weight0 = [w7 w6 w5 w4] [w3 w2 w1 w0]
  // weight1 = [w15 w14 w13 w12] [w11 w10 w9 w8]
  // weight2 = [w23 w22 w21 w20] [w19 w18 w17 w16]
  // weight3 = [w31 w30 w29 w28] [w27 w26 w25 w24]

  const __m256 mul0 = _mm256_mul_ps(inputs256, weight0);
  // [w7.i3 w6.i2 w5.i1 w4.i0] [w3.i3 w2.i2 w1.i1 w0.i0]
  const __m256 mul1 = _mm256_mul_ps(inputs256, weight1);
  // [w15.i3 w14.i2 w13.i1 w12.i0] [w11.i3 w10.i2 w9.i1 w8.i0]
  const __m256 mul2 = _mm256_mul_ps(inputs256, weight2);
  // [w23.i3 w22.i2 w21.i1 w20.i0] [w19.i3 w18.i2 w17.i1 w16.i0]
  const __m256 mul3 = _mm256_mul_ps(inputs256, weight3);
  // [w31.i3 w30.i2 w29.i1 w28.i0] [w27.i3 w26.i2 w25.i1 w24.i0]

  const __m256 hadd0 = _mm256_hadd_ps(mul0, mul1);
  // [w15.i3+w14.i2 w13.i1+w12.i0 w7.i3+w6.i2 w5.i1+w4.i0]
  // [w11.i3+w10.i2 w9.i1+w8.i0 w3.i3+w2.i2 w1.i1+w0.i0]
  const __m256 hadd1 = _mm256_hadd_ps(mul2, mul3);
  // [w31.i3+w30.i2 w29.i1+w28.i0 w23.i3+w22.i2 w21.i1+w20.i0]
  // [w27.i3+w26.i2 w25.i1+w24.i0 w19.i3+w18.i2 w17.i1+w16.i0]

  const __m256 haddhadd = _mm256_hadd_ps(hadd0, hadd1);
  // [w31.i3+w30.i2+w29.i1+w28.i0 w23.i3+w22.i2+w21.i1+w20.i0
  //    w15.i3+w14.i2+w13.i1+w12.i0 w7.i3+w6.i2+w5.i1+w4.i0]
  // [w27.i3+w26.i2+w25.i1+w24.i0 w19.i3+w18.i2+w17.i1+w16.i0
  //    w11.i3+w10.i2+w9.i1+w8.i0 w3.i3+w2.i2+w1.i1+w0.i0]

  // Permute control: [7 6 5 4][3 2 1 0] => [7 3 6 2][5 1 4 0]
  const __m256i ctrl = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
  const __m256 permute = _mm256_permutevar8x32_ps(haddhadd, ctrl);

  *outputs = _mm256_add_ps(*outputs, permute);
}

static void nn_propagate_8to4(const float *const inputs,
                              const float *const weights, __m128 *const outputs,
                              const int num_inputs) {
  const __m256 inputs256 = _mm256_loadu_ps(inputs);
  // [i7 i6 i5 i4] [i3 i2 i1 i0]

  const __m256 weight0 = _mm256_loadu_ps(&weights[0 * num_inputs]);
  // [w7 w6 w5 w4] [w3 w2 w1 w0]
  const __m256 weight1 = _mm256_loadu_ps(&weights[1 * num_inputs]);
  // [w15 w14 w13 w12] [w11 w10 w9 w8]
  const __m256 weight2 = _mm256_loadu_ps(&weights[2 * num_inputs]);
  // [w23 w22 w21 w20] [w19 w18 w17 w16]
  const __m256 weight3 = _mm256_loadu_ps(&weights[3 * num_inputs]);
  // [w31 w30 w29 w28] [w27 w26 w25 w24]

  const __m256 mul0 = _mm256_mul_ps(inputs256, weight0);
  // [w7.i7 w6.i6 w5.i5 w4.i4] [w3.i3 w2.i2 w1.i1 w0.i0]
  const __m256 mul1 = _mm256_mul_ps(inputs256, weight1);
  // [w15.i7 w14.i6 w13.i5 w12.i4] [w11.i3 w10.i2 w9.i1 w8.i0]
  const __m256 mul2 = _mm256_mul_ps(inputs256, weight2);
  // [w23.i7 w22.i6 w21.i5 w20.i4] [w19.i3 w18.i2 w17.i1 w16.i0]
  const __m256 mul3 = _mm256_mul_ps(inputs256, weight3);
  // [w31.i7 w30.i6 w29.i5 w28.i4] [w27.i3 w26.i2 w25.i1 w24.i0]

  const __m256 hadd0 = _mm256_hadd_ps(mul0, mul1);
  // [w15.i7+w14.i6 w13.i5+w12.i4 w7.i7+w6.i6 w5.i5+w4.i4]
  // [w11.i3+w10.i2 w9.i1+w8.i0 w3.i3+w2.i2 w1.i1+w0.i0]
  const __m256 hadd1 = _mm256_hadd_ps(mul2, mul3);
  // [w31.i7+w30.i6 w29.i5+w28.i4 w23.i7+w22.i6 w21.i5+w20.i4]
  // [w27.i3+w26.i2 w25.i1+w24.i0 w19.i3+w18.i2 w17.i1+w16.i0]

  const __m256 haddhadd = _mm256_hadd_ps(hadd0, hadd1);
  // [w31.i7+w30.i6+w29.i5+w28.i4 w23.i7+w22.i6+w21.i5+w20.i4
  //  w15.i7+w14.i6+w13.i5+w12.i4 w7.i7+w6.i6+w5.i5+w4.i4]
  // [w27.i3+w26.i2+w25.i1+w24.i0 w19.i3+w18.i2+w17.i1+w16.i0
  //  w11.i3+w10.i2+w9.i1+w8.i0 w3.i3+w2.i2+w1.i1+w0.i0]

  // Permute control: [7 6 5 4][3 2 1 0] => [7 3 6 2][5 1 4 0]
  const __m256i ctrl = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
  const __m256 permute = _mm256_permutevar8x32_ps(haddhadd, ctrl);
  // [31+30+29+28 27+26+25+24 23+22+21+20 19+18+17+16]
  // [15+14+13+12 11+10+9+8 7+6+5+4 3+2+1+0]

  const __m128 permute_h = _mm256_extractf128_ps(permute, 1);
  const __m128 permute_l = _mm256_castps256_ps128(permute);
  const __m128 hhh = _mm_hadd_ps(permute_l, permute_h);
  // [31+30+29+28+27+26+25+24 23+22+21+20+19+18+17+16
  //  15+14+13+12+11+10+9+8 7+6+5+4+3+2+1+0
  *outputs = _mm_add_ps(*outputs, hhh);
}

static void nn_activate8(__m256 *x) {
  *x = _mm256_max_ps(*x, _mm256_setzero_ps());
}

static void nn_activate4(__m128 *x) { *x = _mm_max_ps(*x, _mm_setzero_ps()); }

// Calculate prediction based on the given input features and neural net config.
// Assume there are no more than NN_MAX_NODES_PER_LAYER nodes in each hidden
// layer.
void av1_nn_predict_avx2(const float *input_nodes,
                         const NN_CONFIG *const nn_config,
                         float *const output) {
  float buf[2][NN_MAX_NODES_PER_LAYER];
  int buf_index = 0;
  int num_inputs = nn_config->num_inputs;

  // Hidden layers, except the final iteration is the output layer.
  for (int layer = 0; layer <= nn_config->num_hidden_layers; layer++) {
    const float *layer_weights = nn_config->weights[layer];
    const float *layer_bias = nn_config->bias[layer];
    bool output_layer = (layer == nn_config->num_hidden_layers);
    float *const output_nodes = output_layer ? output : buf[buf_index];
    const int num_outputs = output_layer ? nn_config->num_outputs
                                         : nn_config->num_hidden_nodes[layer];

    if (num_inputs % 8 == 0 && num_outputs % 8 == 0) {
      for (int out = 0; out < num_outputs; out += 8) {
        __m256 outputs = _mm256_loadu_ps(&layer_bias[out]);
        for (int in = 0; in < num_inputs; in += 8) {
          nn_propagate_8to8(&input_nodes[in],
                            &layer_weights[out * num_inputs + in], &outputs,
                            num_inputs);
        }
        if (!output_layer) nn_activate8(&outputs);
        _mm256_storeu_ps(&output_nodes[out], outputs);
      }
    } else if (num_inputs % 4 == 0 && num_outputs % 8 == 0) {
      for (int out = 0; out < num_outputs; out += 8) {
        __m256 outputs = _mm256_loadu_ps(&layer_bias[out]);
        for (int in = 0; in < num_inputs; in += 4) {
          nn_propagate_4to8(&input_nodes[in],
                            &layer_weights[out * num_inputs + in], &outputs,
                            num_inputs);
        }
        if (!output_layer) nn_activate8(&outputs);
        _mm256_storeu_ps(&output_nodes[out], outputs);
      }
    } else if (num_inputs % 8 == 0 && num_outputs % 4 == 0) {
      for (int out = 0; out < num_outputs; out += 4) {
        __m128 outputs = _mm_loadu_ps(&layer_bias[out]);
        for (int in = 0; in < num_inputs; in += 8) {
          nn_propagate_8to4(&input_nodes[in],
                            &layer_weights[out * num_inputs + in], &outputs,
                            num_inputs);
        }
        if (!output_layer) nn_activate4(&outputs);
        _mm_storeu_ps(&output_nodes[out], outputs);
      }
    } else if (num_inputs % 8 == 0) {
      for (int out = 0; out < num_outputs; out++) {
        float total = layer_bias[out];
        for (int in = 0; in < num_inputs; in += 8) {
          total += nn_propagate_8to1(&input_nodes[in],
                                     &layer_weights[out * num_inputs + in]);
        }
        output_nodes[out] = output_layer ? total : AOMMAX(total, 0.0f);
      }
    } else if (num_inputs % 4 == 0) {
      for (int out = 0; out < num_outputs; out++) {
        float total = layer_bias[out];
        for (int in = 0; in < num_inputs; in += 4) {
          total += nn_propagate_4to1(&input_nodes[in],
                                     &layer_weights[out * num_inputs + in]);
        }
        output_nodes[out] = output_layer ? total : AOMMAX(total, 0.0f);
      }
    } else {
      for (int out = 0; out < num_outputs; out++) {
        float total = layer_bias[out];
        for (int in_node = 0; in_node < num_inputs; in_node++) {
          total +=
              input_nodes[in_node] * layer_weights[num_inputs * out + in_node];
        }
        output_nodes[out] = output_layer ? total : AOMMAX(total, 0.0f);
      }
    }
    input_nodes = output_nodes;
    num_inputs = num_outputs;
    buf_index = 1 - buf_index;
  }
}
