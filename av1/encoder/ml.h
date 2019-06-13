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

#ifndef AOM_AV1_ENCODER_ML_H_
#define AOM_AV1_ENCODER_ML_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "config/av1_rtcd.h"

#define NN_MAX_HIDDEN_LAYERS 10
#define NN_MAX_NODES_PER_LAYER 128

enum Actn{none, relu, sigmoid};      // Activation funciton
enum Loss{softmax_cross_entropy};    // Loss function

struct NN_CONFIG {
  int num_inputs;         // Number of input nodes, i.e. features.
  int num_outputs;        // Number of output nodes.
  int num_hidden_layers;  // Number of hidden layers, maximum 10.
  // Number of nodes for each hidden layer.
  int num_hidden_nodes[NN_MAX_HIDDEN_LAYERS];
  // Weight parameters, indexed by layer.
  const float *weights[NN_MAX_HIDDEN_LAYERS + 1];
  // Bias parameters, indexed by layer.
  const float *bias[NN_MAX_HIDDEN_LAYERS + 1];
};
// Typedef from struct NN_CONFIG to NN_CONFIG is in rtcd_defs

// Fully-connectedly layer configuration
typedef struct {
  const int num_inputs;         // Number of input nodes, i.e. features.
  const int num_outputs;        // Number of output nodes.

  float *weights;               // Weight parameters.
  float *bias;                  // Bias parameters.
  const enum Actn activation;   // Activation function.

  float *output;                // The output array.
  float *dY;                    // Gradient of outputs
  float *dW;                    // Gradient of weights.
  float *db;                    // Gradient of bias
} FC_Layer;

// NN configure structure v2
typedef struct {
  const int num_hidden_layers;             // Number of hidden layers, max = 10.
  FC_Layer layer[NN_MAX_HIDDEN_LAYERS+1];  // The layer array
  const int num_logits;                    // Number of output nodes.
  float *logits;                           // Raw prediction before loss
  const enum Loss loss;                    // Loss function
} NN_CONFIG_v2;


// Calculate prediction based on the given input features and neural net config.
// Assume there are no more than NN_MAX_NODES_PER_LAYER nodes in each hidden
// layer.
void av1_nn_predict_v2(const float *features, NN_CONFIG_v2 *nn_config,
                        float *output);

// Forward prediction in one fc layer, used in function av1_nn_predict_v2
float* av1_nn_fc_forward(const float *input, FC_Layer *layer);

// Applies the ReLu activation to one fc layer
// output[i] = Max(input[i],0.0f)
float* av1_nn_relu(const float *input, FC_Layer *layer);

// Applies the Sigmoid activation to one fc layer
// output[i] = 1/(1+exp(input[i]))
float* av1_nn_sigmoid(const float *input, FC_Layer *layer);

// Applies the softmax normalization function to the input
// to get a valid probability distribution in the output:
// output[i] = exp(input[i]) / sum_{k \in [0,n)}(exp(input[k]))
void av1_nn_softmax(const float *input, float *output, int n);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_ML_H_
