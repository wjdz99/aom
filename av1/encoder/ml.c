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

#include <assert.h>
#include <math.h>

#include "aom_dsp/aom_dsp_common.h"
#include "av1/encoder/ml.h"

// Calculate prediction based on the given input features and neural net config.
// Assume there are no more than NN_MAX_NODES_PER_LAYER nodes in each hidden
// layer.
void av1_nn_predict_c(const float *input_nodes,
                      const NN_CONFIG *const nn_config, float *const output) {
  int num_input_nodes = nn_config->num_inputs;
  int buf_index = 0;
  float buf[2][NN_MAX_NODES_PER_LAYER];

  // Propagate hidden layers.
  const int num_layers = nn_config->num_hidden_layers;
  assert(num_layers <= NN_MAX_HIDDEN_LAYERS);
  for (int layer = 0; layer < num_layers; ++layer) {
    const float *layer_weights = nn_config->weights[layer];
    const float *layer_bias = nn_config->bias[layer];
    float *output_nodes = buf[buf_index];
    const int num_output_nodes = nn_config->num_hidden_nodes[layer];
    assert(num_output_nodes < NN_MAX_NODES_PER_LAYER);
    for (int node = 0; node < num_output_nodes; ++node) {
      float val = layer_bias[node];
      for (int i = 0; i < num_input_nodes; ++i)
        val += layer_weights[node * num_input_nodes + i] * input_nodes[i];
      // ReLU as activation function.
      val = val > 0.0f ? val : 0.0f;  // Could use AOMMAX().
      output_nodes[node] = val;
    }
    num_input_nodes = num_output_nodes;
    input_nodes = output_nodes;
    buf_index = 1 - buf_index;
  }

  // Final output layer.
  const float *layer_weights = nn_config->weights[num_layers];
  const float *layer_bias = nn_config->bias[num_layers];
  for (int node = 0; node < nn_config->num_outputs; ++node) {
    float val = layer_bias[node];
    for (int i = 0; i < num_input_nodes; ++i)
      val += layer_weights[node * num_input_nodes + i] * input_nodes[i];
    output[node] = val;
  }
}

#if CONFIG_NN_V2
/**************************** Forward prediction ******************************/
// Applies the ReLu activation to one fc layer
// output[i] = Max(input[i],0.0f)
static float *nn_relu(const float *input, FC_LAYER *layer) {
  for (int i = 0; i < layer->num_outputs; ++i) {
    layer->output[i] = AOMMAX(input[i], 0.0f);
  }

  return layer->output;
}

// Applies the Sigmoid activation to one fc layer
// output[i] = 1/(1+exp(input[i]))
static float *nn_sigmoid(const float *input, FC_LAYER *layer) {
  for (int i = 0; i < layer->num_outputs; ++i) {
    const float tmp = AOMMIN(AOMMAX(input[i], -10.0f), 10.0f);
    layer->output[i] = 1.0f / (1.0f + expf(-tmp));
  }

  return layer->output;
}

// Forward prediction in one fc layer, used in function av1_nn_predict_V2
static float *nn_fc_forward(const float *input, FC_LAYER *layer) {
  const float *weights = layer->weights;
  const float *bias = layer->bias;
  assert(layer->num_outputs < NN_MAX_NODES_PER_LAYER);
  // fc
  for (int node = 0; node < layer->num_outputs; ++node) {
    float val = bias[node];
    for (int i = 0; i < layer->num_inputs; ++i) val += weights[i] * input[i];
    layer->output[node] = val;
    weights += layer->num_inputs;
  }

  // activation
  switch (layer->activation) {
    case NONE: return layer->output;
    case RELU: return nn_relu(layer->output, layer);
    case SIGMOID: return nn_sigmoid(layer->output, layer);
    case SOFTSIGN:
      assert(0 && "Softsign has not been supported in NN.");  // TO DO
      return NULL;
    default:
      assert(0 && "Unknown activation");  // Unknown activation
      return NULL;
  }
}

void av1_nn_predict_v2(const float *feature, NN_CONFIG_V2 *nn_config,
                       float *output) {
  const float *input_nodes = feature;
  const int num_layers = nn_config->num_hidden_layers;
  assert(num_layers <= NN_MAX_HIDDEN_LAYERS);
  // Copy the input feature to the buffer
  memcpy(nn_config->feature, feature,
         sizeof(*feature) * nn_config->layer[0].num_inputs);

  // Propagate the layers.
  for (int i = 0; i < num_layers; ++i) {
    input_nodes = nn_fc_forward(input_nodes, nn_config->layer + i);
    assert(nn_config->layer[i + 1].num_inputs ==
           nn_config->layer[i].num_outputs);
  }

  // Final layer
  input_nodes = nn_fc_forward(input_nodes, nn_config->layer + num_layers);
  assert(nn_config->layer[num_layers].num_outputs == nn_config->num_logits);
  // Copy the final layer output
  memcpy(output, input_nodes, sizeof(*input_nodes) * nn_config->num_logits);
}

/***************************Backprop for gradient******************************/
// Backprop for ReLU activation
static void nn_relu_back(float *dX_out, FC_LAYER *layer) {
  const float *dY = layer->dY;
  for (int i = 0; i < layer->num_outputs; ++i)
    dX_out[i] = layer->output[i] > 0.0f ? dY[i] : 0.0f;
}

// Backprop for sigmoid activation
static void nn_sigmoid_back(float *dX_out, FC_LAYER *layer) {
  const float *dY = layer->dY;
  for (int i = 0; i < layer->num_outputs; ++i)
    dX_out[i] =
        dY[i] * layer->output[i] * (1 - layer->output[i]);  // dX=dY*sigmoid(X)
}

// Backprop for softmax cross entropy loss
static void nn_softmax_cross_entropy_loss_back(float *dX_out,
                                               const float *logits,
                                               const int num_logits,
                                               const int label) {
  if (num_logits == 1) {
    // sigmoid
    assert(label < 2);  // label [0,1]
    dX_out[0] = 1.0f / (1.0f + expf(-logits[0])) - (float)label;
  } else {
    // softmax
    assert(num_logits > label);  // label [0,1,... num_logits-1]
    av1_nn_softmax(logits, dX_out, num_logits);
    dX_out[label] -= 1;
  }
}

// Backprop in one fc layer, used in function av1_nn_backprop
static void nn_fc_backward(const float *X, float *dX_out, FC_LAYER *layer) {
  // backprop on activation
  float dY_fc[NN_MAX_NODES_PER_LAYER] = { 0 };  // dY for fc
  switch (layer->activation) {
    case NONE:  // no activation, dY_fc <-- dY
      memcpy(dY_fc, layer->dY, sizeof(*(layer->dY)) * layer->num_outputs);
      break;
    case RELU: nn_relu_back(dY_fc, layer); break;
    case SIGMOID: nn_sigmoid_back(dY_fc, layer); break;
    case SOFTSIGN:
      assert(0 && "Softsign has not been supported in NN.");  // TO DO
      break;
    default: assert(0 && "Unknown activation");  // Unknown activation
  }

  // backprop on fc
  // gradient of W, b
  float *dW = layer->dW;
  float *db = layer->db;
  for (int j = 0; j < layer->num_outputs; ++j) {
    for (int i = 0; i < layer->num_inputs; ++i) dW[i] += dY_fc[j] * X[i];
    db[j] += dY_fc[j];
    dW += layer->num_inputs;
  }

  // gradient of the input, i.e., the output of last layer
  if (dX_out) {
    for (int i = 0; i < layer->num_inputs; ++i) {
      float *w = layer->weights + i;
      float val = 0.0f;
      for (int j = 0; j < layer->num_outputs; ++j) {
        val += dY_fc[j] * w[j * layer->num_inputs];
      }
      dX_out[i] = val;
    }
  }
}

void av1_nn_backprop(NN_CONFIG_V2 *nn_config, const int label) {
  const int num_layers = nn_config->num_hidden_layers;

  // loss layer
  switch (nn_config->loss) {
    case SOFTMAX_CROSS_ENTROPY:
      nn_softmax_cross_entropy_loss_back(nn_config->layer[num_layers].dY,
                                         nn_config->logits,
                                         nn_config->num_logits, label);
      break;
  }

  // hidden fc layer
  FC_LAYER *layer_ptr = nn_config->layer + num_layers;
  for (int i = 0; i < num_layers; ++i) {
    nn_fc_backward(layer_ptr[-1].output, layer_ptr[-1].dY, layer_ptr);
    layer_ptr -= 1;
  }

  nn_fc_backward(nn_config->feature, NULL,
                 layer_ptr);  // first layer (no dX for feature)
  ++nn_config->counter;       // increment the counter
}

void av1_nn_outer_product_backprop(NN_CONFIG_V2 *nn_config_hor,
                                   NN_CONFIG_V2 *nn_config_ver,
                                   const int label) {
  assert(nn_config_hor->loss == nn_config_ver->loss);
  const int num_layers1 = nn_config_hor->num_hidden_layers;
  const int num_layers2 = nn_config_ver->num_hidden_layers;
  const float *logits1 = nn_config_hor->logits;
  const int n1 = nn_config_hor->num_logits;
  const float *logits2 = nn_config_ver->logits;
  const int n2 = nn_config_ver->num_logits;

  // outer product
  assert(n1 * n2 <= NN_MAX_NODES_PER_LAYER);
  float out_product[NN_MAX_NODES_PER_LAYER];
  for (int i = 0; i < n2; ++i) {
    for (int j = 0; j < n1; ++j) {
      out_product[i * n1 + j] = logits1[j] * logits2[i];
    }
  }

  // loss layer
  float dY[NN_MAX_NODES_PER_LAYER];
  switch (nn_config_hor->loss) {
    case SOFTMAX_CROSS_ENTROPY:
      nn_softmax_cross_entropy_loss_back(dY, out_product, n1 * n2, label);
      break;
  }

  // hidden fc layer in hor
  FC_LAYER *layer_ptr = nn_config_hor->layer + num_layers1;
  for (int i = 0; i < n1; ++i) {
    layer_ptr->dY[i] = 0;
    for (int j = 0; j < n2; ++j) {
      layer_ptr->dY[i] += dY[i + j * n1] * logits2[j];
    }
  }
  for (int i = 0; i < num_layers1; ++i) {
    nn_fc_backward(layer_ptr[-1].output, layer_ptr[-1].dY, layer_ptr);
    layer_ptr -= 1;
  }
  nn_fc_backward(nn_config_hor->feature, NULL,
                 layer_ptr);  // first layer (no dX of the feature)

  // hidden fc layer in ver
  layer_ptr = nn_config_ver->layer + num_layers2;
  for (int i = 0; i < n2; ++i) {
    layer_ptr->dY[i] = 0;
    for (int j = 0; j < n1; ++j) {
      layer_ptr->dY[i] += dY[i * n1 + j] * logits1[j];
    }
  }
  for (int i = 0; i < num_layers2; ++i) {
    nn_fc_backward(layer_ptr[-1].output, layer_ptr[-1].dY, layer_ptr);
    layer_ptr -= 1;
  }
  nn_fc_backward(nn_config_ver->feature, NULL,
                 layer_ptr);  // first layer (no dX of the feature)

  ++nn_config_hor->counter;  // increment the counter
  ++nn_config_ver->counter;
}

void av1_nn_outer_product_backprop_with_mask(NN_CONFIG_V2 *nn_config_hor,
                                             NN_CONFIG_V2 *nn_config_ver,
                                             const int label, uint16_t mask) {
  assert(nn_config_hor->loss == nn_config_ver->loss);
  const int num_layers1 = nn_config_hor->num_hidden_layers;
  const int num_layers2 = nn_config_ver->num_hidden_layers;
  const float *logits1 = nn_config_hor->logits;
  const int n1 = nn_config_hor->num_logits;
  const float *logits2 = nn_config_ver->logits;
  const int n2 = nn_config_ver->num_logits;

  // outer product
  assert(n1 * n2 <= NN_MAX_NODES_PER_LAYER);
  float out_product[NN_MAX_NODES_PER_LAYER];
  for (int i = 0; i < n2; ++i) {
    for (int j = 0; j < n1; ++j) {
      out_product[i * n1 + j] = logits1[j] * logits2[i];
    }
  }

  // loss layer
  float dY[NN_MAX_NODES_PER_LAYER];
  switch (nn_config_hor->loss) {
    case SOFTMAX_CROSS_ENTROPY:
      nn_softmax_cross_entropy_loss_back(dY, out_product, n1 * n2, label);
      break;
  }

  // mask
  for (int tx_type = 0; tx_type < TX_TYPES; ++tx_type) {
    if (mask & (1 << tx_type)) {
      dY[tx_type] = 0.f;
    }
  }

  // hidden fc layer in hor
  FC_LAYER *layer_ptr = nn_config_hor->layer + num_layers1;
  for (int i = 0; i < n1; ++i) {
    layer_ptr->dY[i] = 0;
    for (int j = 0; j < n2; ++j) {
      layer_ptr->dY[i] += dY[i + j * n1] * logits2[j];
    }
  }
  for (int i = 0; i < num_layers1; ++i) {
    nn_fc_backward(layer_ptr[-1].output, layer_ptr[-1].dY, layer_ptr);
    layer_ptr -= 1;
  }
  nn_fc_backward(nn_config_hor->feature, NULL,
                 layer_ptr);  // first layer (no dX of the feature)

  // hidden fc layer in ver
  layer_ptr = nn_config_ver->layer + num_layers2;
  for (int i = 0; i < n2; ++i) {
    layer_ptr->dY[i] = 0;
    for (int j = 0; j < n1; ++j) {
      layer_ptr->dY[i] += dY[i * n1 + j] * logits1[j];
    }
  }
  for (int i = 0; i < num_layers2; ++i) {
    nn_fc_backward(layer_ptr[-1].output, layer_ptr[-1].dY, layer_ptr);
    layer_ptr -= 1;
  }
  nn_fc_backward(nn_config_ver->feature, NULL,
                 layer_ptr);  // first layer (no dX of the feature)

  ++nn_config_hor->counter;  // increment the counter
  ++nn_config_ver->counter;
}

void av1_nn_update(NN_CONFIG_V2 *nn_config, float mu) {
  const int num_layers = nn_config->num_hidden_layers;

  // Update the weights
  mu /= nn_config->counter;
  for (int i = 0; i <= num_layers; ++i) {
    FC_LAYER layer = nn_config->layer[i];
    for (int j = 0; j < layer.num_inputs * layer.num_outputs; ++j) {
      layer.weights[j] -= mu * layer.dW[j];
      layer.dW[j] = 0.f;
    }
    for (int j = 0; j < layer.num_outputs; ++j) {
      layer.bias[j] -= mu * layer.db[j];
      layer.db[j] = 0.f;
    }
  }
  nn_config->counter = 0;  // reset the counter after update
}
#endif  // CONFIG_NN_V2

void av1_nn_softmax(const float *input, float *output, int n) {
  // Softmax function is invariant to adding the same constant
  // to all input values, so we subtract the maximum input to avoid
  // possible overflow.
  float max_inp = input[0];
  for (int i = 1; i < n; i++) max_inp = AOMMAX(max_inp, input[i]);
  float sum_out = 0.0f;
  for (int i = 0; i < n; i++) {
    // Clamp to range [-10.0, 0.0] to prevent FE_UNDERFLOW errors.
    const float normalized_input = AOMMAX(input[i] - max_inp, -10.0f);
    output[i] = (float)exp(normalized_input);
    sum_out += output[i];
  }
  for (int i = 0; i < n; i++) output[i] /= sum_out;
}
