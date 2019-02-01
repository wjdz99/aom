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

#ifndef AOM_AV1_ENCODER_PARTITION_MODEL_WEIGHTS_H_
#define AOM_AV1_ENCODER_PARTITION_MODEL_WEIGHTS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/ml.h"

#define FEATURE_SIZE 10
#define LABEL_SIZE 16
// nn model for ab partition pruning, 128x128.
extern const float av1_ab_partition_nn_weights_128_layer0[FEATURE_SIZE * 64];
extern const float av1_ab_partition_nn_bias_128_layer0[64];

extern const float av1_ab_partition_nn_weights_128_layer1[64 * LABEL_SIZE];

extern const float av1_ab_partition_nn_bias_128_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_ab_partition_nnconfig_128;

// nn model for ab partition pruning, 64x64.
extern const float av1_ab_partition_nn_weights_64_layer0[FEATURE_SIZE * 64];

extern const float av1_ab_partition_nn_bias_64_layer0[64];

extern const float av1_ab_partition_nn_weights_64_layer1[64 * LABEL_SIZE];

extern const float av1_ab_partition_nn_bias_64_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_ab_partition_nnconfig_64;

// nn model for ab partition pruning, 32x32.
extern const float av1_ab_partition_nn_weights_32_layer0[FEATURE_SIZE * 64];

extern const float av1_ab_partition_nn_bias_32_layer0[64];

extern const float av1_ab_partition_nn_weights_32_layer1[64 * 16];

extern const float av1_ab_partition_nn_bias_32_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_ab_partition_nnconfig_32;

// nn model for ab partition pruning, 16x16.
extern const float av1_ab_partition_nn_weights_16_layer0[FEATURE_SIZE * 64];

extern const float av1_ab_partition_nn_bias_16_layer0[64];

extern const float av1_ab_partition_nn_weights_16_layer1[64 * LABEL_SIZE];

extern const float av1_ab_partition_nn_bias_16_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_ab_partition_nnconfig_16;

#undef FEATURE_SIZE
#undef LABEL_SIZE

#define FEATURE_SIZE 18
#define LABEL_SIZE 4

extern const float av1_4_partition_nn_weights_16_layer0[FEATURE_SIZE * 24];

extern const float av1_4_partition_nn_bias_16_layer0[24];

extern const float av1_4_partition_nn_weights_16_layer1[24 * LABEL_SIZE];

extern const float av1_4_partition_nn_bias_16_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_4_partition_nnconfig_16;

extern const float av1_4_partition_nn_weights_32_layer0[FEATURE_SIZE * 32];

extern const float av1_4_partition_nn_bias_32_layer0[32];

extern const float av1_4_partition_nn_weights_32_layer1[32 * LABEL_SIZE];

extern const float av1_4_partition_nn_bias_32_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_4_partition_nnconfig_32;

extern const float av1_4_partition_nn_weights_64_layer0[FEATURE_SIZE * 24];

extern const float av1_4_partition_nn_bias_64_layer0[24];

extern const float av1_4_partition_nn_weights_64_layer1[24 * LABEL_SIZE];

extern const float av1_4_partition_nn_bias_64_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_4_partition_nnconfig_64;

#undef FEATURE_SIZE
#undef LABEL_SIZE

#define FEATURE_SIZE 4
extern const float
    av1_partition_breakout_nn_weights_128_layer0[FEATURE_SIZE * 32];

extern const float av1_partition_breakout_nn_bias_128_layer0[32];

extern const float av1_partition_breakout_nn_weights_128_layer1[32];

extern const float av1_partition_breakout_nn_bias_128_layer1[1];

extern const NN_CONFIG av1_partition_breakout_nnconfig_128;

extern const float
    av1_partition_breakout_nn_weights_64_layer0[FEATURE_SIZE * 16];

extern const float av1_partition_breakout_nn_bias_64_layer0[16];

extern const float av1_partition_breakout_nn_weights_64_layer1[16];

extern const float av1_partition_breakout_nn_bias_64_layer1[1];

extern const NN_CONFIG av1_partition_breakout_nnconfig_64;

extern const float
    av1_partition_breakout_nn_weights_32_layer0[FEATURE_SIZE * 16];

extern const float av1_partition_breakout_nn_bias_32_layer0[16];

extern const float av1_partition_breakout_nn_weights_32_layer1[16];

extern const float av1_partition_breakout_nn_bias_32_layer1[1];

extern const NN_CONFIG av1_partition_breakout_nnconfig_32;

extern const float
    av1_partition_breakout_nn_weights_16_layer0[FEATURE_SIZE * 16];

extern const float av1_partition_breakout_nn_bias_16_layer0[16];

extern const float av1_partition_breakout_nn_weights_16_layer1[16];

extern const float av1_partition_breakout_nn_bias_16_layer1[1];

extern const NN_CONFIG av1_partition_breakout_nnconfig_16;

extern const float
    av1_partition_breakout_nn_weights_8_layer0[FEATURE_SIZE * 16];

extern const float av1_partition_breakout_nn_bias_8_layer0[16];

extern const float av1_partition_breakout_nn_weights_8_layer1[16];

extern const float av1_partition_breakout_nn_bias_8_layer1[1];

extern const NN_CONFIG av1_partition_breakout_nnconfig_8;
#undef FEATURE_SIZE

#define FEATURE_SIZE 9  // Input layer size
#define NUM_NODES 32    // Hidden layer size
#define LABEL_SIZE 3    // Output layer size

extern const float
    av1_rect_partition_nn_weights_8_layer0[FEATURE_SIZE * NUM_NODES];

extern const float av1_rect_partition_nn_bias_8_layer0[NUM_NODES];

extern const float
    av1_rect_partition_nn_weights_8_layer1[NUM_NODES * LABEL_SIZE];

extern const float av1_rect_partition_nn_bias_8_layer1[LABEL_SIZE];

extern const NN_CONFIG av1_rect_partition_nnconfig_8;

extern const float
    av1_rect_partition_nn_weights_16_layer0[FEATURE_SIZE * NUM_NODES];

extern const float av1_rect_partition_nn_bias_16_layer0[NUM_NODES];

extern const float
    av1_rect_partition_nn_weights_16_layer1[NUM_NODES * LABEL_SIZE];

extern const float av1_rect_partition_nn_bias_16_layer1[3];

extern const NN_CONFIG av1_rect_partition_nnconfig_16;

extern const float
    av1_rect_partition_nn_weights_32_layer0[FEATURE_SIZE * NUM_NODES];

extern const float av1_rect_partition_nn_bias_32_layer0[NUM_NODES];

extern const float
    av1_rect_partition_nn_weights_32_layer1[NUM_NODES * LABEL_SIZE];

extern const float av1_rect_partition_nn_bias_32_layer1[3];

extern const NN_CONFIG av1_rect_partition_nnconfig_32;

extern const float
    av1_rect_partition_nn_weights_64_layer0[FEATURE_SIZE * NUM_NODES];

extern const float av1_rect_partition_nn_bias_64_layer0[NUM_NODES];

extern const float
    av1_rect_partition_nn_weights_64_layer1[NUM_NODES * LABEL_SIZE];

extern const float av1_rect_partition_nn_bias_64_layer1[3];

extern const NN_CONFIG av1_rect_partition_nnconfig_64;

extern const float
    av1_rect_partition_nn_weights_128_layer0[FEATURE_SIZE * NUM_NODES];

extern const float av1_rect_partition_nn_bias_128_layer0[NUM_NODES];

extern const float
    av1_rect_partition_nn_weights_128_layer1[NUM_NODES * LABEL_SIZE];

extern const float av1_rect_partition_nn_bias_128_layer1[3];

extern const NN_CONFIG av1_rect_partition_nnconfig_128;
#undef FEATURE_SIZE
#undef NUM_NODES
#undef LABEL_SIZE

#if CONFIG_ONE_PASS_SVM
#define FEATURE_SIZE 24
extern const float av1_op_svm_early_term_weights_128[FEATURE_SIZE + 1];

extern const float av1_op_svm_early_term_weights_64[FEATURE_SIZE + 1];

extern const float av1_op_svm_early_term_weights_32[FEATURE_SIZE + 1];

extern const float av1_op_svm_early_term_weights_16[FEATURE_SIZE + 1];

extern const float av1_op_svm_early_term_mean_128[FEATURE_SIZE];

extern const float av1_op_svm_early_term_mean_64[FEATURE_SIZE];

extern const float av1_op_svm_early_term_mean_32[FEATURE_SIZE];

extern const float av1_op_svm_early_term_mean_16[FEATURE_SIZE];

extern const float av1_op_svm_early_term_std_128[FEATURE_SIZE];

extern const float av1_op_svm_early_term_std_64[FEATURE_SIZE];

extern const float av1_op_svm_early_term_std_32[FEATURE_SIZE];

extern const float av1_op_svm_early_term_std_16[FEATURE_SIZE];

#undef FEATURE_SIZE
#endif  // CONFIG_ONE_PASS_SVM

// Below are the models used for full_pixel_motion_search_based_split
static float simple_motion_search_prune_rect_thresh_128 = 0.0110f;
static float simple_motion_search_prune_rect_thresh_64 = 0.0474f;
static float simple_motion_search_prune_rect_thresh_32 = 0.0293f;
static float simple_motion_search_prune_rect_thresh_16 = 0.0180f;
static float simple_motion_search_prune_rect_thresh_8 = 0.0f;

// BLOCK_128X128
#define NUM_HIDDEN_LAYERS_128 1
#define NUM_FEATURES_128 6
#define NUM_LAYER_0_UNITS_128 16
#define NUM_LOGITS_128 1

extern const float full_pixel_motion_search_based_split_layer_0_kernel_128[];

extern const float full_pixel_motion_search_based_split_logits_kernel_128[];

extern const float full_pixel_motion_search_based_split_layer_0_bias_128[];

extern const float full_pixel_motion_search_based_split_logits_bias_128[];

extern const NN_CONFIG full_pixel_motion_search_based_split_nn_config_128;

#undef NUM_HIDDEN_LAYERS_128
#undef NUM_FEATURES_128
#undef NUM_LAYER_0_UNITS_128
#undef NUM_LOGITS_128

// BLOCK_64X64
#define NUM_HIDDEN_LAYERS_64 1
#define NUM_FEATURES_64 6
#define NUM_LAYER_0_UNITS_64 16
#define NUM_LOGITS_64 1

extern const float full_pixel_motion_search_based_split_layer_0_kernel_64[];

extern const float full_pixel_motion_search_based_split_logits_kernel_64[];

extern const float full_pixel_motion_search_based_split_layer_0_bias_64[];

extern const float full_pixel_motion_search_based_split_logits_bias_64[];

extern const NN_CONFIG full_pixel_motion_search_based_split_nn_config_64;

#undef NUM_HIDDEN_LAYERS_64
#undef NUM_FEATURES_64
#undef NUM_LAYER_0_UNITS_64
#undef NUM_LOGITS_64

// BLOCK_32X32
#define NUM_HIDDEN_LAYERS_32 1
#define NUM_FEATURES_32 6
#define NUM_LAYER_0_UNITS_32 16
#define NUM_LOGITS_32 1

extern const float full_pixel_motion_search_based_split_layer_0_kernel_32[];

extern const float full_pixel_motion_search_based_split_logits_kernel_32[];

extern const float full_pixel_motion_search_based_split_layer_0_bias_32[];

extern const float full_pixel_motion_search_based_split_logits_bias_32[];

extern const NN_CONFIG full_pixel_motion_search_based_split_nn_config_32;

#undef NUM_HIDDEN_LAYERS_32
#undef NUM_FEATURES_32
#undef NUM_LAYER_0_UNITS_32
#undef NUM_LOGITS_32

// BLOCK_16X16
#define NUM_HIDDEN_LAYERS_16 1
#define NUM_FEATURES_16 6
#define NUM_LAYER_0_UNITS_16 16
#define NUM_LOGITS_16 1

extern const float full_pixel_motion_search_based_split_layer_0_kernel_16[];

extern const float full_pixel_motion_search_based_split_logits_kernel_16[];

extern const float full_pixel_motion_search_based_split_layer_0_bias_16[];

extern const float full_pixel_motion_search_based_split_logits_bias_16[];

extern const NN_CONFIG full_pixel_motion_search_based_split_nn_config_16;

#undef NUM_HIDDEN_LAYERS_16
#undef NUM_FEATURES_16
#undef NUM_LAYER_0_UNITS_16
#undef NUM_LOGITS_16

#if !CONFIG_DISABLE_FULL_PIXEL_SPLIT_8X8
// BLOCK_8X8
#define NUM_HIDDEN_LAYERS_8 1
#define NUM_FEATURES_8 6
#define NUM_LAYER_0_UNITS_8 16
#define NUM_LOGITS_8 1

extern const float full_pixel_motion_search_based_split_layer_0_kernel_8[];

extern const float full_pixel_motion_search_based_split_logits_kernel_8[];

extern const float full_pixel_motion_search_based_split_layer_0_bias_8[];

extern const float full_pixel_motion_search_based_split_logits_bias_8[];

extern const NN_CONFIG full_pixel_motion_search_based_split_nn_config_8;

#endif

// Model based on simple_motion_search

// Thresholds for doing a single type of partition
// TODO(chiyotsai@google.com): Set the thresholds for PARTITION_SPLIT.
extern const float simple_motion_search_prune_part_only_thresh_128[10];
extern const float simple_motion_search_prune_part_only_thresh_64[10];
extern const float simple_motion_search_prune_part_only_thresh_32[10];
extern const float simple_motion_search_prune_part_only_thresh_16[10];
extern const float simple_motion_search_prune_part_only_thresh_8[10];

// Thresholds for pruning a partition type
// TODO(chiyotsai@google.com): Retune the thresholds for rectangular partition.
extern const float simple_motion_search_prune_part_prune_thresh_128[10];
extern const float simple_motion_search_prune_part_prune_thresh_64[10];
extern const float simple_motion_search_prune_part_prune_thresh_32[10];
extern const float simple_motion_search_prune_part_prune_thresh_16[10];
extern const float simple_motion_search_prune_part_prune_thresh_8[10];

// Mean and std
extern const float simple_motion_search_prune_part_mean_128[19];
extern const float simple_motion_search_prune_part_std_128[19];
extern const float simple_motion_search_prune_part_mean_64[19];
extern const float simple_motion_search_prune_part_std_64[19];
extern const float simple_motion_search_prune_part_mean_32[19];
extern const float simple_motion_search_prune_part_std_32[19];
extern const float simple_motion_search_prune_part_mean_16[19];
extern const float simple_motion_search_prune_part_std_16[19];
extern const float simple_motion_search_prune_part_mean_8[19];
extern const float simple_motion_search_prune_part_std_8[19];

// BLOCK_128X128
#define NUM_HIDDEN_LAYERS_128 1
#define NUM_FEATURES_128 19
#define NUM_LAYER_0_UNITS_128 24
#define NUM_LOGITS_128 4

extern const float simple_motion_search_prune_part_logits_kernel_128[];

extern const float simple_motion_search_prune_part_layer_0_kernel_128[];

extern const float simple_motion_search_prune_part_logits_bias_128[];

extern const float simple_motion_search_prune_part_layer_0_bias_128[];

extern const NN_CONFIG simple_motion_search_prune_part_nn_config_128;

#undef NUM_HIDDEN_LAYERS_128
#undef NUM_FEATURES_128
#undef NUM_LAYER_0_UNITS_128
#undef NUM_LOGITS_128

// BLOCK_64X64
#define NUM_HIDDEN_LAYERS_64 1
#define NUM_FEATURES_64 19
#define NUM_LAYER_0_UNITS_64 24
#define NUM_LOGITS_64 10

extern const float simple_motion_search_prune_part_logits_kernel_64[];

extern const float simple_motion_search_prune_part_layer_0_kernel_64[];

extern const float simple_motion_search_prune_part_logits_bias_64[];

extern const float simple_motion_search_prune_part_layer_0_bias_64[];

extern const NN_CONFIG simple_motion_search_prune_part_nn_config_64;

#undef NUM_HIDDEN_LAYERS_64
#undef NUM_FEATURES_64
#undef NUM_LAYER_0_UNITS_64
#undef NUM_LOGITS_64

// BLOCK_32X32
#define NUM_HIDDEN_LAYERS_32 1
#define NUM_FEATURES_32 19
#define NUM_LAYER_0_UNITS_32 24
#define NUM_LOGITS_32 10

extern const float simple_motion_search_prune_part_logits_kernel_32[];

extern const float simple_motion_search_prune_part_layer_0_kernel_32[];

extern const float simple_motion_search_prune_part_logits_bias_32[];

extern const float simple_motion_search_prune_part_layer_0_bias_32[];

extern const NN_CONFIG simple_motion_search_prune_part_nn_config_32;

#undef NUM_HIDDEN_LAYERS_32
#undef NUM_FEATURES_32
#undef NUM_LAYER_0_UNITS_32
#undef NUM_LOGITS_32

// BLOCK_16X16
#define NUM_HIDDEN_LAYERS_16 1
#define NUM_FEATURES_16 19
#define NUM_LAYER_0_UNITS_16 16
#define NUM_LOGITS_16 10

extern const float simple_motion_search_prune_part_logits_kernel_16[];

extern const float simple_motion_search_prune_part_layer_0_kernel_16[];

extern const float simple_motion_search_prune_part_logits_bias_16[];

extern const float simple_motion_search_prune_part_layer_0_bias_16[];

extern const NN_CONFIG simple_motion_search_prune_part_nn_config_16;

#undef NUM_HIDDEN_LAYERS_16
#undef NUM_FEATURES_16
#undef NUM_LAYER_0_UNITS_16
#undef NUM_LOGITS_16

// BLOCK_8X8
#define NUM_HIDDEN_LAYERS_8 1
#define NUM_FEATURES_8 19
#define NUM_LAYER_0_UNITS_8 24
#define NUM_LOGITS_8 4

extern const float simple_motion_search_prune_part_logits_kernel_8[];

extern const float simple_motion_search_prune_part_layer_0_kernel_8[];

extern const float simple_motion_search_prune_part_logits_bias_8[];

extern const float simple_motion_search_prune_part_layer_0_bias_8[];

extern const NN_CONFIG simple_motion_search_prune_part_nn_config_8;

#undef NUM_HIDDEN_LAYERS_8
#undef NUM_FEATURES_8
#undef NUM_LAYER_0_UNITS_8
#undef NUM_LOGITS_8

#define FEATURE_SIZE 19
extern const float two_pass_split_partition_weights_128[FEATURE_SIZE + 1];

extern const float two_pass_split_partition_weights_64[FEATURE_SIZE + 1];

extern const float two_pass_split_partition_weights_32[FEATURE_SIZE + 1];

extern const float two_pass_split_partition_weights_16[FEATURE_SIZE + 1];

extern const float two_pass_split_partition_weights_8[FEATURE_SIZE + 1];

extern const float two_pass_none_partition_weights_128[FEATURE_SIZE + 1];

extern const float two_pass_none_partition_weights_64[FEATURE_SIZE + 1];

extern const float two_pass_none_partition_weights_32[FEATURE_SIZE + 1];

extern const float two_pass_none_partition_weights_16[FEATURE_SIZE + 1];

extern const float two_pass_none_partition_weights_8[FEATURE_SIZE + 1];

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_PARTITION_MODEL_WEIGHTS_H_
