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

#ifndef _AOM_DSP_INTRAPRED_COMMON_H
#define _AOM_DSP_INTRAPRED_COMMON_H

#include "./aom_config.h"

// Weights are quadratic from '1' to '1 / block_size', scaled by
// 2^sm_weight_log2_scale.
static const int sm_weight_log2_scale = 8;

#if CONFIG_TX64X64
// max(block_size_wide[BLOCK_LARGEST], block_size_high[BLOCK_LARGEST])
#define MAX_BLOCK_DIM 64
#else
#define MAX_BLOCK_DIM 32
#endif  // CONFIG_TX64X64

static const uint8_t sm_weight_arrays[2 * MAX_BLOCK_DIM] = {
  // Unused, because we always offset by bs, which is at least 2.
  0, 0,
  // bs = 2
  128, 1,
  // bs = 4
  192, 171, 107, 1,
  // bs = 8
  224, 219, 206, 183, 151, 110, 59, 1,
  // bs = 16
  240, 239, 236, 230, 223, 213, 202, 188, 172, 154, 133, 111, 96, 70, 31, 1,
  // bs = 32
  248, 248, 247, 246, 244, 242, 239, 235, 231, 227, 222, 217, 211, 204, 197,
  190, 182, 173, 164, 155, 145, 134, 123, 111, 99, 87, 74, 60, 46, 31, 16, 1,
#if CONFIG_TX64X64
  // bs = 64
  252, 252, 252, 251, 251, 250, 250, 249, 248, 247, 246, 244, 243, 241, 240,
  238, 236, 234, 231, 229, 227, 224, 221, 218, 215, 212, 209, 206, 202, 199,
  195, 191, 187, 183, 179, 174, 170, 165, 160, 155, 150, 145, 140, 135, 129,
  123, 118, 112, 106, 100, 93, 87, 80, 74, 67, 60, 53, 46, 38, 31, 23, 16, 8, 1
#endif  // CONFIG_TX64X64
};
/*
static const uint8_t sm_weight_arrays[2 * MAX_BLOCK_DIM] = {
  // Unused, because we always offset by bs, which is at least 2.
  0, 0,
  // bs = 2
  255, 128,
  // bs = 4
  255, 149, 85, 64,
  // bs = 8
  255, 197, 146, 105, 73, 50, 37, 32,
  // bs = 16
  255, 225, 196, 170, 145, 123, 102, 84, 68, 54, 43, 33, 26, 20, 17, 16,
  // bs = 32
  255, 240, 225, 210, 196, 182, 169, 157, 145, 133, 122, 111, 101, 92, 83, 74,
  66, 59, 52, 45, 39, 34, 29, 25, 21, 17, 14, 12, 10, 9, 8, 8,
#if CONFIG_TX64X64
  // bs = 64
  255, 248, 240, 233, 225, 218, 210, 203, 196, 189, 182, 176, 169, 163, 156,
  150, 144, 138, 133, 127, 121, 116, 111, 106, 101, 96, 91, 86, 82, 77, 73, 69,
  65, 61, 57, 54, 50, 47, 44, 41, 38, 35, 32, 29, 27, 25, 22, 20, 18, 16, 15,
  13, 12, 10, 9, 8, 7, 6, 6, 5, 5, 4, 4, 4,
#endif  // CONFIG_TX64X64
};
*/
#endif  // _AOM_DSP_INTRAPRED_COMMON_H
