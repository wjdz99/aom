/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <pmmintrin.h>
#include <immintrin.h>

void av1_nn_propagate_4to1(const float *const inputs,
                           const float *const weights, __m128 *const output);

void av1_nn_propagate_4to4(const float *const inputs,
                           const float *const weights, __m128 *const outputs,
                           const int num_inputs);

void av1_nn_propagate_4to8(const float *const inputs,
                           const float *const weights, __m128 *const out_h,
                           __m128 *const out_l, const int num_inputs);
