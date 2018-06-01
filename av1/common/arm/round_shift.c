/*
 *
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>
#include <assert.h>

#include "aom_ports/mem.h"
#include "av1/common/arm/mem_neon.h"
#include "av1/common/arm/transpose_neon.h"

void av1_round_shift_array_neon(int32_t *arr, int size, int bit) {
  int i;
  int32x4_t tmp_q_s32;
  if (bit == 0) {
    return;
  } else {
    const int32x4_t dup_bits_n_32x4 = vdupq_n_s32((int32_t)(-bit));
    if (bit > 0) {
      for (i = 0; i < (size >> 2); i++) {
        tmp_q_s32 = vld1q_s32(arr);
        tmp_q_s32 = vqrshlq_s32(tmp_q_s32, dup_bits_n_32x4);
        vst1q_s32(arr, tmp_q_s32);
        arr = arr + 4;
      }
    } else {
      const int32x4_t one = vdupq_n_s32(1);
      for (i = 0; i < (size >> 2); i++) {
        tmp_q_s32 = vld1q_s32(arr);
        tmp_q_s32 = vmulq_s32(tmp_q_s32, vshlq_s32(one, dup_bits_n_32x4));
        vst1q_s32(arr, tmp_q_s32);
        arr = arr + 4;
      }
    }
  }
}
