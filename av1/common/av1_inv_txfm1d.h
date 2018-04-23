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

#ifndef AV1_INV_TXFM1D_H_
#define AV1_INV_TXFM1D_H_

#include "av1/common/av1_txfm.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE int32_t clamp_value(int32_t value, int8_t bit) {
  if (bit <= 0) return value;  // Do nothing for invalid clamp bit.
  const int64_t max_value = (1LL << (bit - 1)) - 1;
  const int64_t min_value = -(1LL << (bit - 1));
  return (int32_t)clamp64(value, min_value, max_value);
}

static INLINE void clamp_buf(int32_t *buf, int32_t size, int8_t bit) {
  for (int i = 0; i < size; ++i) buf[i] = clamp_value(buf[i], bit);
}

void av1_idct4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void av1_idct8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void av1_idct16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void av1_idct32_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void av1_idct64_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void av1_iadst4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void av1_iadst8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void av1_iadst16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                     const int8_t *stage_range);
void av1_iidentity4_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void av1_iidentity8_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void av1_iidentity16_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);
void av1_iidentity32_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);

#ifdef __cplusplus
}
#endif

#endif  // AV1_INV_TXFM1D_H_
