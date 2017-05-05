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

#ifndef AOM_DSP_BITREADER_BUFFER_H_
#define AOM_DSP_BITREADER_BUFFER_H_

#include <limits.h>

#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*AomRbErrorHandler)(void *data);

struct AomReadBitBuffer {
  const uint8_t *bit_buffer;
  const uint8_t *bit_buffer_end;
  uint32_t bit_offset;

  void *error_handler_data;
  AomRbErrorHandler error_handler;
};

size_t aom_rb_bytes_read(struct AomReadBitBuffer *rb);

int aom_rb_read_bit(struct AomReadBitBuffer *rb);

int aom_rb_read_literal(struct AomReadBitBuffer *rb, int bits);

int aom_rb_read_signed_literal(struct AomReadBitBuffer *rb, int bits);

int aom_rb_read_inv_signed_literal(struct AomReadBitBuffer *rb, int bits);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITREADER_BUFFER_H_
