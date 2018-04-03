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
#include "./tu_reader.h"

#include "./ivfdec.h"

int tu_reader_open(struct AvxInputContext *ctx) {
  if (file_is_ivf(ctx)) {
    ctx->file_type = FILE_TYPE_IVF;
    return 0;
  }
  return -1;
}

int tu_reader_read_temporal_unit(struct AvxInputContext *ctx,
                                 uint8_t **buffer, size_t *bytes_read,
                                 size_t *buffer_size) {
  return 0;
}
