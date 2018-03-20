/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#ifndef OBUDEC_H_
#define OBUDEC_H_

#include "./tools_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ObuDecInputContext {
  struct AvxInputContext *avx_ctx;
  uint8_t *buffer;
  size_t buffer_capacity;
  size_t bytes_buffered;
};

int file_is_obu(struct ObuDecInputContext *obu_ctx);

int obu_read_temporal_unit(struct ObuDecInputContext *obu_ctx, uint8_t **buffer,
                           size_t *bytes_read,
#if CONFIG_SCALABILITY
                           size_t *buffer_size, int last_layer_id
#else
                           size_t *buffer_size
#endif
);

void obudec_free(struct ObuDecInputContext *obu_ctx);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // OBUDEC_H_
