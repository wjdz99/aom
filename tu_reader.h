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
#ifndef TU_READER_H_
#define TU_READER_H_

#include "./tools_common.h"

#ifdef __cplusplus
extern "C" {
#endif

// Returns 0 when file referred to by 'ctx' is in a supported format.
int tu_reader_open(struct AvxInputContext *ctx);

// Reads one Temporal Unit from the input file. Returns 0 when a TU is
// successfully read, 1 when end of file is reached, and less than 0 when an
// error occurs. Stores TU data in 'buffer'. Reallocs buffer to match TU size,
// returns buffer capacity via 'buffer_size', and returns size of buffered data
// via 'bytes_read'.
int tu_reader_read_temporal_unit(struct AvxInputContext *ctx,
                                 uint8_t **buffer, size_t *bytes_read,
                                 size_t *buffer_size);

void tu_reader_close(struct AvxInputContext *ctx);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TU_READER_H_
