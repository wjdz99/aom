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
#ifndef IO_OBUDEC_H_
#define IO_OBUDEC_H_

#include "./tools_common.h"

#include "av1/decoder/obu.h"

#ifdef __cplusplus
extern "C" {
#endif

#define OBU_BUFFER_SIZE (500 * 1024)

#define OBU_HEADER_SIZE 1
#define OBU_EXTENSION_SIZE 1
#define OBU_MAX_LENGTH_FIELD_SIZE 8
#define OBU_DETECTION_SIZE \
  (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + 2 * OBU_MAX_LENGTH_FIELD_SIZE)

struct ObuDecInputContext {
  struct AvxInputContext *avx_ctx;
  uint8_t *buffer;
  size_t buffer_capacity;
  size_t bytes_buffered;
  int is_annexb;
  int last_layer_id;
};

// Returns 1 when file data starts (if Annex B stream, after reading the
// size of the OBU) with what appears to be a Temporal Delimiter
// OBU as defined by Section 5 of the AV1 bitstream specification.
int file_is_obu(struct ObuDecInputContext *obu_ctx);

// Reads unsigned LEB128 integer and returns 0 upon successful read and decode.
// Stores raw bytes in 'value_buffer', length of the number in 'value_length',
// and decoded value in 'value'.
int obudec_read_leb128(FILE *f, uint8_t *value_buffer, size_t *value_length,
                       uint64_t *value);

// Reads OBU header from 'f'. The 'buffer_capacity' passed in must be large
// enough to store an OBU header with extension (2 bytes). Raw OBU data is
// written to 'obu_data', parsed OBU header values are written to 'obu_header',
// and total bytes read from file are written to 'bytes_read'. Returns 0 for
// success, and non-zero on failure. When end of file is reached, the return
// value is 0 and the 'bytes_read' value is set to 0.
int obudec_read_obu_header(FILE *f, size_t buffer_capacity, int is_annexb,
                           uint8_t *obu_data, ObuHeader *obu_header,
                           size_t *bytes_read);

// Reads OBU payload from 'f' and returns 0 for success when all payload bytes
// are read from the file. Payload data is written to 'obu_data', and actual
// bytes read written to 'bytes_read'.
int obudec_read_obu_payload(FILE *f, uint64_t payload_length, uint8_t *obu_data,
                            size_t *bytes_read);

// Reads one Temporal Unit from the input file. Returns 0 when a TU is
// successfully read, 1 when end of file is reached, and less than 0 when an
// error occurs. Stores TU data in 'buffer'. Reallocs buffer to match TU size,
// returns buffer capacity via 'buffer_size', and returns size of buffered data
// via 'bytes_read'.
int obudec_read_temporal_unit(struct ObuDecInputContext *obu_ctx,
                              uint8_t **buffer, size_t *bytes_read,
                              size_t *buffer_size);

void obudec_free(struct ObuDecInputContext *obu_ctx);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // IO_OBUDEC_H_
