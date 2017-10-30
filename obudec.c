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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aom_ports/mem_ops.h"

#include "./obudec.h"

#if CONFIG_OBU_NO_IVF
int obu_read_temporal_unit(FILE *infile, uint8_t **buffer, size_t *bytes_read,
                           size_t *buffer_size) {
  size_t ret;
  const size_t obu_length_header_size = 5;
  uint32_t obu_size = 0;
  uint8_t *data = NULL;

  if (feof(infile)) {
    return 1;
  }

  *buffer_size = 0;
  *bytes_read = 0;
  while (1) {
    // augment the buffer to just contain the next size field
    // and the first byte of the header
    *buffer = realloc(*buffer, (*buffer_size) + obu_length_header_size);
    data = *buffer + (*buffer_size);
    *buffer_size += obu_length_header_size;
    ret = fread(data, 1, obu_length_header_size, infile);
    if (ret == 0) {
      fprintf(stderr, "Found end of stream, ending temporal unit\n");
      break;
    }
    if (ret != obu_length_header_size) {
      warn("Failed to read OBU Header\n");
      return 1;
    }
    *bytes_read += obu_length_header_size;

    if (((data[4] >> 3) & 0xF) == 2) {
      // Stop when a temporal delimiter is found
      fprintf(stderr, "Found temporal delimiter, ending temporal unit\n");
      // prevent decoder to start decoding another frame from this buffer
      *bytes_read -= obu_length_header_size;
      break;
    }

    // otherwise, read the OBU payload into memory
    obu_size = mem_get_le32(data);
    fprintf(stderr, "Found OBU of type %d and size %d\n",
            ((data[4] >> 3) & 0xF), obu_size);
    obu_size--;  // removing the byte of the header already read
    if (obu_size) {
      *buffer = realloc(*buffer, (*buffer_size) + obu_size);
      data = *buffer + (*buffer_size);
      *buffer_size += obu_size;
      ret = fread(data, 1, obu_size, infile);
      if (ret != obu_size) {
        warn("Failed to read OBU Payload\n");
        return 1;
      }
      *bytes_read += obu_size;
    }
  }
  return 0;
}

int file_is_obu(struct AvxInputContext *input_ctx) {
  uint8_t obutd[5];
  int size;

  input_ctx->fourcc = 0;  // Should be AV01
  // The following fields should be read from the OBU SH
  input_ctx->width = 0;
  input_ctx->height = 0;
  input_ctx->framerate.numerator = 1;
  input_ctx->framerate.denominator = 0;

  // Reading the first OBU TD to enable TU end detection at TD start.
  fread(obutd, 1, 5, input_ctx->file);
  size = mem_get_le32(obutd);
  if (size != 1) {
    warn("Expected first OBU size to be 1, got %d", size);
    return 0;
  }
  if (((obutd[4] >> 3) & 0xF) != 2) {
    warn("Expected OBU TD at file start, got %d\n", obutd[4]);
    return 0;
  }
  fprintf(stderr, "Starting to parse OBU stream\n");
  return 1;
}

#endif
