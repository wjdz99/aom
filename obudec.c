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

#include "./obudec.h"

#include "aom_ports/mem_ops.h"
#include "av1/common/common.h"
#include "av1/decoder/obu.h"

#define OBU_HEADER_SIZE 1
#define OBU_EXTENSION_SIZE 1
#define OBU_DETECTION_SIZE (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + 1)

// Unsigned LEB128 OBU length field has maximum size of 8 bytes.
#define OBU_MAX_LENGTH_FIELD_SIZE 8
#define OBU_MAX_HEADER_SIZE \
  (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + OBU_MAX_LENGTH_FIELD_SIZE)

#if CONFIG_OBU_NO_IVF

// header_raw_bytes MUST point to a location that
// can store 2 full bytes.
static int read_obu_header(FILE *f, size_t *consumed, ObuHeader *header,
                           uint8_t *header_raw_bytes) {
  if (!f || !consumed || !header || !header_raw_bytes) return -1;

  memset(header_raw_bytes, 0, OBU_HEADER_SIZE + OBU_EXTENSION_SIZE);
  const size_t bytes_read = fread(header_raw_bytes, 1, OBU_HEADER_SIZE, f);

  aom_codec_err_t parse_result =
      aom_read_obu_header(header_raw_bytes, bytes_read, consumed, header);
  if (parse_result != AOM_CODEC_OK) return parse_result;

  if (header->has_extension) {
    header_raw_bytes[OBU_HEADER_SIZE] = fgetc(f);
    parse_result =
        aom_read_obu_header(header_raw_bytes, bytes_read, consumed, header);
  }

  return parse_result;
}

// Reads OBU size from infile and returns 0 upon success. Returns obu_size via
// output pointer obu_size. Returns -1 when reading or parsing fails. Always
// returns FILE pointer to position at time of call. Returns 0 and sets obu_size
// to 0 when end of file is reached.
static int read_obu_size(FILE *infile, uint64_t *obu_size,
                         size_t *length_field_size) {
  if (!infile || !obu_size) return 1;

  uint8_t read_buffer[OBU_MAX_LENGTH_FIELD_SIZE] = { 0 };
  size_t bytes_read = fread(read_buffer, 1, OBU_MAX_LENGTH_FIELD_SIZE, infile);
  *obu_size = 0;

  if (bytes_read == 0) {
    return 0;
  }

  const int seek_pos = (int)bytes_read;
  if (seek_pos != 0 && fseek(infile, -seek_pos, SEEK_CUR) != 0) return 1;

  if (aom_uleb_decode(read_buffer, bytes_read, obu_size, length_field_size) !=
      0) {
    return 1;
  }

  return 0;
}

int obu_read_temporal_unit(FILE *infile, uint8_t **buffer, size_t *bytes_read,
#if CONFIG_SCALABILITY
                           size_t *buffer_size, int last_layer_id) {
#else
                           size_t *buffer_size) {
#endif
  size_t ret;
  uint64_t obu_size = 0;
  uint8_t *data = NULL;

  if (feof(infile)) {
    return 1;
  }

  *buffer_size = 0;
  *bytes_read = 0;
  while (1) {
    uint8_t obu_header_raw_bytes[OBU_HEADER_SIZE + OBU_DETECTION_SIZE] = { 0 };
    size_t obu_header_size = 0;
    ObuHeader obu_header;
    memset(&obu_header, 0, sizeof(obu_header));

    // BROKEN: Must read bytes into data here.

    if (read_obu_header(infile, &obu_header_size, &obu_header,
                        &obu_header_raw_bytes[0]) != AOM_CODEC_OK) {
      fprintf(stderr, "Failed to read OBU header.\n");
      return 1;
    }

    size_t length_field_size = 0;
    ret = read_obu_size(infile, &obu_size, &length_field_size);
    if (ret != 0) {
      warn("Failed to read OBU size.\n");
      return 1;
    }

    if (obu_header.type == OBU_TEMPORAL_DELIMITER) {
      // Stop when a temporal delimiter is found
      // fprintf(stderr, "Found temporal delimiter, ending temporal unit\n");
      // prevent decoder to start decoding another frame from this buffer
      break;
    } else if (obu_size == 0) {
      fprintf(stderr, "Found end of stream, ending temporal unit\n");
      break;
    }

    const uint64_t obu_payload_size = obu_size + length_field_size;
    obu_size += length_field_size + obu_header.size;

    // Expand the buffer to contain the full OBU and its length.
    const uint8_t *old_buffer = *buffer;
    const uint64_t alloc_size = *buffer_size + obu_size;

#if defined(AOM_MAX_ALLOCABLE_MEMORY)
    *buffer = (alloc_size > AOM_MAX_ALLOCABLE_MEMORY)
                  ? NULL
                  : (uint8_t *)realloc(*buffer, (size_t)alloc_size);
#else
    *buffer = (uint8_t *)realloc(*buffer, (size_t)alloc_size);
#endif

    if (!*buffer) {
      free((void *)old_buffer);
      warn("OBU buffer alloc failed.\n");
      return 1;
    }

    data = *buffer + (*buffer_size);
    *buffer_size += (size_t)obu_size;

    *data = obu_header_raw_bytes[0];
    data += OBU_HEADER_SIZE;

    if (obu_header.has_extension) {
      *data = obu_header_raw_bytes[1];
      data += OBU_EXTENSION_SIZE;
    }

    ret = fread(data, 1, (size_t)obu_payload_size, infile);

    if (ret != obu_payload_size) {
      warn("Failed to read OBU.\n");
      return 1;
    }
    *bytes_read += (size_t)obu_size;

#if CONFIG_SCALABILITY
    // break if obu_extension_flag is found and enhancement_id change
    if (obu_header.has_extension) {
      const int curr_layer_id = obu_header.enhancement_layer_id;
      if (curr_layer_id && (curr_layer_id > last_layer_id)) {
        // new enhancement layer
        *bytes_read -= (size_t)obu_size;
        const int i_obu_size = (int)obu_size;
        fseek(infile, -i_obu_size, SEEK_CUR);
        break;
      }
    }
#endif
  }

  fprintf(stderr, "TU size: %zu\n", *bytes_read);

  return 0;
}

int file_is_obu(struct AvxInputContext *input_ctx) {
  // Parsing of Annex B OBU streams via pipe without framing not supported. This
  // implementation requires seeking backwards in the input stream. Tell caller
  // that this input cannot be processed.
  if (!input_ctx->filename || strcmp(input_ctx->filename, "-") == 0) return 0;

  uint8_t obutd[OBU_DETECTION_SIZE] = { 0 };

  const size_t ret = fread(obutd, 1, OBU_HEADER_SIZE, input_ctx->file);
  if (ret != OBU_HEADER_SIZE) {
    warn("Failed to read OBU Header, not enough data to process header.\n");
    return 0;
  }
  input_ctx->detect.buf_read = ret;

  const int obu_type = (obutd[0] >> 3) & 0xF;
  if (obu_type != OBU_TEMPORAL_DELIMITER) {
    warn("Invalid OBU type found while probing for OBU_TEMPORAL_DELIMITER.\n");
    return 0;
  }

  // TODO(tomfinegan): This block assumes the length field is one byte. It needs
  // to be updated to parse the extension and read/parse the length field.
  const int has_extension = obutd[0] & 0x1;
  if (has_extension) {
    // Skip the length and extension.
    fread(&obutd[1], 1, OBU_EXTENSION_SIZE + 1, input_ctx->file);
  } else {
    // Skip the length.
    fread(&obutd[1], 1, 1, input_ctx->file);
  }

  // fprintf(stderr, "Starting to parse OBU stream\n");
  return 1;
}

#endif
