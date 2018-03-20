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

#define OBU_BUFFER_SIZE (500 * 1024)

#define OBU_HEADER_SIZE 1
#define OBU_EXTENSION_SIZE 1
#define OBU_MAX_LENGTH_FIELD_SIZE 8
#define OBU_DETECTION_SIZE \
  (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + OBU_MAX_LENGTH_FIELD_SIZE)

#if CONFIG_OBU_NO_IVF

// Reads unsigned LEB128 integer and returns 0 upon successful read and decode.
// Stores raw bytes in 'value_buffer', length of the number in 'value_length',
// and decoded value in 'value'.
static int obudec_read_leb128(FILE *f, uint8_t *value_buffer,
                              size_t *value_length, uint64_t *value) {
  if (!f || !value_buffer || !value_length || !value) return -1;
  for (int len = 0; len < OBU_MAX_LENGTH_FIELD_SIZE; ++len) {
    value_buffer[len] = fgetc(f);
    if ((value_buffer[len] >> 7) == 0) {
      *value_length = (size_t)(len + 1);
      break;
    }
  }

  return aom_uleb_decode(value_buffer, OBU_MAX_LENGTH_FIELD_SIZE, value, NULL);
}

// Reads OBU size from infile and returns 0 upon success. Returns obu_size via
// output pointer obu_size. Returns -1 when reading or parsing fails. Always
// returns FILE pointer to position at time of call. Returns 0 and sets obu_size
// to 0 when end of file is reached.
static int obudec_read_obu_size(FILE *infile, uint64_t *obu_size,
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

static int obudec_read_obu_header(FILE *f, size_t buffer_capacity,
                                  uint8_t *obu_data, ObuHeader *obu_header,
                                  size_t *bytes_read) {
  if (!f || buffer_capacity < (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE) ||
      !obu_data || !obu_header) {
    return -1;
  }
  *bytes_read = fread(obu_data, 1, 1, f);
  if (*bytes_read != 1) {
    warn("obudec: Failure reading OBU header.\n");
    return -1;
  }

  const int has_extension = obu_data[0] & 0x1;
  if (has_extension) {
    if (fread(&obu_data[1], 1, 1, f) != 1) {
      warn("obudec: Failure reading OBU extension.");
      return -1;
    }
    ++*bytes_read;
  }

  size_t obu_bytes_parsed = 0;
  aom_codec_err_t parse_result =
      aom_read_obu_header(obu_data, *bytes_read, &obu_bytes_parsed, obu_header);
  if (parse_result != AOM_CODEC_OK || *bytes_read != obu_bytes_parsed) {
    warn("obudec: Error parsing OBU header.\n");
    return -1;
  }

  return 0;
}

static int obudec_read_obu_payload(FILE *f, uint64_t payload_length,
                                   uint8_t *obu_data, size_t *bytes_read) {
  if (!f || payload_length == 0 || !obu_data || !bytes_read) return -1;

  if (fread(obu_data, 1, (size_t)payload_length, f) != payload_length) {
    warn("obudec: Failure reading OBU payload.\n");
    return -1;
  }

  *bytes_read += payload_length;
  return 0;
}

static int obudec_read_one_obu(FILE *f, size_t buffer_capacity,
                               uint8_t *obu_data, uint64_t *obu_length,
                               ObuHeader *obu_header) {
  const size_t kMinimumBufferSize = OBU_DETECTION_SIZE;
  if (!f || !obu_data || !obu_length || !obu_header ||
      buffer_capacity < kMinimumBufferSize) {
    return -1;
  }

  if (feof(f)) {
    *obu_length = 0;
    return 0;
  }

  size_t bytes_read = 0;
  if (obudec_read_obu_header(f, buffer_capacity, obu_data, obu_header,
                             &bytes_read) != 0)
    return -1;

  uint64_t obu_payload_length = 0;
  size_t leb128_length = 0;
  if (obudec_read_leb128(f, &obu_data[bytes_read], &leb128_length,
                         &obu_payload_length) != 0) {
    warn("obudec: Failure reading OBU payload length.\n");
    return -1;
  }
  bytes_read += leb128_length;

  if (bytes_read + obu_payload_length > buffer_capacity) {
    *obu_length = bytes_read + obu_payload_length;
    return -1;
  }

  if (obu_payload_length > 0 &&
      obudec_read_obu_payload(f, obu_payload_length, &obu_data[bytes_read],
                              &bytes_read)) {
    return -1;
  }

  *obu_length = bytes_read + obu_payload_length;
  return 0;
}

int file_is_obu(struct ObuDecInputContext *obu_ctx) {
  if (!obu_ctx || !obu_ctx->avx_ctx) return 0;

  struct AvxInputContext *avx_ctx = obu_ctx->avx_ctx;
  uint8_t detect_buf[OBU_DETECTION_SIZE] = { 0 };
  size_t detect_bytes_read = 0;

  FILE *f = avx_ctx->file;
  uint64_t obu_length = 0;
  ObuHeader obu_header;
  memset(&obu_header, 0, sizeof(obu_header));

  if (obudec_read_one_obu(f, OBU_DETECTION_SIZE, &detect_buf[0], &obu_length,
                   &obu_header) != 0) {
    warn("obudec: Failure reading first OBU.\n");
    return 0;
  }

  if (obu_header.type != OBU_TEMPORAL_DELIMITER) return 0;

  if (obu_header.has_length_field) {
    uint64_t obu_payload_length = 0;
    size_t leb128_length = 0;
    const size_t obu_length_offset = obu_header.has_length_field ? 1 : 2;
    if (aom_uleb_decode(&detect_buf[obu_length_offset], sizeof(leb128_length),
                        &obu_payload_length, &leb128_length) != 0) {
      warn("obudec: Failure decoding OBU payload length.\n");
      return 0;
    }
    if (obu_payload_length != 0) {
      warn("obudec: Invalid OBU_TEMPORAL_DELIMITER payload length (non-zero).");
      return 0;
    }
  }

  // Appears that input is valid Section 5 AV1 stream.
  obu_ctx->buffer = (uint8_t *)calloc(OBU_BUFFER_SIZE, 1);
  if (!obu_ctx->buffer) {
    warn("Out of memory.\n");
    return 0;
  }
  obu_ctx->buffer_capacity = OBU_BUFFER_SIZE;
  memcpy(obu_ctx->buffer, &detect_buf[0], detect_bytes_read);
  obu_ctx->bytes_buffered = detect_bytes_read;

  fprintf(stderr, "Starting to parse OBU stream (file_is_obu happy)\n");
  exit(1);
  return 1;
}


int obu_read_temporal_unit(struct ObuDecInputContext *obu_ctx, uint8_t **buffer,
                           size_t *bytes_read,
#if CONFIG_SCALABILITY
                           size_t *buffer_size, int last_layer_id
#else
                           size_t *buffer_size
#endif
) {
  size_t ret = 0;
  uint64_t obu_size = 0;
  uint8_t *data = NULL;
  FILE *infile = obu_ctx->avx_ctx->file;

  if (feof(infile)) {
    return 1;
  }

  *buffer_size = 0;
  *bytes_read = 0;
  while (1) {
    uint8_t obu_header_raw_bytes[OBU_HEADER_SIZE + OBU_EXTENSION_SIZE] = { 0 };
    //    size_t obu_header_size = 0;
    ObuHeader obu_header;
    memset(&obu_header, 0, sizeof(obu_header));

    // BROKEN: Must read bytes into data here.

    //    if (read_obu_header(infile, &obu_header_size, &obu_header,
    //                        &obu_header_raw_bytes[0]) != AOM_CODEC_OK) {
    //      fprintf(stderr, "Failed to read OBU header.\n");
    //      return 1;
    //    }

    size_t length_field_size = 0;
    //    ret = read_obu_size(infile, &obu_size, &length_field_size);
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


void obudec_free(struct ObuDecInputContext *obu_ctx) { free(obu_ctx->buffer); }

#endif  // CONFIG_OBU_NO_IVF
