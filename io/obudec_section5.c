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
#include "io/obudec_section5.h"

#include "av1/decoder/obu.h"

static int obudec_read_one_obu(FILE *f, size_t buffer_capacity,
                               uint8_t *obu_data, uint64_t *obu_length,
                               ObuHeader *obu_header) {
  const size_t kMinimumBufferSize = OBU_DETECTION_SIZE;
  if (!f || !obu_data || !obu_length || !obu_header ||
      buffer_capacity < kMinimumBufferSize) {
    return -1;
  }

  uint64_t obu_payload_length = 0;
  size_t leb128_length = 0;

  size_t bytes_read = 0;
  if (obudec_read_obu_header(f, buffer_capacity, 0 /* is_annexb */,
                             obu_data + leb128_length, obu_header,
                             &bytes_read) != 0) {
    return -1;
  } else if (bytes_read == 0) {
    *obu_length = 0;
    return 0;
  }

  if (obudec_read_leb128(f, &obu_data[bytes_read], &leb128_length,
                         &obu_payload_length) != 0) {
    fprintf(stderr, "obudec: Failure reading OBU payload length.\n");
    return -1;
  }

  bytes_read += leb128_length;

  if (UINT64_MAX - bytes_read < obu_payload_length) return -1;
  if (bytes_read + obu_payload_length > buffer_capacity) {
    *obu_length = bytes_read + obu_payload_length;
    return -1;
  }

  if (obu_payload_length > 0 &&
      obudec_read_obu_payload(f, obu_payload_length, &obu_data[bytes_read],
                              &bytes_read) != 0) {
    return -1;
  }

  *obu_length = bytes_read;
  return 0;
}

int file_is_obu_section5(struct ObuDecInputContext *obu_ctx) {
  if (!obu_ctx || !obu_ctx->avx_ctx) return 0;

  struct AvxInputContext *avx_ctx = obu_ctx->avx_ctx;
  uint8_t detect_buf[OBU_DETECTION_SIZE] = { 0 };
  FILE *f = avx_ctx->file;
  uint64_t obu_length = 0;
  ObuHeader obu_header;
  memset(&obu_header, 0, sizeof(obu_header));

  if (obudec_read_one_obu(f, OBU_DETECTION_SIZE, &detect_buf[0], &obu_length,
                          &obu_header) != 0) {
    fprintf(stderr, "obudec: Failure reading first OBU.\n");
    rewind(f);
    return 0;
  }

  if (obu_header.type != OBU_TEMPORAL_DELIMITER) return 0;

  if (obu_header.has_length_field) {
    uint64_t obu_payload_length = 0;
    size_t leb128_length = 0;
    const size_t obu_length_offset = obu_header.has_length_field ? 1 : 2;
    if (aom_uleb_decode(&detect_buf[obu_length_offset], sizeof(leb128_length),
                        &obu_payload_length, &leb128_length) != 0) {
      fprintf(stderr, "obudec: Failure decoding OBU payload length.\n");
      rewind(f);
      return 0;
    }
    if (obu_payload_length != 0) {
      fprintf(
          stderr,
          "obudec: Invalid OBU_TEMPORAL_DELIMITER payload length (non-zero).");
      rewind(f);
      return 0;
    }
  } else {
    fprintf(stderr, "obudec: OBU size fields required, cannot decode input.\n");
    rewind(f);
    return (0);
  }

  // Appears that input is valid Section 5 AV1 stream.
  obu_ctx->buffer = (uint8_t *)calloc(OBU_BUFFER_SIZE, 1);
  if (!obu_ctx->buffer) {
    fprintf(stderr, "Out of memory.\n");
    rewind(f);
    return 0;
  }
  obu_ctx->buffer_capacity = OBU_BUFFER_SIZE;

  memcpy(obu_ctx->buffer, &detect_buf[0], (size_t)obu_length);
  obu_ctx->bytes_buffered = (size_t)obu_length;

  return 1;
}


