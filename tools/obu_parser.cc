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

#include <cstdio>

// TODO(tomfinegan): Remove the aom_config.h include as soon as possible. At
// present it's required because w/out aom_config.h it's impossible to know how
// Obus are being written out by the library (because
// CONFIG_ADD_4BYTES_OBUSIZE).
#include "./aom_config.h"
#include "aom_ports/mem_ops.h"
#include "tools/obu_parser.h"

namespace {
  const aom_tools::ObuExtensionHeader kEmptyObuExt = { 0, 0, 0, false };
  const aom_tools::ObuHeader kEmptyObu = { 0, 0, false, kEmptyObuExt };
}  // namespace

namespace aom_tools {

// Basic OBU syntax
// 4 bytes: length
// 8 bits: Header
//   7
//     forbidden bit
//   6,5,4,3
//     type bits
//   2,1
//     reserved bits
//   0
//     extension bit
const uint32_t kObuForbiddenBit = 0x80;
const uint32_t kObuForbiddenBitShift = 7;
const uint32_t kObuTypeBits = 0x78;
const uint32_t kObuTypeBitsShift = 3;
const uint32_t kObuReservedBits = 0x6;
const uint32_t kObuReservedBitsShift = 1;
const uint32_t kObuExtensionFlagBit = 0x1;

// When extension bit is set:
// 8 bits: extension header
// 7,6,5
//   temporal ID
// 4,3
//   spatial ID
// 2,1
//   quality ID
// 0
//   reserved bit
const uint32_t kObuExtTemporalIdBits = 0xE0;
const uint32_t kObuExtTemporalIdBitsShift = 5;
const uint32_t kObuExtSpatialIdBits = 0x18;
const uint32_t kObuExtSpatialIdBitsShift = 3;
const uint32_t kObuExtQualityIdBits = 0x6;
const uint32_t kObuExtQualityIdBitsShift = 1;
const uint32_t kObuExtReservedFlagBit = 0x1;

bool valid_obu_type(int obu_type) {
  switch(obu_type) {
    case OBU_SEQUENCE_HEADER:
    case OBU_TEMPORAL_DELIMITER:
    case OBU_FRAME_HEADER:
    case OBU_TILE_GROUP:
    case OBU_METADATA:
    case OBU_PADDING:
      return true;
  }
  return false;
}

bool parse_obu_header(uint8_t obu_header_byte, ObuHeader *obu_header) {
  bool forbidden_bit =
      (obu_header_byte & kObuForbiddenBit) >> kObuForbiddenBitShift;
  if (forbidden_bit) {
    fprintf(stderr, "Invalid OBU, forbidden bit set.\n");
    return false;
  }

  obu_header->type = (obu_header_byte & kObuTypeBits) >> kObuTypeBitsShift;
  if (!valid_obu_type(obu_header->type)) {
    fprintf(stderr, "Invalid OBU type: %d.\n", obu_header->type);
    return false;
  }

  obu_header->reserved =
      (obu_header_byte & kObuReservedBits) >> kObuReservedBitsShift;
  if (obu_header->reserved != 0) {
    fprintf(stderr, "Invalid OBU type: reserved bit(s) set.\n");
    return false;
  }

  obu_header->has_extension = obu_header_byte & kObuExtensionFlagBit;
  return true;
}

bool parse_obu_extension_header(uint8_t ext_header_byte,
                                ObuExtensionHeader *ext_header) {
  return false;
}
  
void dump_obu(const uint8_t *data, int length) {
  const int kMinimumBytesRequired = 5;  // Assumes CONFIG_ADD_4BYTES_OBUSIZE.
  int consumed = 0;
  ObuHeader obu_header;
  while (consumed < length) {
    const int remaining = length - consumed;
    if (remaining < kMinimumBytesRequired) {
      if (remaining > 0) {
        fprintf(stderr, "OBU parse error. Did not consume all data, %d bytes "
                "remain.\n", remaining);
      }
      break;
    }

    const int current_obu_length = mem_get_le32(data);
    if (current_obu_length > remaining) {
      fprintf(stderr, "OBU parsing failed at offset %d with bad length of %d "
              "and %d bytes left.\n", consumed, current_obu_length, remaining);
      break;
    }
    consumed += 4;

    obu_header = kEmptyObu;
    const uint8_t obu_header_byte = *(data + consumed);
    if (!parse_obu_header(obu_header_byte, &obu_header)) {
      fprintf(stderr, "OBU parsing failed at offset %d.\n", consumed);
      break;
    }

    if (obu_header.has_extension) {
      const uint8_t obu_ext_header_byte = *(data + consumed + 1);
      if (!parse_obu_extension_header(obu_ext_header_byte,
                                      &obu_header.ext_header)) {
        fprintf(stderr, "OBU extension parsing failed at offset %d.\n",
                consumed);
        break;
      }
    }

    // TODO(tomfinegan): Parse OBU payload. For now just consume it.
    consumed += current_obu_length;
  }
}

}  // namespace aom_tools

