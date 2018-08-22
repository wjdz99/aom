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
#ifndef COMMON_AV1C_H_
#define COMMON_AV1C_H_

#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

// Struct representing ISOBMFF/Matroska AV1 config. The AV1 config has the
// following format:
//
// unsigned int (1) marker = 1;
// unsigned int (7) version = 1;
// unsigned int (3) seq_profile;
// unsigned int (5) seq_level_idx_0;
// unsigned int (1) seq_tier_0;
// unsigned int (1) high_bitdepth;
// unsigned int (1) twelve_bit;
// unsigned int (1) monochrome;
// unsigned int (1) chroma_subsampling_x;
// unsigned int (1) chroma_subsampling_y;
// unsigned int (2) chroma_sample_position;
// unsigned int (3) reserved = 0;
//
// unsigned int (1) initial_presentation_delay_present;
// if (initial_presentation_delay_present) {
//   unsigned int (4) initial_presentation_delay_minus_one;
// } else {
//   unsigned int (4) reserved = 0;
// }
//
// unsigned int (8)[] configOBUs;
//
// Note: get_av1c() does not currently store 'configOBUs' data, so the field is
// omitted.
struct av1c {
  uint8_t marker;
  uint8_t version;
  uint8_t seq_profile;
  uint8_t seq_level_idx_0;
  uint8_t seq_tier_0;
  uint8_t high_bitdepth;
  uint8_t twelve_bit;
  uint8_t monochrome;
  uint8_t chroma_subsampling_x;
  uint8_t chroma_subsampling_y;
  uint8_t chroma_sample_position;
  uint8_t initial_presentation_delay_present;
  uint8_t initial_presentation_delay_minus_one;
};

// Attempts to parse a Sequence Header OBU and set the paramenters of 'config'.
// Returns 0 upon success, and -1 upon failure. 'buffer' can contain multiple
// OBUs, but the Sequence Header OBU must be the first OBU within the buffer.
int get_av1c(const uint8_t *const buffer, size_t length, int is_annexb,
             struct av1c *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // COMMON_AV1C_H_
