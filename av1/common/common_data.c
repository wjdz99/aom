/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/common/common_data.h"

// Width/height lookup tables in units of mode info (4x4).
// The Num_4x4_Blocks_Wide table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_wide[BLOCK_SIZES_ALL] = { 1, 1, 2, 2,  2,  4,  4,  4,
                                                8, 8, 8, 16, 16, 16, 32, 32,
                                                1, 4, 2, 8,  4,  16 };

// The Num_4x4_Blocks_High table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_high[BLOCK_SIZES_ALL] = { 1, 2, 1,  2, 4,  2,  4,  8,
                                                4, 8, 16, 8, 16, 32, 16, 32,
                                                4, 1, 8,  2, 16, 4 };

// The Subsampled_Size table in the spec (Section 5.11.38. Get plane residual
// size function).
/* clang-format off */
const BLOCK_SIZE av1_ss_size_lookup[BLOCK_SIZES_ALL][2][2] = {
  //  ss_x == 0      ss_x == 0          ss_x == 1      ss_x == 1
  //  ss_y == 0      ss_y == 1          ss_y == 0      ss_y == 1
  { { BLOCK_4X4,     BLOCK_4X4 },     { BLOCK_4X4,     BLOCK_4X4 } },
  { { BLOCK_4X8,     BLOCK_4X4 },     { BLOCK_INVALID, BLOCK_4X4 } },
  { { BLOCK_8X4,     BLOCK_INVALID }, { BLOCK_4X4,     BLOCK_4X4 } },
  { { BLOCK_8X8,     BLOCK_8X4 },     { BLOCK_4X8,     BLOCK_4X4 } },
  { { BLOCK_8X16,    BLOCK_8X8 },     { BLOCK_INVALID, BLOCK_4X8 } },
  { { BLOCK_16X8,    BLOCK_INVALID }, { BLOCK_8X8,     BLOCK_8X4 } },
  { { BLOCK_16X16,   BLOCK_16X8 },    { BLOCK_8X16,    BLOCK_8X8 } },
  { { BLOCK_16X32,   BLOCK_16X16 },   { BLOCK_INVALID, BLOCK_8X16 } },
  { { BLOCK_32X16,   BLOCK_INVALID }, { BLOCK_16X16,   BLOCK_16X8 } },
  { { BLOCK_32X32,   BLOCK_32X16 },   { BLOCK_16X32,   BLOCK_16X16 } },
  { { BLOCK_32X64,   BLOCK_32X32 },   { BLOCK_INVALID, BLOCK_16X32 } },
  { { BLOCK_64X32,   BLOCK_INVALID }, { BLOCK_32X32,   BLOCK_32X16 } },
  { { BLOCK_64X64,   BLOCK_64X32 },   { BLOCK_32X64,   BLOCK_32X32 } },
  { { BLOCK_64X128,  BLOCK_64X64 },   { BLOCK_INVALID, BLOCK_32X64 } },
  { { BLOCK_128X64,  BLOCK_INVALID }, { BLOCK_64X64,   BLOCK_64X32 } },
  { { BLOCK_128X128, BLOCK_128X64 },  { BLOCK_64X128,  BLOCK_64X64 } },
  { { BLOCK_4X16,    BLOCK_4X8 },     { BLOCK_INVALID, BLOCK_4X8 } },
  { { BLOCK_16X4,    BLOCK_INVALID }, { BLOCK_8X4,     BLOCK_8X4 } },
  { { BLOCK_8X32,    BLOCK_8X16 },    { BLOCK_INVALID, BLOCK_4X16 } },
  { { BLOCK_32X8,    BLOCK_INVALID }, { BLOCK_16X8,    BLOCK_16X4 } },
  { { BLOCK_16X64,   BLOCK_16X32 },   { BLOCK_INVALID, BLOCK_8X32 } },
  { { BLOCK_64X16,   BLOCK_INVALID }, { BLOCK_32X16,   BLOCK_32X8 } }
};
/* clang-format on */
