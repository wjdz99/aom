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

// Log 2 conversion lookup tables in units of mode info (4x4).
// The Mi_Width_Log2 table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_wide_log2[BLOCK_SIZES_ALL] = {
  0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 0, 2, 1, 3, 2, 4
};
// The Mi_Height_Log2 table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_high_log2[BLOCK_SIZES_ALL] = {
  0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 2, 0, 3, 1, 4, 2
};

// Width/height lookup tables in units of mode info (4x4).
// The Num_4x4_Blocks_Wide table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_wide[BLOCK_SIZES_ALL] = {
  1, 1, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 1, 4, 2, 8, 4, 16
};

// The Num_4x4_Blocks_High table in the spec (Section 9.3. Conversion tables).
const uint8_t mi_size_high[BLOCK_SIZES_ALL] = {
  1, 2, 1, 2, 4, 2, 4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 4, 1, 8, 2, 16, 4
};

// Width/height lookup tables in units of samples.
// The Block_Width table in the spec (Section 9.3. Conversion tables).
const uint8_t block_size_wide[BLOCK_SIZES_ALL] = {
  4,  4,  8,  8,   8,   16, 16, 16, 32, 32, 32,
  64, 64, 64, 128, 128, 4,  16, 8,  32, 16, 64
};

// The Block_Height table in the spec (Section 9.3. Conversion tables).
const uint8_t block_size_high[BLOCK_SIZES_ALL] = {
  4,  8,  4,   8,  16,  8,  16, 32, 16, 32, 64,
  32, 64, 128, 64, 128, 16, 4,  32, 8,  64, 16
};

/* clang-format off */
const TX_SIZE max_txsize_rect_lookup[BLOCK_SIZES_ALL] = {
  // 4X4
  TX_4X4,
  // 4X8,    8X4,      8X8
  TX_4X8, TX_8X4, TX_8X8,
  // 8X16,   16X8,     16X16
  TX_8X16, TX_16X8, TX_16X16,
  // 16X32,  32X16,    32X32
  TX_16X32, TX_32X16, TX_32X32,
  // 32X64,  64X32,
  TX_32X64, TX_64X32,
  // 64X64
  TX_64X64,
  // 64x128, 128x64,   128x128
  TX_64X64, TX_64X64, TX_64X64,
  // 4x16,   16x4,
  TX_4X16, TX_16X4,
  // 8x32,   32x8
  TX_8X32, TX_32X8,
  // 16x64,  64x16
  TX_16X64, TX_64X16
};

const TX_TYPE_1D vtx_tab[TX_TYPES] = {
  DCT_1D,      ADST_1D, DCT_1D,      ADST_1D, FLIPADST_1D, DCT_1D,
  FLIPADST_1D, ADST_1D, FLIPADST_1D, IDTX_1D, DCT_1D,      IDTX_1D,
  ADST_1D,     IDTX_1D, FLIPADST_1D, IDTX_1D,
};

const TX_TYPE_1D htx_tab[TX_TYPES] = {
  DCT_1D,      DCT_1D,      ADST_1D, ADST_1D,     DCT_1D,  FLIPADST_1D,
  FLIPADST_1D, FLIPADST_1D, ADST_1D, IDTX_1D,     IDTX_1D, DCT_1D,
  IDTX_1D,     ADST_1D,     IDTX_1D, FLIPADST_1D,
};

/* clang-format on */

const TX_SIZE sub_tx_size_map[TX_SIZES_ALL] = {
  TX_4X4,    // TX_4X4
  TX_4X4,    // TX_8X8
  TX_8X8,    // TX_16X16
  TX_16X16,  // TX_32X32
  TX_32X32,  // TX_64X64
  TX_4X4,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_8X8,    // TX_16X8
  TX_16X16,  // TX_16X32
  TX_16X16,  // TX_32X16
  TX_32X32,  // TX_32X64
  TX_32X32,  // TX_64X32
  TX_4X8,    // TX_4X16
  TX_8X4,    // TX_16X4
  TX_8X16,   // TX_8X32
  TX_16X8,   // TX_32X8
  TX_16X32,  // TX_16X64
  TX_32X16,  // TX_64X16
};

const TX_SIZE txsize_horz_map[TX_SIZES_ALL] = {
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_64X64,  // TX_64X64
  TX_4X4,    // TX_4X8
  TX_8X8,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_16X16,  // TX_16X8
  TX_16X16,  // TX_16X32
  TX_32X32,  // TX_32X16
  TX_32X32,  // TX_32X64
  TX_64X64,  // TX_64X32
  TX_4X4,    // TX_4X16
  TX_16X16,  // TX_16X4
  TX_8X8,    // TX_8X32
  TX_32X32,  // TX_32X8
  TX_16X16,  // TX_16X64
  TX_64X64,  // TX_64X16
};

const TX_SIZE txsize_vert_map[TX_SIZES_ALL] = {
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_64X64,  // TX_64X64
  TX_8X8,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_16X16,  // TX_8X16
  TX_8X8,    // TX_16X8
  TX_32X32,  // TX_16X32
  TX_16X16,  // TX_32X16
  TX_64X64,  // TX_32X64
  TX_32X32,  // TX_64X32
  TX_16X16,  // TX_4X16
  TX_4X4,    // TX_16X4
  TX_32X32,  // TX_8X32
  TX_8X8,    // TX_32X8
  TX_64X64,  // TX_16X64
  TX_16X16,  // TX_64X16
};

// Transform block width in pixels
const int tx_size_wide[TX_SIZES_ALL] = {
  4, 8, 16, 32, 64, 4, 8, 8, 16, 16, 32, 32, 64, 4, 16, 8, 32, 16, 64,
};

// Transform block height in pixels
const int tx_size_high[TX_SIZES_ALL] = {
  4, 8, 16, 32, 64, 8, 4, 16, 8, 32, 16, 64, 32, 16, 4, 32, 8, 64, 16,
};

// Transform block width in unit
const int tx_size_wide_unit[TX_SIZES_ALL] = {
  1, 2, 4, 8, 16, 1, 2, 2, 4, 4, 8, 8, 16, 1, 4, 2, 8, 4, 16,
};

// Transform block height in unit
const int tx_size_high_unit[TX_SIZES_ALL] = {
  1, 2, 4, 8, 16, 2, 1, 4, 2, 8, 4, 16, 8, 4, 1, 8, 2, 16, 4,
};

// Transform block width in log2
const int tx_size_wide_log2[TX_SIZES_ALL] = {
  2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 5, 6, 2, 4, 3, 5, 4, 6,
};

// Transform block width in log2 unit
const int tx_size_wide_unit_log2[TX_SIZES_ALL] = {
  0, 1, 2, 3, 4, 0, 1, 1, 2, 2, 3, 3, 4, 0, 2, 1, 3, 2, 4,
};

// Transform block height in log2
const int tx_size_high_log2[TX_SIZES_ALL] = {
  2, 3, 4, 5, 6, 3, 2, 4, 3, 5, 4, 6, 5, 4, 2, 5, 3, 6, 4,
};

// Transform block height in log2 unit
const int tx_size_high_unit_log2[TX_SIZES_ALL] = {
  0, 1, 2, 3, 4, 1, 0, 2, 1, 3, 2, 4, 3, 2, 0, 3, 1, 4, 2,
};

const int tx_size_2d[TX_SIZES_ALL + 1] = {
  16,  64,   256,  1024, 4096, 32,  32,  128,  128,  512,
  512, 2048, 2048, 64,   64,   256, 256, 1024, 1024,
};

const TX_SIZE txsize_sqr_map[TX_SIZES_ALL] = {
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_64X64,  // TX_64X64
  TX_4X4,    // TX_4X8
  TX_4X4,    // TX_8X4
  TX_8X8,    // TX_8X16
  TX_8X8,    // TX_16X8
  TX_16X16,  // TX_16X32
  TX_16X16,  // TX_32X16
  TX_32X32,  // TX_32X64
  TX_32X32,  // TX_64X32
  TX_4X4,    // TX_4X16
  TX_4X4,    // TX_16X4
  TX_8X8,    // TX_8X32
  TX_8X8,    // TX_32X8
  TX_16X16,  // TX_16X64
  TX_16X16,  // TX_64X16
};

const TX_SIZE txsize_sqr_up_map[TX_SIZES_ALL] = {
  TX_4X4,    // TX_4X4
  TX_8X8,    // TX_8X8
  TX_16X16,  // TX_16X16
  TX_32X32,  // TX_32X32
  TX_64X64,  // TX_64X64
  TX_8X8,    // TX_4X8
  TX_8X8,    // TX_8X4
  TX_16X16,  // TX_8X16
  TX_16X16,  // TX_16X8
  TX_32X32,  // TX_16X32
  TX_32X32,  // TX_32X16
  TX_64X64,  // TX_32X64
  TX_64X64,  // TX_64X32
  TX_16X16,  // TX_4X16
  TX_16X16,  // TX_16X4
  TX_32X32,  // TX_8X32
  TX_32X32,  // TX_32X8
  TX_64X64,  // TX_16X64
  TX_64X64,  // TX_64X16
};

const int8_t txsize_log2_minus4[TX_SIZES_ALL] = {
  0,  // TX_4X4
  2,  // TX_8X8
  4,  // TX_16X16
  6,  // TX_32X32
  6,  // TX_64X64
  1,  // TX_4X8
  1,  // TX_8X4
  3,  // TX_8X16
  3,  // TX_16X8
  5,  // TX_16X32
  5,  // TX_32X16
  6,  // TX_32X64
  6,  // TX_64X32
  2,  // TX_4X16
  2,  // TX_16X4
  4,  // TX_8X32
  4,  // TX_32X8
  5,  // TX_16X64
  5,  // TX_64X16
};

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

const int intra_mode_context[INTRA_MODES] = {
  0, 1, 2, 3, 4, 4, 4, 4, 3, 0, 1, 2, 0,
};

// Note: this is also used in unit tests. So whenever one changes the table,
// the unit tests need to be changed accordingly.
const int quant_dist_weight[4][2] = {
  { 2, 3 }, { 2, 5 }, { 2, 7 }, { 1, MAX_FRAME_DISTANCE }
};

const int quant_dist_lookup_table[4][2] = {
  { 9, 7 },
  { 11, 5 },
  { 12, 4 },
  { 13, 3 },
};
