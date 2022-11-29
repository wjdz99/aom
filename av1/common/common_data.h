/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_COMMON_COMMON_DATA_H_
#define AOM_AV1_COMMON_COMMON_DATA_H_

#include "av1/common/enums.h"
#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"

#ifdef __cplusplus
extern "C" {
#endif

// Log 2 conversion lookup tables in units of mode info (4x4).
// The Mi_Width_Log2 table in the spec (Section 9.3. Conversion tables).
extern const uint8_t mi_size_wide_log2[BLOCK_SIZES_ALL];

// The Mi_Height_Log2 table in the spec (Section 9.3. Conversion tables).
extern const uint8_t mi_size_high_log2[BLOCK_SIZES_ALL];

// Width/height lookup tables in units of mode info (4x4).
// The Num_4x4_Blocks_Wide table in the spec (Section 9.3. Conversion tables).
extern const uint8_t mi_size_wide[BLOCK_SIZES_ALL];

// The Num_4x4_Blocks_High table in the spec (Section 9.3. Conversion tables).
extern const uint8_t mi_size_high[BLOCK_SIZES_ALL];

// Width/height lookup tables in units of samples.
// The Block_Width table in the spec (Section 9.3. Conversion tables).
extern const uint8_t block_size_wide[BLOCK_SIZES_ALL];

// The Block_Height table in the spec (Section 9.3. Conversion tables).
extern const uint8_t block_size_high[BLOCK_SIZES_ALL];

// Maps a block size to a context.
// The Size_Group table in the spec (Section 9.3. Conversion tables).
// AOMMIN(3, AOMMIN(mi_size_wide_log2(bsize), mi_size_high_log2(bsize)))
static const uint8_t size_group_lookup[BLOCK_SIZES_ALL] = {
  0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 2, 2
};

static const uint8_t num_pels_log2_lookup[BLOCK_SIZES_ALL] = {
  4, 5, 5, 6, 7, 7, 8, 9, 9, 10, 11, 11, 12, 13, 13, 14, 6, 6, 8, 8, 10, 10
};

// A compressed version of the Partition_Subsize table in the spec (9.3.
// Conversion tables), for square block sizes only.
/* clang-format off */
static const BLOCK_SIZE subsize_lookup[EXT_PARTITION_TYPES][SQR_BLOCK_SIZES] = {
  {     // PARTITION_NONE
    BLOCK_4X4, BLOCK_8X8, BLOCK_16X16,
    BLOCK_32X32, BLOCK_64X64, BLOCK_128X128
  }, {  // PARTITION_HORZ
    BLOCK_INVALID, BLOCK_8X4, BLOCK_16X8,
    BLOCK_32X16, BLOCK_64X32, BLOCK_128X64
  }, {  // PARTITION_VERT
    BLOCK_INVALID, BLOCK_4X8, BLOCK_8X16,
    BLOCK_16X32, BLOCK_32X64, BLOCK_64X128
  }, {  // PARTITION_SPLIT
    BLOCK_INVALID, BLOCK_4X4, BLOCK_8X8,
    BLOCK_16X16, BLOCK_32X32, BLOCK_64X64
  }, {  // PARTITION_HORZ_A
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    BLOCK_32X16, BLOCK_64X32, BLOCK_128X64
  }, {  // PARTITION_HORZ_B
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    BLOCK_32X16, BLOCK_64X32, BLOCK_128X64
  }, {  // PARTITION_VERT_A
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
    BLOCK_16X32, BLOCK_32X64, BLOCK_64X128
  }, {  // PARTITION_VERT_B
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
    BLOCK_16X32, BLOCK_32X64, BLOCK_64X128
  }, {  // PARTITION_HORZ_4
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X4,
    BLOCK_32X8, BLOCK_64X16, BLOCK_INVALID
  }, {  // PARTITION_VERT_4
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X16,
    BLOCK_8X32, BLOCK_16X64, BLOCK_INVALID
  }
};

static const TX_SIZE max_txsize_lookup[BLOCK_SIZES_ALL] = {
  //                   4X4
                       TX_4X4,
  // 4X8,    8X4,      8X8
  TX_4X4,    TX_4X4,   TX_8X8,
  // 8X16,   16X8,     16X16
  TX_8X8,    TX_8X8,   TX_16X16,
  // 16X32,  32X16,    32X32
  TX_16X16,  TX_16X16, TX_32X32,
  // 32X64,  64X32,
  TX_32X32,  TX_32X32,
  // 64X64
  TX_64X64,
  // 64x128, 128x64,   128x128
  TX_64X64,  TX_64X64, TX_64X64,
  // 4x16,   16x4,     8x32
  TX_4X4,    TX_4X4,   TX_8X8,
  // 32x8,   16x64     64x16
  TX_8X8,    TX_16X16, TX_16X16
};

extern const TX_SIZE max_txsize_rect_lookup[BLOCK_SIZES_ALL];

extern const TX_TYPE_1D vtx_tab[TX_TYPES];

extern const TX_TYPE_1D htx_tab[TX_TYPES];

#define TXSIZE_CAT_INVALID (-1)

/* clang-format on */

extern const TX_SIZE sub_tx_size_map[TX_SIZES_ALL];

extern const TX_SIZE txsize_horz_map[TX_SIZES_ALL];

extern const TX_SIZE txsize_vert_map[TX_SIZES_ALL];

#define TX_SIZE_W_MIN 4

// Transform block width in pixels
extern const int tx_size_wide[TX_SIZES_ALL];

#define TX_SIZE_H_MIN 4

// Transform block height in pixels
extern const int tx_size_high[TX_SIZES_ALL];

// Transform block width in unit
extern const int tx_size_wide_unit[TX_SIZES_ALL];

// Transform block height in unit
extern const int tx_size_high_unit[TX_SIZES_ALL];

// Transform block width in log2
extern const int tx_size_wide_log2[TX_SIZES_ALL];

// Transform block width in log2 unit
extern const int tx_size_wide_unit_log2[TX_SIZES_ALL];

// Transform block height in log2
extern const int tx_size_high_log2[TX_SIZES_ALL];

// Transform block height in log2 unit
extern const int tx_size_high_unit_log2[TX_SIZES_ALL];

extern const int tx_size_2d[TX_SIZES_ALL + 1];

static const BLOCK_SIZE txsize_to_bsize[TX_SIZES_ALL] = {
  BLOCK_4X4,    // TX_4X4
  BLOCK_8X8,    // TX_8X8
  BLOCK_16X16,  // TX_16X16
  BLOCK_32X32,  // TX_32X32
  BLOCK_64X64,  // TX_64X64
  BLOCK_4X8,    // TX_4X8
  BLOCK_8X4,    // TX_8X4
  BLOCK_8X16,   // TX_8X16
  BLOCK_16X8,   // TX_16X8
  BLOCK_16X32,  // TX_16X32
  BLOCK_32X16,  // TX_32X16
  BLOCK_32X64,  // TX_32X64
  BLOCK_64X32,  // TX_64X32
  BLOCK_4X16,   // TX_4X16
  BLOCK_16X4,   // TX_16X4
  BLOCK_8X32,   // TX_8X32
  BLOCK_32X8,   // TX_32X8
  BLOCK_16X64,  // TX_16X64
  BLOCK_64X16,  // TX_64X16
};

extern const TX_SIZE txsize_sqr_map[TX_SIZES_ALL];

extern const TX_SIZE txsize_sqr_up_map[TX_SIZES_ALL];

extern const int8_t txsize_log2_minus4[TX_SIZES_ALL];

static const TX_SIZE tx_mode_to_biggest_tx_size[TX_MODES] = {
  TX_4X4,    // ONLY_4X4
  TX_64X64,  // TX_MODE_LARGEST
  TX_64X64,  // TX_MODE_SELECT
};

// The Subsampled_Size table in the spec (Section 5.11.38. Get plane residual
// size function).
extern const BLOCK_SIZE av1_ss_size_lookup[BLOCK_SIZES_ALL][2][2];

// Generates 5 bit field in which each bit set to 1 represents
// a blocksize partition  11111 means we split 128x128, 64x64, 32x32, 16x16
// and 8x8.  10000 means we just split the 128x128 to 64x64
/* clang-format off */
static const struct {
  PARTITION_CONTEXT above;
  PARTITION_CONTEXT left;
} partition_context_lookup[BLOCK_SIZES_ALL] = {
  { 31, 31 },  // 4X4   - {0b11111, 0b11111}
  { 31, 30 },  // 4X8   - {0b11111, 0b11110}
  { 30, 31 },  // 8X4   - {0b11110, 0b11111}
  { 30, 30 },  // 8X8   - {0b11110, 0b11110}
  { 30, 28 },  // 8X16  - {0b11110, 0b11100}
  { 28, 30 },  // 16X8  - {0b11100, 0b11110}
  { 28, 28 },  // 16X16 - {0b11100, 0b11100}
  { 28, 24 },  // 16X32 - {0b11100, 0b11000}
  { 24, 28 },  // 32X16 - {0b11000, 0b11100}
  { 24, 24 },  // 32X32 - {0b11000, 0b11000}
  { 24, 16 },  // 32X64 - {0b11000, 0b10000}
  { 16, 24 },  // 64X32 - {0b10000, 0b11000}
  { 16, 16 },  // 64X64 - {0b10000, 0b10000}
  { 16, 0 },   // 64X128- {0b10000, 0b00000}
  { 0, 16 },   // 128X64- {0b00000, 0b10000}
  { 0, 0 },    // 128X128-{0b00000, 0b00000}
  { 31, 28 },  // 4X16  - {0b11111, 0b11100}
  { 28, 31 },  // 16X4  - {0b11100, 0b11111}
  { 30, 24 },  // 8X32  - {0b11110, 0b11000}
  { 24, 30 },  // 32X8  - {0b11000, 0b11110}
  { 28, 16 },  // 16X64 - {0b11100, 0b10000}
  { 16, 28 },  // 64X16 - {0b10000, 0b11100}
};
/* clang-format on */

extern const int intra_mode_context[INTRA_MODES];

// Note: this is also used in unit tests. So whenever one changes the table,
// the unit tests need to be changed accordingly.
extern const int quant_dist_weight[4][2];

extern const int quant_dist_lookup_table[4][2];

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_COMMON_DATA_H_
