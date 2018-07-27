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

// Header used by lightfield example.

#ifndef EXAMPLES_LIGHTFIELD_H_
#define EXAMPLES_LIGHTFIELD_H_

#define MAX_TILES 512

// Output frame size
const int output_frame_width = 512;
const int output_frame_height = 512;

// Spec:
// typedef struct {
//   uint8_t anchor_frame_idx;
//   uint8_t tile_row;
//   uint8_t tile_col;
//   uint16_t coded_tile_data_size_minus_1;
//   uint8_t *coded_tile_data;
// } TILE_LIST_ENTRY;

// Tile list entry provided by the application
typedef struct {
  int image_idx;
  int reference_idx;
  int tile_col;
  int tile_row;
} TILE_LIST_INFO;

// M references: 0 - M-1; N images(including references): 0 - N-1;
// Note: order the image index incrementally, so that we only go through the
// bitstream once to construct the tile list.
const int num_tile_lists = 2;
const uint16_t tile_count_minus_1[2] = { 9 - 1, 66 - 1 };
const TILE_LIST_INFO tile_list[2][MAX_TILES] = {
  // tile list 0
  { { 16, 0, 4, 5 },
    { 83, 3, 13, 2 },
    { 57, 2, 2, 6 },
    { 31, 1, 11, 5 },
    { 2, 0, 7, 4 },
    { 77, 3, 9, 9 },
    { 49, 1, 0, 1 },
    { 6, 0, 3, 10 },
    { 63, 2, 5, 8 } },
  // tile list 1
  { { 65, 2, 11, 1 }, { 42, 1, 3, 7 },  { 88, 3, 8, 4 },  { 76, 3, 1, 15 },
    { 1, 0, 2, 2 },   { 19, 0, 5, 6 },  { 60, 2, 4, 0 },  { 25, 1, 11, 15 },
    { 50, 2, 5, 4 },  { 3, 0, 1, 0 },   { 92, 3, 0, 11 }, { 65, 1, 13, 4 },
    { 81, 2, 2, 4 },  { 9, 2, 7, 14 },  { 23, 2, 6, 10 }, { 0, 3, 0, 0 },
    { 16, 3, 6, 2 },  { 35, 0, 8, 9 },  { 20, 1, 6, 3 },  { 64, 2, 15, 15 },
    { 13, 3, 0, 6 },  { 47, 0, 8, 5 },  { 82, 3, 12, 1 }, { 19, 2, 6, 15 },
    { 72, 1, 14, 2 }, { 67, 2, 5, 11 }, { 1, 1, 9, 3 },   { 5, 3, 10, 7 },
    { 14, 2, 11, 7 }, { 53, 2, 0, 14 }, { 95, 0, 9, 7 },  { 59, 1, 6, 6 },
    { 65, 2, 11, 1 }, { 42, 1, 3, 7 },  { 88, 3, 8, 4 },  { 76, 3, 1, 15 },
    { 1, 0, 2, 2 },   { 19, 0, 5, 6 },  { 60, 2, 4, 0 },  { 25, 1, 11, 15 },
    { 50, 2, 5, 4 },  { 3, 0, 1, 0 },   { 92, 3, 0, 11 }, { 65, 1, 13, 4 },
    { 81, 2, 2, 4 },  { 9, 2, 7, 14 },  { 23, 2, 6, 10 }, { 0, 3, 0, 0 },
    { 16, 3, 6, 2 },  { 35, 0, 8, 9 },  { 20, 1, 6, 3 },  { 64, 2, 15, 15 },
    { 13, 3, 0, 6 },  { 47, 0, 8, 5 },  { 82, 3, 12, 1 }, { 19, 2, 6, 15 },
    { 72, 1, 14, 2 }, { 67, 2, 5, 11 }, { 1, 1, 9, 3 },   { 5, 3, 10, 7 },
    { 14, 2, 11, 7 }, { 53, 2, 0, 14 }, { 95, 0, 9, 7 },  { 59, 1, 6, 6 },
    { 8, 0, 12, 9 },  { 37, 3, 0, 12 } },
};

#endif  // EXAMPLES_LIGHTFIELD_H_
