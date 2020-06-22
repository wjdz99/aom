/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom_ports/system_state.h"

#include "config/aom_dsp_rtcd.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/encoder_utils.h"
#include "av1/encoder/rc_utils.h"

#define MIN_BOOST_COMBINE_FACTOR 4.0
#define MAX_BOOST_COMBINE_FACTOR 12.0

const int default_tx_type_probs[FRAME_UPDATE_TYPES][TX_SIZES_ALL][TX_TYPES] = {
  { { 221, 189, 214, 292, 0, 0, 0, 0, 0, 2, 38, 68, 0, 0, 0, 0 },
    { 262, 203, 216, 239, 0, 0, 0, 0, 0, 1, 37, 66, 0, 0, 0, 0 },
    { 315, 231, 239, 226, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 222, 188, 214, 287, 0, 0, 0, 0, 0, 2, 50, 61, 0, 0, 0, 0 },
    { 256, 182, 205, 282, 0, 0, 0, 0, 0, 2, 21, 76, 0, 0, 0, 0 },
    { 281, 214, 217, 222, 0, 0, 0, 0, 0, 1, 48, 41, 0, 0, 0, 0 },
    { 263, 194, 225, 225, 0, 0, 0, 0, 0, 2, 15, 100, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 170, 192, 242, 293, 0, 0, 0, 0, 0, 1, 68, 58, 0, 0, 0, 0 },
    { 199, 210, 213, 291, 0, 0, 0, 0, 0, 1, 14, 96, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 106, 69, 107, 278, 9, 15, 20, 45, 49, 23, 23, 88, 36, 74, 25, 57 },
    { 105, 72, 81, 98, 45, 49, 47, 50, 56, 72, 30, 81, 33, 95, 27, 83 },
    { 211, 105, 109, 120, 57, 62, 43, 49, 52, 58, 42, 116, 0, 0, 0, 0 },
    { 1008, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 131, 57, 98, 172, 19, 40, 37, 64, 69, 22, 41, 52, 51, 77, 35, 59 },
    { 176, 83, 93, 202, 22, 24, 28, 47, 50, 16, 12, 93, 26, 76, 17, 59 },
    { 136, 72, 89, 95, 46, 59, 47, 56, 61, 68, 35, 51, 32, 82, 26, 69 },
    { 122, 80, 87, 105, 49, 47, 46, 46, 57, 52, 13, 90, 19, 103, 15, 93 },
    { 1009, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0 },
    { 1011, 0, 0, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 202, 20, 84, 114, 14, 60, 41, 79, 99, 21, 41, 15, 50, 84, 34, 66 },
    { 196, 44, 23, 72, 30, 22, 28, 57, 67, 13, 4, 165, 15, 148, 9, 131 },
    { 882, 0, 0, 0, 0, 0, 0, 0, 0, 142, 0, 0, 0, 0, 0, 0 },
    { 840, 0, 0, 0, 0, 0, 0, 0, 0, 184, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 } },
  { { 213, 110, 141, 269, 12, 16, 15, 19, 21, 11, 38, 68, 22, 29, 16, 24 },
    { 216, 119, 128, 143, 38, 41, 26, 30, 31, 30, 42, 70, 23, 36, 19, 32 },
    { 367, 149, 154, 154, 38, 35, 17, 21, 21, 10, 22, 36, 0, 0, 0, 0 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 219, 96, 127, 191, 21, 40, 25, 32, 34, 18, 45, 45, 33, 39, 26, 33 },
    { 296, 99, 122, 198, 23, 21, 19, 24, 25, 13, 20, 64, 23, 32, 18, 27 },
    { 275, 128, 142, 143, 35, 48, 23, 30, 29, 18, 42, 36, 18, 23, 14, 20 },
    { 239, 132, 166, 175, 36, 27, 19, 21, 24, 14, 13, 85, 9, 31, 8, 25 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 309, 25, 79, 59, 25, 80, 34, 53, 61, 25, 49, 23, 43, 64, 36, 59 },
    { 270, 57, 40, 54, 50, 42, 41, 53, 56, 28, 17, 81, 45, 86, 34, 70 },
    { 1005, 0, 0, 0, 0, 0, 0, 0, 0, 19, 0, 0, 0, 0, 0, 0 },
    { 992, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 133, 63, 55, 83, 57, 87, 58, 72, 68, 16, 24, 35, 29, 105, 25, 114 },
    { 131, 75, 74, 60, 71, 77, 65, 66, 73, 33, 21, 79, 20, 83, 18, 78 },
    { 276, 95, 82, 58, 86, 93, 63, 60, 64, 17, 38, 92, 0, 0, 0, 0 },
    { 1006, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 147, 49, 75, 78, 50, 97, 60, 67, 76, 17, 42, 35, 31, 93, 27, 80 },
    { 157, 49, 58, 75, 61, 52, 56, 67, 69, 12, 15, 79, 24, 119, 11, 120 },
    { 178, 69, 83, 77, 69, 85, 72, 77, 77, 20, 35, 40, 25, 48, 23, 46 },
    { 174, 55, 64, 57, 73, 68, 62, 61, 75, 15, 12, 90, 17, 99, 16, 86 },
    { 1008, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0 },
    { 1018, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 266, 31, 63, 64, 21, 52, 39, 54, 63, 30, 52, 31, 48, 89, 46, 75 },
    { 272, 26, 32, 44, 29, 31, 32, 53, 51, 13, 13, 88, 22, 153, 16, 149 },
    { 923, 0, 0, 0, 0, 0, 0, 0, 0, 101, 0, 0, 0, 0, 0, 0 },
    { 969, 0, 0, 0, 0, 0, 0, 0, 0, 55, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 } },
  { { 158, 92, 125, 298, 12, 15, 20, 29, 31, 12, 29, 67, 34, 44, 23, 35 },
    { 147, 94, 103, 123, 45, 48, 38, 41, 46, 48, 37, 78, 33, 63, 27, 53 },
    { 268, 126, 125, 136, 54, 53, 31, 38, 38, 33, 35, 87, 0, 0, 0, 0 },
    { 1018, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 159, 72, 103, 194, 20, 35, 37, 50, 56, 21, 39, 40, 51, 61, 38, 48 },
    { 259, 86, 95, 188, 32, 20, 25, 34, 37, 13, 12, 85, 25, 53, 17, 43 },
    { 189, 99, 113, 123, 45, 59, 37, 46, 48, 44, 39, 41, 31, 47, 26, 37 },
    { 175, 110, 113, 128, 58, 38, 33, 33, 43, 29, 13, 100, 14, 68, 12, 57 },
    { 1017, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0 },
    { 1019, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 208, 22, 84, 101, 21, 59, 44, 70, 90, 25, 59, 13, 64, 67, 49, 48 },
    { 277, 52, 32, 63, 43, 26, 33, 48, 54, 11, 6, 130, 18, 119, 11, 101 },
    { 963, 0, 0, 0, 0, 0, 0, 0, 0, 61, 0, 0, 0, 0, 0, 0 },
    { 979, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } }
};

const int default_obmc_probs[FRAME_UPDATE_TYPES][BLOCK_SIZES_ALL] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  106, 90, 90, 97, 67, 59, 70, 28,
    30, 38, 16, 16,  16, 0,  0,  44, 50, 26, 25 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  98, 93, 97, 68, 82, 85, 33, 30,
    33, 16, 16, 16, 16, 0,  0,  43, 37, 26, 16 },
  { 0,  0,  0,  91, 80, 76, 78, 55, 49, 24, 16,
    16, 16, 16, 16, 16, 0,  0,  29, 45, 16, 38 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  103, 89, 89, 89, 62, 63, 76, 34,
    35, 32, 19, 16,  16, 0,  0,  49, 55, 29, 19 }
};

const int default_warped_probs[FRAME_UPDATE_TYPES] = { 64, 64, 64, 64,
                                                       64, 64, 64 };

// TODO(yunqing): the default probs can be trained later from better
// performance.
const int default_switchable_interp_probs[FRAME_UPDATE_TYPES]
                                         [SWITCHABLE_FILTER_CONTEXTS]
                                         [SWITCHABLE_FILTERS] = {
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } }
                                         };

#define HIGHBD_BFP(BT, SDF, SDAF, VF, SVF, SVAF, SDX4DF, JSDAF, JSVAF) \
  cpi->fn_ptr[BT].sdf = SDF;                                           \
  cpi->fn_ptr[BT].sdaf = SDAF;                                         \
  cpi->fn_ptr[BT].vf = VF;                                             \
  cpi->fn_ptr[BT].svf = SVF;                                           \
  cpi->fn_ptr[BT].svaf = SVAF;                                         \
  cpi->fn_ptr[BT].sdx4df = SDX4DF;                                     \
  cpi->fn_ptr[BT].jsdaf = JSDAF;                                       \
  cpi->fn_ptr[BT].jsvaf = JSVAF;

#define HIGHBD_BFP_WRAPPER(WIDTH, HEIGHT, BD)                                \
  HIGHBD_BFP(                                                                \
      BLOCK_##WIDTH##X##HEIGHT, aom_highbd_sad##WIDTH##x##HEIGHT##_bits##BD, \
      aom_highbd_sad##WIDTH##x##HEIGHT##_avg_bits##BD,                       \
      aom_highbd_##BD##_variance##WIDTH##x##HEIGHT,                          \
      aom_highbd_##BD##_sub_pixel_variance##WIDTH##x##HEIGHT,                \
      aom_highbd_##BD##_sub_pixel_avg_variance##WIDTH##x##HEIGHT,            \
      aom_highbd_sad##WIDTH##x##HEIGHT##x4d_bits##BD,                        \
      aom_highbd_dist_wtd_sad##WIDTH##x##HEIGHT##_avg_bits##BD,              \
      aom_highbd_##BD##_dist_wtd_sub_pixel_avg_variance##WIDTH##x##HEIGHT)

#define MAKE_BFP_SAD_WRAPPER(fnname)                                           \
  static unsigned int fnname##_bits8(const uint8_t *src_ptr,                   \
                                     int source_stride,                        \
                                     const uint8_t *ref_ptr, int ref_stride) { \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride);                \
  }                                                                            \
  static unsigned int fnname##_bits10(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride) {                                                        \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride) >> 2;           \
  }                                                                            \
  static unsigned int fnname##_bits12(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride) {                                                        \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride) >> 4;           \
  }

#define MAKE_BFP_SADAVG_WRAPPER(fnname)                                        \
  static unsigned int fnname##_bits8(                                          \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred);   \
  }                                                                            \
  static unsigned int fnname##_bits10(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred) >> \
           2;                                                                  \
  }                                                                            \
  static unsigned int fnname##_bits12(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred) >> \
           4;                                                                  \
  }

#define MAKE_BFP_SAD4D_WRAPPER(fnname)                                        \
  static void fnname##_bits8(const uint8_t *src_ptr, int source_stride,       \
                             const uint8_t *const ref_ptr[], int ref_stride,  \
                             unsigned int *sad_array) {                       \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
  }                                                                           \
  static void fnname##_bits10(const uint8_t *src_ptr, int source_stride,      \
                              const uint8_t *const ref_ptr[], int ref_stride, \
                              unsigned int *sad_array) {                      \
    int i;                                                                    \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
    for (i = 0; i < 4; i++) sad_array[i] >>= 2;                               \
  }                                                                           \
  static void fnname##_bits12(const uint8_t *src_ptr, int source_stride,      \
                              const uint8_t *const ref_ptr[], int ref_stride, \
                              unsigned int *sad_array) {                      \
    int i;                                                                    \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
    for (i = 0; i < 4; i++) sad_array[i] >>= 4;                               \
  }

#define MAKE_BFP_JSADAVG_WRAPPER(fnname)                                    \
  static unsigned int fnname##_bits8(                                       \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param);                                               \
  }                                                                         \
  static unsigned int fnname##_bits10(                                      \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param) >>                                             \
           2;                                                               \
  }                                                                         \
  static unsigned int fnname##_bits12(                                      \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param) >>                                             \
           4;                                                               \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad128x128)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad128x128_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad128x128x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad128x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad128x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad128x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x128)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x128_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x128x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x4x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x4x4d)

MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x4x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x16x4d)

MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad128x128_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad128x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x128_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x16_avg)
#endif  // CONFIG_AV1_HIGHBITDEPTH

#define HIGHBD_MBFP(BT, MCSDF, MCSVF) \
  cpi->fn_ptr[BT].msdf = MCSDF;       \
  cpi->fn_ptr[BT].msvf = MCSVF;

#define HIGHBD_MBFP_WRAPPER(WIDTH, HEIGHT, BD)                    \
  HIGHBD_MBFP(BLOCK_##WIDTH##X##HEIGHT,                           \
              aom_highbd_masked_sad##WIDTH##x##HEIGHT##_bits##BD, \
              aom_highbd_##BD##_masked_sub_pixel_variance##WIDTH##x##HEIGHT)

#define MAKE_MBFP_COMPOUND_SAD_WRAPPER(fnname)                           \
  static unsigned int fnname##_bits8(                                    \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask);            \
  }                                                                      \
  static unsigned int fnname##_bits10(                                   \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask) >>          \
           2;                                                            \
  }                                                                      \
  static unsigned int fnname##_bits12(                                   \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask) >>          \
           4;                                                            \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad128x128)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad128x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x128)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x16)
#endif

#define HIGHBD_OBFP(BT, OSDF, OVF, OSVF) \
  cpi->fn_ptr[BT].osdf = OSDF;           \
  cpi->fn_ptr[BT].ovf = OVF;             \
  cpi->fn_ptr[BT].osvf = OSVF;

#define HIGHBD_OBFP_WRAPPER(WIDTH, HEIGHT, BD)                   \
  HIGHBD_OBFP(BLOCK_##WIDTH##X##HEIGHT,                          \
              aom_highbd_obmc_sad##WIDTH##x##HEIGHT##_bits##BD,  \
              aom_highbd_##BD##_obmc_variance##WIDTH##x##HEIGHT, \
              aom_highbd_##BD##_obmc_sub_pixel_variance##WIDTH##x##HEIGHT)

#define LOWBD_OBFP_WRAPPER(WIDTH, HEIGHT)                    \
  HIGHBD_OBFP(BLOCK_##WIDTH##X##HEIGHT,                      \
              aom_highbd_obmc_sad##WIDTH##x##HEIGHT##_bits8, \
              aom_highbd_obmc_variance##WIDTH##x##HEIGHT,    \
              aom_highbd_obmc_sub_pixel_variance##WIDTH##x##HEIGHT)

#define MAKE_OBFP_SAD_WRAPPER(fnname)                                     \
  static unsigned int fnname##_bits8(const uint8_t *ref, int ref_stride,  \
                                     const int32_t *wsrc,                 \
                                     const int32_t *msk) {                \
    return fnname(ref, ref_stride, wsrc, msk);                            \
  }                                                                       \
  static unsigned int fnname##_bits10(const uint8_t *ref, int ref_stride, \
                                      const int32_t *wsrc,                \
                                      const int32_t *msk) {               \
    return fnname(ref, ref_stride, wsrc, msk) >> 2;                       \
  }                                                                       \
  static unsigned int fnname##_bits12(const uint8_t *ref, int ref_stride, \
                                      const int32_t *wsrc,                \
                                      const int32_t *msk) {               \
    return fnname(ref, ref_stride, wsrc, msk) >> 4;                       \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad128x128)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad128x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x128)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x16)

void highbd_set_var_fns(AV1_COMP *const cpi) {
  AV1_COMMON *const cm = &cpi->common;
  if (cm->seq_params.use_highbitdepth) {
    switch (cm->seq_params.bit_depth) {
      case AOM_BITS_8:
        HIGHBD_BFP_WRAPPER(64, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 64, 8)
        HIGHBD_BFP_WRAPPER(32, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 32, 8)
        HIGHBD_BFP_WRAPPER(16, 4, 8)
        HIGHBD_BFP_WRAPPER(4, 16, 8)
        HIGHBD_BFP_WRAPPER(32, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 32, 8)
        HIGHBD_BFP_WRAPPER(64, 32, 8)
        HIGHBD_BFP_WRAPPER(32, 64, 8)
        HIGHBD_BFP_WRAPPER(32, 32, 8)
        HIGHBD_BFP_WRAPPER(64, 64, 8)
        HIGHBD_BFP_WRAPPER(16, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 16, 8)
        HIGHBD_BFP_WRAPPER(8, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 4, 8)
        HIGHBD_BFP_WRAPPER(4, 8, 8)
        HIGHBD_BFP_WRAPPER(4, 4, 8)
        HIGHBD_BFP_WRAPPER(128, 128, 8)
        HIGHBD_BFP_WRAPPER(128, 64, 8)
        HIGHBD_BFP_WRAPPER(64, 128, 8)

        HIGHBD_MBFP_WRAPPER(128, 128, 8)
        HIGHBD_MBFP_WRAPPER(128, 64, 8)
        HIGHBD_MBFP_WRAPPER(64, 128, 8)
        HIGHBD_MBFP_WRAPPER(64, 64, 8)
        HIGHBD_MBFP_WRAPPER(64, 32, 8)
        HIGHBD_MBFP_WRAPPER(32, 64, 8)
        HIGHBD_MBFP_WRAPPER(32, 32, 8)
        HIGHBD_MBFP_WRAPPER(32, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 32, 8)
        HIGHBD_MBFP_WRAPPER(16, 16, 8)
        HIGHBD_MBFP_WRAPPER(8, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 8, 8)
        HIGHBD_MBFP_WRAPPER(4, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 4, 8)
        HIGHBD_MBFP_WRAPPER(4, 4, 8)
        HIGHBD_MBFP_WRAPPER(64, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 64, 8)
        HIGHBD_MBFP_WRAPPER(32, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 32, 8)
        HIGHBD_MBFP_WRAPPER(16, 4, 8)
        HIGHBD_MBFP_WRAPPER(4, 16, 8)

        LOWBD_OBFP_WRAPPER(128, 128)
        LOWBD_OBFP_WRAPPER(128, 64)
        LOWBD_OBFP_WRAPPER(64, 128)
        LOWBD_OBFP_WRAPPER(64, 64)
        LOWBD_OBFP_WRAPPER(64, 32)
        LOWBD_OBFP_WRAPPER(32, 64)
        LOWBD_OBFP_WRAPPER(32, 32)
        LOWBD_OBFP_WRAPPER(32, 16)
        LOWBD_OBFP_WRAPPER(16, 32)
        LOWBD_OBFP_WRAPPER(16, 16)
        LOWBD_OBFP_WRAPPER(8, 16)
        LOWBD_OBFP_WRAPPER(16, 8)
        LOWBD_OBFP_WRAPPER(8, 8)
        LOWBD_OBFP_WRAPPER(4, 8)
        LOWBD_OBFP_WRAPPER(8, 4)
        LOWBD_OBFP_WRAPPER(4, 4)
        LOWBD_OBFP_WRAPPER(64, 16)
        LOWBD_OBFP_WRAPPER(16, 64)
        LOWBD_OBFP_WRAPPER(32, 8)
        LOWBD_OBFP_WRAPPER(8, 32)
        LOWBD_OBFP_WRAPPER(16, 4)
        LOWBD_OBFP_WRAPPER(4, 16)
        break;

      case AOM_BITS_10:
        HIGHBD_BFP_WRAPPER(64, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 64, 10)
        HIGHBD_BFP_WRAPPER(32, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 32, 10)
        HIGHBD_BFP_WRAPPER(16, 4, 10)
        HIGHBD_BFP_WRAPPER(4, 16, 10)
        HIGHBD_BFP_WRAPPER(32, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 32, 10)
        HIGHBD_BFP_WRAPPER(64, 32, 10)
        HIGHBD_BFP_WRAPPER(32, 64, 10)
        HIGHBD_BFP_WRAPPER(32, 32, 10)
        HIGHBD_BFP_WRAPPER(64, 64, 10)
        HIGHBD_BFP_WRAPPER(16, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 16, 10)
        HIGHBD_BFP_WRAPPER(8, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 4, 10)
        HIGHBD_BFP_WRAPPER(4, 8, 10)
        HIGHBD_BFP_WRAPPER(4, 4, 10)
        HIGHBD_BFP_WRAPPER(128, 128, 10)
        HIGHBD_BFP_WRAPPER(128, 64, 10)
        HIGHBD_BFP_WRAPPER(64, 128, 10)

        HIGHBD_MBFP_WRAPPER(128, 128, 10)
        HIGHBD_MBFP_WRAPPER(128, 64, 10)
        HIGHBD_MBFP_WRAPPER(64, 128, 10)
        HIGHBD_MBFP_WRAPPER(64, 64, 10)
        HIGHBD_MBFP_WRAPPER(64, 32, 10)
        HIGHBD_MBFP_WRAPPER(32, 64, 10)
        HIGHBD_MBFP_WRAPPER(32, 32, 10)
        HIGHBD_MBFP_WRAPPER(32, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 32, 10)
        HIGHBD_MBFP_WRAPPER(16, 16, 10)
        HIGHBD_MBFP_WRAPPER(8, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 8, 10)
        HIGHBD_MBFP_WRAPPER(4, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 4, 10)
        HIGHBD_MBFP_WRAPPER(4, 4, 10)
        HIGHBD_MBFP_WRAPPER(64, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 64, 10)
        HIGHBD_MBFP_WRAPPER(32, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 32, 10)
        HIGHBD_MBFP_WRAPPER(16, 4, 10)
        HIGHBD_MBFP_WRAPPER(4, 16, 10)

        HIGHBD_OBFP_WRAPPER(128, 128, 10)
        HIGHBD_OBFP_WRAPPER(128, 64, 10)
        HIGHBD_OBFP_WRAPPER(64, 128, 10)
        HIGHBD_OBFP_WRAPPER(64, 64, 10)
        HIGHBD_OBFP_WRAPPER(64, 32, 10)
        HIGHBD_OBFP_WRAPPER(32, 64, 10)
        HIGHBD_OBFP_WRAPPER(32, 32, 10)
        HIGHBD_OBFP_WRAPPER(32, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 32, 10)
        HIGHBD_OBFP_WRAPPER(16, 16, 10)
        HIGHBD_OBFP_WRAPPER(8, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 8, 10)
        HIGHBD_OBFP_WRAPPER(4, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 4, 10)
        HIGHBD_OBFP_WRAPPER(4, 4, 10)
        HIGHBD_OBFP_WRAPPER(64, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 64, 10)
        HIGHBD_OBFP_WRAPPER(32, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 32, 10)
        HIGHBD_OBFP_WRAPPER(16, 4, 10)
        HIGHBD_OBFP_WRAPPER(4, 16, 10)
        break;

      case AOM_BITS_12:
        HIGHBD_BFP_WRAPPER(64, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 64, 12)
        HIGHBD_BFP_WRAPPER(32, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 32, 12)
        HIGHBD_BFP_WRAPPER(16, 4, 12)
        HIGHBD_BFP_WRAPPER(4, 16, 12)
        HIGHBD_BFP_WRAPPER(32, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 32, 12)
        HIGHBD_BFP_WRAPPER(64, 32, 12)
        HIGHBD_BFP_WRAPPER(32, 64, 12)
        HIGHBD_BFP_WRAPPER(32, 32, 12)
        HIGHBD_BFP_WRAPPER(64, 64, 12)
        HIGHBD_BFP_WRAPPER(16, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 16, 12)
        HIGHBD_BFP_WRAPPER(8, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 4, 12)
        HIGHBD_BFP_WRAPPER(4, 8, 12)
        HIGHBD_BFP_WRAPPER(4, 4, 12)
        HIGHBD_BFP_WRAPPER(128, 128, 12)
        HIGHBD_BFP_WRAPPER(128, 64, 12)
        HIGHBD_BFP_WRAPPER(64, 128, 12)

        HIGHBD_MBFP_WRAPPER(128, 128, 12)
        HIGHBD_MBFP_WRAPPER(128, 64, 12)
        HIGHBD_MBFP_WRAPPER(64, 128, 12)
        HIGHBD_MBFP_WRAPPER(64, 64, 12)
        HIGHBD_MBFP_WRAPPER(64, 32, 12)
        HIGHBD_MBFP_WRAPPER(32, 64, 12)
        HIGHBD_MBFP_WRAPPER(32, 32, 12)
        HIGHBD_MBFP_WRAPPER(32, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 32, 12)
        HIGHBD_MBFP_WRAPPER(16, 16, 12)
        HIGHBD_MBFP_WRAPPER(8, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 8, 12)
        HIGHBD_MBFP_WRAPPER(4, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 4, 12)
        HIGHBD_MBFP_WRAPPER(4, 4, 12)
        HIGHBD_MBFP_WRAPPER(64, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 64, 12)
        HIGHBD_MBFP_WRAPPER(32, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 32, 12)
        HIGHBD_MBFP_WRAPPER(16, 4, 12)
        HIGHBD_MBFP_WRAPPER(4, 16, 12)

        HIGHBD_OBFP_WRAPPER(128, 128, 12)
        HIGHBD_OBFP_WRAPPER(128, 64, 12)
        HIGHBD_OBFP_WRAPPER(64, 128, 12)
        HIGHBD_OBFP_WRAPPER(64, 64, 12)
        HIGHBD_OBFP_WRAPPER(64, 32, 12)
        HIGHBD_OBFP_WRAPPER(32, 64, 12)
        HIGHBD_OBFP_WRAPPER(32, 32, 12)
        HIGHBD_OBFP_WRAPPER(32, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 32, 12)
        HIGHBD_OBFP_WRAPPER(16, 16, 12)
        HIGHBD_OBFP_WRAPPER(8, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 8, 12)
        HIGHBD_OBFP_WRAPPER(4, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 4, 12)
        HIGHBD_OBFP_WRAPPER(4, 4, 12)
        HIGHBD_OBFP_WRAPPER(64, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 64, 12)
        HIGHBD_OBFP_WRAPPER(32, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 32, 12)
        HIGHBD_OBFP_WRAPPER(16, 4, 12)
        HIGHBD_OBFP_WRAPPER(4, 16, 12)
        break;

      default:
        assert(0 &&
               "cm->seq_params.bit_depth should be AOM_BITS_8, "
               "AOM_BITS_10 or AOM_BITS_12");
    }
  }
}
#endif  // CONFIG_AV1_HIGHBITDEPTH

#if !CONFIG_REALTIME_ONLY
void process_tpl_stats_frame(AV1_COMP *cpi) {
  const GF_GROUP *const gf_group = &cpi->gf_group;
  AV1_COMMON *const cm = &cpi->common;

  assert(IMPLIES(gf_group->size > 0, gf_group->index < gf_group->size));

  const int tpl_idx = gf_group->index;
  TplParams *const tpl_data = &cpi->tpl_data;
  TplDepFrame *tpl_frame = &tpl_data->tpl_frame[tpl_idx];
  TplDepStats *tpl_stats = tpl_frame->tpl_stats_ptr;

  if (tpl_frame->is_valid) {
    int tpl_stride = tpl_frame->stride;
    int64_t intra_cost_base = 0;
    int64_t mc_dep_cost_base = 0;
    int64_t mc_saved_base = 0;
    int64_t mc_count_base = 0;
    const int step = 1 << tpl_data->tpl_stats_block_mis_log2;
    const int mi_cols_sr = av1_pixels_to_mi(cm->superres_upscaled_width);

    for (int row = 0; row < cm->mi_params.mi_rows; row += step) {
      for (int col = 0; col < mi_cols_sr; col += step) {
        TplDepStats *this_stats = &tpl_stats[av1_tpl_ptr_pos(
            row, col, tpl_stride, tpl_data->tpl_stats_block_mis_log2)];
        int64_t mc_dep_delta =
            RDCOST(tpl_frame->base_rdmult, this_stats->mc_dep_rate,
                   this_stats->mc_dep_dist);
        intra_cost_base += (this_stats->recrf_dist << RDDIV_BITS);
        mc_dep_cost_base +=
            (this_stats->recrf_dist << RDDIV_BITS) + mc_dep_delta;
        mc_count_base += this_stats->mc_count;
        mc_saved_base += this_stats->mc_saved;
      }
    }

    if (mc_dep_cost_base == 0) {
      tpl_frame->is_valid = 0;
    } else {
      aom_clear_system_state();
      cpi->rd.r0 = (double)intra_cost_base / mc_dep_cost_base;
      if (is_frame_arf_and_tpl_eligible(gf_group)) {
        cpi->rd.arf_r0 = cpi->rd.r0;
        if (cpi->lap_enabled) {
          double min_boost_factor = sqrt(cpi->rc.baseline_gf_interval);
          const int gfu_boost = get_gfu_boost_from_r0_lap(
              min_boost_factor, MAX_GFUBOOST_FACTOR, cpi->rd.arf_r0,
              cpi->rc.num_stats_required_for_gfu_boost);
          // printf("old boost %d new boost %d\n", cpi->rc.gfu_boost,
          //        gfu_boost);
          cpi->rc.gfu_boost = combine_prior_with_tpl_boost(
              min_boost_factor, MAX_BOOST_COMBINE_FACTOR, cpi->rc.gfu_boost,
              gfu_boost, cpi->rc.num_stats_used_for_gfu_boost);
        } else {
          const int gfu_boost = (int)(200.0 / cpi->rd.r0);
          cpi->rc.gfu_boost = combine_prior_with_tpl_boost(
              MIN_BOOST_COMBINE_FACTOR, MAX_BOOST_COMBINE_FACTOR,
              cpi->rc.gfu_boost, gfu_boost, cpi->rc.frames_to_key);
        }
      } else if (frame_is_intra_only(cm)) {
        // TODO(debargha): Turn off q adjustment for kf temporarily to
        // reduce impact on speed of encoding. Need to investigate how
        // to mitigate the issue.
        if (cpi->oxcf.rc_cfg.mode == AOM_Q) {
          const int kf_boost =
              get_kf_boost_from_r0(cpi->rd.r0, cpi->rc.frames_to_key);
          if (cpi->lap_enabled) {
            cpi->rc.kf_boost = combine_prior_with_tpl_boost(
                MIN_BOOST_COMBINE_FACTOR, MAX_BOOST_COMBINE_FACTOR,
                cpi->rc.kf_boost, kf_boost,
                cpi->rc.num_stats_used_for_kf_boost);
          } else {
            cpi->rc.kf_boost = combine_prior_with_tpl_boost(
                MIN_BOOST_COMBINE_FACTOR, MAX_BOOST_COMBINE_FACTOR,
                cpi->rc.kf_boost, kf_boost, cpi->rc.frames_to_key);
          }
        }
      }
      cpi->rd.mc_count_base = (double)mc_count_base /
                              (cm->mi_params.mi_rows * cm->mi_params.mi_cols);
      cpi->rd.mc_saved_base = (double)mc_saved_base /
                              (cm->mi_params.mi_rows * cm->mi_params.mi_cols);
      aom_clear_system_state();
    }
  }
}
#endif  // !CONFIG_REALTIME_ONLY
