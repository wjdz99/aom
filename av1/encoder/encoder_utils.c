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

#include "av1/encoder/aq_complexity.h"
#include "av1/encoder/aq_variance.h"
#include "av1/encoder/bitstream.h"
#include "av1/encoder/encodeframe.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encoder_alloc.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/encoder_utils.h"
#include "av1/encoder/grain_test_vectors.h"
#include "av1/encoder/mv_prec.h"
#include "av1/encoder/pickcdef.h"
#include "av1/encoder/picklpf.h"
#include "av1/encoder/pickrst.h"
#include "av1/encoder/rc_utils.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/superres_scale.h"
#include "av1/encoder/var_based_part.h"

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

static void configure_static_seg_features(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const RATE_CONTROL *const rc = &cpi->rc;
  struct segmentation *const seg = &cm->seg;

  int high_q = (int)(rc->avg_q > 48.0);
  int qi_delta;

  // Disable and clear down for KF
  if (cm->current_frame.frame_type == KEY_FRAME) {
    // Clear down the global segmentation map
    memset(cpi->enc_seg.map, 0, cm->mi_params.mi_rows * cm->mi_params.mi_cols);
    seg->update_map = 0;
    seg->update_data = 0;

    // Disable segmentation
    av1_disable_segmentation(seg);

    // Clear down the segment features.
    av1_clearall_segfeatures(seg);
  } else if (cpi->refresh_frame.alt_ref_frame) {
    // If this is an alt ref frame
    // Clear down the global segmentation map
    memset(cpi->enc_seg.map, 0, cm->mi_params.mi_rows * cm->mi_params.mi_cols);
    seg->update_map = 0;
    seg->update_data = 0;

    // Disable segmentation and individual segment features by default
    av1_disable_segmentation(seg);
    av1_clearall_segfeatures(seg);

    // If segmentation was enabled set those features needed for the
    // arf itself.
    if (seg->enabled) {
      seg->update_map = 1;
      seg->update_data = 1;

      qi_delta = av1_compute_qdelta(rc, rc->avg_q, rc->avg_q * 0.875,
                                    cm->seq_params.bit_depth);
      av1_set_segdata(seg, 1, SEG_LVL_ALT_Q, qi_delta - 2);
      av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_Y_H, -2);
      av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_Y_V, -2);
      av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_U, -2);
      av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_V, -2);

      av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_Y_H);
      av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_Y_V);
      av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_U);
      av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_V);

      av1_enable_segfeature(seg, 1, SEG_LVL_ALT_Q);
    }
  } else if (seg->enabled) {
    // All other frames if segmentation has been enabled

    // First normal frame in a valid gf or alt ref group
    if (rc->frames_since_golden == 0) {
      // Set up segment features for normal frames in an arf group
      if (rc->source_alt_ref_active) {
        seg->update_map = 0;
        seg->update_data = 1;

        qi_delta = av1_compute_qdelta(rc, rc->avg_q, rc->avg_q * 1.125,
                                      cm->seq_params.bit_depth);
        av1_set_segdata(seg, 1, SEG_LVL_ALT_Q, qi_delta + 2);
        av1_enable_segfeature(seg, 1, SEG_LVL_ALT_Q);

        av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_Y_H, -2);
        av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_Y_V, -2);
        av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_U, -2);
        av1_set_segdata(seg, 1, SEG_LVL_ALT_LF_V, -2);

        av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_Y_H);
        av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_Y_V);
        av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_U);
        av1_enable_segfeature(seg, 1, SEG_LVL_ALT_LF_V);

        // Segment coding disabled for compred testing
        if (high_q) {
          av1_set_segdata(seg, 1, SEG_LVL_REF_FRAME, ALTREF_FRAME);
          av1_enable_segfeature(seg, 1, SEG_LVL_REF_FRAME);
          av1_enable_segfeature(seg, 1, SEG_LVL_SKIP);
        }
      } else {
        // Disable segmentation and clear down features if alt ref
        // is not active for this group

        av1_disable_segmentation(seg);

        memset(cpi->enc_seg.map, 0,
               cm->mi_params.mi_rows * cm->mi_params.mi_cols);

        seg->update_map = 0;
        seg->update_data = 0;

        av1_clearall_segfeatures(seg);
      }
    } else if (rc->is_src_frame_alt_ref) {
      // Special case where we are coding over the top of a previous
      // alt ref frame.
      // Segment coding disabled for compred testing

      // Enable ref frame features for segment 0 as well
      av1_enable_segfeature(seg, 0, SEG_LVL_REF_FRAME);
      av1_enable_segfeature(seg, 1, SEG_LVL_REF_FRAME);

      // All mbs should use ALTREF_FRAME
      av1_clear_segdata(seg, 0, SEG_LVL_REF_FRAME);
      av1_set_segdata(seg, 0, SEG_LVL_REF_FRAME, ALTREF_FRAME);
      av1_clear_segdata(seg, 1, SEG_LVL_REF_FRAME);
      av1_set_segdata(seg, 1, SEG_LVL_REF_FRAME, ALTREF_FRAME);

      // Skip all MBs if high Q (0,0 mv and skip coeffs)
      if (high_q) {
        av1_enable_segfeature(seg, 0, SEG_LVL_SKIP);
        av1_enable_segfeature(seg, 1, SEG_LVL_SKIP);
      }
      // Enable data update
      seg->update_data = 1;
    } else {
      // All other frames.

      // No updates.. leave things as they are.
      seg->update_map = 0;
      seg->update_data = 0;
    }
  }
}

static void apply_active_map(AV1_COMP *cpi) {
  struct segmentation *const seg = &cpi->common.seg;
  unsigned char *const seg_map = cpi->enc_seg.map;
  const unsigned char *const active_map = cpi->active_map.map;
  int i;

  assert(AM_SEGMENT_ID_ACTIVE == CR_SEGMENT_ID_BASE);

  if (frame_is_intra_only(&cpi->common)) {
    cpi->active_map.enabled = 0;
    cpi->active_map.update = 1;
  }

  if (cpi->active_map.update) {
    if (cpi->active_map.enabled) {
      for (i = 0;
           i < cpi->common.mi_params.mi_rows * cpi->common.mi_params.mi_cols;
           ++i)
        if (seg_map[i] == AM_SEGMENT_ID_ACTIVE) seg_map[i] = active_map[i];
      av1_enable_segmentation(seg);
      av1_enable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_SKIP);
      av1_enable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_H);
      av1_enable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_V);
      av1_enable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_U);
      av1_enable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_V);

      av1_set_segdata(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_H,
                      -MAX_LOOP_FILTER);
      av1_set_segdata(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_V,
                      -MAX_LOOP_FILTER);
      av1_set_segdata(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_U,
                      -MAX_LOOP_FILTER);
      av1_set_segdata(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_V,
                      -MAX_LOOP_FILTER);
    } else {
      av1_disable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_SKIP);
      av1_disable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_H);
      av1_disable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_Y_V);
      av1_disable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_U);
      av1_disable_segfeature(seg, AM_SEGMENT_ID_INACTIVE, SEG_LVL_ALT_LF_V);
      if (seg->enabled) {
        seg->update_data = 1;
        seg->update_map = 1;
      }
    }
    cpi->active_map.update = 0;
  }
}

static void set_size_dependent_vars(AV1_COMP *cpi, int *q, int *bottom_index,
                                    int *top_index) {
  AV1_COMMON *const cm = &cpi->common;

  // Setup variables that depend on the dimensions of the frame.
  av1_set_speed_features_framesize_dependent(cpi, cpi->speed);

#if !CONFIG_REALTIME_ONLY
  if (cpi->oxcf.enable_tpl_model && is_frame_tpl_eligible(cpi)) {
    process_tpl_stats_frame(cpi);
    av1_tpl_rdmult_setup(cpi);
  }
#endif

  // Decide q and q bounds.
  *q = av1_rc_pick_q_and_bounds(cpi, &cpi->rc, cm->width, cm->height,
                                cpi->gf_group.index, bottom_index, top_index);

  // Configure experimental use of segmentation for enhanced coding of
  // static regions if indicated.
  // Only allowed in the second pass of a two pass encode, as it requires
  // lagged coding, and if the relevant speed feature flag is set.
  if (is_stat_consumption_stage_twopass(cpi) &&
      cpi->sf.hl_sf.static_segmentation)
    configure_static_seg_features(cpi);
}

static void reset_film_grain_chroma_params(aom_film_grain_t *pars) {
  pars->num_cr_points = 0;
  pars->cr_mult = 0;
  pars->cr_luma_mult = 0;
  memset(pars->scaling_points_cr, 0, sizeof(pars->scaling_points_cr));
  memset(pars->ar_coeffs_cr, 0, sizeof(pars->ar_coeffs_cr));
  pars->num_cb_points = 0;
  pars->cb_mult = 0;
  pars->cb_luma_mult = 0;
  pars->chroma_scaling_from_luma = 0;
  memset(pars->scaling_points_cb, 0, sizeof(pars->scaling_points_cb));
  memset(pars->ar_coeffs_cb, 0, sizeof(pars->ar_coeffs_cb));
}

void update_film_grain_parameters(struct AV1_COMP *cpi,
                                  const AV1EncoderConfig *oxcf) {
  AV1_COMMON *const cm = &cpi->common;
  cpi->oxcf = *oxcf;

  if (cpi->film_grain_table) {
    aom_film_grain_table_free(cpi->film_grain_table);
    aom_free(cpi->film_grain_table);
    cpi->film_grain_table = NULL;
  }

  if (oxcf->film_grain_test_vector) {
    cm->seq_params.film_grain_params_present = 1;
    if (cm->current_frame.frame_type == KEY_FRAME) {
      memcpy(&cm->film_grain_params,
             film_grain_test_vectors + oxcf->film_grain_test_vector - 1,
             sizeof(cm->film_grain_params));
      if (oxcf->monochrome)
        reset_film_grain_chroma_params(&cm->film_grain_params);
      cm->film_grain_params.bit_depth = cm->seq_params.bit_depth;
      if (cm->seq_params.color_range == AOM_CR_FULL_RANGE) {
        cm->film_grain_params.clip_to_restricted_range = 0;
      }
    }
  } else if (oxcf->film_grain_table_filename) {
    cm->seq_params.film_grain_params_present = 1;

    cpi->film_grain_table = aom_malloc(sizeof(*cpi->film_grain_table));
    memset(cpi->film_grain_table, 0, sizeof(aom_film_grain_table_t));

    aom_film_grain_table_read(cpi->film_grain_table,
                              oxcf->film_grain_table_filename, &cm->error);
  } else {
#if CONFIG_DENOISE
    cm->seq_params.film_grain_params_present = (cpi->oxcf.noise_level > 0);
#else
    cm->seq_params.film_grain_params_present = 0;
#endif
    memset(&cm->film_grain_params, 0, sizeof(cm->film_grain_params));
  }
}

static void set_size_independent_vars(AV1_COMP *cpi) {
  int i;
  AV1_COMMON *const cm = &cpi->common;
  for (i = LAST_FRAME; i <= ALTREF_FRAME; ++i) {
    cm->global_motion[i] = default_warp_params;
  }
  cpi->gm_info.search_done = 0;

  av1_set_speed_features_framesize_independent(cpi, cpi->speed);
  av1_set_rd_speed_thresholds(cpi);
  cm->features.interp_filter = SWITCHABLE;
  cm->features.switchable_motion_mode = 1;
}

static void scale_references(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;
  const int num_planes = av1_num_planes(cm);
  MV_REFERENCE_FRAME ref_frame;

  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    // Need to convert from AOM_REFFRAME to index into ref_mask (subtract 1).
    if (cpi->ref_frame_flags & av1_ref_frame_flag_list[ref_frame]) {
      BufferPool *const pool = cm->buffer_pool;
      const YV12_BUFFER_CONFIG *const ref =
          get_ref_frame_yv12_buf(cm, ref_frame);

      if (ref == NULL) {
        cpi->scaled_ref_buf[ref_frame - 1] = NULL;
        continue;
      }

      if (ref->y_crop_width != cm->width || ref->y_crop_height != cm->height) {
        // Replace the reference buffer with a copy having a thicker border,
        // if the reference buffer is higher resolution than the current
        // frame, and the border is thin.
        if ((ref->y_crop_width > cm->width ||
             ref->y_crop_height > cm->height) &&
            ref->border < AOM_BORDER_IN_PIXELS) {
          RefCntBuffer *ref_fb = get_ref_frame_buf(cm, ref_frame);
          if (aom_yv12_realloc_with_new_border(
                  &ref_fb->buf, AOM_BORDER_IN_PIXELS,
                  cm->features.byte_alignment, num_planes) != 0) {
            aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                               "Failed to allocate frame buffer");
          }
        }
        int force_scaling = 0;
        RefCntBuffer *new_fb = cpi->scaled_ref_buf[ref_frame - 1];
        if (new_fb == NULL) {
          const int new_fb_idx = get_free_fb(cm);
          if (new_fb_idx == INVALID_IDX) {
            aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                               "Unable to find free frame buffer");
          }
          force_scaling = 1;
          new_fb = &pool->frame_bufs[new_fb_idx];
        }

        if (force_scaling || new_fb->buf.y_crop_width != cm->width ||
            new_fb->buf.y_crop_height != cm->height) {
          if (aom_realloc_frame_buffer(
                  &new_fb->buf, cm->width, cm->height,
                  cm->seq_params.subsampling_x, cm->seq_params.subsampling_y,
                  cm->seq_params.use_highbitdepth, AOM_BORDER_IN_PIXELS,
                  cm->features.byte_alignment, NULL, NULL, NULL)) {
            if (force_scaling) {
              // Release the reference acquired in the get_free_fb() call above.
              --new_fb->ref_count;
            }
            aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                               "Failed to allocate frame buffer");
          }
          av1_resize_and_extend_frame(
              ref, &new_fb->buf, (int)cm->seq_params.bit_depth, num_planes);
          cpi->scaled_ref_buf[ref_frame - 1] = new_fb;
          alloc_frame_mvs(cm, new_fb);
        }
      } else {
        RefCntBuffer *buf = get_ref_frame_buf(cm, ref_frame);
        buf->buf.y_crop_width = ref->y_crop_width;
        buf->buf.y_crop_height = ref->y_crop_height;
        cpi->scaled_ref_buf[ref_frame - 1] = buf;
        ++buf->ref_count;
      }
    } else {
      if (!has_no_stats_stage(cpi)) cpi->scaled_ref_buf[ref_frame - 1] = NULL;
    }
  }
}

#if !CONFIG_REALTIME_ONLY
static int get_interp_filter_selected(const AV1_COMMON *const cm,
                                      MV_REFERENCE_FRAME ref,
                                      InterpFilter ifilter) {
  const RefCntBuffer *const buf = get_ref_frame_buf(cm, ref);
  if (buf == NULL) return 0;
  return buf->interp_filter_selected[ifilter];
}

static uint16_t setup_interp_filter_search_mask(AV1_COMP *cpi) {
  const AV1_COMMON *const cm = &cpi->common;
  int ref_total[REF_FRAMES] = { 0 };
  uint16_t mask = ALLOW_ALL_INTERP_FILT_MASK;

  if (cpi->last_frame_type == KEY_FRAME || cpi->refresh_frame.alt_ref_frame)
    return mask;

  for (MV_REFERENCE_FRAME ref = LAST_FRAME; ref <= ALTREF_FRAME; ++ref) {
    for (InterpFilter ifilter = EIGHTTAP_REGULAR; ifilter <= MULTITAP_SHARP;
         ++ifilter) {
      ref_total[ref] += get_interp_filter_selected(cm, ref, ifilter);
    }
  }
  int ref_total_total = (ref_total[LAST2_FRAME] + ref_total[LAST3_FRAME] +
                         ref_total[GOLDEN_FRAME] + ref_total[BWDREF_FRAME] +
                         ref_total[ALTREF2_FRAME] + ref_total[ALTREF_FRAME]);

  for (InterpFilter ifilter = EIGHTTAP_REGULAR; ifilter <= MULTITAP_SHARP;
       ++ifilter) {
    int last_score = get_interp_filter_selected(cm, LAST_FRAME, ifilter) * 30;
    if (ref_total[LAST_FRAME] && last_score <= ref_total[LAST_FRAME]) {
      int filter_score =
          get_interp_filter_selected(cm, LAST2_FRAME, ifilter) * 20 +
          get_interp_filter_selected(cm, LAST3_FRAME, ifilter) * 20 +
          get_interp_filter_selected(cm, GOLDEN_FRAME, ifilter) * 20 +
          get_interp_filter_selected(cm, BWDREF_FRAME, ifilter) * 10 +
          get_interp_filter_selected(cm, ALTREF2_FRAME, ifilter) * 10 +
          get_interp_filter_selected(cm, ALTREF_FRAME, ifilter) * 10;
      if (filter_score < ref_total_total) {
        DUAL_FILTER_TYPE filt_type = ifilter + SWITCHABLE_FILTERS * ifilter;
        reset_interp_filter_allowed_mask(&mask, filt_type);
      }
    }
  }
  return mask;
}

#define STRICT_PSNR_DIFF_THRESH 0.9
// Encode key frame with/without screen content tools to determine whether
// screen content tools should be enabled for this key frame group or not.
// The first encoding is without screen content tools.
// The second encoding is with screen content tools.
// We compare the psnr and frame size to make the decision.
static void screen_content_tools_determination(
    AV1_COMP *cpi, const int allow_screen_content_tools_orig_decision,
    const int allow_intrabc_orig_decision,
    const int is_screen_content_type_orig_decision, const int pass,
    int *projected_size_pass, PSNR_STATS *psnr) {
  AV1_COMMON *const cm = &cpi->common;
  FeatureFlags *const features = &cm->features;
  projected_size_pass[pass] = cpi->rc.projected_frame_size;
#if CONFIG_AV1_HIGHBITDEPTH
  const uint32_t in_bit_depth = cpi->oxcf.input_bit_depth;
  const uint32_t bit_depth = cpi->td.mb.e_mbd.bd;
  aom_calc_highbd_psnr(cpi->source, &cpi->common.cur_frame->buf, &psnr[pass],
                       bit_depth, in_bit_depth);
#else
  aom_calc_psnr(cpi->source, &cpi->common.cur_frame->buf, &psnr[pass]);
#endif
  if (pass != 1) return;

  const double psnr_diff = psnr[1].psnr[0] - psnr[0].psnr[0];
  const int is_sc_encoding_much_better = psnr_diff > STRICT_PSNR_DIFF_THRESH;
  if (is_sc_encoding_much_better) {
    // Use screen content tools, if we get coding gain.
    features->allow_screen_content_tools = 1;
    features->allow_intrabc = cpi->intrabc_used;
    cpi->is_screen_content_type = 1;
  } else {
    // Use original screen content decision.
    features->allow_screen_content_tools =
        allow_screen_content_tools_orig_decision;
    features->allow_intrabc = allow_intrabc_orig_decision;
    cpi->is_screen_content_type = is_screen_content_type_orig_decision;
  }
}

// Set some encoding parameters to make the encoding process fast.
// A fixed block partition size, and a large q is used.
static void set_encoding_params_for_screen_content(AV1_COMP *cpi,
                                                   const int pass) {
  AV1_COMMON *const cm = &cpi->common;
  if (pass == 0) {
    // In the first pass, encode without screen content tools.
    // Use a high q, and a fixed block size for fast encoding.
    cm->features.allow_screen_content_tools = 0;
    cm->features.allow_intrabc = 0;
    cpi->is_screen_content_type = 0;
    cpi->sf.part_sf.partition_search_type = FIXED_PARTITION;
    cpi->sf.part_sf.always_this_block_size = BLOCK_32X32;
    return;
  }
  assert(pass == 1);
  // In the second pass, encode with screen content tools.
  // Use a high q, and a fixed block size for fast encoding.
  cm->features.allow_screen_content_tools = 1;
  // TODO(chengchen): turn intrabc on could lead to data race issue.
  // cm->allow_intrabc = 1;
  cpi->is_screen_content_type = 1;
  cpi->sf.part_sf.partition_search_type = FIXED_PARTITION;
  cpi->sf.part_sf.always_this_block_size = BLOCK_32X32;
}

BLOCK_SIZE select_sb_size(const AV1_COMP *const cpi) {
  const AV1_COMMON *const cm = &cpi->common;
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;

  if (oxcf->superblock_size == AOM_SUPERBLOCK_SIZE_64X64) return BLOCK_64X64;
  if (oxcf->superblock_size == AOM_SUPERBLOCK_SIZE_128X128)
    return BLOCK_128X128;

  assert(oxcf->superblock_size == AOM_SUPERBLOCK_SIZE_DYNAMIC);

  if (cpi->svc.number_spatial_layers > 1) {
    // Use the configured size (top resolution) for spatial layers.
    return AOMMIN(oxcf->frm_dim_cfg.width, oxcf->frm_dim_cfg.height) > 480
               ? BLOCK_128X128
               : BLOCK_64X64;
  }

  // TODO(any): Possibly could improve this with a heuristic.
  // When superres / resize is on, 'cm->width / height' can change between
  // calls, so we don't apply this heuristic there.
  // Things break if superblock size changes between the first pass and second
  // pass encoding, which is why this heuristic is not configured as a
  // speed-feature.
  if (oxcf->superres_cfg.superres_mode == AOM_SUPERRES_NONE &&
      oxcf->resize_cfg.resize_mode == RESIZE_NONE && oxcf->speed >= 1) {
    return AOMMIN(cm->width, cm->height) > 480 ? BLOCK_128X128 : BLOCK_64X64;
  }

  return BLOCK_128X128;
}

static void setup_frame(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  // Set up entropy context depending on frame type. The decoder mandates
  // the use of the default context, index 0, for keyframes and inter
  // frames where the error_resilient_mode or intra_only flag is set. For
  // other inter-frames the encoder currently uses only two contexts;
  // context 1 for ALTREF frames and context 0 for the others.

  if (frame_is_intra_only(cm) || cm->features.error_resilient_mode ||
      cpi->ext_flags.use_primary_ref_none) {
    av1_setup_past_independence(cm);
  }

  if ((cm->current_frame.frame_type == KEY_FRAME && cm->show_frame) ||
      frame_is_sframe(cm)) {
    if (!cpi->seq_params_locked) {
      set_sb_size(&cm->seq_params, select_sb_size(cpi));
    }
  } else {
    const RefCntBuffer *const primary_ref_buf = get_primary_ref_frame_buf(cm);
    if (primary_ref_buf == NULL) {
      av1_setup_past_independence(cm);
      cm->seg.update_map = 1;
      cm->seg.update_data = 1;
    } else {
      *cm->fc = primary_ref_buf->frame_context;
    }
  }

  av1_zero(cm->cur_frame->interp_filter_selected);
  cm->prev_frame = get_primary_ref_frame_buf(cm);
  cpi->vaq_refresh = 0;
}

// Determines whether to use screen content tools for the key frame group.
// This function modifies "cm->features.allow_screen_content_tools",
// "cm->features.allow_intrabc" and "cpi->is_screen_content_type".
static void determine_sc_tools_with_encoding(AV1_COMP *cpi, const int q_orig) {
  AV1_COMMON *const cm = &cpi->common;
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;
  const QuantizationCfg *const q_cfg = &oxcf->q_cfg;
  // Variables to help determine if we should allow screen content tools.
  int projected_size_pass[3] = { 0 };
  PSNR_STATS psnr[3];
  const int is_key_frame = cm->current_frame.frame_type == KEY_FRAME;
  const int allow_screen_content_tools_orig_decision =
      cm->features.allow_screen_content_tools;
  const int allow_intrabc_orig_decision = cm->features.allow_intrabc;
  const int is_screen_content_type_orig_decision = cpi->is_screen_content_type;
  // Turn off the encoding trial for forward key frame and superres.
  if (cpi->sf.rt_sf.use_nonrd_pick_mode || oxcf->kf_cfg.fwd_kf_enabled ||
      cpi->superres_mode != AOM_SUPERRES_NONE || oxcf->mode == REALTIME ||
      is_screen_content_type_orig_decision || !is_key_frame) {
    return;
  }

  // TODO(chengchen): multiple encoding for the lossless mode is time consuming.
  // Find a better way to determine whether screen content tools should be used
  // for lossless coding.
  // Use a high q and a fixed partition to do quick encoding.
  const int q_for_screen_content_quick_run =
      is_lossless_requested(&oxcf->rc_cfg) ? q_orig : AOMMAX(q_orig, 244);
  const int partition_search_type_orig = cpi->sf.part_sf.partition_search_type;
  const BLOCK_SIZE fixed_partition_block_size_orig =
      cpi->sf.part_sf.always_this_block_size;

  // Setup necessary params for encoding, including frame source, etc.
  aom_clear_system_state();

  cpi->source =
      av1_scale_if_required(cm, cpi->unscaled_source, &cpi->scaled_source);
  if (cpi->unscaled_last_source != NULL) {
    cpi->last_source = av1_scale_if_required(cm, cpi->unscaled_last_source,
                                             &cpi->scaled_last_source);
  }

  setup_frame(cpi);

  if (cm->seg.enabled) {
    if (!cm->seg.update_data && cm->prev_frame) {
      segfeatures_copy(&cm->seg, &cm->prev_frame->seg);
      cm->seg.enabled = cm->prev_frame->seg.enabled;
    } else {
      av1_calculate_segdata(&cm->seg);
    }
  } else {
    memset(&cm->seg, 0, sizeof(cm->seg));
  }
  segfeatures_copy(&cm->cur_frame->seg, &cm->seg);
  cm->cur_frame->seg.enabled = cm->seg.enabled;

  // The two encoding passes aim to help determine whether to use screen
  // content tools, with a high q and fixed partition.
  for (int pass = 0; pass < 2; ++pass) {
    set_encoding_params_for_screen_content(cpi, pass);
#if CONFIG_TUNE_VMAF
    if (oxcf->tuning == AOM_TUNE_VMAF_WITH_PREPROCESSING ||
        oxcf->tuning == AOM_TUNE_VMAF_WITHOUT_PREPROCESSING ||
        oxcf->tuning == AOM_TUNE_VMAF_MAX_GAIN) {
      av1_set_quantizer(
          cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel,
          av1_get_vmaf_base_qindex(cpi, q_for_screen_content_quick_run),
          q_cfg->enable_chroma_deltaq);
    } else {
#endif
      av1_set_quantizer(cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel,
                        q_for_screen_content_quick_run,
                        q_cfg->enable_chroma_deltaq);
#if CONFIG_TUNE_VMAF
    }
#endif
    av1_set_speed_features_qindex_dependent(cpi, oxcf->speed);
    if (q_cfg->deltaq_mode != NO_DELTA_Q)
      av1_init_quantizer(&cpi->enc_quant_dequant_params, &cm->quant_params,
                         cm->seq_params.bit_depth);

    av1_set_variance_partition_thresholds(cpi, q_for_screen_content_quick_run,
                                          0);
    // transform / motion compensation build reconstruction frame
    av1_encode_frame(cpi);
    // Screen content decision
    screen_content_tools_determination(
        cpi, allow_screen_content_tools_orig_decision,
        allow_intrabc_orig_decision, is_screen_content_type_orig_decision, pass,
        projected_size_pass, psnr);
  }

  // Set partition speed feature back.
  cpi->sf.part_sf.partition_search_type = partition_search_type_orig;
  cpi->sf.part_sf.always_this_block_size = fixed_partition_block_size_orig;
}
#endif  // CONFIG_REALTIME_ONLY

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

#if !CONFIG_REALTIME_ONLY
#define GM_RECODE_LOOP_NUM4X4_FACTOR 192
static int recode_loop_test_global_motion(
    WarpedMotionParams *const global_motion,
    const int *const global_motion_used, int *const gm_params_cost) {
  int i;
  int recode = 0;
  for (i = LAST_FRAME; i <= ALTREF_FRAME; ++i) {
    if (global_motion[i].wmtype != IDENTITY &&
        global_motion_used[i] * GM_RECODE_LOOP_NUM4X4_FACTOR <
            gm_params_cost[i]) {
      global_motion[i] = default_warp_params;
      assert(global_motion[i].wmtype == IDENTITY);
      gm_params_cost[i] = 0;
      recode = 1;
      // TODO(sarahparker): The earlier condition for recoding here was:
      // "recode |= (rdc->global_motion_used[i] > 0);". Can we bring something
      // similar to that back to speed up global motion?
    }
  }
  return recode;
}
#endif  // !CONFIG_REALTIME_ONLY

static void release_scaled_references(AV1_COMP *cpi) {
  // TODO(isbs): only refresh the necessary frames, rather than all of them
  for (int i = 0; i < INTER_REFS_PER_FRAME; ++i) {
    RefCntBuffer *const buf = cpi->scaled_ref_buf[i];
    if (buf != NULL) {
      --buf->ref_count;
      cpi->scaled_ref_buf[i] = NULL;
    }
  }
}

static void fix_interp_filter(InterpFilter *const interp_filter,
                              const FRAME_COUNTS *const counts) {
  if (*interp_filter == SWITCHABLE) {
    // Check to see if only one of the filters is actually used
    int count[SWITCHABLE_FILTERS] = { 0 };
    int num_filters_used = 0;
    for (int i = 0; i < SWITCHABLE_FILTERS; ++i) {
      for (int j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
        count[i] += counts->switchable_interp[j][i];
      num_filters_used += (count[i] > 0);
    }
    if (num_filters_used == 1) {
      // Only one filter is used. So set the filter at frame level
      for (int i = 0; i < SWITCHABLE_FILTERS; ++i) {
        if (count[i]) {
          if (i == EIGHTTAP_REGULAR) *interp_filter = i;
          break;
        }
      }
    }
  }
}

static void finalize_encoded_frame(AV1_COMP *const cpi) {
  AV1_COMMON *const cm = &cpi->common;
  CurrentFrame *const current_frame = &cm->current_frame;

  if (!cm->seq_params.reduced_still_picture_hdr &&
      encode_show_existing_frame(cm)) {
    RefCntBuffer *const frame_to_show =
        cm->ref_frame_map[cpi->existing_fb_idx_to_show];

    if (frame_to_show == NULL) {
      aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
                         "Buffer does not contain a reconstructed frame");
    }
    assert(frame_to_show->ref_count > 0);
    assign_frame_buffer_p(&cm->cur_frame, frame_to_show);
  }

  if (!encode_show_existing_frame(cm) &&
      cm->seq_params.film_grain_params_present &&
      (cm->show_frame || cm->showable_frame)) {
    // Copy the current frame's film grain params to the its corresponding
    // RefCntBuffer slot.
    cm->cur_frame->film_grain_params = cm->film_grain_params;

    // We must update the parameters if this is not an INTER_FRAME
    if (current_frame->frame_type != INTER_FRAME)
      cm->cur_frame->film_grain_params.update_parameters = 1;

    // Iterate the random seed for the next frame.
    cm->film_grain_params.random_seed += 3381;
    if (cm->film_grain_params.random_seed == 0)
      cm->film_grain_params.random_seed = 7391;
  }

  // Initialise all tiles' contexts from the global frame context
  for (int tile_col = 0; tile_col < cm->tiles.cols; tile_col++) {
    for (int tile_row = 0; tile_row < cm->tiles.rows; tile_row++) {
      const int tile_idx = tile_row * cm->tiles.cols + tile_col;
      cpi->tile_data[tile_idx].tctx = *cm->fc;
    }
  }

  fix_interp_filter(&cm->features.interp_filter, cpi->td.counts);
}

// Refresh reference frame buffers according to refresh_frame_flags.
static void refresh_reference_frames(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  // All buffers are refreshed for shown keyframes and S-frames.

  for (int ref_frame = 0; ref_frame < REF_FRAMES; ref_frame++) {
    if (((cm->current_frame.refresh_frame_flags >> ref_frame) & 1) == 1) {
      assign_frame_buffer_p(&cm->ref_frame_map[ref_frame], cm->cur_frame);
    }
  }
}

/*!\brief Encode a frame without the recode loop, usually used in one-pass
 * encoding and realtime coding.
 *
 * \ingroup high_level_algo
 *
 * \param[in]    cpi             Top-level encoder structure
 *
 * \return Returns a value to indicate if the encoding is done successfully.
 * \retval #AOM_CODEC_OK
 * \retval #AOM_CODEC_ERROR
 */
static int encode_without_recode(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const QuantizationCfg *const q_cfg = &cpi->oxcf.q_cfg;
  int top_index = 0, bottom_index = 0, q = 0;

  set_size_independent_vars(cpi);
  av1_setup_frame_size(cpi);
  set_size_dependent_vars(cpi, &q, &bottom_index, &top_index);

  if (cpi->sf.part_sf.partition_search_type == VAR_BASED_PARTITION)
    variance_partition_alloc(cpi);

  if (cm->current_frame.frame_type == KEY_FRAME) copy_frame_prob_info(cpi);

#if CONFIG_COLLECT_COMPONENT_TIMING
  printf("\n Encoding a frame:");
#endif

  aom_clear_system_state();

  cpi->source =
      av1_scale_if_required(cm, cpi->unscaled_source, &cpi->scaled_source);
  if (cpi->unscaled_last_source != NULL) {
    cpi->last_source = av1_scale_if_required(cm, cpi->unscaled_last_source,
                                             &cpi->scaled_last_source);
  }
  if (!frame_is_intra_only(cm)) scale_references(cpi);

  av1_set_quantizer(cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel, q,
                    q_cfg->enable_chroma_deltaq);
  av1_set_speed_features_qindex_dependent(cpi, cpi->oxcf.speed);
  if (q_cfg->deltaq_mode != NO_DELTA_Q)
    av1_init_quantizer(&cpi->enc_quant_dequant_params, &cm->quant_params,
                       cm->seq_params.bit_depth);
  av1_set_variance_partition_thresholds(cpi, q, 0);
  setup_frame(cpi);

  // Check if this high_source_sad (scene/slide change) frame should be
  // encoded at high/max QP, and if so, set the q and adjust some rate
  // control parameters.
  if (cpi->sf.rt_sf.overshoot_detection_cbr == FAST_DETECTION_MAXQ &&
      cpi->rc.high_source_sad) {
    if (av1_encodedframe_overshoot(cpi, &q)) {
      av1_set_quantizer(cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel, q,
                        q_cfg->enable_chroma_deltaq);
      av1_set_speed_features_qindex_dependent(cpi, cpi->oxcf.speed);
      if (q_cfg->deltaq_mode != NO_DELTA_Q)
        av1_init_quantizer(&cpi->enc_quant_dequant_params, &cm->quant_params,
                           cm->seq_params.bit_depth);
      av1_set_variance_partition_thresholds(cpi, q, 0);
    }
  }

  if (q_cfg->aq_mode == CYCLIC_REFRESH_AQ) {
    suppress_active_map(cpi);
    av1_cyclic_refresh_setup(cpi);
    apply_active_map(cpi);
  }
  if (cm->seg.enabled) {
    if (!cm->seg.update_data && cm->prev_frame) {
      segfeatures_copy(&cm->seg, &cm->prev_frame->seg);
      cm->seg.enabled = cm->prev_frame->seg.enabled;
    } else {
      av1_calculate_segdata(&cm->seg);
    }
  } else {
    memset(&cm->seg, 0, sizeof(cm->seg));
  }
  segfeatures_copy(&cm->cur_frame->seg, &cm->seg);
  cm->cur_frame->seg.enabled = cm->seg.enabled;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_encode_frame_time);
#endif

  // Set the motion vector precision based on mv stats from the last coded
  // frame.
  if (!frame_is_intra_only(cm)) av1_pick_and_set_high_precision_mv(cpi, q);

  // transform / motion compensation build reconstruction frame
  av1_encode_frame(cpi);

  // Update some stats from cyclic refresh.
  if (q_cfg->aq_mode == CYCLIC_REFRESH_AQ && !frame_is_intra_only(cm))
    av1_cyclic_refresh_postencode(cpi);

#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_encode_frame_time);
#endif
#if CONFIG_INTERNAL_STATS
  ++cpi->tot_recode_hits;
#endif

  aom_clear_system_state();

  return AOM_CODEC_OK;
}

static int is_integer_mv(const YV12_BUFFER_CONFIG *cur_picture,
                         const YV12_BUFFER_CONFIG *last_picture,
                         ForceIntegerMVInfo *const force_intpel_info) {
  aom_clear_system_state();
  // check use hash ME
  int k;

  const int block_size = FORCE_INT_MV_DECISION_BLOCK_SIZE;
  const double threshold_current = 0.8;
  const double threshold_average = 0.95;
  const int max_history_size = 32;
  int T = 0;  // total block
  int C = 0;  // match with collocated block
  int S = 0;  // smooth region but not match with collocated block

  const int pic_width = cur_picture->y_width;
  const int pic_height = cur_picture->y_height;
  for (int i = 0; i + block_size <= pic_height; i += block_size) {
    for (int j = 0; j + block_size <= pic_width; j += block_size) {
      const int x_pos = j;
      const int y_pos = i;
      int match = 1;
      T++;

      // check whether collocated block match with current
      uint8_t *p_cur = cur_picture->y_buffer;
      uint8_t *p_ref = last_picture->y_buffer;
      int stride_cur = cur_picture->y_stride;
      int stride_ref = last_picture->y_stride;
      p_cur += (y_pos * stride_cur + x_pos);
      p_ref += (y_pos * stride_ref + x_pos);

      if (cur_picture->flags & YV12_FLAG_HIGHBITDEPTH) {
        uint16_t *p16_cur = CONVERT_TO_SHORTPTR(p_cur);
        uint16_t *p16_ref = CONVERT_TO_SHORTPTR(p_ref);
        for (int tmpY = 0; tmpY < block_size && match; tmpY++) {
          for (int tmpX = 0; tmpX < block_size && match; tmpX++) {
            if (p16_cur[tmpX] != p16_ref[tmpX]) {
              match = 0;
            }
          }
          p16_cur += stride_cur;
          p16_ref += stride_ref;
        }
      } else {
        for (int tmpY = 0; tmpY < block_size && match; tmpY++) {
          for (int tmpX = 0; tmpX < block_size && match; tmpX++) {
            if (p_cur[tmpX] != p_ref[tmpX]) {
              match = 0;
            }
          }
          p_cur += stride_cur;
          p_ref += stride_ref;
        }
      }

      if (match) {
        C++;
        continue;
      }

      if (av1_hash_is_horizontal_perfect(cur_picture, block_size, x_pos,
                                         y_pos) ||
          av1_hash_is_vertical_perfect(cur_picture, block_size, x_pos, y_pos)) {
        S++;
        continue;
      }
    }
  }

  assert(T > 0);
  double cs_rate = ((double)(C + S)) / ((double)(T));

  force_intpel_info->cs_rate_array[force_intpel_info->rate_index] = cs_rate;

  force_intpel_info->rate_index =
      (force_intpel_info->rate_index + 1) % max_history_size;
  force_intpel_info->rate_size++;
  force_intpel_info->rate_size =
      AOMMIN(force_intpel_info->rate_size, max_history_size);

  if (cs_rate < threshold_current) {
    return 0;
  }

  if (C == T) {
    return 1;
  }

  double cs_average = 0.0;

  for (k = 0; k < force_intpel_info->rate_size; k++) {
    cs_average += force_intpel_info->cs_rate_array[k];
  }
  cs_average /= force_intpel_info->rate_size;

  if (cs_average < threshold_average) {
    return 0;
  }

  if ((T - C - S) < 0) {
    return 1;
  }

  if (cs_average > 1.01) {
    return 1;
  }

  return 0;
}

static void set_mb_ssim_rdmult_scaling(AV1_COMP *cpi) {
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  ThreadData *td = &cpi->td;
  MACROBLOCK *x = &td->mb;
  MACROBLOCKD *xd = &x->e_mbd;
  uint8_t *y_buffer = cpi->source->y_buffer;
  const int y_stride = cpi->source->y_stride;
  const int block_size = BLOCK_16X16;

  const int num_mi_w = mi_size_wide[block_size];
  const int num_mi_h = mi_size_high[block_size];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (mi_params->mi_rows + num_mi_h - 1) / num_mi_h;
  double log_sum = 0.0;
  const int use_hbd = cpi->source->flags & YV12_FLAG_HIGHBITDEPTH;

  // Loop through each 16x16 block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double var = 0.0, num_of_var = 0.0;
      const int index = row * num_cols + col;

      // Loop through each 8x8 block.
      for (int mi_row = row * num_mi_h;
           mi_row < mi_params->mi_rows && mi_row < (row + 1) * num_mi_h;
           mi_row += 2) {
        for (int mi_col = col * num_mi_w;
             mi_col < mi_params->mi_cols && mi_col < (col + 1) * num_mi_w;
             mi_col += 2) {
          struct buf_2d buf;
          const int row_offset_y = mi_row << 2;
          const int col_offset_y = mi_col << 2;

          buf.buf = y_buffer + row_offset_y * y_stride + col_offset_y;
          buf.stride = y_stride;

          if (use_hbd) {
            var += av1_high_get_sby_perpixel_variance(cpi, &buf, BLOCK_8X8,
                                                      xd->bd);
          } else {
            var += av1_get_sby_perpixel_variance(cpi, &buf, BLOCK_8X8);
          }

          num_of_var += 1.0;
        }
      }
      var = var / num_of_var;

      // Curve fitting with an exponential model on all 16x16 blocks from the
      // midres dataset.
      var = 67.035434 * (1 - exp(-0.0021489 * var)) + 17.492222;
      cpi->ssim_rdmult_scaling_factors[index] = var;
      log_sum += log(var);
    }
  }
  log_sum = exp(log_sum / (double)(num_rows * num_cols));

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      cpi->ssim_rdmult_scaling_factors[index] /= log_sum;
    }
  }
}

#if CONFIG_SUPERRES_IN_RECODE

static void save_cur_buf(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;
  const YV12_BUFFER_CONFIG *ybf = &cm->cur_frame->buf;
  memset(&cc->copy_buffer, 0, sizeof(cc->copy_buffer));
  if (aom_alloc_frame_buffer(&cc->copy_buffer, ybf->y_crop_width,
                             ybf->y_crop_height, ybf->subsampling_x,
                             ybf->subsampling_y,
                             ybf->flags & YV12_FLAG_HIGHBITDEPTH, ybf->border,
                             cm->features.byte_alignment) != AOM_CODEC_OK) {
    aom_internal_error(
        &cm->error, AOM_CODEC_MEM_ERROR,
        "Failed to allocate copy buffer for saving coding context");
  }
  aom_yv12_copy_frame(ybf, &cc->copy_buffer, av1_num_planes(cm));
}

// Coding context that only needs to be saved when recode loop includes
// filtering (deblocking, CDEF, superres post-encode upscale and/or loop
// restoraton).
static void save_extra_coding_context(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;

  cc->lf = cm->lf;
  cc->cdef_info = cm->cdef_info;
  cc->rc = cpi->rc;
}

static void save_all_coding_context(AV1_COMP *cpi) {
  save_cur_buf(cpi);
  save_extra_coding_context(cpi);
  if (!frame_is_intra_only(&cpi->common)) release_scaled_references(cpi);
}

static void restore_all_coding_context(AV1_COMP *cpi) {
  restore_cur_buf(cpi);
  restore_extra_coding_context(cpi);
  if (!frame_is_intra_only(&cpi->common)) release_scaled_references(cpi);
}

#if !CONFIG_REALTIME_ONLY

/*!\brief Recode loop for encoding one frame. the purpose of encoding one frame
 * for multiple times can be approaching a target bitrate or adjusting the usage
 * of global motions.
 *
 * \ingroup high_level_algo
 *
 * \param[in]    cpi             Top-level encoder structure
 * \param[in]    size            Bitstream size
 * \param[in]    dest            Bitstream output
 *
 * \return Returns a value to indicate if the encoding is done successfully.
 * \retval #AOM_CODEC_OK
 * \retval -1
 * \retval #AOM_CODEC_ERROR
 */
static int encode_with_recode_loop(AV1_COMP *cpi, size_t *size, uint8_t *dest) {
  AV1_COMMON *const cm = &cpi->common;
  RATE_CONTROL *const rc = &cpi->rc;
  GlobalMotionInfo *const gm_info = &cpi->gm_info;
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;
  const QuantizationCfg *const q_cfg = &oxcf->q_cfg;
  const int allow_recode = (cpi->sf.hl_sf.recode_loop != DISALLOW_RECODE);
  // Must allow recode if minimum compression ratio is set.
  assert(IMPLIES(oxcf->rc_cfg.min_cr > 0, allow_recode));

  set_size_independent_vars(cpi);
  if (is_stat_consumption_stage_twopass(cpi) &&
      cpi->sf.interp_sf.adaptive_interp_filter_search)
    cpi->interp_search_flags.interp_filter_search_mask =
        setup_interp_filter_search_mask(cpi);
  cpi->source->buf_8bit_valid = 0;

  av1_setup_frame_size(cpi);

#if CONFIG_SUPERRES_IN_RECODE
  if (av1_superres_in_recode_allowed(cpi) &&
      cpi->superres_mode != AOM_SUPERRES_NONE &&
      cm->superres_scale_denominator == SCALE_NUMERATOR) {
    // Superres mode is currently enabled, but the denominator selected will
    // disable superres. So no need to continue, as we will go through another
    // recode loop for full-resolution after this anyway.
    return -1;
  }
#endif  // CONFIG_SUPERRES_IN_RECODE

  int top_index = 0, bottom_index = 0;
  int q = 0, q_low = 0, q_high = 0;
  set_size_dependent_vars(cpi, &q, &bottom_index, &top_index);
  q_low = bottom_index;
  q_high = top_index;

  if (cpi->sf.part_sf.partition_search_type == VAR_BASED_PARTITION)
    variance_partition_alloc(cpi);

  if (cm->current_frame.frame_type == KEY_FRAME) copy_frame_prob_info(cpi);

#if CONFIG_COLLECT_COMPONENT_TIMING
  printf("\n Encoding a frame:");
#endif

  // Determine whether to use screen content tools using two fast encoding.
  determine_sc_tools_with_encoding(cpi, q);

  // Loop variables
  int loop = 0;
  int loop_count = 0;
  int loop_at_this_size = 0;
  int overshoot_seen = 0;
  int undershoot_seen = 0;
  int low_cr_seen = 0;
  int last_loop_allow_hp = 0;

  do {
    loop = 0;
    aom_clear_system_state();

    // if frame was scaled calculate global_motion_search again if already
    // done
    if (loop_count > 0 && cpi->source && gm_info->search_done) {
      if (cpi->source->y_crop_width != cm->width ||
          cpi->source->y_crop_height != cm->height) {
        gm_info->search_done = 0;
      }
    }
    cpi->source =
        av1_scale_if_required(cm, cpi->unscaled_source, &cpi->scaled_source);
    if (cpi->unscaled_last_source != NULL) {
      cpi->last_source = av1_scale_if_required(cm, cpi->unscaled_last_source,
                                               &cpi->scaled_last_source);
    }

    if (!frame_is_intra_only(cm)) {
      if (loop_count > 0) {
        release_scaled_references(cpi);
      }
      scale_references(cpi);
    }
#if CONFIG_TUNE_VMAF
    if (oxcf->tuning == AOM_TUNE_VMAF_WITH_PREPROCESSING ||
        oxcf->tuning == AOM_TUNE_VMAF_WITHOUT_PREPROCESSING ||
        oxcf->tuning == AOM_TUNE_VMAF_MAX_GAIN) {
      av1_set_quantizer(cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel,
                        av1_get_vmaf_base_qindex(cpi, q),
                        q_cfg->enable_chroma_deltaq);
    } else {
#endif
      av1_set_quantizer(cm, q_cfg->qm_minlevel, q_cfg->qm_maxlevel, q,
                        q_cfg->enable_chroma_deltaq);
#if CONFIG_TUNE_VMAF
    }
#endif
    av1_set_speed_features_qindex_dependent(cpi, oxcf->speed);

    if (q_cfg->deltaq_mode != NO_DELTA_Q)
      av1_init_quantizer(&cpi->enc_quant_dequant_params, &cm->quant_params,
                         cm->seq_params.bit_depth);

    av1_set_variance_partition_thresholds(cpi, q, 0);

    // printf("Frame %d/%d: q = %d, frame_type = %d superres_denom = %d\n",
    //        cm->current_frame.frame_number, cm->show_frame, q,
    //        cm->current_frame.frame_type, cm->superres_scale_denominator);

    if (loop_count == 0) {
      setup_frame(cpi);
    } else if (get_primary_ref_frame_buf(cm) == NULL) {
      // Base q-index may have changed, so we need to assign proper default coef
      // probs before every iteration.
      av1_default_coef_probs(cm);
      av1_setup_frame_contexts(cm);
    }

    if (q_cfg->aq_mode == VARIANCE_AQ) {
      av1_vaq_frame_setup(cpi);
    } else if (q_cfg->aq_mode == COMPLEXITY_AQ) {
      av1_setup_in_frame_q_adj(cpi);
    }

    if (cm->seg.enabled) {
      if (!cm->seg.update_data && cm->prev_frame) {
        segfeatures_copy(&cm->seg, &cm->prev_frame->seg);
        cm->seg.enabled = cm->prev_frame->seg.enabled;
      } else {
        av1_calculate_segdata(&cm->seg);
      }
    } else {
      memset(&cm->seg, 0, sizeof(cm->seg));
    }
    segfeatures_copy(&cm->cur_frame->seg, &cm->seg);
    cm->cur_frame->seg.enabled = cm->seg.enabled;

#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, av1_encode_frame_time);
#endif
    // Set the motion vector precision based on mv stats from the last coded
    // frame.
    if (!frame_is_intra_only(cm)) {
      av1_pick_and_set_high_precision_mv(cpi, q);

      // If the precision has changed during different iteration of the loop,
      // then we need to reset the global motion vectors
      if (loop_count > 0 &&
          cm->features.allow_high_precision_mv != last_loop_allow_hp) {
        gm_info->search_done = 0;
      }
      last_loop_allow_hp = cm->features.allow_high_precision_mv;
    }

    // transform / motion compensation build reconstruction frame
    av1_encode_frame(cpi);

    // Reset the mv_stats in case we are interrupted by an intraframe or an
    // overlay frame.
    if (cpi->mv_stats.valid) {
      av1_zero(cpi->mv_stats);
    }
    // Gather the mv_stats for the next frame
    if (cpi->sf.hl_sf.high_precision_mv_usage == LAST_MV_DATA &&
        av1_frame_allows_smart_mv(cpi)) {
      av1_collect_mv_stats(cpi, q);
    }

#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, av1_encode_frame_time);
#endif

    aom_clear_system_state();

    // Dummy pack of the bitstream using up to date stats to get an
    // accurate estimate of output frame size to determine if we need
    // to recode.
    const int do_dummy_pack =
        (cpi->sf.hl_sf.recode_loop >= ALLOW_RECODE_KFARFGF &&
         oxcf->rc_cfg.mode != AOM_Q) ||
        oxcf->rc_cfg.min_cr > 0;
    if (do_dummy_pack) {
      finalize_encoded_frame(cpi);
      int largest_tile_id = 0;  // Output from bitstream: unused here
      if (av1_pack_bitstream(cpi, dest, size, &largest_tile_id) !=
          AOM_CODEC_OK) {
        return AOM_CODEC_ERROR;
      }

      rc->projected_frame_size = (int)(*size) << 3;
    }

    if (allow_recode) {
      // Update q and decide whether to do a recode loop
      recode_loop_update_q(cpi, &loop, &q, &q_low, &q_high, top_index,
                           bottom_index, &undershoot_seen, &overshoot_seen,
                           &low_cr_seen, loop_at_this_size);
    }

    // Special case for overlay frame.
    if (loop && rc->is_src_frame_alt_ref &&
        rc->projected_frame_size < rc->max_frame_bandwidth) {
      loop = 0;
    }

    if (allow_recode && !cpi->sf.gm_sf.gm_disable_recode &&
        recode_loop_test_global_motion(cm->global_motion,
                                       cpi->td.rd_counts.global_motion_used,
                                       gm_info->params_cost)) {
      loop = 1;
    }

    if (loop) {
      ++loop_count;
      ++loop_at_this_size;

#if CONFIG_INTERNAL_STATS
      ++cpi->tot_recode_hits;
#endif
    }
#if CONFIG_COLLECT_COMPONENT_TIMING
    if (loop) printf("\n Recoding:");
#endif
  } while (loop);

  return AOM_CODEC_OK;
}
#endif  // !CONFIG_REALTIME_ONLY

/*!\brief Select and apply cdef filters and switchable restoration filters
 *
 * \ingroup high_level_algo
 */
static void cdef_restoration_frame(AV1_COMP *cpi, AV1_COMMON *cm,
                                   MACROBLOCKD *xd, int use_restoration,
                                   int use_cdef) {
  MultiThreadInfo *const mt_info = &cpi->mt_info;
  const int num_workers = mt_info->num_workers;
  if (use_restoration)
    av1_loop_restoration_save_boundary_lines(&cm->cur_frame->buf, cm, 0);

  if (use_cdef) {
#if CONFIG_COLLECT_COMPONENT_TIMING
    start_timing(cpi, cdef_time);
#endif
    // Find CDEF parameters
    av1_cdef_search(&cm->cur_frame->buf, cpi->source, cm, xd,
                    cpi->sf.lpf_sf.cdef_pick_method, cpi->td.mb.rdmult);

    // Apply the filter
    av1_cdef_frame(&cm->cur_frame->buf, cm, xd);
#if CONFIG_COLLECT_COMPONENT_TIMING
    end_timing(cpi, cdef_time);
#endif
  } else {
    cm->cdef_info.cdef_bits = 0;
    cm->cdef_info.cdef_strengths[0] = 0;
    cm->cdef_info.nb_cdef_strengths = 1;
    cm->cdef_info.cdef_uv_strengths[0] = 0;
  }

  av1_superres_post_encode(cpi);

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, loop_restoration_time);
#endif
  if (use_restoration) {
    av1_loop_restoration_save_boundary_lines(&cm->cur_frame->buf, cm, 1);
    av1_pick_filter_restoration(cpi->source, cpi);
    if (cm->rst_info[0].frame_restoration_type != RESTORE_NONE ||
        cm->rst_info[1].frame_restoration_type != RESTORE_NONE ||
        cm->rst_info[2].frame_restoration_type != RESTORE_NONE) {
      if (num_workers > 1)
        av1_loop_restoration_filter_frame_mt(
            &cm->cur_frame->buf, cm, 0, mt_info->workers, num_workers,
            &mt_info->lr_row_sync, &cpi->lr_ctxt);
      else
        av1_loop_restoration_filter_frame(&cm->cur_frame->buf, cm, 0,
                                          &cpi->lr_ctxt);
    }
  } else {
    cm->rst_info[0].frame_restoration_type = RESTORE_NONE;
    cm->rst_info[1].frame_restoration_type = RESTORE_NONE;
    cm->rst_info[2].frame_restoration_type = RESTORE_NONE;
  }
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, loop_restoration_time);
#endif
}

/*!\brief Select and apply in-loop deblocking filters, cdef filters, and
 * restoration filters
 *
 * \ingroup high_level_algo
 */
static void loopfilter_frame(AV1_COMP *cpi, AV1_COMMON *cm) {
  MultiThreadInfo *const mt_info = &cpi->mt_info;
  const int num_workers = mt_info->num_workers;
  const int num_planes = av1_num_planes(cm);
  MACROBLOCKD *xd = &cpi->td.mb.e_mbd;

  assert(IMPLIES(is_lossless_requested(&cpi->oxcf.rc_cfg),
                 cm->features.coded_lossless && cm->features.all_lossless));

  const int use_loopfilter =
      !cm->features.coded_lossless && !cm->tiles.large_scale;
  const int use_cdef = cm->seq_params.enable_cdef &&
                       !cm->features.coded_lossless && !cm->tiles.large_scale;
  const int use_restoration = cm->seq_params.enable_restoration &&
                              !cm->features.all_lossless &&
                              !cm->tiles.large_scale;

  struct loopfilter *lf = &cm->lf;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, loop_filter_time);
#endif
  if (use_loopfilter) {
    aom_clear_system_state();
    av1_pick_filter_level(cpi->source, cpi, cpi->sf.lpf_sf.lpf_pick);
  } else {
    lf->filter_level[0] = 0;
    lf->filter_level[1] = 0;
  }

  if (lf->filter_level[0] || lf->filter_level[1]) {
    if (num_workers > 1)
      av1_loop_filter_frame_mt(&cm->cur_frame->buf, cm, xd, 0, num_planes, 0,
#if CONFIG_LPF_MASK
                               0,
#endif
                               mt_info->workers, num_workers,
                               &mt_info->lf_row_sync);
    else
      av1_loop_filter_frame(&cm->cur_frame->buf, cm, xd,
#if CONFIG_LPF_MASK
                            0,
#endif
                            0, num_planes, 0);
  }
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, loop_filter_time);
#endif

  cdef_restoration_frame(cpi, cm, xd, use_restoration, use_cdef);
}

#ifdef OUTPUT_YUV_REC
void aom_write_one_yuv_frame(AV1_COMMON *cm, YV12_BUFFER_CONFIG *s) {
  uint8_t *src = s->y_buffer;
  int h = cm->height;
  if (yuv_rec_file == NULL) return;
  if (s->flags & YV12_FLAG_HIGHBITDEPTH) {
    uint16_t *src16 = CONVERT_TO_SHORTPTR(s->y_buffer);

    do {
      fwrite(src16, s->y_width, 2, yuv_rec_file);
      src16 += s->y_stride;
    } while (--h);

    src16 = CONVERT_TO_SHORTPTR(s->u_buffer);
    h = s->uv_height;

    do {
      fwrite(src16, s->uv_width, 2, yuv_rec_file);
      src16 += s->uv_stride;
    } while (--h);

    src16 = CONVERT_TO_SHORTPTR(s->v_buffer);
    h = s->uv_height;

    do {
      fwrite(src16, s->uv_width, 2, yuv_rec_file);
      src16 += s->uv_stride;
    } while (--h);

    fflush(yuv_rec_file);
    return;
  }

  do {
    fwrite(src, s->y_width, 1, yuv_rec_file);
    src += s->y_stride;
  } while (--h);

  src = s->u_buffer;
  h = s->uv_height;

  do {
    fwrite(src, s->uv_width, 1, yuv_rec_file);
    src += s->uv_stride;
  } while (--h);

  src = s->v_buffer;
  h = s->uv_height;

  do {
    fwrite(src, s->uv_width, 1, yuv_rec_file);
    src += s->uv_stride;
  } while (--h);

  fflush(yuv_rec_file);
}
#endif  // OUTPUT_YUV_REC

/*!\brief Recode loop or a single loop for encoding one frame, followed by
 * in-loop deblocking filters, CDEF filters, and restoration filters.
 *
 * \ingroup high_level_algo
 * \callgraph
 * \callergraph
 *
 * \param[in]    cpi             Top-level encoder structure
 * \param[in]    size            Bitstream size
 * \param[in]    dest            Bitstream output
 * \param[in]    sse             Total distortion of the frame
 * \param[in]    rate            Total rate of the frame
 * \param[in]    largest_tile_id Tile id of the last tile
 *
 * \return Returns a value to indicate if the encoding is done successfully.
 * \retval #AOM_CODEC_OK
 * \retval #AOM_CODEC_ERROR
 */
static int encode_with_recode_loop_and_filter(AV1_COMP *cpi, size_t *size,
                                              uint8_t *dest, int64_t *sse,
                                              int64_t *rate,
                                              int *largest_tile_id) {
#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, encode_with_recode_loop_time);
#endif
  int err;
#if CONFIG_REALTIME_ONLY
  err = encode_without_recode(cpi);
#else
  if (cpi->sf.hl_sf.recode_loop == DISALLOW_RECODE &&
      cpi->oxcf.mode == REALTIME)
    err = encode_without_recode(cpi);
  else
    err = encode_with_recode_loop(cpi, size, dest);
#endif
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, encode_with_recode_loop_time);
#endif
  if (err != AOM_CODEC_OK) {
    if (err == -1) {
      // special case as described in encode_with_recode_loop().
      // Encoding was skipped.
      err = AOM_CODEC_OK;
      if (sse != NULL) *sse = INT64_MAX;
      if (rate != NULL) *rate = INT64_MAX;
      *largest_tile_id = 0;
    }
    return err;
  }

  AV1_COMMON *const cm = &cpi->common;
  SequenceHeader *const seq_params = &cm->seq_params;

  // Special case code to reduce pulsing when key frames are forced at a
  // fixed interval. Note the reconstruction error if it is the frame before
  // the force key frame
  if (cpi->rc.next_key_frame_forced && cpi->rc.frames_to_key == 1) {
#if CONFIG_AV1_HIGHBITDEPTH
    if (seq_params->use_highbitdepth) {
      cpi->ambient_err = aom_highbd_get_y_sse(cpi->source, &cm->cur_frame->buf);
    } else {
      cpi->ambient_err = aom_get_y_sse(cpi->source, &cm->cur_frame->buf);
    }
#else
    cpi->ambient_err = aom_get_y_sse(cpi->source, &cm->cur_frame->buf);
#endif
  }

  cm->cur_frame->buf.color_primaries = seq_params->color_primaries;
  cm->cur_frame->buf.transfer_characteristics =
      seq_params->transfer_characteristics;
  cm->cur_frame->buf.matrix_coefficients = seq_params->matrix_coefficients;
  cm->cur_frame->buf.monochrome = seq_params->monochrome;
  cm->cur_frame->buf.chroma_sample_position =
      seq_params->chroma_sample_position;
  cm->cur_frame->buf.color_range = seq_params->color_range;
  cm->cur_frame->buf.render_width = cm->render_width;
  cm->cur_frame->buf.render_height = cm->render_height;

  // TODO(zoeliu): For non-ref frames, loop filtering may need to be turned
  // off.

  // Pick the loop filter level for the frame.
  if (!cm->features.allow_intrabc) {
    loopfilter_frame(cpi, cm);
  } else {
    cm->lf.filter_level[0] = 0;
    cm->lf.filter_level[1] = 0;
    cm->cdef_info.cdef_bits = 0;
    cm->cdef_info.cdef_strengths[0] = 0;
    cm->cdef_info.nb_cdef_strengths = 1;
    cm->cdef_info.cdef_uv_strengths[0] = 0;
    cm->rst_info[0].frame_restoration_type = RESTORE_NONE;
    cm->rst_info[1].frame_restoration_type = RESTORE_NONE;
    cm->rst_info[2].frame_restoration_type = RESTORE_NONE;
  }

  // TODO(debargha): Fix mv search range on encoder side
  // aom_extend_frame_inner_borders(&cm->cur_frame->buf, av1_num_planes(cm));
  aom_extend_frame_borders(&cm->cur_frame->buf, av1_num_planes(cm));

#ifdef OUTPUT_YUV_REC
  aom_write_one_yuv_frame(cm, &cm->cur_frame->buf);
#endif

  finalize_encoded_frame(cpi);
  // Build the bitstream
#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, av1_pack_bitstream_final_time);
#endif
  if (av1_pack_bitstream(cpi, dest, size, largest_tile_id) != AOM_CODEC_OK)
    return AOM_CODEC_ERROR;
#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, av1_pack_bitstream_final_time);
#endif

  // Compute sse and rate.
  if (sse != NULL) {
#if CONFIG_AV1_HIGHBITDEPTH
    *sse = (seq_params->use_highbitdepth)
               ? aom_highbd_get_y_sse(cpi->source, &cm->cur_frame->buf)
               : aom_get_y_sse(cpi->source, &cm->cur_frame->buf);
#else
    *sse = aom_get_y_sse(cpi->source, &cm->cur_frame->buf);
#endif
  }
  if (rate != NULL) {
    const int64_t bits = (*size << 3);
    *rate = (bits << 5);  // To match scale.
  }
  return AOM_CODEC_OK;
}

static int encode_with_and_without_superres(AV1_COMP *cpi, size_t *size,
                                            uint8_t *dest,
                                            int *largest_tile_id) {
  const AV1_COMMON *const cm = &cpi->common;
  assert(cm->seq_params.enable_superres);
  assert(av1_superres_in_recode_allowed(cpi));
  aom_codec_err_t err = AOM_CODEC_OK;
  save_all_coding_context(cpi);

  // Encode with superres.
#if SUPERRES_RECODE_ALL_RATIOS
  SuperResCfg *const superres_cfg = &cpi->oxcf.superres_cfg;
  int64_t superres_sses[SCALE_NUMERATOR];
  int64_t superres_rates[SCALE_NUMERATOR];
  int superres_largest_tile_ids[SCALE_NUMERATOR];
  // Use superres for Key-frames and Alt-ref frames only.
  const GF_GROUP *const gf_group = &cpi->gf_group;
  if (gf_group->update_type[gf_group->index] != OVERLAY_UPDATE &&
      gf_group->update_type[gf_group->index] != INTNL_OVERLAY_UPDATE) {
    for (int denom = SCALE_NUMERATOR + 1; denom <= 2 * SCALE_NUMERATOR;
         ++denom) {
      superres_cfg->superres_scale_denominator = denom;
      superres_cfg->superres_kf_scale_denominator = denom;
      const int this_index = denom - (SCALE_NUMERATOR + 1);
      err = encode_with_recode_loop_and_filter(
          cpi, size, dest, &superres_sses[this_index],
          &superres_rates[this_index], &superres_largest_tile_ids[this_index]);
      if (err != AOM_CODEC_OK) return err;
      restore_all_coding_context(cpi);
    }
    // Reset.
    superres_cfg->superres_scale_denominator = SCALE_NUMERATOR;
    superres_cfg->superres_kf_scale_denominator = SCALE_NUMERATOR;
  } else {
    for (int denom = SCALE_NUMERATOR + 1; denom <= 2 * SCALE_NUMERATOR;
         ++denom) {
      const int this_index = denom - (SCALE_NUMERATOR + 1);
      superres_sses[this_index] = INT64_MAX;
      superres_rates[this_index] = INT64_MAX;
    }
  }
#else
  int64_t sse1 = INT64_MAX;
  int64_t rate1 = INT64_MAX;
  int largest_tile_id1;
  err = encode_with_recode_loop_and_filter(cpi, size, dest, &sse1, &rate1,
                                           &largest_tile_id1);
  if (err != AOM_CODEC_OK) return err;
  restore_all_coding_context(cpi);
#endif  // SUPERRES_RECODE_ALL_RATIOS

  // Encode without superres.
  int64_t sse2 = INT64_MAX;
  int64_t rate2 = INT64_MAX;
  int largest_tile_id2;
  cpi->superres_mode = AOM_SUPERRES_NONE;  // To force full-res.
  err = encode_with_recode_loop_and_filter(cpi, size, dest, &sse2, &rate2,
                                           &largest_tile_id2);
  cpi->superres_mode = cpi->oxcf.superres_cfg.superres_mode;  // Reset.
  assert(cpi->oxcf.superres_cfg.superres_mode == AOM_SUPERRES_AUTO);
  if (err != AOM_CODEC_OK) return err;

  // Note: Both use common rdmult based on base qindex of fullres.
  const int64_t rdmult =
      av1_compute_rd_mult_based_on_qindex(cpi, cm->quant_params.base_qindex);

#if SUPERRES_RECODE_ALL_RATIOS
  // Find the best rdcost among all superres denoms.
  double proj_rdcost1 = DBL_MAX;
  int64_t sse1 = INT64_MAX;
  int64_t rate1 = INT64_MAX;
  int largest_tile_id1 = 0;
  (void)sse1;
  (void)rate1;
  (void)largest_tile_id1;
  int best_denom = -1;
  for (int denom = SCALE_NUMERATOR + 1; denom <= 2 * SCALE_NUMERATOR; ++denom) {
    const int this_index = denom - (SCALE_NUMERATOR + 1);
    const int64_t this_sse = superres_sses[this_index];
    const int64_t this_rate = superres_rates[this_index];
    const int this_largest_tile_id = superres_largest_tile_ids[this_index];
    const double this_rdcost = RDCOST_DBL(rdmult, this_rate, this_sse);
    if (this_rdcost < proj_rdcost1) {
      sse1 = this_sse;
      rate1 = this_rate;
      largest_tile_id1 = this_largest_tile_id;
      proj_rdcost1 = this_rdcost;
      best_denom = denom;
    }
  }
#else
  const double proj_rdcost1 = RDCOST_DBL(rdmult, rate1, sse1);
#endif  // SUPERRES_RECODE_ALL_RATIOS
  const double proj_rdcost2 = RDCOST_DBL(rdmult, rate2, sse2);

  // Re-encode with superres if it's better.
  if (proj_rdcost1 < proj_rdcost2) {
    restore_all_coding_context(cpi);
    // TODO(urvang): We should avoid rerunning the recode loop by saving
    // previous output+state, or running encode only for the selected 'q' in
    // previous step.
#if SUPERRES_RECODE_ALL_RATIOS
    // Again, temporarily force the best denom.
    superres_cfg->superres_scale_denominator = best_denom;
    superres_cfg->superres_kf_scale_denominator = best_denom;
#endif  // SUPERRES_RECODE_ALL_RATIOS
    int64_t sse3 = INT64_MAX;
    int64_t rate3 = INT64_MAX;
    err = encode_with_recode_loop_and_filter(cpi, size, dest, &sse3, &rate3,
                                             largest_tile_id);
    assert(sse1 == sse3);
    assert(rate1 == rate3);
    assert(largest_tile_id1 == *largest_tile_id);
#if SUPERRES_RECODE_ALL_RATIOS
    // Reset.
    superres_cfg->superres_scale_denominator = SCALE_NUMERATOR;
    superres_cfg->superres_kf_scale_denominator = SCALE_NUMERATOR;
#endif  // SUPERRES_RECODE_ALL_RATIOS
  } else {
    *largest_tile_id = largest_tile_id2;
  }

  release_copy_buffer(&cpi->coding_context);

  return err;
}
#endif  // CONFIG_SUPERRES_IN_RECODE

static void update_reference_segmentation_map(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  MB_MODE_INFO **mi_4x4_ptr = mi_params->mi_grid_base;
  uint8_t *cache_ptr = cm->cur_frame->seg_map;

  for (int row = 0; row < mi_params->mi_rows; row++) {
    MB_MODE_INFO **mi_4x4 = mi_4x4_ptr;
    uint8_t *cache = cache_ptr;
    for (int col = 0; col < mi_params->mi_cols; col++, mi_4x4++, cache++)
      cache[0] = mi_4x4[0]->segment_id;
    mi_4x4_ptr += mi_params->mi_stride;
    cache_ptr += mi_params->mi_cols;
  }
}

#if DUMP_RECON_FRAMES == 1
// NOTE(zoeliu): For debug - Output the filtered reconstructed video.
static void dump_filtered_recon_frames(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const CurrentFrame *const current_frame = &cm->current_frame;
  const YV12_BUFFER_CONFIG *recon_buf = &cm->cur_frame->buf;

  if (recon_buf == NULL) {
    printf("Frame %d is not ready.\n", current_frame->frame_number);
    return;
  }

  static const int flag_list[REF_FRAMES] = { 0,
                                             AOM_LAST_FLAG,
                                             AOM_LAST2_FLAG,
                                             AOM_LAST3_FLAG,
                                             AOM_GOLD_FLAG,
                                             AOM_BWD_FLAG,
                                             AOM_ALT2_FLAG,
                                             AOM_ALT_FLAG };
  printf(
      "\n***Frame=%d (frame_offset=%d, show_frame=%d, "
      "show_existing_frame=%d) "
      "[LAST LAST2 LAST3 GOLDEN BWD ALT2 ALT]=[",
      current_frame->frame_number, current_frame->order_hint, cm->show_frame,
      cm->show_existing_frame);
  for (int ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
    const RefCntBuffer *const buf = get_ref_frame_buf(cm, ref_frame);
    const int ref_offset = buf != NULL ? (int)buf->order_hint : -1;
    printf(" %d(%c)", ref_offset,
           (cpi->ref_frame_flags & flag_list[ref_frame]) ? 'Y' : 'N');
  }
  printf(" ]\n");

  if (!cm->show_frame) {
    printf("Frame %d is a no show frame, so no image dump.\n",
           current_frame->frame_number);
    return;
  }

  int h;
  char file_name[256] = "/tmp/enc_filtered_recon.yuv";
  FILE *f_recon = NULL;

  if (current_frame->frame_number == 0) {
    if ((f_recon = fopen(file_name, "wb")) == NULL) {
      printf("Unable to open file %s to write.\n", file_name);
      return;
    }
  } else {
    if ((f_recon = fopen(file_name, "ab")) == NULL) {
      printf("Unable to open file %s to append.\n", file_name);
      return;
    }
  }
  printf(
      "\nFrame=%5d, encode_update_type[%5d]=%1d, frame_offset=%d, "
      "show_frame=%d, show_existing_frame=%d, source_alt_ref_active=%d, "
      "refresh_alt_ref_frame=%d, "
      "y_stride=%4d, uv_stride=%4d, cm->width=%4d, cm->height=%4d\n\n",
      current_frame->frame_number, cpi->gf_group.index,
      cpi->gf_group.update_type[cpi->gf_group.index], current_frame->order_hint,
      cm->show_frame, cm->show_existing_frame, cpi->rc.source_alt_ref_active,
      cpi->refresh_frame.alt_ref_frame, recon_buf->y_stride,
      recon_buf->uv_stride, cm->width, cm->height);
#if 0
  int ref_frame;
  printf("get_ref_frame_map_idx: [");
  for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame)
    printf(" %d", get_ref_frame_map_idx(cm, ref_frame));
  printf(" ]\n");
#endif  // 0

  // --- Y ---
  for (h = 0; h < cm->height; ++h) {
    fwrite(&recon_buf->y_buffer[h * recon_buf->y_stride], 1, cm->width,
           f_recon);
  }
  // --- U ---
  for (h = 0; h < (cm->height >> 1); ++h) {
    fwrite(&recon_buf->u_buffer[h * recon_buf->uv_stride], 1, (cm->width >> 1),
           f_recon);
  }
  // --- V ---
  for (h = 0; h < (cm->height >> 1); ++h) {
    fwrite(&recon_buf->v_buffer[h * recon_buf->uv_stride], 1, (cm->width >> 1),
           f_recon);
  }

  fclose(f_recon);
}
#endif  // DUMP_RECON_FRAMES

extern void av1_print_frame_contexts(const FRAME_CONTEXT *fc,
                                     const char *filename);

/*!\brief Run the final pass encoding for 1-pass/2-pass encoding mode, and pack
 * the bitstream
 *
 * \ingroup high_level_algo
 * \callgraph
 * \callergraph
 *
 * \param[in]    cpi             Top-level encoder structure
 * \param[in]    size            Bitstream size
 * \param[in]    dest            Bitstream output
 *
 * \return Returns a value to indicate if the encoding is done successfully.
 * \retval #AOM_CODEC_OK
 * \retval #AOM_CODEC_ERROR
 */
int encode_frame_to_data_rate(AV1_COMP *cpi, size_t *size, uint8_t *dest) {
  AV1_COMMON *const cm = &cpi->common;
  SequenceHeader *const seq_params = &cm->seq_params;
  CurrentFrame *const current_frame = &cm->current_frame;
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;
  struct segmentation *const seg = &cm->seg;
  FeatureFlags *const features = &cm->features;
  const TileConfig *const tile_cfg = &oxcf->tile_cfg;

#if CONFIG_COLLECT_COMPONENT_TIMING
  start_timing(cpi, encode_frame_to_data_rate_time);
#endif

  if (frame_is_intra_only(cm)) {
    av1_set_screen_content_options(cpi, features);
    cpi->is_screen_content_type = features->allow_screen_content_tools;
  }

  // frame type has been decided outside of this function call
  cm->cur_frame->frame_type = current_frame->frame_type;

  cm->tiles.large_scale = tile_cfg->enable_large_scale_tile;
  cm->tiles.single_tile_decoding = tile_cfg->enable_single_tile_decoding;

  features->allow_ref_frame_mvs &= frame_might_allow_ref_frame_mvs(cm);
  // features->allow_ref_frame_mvs needs to be written into the frame header
  // while cm->tiles.large_scale is 1, therefore, "cm->tiles.large_scale=1" case
  // is separated from frame_might_allow_ref_frame_mvs().
  features->allow_ref_frame_mvs &= !cm->tiles.large_scale;

  features->allow_warped_motion = oxcf->motion_mode_cfg.allow_warped_motion &&
                                  frame_might_allow_warped_motion(cm);

  cpi->last_frame_type = current_frame->frame_type;

  if (frame_is_sframe(cm)) {
    GF_GROUP *gf_group = &cpi->gf_group;
    RATE_CONTROL *const rc = &cpi->rc;
    // S frame will wipe out any previously encoded altref so we cannot place
    // an overlay frame
    gf_group->update_type[gf_group->size] = GF_UPDATE;
    rc->source_alt_ref_active = 0;
  }

  if (encode_show_existing_frame(cm)) {
    finalize_encoded_frame(cpi);
    // Build the bitstream
    int largest_tile_id = 0;  // Output from bitstream: unused here
    if (av1_pack_bitstream(cpi, dest, size, &largest_tile_id) != AOM_CODEC_OK)
      return AOM_CODEC_ERROR;

    if (seq_params->frame_id_numbers_present_flag &&
        current_frame->frame_type == KEY_FRAME) {
      // Displaying a forward key-frame, so reset the ref buffer IDs
      int display_frame_id = cm->ref_frame_id[cpi->existing_fb_idx_to_show];
      for (int i = 0; i < REF_FRAMES; i++)
        cm->ref_frame_id[i] = display_frame_id;
    }

    cpi->seq_params_locked = 1;

#if DUMP_RECON_FRAMES == 1
    // NOTE(zoeliu): For debug - Output the filtered reconstructed video.
    dump_filtered_recon_frames(cpi);
#endif  // DUMP_RECON_FRAMES

    // NOTE: Save the new show frame buffer index for --test-code=warn, i.e.,
    //       for the purpose to verify no mismatch between encoder and decoder.
    if (cm->show_frame) cpi->last_show_frame_buf = cm->cur_frame;

    refresh_reference_frames(cpi);

    // Since we allocate a spot for the OVERLAY frame in the gf group, we need
    // to do post-encoding update accordingly.
    if (cpi->rc.is_src_frame_alt_ref) {
      av1_set_target_rate(cpi, cm->width, cm->height);
      av1_rc_postencode_update(cpi, *size);
    }

    if (is_psnr_calc_enabled(cpi)) {
      cpi->source =
          realloc_and_scale_source(cpi, cm->cur_frame->buf.y_crop_width,
                                   cm->cur_frame->buf.y_crop_height);
    }
    ++current_frame->frame_number;

    return AOM_CODEC_OK;
  }

  // Work out whether to force_integer_mv this frame
  if (!is_stat_generation_stage(cpi) &&
      cpi->common.features.allow_screen_content_tools &&
      !frame_is_intra_only(cm)) {
    if (cpi->common.seq_params.force_integer_mv == 2) {
      // Adaptive mode: see what previous frame encoded did
      if (cpi->unscaled_last_source != NULL) {
        features->cur_frame_force_integer_mv = is_integer_mv(
            cpi->source, cpi->unscaled_last_source, &cpi->force_intpel_info);
      } else {
        cpi->common.features.cur_frame_force_integer_mv = 0;
      }
    } else {
      cpi->common.features.cur_frame_force_integer_mv =
          cpi->common.seq_params.force_integer_mv;
    }
  } else {
    cpi->common.features.cur_frame_force_integer_mv = 0;
  }

  // Set default state for segment based loop filter update flags.
  cm->lf.mode_ref_delta_update = 0;

  // Set various flags etc to special state if it is a key frame.
  if (frame_is_intra_only(cm) || frame_is_sframe(cm)) {
    // Reset the loop filter deltas and segmentation map.
    av1_reset_segment_features(cm);

    // If segmentation is enabled force a map update for key frames.
    if (seg->enabled) {
      seg->update_map = 1;
      seg->update_data = 1;
    }

    // The alternate reference frame cannot be active for a key frame.
    cpi->rc.source_alt_ref_active = 0;
  }
  if (tile_cfg->mtu == 0) {
    cpi->num_tg = tile_cfg->num_tile_groups;
  } else {
    // Use a default value for the purposes of weighting costs in probability
    // updates
    cpi->num_tg = DEFAULT_MAX_NUM_TG;
  }

  // For 1 pass CBR, check if we are dropping this frame.
  // Never drop on key frame.
  if (has_no_stats_stage(cpi) && oxcf->rc_cfg.mode == AOM_CBR &&
      current_frame->frame_type != KEY_FRAME) {
    if (av1_rc_drop_frame(cpi)) {
      av1_rc_postencode_update_drop_frame(cpi);
      release_scaled_references(cpi);
      return AOM_CODEC_OK;
    }
  }

  if (oxcf->tuning == AOM_TUNE_SSIM) set_mb_ssim_rdmult_scaling(cpi);

#if CONFIG_TUNE_VMAF
  if (oxcf->tuning == AOM_TUNE_VMAF_WITHOUT_PREPROCESSING ||
      oxcf->tuning == AOM_TUNE_VMAF_MAX_GAIN) {
    av1_set_mb_vmaf_rdmult_scaling(cpi);
  }
#endif

  aom_clear_system_state();

#if CONFIG_INTERNAL_STATS
  memset(cpi->mode_chosen_counts, 0,
         MAX_MODES * sizeof(*cpi->mode_chosen_counts));
#endif

  if (seq_params->frame_id_numbers_present_flag) {
    /* Non-normative definition of current_frame_id ("frame counter" with
     * wraparound) */
    if (cm->current_frame_id == -1) {
      int lsb, msb;
      /* quasi-random initialization of current_frame_id for a key frame */
      if (cpi->source->flags & YV12_FLAG_HIGHBITDEPTH) {
        lsb = CONVERT_TO_SHORTPTR(cpi->source->y_buffer)[0] & 0xff;
        msb = CONVERT_TO_SHORTPTR(cpi->source->y_buffer)[1] & 0xff;
      } else {
        lsb = cpi->source->y_buffer[0] & 0xff;
        msb = cpi->source->y_buffer[1] & 0xff;
      }
      cm->current_frame_id =
          ((msb << 8) + lsb) % (1 << seq_params->frame_id_length);

      // S_frame is meant for stitching different streams of different
      // resolutions together, so current_frame_id must be the
      // same across different streams of the same content current_frame_id
      // should be the same and not random. 0x37 is a chosen number as start
      // point
      if (oxcf->kf_cfg.sframe_dist != 0) cm->current_frame_id = 0x37;
    } else {
      cm->current_frame_id =
          (cm->current_frame_id + 1 + (1 << seq_params->frame_id_length)) %
          (1 << seq_params->frame_id_length);
    }
  }

  switch (oxcf->cdf_update_mode) {
    case 0:  // No CDF update for any frames(4~6% compression loss).
      features->disable_cdf_update = 1;
      break;
    case 1:  // Enable CDF update for all frames.
      features->disable_cdf_update = 0;
      break;
    case 2:
      // Strategically determine at which frames to do CDF update.
      // Currently only enable CDF update for all-intra and no-show frames(1.5%
      // compression loss).
      // TODO(huisu@google.com): design schemes for various trade-offs between
      // compression quality and decoding speed.
      features->disable_cdf_update =
          (frame_is_intra_only(cm) || !cm->show_frame) ? 0 : 1;
      break;
  }
  seq_params->timing_info_present &= !seq_params->reduced_still_picture_hdr;

  int largest_tile_id = 0;
#if CONFIG_SUPERRES_IN_RECODE
  if (av1_superres_in_recode_allowed(cpi)) {
    if (encode_with_and_without_superres(cpi, size, dest, &largest_tile_id) !=
        AOM_CODEC_OK) {
      return AOM_CODEC_ERROR;
    }
  } else {
#endif  // CONFIG_SUPERRES_IN_RECODE
    if (encode_with_recode_loop_and_filter(cpi, size, dest, NULL, NULL,
                                           &largest_tile_id) != AOM_CODEC_OK) {
      return AOM_CODEC_ERROR;
    }
#if CONFIG_SUPERRES_IN_RECODE
  }
#endif  // CONFIG_SUPERRES_IN_RECODE

  cpi->seq_params_locked = 1;

  // Update reference frame ids for reference frames this frame will overwrite
  if (seq_params->frame_id_numbers_present_flag) {
    for (int i = 0; i < REF_FRAMES; i++) {
      if ((current_frame->refresh_frame_flags >> i) & 1) {
        cm->ref_frame_id[i] = cm->current_frame_id;
      }
    }
  }

#if DUMP_RECON_FRAMES == 1
  // NOTE(zoeliu): For debug - Output the filtered reconstructed video.
  dump_filtered_recon_frames(cpi);
#endif  // DUMP_RECON_FRAMES

  if (cm->seg.enabled) {
    if (cm->seg.update_map) {
      update_reference_segmentation_map(cpi);
    } else if (cm->last_frame_seg_map) {
      memcpy(cm->cur_frame->seg_map, cm->last_frame_seg_map,
             cm->mi_params.mi_cols * cm->mi_params.mi_rows * sizeof(uint8_t));
    }
  }

  if (frame_is_intra_only(cm) == 0) {
    release_scaled_references(cpi);
  }

  // NOTE: Save the new show frame buffer index for --test-code=warn, i.e.,
  //       for the purpose to verify no mismatch between encoder and decoder.
  if (cm->show_frame) cpi->last_show_frame_buf = cm->cur_frame;

  refresh_reference_frames(cpi);

#if CONFIG_ENTROPY_STATS
  av1_accumulate_frame_counts(&aggregate_fc, &cpi->counts);
#endif  // CONFIG_ENTROPY_STATS

  if (features->refresh_frame_context == REFRESH_FRAME_CONTEXT_BACKWARD) {
    *cm->fc = cpi->tile_data[largest_tile_id].tctx;
    av1_reset_cdf_symbol_counters(cm->fc);
  }
  if (!cm->tiles.large_scale) {
    cm->cur_frame->frame_context = *cm->fc;
  }

  if (tile_cfg->enable_ext_tile_debug) {
    // (yunqing) This test ensures the correctness of large scale tile coding.
    if (cm->tiles.large_scale && is_stat_consumption_stage(cpi)) {
      char fn[20] = "./fc";
      fn[4] = current_frame->frame_number / 100 + '0';
      fn[5] = (current_frame->frame_number % 100) / 10 + '0';
      fn[6] = (current_frame->frame_number % 10) + '0';
      fn[7] = '\0';
      av1_print_frame_contexts(cm->fc, fn);
    }
  }

#if CONFIG_COLLECT_COMPONENT_TIMING
  end_timing(cpi, encode_frame_to_data_rate_time);

  // Print out timing information.
  int i;
  fprintf(stderr, "\n Frame number: %d, Frame type: %s, Show Frame: %d\n",
          cm->current_frame.frame_number,
          get_frame_type_enum(cm->current_frame.frame_type), cm->show_frame);
  for (i = 0; i < kTimingComponents; i++) {
    cpi->component_time[i] += cpi->frame_component_time[i];
    fprintf(stderr, " %s:  %" PRId64 " us (total: %" PRId64 " us)\n",
            get_component_name(i), cpi->frame_component_time[i],
            cpi->component_time[i]);
    cpi->frame_component_time[i] = 0;
  }
#endif

  cpi->last_frame_type = current_frame->frame_type;

  av1_rc_postencode_update(cpi, *size);

  // Clear the one shot update flags for segmentation map and mode/ref loop
  // filter deltas.
  cm->seg.update_map = 0;
  cm->seg.update_data = 0;
  cm->lf.mode_ref_delta_update = 0;

  // A droppable frame might not be shown but it always
  // takes a space in the gf group. Therefore, even when
  // it is not shown, we still need update the count down.

  if (cm->show_frame) {
    // Don't increment frame counters if this was an altref buffer
    // update not a real frame
    ++current_frame->frame_number;
  }

  return AOM_CODEC_OK;
}
