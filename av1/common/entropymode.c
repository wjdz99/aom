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

#include "aom_mem/aom_mem.h"

#include "av1/common/reconinter.h"
#include "av1/common/scan.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/seg_common.h"

#if CONFIG_EC_MULTISYMBOL
aom_cdf_prob av1_kf_y_mode_cdf[INTRA_MODES][INTRA_MODES][INTRA_MODES];
#endif

#if CONFIG_ALT_INTRA

const aom_prob av1_kf_y_mode_prob[INTRA_MODES][INTRA_MODES][INTRA_MODES - 1] = {
  {
      // above = dc
      { 111, 22, 46, 107, 169, 150, 31, 33, 43, 43, 43 },  // left = dc
      { 82, 27, 103, 82, 118, 133, 35, 39, 39, 29, 50 },   // left = v
      { 46, 19, 16, 157, 233, 182, 17, 13, 47, 44, 69 },   // left = h
      { 99, 14, 41, 82, 139, 157, 55, 40, 46, 52, 46 },    // left = d45
      { 60, 9, 27, 129, 48, 161, 32, 26, 80, 40, 33 },     // left = d135
      { 84, 15, 53, 153, 40, 103, 63, 94, 101, 99, 158 },  // left = d117
      { 55, 8, 20, 114, 109, 197, 20, 16, 97, 67, 40 },    // left = d153
      { 67, 6, 25, 84, 176, 168, 25, 19, 39, 95, 54 },     // left = d207
      { 84, 15, 36, 69, 108, 141, 42, 67, 44, 39, 47 },    // left = d63
      { 53, 36, 48, 110, 198, 165, 26, 31, 43, 39, 73 },   // left = paeth
      { 65, 13, 37, 70, 166, 223, 19, 29, 36, 33, 27 },    // left = smooth
      { 50, 53, 49, 114, 204, 163, 28, 28, 46, 45, 76 },   // left = tm
  },
  {
      // above = v
      { 35, 17, 133, 87, 109, 88, 15, 44, 25, 16, 61 },    // left = dc
      { 27, 21, 175, 59, 91, 83, 18, 47, 17, 14, 65 },     // left = v
      { 26, 12, 67, 94, 208, 127, 15, 26, 21, 18, 86 },    // left = h
      { 58, 11, 91, 68, 111, 96, 42, 53, 26, 40, 44 },     // left = d45
      { 27, 9, 94, 101, 41, 119, 24, 45, 40, 17, 58 },     // left = d135
      { 36, 11, 127, 139, 33, 71, 37, 118, 55, 60, 171 },  // left = d117
      { 29, 8, 85, 83, 74, 146, 19, 34, 62, 29, 67 },      // left = d153
      { 40, 7, 83, 66, 143, 118, 25, 37, 27, 47, 65 },     // left = d207
      { 45, 12, 106, 64, 71, 93, 26, 81, 21, 22, 61 },     // left = d63
      { 24, 20, 137, 68, 169, 101, 20, 36, 16, 16, 85 },   // left = paeth
      { 28, 12, 120, 49, 147, 189, 15, 37, 20, 15, 44 },   // left = smooth
      { 22, 21, 136, 60, 153, 99, 22, 44, 16, 19, 85 },    // left = tm
  },
  {
      // above = h
      { 88, 21, 22, 140, 219, 174, 20, 22, 48, 56, 53 },    // left = dc
      { 66, 36, 66, 117, 189, 156, 24, 30, 36, 40, 73 },    // left = v
      { 39, 16, 7, 176, 244, 209, 14, 11, 46, 58, 88 },     // left = h
      { 87, 22, 23, 87, 196, 189, 44, 24, 50, 69, 67 },     // left = d45
      { 67, 17, 20, 122, 114, 176, 26, 27, 70, 52, 52 },    // left = d135
      { 84, 22, 43, 145, 111, 121, 59, 56, 90, 113, 157 },  // left = d117
      { 55, 11, 13, 115, 175, 207, 18, 11, 92, 78, 47 },    // left = d153
      { 58, 8, 10, 109, 224, 196, 21, 13, 39, 110, 66 },    // left = d207
      { 67, 24, 18, 93, 167, 175, 27, 48, 42, 47, 60 },     // left = d63
      { 38, 33, 22, 142, 232, 182, 21, 19, 36, 40, 86 },    // left = paeth
      { 53, 19, 21, 103, 217, 236, 16, 19, 40, 40, 36 },    // left = smooth
      { 36, 40, 24, 141, 230, 191, 24, 22, 39, 41, 77 },    // left = tm
  },
  {
      // above = d45
      { 82, 9, 34, 80, 140, 146, 53, 56, 42, 54, 46 },   // left = dc
      { 63, 13, 86, 57, 111, 118, 71, 70, 37, 30, 47 },  // left = v
      { 48, 9, 15, 117, 217, 164, 38, 33, 42, 46, 60 },  // left = h
      { 91, 5, 16, 42, 168, 158, 120, 62, 30, 63, 64 },  // left = d45
      { 54, 5, 27, 98, 56, 135, 63, 54, 62, 46, 45 },    // left = d135
      { 63, 8, 41, 124, 47, 96, 83, 121, 89, 95, 169 },  // left = d117
      { 50, 5, 16, 99, 111, 175, 42, 37, 74, 79, 42 },   // left = d153
      { 82, 4, 27, 73, 141, 148, 54, 37, 33, 102, 62 },  // left = d207
      { 75, 7, 20, 48, 144, 130, 86, 116, 36, 46, 46 },  // left = d63
      { 46, 17, 42, 75, 168, 135, 77, 56, 42, 44, 80 },  // left = paeth
      { 56, 8, 34, 56, 150, 220, 49, 52, 38, 39, 35 },   // left = smooth
      { 42, 16, 45, 72, 171, 140, 81, 65, 31, 36, 75 },  // left = tm
  },
  {
      // above = d135
      { 58, 8, 23, 123, 68, 176, 18, 23, 82, 50, 44 },     // left = dc
      { 57, 14, 49, 114, 58, 149, 22, 38, 64, 33, 39 },    // left = v
      { 32, 7, 11, 133, 176, 179, 18, 17, 64, 52, 57 },    // left = h
      { 63, 8, 19, 81, 80, 178, 69, 29, 52, 70, 37 },      // left = d45
      { 38, 4, 15, 157, 22, 186, 23, 28, 90, 38, 33 },     // left = d135
      { 55, 8, 29, 171, 27, 134, 41, 72, 141, 116, 141 },  // left = d117
      { 35, 5, 12, 121, 54, 201, 18, 22, 111, 60, 47 },    // left = d153
      { 56, 6, 14, 95, 114, 173, 27, 23, 50, 87, 52 },     // left = d207
      { 59, 11, 24, 88, 61, 146, 34, 41, 56, 48, 38 },     // left = d63
      { 39, 13, 29, 108, 134, 160, 23, 34, 49, 45, 49 },   // left = paeth
      { 47, 9, 26, 90, 93, 221, 19, 29, 47, 40, 23 },      // left = smooth
      { 41, 15, 30, 105, 128, 164, 27, 36, 51, 41, 61 },   // left = tm
  },
  {
      // above = d117
      { 55, 10, 56, 184, 45, 98, 37, 73, 83, 86, 153 },    // left = dc
      { 53, 13, 84, 175, 32, 83, 40, 94, 56, 70, 147 },    // left = v
      { 35, 10, 31, 162, 125, 126, 35, 40, 55, 75, 156 },  // left = h
      { 75, 10, 48, 120, 58, 139, 67, 85, 60, 126, 165 },  // left = d45
      { 38, 6, 37, 190, 17, 119, 45, 72, 95, 97, 156 },    // left = d135
      { 42, 6, 58, 196, 11, 74, 47, 89, 69, 87, 127 },     // left = d117
      { 36, 5, 35, 153, 52, 139, 29, 60, 109, 108, 164 },  // left = d153
      { 45, 5, 39, 115, 83, 142, 42, 68, 46, 118, 145 },   // left = d207
      { 66, 11, 56, 137, 40, 117, 52, 115, 52, 98, 156 },  // left = d63
      { 39, 17, 62, 147, 87, 94, 43, 75, 53, 69, 169 },    // left = paeth
      { 27, 10, 39, 113, 51, 102, 23, 40, 35, 40, 43 },    // left = smooth
      { 34, 17, 64, 151, 94, 105, 48, 63, 46, 74, 165 },   // left = tm
  },
  {
      // above = d153
      { 75, 10, 20, 109, 154, 178, 17, 17, 94, 80, 41 },   // left = dc
      { 70, 21, 49, 92, 134, 149, 23, 30, 75, 43, 53 },    // left = v
      { 34, 6, 6, 134, 212, 207, 13, 10, 79, 68, 71 },     // left = h
      { 81, 12, 17, 86, 135, 201, 41, 29, 71, 73, 49 },    // left = d45
      { 56, 6, 14, 122, 60, 176, 27, 20, 110, 58, 37 },    // left = d135
      { 84, 8, 23, 137, 51, 129, 41, 54, 148, 125, 144 },  // left = d117
      { 47, 4, 8, 93, 111, 211, 11, 9, 136, 92, 41 },      // left = d153
      { 51, 4, 9, 89, 163, 195, 22, 11, 80, 94, 52 },      // left = d207
      { 87, 14, 21, 89, 95, 181, 28, 41, 66, 48, 48 },     // left = d63
      { 54, 15, 23, 109, 188, 173, 20, 21, 72, 62, 66 },   // left = paeth
      { 54, 8, 18, 78, 151, 232, 18, 24, 67, 55, 33 },     // left = smooth
      { 50, 18, 25, 109, 180, 147, 19, 21, 76, 46, 65 },   // left = tm
  },
  {
      // above = d207
      { 89, 12, 25, 92, 170, 178, 28, 30, 42, 69, 46 },    // left = dc
      { 72, 23, 60, 73, 127, 137, 31, 39, 33, 52, 59 },    // left = v
      { 55, 10, 11, 128, 233, 195, 26, 19, 43, 74, 75 },   // left = h
      { 85, 14, 17, 64, 138, 185, 59, 44, 36, 87, 62 },    // left = d45
      { 63, 11, 19, 85, 88, 174, 49, 36, 68, 62, 42 },     // left = d135
      { 77, 11, 37, 106, 81, 128, 68, 71, 83, 128, 161 },  // left = d117
      { 61, 7, 12, 80, 133, 202, 33, 19, 76, 84, 57 },     // left = d153
      { 83, 3, 13, 51, 188, 167, 33, 16, 19, 155, 67 },    // left = d207
      { 72, 17, 21, 64, 120, 164, 43, 56, 38, 64, 42 },    // left = d63
      { 51, 20, 25, 88, 193, 166, 42, 34, 40, 63, 63 },    // left = paeth
      { 62, 11, 21, 71, 179, 232, 27, 23, 35, 59, 31 },    // left = smooth
      { 50, 19, 28, 79, 200, 174, 52, 34, 37, 73, 69 },    // left = tm
  },
  {
      // above = d63
      { 57, 8, 47, 67, 112, 110, 41, 91, 35, 31, 33 },    // left = dc
      { 48, 10, 77, 48, 78, 99, 51, 123, 27, 23, 39 },    // left = v
      { 38, 7, 26, 103, 212, 134, 35, 57, 30, 39, 57 },   // left = h
      { 77, 6, 33, 42, 136, 163, 78, 95, 33, 57, 38 },    // left = d45
      { 38, 6, 31, 76, 57, 134, 57, 84, 46, 34, 32 },     // left = d135
      { 47, 6, 58, 108, 36, 105, 62, 182, 63, 81, 142 },  // left = d117
      { 41, 6, 24, 67, 87, 165, 40, 72, 59, 59, 38 },     // left = d153
      { 62, 4, 29, 53, 151, 152, 49, 76, 24, 72, 44 },    // left = d207
      { 63, 6, 31, 32, 106, 141, 63, 149, 27, 39, 42 },   // left = d63
      { 39, 9, 53, 63, 163, 113, 54, 99, 24, 32, 71 },    // left = paeth
      { 46, 7, 41, 47, 148, 210, 35, 80, 29, 31, 24 },    // left = smooth
      { 36, 9, 57, 44, 147, 111, 57, 116, 23, 34, 59 },   // left = tm
  },
  {
      // above = paeth
      { 56, 29, 74, 108, 182, 120, 18, 31, 33, 36, 76 },   // left = dc
      { 37, 31, 129, 79, 155, 124, 19, 39, 21, 20, 84 },   // left = v
      { 27, 18, 26, 143, 239, 152, 17, 14, 31, 36, 101 },  // left = h
      { 67, 16, 43, 78, 158, 148, 45, 39, 43, 71, 59 },    // left = d45
      { 42, 14, 51, 94, 109, 149, 28, 33, 43, 43, 68 },    // left = d135
      { 57, 21, 83, 132, 93, 98, 53, 106, 53, 55, 167 },   // left = d117
      { 39, 13, 45, 100, 143, 195, 18, 21, 70, 60, 78 },   // left = d153
      { 50, 9, 31, 89, 213, 168, 23, 19, 29, 95, 89 },     // left = d207
      { 48, 19, 57, 67, 128, 131, 27, 55, 34, 30, 62 },    // left = d63
      { 31, 36, 64, 101, 217, 133, 21, 31, 20, 37, 116 },  // left = paeth
      { 39, 19, 55, 80, 198, 226, 16, 26, 27, 25, 54 },    // left = smooth
      { 29, 39, 67, 93, 207, 133, 24, 30, 25, 32, 104 },   // left = tm
  },
  {
      // above = smooth
      { 65, 10, 38, 80, 178, 216, 16, 26, 35, 31, 25 },   // left = dc
      { 47, 20, 89, 60, 146, 214, 17, 35, 28, 18, 34 },   // left = v
      { 36, 10, 18, 126, 235, 228, 16, 18, 33, 35, 48 },  // left = h
      { 65, 8, 29, 58, 154, 233, 40, 37, 37, 48, 35 },    // left = d45
      { 49, 7, 30, 83, 86, 220, 24, 34, 47, 34, 26 },     // left = d135
      { 33, 11, 41, 101, 58, 112, 29, 52, 42, 38, 41 },   // left = d117
      { 50, 7, 23, 77, 144, 230, 19, 25, 63, 49, 34 },    // left = d153
      { 53, 5, 24, 72, 196, 233, 20, 22, 31, 70, 43 },    // left = d207
      { 57, 12, 38, 54, 129, 223, 25, 54, 33, 32, 27 },   // left = d63
      { 38, 18, 47, 86, 209, 223, 22, 26, 30, 28, 42 },   // left = paeth
      { 50, 10, 37, 62, 167, 224, 16, 27, 29, 28, 22 },   // left = smooth
      { 36, 22, 53, 79, 202, 218, 23, 32, 27, 26, 42 },   // left = tm
  },
  {
      // above = tm
      { 50, 35, 74, 104, 179, 111, 20, 39, 38, 38, 73 },  // left = dc
      { 37, 36, 125, 75, 159, 121, 23, 38, 22, 26, 71 },  // left = v
      { 23, 23, 21, 143, 239, 161, 17, 15, 25, 41, 99 },  // left = h
      { 61, 22, 38, 70, 174, 121, 58, 40, 36, 75, 60 },   // left = d45
      { 40, 16, 53, 105, 102, 154, 28, 37, 49, 44, 54 },  // left = d135
      { 52, 25, 82, 131, 97, 98, 45, 92, 58, 63, 155 },   // left = d117
      { 35, 18, 44, 98, 145, 191, 18, 23, 68, 63, 74 },   // left = d153
      { 44, 9, 22, 87, 219, 175, 23, 20, 23, 107, 87 },   // left = d207
      { 49, 20, 50, 74, 147, 135, 33, 58, 31, 36, 67 },   // left = d63
      { 29, 40, 61, 101, 215, 154, 24, 25, 20, 37, 92 },  // left = paeth
      { 40, 24, 53, 81, 204, 222, 16, 25, 25, 30, 45 },   // left = smooth
      { 30, 42, 60, 96, 202, 160, 29, 31, 24, 43, 85 },   // left = tm
  },
};

static const aom_prob default_if_y_probs[BLOCK_SIZE_GROUPS][INTRA_MODES - 1] = {
  { 57, 10, 33, 107, 90, 140, 64, 55, 91, 80, 95 },   // block_size < 8x8
  { 103, 97, 51, 97, 125, 141, 20, 43, 41, 44, 61 },  // block_size < 16x16
  { 127, 91, 51, 98, 189, 136, 15, 31, 27, 27, 39 },  // block_size < 32x32
  { 172, 151, 59, 110, 240, 128, 11, 9, 9, 12, 32 },  // block_size >= 32x32
};

static const aom_prob default_uv_probs[INTRA_MODES][INTRA_MODES - 1] = {
  { 183, 2, 68, 141, 192, 108, 39, 49, 61, 81, 47 },    // y = dc
  { 24, 2, 203, 112, 113, 54, 16, 80, 42, 17, 45 },     // y = v
  { 35, 1, 7, 221, 248, 139, 15, 16, 82, 103, 85 },     // y = h
  { 139, 2, 46, 54, 159, 99, 117, 120, 45, 86, 65 },    // y = d45
  { 130, 1, 28, 166, 50, 162, 18, 28, 114, 61, 37 },    // y = d135
  { 110, 2, 55, 217, 22, 44, 36, 75, 83, 65, 92 },      // y = d117
  { 132, 1, 20, 131, 155, 186, 11, 13, 168, 112, 55 },  // y = d153
  { 129, 2, 12, 118, 226, 127, 36, 23, 37, 167, 97 },   // y = d207
  { 116, 3, 65, 50, 131, 80, 40, 177, 44, 47, 48 },     // y = d63
  { 127, 8, 84, 136, 213, 71, 15, 44, 45, 49, 85 },     // y = paeth
  { 133, 3, 54, 73, 175, 201, 20, 30, 55, 38, 27 },     // y = smooth
  { 147, 8, 89, 146, 208, 66, 15, 53, 62, 53, 80 },     // y = tm
};

#else

const aom_prob av1_kf_y_mode_prob[INTRA_MODES][INTRA_MODES][INTRA_MODES - 1] = {
  {
      // above = dc
      { 137, 30, 42, 148, 151, 207, 70, 52, 91 },   // left = dc
      { 92, 45, 102, 136, 116, 180, 74, 90, 100 },  // left = v
      { 73, 32, 19, 187, 222, 215, 46, 34, 100 },   // left = h
      { 91, 30, 32, 116, 121, 186, 93, 86, 94 },    // left = d45
      { 72, 35, 36, 149, 68, 206, 68, 63, 105 },    // left = d135
      { 73, 31, 28, 138, 57, 124, 55, 122, 151 },   // left = d117
      { 67, 23, 21, 140, 126, 197, 40, 37, 171 },   // left = d153
      { 86, 27, 28, 128, 154, 212, 45, 43, 53 },    // left = d207
      { 74, 32, 27, 107, 86, 160, 63, 134, 102 },   // left = d63
      { 59, 67, 44, 140, 161, 202, 78, 67, 119 }    // left = tm
  },
  {
      // above = v
      { 63, 36, 126, 146, 123, 158, 60, 90, 96 },   // left = dc
      { 43, 46, 168, 134, 107, 128, 69, 142, 92 },  // left = v
      { 44, 29, 68, 159, 201, 177, 50, 57, 77 },    // left = h
      { 58, 38, 76, 114, 97, 172, 78, 133, 92 },    // left = d45
      { 46, 41, 76, 140, 63, 184, 69, 112, 57 },    // left = d135
      { 38, 32, 85, 140, 46, 112, 54, 151, 133 },   // left = d117
      { 39, 27, 61, 131, 110, 175, 44, 75, 136 },   // left = d153
      { 52, 30, 74, 113, 130, 175, 51, 64, 58 },    // left = d207
      { 47, 35, 80, 100, 74, 143, 64, 163, 74 },    // left = d63
      { 36, 61, 116, 114, 128, 162, 80, 125, 82 }   // left = tm
  },
  {
      // above = h
      { 82, 26, 26, 171, 208, 204, 44, 32, 105 },  // left = dc
      { 55, 44, 68, 166, 179, 192, 57, 57, 108 },  // left = v
      { 42, 26, 11, 199, 241, 228, 23, 15, 85 },   // left = h
      { 68, 42, 19, 131, 160, 199, 55, 52, 83 },   // left = d45
      { 58, 50, 25, 139, 115, 232, 39, 52, 118 },  // left = d135
      { 50, 35, 33, 153, 104, 162, 64, 59, 131 },  // left = d117
      { 44, 24, 16, 150, 177, 202, 33, 19, 156 },  // left = d153
      { 55, 27, 12, 153, 203, 218, 26, 27, 49 },   // left = d207
      { 53, 49, 21, 110, 116, 168, 59, 80, 76 },   // left = d63
      { 38, 72, 19, 168, 203, 212, 50, 50, 107 }   // left = tm
  },
  {
      // above = d45
      { 103, 26, 36, 129, 132, 201, 83, 80, 93 },  // left = dc
      { 59, 38, 83, 112, 103, 162, 98, 136, 90 },  // left = v
      { 62, 30, 23, 158, 200, 207, 59, 57, 50 },   // left = h
      { 67, 30, 29, 84, 86, 191, 102, 91, 59 },    // left = d45
      { 60, 32, 33, 112, 71, 220, 64, 89, 104 },   // left = d135
      { 53, 26, 34, 130, 56, 149, 84, 120, 103 },  // left = d117
      { 53, 21, 23, 133, 109, 210, 56, 77, 172 },  // left = d153
      { 77, 19, 29, 112, 142, 228, 55, 66, 36 },   // left = d207
      { 61, 29, 29, 93, 97, 165, 83, 175, 162 },   // left = d63
      { 47, 47, 43, 114, 137, 181, 100, 99, 95 }   // left = tm
  },
  {
      // above = d135
      { 69, 23, 29, 128, 83, 199, 46, 44, 101 },   // left = dc
      { 53, 40, 55, 139, 69, 183, 61, 80, 110 },   // left = v
      { 40, 29, 19, 161, 180, 207, 43, 24, 91 },   // left = h
      { 60, 34, 19, 105, 61, 198, 53, 64, 89 },    // left = d45
      { 52, 31, 22, 158, 40, 209, 58, 62, 89 },    // left = d135
      { 44, 31, 29, 147, 46, 158, 56, 102, 198 },  // left = d117
      { 35, 19, 12, 135, 87, 209, 41, 45, 167 },   // left = d153
      { 55, 25, 21, 118, 95, 215, 38, 39, 66 },    // left = d207
      { 51, 38, 25, 113, 58, 164, 70, 93, 97 },    // left = d63
      { 47, 54, 34, 146, 108, 203, 72, 103, 151 }  // left = tm
  },
  {
      // above = d117
      { 64, 19, 37, 156, 66, 138, 49, 95, 133 },   // left = dc
      { 46, 27, 80, 150, 55, 124, 55, 121, 135 },  // left = v
      { 36, 23, 27, 165, 149, 166, 54, 64, 118 },  // left = h
      { 53, 21, 36, 131, 63, 163, 60, 109, 81 },   // left = d45
      { 40, 26, 35, 154, 40, 185, 51, 97, 123 },   // left = d135
      { 35, 19, 34, 179, 19, 97, 48, 129, 124 },   // left = d117
      { 36, 20, 26, 136, 62, 164, 33, 77, 154 },   // left = d153
      { 45, 18, 32, 130, 90, 157, 40, 79, 91 },    // left = d207
      { 45, 26, 28, 129, 45, 129, 49, 147, 123 },  // left = d63
      { 38, 44, 51, 136, 74, 162, 57, 97, 121 }    // left = tm
  },
  {
      // above = d153
      { 75, 17, 22, 136, 138, 185, 32, 34, 166 },  // left = dc
      { 56, 39, 58, 133, 117, 173, 48, 53, 187 },  // left = v
      { 35, 21, 12, 161, 212, 207, 20, 23, 145 },  // left = h
      { 56, 29, 19, 117, 109, 181, 55, 68, 112 },  // left = d45
      { 47, 29, 17, 153, 64, 220, 59, 51, 114 },   // left = d135
      { 46, 16, 24, 136, 76, 147, 41, 64, 172 },   // left = d117
      { 34, 17, 11, 108, 152, 187, 13, 15, 209 },  // left = d153
      { 51, 24, 14, 115, 133, 209, 32, 26, 104 },  // left = d207
      { 55, 30, 18, 122, 79, 179, 44, 88, 116 },   // left = d63
      { 37, 49, 25, 129, 168, 164, 41, 54, 148 }   // left = tm
  },
  {
      // above = d207
      { 82, 22, 32, 127, 143, 213, 39, 41, 70 },   // left = dc
      { 62, 44, 61, 123, 105, 189, 48, 57, 64 },   // left = v
      { 47, 25, 17, 175, 222, 220, 24, 30, 86 },   // left = h
      { 68, 36, 17, 106, 102, 206, 59, 74, 74 },   // left = d45
      { 57, 39, 23, 151, 68, 216, 55, 63, 58 },    // left = d135
      { 49, 30, 35, 141, 70, 168, 82, 40, 115 },   // left = d117
      { 51, 25, 15, 136, 129, 202, 38, 35, 139 },  // left = d153
      { 68, 26, 16, 111, 141, 215, 29, 28, 28 },   // left = d207
      { 59, 39, 19, 114, 75, 180, 77, 104, 42 },   // left = d63
      { 40, 61, 26, 126, 152, 206, 61, 59, 93 }    // left = tm
  },
  {
      // above = d63
      { 78, 23, 39, 111, 117, 170, 74, 124, 94 },   // left = dc
      { 48, 34, 86, 101, 92, 146, 78, 179, 134 },   // left = v
      { 47, 22, 24, 138, 187, 178, 68, 69, 59 },    // left = h
      { 56, 25, 33, 105, 112, 187, 95, 177, 129 },  // left = d45
      { 48, 31, 27, 114, 63, 183, 82, 116, 56 },    // left = d135
      { 43, 28, 37, 121, 63, 123, 61, 192, 169 },   // left = d117
      { 42, 17, 24, 109, 97, 177, 56, 76, 122 },    // left = d153
      { 58, 18, 28, 105, 139, 182, 70, 92, 63 },    // left = d207
      { 46, 23, 32, 74, 86, 150, 67, 183, 88 },     // left = d63
      { 36, 38, 48, 92, 122, 165, 88, 137, 91 }     // left = tm
  },
  {
      // above = tm
      { 65, 70, 60, 155, 159, 199, 61, 60, 81 },    // left = dc
      { 44, 78, 115, 132, 119, 173, 71, 112, 93 },  // left = v
      { 39, 38, 21, 184, 227, 206, 42, 32, 64 },    // left = h
      { 58, 47, 36, 124, 137, 193, 80, 82, 78 },    // left = d45
      { 49, 50, 35, 144, 95, 205, 63, 78, 59 },     // left = d135
      { 41, 53, 52, 148, 71, 142, 65, 128, 51 },    // left = d117
      { 40, 36, 28, 143, 143, 202, 40, 55, 137 },   // left = d153
      { 52, 34, 29, 129, 183, 227, 42, 35, 43 },    // left = d207
      { 42, 44, 44, 104, 105, 164, 64, 130, 80 },   // left = d63
      { 43, 81, 53, 140, 169, 204, 68, 84, 72 }     // left = tm
  }
};

// Default probabilities for signaling Intra mode for Y plane -- used only for
// inter frames. ('av1_kf_y_mode_prob' is used for intra-only frames).
// Context used: block size group.
static const aom_prob default_if_y_probs[BLOCK_SIZE_GROUPS][INTRA_MODES - 1] = {
  { 65, 32, 18, 144, 162, 194, 41, 51, 98 },   // block_size < 8x8
  { 132, 68, 18, 165, 217, 196, 45, 40, 78 },  // block_size < 16x16
  { 173, 80, 19, 176, 240, 193, 64, 35, 46 },  // block_size < 32x32
  { 221, 135, 38, 194, 248, 121, 96, 85, 29 }  // block_size >= 32x32
};

// Default probabilities for signaling Intra mode for UV plane -- common for
// both intra and inter frames.
// Context used: Intra mode used by Y plane of the same block.
static const aom_prob default_uv_probs[INTRA_MODES][INTRA_MODES - 1] = {
  { 120, 7, 76, 176, 208, 126, 28, 54, 103 },   // y = dc
  { 48, 12, 154, 155, 139, 90, 34, 117, 119 },  // y = v
  { 67, 6, 25, 204, 243, 158, 13, 21, 96 },     // y = h
  { 97, 5, 44, 131, 176, 139, 48, 68, 97 },     // y = d45
  { 83, 5, 42, 156, 111, 152, 26, 49, 152 },    // y = d135
  { 80, 5, 58, 178, 74, 83, 33, 62, 145 },      // y = d117
  { 86, 5, 32, 154, 192, 168, 14, 22, 163 },    // y = d153
  { 85, 5, 32, 156, 216, 148, 19, 29, 73 },     // y = d207
  { 77, 7, 64, 116, 132, 122, 37, 126, 120 },   // y = d63
  { 101, 21, 107, 181, 192, 103, 19, 67, 125 }  // y = tm
};

#endif  // CONFIG_ALT_INTRA

#if CONFIG_EXT_PARTITION_TYPES
static const aom_prob
    default_partition_probs[PARTITION_CONTEXTS][EXT_PARTITION_TYPES - 1] = {
      // 8x8 -> 4x4
      { 199, 122, 141, 128, 128, 128, 128 },  // a/l both not split
      { 147, 63, 159, 128, 128, 128, 128 },   // a split, l not split
      { 148, 133, 118, 128, 128, 128, 128 },  // l split, a not split
      { 121, 104, 114, 128, 128, 128, 128 },  // a/l both split
      // 16x16 -> 8x8
      { 174, 73, 87, 128, 128, 128, 128 },  // a/l both not split
      { 92, 41, 83, 128, 128, 128, 128 },   // a split, l not split
      { 82, 99, 50, 128, 128, 128, 128 },   // l split, a not split
      { 53, 39, 39, 128, 128, 128, 128 },   // a/l both split
      // 32x32 -> 16x16
      { 177, 58, 59, 128, 128, 128, 128 },  // a/l both not split
      { 68, 26, 63, 128, 128, 128, 128 },   // a split, l not split
      { 52, 79, 25, 128, 128, 128, 128 },   // l split, a not split
      { 17, 14, 12, 128, 128, 128, 128 },   // a/l both split
      // 64x64 -> 32x32
      { 222, 34, 30, 128, 128, 128, 128 },  // a/l both not split
      { 72, 16, 44, 128, 128, 128, 128 },   // a split, l not split
      { 58, 32, 12, 128, 128, 128, 128 },   // l split, a not split
      { 10, 7, 6, 128, 128, 128, 128 },     // a/l both split
#if CONFIG_EXT_PARTITION
      // 128x128 -> 64x64
      { 222, 34, 30, 128, 128, 128, 128 },  // a/l both not split
      { 72, 16, 44, 128, 128, 128, 128 },   // a split, l not split
      { 58, 32, 12, 128, 128, 128, 128 },   // l split, a not split
      { 10, 7, 6, 128, 128, 128, 128 },     // a/l both split
#endif                                      // CONFIG_EXT_PARTITION
#if CONFIG_UNPOISON_PARTITION_CTX
      { 0, 0, 141, 0, 0, 0, 0 },  // 8x8 -> 4x4
      { 0, 0, 87, 0, 0, 0, 0 },   // 16x16 -> 8x8
      { 0, 0, 59, 0, 0, 0, 0 },   // 32x32 -> 16x16
      { 0, 0, 30, 0, 0, 0, 0 },   // 64x64 -> 32x32
#if CONFIG_EXT_PARTITION
      { 0, 0, 30, 0, 0, 0, 0 },   // 128x128 -> 64x64
#endif                            // CONFIG_EXT_PARTITION
      { 0, 122, 0, 0, 0, 0, 0 },  // 8x8 -> 4x4
      { 0, 73, 0, 0, 0, 0, 0 },   // 16x16 -> 8x8
      { 0, 58, 0, 0, 0, 0, 0 },   // 32x32 -> 16x16
      { 0, 34, 0, 0, 0, 0, 0 },   // 64x64 -> 32x32
#if CONFIG_EXT_PARTITION
      { 0, 34, 0, 0, 0, 0, 0 },  // 128x128 -> 64x64
#endif                           // CONFIG_EXT_PARTITION
#endif                           // CONFIG_UNPOISON_PARTITION_CTX
    };
#else
static const aom_prob
    default_partition_probs[PARTITION_CONTEXTS][PARTITION_TYPES - 1] = {
      // 8x8 -> 4x4
      { 199, 122, 141 },  // a/l both not split
      { 147, 63, 159 },   // a split, l not split
      { 148, 133, 118 },  // l split, a not split
      { 121, 104, 114 },  // a/l both split
      // 16x16 -> 8x8
      { 174, 73, 87 },  // a/l both not split
      { 92, 41, 83 },   // a split, l not split
      { 82, 99, 50 },   // l split, a not split
      { 53, 39, 39 },   // a/l both split
      // 32x32 -> 16x16
      { 177, 58, 59 },  // a/l both not split
      { 68, 26, 63 },   // a split, l not split
      { 52, 79, 25 },   // l split, a not split
      { 17, 14, 12 },   // a/l both split
      // 64x64 -> 32x32
      { 222, 34, 30 },  // a/l both not split
      { 72, 16, 44 },   // a split, l not split
      { 58, 32, 12 },   // l split, a not split
      { 10, 7, 6 },     // a/l both split
#if CONFIG_EXT_PARTITION
      // 128x128 -> 64x64
      { 222, 34, 30 },  // a/l both not split
      { 72, 16, 44 },   // a split, l not split
      { 58, 32, 12 },   // l split, a not split
      { 10, 7, 6 },     // a/l both split
#endif  // CONFIG_EXT_PARTITION
#if CONFIG_UNPOISON_PARTITION_CTX
      { 0, 0, 141 },    // 8x8 -> 4x4
      { 0, 0, 87 },     // 16x16 -> 8x8
      { 0, 0, 59 },     // 32x32 -> 16x16
      { 0, 0, 30 },     // 64x64 -> 32x32
#if CONFIG_EXT_PARTITION
      { 0, 0, 30 },     // 128x128 -> 64x64
#endif  // CONFIG_EXT_PARTITION
      { 0, 122, 0 },    // 8x8 -> 4x4
      { 0, 73, 0 },     // 16x16 -> 8x8
      { 0, 58, 0 },     // 32x32 -> 16x16
      { 0, 34, 0 },     // 64x64 -> 32x32
#if CONFIG_EXT_PARTITION
      { 0, 34, 0 },     // 128x128 -> 64x64
#endif  // CONFIG_EXT_PARTITION
#endif  // CONFIG_UNPOISON_PARTITION_CTX
    };
#endif  // CONFIG_EXT_PARTITION_TYPES

#if CONFIG_REF_MV
static const aom_prob default_newmv_prob[NEWMV_MODE_CONTEXTS] = {
  200, 180, 150, 150, 110, 70, 60,
};

static const aom_prob default_zeromv_prob[ZEROMV_MODE_CONTEXTS] = {
  192, 64,
};

static const aom_prob default_refmv_prob[REFMV_MODE_CONTEXTS] = {
  220, 220, 200, 200, 180, 128, 30, 220, 30,
};

static const aom_prob default_drl_prob[DRL_MODE_CONTEXTS] = { 128, 160, 180,
                                                              128, 160 };

#if CONFIG_EXT_INTER
static const aom_prob default_new2mv_prob = 180;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV

static const aom_prob
    default_inter_mode_probs[INTER_MODE_CONTEXTS][INTER_MODES - 1] = {
#if CONFIG_EXT_INTER
      // TODO(zoeliu): To adjust the initial default probs
      { 2, 173, 34, 173 },  // 0 = both zero mv
      { 7, 145, 85, 145 },  // 1 = one zero mv + one a predicted mv
      { 7, 166, 63, 166 },  // 2 = two predicted mvs
      { 7, 94, 66, 128 },   // 3 = one predicted/zero and one new mv
      { 8, 64, 46, 128 },   // 4 = two new mvs
      { 17, 81, 31, 128 },  // 5 = one intra neighbour + x
      { 25, 29, 30, 96 },   // 6 = two intra neighbours
#else
      { 2, 173, 34 },  // 0 = both zero mv
      { 7, 145, 85 },  // 1 = one zero mv + one a predicted mv
      { 7, 166, 63 },  // 2 = two predicted mvs
      { 7, 94, 66 },   // 3 = one predicted/zero and one new mv
      { 8, 64, 46 },   // 4 = two new mvs
      { 17, 81, 31 },  // 5 = one intra neighbour + x
      { 25, 29, 30 },  // 6 = two intra neighbours
#endif  // CONFIG_EXT_INTER
    };

#if CONFIG_EXT_INTER
static const aom_prob default_inter_compound_mode_probs
    [INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES - 1] = {
      { 2, 173, 68, 192, 64, 192, 128, 180, 180 },   // 0 = both zero mv
      { 7, 145, 160, 192, 64, 192, 128, 180, 180 },  // 1 = 1 zero + 1 predicted
      { 7, 166, 126, 192, 64, 192, 128, 180, 180 },  // 2 = two predicted mvs
      { 7, 94, 132, 192, 64, 192, 128, 180, 180 },   // 3 = 1 pred/zero, 1 new
      { 8, 64, 64, 192, 64, 192, 128, 180, 180 },    // 4 = two new mvs
      { 17, 81, 52, 192, 64, 192, 128, 180, 180 },   // 5 = one intra neighbour
      { 25, 29, 50, 192, 64, 192, 128, 180, 180 },   // 6 = two intra neighbours
    };

#if CONFIG_COMPOUND_SEGMENT
static const aom_prob
    default_compound_type_probs[BLOCK_SIZES][COMPOUND_TYPES - 1] = {
#if CONFIG_CB4X4
      { 255, 255 }, { 255, 255 }, { 255, 255 },
#endif
      { 208, 200 }, { 208, 200 }, { 208, 200 }, { 208, 200 }, { 208, 200 },
      { 208, 200 }, { 216, 200 }, { 216, 200 }, { 216, 200 }, { 224, 200 },
      { 224, 200 }, { 240, 200 }, { 240, 200 },
#if CONFIG_EXT_PARTITION
      { 255, 200 }, { 255, 200 }, { 255, 200 },
#endif  // CONFIG_EXT_PARTITION
    };
#else  // !CONFIG_COMPOUND_SEGMENT
static const aom_prob default_compound_type_probs[BLOCK_SIZES]
                                                 [COMPOUND_TYPES - 1] = {
#if CONFIG_CB4X4
                                                   { 208 }, { 208 }, { 208 },
#endif
                                                   { 208 }, { 208 }, { 208 },
                                                   { 208 }, { 208 }, { 208 },
                                                   { 216 }, { 216 }, { 216 },
                                                   { 224 }, { 224 }, { 240 },
                                                   { 240 },
#if CONFIG_EXT_PARTITION
                                                   { 255 }, { 255 }, { 255 },
#endif  // CONFIG_EXT_PARTITION
                                                 };
#endif  // CONFIG_COMPOUND_SEGMENT

static const aom_prob default_interintra_prob[BLOCK_SIZE_GROUPS] = {
  208, 208, 208, 208,
};

static const aom_prob
    default_interintra_mode_prob[BLOCK_SIZE_GROUPS][INTERINTRA_MODES - 1] = {
      { 65, 32, 18, 144, 162, 194, 41, 51, 98 },   // block_size < 8x8
      { 132, 68, 18, 165, 217, 196, 45, 40, 78 },  // block_size < 16x16
      { 173, 80, 19, 176, 240, 193, 64, 35, 46 },  // block_size < 32x32
      { 221, 135, 38, 194, 248, 121, 96, 85, 29 }  // block_size >= 32x32
    };

static const aom_prob default_wedge_interintra_prob[BLOCK_SIZES] = {
#if CONFIG_CB4X4
  208, 208, 208,
#endif
  208, 208, 208, 208, 208, 208, 216, 216, 216, 224, 224, 224, 240,
#if CONFIG_EXT_PARTITION
  208, 208, 208
#endif  // CONFIG_EXT_PARTITION
};
#endif  // CONFIG_EXT_INTER

// Change this section appropriately once warped motion is supported
#if CONFIG_MOTION_VAR && !CONFIG_WARPED_MOTION
const aom_tree_index av1_motion_mode_tree[TREE_SIZE(MOTION_MODES)] = {
  -SIMPLE_TRANSLATION, -OBMC_CAUSAL
};
static const aom_prob default_motion_mode_prob[BLOCK_SIZES]
                                              [MOTION_MODES - 1] = {
#if CONFIG_CB4X4
                                                { 255 }, { 255 }, { 255 },
#endif
                                                { 255 }, { 255 }, { 255 },
                                                { 151 }, { 153 }, { 144 },
                                                { 178 }, { 165 }, { 160 },
                                                { 207 }, { 195 }, { 168 },
                                                { 244 },
#if CONFIG_EXT_PARTITION
                                                { 252 }, { 252 }, { 252 },
#endif  // CONFIG_EXT_PARTITION
                                              };

#elif !CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION

const aom_tree_index av1_motion_mode_tree[TREE_SIZE(MOTION_MODES)] = {
  -SIMPLE_TRANSLATION, -WARPED_CAUSAL
};

static const aom_prob default_motion_mode_prob[BLOCK_SIZES]
                                              [MOTION_MODES - 1] = {
#if CONFIG_CB4X4
                                                { 255 }, { 255 }, { 255 },
#endif
                                                { 255 }, { 255 }, { 255 },
                                                { 151 }, { 153 }, { 144 },
                                                { 178 }, { 165 }, { 160 },
                                                { 207 }, { 195 }, { 168 },
                                                { 244 },
#if CONFIG_EXT_PARTITION
                                                { 252 }, { 252 }, { 252 },
#endif  // CONFIG_EXT_PARTITION
                                              };

#elif CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION

const aom_tree_index av1_motion_mode_tree[TREE_SIZE(MOTION_MODES)] = {
  -SIMPLE_TRANSLATION, 2, -OBMC_CAUSAL, -WARPED_CAUSAL,
};
static const aom_prob default_motion_mode_prob[BLOCK_SIZES][MOTION_MODES - 1] =
    {
#if CONFIG_CB4X4
      { 255, 200 }, { 255, 200 }, { 255, 200 },
#endif
      { 255, 200 }, { 255, 200 }, { 255, 200 }, { 151, 200 }, { 153, 200 },
      { 144, 200 }, { 178, 200 }, { 165, 200 }, { 160, 200 }, { 207, 200 },
      { 195, 200 }, { 168, 200 }, { 244, 200 },
#if CONFIG_EXT_PARTITION
      { 252, 200 }, { 252, 200 }, { 252, 200 },
#endif  // CONFIG_EXT_PARTITION
    };

// Probability for the case that only 1 additional motion mode is allowed
static const aom_prob default_obmc_prob[BLOCK_SIZES] = {
#if CONFIG_CB4X4
  255, 255, 255,
#endif
  255, 255, 255, 151, 153, 144, 178, 165, 160, 207, 195, 168, 244,
#if CONFIG_EXT_PARTITION
  252, 252, 252,
#endif  // CONFIG_EXT_PARTITION
};
#endif

#if CONFIG_DELTA_Q
static const aom_prob default_delta_q_probs[DELTA_Q_CONTEXTS] = { 220, 220,
                                                                  220 };
#endif
#if CONFIG_EC_MULTISYMBOL
int av1_intra_mode_ind[INTRA_MODES];
int av1_intra_mode_inv[INTRA_MODES];
int av1_inter_mode_ind[INTER_MODES];
int av1_inter_mode_inv[INTER_MODES];
#endif

/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */
#if CONFIG_ALT_INTRA
const aom_tree_index av1_intra_mode_tree[TREE_SIZE(INTRA_MODES)] = {
  -DC_PRED,    2,            /* 0 = DC_NODE */
  -TM_PRED,    4,            /* 1 = TM_NODE */
  -V_PRED,     6,            /* 2 = V_NODE */
  8,           12,           /* 3 = COM_NODE */
  -H_PRED,     10,           /* 4 = H_NODE */
  -D135_PRED,  -D117_PRED,   /* 5 = D135_NODE */
  -D45_PRED,   14,           /* 6 = D45_NODE */
  -D63_PRED,   16,           /* 7 = D63_NODE */
  -D153_PRED,  18,           /* 8 = D153_NODE */
  -D207_PRED,  20,           /* 9 = D207_NODE */
  -PAETH_PRED, -SMOOTH_PRED, /* 10 = PAETH_NODE */
};
#else
const aom_tree_index av1_intra_mode_tree[TREE_SIZE(INTRA_MODES)] = {
  -DC_PRED,   2,          /* 0 = DC_NODE */
  -TM_PRED,   4,          /* 1 = TM_NODE */
  -V_PRED,    6,          /* 2 = V_NODE */
  8,          12,         /* 3 = COM_NODE */
  -H_PRED,    10,         /* 4 = H_NODE */
  -D135_PRED, -D117_PRED, /* 5 = D135_NODE */
  -D45_PRED,  14,         /* 6 = D45_NODE */
  -D63_PRED,  16,         /* 7 = D63_NODE */
  -D153_PRED, -D207_PRED  /* 8 = D153_NODE */
};
#endif  // CONFIG_ALT_INTRA

const aom_tree_index av1_inter_mode_tree[TREE_SIZE(INTER_MODES)] = {
  -INTER_OFFSET(ZEROMV),    2,
  -INTER_OFFSET(NEARESTMV), 4,
#if CONFIG_EXT_INTER
  -INTER_OFFSET(NEARMV),    6,
  -INTER_OFFSET(NEWMV),     -INTER_OFFSET(NEWFROMNEARMV)
#else
  -INTER_OFFSET(NEARMV),    -INTER_OFFSET(NEWMV)
#endif  // CONFIG_EXT_INTER
};

#if CONFIG_EXT_INTER
/* clang-format off */
const aom_tree_index av1_interintra_mode_tree[TREE_SIZE(INTERINTRA_MODES)] = {
  -II_DC_PRED, 2,                   /* 0 = II_DC_NODE     */
  -II_TM_PRED, 4,                   /* 1 = II_TM_NODE     */
  -II_V_PRED, 6,                    /* 2 = II_V_NODE      */
  8, 12,                            /* 3 = II_COM_NODE    */
  -II_H_PRED, 10,                   /* 4 = II_H_NODE      */
  -II_D135_PRED, -II_D117_PRED,     /* 5 = II_D135_NODE   */
  -II_D45_PRED, 14,                 /* 6 = II_D45_NODE    */
  -II_D63_PRED, 16,                 /* 7 = II_D63_NODE    */
  -II_D153_PRED, -II_D207_PRED      /* 8 = II_D153_NODE   */
};

const aom_tree_index av1_inter_compound_mode_tree
    [TREE_SIZE(INTER_COMPOUND_MODES)] = {
  -INTER_COMPOUND_OFFSET(ZERO_ZEROMV), 2,
  -INTER_COMPOUND_OFFSET(NEAREST_NEARESTMV), 4,
  6, -INTER_COMPOUND_OFFSET(NEW_NEWMV),
  8, 12,
  -INTER_COMPOUND_OFFSET(NEAR_NEARMV), 10,
  -INTER_COMPOUND_OFFSET(NEAREST_NEARMV),
      -INTER_COMPOUND_OFFSET(NEAR_NEARESTMV),
  14, 16,
  -INTER_COMPOUND_OFFSET(NEAREST_NEWMV), -INTER_COMPOUND_OFFSET(NEW_NEARESTMV),
  -INTER_COMPOUND_OFFSET(NEAR_NEWMV), -INTER_COMPOUND_OFFSET(NEW_NEARMV)
};

#if CONFIG_COMPOUND_SEGMENT
const aom_tree_index av1_compound_type_tree[TREE_SIZE(COMPOUND_TYPES)] = {
  -COMPOUND_AVERAGE, 2, -COMPOUND_WEDGE, -COMPOUND_SEG
};
#else  // !CONFIG_COMPOUND_SEGMENT
const aom_tree_index av1_compound_type_tree[TREE_SIZE(COMPOUND_TYPES)] = {
  -COMPOUND_AVERAGE, -COMPOUND_WEDGE
};
#endif  // CONFIG_COMPOUND_SEGMENT
/* clang-format on */
#endif  // CONFIG_EXT_INTER

const aom_tree_index av1_partition_tree[TREE_SIZE(PARTITION_TYPES)] = {
  -PARTITION_NONE, 2, -PARTITION_HORZ, 4, -PARTITION_VERT, -PARTITION_SPLIT
};

#if CONFIG_EXT_PARTITION_TYPES
/* clang-format off */
const aom_tree_index av1_ext_partition_tree[TREE_SIZE(EXT_PARTITION_TYPES)] = {
  -PARTITION_NONE, 2,
  6, 4,
  8, -PARTITION_SPLIT,
  -PARTITION_HORZ, 10,
  -PARTITION_VERT, 12,
  -PARTITION_HORZ_A, -PARTITION_HORZ_B,
  -PARTITION_VERT_A, -PARTITION_VERT_B
};
/* clang-format on */
#endif  // CONFIG_EXT_PARTITION_TYPES

static const aom_prob default_intra_inter_p[INTRA_INTER_CONTEXTS] = {
  9, 102, 187, 225
};

static const aom_prob default_comp_inter_p[COMP_INTER_CONTEXTS] = {
  239, 183, 119, 96, 41
};

#if CONFIG_EXT_REFS
static const aom_prob default_comp_ref_p[REF_CONTEXTS][FWD_REFS - 1] = {
  // TODO(zoeliu): To adjust the initial prob values.
  { 33, 16, 16 },
  { 77, 74, 74 },
  { 142, 142, 142 },
  { 172, 170, 170 },
  { 238, 247, 247 }
};
static const aom_prob default_comp_bwdref_p[REF_CONTEXTS][BWD_REFS - 1] = {
  { 16 }, { 74 }, { 142 }, { 170 }, { 247 }
};
#else
static const aom_prob default_comp_ref_p[REF_CONTEXTS][COMP_REFS - 1] = {
  { 50 }, { 126 }, { 123 }, { 221 }, { 226 }
};
#endif  // CONFIG_EXT_REFS

static const aom_prob default_single_ref_p[REF_CONTEXTS][SINGLE_REFS - 1] = {
#if CONFIG_EXT_REFS
  { 33, 16, 16, 16, 16 },
  { 77, 74, 74, 74, 74 },
  { 142, 142, 142, 142, 142 },
  { 172, 170, 170, 170, 170 },
  { 238, 247, 247, 247, 247 }
#else
  { 33, 16 }, { 77, 74 }, { 142, 142 }, { 172, 170 }, { 238, 247 }
#endif  // CONFIG_EXT_REFS
};

#if CONFIG_PALETTE

// Tree to code palette size (number of colors in a palette) and the
// corresponding probabilities for Y and UV planes.
const aom_tree_index av1_palette_size_tree[TREE_SIZE(PALETTE_SIZES)] = {
  -TWO_COLORS,  2, -THREE_COLORS, 4,  -FOUR_COLORS,  6,
  -FIVE_COLORS, 8, -SIX_COLORS,   10, -SEVEN_COLORS, -EIGHT_COLORS,
};

// TODO(huisu): tune these probs
const aom_prob
    av1_default_palette_y_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1] = {
      { 96, 89, 100, 64, 77, 130 },   { 22, 15, 44, 16, 34, 82 },
      { 30, 19, 57, 18, 38, 86 },     { 94, 36, 104, 23, 43, 92 },
      { 116, 76, 107, 46, 65, 105 },  { 112, 82, 94, 40, 70, 112 },
      { 147, 124, 123, 58, 69, 103 }, { 180, 113, 136, 49, 45, 114 },
      { 107, 70, 87, 49, 154, 156 },  { 98, 105, 142, 63, 64, 152 },
#if CONFIG_EXT_PARTITION
      { 98, 105, 142, 63, 64, 152 },  { 98, 105, 142, 63, 64, 152 },
      { 98, 105, 142, 63, 64, 152 },
#endif  // CONFIG_EXT_PARTITION
    };

const aom_prob
    av1_default_palette_uv_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1] = {
      { 160, 196, 228, 213, 175, 230 }, { 87, 148, 208, 141, 166, 163 },
      { 72, 151, 204, 139, 155, 161 },  { 78, 135, 171, 104, 120, 173 },
      { 59, 92, 131, 78, 92, 142 },     { 75, 118, 149, 84, 90, 128 },
      { 89, 87, 92, 66, 66, 128 },      { 67, 53, 54, 55, 66, 93 },
      { 120, 130, 83, 171, 75, 214 },   { 72, 55, 66, 68, 79, 107 },
#if CONFIG_EXT_PARTITION
      { 72, 55, 66, 68, 79, 107 },      { 72, 55, 66, 68, 79, 107 },
      { 72, 55, 66, 68, 79, 107 },
#endif  // CONFIG_EXT_PARTITION
    };

// When palette mode is enabled, following probability tables indicate the
// probabilities to code the "is_palette" bit (i.e. the bit that indicates
// if this block uses palette mode or DC_PRED mode).
const aom_prob
    av1_default_palette_y_mode_prob[PALETTE_BLOCK_SIZES]
                                   [PALETTE_Y_MODE_CONTEXTS] = {
                                     { 240, 180, 100 }, { 240, 180, 100 },
                                     { 240, 180, 100 }, { 240, 180, 100 },
                                     { 240, 180, 100 }, { 240, 180, 100 },
                                     { 240, 180, 100 }, { 240, 180, 100 },
                                     { 240, 180, 100 }, { 240, 180, 100 },
#if CONFIG_EXT_PARTITION
                                     { 240, 180, 100 }, { 240, 180, 100 },
                                     { 240, 180, 100 },
#endif  // CONFIG_EXT_PARTITION
                                   };

const aom_prob av1_default_palette_uv_mode_prob[PALETTE_UV_MODE_CONTEXTS] = {
  253, 229
};

// Trees to code palette color indices (for various palette sizes), and the
// corresponding probability tables for Y and UV planes.
const aom_tree_index
    av1_palette_color_index_tree[PALETTE_MAX_SIZE -
                                 1][TREE_SIZE(PALETTE_COLORS)] = {
      { // 2 colors
        -PALETTE_COLOR_ONE, -PALETTE_COLOR_TWO },
      { // 3 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, -PALETTE_COLOR_THREE },
      { // 4 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, 4, -PALETTE_COLOR_THREE,
        -PALETTE_COLOR_FOUR },
      { // 5 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, 4, -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, -PALETTE_COLOR_FIVE },
      { // 6 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, 4, -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8, -PALETTE_COLOR_FIVE, -PALETTE_COLOR_SIX },
      { // 7 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, 4, -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8, -PALETTE_COLOR_FIVE, 10, -PALETTE_COLOR_SIX,
        -PALETTE_COLOR_SEVEN },
      { // 8 colors
        -PALETTE_COLOR_ONE, 2, -PALETTE_COLOR_TWO, 4, -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8, -PALETTE_COLOR_FIVE, 10, -PALETTE_COLOR_SIX, 12,
        -PALETTE_COLOR_SEVEN, -PALETTE_COLOR_EIGHT },
    };

// Note: Has to be non-zero to avoid any asserts triggering.
#define UNUSED_PROB 128

const aom_prob av1_default_palette_y_color_index_prob
    [PALETTE_MAX_SIZE - 1][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS - 1] = {
      {
          // 2 colors
          { 231, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { 69, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 224, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { 249, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 3 colors
          { 219, 124, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 91, 191, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 34, 237, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 184, 118, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 252, 124, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
      },
      {
          // 4 colors
          { 204, 87, 97, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 74, 144, 129, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 52, 191, 134, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 151, 85, 147, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 248, 60, 115, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 5 colors
          { 218, 69, 62, 106, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 76, 143, 89, 127, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 21, 233, 94, 131, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 172, 72, 89, 112, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 253, 66, 65, 128, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 6 colors
          { 190, 60, 47, 54, 74, UNUSED_PROB, UNUSED_PROB },
          { 62, 106, 51, 95, 110, UNUSED_PROB, UNUSED_PROB },
          { 52, 180, 69, 72, 107, UNUSED_PROB, UNUSED_PROB },
          { 156, 83, 72, 83, 101, UNUSED_PROB, UNUSED_PROB },
          { 245, 45, 37, 52, 91, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 7 colors
          { 206, 56, 42, 42, 53, 85, UNUSED_PROB },
          { 70, 100, 45, 68, 77, 94, UNUSED_PROB },
          { 57, 169, 51, 62, 74, 119, UNUSED_PROB },
          { 172, 76, 71, 40, 59, 76, UNUSED_PROB },
          { 248, 47, 36, 53, 61, 110, UNUSED_PROB },
      },
      {
          // 8 colors
          { 208, 52, 38, 34, 34, 44, 66 },
          { 52, 107, 34, 73, 69, 82, 87 },
          { 28, 208, 53, 43, 62, 70, 102 },
          { 184, 64, 45, 37, 37, 69, 105 },
          { 251, 18, 31, 45, 47, 61, 104 },
      },
    };

const aom_prob av1_default_palette_uv_color_index_prob
    [PALETTE_MAX_SIZE - 1][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS - 1] = {
      {
          // 2 colors
          { 233, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { 69, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 240, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
          { 248, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 3 colors
          { 216, 128, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 110, 171, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 40, 239, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 191, 104, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
          { 247, 134, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB,
            UNUSED_PROB },
      },
      {
          // 4 colors
          { 202, 89, 132, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 90, 132, 136, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 63, 195, 149, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 152, 84, 152, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 241, 87, 136, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 5 colors
          { 209, 54, 82, 134, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 94, 173, 180, 93, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 10, 251, 127, 84, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 183, 20, 150, 47, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
          { 252, 73, 111, 150, UNUSED_PROB, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 6 colors
          { 192, 67, 59, 46, 184, UNUSED_PROB, UNUSED_PROB },
          { 59, 92, 61, 100, 130, UNUSED_PROB, UNUSED_PROB },
          { 49, 162, 68, 91, 150, UNUSED_PROB, UNUSED_PROB },
          { 133, 29, 36, 153, 101, UNUSED_PROB, UNUSED_PROB },
          { 247, 71, 44, 90, 129, UNUSED_PROB, UNUSED_PROB },
      },
      {
          // 7 colors
          { 182, 62, 80, 78, 46, 116, UNUSED_PROB },
          { 59, 62, 39, 81, 65, 99, UNUSED_PROB },
          { 54, 177, 48, 58, 93, 104, UNUSED_PROB },
          { 137, 79, 54, 55, 44, 134, UNUSED_PROB },
          { 239, 82, 79, 44, 69, 71, UNUSED_PROB },
      },
      {
          // 8 colors
          { 172, 53, 27, 67, 30, 79, 113 },
          { 63, 57, 45, 81, 62, 35, 47 },
          { 51, 200, 36, 47, 82, 165, 129 },
          { 141, 100, 47, 29, 33, 37, 129 },
          { 236, 42, 50, 91, 24, 154, 65 },
      },
    };

#undef UNUSED_PROB

#define MAX_COLOR_CONTEXT_HASH 8
// Negative values are invalid
static const int palette_color_index_context_lookup[MAX_COLOR_CONTEXT_HASH +
                                                    1] = { -1, -1, 0, -1, -1,
                                                           4,  3,  2, 1 };

#endif  // CONFIG_PALETTE

// The transform size is coded as an offset to the smallest transform
// block size.
const aom_tree_index av1_tx_size_tree[MAX_TX_DEPTH][TREE_SIZE(TX_SIZES)] = {
  {
      // Max tx_size is 8X8
      -0, -1,
  },
  {
      // Max tx_size is 16X16
      -0, 2, -1, -2,
  },
  {
      // Max tx_size is 32X32
      -0, 2, -1, 4, -2, -3,
  },
#if CONFIG_TX64X64
  {
      // Max tx_size is 64X64
      -0, 2, -1, 4, -2, 6, -3, -4,
  },
#endif  // CONFIG_TX64X64
};

static const aom_prob default_tx_size_prob[MAX_TX_DEPTH][TX_SIZE_CONTEXTS]
                                          [MAX_TX_DEPTH] = {
                                            {
                                                // Max tx_size is 8X8
                                                { 100 },
                                                { 66 },
                                            },
                                            {
                                                // Max tx_size is 16X16
                                                { 20, 152 },
                                                { 15, 101 },
                                            },
                                            {
                                                // Max tx_size is 32X32
                                                { 3, 136, 37 },
                                                { 5, 52, 13 },
                                            },
#if CONFIG_TX64X64
                                            {
                                                // Max tx_size is 64X64
                                                { 1, 64, 136, 127 },
                                                { 1, 32, 52, 67 },
                                            },
#endif  // CONFIG_TX64X64
                                          };

#if CONFIG_LOOP_RESTORATION
const aom_tree_index av1_switchable_restore_tree[TREE_SIZE(
    RESTORE_SWITCHABLE_TYPES)] = {
  -RESTORE_NONE, 2, -RESTORE_WIENER, 4, -RESTORE_SGRPROJ, -RESTORE_DOMAINTXFMRF,
};

static const aom_prob
    default_switchable_restore_prob[RESTORE_SWITCHABLE_TYPES - 1] = {
      32, 128, 128,
    };
#endif  // CONFIG_LOOP_RESTORATION

#if CONFIG_PALETTE
#define NUM_PALETTE_NEIGHBORS 3  // left, top-left and top.
int av1_get_palette_color_index_context(const uint8_t *color_map, int stride,
                                        int r, int c, int palette_size,
                                        uint8_t *color_order, int *color_idx) {
  int i;
  // The +10 below should not be needed. But we get a warning "array subscript
  // is above array bounds [-Werror=array-bounds]" without it, possibly due to
  // this (or similar) bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59124
  int scores[PALETTE_MAX_SIZE + 10];
  const int weights[NUM_PALETTE_NEIGHBORS] = { 2, 1, 2 };
  const int hash_multipliers[NUM_PALETTE_NEIGHBORS] = { 1, 2, 2 };
  int color_index_ctx_hash;
  int color_index_ctx;
  int color_neighbors[NUM_PALETTE_NEIGHBORS];
  int inverse_color_order[PALETTE_MAX_SIZE];
  assert(palette_size <= PALETTE_MAX_SIZE);
  assert(r > 0 || c > 0);

  // Get color indices of neighbors.
  color_neighbors[0] = (c - 1 >= 0) ? color_map[r * stride + c - 1] : -1;
  color_neighbors[1] =
      (c - 1 >= 0 && r - 1 >= 0) ? color_map[(r - 1) * stride + c - 1] : -1;
  color_neighbors[2] = (r - 1 >= 0) ? color_map[(r - 1) * stride + c] : -1;

  for (i = 0; i < PALETTE_MAX_SIZE; ++i) {
    color_order[i] = i;
    inverse_color_order[i] = i;
  }
  memset(scores, 0, PALETTE_MAX_SIZE * sizeof(scores[0]));
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    if (color_neighbors[i] >= 0) {
      scores[color_neighbors[i]] += weights[i];
    }
  }

  // Get the top NUM_PALETTE_NEIGHBORS scores (sorted from large to small).
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    int max = scores[i];
    int max_idx = i;
    int j;
    for (j = i + 1; j < palette_size; ++j) {
      if (scores[j] > max) {
        max = scores[j];
        max_idx = j;
      }
    }
    if (max_idx != i) {
      // Move the score at index 'max_idx' to index 'i', and shift the scores
      // from 'i' to 'max_idx - 1' by 1.
      const int max_score = scores[max_idx];
      const uint8_t max_color_order = color_order[max_idx];
      int k;
      for (k = max_idx; k > i; --k) {
        scores[k] = scores[k - 1];
        color_order[k] = color_order[k - 1];
        inverse_color_order[color_order[k]] = k;
      }
      scores[i] = max_score;
      color_order[i] = max_color_order;
      inverse_color_order[color_order[i]] = i;
    }
  }

  // Get hash value of context.
  color_index_ctx_hash = 0;
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    color_index_ctx_hash += scores[i] * hash_multipliers[i];
  }
  assert(color_index_ctx_hash > 0);
  assert(color_index_ctx_hash <= MAX_COLOR_CONTEXT_HASH);

  // Lookup context from hash.
  color_index_ctx = palette_color_index_context_lookup[color_index_ctx_hash];
  assert(color_index_ctx >= 0);
  assert(color_index_ctx < PALETTE_COLOR_INDEX_CONTEXTS);

  if (color_idx != NULL) {
    *color_idx = inverse_color_order[color_map[r * stride + c]];
  }
  return color_index_ctx;
}
#undef NUM_PALETTE_NEIGHBORS
#undef MAX_COLOR_CONTEXT_HASH

#endif  // CONFIG_PALETTE

#if CONFIG_VAR_TX
static const aom_prob default_txfm_partition_probs[TXFM_PARTITION_CONTEXTS] = {
  250, 231, 212, 241, 166, 66, 241, 230, 135, 243, 154, 64, 248, 161, 63, 128,
};
#endif

static const aom_prob default_skip_probs[SKIP_CONTEXTS] = { 192, 128, 64 };

#if CONFIG_DUAL_FILTER
static const aom_prob
    default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                  [SWITCHABLE_FILTERS - 1] = {
                                    { 235, 192, 128 }, { 36, 243, 48 },
                                    { 34, 16, 128 },   { 34, 16, 128 },
                                    { 149, 160, 128 }, { 235, 192, 128 },
                                    { 36, 243, 48 },   { 34, 16, 128 },
                                    { 34, 16, 128 },   { 149, 160, 128 },
                                    { 235, 192, 128 }, { 36, 243, 48 },
                                    { 34, 16, 128 },   { 34, 16, 128 },
                                    { 149, 160, 128 }, { 235, 192, 128 },
                                    { 36, 243, 48 },   { 34, 16, 128 },
                                    { 34, 16, 128 },   { 149, 160, 128 },
                                  };
#else   // CONFIG_DUAL_FILTER
static const aom_prob default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                                    [SWITCHABLE_FILTERS - 1] = {
                                                      { 235, 162 },
                                                      { 36, 255 },
                                                      { 34, 3 },
                                                      { 149, 144 },
                                                    };
#endif  // CONFIG_DUAL_FILTER

#if CONFIG_EXT_TX
/* clang-format off */
const aom_tree_index av1_ext_tx_inter_tree[EXT_TX_SETS_INTER]
                                           [TREE_SIZE(TX_TYPES)] = {
  { // ToDo(yaowu): remove used entry 0.
    0
  }, {
    -IDTX, 2,
    4, 14,
    6, 8,
    -V_DCT, -H_DCT,
    10, 12,
    -V_ADST, -H_ADST,
    -V_FLIPADST, -H_FLIPADST,
    -DCT_DCT, 16,
    18, 24,
    20, 22,
    -ADST_DCT, -DCT_ADST,
    -FLIPADST_DCT, -DCT_FLIPADST,
    26, 28,
    -ADST_ADST, -FLIPADST_FLIPADST,
    -ADST_FLIPADST, -FLIPADST_ADST
  }, {
    -IDTX, 2,
    4, 6,
    -V_DCT, -H_DCT,
    -DCT_DCT, 8,
    10, 16,
    12, 14,
    -ADST_DCT, -DCT_ADST,
    -FLIPADST_DCT, -DCT_FLIPADST,
    18, 20,
    -ADST_ADST, -FLIPADST_FLIPADST,
    -ADST_FLIPADST, -FLIPADST_ADST
  }, {
    -IDTX, -DCT_DCT,
  }
};

const aom_tree_index av1_ext_tx_intra_tree[EXT_TX_SETS_INTRA]
                                           [TREE_SIZE(TX_TYPES)] = {
  {  // ToDo(yaowu): remove unused entry 0.
    0
  }, {
    -IDTX, 2,
    -DCT_DCT, 4,
    6, 8,
    -V_DCT, -H_DCT,
    -ADST_ADST, 10,
    -ADST_DCT, -DCT_ADST,
  }, {
    -IDTX, 2,
    -DCT_DCT, 4,
    -ADST_ADST, 6,
    -ADST_DCT, -DCT_ADST,
  }
};
/* clang-format on */

static const aom_prob
    default_inter_ext_tx_prob[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES - 1] = {
      {
// ToDo(yaowu): remove unused entry 0.
#if CONFIG_CB4X4
          { 0 },
#endif
          { 0 },
          { 0 },
          { 0 },
          { 0 },
      },
      {
#if CONFIG_CB4X4
          { 0 },
#endif
          { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128,
            128 },
          { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128,
            128 },
          { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128,
            128 },
          { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128,
            128 },
      },
      {
#if CONFIG_CB4X4
          { 0 },
#endif
          { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
          { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
          { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
          { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
      },
      {
#if CONFIG_CB4X4
          { 0 },
#endif
          { 12 },
          { 12 },
          { 12 },
          { 12 },
      }
    };

static const aom_prob
    default_intra_ext_tx_prob[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                             [TX_TYPES - 1] = {
                               {
// ToDo(yaowu): remove unused entry 0.
#if CONFIG_CB4X4
                                   {
                                       { 0 },
                                   },
#endif
                                   {
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                   },
                                   {
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                   },
                                   {
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                   },
                                   {
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                       { 0 },
                                   },
                               },
                               {
#if CONFIG_CB4X4
                                   {
                                       { 0 },
                                   },
#endif
                                   {
                                       { 8, 224, 32, 128, 64, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 9, 200, 32, 128, 64, 128 },
                                       { 8, 8, 32, 128, 224, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 32, 32, 128, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 32, 128, 64, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 9, 200, 32, 128, 64, 128 },
                                       { 8, 8, 32, 128, 224, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 32, 32, 128, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 32, 128, 64, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 9, 200, 32, 128, 64, 128 },
                                       { 8, 8, 32, 128, 224, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 32, 32, 128, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 32, 128, 64, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 9, 200, 32, 128, 64, 128 },
                                       { 8, 8, 32, 128, 224, 128 },
                                       { 10, 32, 32, 128, 16, 192 },
                                       { 10, 32, 32, 128, 16, 64 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 23, 32, 128, 80, 176 },
                                       { 10, 32, 32, 128, 16, 64 },
                                   },
                               },
                               {
#if CONFIG_CB4X4
                                   {
                                       { 0 },
                                   },
#endif
                                   {
                                       { 8, 224, 64, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 9, 200, 64, 128 },
                                       { 8, 8, 224, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 10, 23, 80, 176 },
                                       { 10, 23, 80, 176 },
                                       { 10, 32, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 64, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 9, 200, 64, 128 },
                                       { 8, 8, 224, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 10, 23, 80, 176 },
                                       { 10, 23, 80, 176 },
                                       { 10, 32, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 64, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 9, 200, 64, 128 },
                                       { 8, 8, 224, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 10, 23, 80, 176 },
                                       { 10, 23, 80, 176 },
                                       { 10, 32, 16, 64 },
                                   },
                                   {
                                       { 8, 224, 64, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 9, 200, 64, 128 },
                                       { 8, 8, 224, 128 },
                                       { 10, 32, 16, 192 },
                                       { 10, 32, 16, 64 },
                                       { 10, 23, 80, 176 },
                                       { 10, 23, 80, 176 },
                                       { 10, 32, 16, 64 },
                                   },
                               },
                             };
#else

/* clang-format off */
const aom_tree_index av1_ext_tx_tree[TREE_SIZE(TX_TYPES)] = {
  -DCT_DCT, 2,
  -ADST_ADST, 4,
  -ADST_DCT, -DCT_ADST
};
/* clang-format on */

int av1_ext_tx_ind[TX_TYPES];
int av1_ext_tx_inv[TX_TYPES];

static const aom_prob
    default_intra_ext_tx_prob[EXT_TX_SIZES][TX_TYPES][TX_TYPES - 1] = {
#if CONFIG_CB4X4
      { { 240, 85, 128 }, { 4, 1, 248 }, { 4, 1, 8 }, { 4, 248, 128 } },
#endif
      { { 240, 85, 128 }, { 4, 1, 248 }, { 4, 1, 8 }, { 4, 248, 128 } },
      { { 244, 85, 128 }, { 8, 2, 248 }, { 8, 2, 8 }, { 8, 248, 128 } },
      { { 248, 85, 128 }, { 16, 4, 248 }, { 16, 4, 8 }, { 16, 248, 128 } },
    };

static const aom_prob default_inter_ext_tx_prob[EXT_TX_SIZES][TX_TYPES - 1] = {
#if CONFIG_CB4X4
  { 160, 85, 128 },
#endif
  { 160, 85, 128 },
  { 176, 85, 128 },
  { 192, 85, 128 },
};
#endif  // CONFIG_EXT_TX

#if CONFIG_EXT_INTRA
#if CONFIG_INTRA_INTERP
static const aom_prob default_intra_filter_probs[INTRA_FILTERS + 1]
                                                [INTRA_FILTERS - 1] = {
                                                  { 98, 63, 60 },
                                                  { 98, 82, 80 },
                                                  { 94, 65, 103 },
                                                  { 49, 25, 24 },
                                                  { 72, 38, 50 },
                                                };
const aom_tree_index av1_intra_filter_tree[TREE_SIZE(INTRA_FILTERS)] = {
  -INTRA_FILTER_LINEAR,      2, -INTRA_FILTER_8TAP, 4, -INTRA_FILTER_8TAP_SHARP,
  -INTRA_FILTER_8TAP_SMOOTH,
};
#endif  // CONFIG_INTRA_INTERP
#endif  // CONFIG_EXT_INTRA

#if CONFIG_FILTER_INTRA
static const aom_prob default_filter_intra_probs[2] = { 230, 230 };
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_SUPERTX
static const aom_prob default_supertx_prob[PARTITION_SUPERTX_CONTEXTS]
                                          [TX_SIZES] = {
#if CONFIG_CB4X4
#if CONFIG_TX64X64
                                            { 1, 1, 160, 160, 170, 180 },
                                            { 1, 1, 200, 200, 210, 220 },
#else
                                            { 1, 1, 160, 160, 170 },
                                            { 1, 1, 200, 200, 210 },
#endif  // CONFIG_TX64X64
#else
#if CONFIG_TX64X64
                                            { 1, 160, 160, 170, 180 },
                                            { 1, 200, 200, 210, 220 },
#else
                                            { 1, 160, 160, 170 },
                                            { 1, 200, 200, 210 },
#endif  // CONFIG_CB4X4
#endif  // CONFIG_TX64X64
                                          };
#endif  // CONFIG_SUPERTX

// FIXME(someone) need real defaults here
static const aom_prob default_segment_tree_probs[SEG_TREE_PROBS] = {
  128, 128, 128, 128, 128, 128, 128
};
// clang-format off
static const aom_prob default_segment_pred_probs[PREDICTION_PROBS] = {
  128, 128, 128
};
// clang-format on

static void init_mode_probs(FRAME_CONTEXT *fc) {
  av1_copy(fc->uv_mode_prob, default_uv_probs);
  av1_copy(fc->y_mode_prob, default_if_y_probs);
  av1_copy(fc->switchable_interp_prob, default_switchable_interp_prob);
  av1_copy(fc->partition_prob, default_partition_probs);
  av1_copy(fc->intra_inter_prob, default_intra_inter_p);
  av1_copy(fc->comp_inter_prob, default_comp_inter_p);
  av1_copy(fc->comp_ref_prob, default_comp_ref_p);
#if CONFIG_EXT_REFS
  av1_copy(fc->comp_bwdref_prob, default_comp_bwdref_p);
#endif  // CONFIG_EXT_REFS
  av1_copy(fc->single_ref_prob, default_single_ref_p);
  av1_copy(fc->tx_size_probs, default_tx_size_prob);
#if CONFIG_VAR_TX
  av1_copy(fc->txfm_partition_prob, default_txfm_partition_probs);
#endif
  av1_copy(fc->skip_probs, default_skip_probs);
#if CONFIG_REF_MV
  av1_copy(fc->newmv_prob, default_newmv_prob);
  av1_copy(fc->zeromv_prob, default_zeromv_prob);
  av1_copy(fc->refmv_prob, default_refmv_prob);
  av1_copy(fc->drl_prob, default_drl_prob);
#if CONFIG_EXT_INTER
  fc->new2mv_prob = default_new2mv_prob;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV
  av1_copy(fc->inter_mode_probs, default_inter_mode_probs);
#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  av1_copy(fc->motion_mode_prob, default_motion_mode_prob);
#if CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
  av1_copy(fc->obmc_prob, default_obmc_prob);
#endif  // CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
#if CONFIG_EXT_INTER
  av1_copy(fc->inter_compound_mode_probs, default_inter_compound_mode_probs);
  av1_copy(fc->compound_type_prob, default_compound_type_probs);
  av1_copy(fc->interintra_prob, default_interintra_prob);
  av1_copy(fc->interintra_mode_prob, default_interintra_mode_prob);
  av1_copy(fc->wedge_interintra_prob, default_wedge_interintra_prob);
#endif  // CONFIG_EXT_INTER
#if CONFIG_SUPERTX
  av1_copy(fc->supertx_prob, default_supertx_prob);
#endif  // CONFIG_SUPERTX
  av1_copy(fc->seg.tree_probs, default_segment_tree_probs);
  av1_copy(fc->seg.pred_probs, default_segment_pred_probs);
#if CONFIG_EXT_INTRA
#if CONFIG_INTRA_INTERP
  av1_copy(fc->intra_filter_probs, default_intra_filter_probs);
#endif  // CONFIG_INTRA_INTERP
#endif  // CONFIG_EXT_INTRA
#if CONFIG_FILTER_INTRA
  av1_copy(fc->filter_intra_probs, default_filter_intra_probs);
#endif  // CONFIG_FILTER_INTRA
  av1_copy(fc->inter_ext_tx_prob, default_inter_ext_tx_prob);
  av1_copy(fc->intra_ext_tx_prob, default_intra_ext_tx_prob);
#if CONFIG_LOOP_RESTORATION
  av1_copy(fc->switchable_restore_prob, default_switchable_restore_prob);
#endif  // CONFIG_LOOP_RESTORATION
#if CONFIG_EC_MULTISYMBOL
  av1_tree_to_cdf_1D(av1_intra_mode_tree, fc->y_mode_prob, fc->y_mode_cdf,
                     BLOCK_SIZE_GROUPS);
  av1_tree_to_cdf_1D(av1_intra_mode_tree, fc->uv_mode_prob, fc->uv_mode_cdf,
                     INTRA_MODES);
  av1_tree_to_cdf_1D(av1_switchable_interp_tree, fc->switchable_interp_prob,
                     fc->switchable_interp_cdf, SWITCHABLE_FILTER_CONTEXTS);
  av1_tree_to_cdf_1D(av1_partition_tree, fc->partition_prob, fc->partition_cdf,
                     PARTITION_CONTEXTS);
  av1_tree_to_cdf_1D(av1_inter_mode_tree, fc->inter_mode_probs,
                     fc->inter_mode_cdf, INTER_MODE_CONTEXTS);
#if !CONFIG_EXT_TX
  av1_tree_to_cdf_2D(av1_ext_tx_tree, fc->intra_ext_tx_prob,
                     fc->intra_ext_tx_cdf, EXT_TX_SIZES, TX_TYPES);
  av1_tree_to_cdf_1D(av1_ext_tx_tree, fc->inter_ext_tx_prob,
                     fc->inter_ext_tx_cdf, EXT_TX_SIZES);
#endif
  av1_tree_to_cdf_2D(av1_intra_mode_tree, av1_kf_y_mode_prob, av1_kf_y_mode_cdf,
                     INTRA_MODES, INTRA_MODES);
  av1_tree_to_cdf(av1_segment_tree, fc->seg.tree_probs, fc->seg.tree_cdf);
#endif
#if CONFIG_DELTA_Q
  av1_copy(fc->delta_q_prob, default_delta_q_probs);
#endif
}

#if CONFIG_EC_MULTISYMBOL
int av1_switchable_interp_ind[SWITCHABLE_FILTERS];
int av1_switchable_interp_inv[SWITCHABLE_FILTERS];

void av1_set_mode_cdfs(struct AV1Common *cm) {
  FRAME_CONTEXT *fc = cm->fc;
  int i, j;
  if (cm->seg.enabled && cm->seg.update_map) {
    av1_tree_to_cdf(av1_segment_tree, cm->fc->seg.tree_probs,
                    cm->fc->seg.tree_cdf);
  }

  for (i = 0; i < INTRA_MODES; ++i)
    av1_tree_to_cdf(av1_intra_mode_tree, fc->uv_mode_prob[i],
                    fc->uv_mode_cdf[i]);

  for (i = 0; i < PARTITION_CONTEXTS; ++i)
    av1_tree_to_cdf(av1_partition_tree, fc->partition_prob[i],
                    fc->partition_cdf[i]);

  for (i = 0; i < INTRA_MODES; ++i)
    for (j = 0; j < INTRA_MODES; ++j)
      av1_tree_to_cdf(av1_intra_mode_tree, cm->kf_y_prob[i][j],
                      cm->fc->kf_y_cdf[i][j]);

  for (j = 0; j < SWITCHABLE_FILTER_CONTEXTS; ++j)
    av1_tree_to_cdf(av1_switchable_interp_tree, fc->switchable_interp_prob[j],
                    fc->switchable_interp_cdf[j]);

  for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
    av1_tree_to_cdf(av1_inter_mode_tree, fc->inter_mode_probs[i],
                    fc->inter_mode_cdf[i]);

  for (i = 0; i < BLOCK_SIZE_GROUPS; ++i)
    av1_tree_to_cdf(av1_intra_mode_tree, fc->y_mode_prob[i], fc->y_mode_cdf[i]);

#if !CONFIG_EXT_TX
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i)
    for (j = 0; j < TX_TYPES; ++j)
      av1_tree_to_cdf(av1_ext_tx_tree, fc->intra_ext_tx_prob[i][j],
                      fc->intra_ext_tx_cdf[i][j]);

  for (i = TX_4X4; i < EXT_TX_SIZES; ++i)
    av1_tree_to_cdf(av1_ext_tx_tree, fc->inter_ext_tx_prob[i],
                    fc->inter_ext_tx_cdf[i]);
#endif
}
#endif

#if CONFIG_DUAL_FILTER
const aom_tree_index av1_switchable_interp_tree[TREE_SIZE(SWITCHABLE_FILTERS)] =
    {
      -EIGHTTAP_REGULAR, 2, 4, -MULTITAP_SHARP, -EIGHTTAP_SMOOTH,
      -EIGHTTAP_SMOOTH2,
    };
#else
const aom_tree_index av1_switchable_interp_tree[TREE_SIZE(SWITCHABLE_FILTERS)] =
    { -EIGHTTAP_REGULAR, 2, -EIGHTTAP_SMOOTH, -MULTITAP_SHARP };
#endif  // CONFIG_DUAL_FILTER

void av1_adapt_inter_frame_probs(AV1_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = cm->fc;
  const FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];
  const FRAME_COUNTS *counts = &cm->counts;

  for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
    fc->intra_inter_prob[i] = av1_mode_mv_merge_probs(
        pre_fc->intra_inter_prob[i], counts->intra_inter[i]);
  for (i = 0; i < COMP_INTER_CONTEXTS; i++)
    fc->comp_inter_prob[i] = av1_mode_mv_merge_probs(pre_fc->comp_inter_prob[i],
                                                     counts->comp_inter[i]);

#if CONFIG_EXT_REFS
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (FWD_REFS - 1); j++)
      fc->comp_ref_prob[i][j] = mode_mv_merge_probs(pre_fc->comp_ref_prob[i][j],
                                                    counts->comp_ref[i][j]);
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (BWD_REFS - 1); j++)
      fc->comp_bwdref_prob[i][j] = mode_mv_merge_probs(
          pre_fc->comp_bwdref_prob[i][j], counts->comp_bwdref[i][j]);
#else
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (COMP_REFS - 1); j++)
      fc->comp_ref_prob[i][j] = mode_mv_merge_probs(pre_fc->comp_ref_prob[i][j],
                                                    counts->comp_ref[i][j]);
#endif  // CONFIG_EXT_REFS

  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (SINGLE_REFS - 1); j++)
      fc->single_ref_prob[i][j] = av1_mode_mv_merge_probs(
          pre_fc->single_ref_prob[i][j], counts->single_ref[i][j]);

#if CONFIG_REF_MV
  for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i)
    fc->newmv_prob[i] =
        av1_mode_mv_merge_probs(pre_fc->newmv_prob[i], counts->newmv_mode[i]);
  for (i = 0; i < ZEROMV_MODE_CONTEXTS; ++i)
    fc->zeromv_prob[i] =
        av1_mode_mv_merge_probs(pre_fc->zeromv_prob[i], counts->zeromv_mode[i]);
  for (i = 0; i < REFMV_MODE_CONTEXTS; ++i)
    fc->refmv_prob[i] =
        av1_mode_mv_merge_probs(pre_fc->refmv_prob[i], counts->refmv_mode[i]);

  for (i = 0; i < DRL_MODE_CONTEXTS; ++i)
    fc->drl_prob[i] =
        av1_mode_mv_merge_probs(pre_fc->drl_prob[i], counts->drl_mode[i]);
#if CONFIG_EXT_INTER
  fc->new2mv_prob =
      av1_mode_mv_merge_probs(pre_fc->new2mv_prob, counts->new2mv_mode);
#endif  // CONFIG_EXT_INTER
#else
  for (i = 0; i < INTER_MODE_CONTEXTS; i++)
    aom_tree_merge_probs(av1_inter_mode_tree, pre_fc->inter_mode_probs[i],
                         counts->inter_mode[i], fc->inter_mode_probs[i]);
#endif

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  for (i = BLOCK_8X8; i < BLOCK_SIZES; ++i)
    aom_tree_merge_probs(av1_motion_mode_tree, pre_fc->motion_mode_prob[i],
                         counts->motion_mode[i], fc->motion_mode_prob[i]);
#if CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
  for (i = BLOCK_8X8; i < BLOCK_SIZES; ++i)
    fc->obmc_prob[i] =
        av1_mode_mv_merge_probs(pre_fc->obmc_prob[i], counts->obmc[i]);
#endif  // CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

#if CONFIG_SUPERTX
  for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
    for (j = TX_8X8; j < TX_SIZES; ++j) {
      fc->supertx_prob[i][j] = av1_mode_mv_merge_probs(
          pre_fc->supertx_prob[i][j], counts->supertx[i][j]);
    }
  }
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_INTER
  for (i = 0; i < INTER_MODE_CONTEXTS; i++)
    aom_tree_merge_probs(
        av1_inter_compound_mode_tree, pre_fc->inter_compound_mode_probs[i],
        counts->inter_compound_mode[i], fc->inter_compound_mode_probs[i]);
  for (i = 0; i < BLOCK_SIZE_GROUPS; ++i) {
    if (is_interintra_allowed_bsize_group(i))
      fc->interintra_prob[i] = av1_mode_mv_merge_probs(
          pre_fc->interintra_prob[i], counts->interintra[i]);
  }
  for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
    aom_tree_merge_probs(
        av1_interintra_mode_tree, pre_fc->interintra_mode_prob[i],
        counts->interintra_mode[i], fc->interintra_mode_prob[i]);
  }
  for (i = 0; i < BLOCK_SIZES; ++i) {
    if (is_interintra_allowed_bsize(i) && is_interintra_wedge_used(i))
      fc->wedge_interintra_prob[i] = av1_mode_mv_merge_probs(
          pre_fc->wedge_interintra_prob[i], counts->wedge_interintra[i]);
  }

  for (i = 0; i < BLOCK_SIZES; ++i) {
    aom_tree_merge_probs(av1_compound_type_tree, pre_fc->compound_type_prob[i],
                         counts->compound_interinter[i],
                         fc->compound_type_prob[i]);
  }
#endif  // CONFIG_EXT_INTER

  for (i = 0; i < BLOCK_SIZE_GROUPS; i++)
    aom_tree_merge_probs(av1_intra_mode_tree, pre_fc->y_mode_prob[i],
                         counts->y_mode[i], fc->y_mode_prob[i]);

  if (cm->interp_filter == SWITCHABLE) {
    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; i++)
      aom_tree_merge_probs(
          av1_switchable_interp_tree, pre_fc->switchable_interp_prob[i],
          counts->switchable_interp[i], fc->switchable_interp_prob[i]);
  }

#if CONFIG_DELTA_Q
  for (i = 0; i < DELTA_Q_CONTEXTS; ++i)
    fc->delta_q_prob[i] =
        mode_mv_merge_probs(pre_fc->delta_q_prob[i], counts->delta_q[i]);
#endif
}

void av1_adapt_intra_frame_probs(AV1_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = cm->fc;
  const FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];
  const FRAME_COUNTS *counts = &cm->counts;

  if (cm->tx_mode == TX_MODE_SELECT) {
    for (i = 0; i < MAX_TX_DEPTH; ++i) {
      for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
        aom_tree_merge_probs(av1_tx_size_tree[i], pre_fc->tx_size_probs[i][j],
                             counts->tx_size[i][j], fc->tx_size_probs[i][j]);
    }
  }

#if CONFIG_VAR_TX
  if (cm->tx_mode == TX_MODE_SELECT) {
    for (i = 0; i < TXFM_PARTITION_CONTEXTS; ++i)
      fc->txfm_partition_prob[i] = av1_mode_mv_merge_probs(
          pre_fc->txfm_partition_prob[i], counts->txfm_partition[i]);
  }
#endif

  for (i = 0; i < SKIP_CONTEXTS; ++i)
    fc->skip_probs[i] =
        av1_mode_mv_merge_probs(pre_fc->skip_probs[i], counts->skip[i]);

#if CONFIG_EXT_TX
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    int s;
    for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
      if (use_inter_ext_tx_for_txsize[s][i]) {
        aom_tree_merge_probs(
            av1_ext_tx_inter_tree[s], pre_fc->inter_ext_tx_prob[s][i],
            counts->inter_ext_tx[s][i], fc->inter_ext_tx_prob[s][i]);
      }
    }
    for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
      if (use_intra_ext_tx_for_txsize[s][i]) {
        for (j = 0; j < INTRA_MODES; ++j)
          aom_tree_merge_probs(
              av1_ext_tx_intra_tree[s], pre_fc->intra_ext_tx_prob[s][i][j],
              counts->intra_ext_tx[s][i][j], fc->intra_ext_tx_prob[s][i][j]);
      }
    }
  }
#else
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    for (j = 0; j < TX_TYPES; ++j) {
      aom_tree_merge_probs(av1_ext_tx_tree, pre_fc->intra_ext_tx_prob[i][j],
                           counts->intra_ext_tx[i][j],
                           fc->intra_ext_tx_prob[i][j]);
    }
  }
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    aom_tree_merge_probs(av1_ext_tx_tree, pre_fc->inter_ext_tx_prob[i],
                         counts->inter_ext_tx[i], fc->inter_ext_tx_prob[i]);
  }
#endif  // CONFIG_EXT_TX

  if (cm->seg.temporal_update) {
    for (i = 0; i < PREDICTION_PROBS; i++)
      fc->seg.pred_probs[i] = av1_mode_mv_merge_probs(pre_fc->seg.pred_probs[i],
                                                      counts->seg.pred[i]);

    aom_tree_merge_probs(av1_segment_tree, pre_fc->seg.tree_probs,
                         counts->seg.tree_mispred, fc->seg.tree_probs);
  } else {
    aom_tree_merge_probs(av1_segment_tree, pre_fc->seg.tree_probs,
                         counts->seg.tree_total, fc->seg.tree_probs);
  }

  for (i = 0; i < INTRA_MODES; ++i)
    aom_tree_merge_probs(av1_intra_mode_tree, pre_fc->uv_mode_prob[i],
                         counts->uv_mode[i], fc->uv_mode_prob[i]);

#if CONFIG_EXT_PARTITION_TYPES
  aom_tree_merge_probs(av1_partition_tree, pre_fc->partition_prob[0],
                       counts->partition[0], fc->partition_prob[0]);
  for (i = 1; i < PARTITION_CONTEXTS_PRIMARY; i++)
    aom_tree_merge_probs(av1_ext_partition_tree, pre_fc->partition_prob[i],
                         counts->partition[i], fc->partition_prob[i]);
#else
  for (i = 0; i < PARTITION_CONTEXTS_PRIMARY; i++) {
    aom_tree_merge_probs(av1_partition_tree, pre_fc->partition_prob[i],
                         counts->partition[i], fc->partition_prob[i]);
  }
#endif  // CONFIG_EXT_PARTITION_TYPES
#if CONFIG_UNPOISON_PARTITION_CTX
  for (i = PARTITION_CONTEXTS_PRIMARY;
       i < PARTITION_CONTEXTS_PRIMARY + PARTITION_BLOCK_SIZES; ++i) {
    unsigned int ct[2] = { counts->partition[i][PARTITION_VERT],
                           counts->partition[i][PARTITION_SPLIT] };
    assert(counts->partition[i][PARTITION_NONE] == 0);
    assert(counts->partition[i][PARTITION_HORZ] == 0);
    assert(fc->partition_prob[i][PARTITION_NONE] == 0);
    assert(fc->partition_prob[i][PARTITION_HORZ] == 0);
    fc->partition_prob[i][PARTITION_VERT] =
        av1_mode_mv_merge_probs(pre_fc->partition_prob[i][PARTITION_VERT], ct);
  }
  for (i = PARTITION_CONTEXTS_PRIMARY + PARTITION_BLOCK_SIZES;
       i < PARTITION_CONTEXTS_PRIMARY + 2 * PARTITION_BLOCK_SIZES; ++i) {
    unsigned int ct[2] = { counts->partition[i][PARTITION_HORZ],
                           counts->partition[i][PARTITION_SPLIT] };
    assert(counts->partition[i][PARTITION_NONE] == 0);
    assert(counts->partition[i][PARTITION_VERT] == 0);
    assert(fc->partition_prob[i][PARTITION_NONE] == 0);
    assert(fc->partition_prob[i][PARTITION_VERT] == 0);
    fc->partition_prob[i][PARTITION_HORZ] =
        av1_mode_mv_merge_probs(pre_fc->partition_prob[i][PARTITION_HORZ], ct);
  }
#endif
#if CONFIG_DELTA_Q
  for (i = 0; i < DELTA_Q_CONTEXTS; ++i)
    fc->delta_q_prob[i] =
        mode_mv_merge_probs(pre_fc->delta_q_prob[i], counts->delta_q[i]);
#endif
#if CONFIG_EXT_INTRA
#if CONFIG_INTRA_INTERP
  for (i = 0; i < INTRA_FILTERS + 1; ++i) {
    aom_tree_merge_probs(av1_intra_filter_tree, pre_fc->intra_filter_probs[i],
                         counts->intra_filter[i], fc->intra_filter_probs[i]);
  }
#endif  // CONFIG_INTRA_INTERP
#endif  // CONFIG_EXT_INTRA
#if CONFIG_FILTER_INTRA
  for (i = 0; i < PLANE_TYPES; ++i) {
    fc->filter_intra_probs[i] = av1_mode_mv_merge_probs(
        pre_fc->filter_intra_probs[i], counts->filter_intra[i]);
  }
#endif  // CONFIG_FILTER_INTRA
}

static void set_default_lf_deltas(struct loopfilter *lf) {
  lf->mode_ref_delta_enabled = 1;
  lf->mode_ref_delta_update = 1;

  lf->ref_deltas[INTRA_FRAME] = 1;
  lf->ref_deltas[LAST_FRAME] = 0;
#if CONFIG_EXT_REFS
  lf->ref_deltas[LAST2_FRAME] = lf->ref_deltas[LAST_FRAME];
  lf->ref_deltas[LAST3_FRAME] = lf->ref_deltas[LAST_FRAME];
  lf->ref_deltas[BWDREF_FRAME] = lf->ref_deltas[LAST_FRAME];
#endif  // CONFIG_EXT_REFS
  lf->ref_deltas[GOLDEN_FRAME] = -1;
  lf->ref_deltas[ALTREF_FRAME] = -1;

  lf->mode_deltas[0] = 0;
  lf->mode_deltas[1] = 0;
}

void av1_setup_past_independence(AV1_COMMON *cm) {
  // Reset the segment feature data to the default stats:
  // Features disabled, 0, with delta coding (Default state).
  struct loopfilter *const lf = &cm->lf;

  int i;
  av1_clearall_segfeatures(&cm->seg);
  cm->seg.abs_delta = SEGMENT_DELTADATA;

  if (cm->last_frame_seg_map && !cm->frame_parallel_decode)
    memset(cm->last_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  if (cm->current_frame_seg_map)
    memset(cm->current_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  // Reset the mode ref deltas for loop filter
  av1_zero(lf->last_ref_deltas);
  av1_zero(lf->last_mode_deltas);
  set_default_lf_deltas(lf);

  // To force update of the sharpness
  lf->last_sharpness_level = -1;

  av1_default_coef_probs(cm);
  init_mode_probs(cm->fc);
  av1_init_mv_probs(cm);
#if CONFIG_ADAPT_SCAN
  av1_init_scan_order(cm);
#endif
  av1_convolve_init();
  cm->fc->initialized = 1;

  if (cm->frame_type == KEY_FRAME || cm->error_resilient_mode ||
      cm->reset_frame_context == RESET_FRAME_CONTEXT_ALL) {
    // Reset all frame contexts.
    for (i = 0; i < FRAME_CONTEXTS; ++i) cm->frame_contexts[i] = *cm->fc;
  } else if (cm->reset_frame_context == RESET_FRAME_CONTEXT_CURRENT) {
    // Reset only the frame context specified in the frame header.
    cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
  }

  // prev_mip will only be allocated in encoder.
  if (frame_is_intra_only(cm) && cm->prev_mip && !cm->frame_parallel_decode)
    memset(cm->prev_mip, 0,
           cm->mi_stride * (cm->mi_rows + 1) * sizeof(*cm->prev_mip));

  cm->frame_context_idx = 0;
}
