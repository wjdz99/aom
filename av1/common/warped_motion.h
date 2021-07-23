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

#ifndef AOM_AV1_COMMON_WARPED_MOTION_H_
#define AOM_AV1_COMMON_WARPED_MOTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "config/aom_config.h"

#include "aom_ports/mem.h"
#include "aom_dsp/aom_dsp_common.h"
#include "av1/common/mv.h"
#include "av1/common/convolve.h"

#define MAX_PARAMDIM 9
#define LEAST_SQUARES_SAMPLES_MAX_BITS 3
#define LEAST_SQUARES_SAMPLES_MAX (1 << LEAST_SQUARES_SAMPLES_MAX_BITS)
#define SAMPLES_ARRAY_SIZE (LEAST_SQUARES_SAMPLES_MAX * 2)
#define WARPED_MOTION_DEBUG 0
#define DEFAULT_WMTYPE AFFINE
#define WARP_ERROR_BLOCK_LOG 5
#define WARP_ERROR_BLOCK (1 << WARP_ERROR_BLOCK_LOG)

extern const int16_t av1_warped_filter[WARPEDPIXEL_PREC_SHIFTS * 3 + 1][8];

DECLARE_ALIGNED(8, extern const int8_t,
                av1_filter_8bit[WARPEDPIXEL_PREC_SHIFTS * 3 + 1][8]);

/* clang-format off */
static const int error_measure_lut[512] = {
    // pow 0.7
    16384, 16339, 16294, 16249, 16204, 16158, 16113, 16068,
    16022, 15977, 15932, 15886, 15840, 15795, 15749, 15703,
    15657, 15612, 15566, 15520, 15474, 15427, 15381, 15335,
    15289, 15242, 15196, 15149, 15103, 15056, 15010, 14963,
    14916, 14869, 14822, 14775, 14728, 14681, 14634, 14587,
    14539, 14492, 14445, 14397, 14350, 14302, 14254, 14206,
    14159, 14111, 14063, 14015, 13967, 13918, 13870, 13822,
    13773, 13725, 13676, 13628, 13579, 13530, 13481, 13432,
    13383, 13334, 13285, 13236, 13187, 13137, 13088, 13038,
    12988, 12939, 12889, 12839, 12789, 12739, 12689, 12639,
    12588, 12538, 12487, 12437, 12386, 12335, 12285, 12234,
    12183, 12132, 12080, 12029, 11978, 11926, 11875, 11823,
    11771, 11719, 11667, 11615, 11563, 11511, 11458, 11406,
    11353, 11301, 11248, 11195, 11142, 11089, 11036, 10982,
    10929, 10875, 10822, 10768, 10714, 10660, 10606, 10552,
    10497, 10443, 10388, 10333, 10279, 10224, 10168, 10113,
    10058, 10002,  9947,  9891,  9835,  9779,  9723,  9666,
    9610, 9553, 9497, 9440, 9383, 9326, 9268, 9211,
    9153, 9095, 9037, 8979, 8921, 8862, 8804, 8745,
    8686, 8627, 8568, 8508, 8449, 8389, 8329, 8269,
    8208, 8148, 8087, 8026, 7965, 7903, 7842, 7780,
    7718, 7656, 7593, 7531, 7468, 7405, 7341, 7278,
    7214, 7150, 7086, 7021, 6956, 6891, 6826, 6760,
    6695, 6628, 6562, 6495, 6428, 6361, 6293, 6225,
    6157, 6089, 6020, 5950, 5881, 5811, 5741, 5670,
    5599, 5527, 5456, 5383, 5311, 5237, 5164, 5090,
    5015, 4941, 4865, 4789, 4713, 4636, 4558, 4480,
    4401, 4322, 4242, 4162, 4080, 3998, 3916, 3832,
    3748, 3663, 3577, 3490, 3402, 3314, 3224, 3133,
    3041, 2948, 2854, 2758, 2661, 2562, 2461, 2359,
    2255, 2148, 2040, 1929, 1815, 1698, 1577, 1452,
    1323, 1187, 1045,  894,  731,  550,  339,    0,
    339,  550,  731,  894, 1045, 1187, 1323, 1452,
    1577, 1698, 1815, 1929, 2040, 2148, 2255, 2359,
    2461, 2562, 2661, 2758, 2854, 2948, 3041, 3133,
    3224, 3314, 3402, 3490, 3577, 3663, 3748, 3832,
    3916, 3998, 4080, 4162, 4242, 4322, 4401, 4480,
    4558, 4636, 4713, 4789, 4865, 4941, 5015, 5090,
    5164, 5237, 5311, 5383, 5456, 5527, 5599, 5670,
    5741, 5811, 5881, 5950, 6020, 6089, 6157, 6225,
    6293, 6361, 6428, 6495, 6562, 6628, 6695, 6760,
    6826, 6891, 6956, 7021, 7086, 7150, 7214, 7278,
    7341, 7405, 7468, 7531, 7593, 7656, 7718, 7780,
    7842, 7903, 7965, 8026, 8087, 8148, 8208, 8269,
    8329, 8389, 8449, 8508, 8568, 8627, 8686, 8745,
    8804, 8862, 8921, 8979, 9037, 9095, 9153, 9211,
    9268, 9326, 9383, 9440, 9497, 9553, 9610, 9666,
    9723,  9779,  9835,  9891,  9947, 10002, 10058, 10113,
    10168, 10224, 10279, 10333, 10388, 10443, 10497, 10552,
    10606, 10660, 10714, 10768, 10822, 10875, 10929, 10982,
    11036, 11089, 11142, 11195, 11248, 11301, 11353, 11406,
    11458, 11511, 11563, 11615, 11667, 11719, 11771, 11823,
    11875, 11926, 11978, 12029, 12080, 12132, 12183, 12234,
    12285, 12335, 12386, 12437, 12487, 12538, 12588, 12639,
    12689, 12739, 12789, 12839, 12889, 12939, 12988, 13038,
    13088, 13137, 13187, 13236, 13285, 13334, 13383, 13432,
    13481, 13530, 13579, 13628, 13676, 13725, 13773, 13822,
    13870, 13918, 13967, 14015, 14063, 14111, 14159, 14206,
    14254, 14302, 14350, 14397, 14445, 14492, 14539, 14587,
    14634, 14681, 14728, 14775, 14822, 14869, 14916, 14963,
    15010, 15056, 15103, 15149, 15196, 15242, 15289, 15335,
    15381, 15427, 15474, 15520, 15566, 15612, 15657, 15703,
    15749, 15795, 15840, 15886, 15932, 15977, 16022, 16068,
    16113, 16158, 16204, 16249, 16294, 16339, 16384, 16384,
};
/* clang-format on */

static const uint8_t warp_pad_left[14][16] = {
  { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 4, 4, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 5, 5, 5, 5, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 6, 6, 6, 6, 6, 6, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 7, 7, 7, 7, 7, 7, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 10, 11, 12, 13, 14, 15 },
  { 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 11, 12, 13, 14, 15 },
  { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 12, 13, 14, 15 },
  { 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 13, 14, 15 },
  { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 14, 15 },
  { 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 15 },
  { 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15 },
};

static const uint8_t warp_pad_right[14][16] = {
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 11, 11 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 10 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9, 9, 9, 9 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7 },
  { 0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
  { 0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 },
  { 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 },
  { 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 },
  { 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
  { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
};

static INLINE int error_measure(int err) {
  return error_measure_lut[255 + err];
}

#define DIV_LUT_PREC_BITS 14
#define DIV_LUT_BITS 8
#define DIV_LUT_NUM (1 << DIV_LUT_BITS)

static const uint16_t div_lut[DIV_LUT_NUM + 1] = {
  16384, 16320, 16257, 16194, 16132, 16070, 16009, 15948, 15888, 15828, 15768,
  15709, 15650, 15592, 15534, 15477, 15420, 15364, 15308, 15252, 15197, 15142,
  15087, 15033, 14980, 14926, 14873, 14821, 14769, 14717, 14665, 14614, 14564,
  14513, 14463, 14413, 14364, 14315, 14266, 14218, 14170, 14122, 14075, 14028,
  13981, 13935, 13888, 13843, 13797, 13752, 13707, 13662, 13618, 13574, 13530,
  13487, 13443, 13400, 13358, 13315, 13273, 13231, 13190, 13148, 13107, 13066,
  13026, 12985, 12945, 12906, 12866, 12827, 12788, 12749, 12710, 12672, 12633,
  12596, 12558, 12520, 12483, 12446, 12409, 12373, 12336, 12300, 12264, 12228,
  12193, 12157, 12122, 12087, 12053, 12018, 11984, 11950, 11916, 11882, 11848,
  11815, 11782, 11749, 11716, 11683, 11651, 11619, 11586, 11555, 11523, 11491,
  11460, 11429, 11398, 11367, 11336, 11305, 11275, 11245, 11215, 11185, 11155,
  11125, 11096, 11067, 11038, 11009, 10980, 10951, 10923, 10894, 10866, 10838,
  10810, 10782, 10755, 10727, 10700, 10673, 10645, 10618, 10592, 10565, 10538,
  10512, 10486, 10460, 10434, 10408, 10382, 10356, 10331, 10305, 10280, 10255,
  10230, 10205, 10180, 10156, 10131, 10107, 10082, 10058, 10034, 10010, 9986,
  9963,  9939,  9916,  9892,  9869,  9846,  9823,  9800,  9777,  9754,  9732,
  9709,  9687,  9664,  9642,  9620,  9598,  9576,  9554,  9533,  9511,  9489,
  9468,  9447,  9425,  9404,  9383,  9362,  9341,  9321,  9300,  9279,  9259,
  9239,  9218,  9198,  9178,  9158,  9138,  9118,  9098,  9079,  9059,  9039,
  9020,  9001,  8981,  8962,  8943,  8924,  8905,  8886,  8867,  8849,  8830,
  8812,  8793,  8775,  8756,  8738,  8720,  8702,  8684,  8666,  8648,  8630,
  8613,  8595,  8577,  8560,  8542,  8525,  8508,  8490,  8473,  8456,  8439,
  8422,  8405,  8389,  8372,  8355,  8339,  8322,  8306,  8289,  8273,  8257,
  8240,  8224,  8208,  8192,
};

// Decomposes a divisor D such that 1/D = y/2^shift, where y is returned
// at precision of DIV_LUT_PREC_BITS along with the shift.
static int16_t resolve_divisor_64(uint64_t D, int16_t *shift) {
  int64_t f;
  *shift = (int16_t)((D >> 32) ? get_msb((unsigned int)(D >> 32)) + 32
                               : get_msb((unsigned int)D));
  // e is obtained from D after resetting the most significant 1 bit.
  const int64_t e = D - ((uint64_t)1 << *shift);
  // Get the most significant DIV_LUT_BITS (8) bits of e into f
  if (*shift > DIV_LUT_BITS)
    f = ROUND_POWER_OF_TWO_64(e, *shift - DIV_LUT_BITS);
  else
    f = e << (DIV_LUT_BITS - *shift);
  assert(f <= DIV_LUT_NUM);
  *shift += DIV_LUT_PREC_BITS;
  // Use f as lookup into the precomputed table of multipliers
  return div_lut[f];
}

static int16_t resolve_divisor_32(uint32_t D, int16_t *shift) {
  int32_t f;
  *shift = get_msb(D);
  // e is obtained from D after resetting the most significant 1 bit.
  const int32_t e = D - ((uint32_t)1 << *shift);
  // Get the most significant DIV_LUT_BITS (8) bits of e into f
  if (*shift > DIV_LUT_BITS)
    f = ROUND_POWER_OF_TWO(e, *shift - DIV_LUT_BITS);
  else
    f = e << (DIV_LUT_BITS - *shift);
  assert(f <= DIV_LUT_NUM);
  *shift += DIV_LUT_PREC_BITS;
  // Use f as lookup into the precomputed table of multipliers
  return div_lut[f];
}

// Returns the error between the frame described by 'ref' and the frame
// described by 'dst'.
int64_t av1_frame_error(int use_hbd, int bd, const uint8_t *ref, int stride,
                        uint8_t *dst, int p_width, int p_height, int p_stride);

int64_t av1_segmented_frame_error(int use_hbd, int bd, const uint8_t *ref,
                                  int stride, uint8_t *dst, int p_width,
                                  int p_height, int p_stride,
                                  uint8_t *segment_map, int segment_map_stride);

int64_t av1_calc_highbd_frame_error(const uint16_t *const ref, int stride,
                                    const uint16_t *const dst, int p_width,
                                    int p_height, int p_stride, int bd);

void highbd_warp_plane(WarpedMotionParams *wm, const uint16_t *const ref,
                       int width, int height, int stride, uint16_t *const pred,
                       int p_col, int p_row, int p_width, int p_height,
                       int p_stride, int subsampling_x, int subsampling_y,
                       int bd, ConvolveParams *conv_params);

void warp_plane(WarpedMotionParams *wm, const uint8_t *const ref, int width,
                int height, int stride, uint8_t *pred, int p_col, int p_row,
                int p_width, int p_height, int p_stride, int subsampling_x,
                int subsampling_y, ConvolveParams *conv_params);

void av1_warp_plane(WarpedMotionParams *wm, int use_hbd, int bd,
                    const uint8_t *ref, int width, int height, int stride,
                    uint8_t *pred, int p_col, int p_row, int p_width,
                    int p_height, int p_stride, int subsampling_x,
                    int subsampling_y, ConvolveParams *conv_params);

int av1_find_projection(int np, const int *pts1, const int *pts2,
                        BLOCK_SIZE bsize, int mvy, int mvx,
                        WarpedMotionParams *wm_params, int mi_row, int mi_col);

int av1_get_shear_params(WarpedMotionParams *wm);
#endif  // AOM_AV1_COMMON_WARPED_MOTION_H_
