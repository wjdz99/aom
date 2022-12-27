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

#ifndef AOM_AV1_ENCODER_NONRD_OPT_H_
#define AOM_AV1_ENCODER_NONRD_OPT_H_

#include "av1/encoder/rdopt_utils.h"

#define RTC_INTER_MODES (4)
#define RTC_INTRA_MODES (4)
#define RTC_MODES (AOMMAX(RTC_INTER_MODES, RTC_INTRA_MODES))
#define CALC_BIASED_RDCOST(rdcost) (7 * (rdcost) >> 3)
#define NUM_COMP_INTER_MODES_RT (6)
#define NUM_INTER_MODES 12
#define CAP_TX_SIZE_FOR_BSIZE_GT32(tx_mode_search_type, bsize) \
  (((tx_mode_search_type) != ONLY_4X4 && (bsize) > BLOCK_32X32) ? true : false)
#define TX_SIZE_FOR_BSIZE_GT32 (TX_16X16)
#define FILTER_SEARCH_SIZE 2
#if !CONFIG_REALTIME_ONLY
#define MOTION_MODE_SEARCH_SIZE 2
#endif

extern int g_pick_inter_mode_cnt;
/*!\cond */
typedef struct {
  uint8_t *data;
  int stride;
  int in_use;
} PRED_BUFFER;

typedef struct {
  PRED_BUFFER *best_pred;
  PREDICTION_MODE best_mode;
  TX_SIZE best_tx_size;
  TX_TYPE tx_type;
  MV_REFERENCE_FRAME best_ref_frame;
  MV_REFERENCE_FRAME best_second_ref_frame;
  uint8_t best_mode_skip_txfm;
  uint8_t best_mode_initial_skip_flag;
  int_interpfilters best_pred_filter;
  MOTION_MODE best_motion_mode;
  WarpedMotionParams wm_params;
  int num_proj_ref;
  uint8_t blk_skip[MAX_MIB_SIZE * MAX_MIB_SIZE / 4];
  PALETTE_MODE_INFO pmi;
  int64_t best_sse;
} BEST_PICKMODE;

typedef struct {
  MV_REFERENCE_FRAME ref_frame;
  PREDICTION_MODE pred_mode;
} REF_MODE;

typedef struct {
  MV_REFERENCE_FRAME ref_frame[2];
  PREDICTION_MODE pred_mode;
} COMP_REF_MODE;

struct estimate_block_intra_args {
  AV1_COMP *cpi;
  MACROBLOCK *x;
  PREDICTION_MODE mode;
  int skippable;
  RD_STATS *rdc;
};
/*!\endcond */

/*!\brief Structure to store parameters and statistics used in non-rd inter mode
 * evaluation.
 */
typedef struct {
  //! Structure to hold best inter mode data
  BEST_PICKMODE best_pickmode;
  //! Structure to RD cost of current mode
  RD_STATS this_rdc;
  //! Pointer to the RD Cost for the best mode found so far
  RD_STATS best_rdc;
  //! Distortion of chroma planes for all modes and reference frames
  int64_t uv_dist[RTC_INTER_MODES][REF_FRAMES];
  //! Buffer to hold predicted block for all reference frames and planes
  struct buf_2d yv12_mb[REF_FRAMES][MAX_MB_PLANE];
  //! Array to hold variance of all modes and reference frames
  unsigned int vars[RTC_INTER_MODES][REF_FRAMES];
  //! Array to hold ref cost of single reference mode for all ref frames
  unsigned int ref_costs_single[REF_FRAMES];
  //! Array to hold motion vector for all modes and reference frames
  int_mv frame_mv[MB_MODE_COUNT][REF_FRAMES];
  //! Array to hold best mv for all modes and reference frames
  int_mv frame_mv_best[MB_MODE_COUNT][REF_FRAMES];
  //! Array to hold inter mode cost of single ref mode for all ref frames
  int single_inter_mode_costs[RTC_INTER_MODES][REF_FRAMES];
  //! Array to hold use reference frame mask for each reference frame
  int use_ref_frame_mask[REF_FRAMES];
  //! Array to hold flags of evaluated modes for each reference frame
  uint8_t mode_checked[MB_MODE_COUNT][REF_FRAMES];
} InterModeSearchStateNonrd;

static const uint8_t b_width_log2_lookup[BLOCK_SIZES] = { 0, 0, 1, 1, 1, 2,
                                                          2, 2, 3, 3, 3, 4,
                                                          4, 4, 5, 5 };
static const uint8_t b_height_log2_lookup[BLOCK_SIZES] = { 0, 1, 0, 1, 2, 1,
                                                           2, 3, 2, 3, 4, 3,
                                                           4, 5, 4, 5 };

static const PREDICTION_MODE intra_mode_list[] = { DC_PRED, V_PRED, H_PRED,
                                                   SMOOTH_PRED };

static const PREDICTION_MODE inter_mode_list[] = { NEARESTMV, NEARMV, GLOBALMV,
                                                   NEWMV };

static const THR_MODES mode_idx[REF_FRAMES][RTC_MODES] = {
  { THR_DC, THR_V_PRED, THR_H_PRED, THR_SMOOTH },
  { THR_NEARESTMV, THR_NEARMV, THR_GLOBALMV, THR_NEWMV },
  { THR_NEARESTL2, THR_NEARL2, THR_GLOBALL2, THR_NEWL2 },
  { THR_NEARESTL3, THR_NEARL3, THR_GLOBALL3, THR_NEWL3 },
  { THR_NEARESTG, THR_NEARG, THR_GLOBALG, THR_NEWG },
  { THR_NEARESTB, THR_NEARB, THR_GLOBALB, THR_NEWB },
  { THR_NEARESTA2, THR_NEARA2, THR_GLOBALA2, THR_NEWA2 },
  { THR_NEARESTA, THR_NEARA, THR_GLOBALA, THR_NEWA },
};

// GLOBALMV in the set below is in fact ZEROMV as we don't do global ME in RT
// mode
static const REF_MODE ref_mode_set[NUM_INTER_MODES] = {
  { LAST_FRAME, NEARESTMV },   { LAST_FRAME, NEARMV },
  { LAST_FRAME, GLOBALMV },    { LAST_FRAME, NEWMV },
  { GOLDEN_FRAME, NEARESTMV }, { GOLDEN_FRAME, NEARMV },
  { GOLDEN_FRAME, GLOBALMV },  { GOLDEN_FRAME, NEWMV },
  { ALTREF_FRAME, NEARESTMV }, { ALTREF_FRAME, NEARMV },
  { ALTREF_FRAME, GLOBALMV },  { ALTREF_FRAME, NEWMV },
};

static const COMP_REF_MODE comp_ref_mode_set[NUM_COMP_INTER_MODES_RT] = {
  { { LAST_FRAME, GOLDEN_FRAME }, GLOBAL_GLOBALMV },
  { { LAST_FRAME, GOLDEN_FRAME }, NEAREST_NEARESTMV },
  { { LAST_FRAME, LAST2_FRAME }, GLOBAL_GLOBALMV },
  { { LAST_FRAME, LAST2_FRAME }, NEAREST_NEARESTMV },
  { { LAST_FRAME, ALTREF_FRAME }, GLOBAL_GLOBALMV },
  { { LAST_FRAME, ALTREF_FRAME }, NEAREST_NEARESTMV },
};

static const int_interpfilters filters_ref_set[9] = {
  [0].as_filters = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR },
  [1].as_filters = { EIGHTTAP_SMOOTH, EIGHTTAP_SMOOTH },
  [2].as_filters = { EIGHTTAP_REGULAR, EIGHTTAP_SMOOTH },
  [3].as_filters = { EIGHTTAP_SMOOTH, EIGHTTAP_REGULAR },
  [4].as_filters = { MULTITAP_SHARP, MULTITAP_SHARP },
  [5].as_filters = { EIGHTTAP_REGULAR, MULTITAP_SHARP },
  [6].as_filters = { MULTITAP_SHARP, EIGHTTAP_REGULAR },
  [7].as_filters = { EIGHTTAP_SMOOTH, MULTITAP_SHARP },
  [8].as_filters = { MULTITAP_SHARP, EIGHTTAP_SMOOTH }
};

enum {
  //  INTER_ALL = (1 << NEARESTMV) | (1 << NEARMV) | (1 << NEWMV),
  INTER_NEAREST = (1 << NEARESTMV),
  INTER_NEAREST_NEW = (1 << NEARESTMV) | (1 << NEWMV),
  INTER_NEAREST_NEAR = (1 << NEARESTMV) | (1 << NEARMV),
  INTER_NEAR_NEW = (1 << NEARMV) | (1 << NEWMV),
};

// The original scan order (default_scan_8x8) is modified according to the extra
// transpose in hadamard c implementation, i.e., aom_hadamard_lp_8x8_c and
// aom_hadamard_8x8_c.
static const int16_t default_scan_8x8_transpose[64] = {
  0,  8,  1,  2,  9,  16, 24, 17, 10, 3,  4,  11, 18, 25, 32, 40,
  33, 26, 19, 12, 5,  6,  13, 20, 27, 34, 41, 48, 56, 49, 42, 35,
  28, 21, 14, 7,  15, 22, 29, 36, 43, 50, 57, 58, 51, 44, 37, 30,
  23, 31, 38, 45, 52, 59, 60, 53, 46, 39, 47, 54, 61, 62, 55, 63
};

// The original scan order (av1_default_iscan_8x8) is modified to match
// hadamard AVX2 implementation, i.e., aom_hadamard_lp_8x8_avx2 and
// aom_hadamard_8x8_avx2. Since hadamard AVX2 implementation will modify the
// order of coefficients, such that the normal scan order is no longer
// guaranteed to scan low coefficients first, therefore we modify the scan order
// accordingly.
// Note that this one has to be used together with default_scan_8x8_transpose.
static const int16_t av1_default_iscan_8x8_transpose[64] = {
  0,  2,  3,  9,  10, 20, 21, 35, 1,  4,  8,  11, 19, 22, 34, 36,
  5,  7,  12, 18, 23, 33, 37, 48, 6,  13, 17, 24, 32, 38, 47, 49,
  14, 16, 25, 31, 39, 46, 50, 57, 15, 26, 30, 40, 45, 51, 56, 58,
  27, 29, 41, 44, 52, 55, 59, 62, 28, 42, 43, 53, 54, 60, 61, 63
};

// The original scan order (default_scan_16x16) is modified according to the
// extra transpose in hadamard c implementation in lp case, i.e.,
// aom_hadamard_lp_16x16_c.
static const int16_t default_scan_lp_16x16_transpose[256] = {
  0,   8,   2,   4,   10,  16,  24,  18,  12,  6,   64,  14,  20,  26,  32,
  40,  34,  28,  22,  72,  66,  68,  74,  80,  30,  36,  42,  48,  56,  50,
  44,  38,  88,  82,  76,  70,  128, 78,  84,  90,  96,  46,  52,  58,  1,
  9,   3,   60,  54,  104, 98,  92,  86,  136, 130, 132, 138, 144, 94,  100,
  106, 112, 62,  5,   11,  17,  25,  19,  13,  7,   120, 114, 108, 102, 152,
  146, 140, 134, 192, 142, 148, 154, 160, 110, 116, 122, 65,  15,  21,  27,
  33,  41,  35,  29,  23,  73,  67,  124, 118, 168, 162, 156, 150, 200, 194,
  196, 202, 208, 158, 164, 170, 176, 126, 69,  75,  81,  31,  37,  43,  49,
  57,  51,  45,  39,  89,  83,  77,  71,  184, 178, 172, 166, 216, 210, 204,
  198, 206, 212, 218, 224, 174, 180, 186, 129, 79,  85,  91,  97,  47,  53,
  59,  61,  55,  105, 99,  93,  87,  137, 131, 188, 182, 232, 226, 220, 214,
  222, 228, 234, 240, 190, 133, 139, 145, 95,  101, 107, 113, 63,  121, 115,
  109, 103, 153, 147, 141, 135, 248, 242, 236, 230, 238, 244, 250, 193, 143,
  149, 155, 161, 111, 117, 123, 125, 119, 169, 163, 157, 151, 201, 195, 252,
  246, 254, 197, 203, 209, 159, 165, 171, 177, 127, 185, 179, 173, 167, 217,
  211, 205, 199, 207, 213, 219, 225, 175, 181, 187, 189, 183, 233, 227, 221,
  215, 223, 229, 235, 241, 191, 249, 243, 237, 231, 239, 245, 251, 253, 247,
  255
};

#if CONFIG_AV1_HIGHBITDEPTH
// The original scan order (default_scan_16x16) is modified according to the
// extra shift in hadamard c implementation in fp case, i.e.,
// aom_hadamard_16x16_c. Note that 16x16 lp and fp hadamard generate different
// outputs, so we handle them separately.
static const int16_t default_scan_fp_16x16_transpose[256] = {
  0,   4,   2,   8,   6,   16,  20,  18,  12,  10,  64,  14,  24,  22,  32,
  36,  34,  28,  26,  68,  66,  72,  70,  80,  30,  40,  38,  48,  52,  50,
  44,  42,  84,  82,  76,  74,  128, 78,  88,  86,  96,  46,  56,  54,  1,
  5,   3,   60,  58,  100, 98,  92,  90,  132, 130, 136, 134, 144, 94,  104,
  102, 112, 62,  9,   7,   17,  21,  19,  13,  11,  116, 114, 108, 106, 148,
  146, 140, 138, 192, 142, 152, 150, 160, 110, 120, 118, 65,  15,  25,  23,
  33,  37,  35,  29,  27,  69,  67,  124, 122, 164, 162, 156, 154, 196, 194,
  200, 198, 208, 158, 168, 166, 176, 126, 73,  71,  81,  31,  41,  39,  49,
  53,  51,  45,  43,  85,  83,  77,  75,  180, 178, 172, 170, 212, 210, 204,
  202, 206, 216, 214, 224, 174, 184, 182, 129, 79,  89,  87,  97,  47,  57,
  55,  61,  59,  101, 99,  93,  91,  133, 131, 188, 186, 228, 226, 220, 218,
  222, 232, 230, 240, 190, 137, 135, 145, 95,  105, 103, 113, 63,  117, 115,
  109, 107, 149, 147, 141, 139, 244, 242, 236, 234, 238, 248, 246, 193, 143,
  153, 151, 161, 111, 121, 119, 125, 123, 165, 163, 157, 155, 197, 195, 252,
  250, 254, 201, 199, 209, 159, 169, 167, 177, 127, 181, 179, 173, 171, 213,
  211, 205, 203, 207, 217, 215, 225, 175, 185, 183, 189, 187, 229, 227, 221,
  219, 223, 233, 231, 241, 191, 245, 243, 237, 235, 239, 249, 247, 253, 251,
  255
};
#endif

// The original scan order (av1_default_iscan_16x16) is modified to match
// hadamard AVX2 implementation, i.e., aom_hadamard_lp_16x16_avx2.
// Since hadamard AVX2 implementation will modify the order of coefficients,
// such that the normal scan order is no longer guaranteed to scan low
// coefficients first, therefore we modify the scan order accordingly. Note that
// this one has to be used together with default_scan_lp_16x16_transpose.
static const int16_t av1_default_iscan_lp_16x16_transpose[256] = {
  0,   44,  2,   46,  3,   63,  9,   69,  1,   45,  4,   64,  8,   68,  11,
  87,  5,   65,  7,   67,  12,  88,  18,  94,  6,   66,  13,  89,  17,  93,
  24,  116, 14,  90,  16,  92,  25,  117, 31,  123, 15,  91,  26,  118, 30,
  122, 41,  148, 27,  119, 29,  121, 42,  149, 48,  152, 28,  120, 43,  150,
  47,  151, 62,  177, 10,  86,  20,  96,  21,  113, 35,  127, 19,  95,  22,
  114, 34,  126, 37,  144, 23,  115, 33,  125, 38,  145, 52,  156, 32,  124,
  39,  146, 51,  155, 58,  173, 40,  147, 50,  154, 59,  174, 73,  181, 49,
  153, 60,  175, 72,  180, 83,  198, 61,  176, 71,  179, 84,  199, 98,  202,
  70,  178, 85,  200, 97,  201, 112, 219, 36,  143, 54,  158, 55,  170, 77,
  185, 53,  157, 56,  171, 76,  184, 79,  194, 57,  172, 75,  183, 80,  195,
  102, 206, 74,  182, 81,  196, 101, 205, 108, 215, 82,  197, 100, 204, 109,
  216, 131, 223, 99,  203, 110, 217, 130, 222, 140, 232, 111, 218, 129, 221,
  141, 233, 160, 236, 128, 220, 142, 234, 159, 235, 169, 245, 78,  193, 104,
  208, 105, 212, 135, 227, 103, 207, 106, 213, 134, 226, 136, 228, 107, 214,
  133, 225, 137, 229, 164, 240, 132, 224, 138, 230, 163, 239, 165, 241, 139,
  231, 162, 238, 166, 242, 189, 249, 161, 237, 167, 243, 188, 248, 190, 250,
  168, 244, 187, 247, 191, 251, 210, 254, 186, 246, 192, 252, 209, 253, 211,
  255
};

#if CONFIG_AV1_HIGHBITDEPTH
// The original scan order (av1_default_iscan_16x16) is modified to match
// hadamard AVX2 implementation, i.e., aom_hadamard_16x16_avx2.
// Since hadamard AVX2 implementation will modify the order of coefficients,
// such that the normal scan order is no longer guaranteed to scan low
// coefficients first, therefore we modify the scan order accordingly. Note that
// this one has to be used together with default_scan_fp_16x16_transpose.
static const int16_t av1_default_iscan_fp_16x16_transpose[256] = {
  0,   44,  2,   46,  1,   45,  4,   64,  3,   63,  9,   69,  8,   68,  11,
  87,  5,   65,  7,   67,  6,   66,  13,  89,  12,  88,  18,  94,  17,  93,
  24,  116, 14,  90,  16,  92,  15,  91,  26,  118, 25,  117, 31,  123, 30,
  122, 41,  148, 27,  119, 29,  121, 28,  120, 43,  150, 42,  149, 48,  152,
  47,  151, 62,  177, 10,  86,  20,  96,  19,  95,  22,  114, 21,  113, 35,
  127, 34,  126, 37,  144, 23,  115, 33,  125, 32,  124, 39,  146, 38,  145,
  52,  156, 51,  155, 58,  173, 40,  147, 50,  154, 49,  153, 60,  175, 59,
  174, 73,  181, 72,  180, 83,  198, 61,  176, 71,  179, 70,  178, 85,  200,
  84,  199, 98,  202, 97,  201, 112, 219, 36,  143, 54,  158, 53,  157, 56,
  171, 55,  170, 77,  185, 76,  184, 79,  194, 57,  172, 75,  183, 74,  182,
  81,  196, 80,  195, 102, 206, 101, 205, 108, 215, 82,  197, 100, 204, 99,
  203, 110, 217, 109, 216, 131, 223, 130, 222, 140, 232, 111, 218, 129, 221,
  128, 220, 142, 234, 141, 233, 160, 236, 159, 235, 169, 245, 78,  193, 104,
  208, 103, 207, 106, 213, 105, 212, 135, 227, 134, 226, 136, 228, 107, 214,
  133, 225, 132, 224, 138, 230, 137, 229, 164, 240, 163, 239, 165, 241, 139,
  231, 162, 238, 161, 237, 167, 243, 166, 242, 189, 249, 188, 248, 190, 250,
  168, 244, 187, 247, 186, 246, 192, 252, 191, 251, 210, 254, 209, 253, 211,
  255
};
#endif

// For entropy coding, IDTX shares the scan orders of the other 2D-transforms,
// but the fastest way to calculate the IDTX transform (i.e. no transposes)
// results in coefficients that are a transposition of the entropy coding
// versions. These tables are used as substitute for the scan order for the
// faster version of IDTX.

// Must be used together with av1_fast_idtx_iscan_4x4
static const int16_t av1_fast_idtx_scan_4x4[16] = {
  0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15
};

// Must be used together with av1_fast_idtx_scan_4x4
static const int16_t av1_fast_idtx_iscan_4x4[16] = { 0, 1,  5,  6, 2,  4,
                                                     7, 12, 3,  8, 11, 13,
                                                     9, 10, 14, 15 };

static const SCAN_ORDER av1_fast_idtx_scan_order_4x4 = {
  av1_fast_idtx_scan_4x4, av1_fast_idtx_iscan_4x4
};

// Must be used together with av1_fast_idtx_iscan_8x8
static const int16_t av1_fast_idtx_scan_8x8[64] = {
  0,  1,  8,  16, 9,  2,  3,  10, 17, 24, 32, 25, 18, 11, 4,  5,
  12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6,  7,  14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};

// Must be used together with av1_fast_idtx_scan_8x8
static const int16_t av1_fast_idtx_iscan_8x8[64] = {
  0,  1,  5,  6,  14, 15, 27, 28, 2,  4,  7,  13, 16, 26, 29, 42,
  3,  8,  12, 17, 25, 30, 41, 43, 9,  11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54, 20, 22, 33, 38, 46, 51, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61, 35, 36, 48, 49, 57, 58, 62, 63
};

static const SCAN_ORDER av1_fast_idtx_scan_order_8x8 = {
  av1_fast_idtx_scan_8x8, av1_fast_idtx_iscan_8x8
};

// Must be used together with av1_fast_idtx_iscan_16x16
static const int16_t av1_fast_idtx_scan_16x16[256] = {
  0,   1,   16,  32,  17,  2,   3,   18,  33,  48,  64,  49,  34,  19,  4,
  5,   20,  35,  50,  65,  80,  96,  81,  66,  51,  36,  21,  6,   7,   22,
  37,  52,  67,  82,  97,  112, 128, 113, 98,  83,  68,  53,  38,  23,  8,
  9,   24,  39,  54,  69,  84,  99,  114, 129, 144, 160, 145, 130, 115, 100,
  85,  70,  55,  40,  25,  10,  11,  26,  41,  56,  71,  86,  101, 116, 131,
  146, 161, 176, 192, 177, 162, 147, 132, 117, 102, 87,  72,  57,  42,  27,
  12,  13,  28,  43,  58,  73,  88,  103, 118, 133, 148, 163, 178, 193, 208,
  224, 209, 194, 179, 164, 149, 134, 119, 104, 89,  74,  59,  44,  29,  14,
  15,  30,  45,  60,  75,  90,  105, 120, 135, 150, 165, 180, 195, 210, 225,
  240, 241, 226, 211, 196, 181, 166, 151, 136, 121, 106, 91,  76,  61,  46,
  31,  47,  62,  77,  92,  107, 122, 137, 152, 167, 182, 197, 212, 227, 242,
  243, 228, 213, 198, 183, 168, 153, 138, 123, 108, 93,  78,  63,  79,  94,
  109, 124, 139, 154, 169, 184, 199, 214, 229, 244, 245, 230, 215, 200, 185,
  170, 155, 140, 125, 110, 95,  111, 126, 141, 156, 171, 186, 201, 216, 231,
  246, 247, 232, 217, 202, 187, 172, 157, 142, 127, 143, 158, 173, 188, 203,
  218, 233, 248, 249, 234, 219, 204, 189, 174, 159, 175, 190, 205, 220, 235,
  250, 251, 236, 221, 206, 191, 207, 222, 237, 252, 253, 238, 223, 239, 254,
  255
};

// Must be used together with av1_fast_idtx_scan_16x16
static const int16_t av1_fast_idtx_iscan_16x16[256] = {
  0,   1,   5,   6,   14,  15,  27,  28,  44,  45,  65,  66,  90,  91,  119,
  120, 2,   4,   7,   13,  16,  26,  29,  43,  46,  64,  67,  89,  92,  118,
  121, 150, 3,   8,   12,  17,  25,  30,  42,  47,  63,  68,  88,  93,  117,
  122, 149, 151, 9,   11,  18,  24,  31,  41,  48,  62,  69,  87,  94,  116,
  123, 148, 152, 177, 10,  19,  23,  32,  40,  49,  61,  70,  86,  95,  115,
  124, 147, 153, 176, 178, 20,  22,  33,  39,  50,  60,  71,  85,  96,  114,
  125, 146, 154, 175, 179, 200, 21,  34,  38,  51,  59,  72,  84,  97,  113,
  126, 145, 155, 174, 180, 199, 201, 35,  37,  52,  58,  73,  83,  98,  112,
  127, 144, 156, 173, 181, 198, 202, 219, 36,  53,  57,  74,  82,  99,  111,
  128, 143, 157, 172, 182, 197, 203, 218, 220, 54,  56,  75,  81,  100, 110,
  129, 142, 158, 171, 183, 196, 204, 217, 221, 234, 55,  76,  80,  101, 109,
  130, 141, 159, 170, 184, 195, 205, 216, 222, 233, 235, 77,  79,  102, 108,
  131, 140, 160, 169, 185, 194, 206, 215, 223, 232, 236, 245, 78,  103, 107,
  132, 139, 161, 168, 186, 193, 207, 214, 224, 231, 237, 244, 246, 104, 106,
  133, 138, 162, 167, 187, 192, 208, 213, 225, 230, 238, 243, 247, 252, 105,
  134, 137, 163, 166, 188, 191, 209, 212, 226, 229, 239, 242, 248, 251, 253,
  135, 136, 164, 165, 189, 190, 210, 211, 227, 228, 240, 241, 249, 250, 254,
  255
};

// Indicates the blocks for which RD model should be based on special logic
static INLINE int get_model_rd_flag(const AV1_COMP *cpi, const MACROBLOCKD *xd,
                                    BLOCK_SIZE bsize) {
  const int large_block = bsize >= BLOCK_32X32;
  const AV1_COMMON *const cm = &cpi->common;
  return cpi->oxcf.rc_cfg.mode == AOM_CBR && large_block &&
         !cyclic_refresh_segment_id_boosted(xd->mi[0]->segment_id) &&
         cm->quant_params.base_qindex &&
         cm->seq_params->bit_depth == AOM_BITS_8;
}
/*!\brief Finds predicted motion vectors for a block.
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Finds predicted motion vectors for a block from a certain reference frame.
 * First, it fills reference MV stack, then picks the test from the stack and
 * predicts the final MV for a block for each mode.
 * \param[in]    cpi                      Top-level encoder structure
 * \param[in]    x                        Pointer to structure holding all the
 *                                        data for the current macroblock
 * \param[in]    ref_frame                Reference frame for which to find
 *                                        ref MVs
 * \param[in]    frame_mv                 Predicted MVs for a block
 * \param[in]    yv12_mb                  Buffer to hold predicted block
 * \param[in]    bsize                    Current block size
 * \param[in]    force_skip_low_temp_var  Flag indicating possible mode search
 *                                        prune for low temporal variance block
 * \param[in]    skip_pred_mv             Flag indicating to skip av1_mv_pred
 *
 * \remark Nothing is returned. Instead, predicted MVs are placed into
 * \c frame_mv array
 */
static INLINE void find_predictors(AV1_COMP *cpi, MACROBLOCK *x,
                                   MV_REFERENCE_FRAME ref_frame,
                                   int_mv frame_mv[MB_MODE_COUNT][REF_FRAMES],
                                   struct buf_2d yv12_mb[8][MAX_MB_PLANE],
                                   BLOCK_SIZE bsize,
                                   int force_skip_low_temp_var,
                                   int skip_pred_mv) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = xd->mi[0];
  MB_MODE_INFO_EXT *const mbmi_ext = &x->mbmi_ext;
  const YV12_BUFFER_CONFIG *yv12 = get_ref_frame_yv12_buf(cm, ref_frame);
  const int num_planes = av1_num_planes(cm);

  x->pred_mv_sad[ref_frame] = INT_MAX;
  x->pred_mv0_sad[ref_frame] = INT_MAX;
  x->pred_mv1_sad[ref_frame] = INT_MAX;
  frame_mv[NEWMV][ref_frame].as_int = INVALID_MV;
  // TODO(kyslov) this needs various further optimizations. to be continued..
  assert(yv12 != NULL);
  if (yv12 != NULL) {
    const struct scale_factors *const sf =
        get_ref_scale_factors_const(cm, ref_frame);
    av1_setup_pred_block(xd, yv12_mb[ref_frame], yv12, sf, sf, num_planes);
    av1_find_mv_refs(cm, xd, mbmi, ref_frame, mbmi_ext->ref_mv_count,
                     xd->ref_mv_stack, xd->weight, NULL, mbmi_ext->global_mvs,
                     mbmi_ext->mode_context);
    // TODO(Ravi): Populate mbmi_ext->ref_mv_stack[ref_frame][4] and
    // mbmi_ext->weight[ref_frame][4] inside av1_find_mv_refs.
    av1_copy_usable_ref_mv_stack_and_weight(xd, mbmi_ext, ref_frame);
    av1_find_best_ref_mvs_from_stack(
        cm->features.allow_high_precision_mv, mbmi_ext, ref_frame,
        &frame_mv[NEARESTMV][ref_frame], &frame_mv[NEARMV][ref_frame], 0);
    frame_mv[GLOBALMV][ref_frame] = mbmi_ext->global_mvs[ref_frame];
    // Early exit for non-LAST frame if force_skip_low_temp_var is set.
    if (!av1_is_scaled(sf) && bsize >= BLOCK_8X8 && !skip_pred_mv &&
        !(force_skip_low_temp_var && ref_frame != LAST_FRAME)) {
      av1_mv_pred(cpi, x, yv12_mb[ref_frame][0].buf, yv12->y_stride, ref_frame,
                  bsize);
    }
  }
  if (cm->features.switchable_motion_mode) {
    av1_count_overlappable_neighbors(cm, xd);
  }
  mbmi->num_proj_ref = 1;
}

static INLINE int early_term_inter_search_with_sse(int early_term_idx,
                                                   BLOCK_SIZE bsize,
                                                   int64_t this_sse,
                                                   int64_t best_sse,
                                                   PREDICTION_MODE this_mode) {
  // Aggressiveness to terminate inter mode search early is adjusted based on
  // speed and block size.
  static const double early_term_thresh[4][4] = { { 0.65, 0.65, 0.65, 0.7 },
                                                  { 0.6, 0.65, 0.85, 0.9 },
                                                  { 0.5, 0.5, 0.55, 0.6 },
                                                  { 0.6, 0.75, 0.85, 0.85 } };
  static const double early_term_thresh_newmv_nearestmv[4] = { 0.3, 0.3, 0.3,
                                                               0.3 };

  const int size_group = size_group_lookup[bsize];
  assert(size_group < 4);
  assert((early_term_idx > 0) && (early_term_idx < EARLY_TERM_INDICES));
  const double threshold =
      ((early_term_idx == EARLY_TERM_IDX_4) &&
       (this_mode == NEWMV || this_mode == NEARESTMV))
          ? early_term_thresh_newmv_nearestmv[size_group]
          : early_term_thresh[early_term_idx - 1][size_group];

  // Terminate inter mode search early based on best sse so far.
  if ((early_term_idx > 0) && (threshold * this_sse > best_sse)) {
    return 1;
  }
  return 0;
}

static INLINE void init_best_pickmode(BEST_PICKMODE *bp) {
  bp->best_sse = INT64_MAX;
  bp->best_mode = NEARESTMV;
  bp->best_ref_frame = LAST_FRAME;
  bp->best_second_ref_frame = NONE_FRAME;
  bp->best_tx_size = TX_8X8;
  bp->tx_type = DCT_DCT;
  bp->best_pred_filter = av1_broadcast_interp_filter(EIGHTTAP_REGULAR);
  bp->best_mode_skip_txfm = 0;
  bp->best_mode_initial_skip_flag = 0;
  bp->best_pred = NULL;
  bp->best_motion_mode = SIMPLE_TRANSLATION;
  bp->num_proj_ref = 0;
  av1_zero(bp->wm_params);
  av1_zero(bp->blk_skip);
  av1_zero(bp->pmi);
}

// Copy best inter mode parameters to best_pickmode
static INLINE void update_search_state_nonrd(
    InterModeSearchStateNonrd *search_state, MB_MODE_INFO *const mi,
    TxfmSearchInfo *txfm_info, RD_STATS *nonskip_rdc, PICK_MODE_CONTEXT *ctx,
    PREDICTION_MODE this_best_mode, const int64_t sse_y) {
  BEST_PICKMODE *const best_pickmode = &search_state->best_pickmode;
  const int num_8x8_blocks = ctx->num_4x4_blk / 4;

  best_pickmode->best_sse = sse_y;
  best_pickmode->best_mode = this_best_mode;
  best_pickmode->best_motion_mode = mi->motion_mode;
  best_pickmode->wm_params = mi->wm_params;
  best_pickmode->num_proj_ref = mi->num_proj_ref;
  best_pickmode->best_pred_filter = mi->interp_filters;
  best_pickmode->best_tx_size = mi->tx_size;
  best_pickmode->best_ref_frame = mi->ref_frame[0];
  best_pickmode->best_second_ref_frame = mi->ref_frame[1];
  best_pickmode->best_mode_skip_txfm = search_state->this_rdc.skip_txfm;
  best_pickmode->best_mode_initial_skip_flag =
      (nonskip_rdc->rate == INT_MAX && search_state->this_rdc.skip_txfm);
  if (!best_pickmode->best_mode_skip_txfm) {
    memcpy(best_pickmode->blk_skip, txfm_info->blk_skip,
           sizeof(txfm_info->blk_skip[0]) * num_8x8_blocks);
  }
}

static INLINE int subpel_select(AV1_COMP *cpi, MACROBLOCK *x, BLOCK_SIZE bsize,
                                int_mv *mv, MV ref_mv, FULLPEL_MV start_mv,
                                bool fullpel_performed_well) {
  const int frame_lowmotion = cpi->rc.avg_frame_low_motion;
  const int reduce_mv_pel_precision_highmotion =
      cpi->sf.rt_sf.reduce_mv_pel_precision_highmotion;

  // Reduce MV precision for higher int MV value & frame-level motion
  if (reduce_mv_pel_precision_highmotion >= 3) {
    int mv_thresh = 4;
    const int is_low_resoln =
        (cpi->common.width * cpi->common.height <= 320 * 240);
    mv_thresh = (bsize > BLOCK_32X32) ? 2 : (bsize > BLOCK_16X16) ? 4 : 6;
    if (frame_lowmotion > 0 && frame_lowmotion < 40) mv_thresh = 12;
    mv_thresh = (is_low_resoln) ? mv_thresh >> 1 : mv_thresh;
    if (abs(mv->as_fullmv.row) >= mv_thresh ||
        abs(mv->as_fullmv.col) >= mv_thresh)
      return HALF_PEL;
  } else if (reduce_mv_pel_precision_highmotion >= 1) {
    int mv_thresh;
    const int th_vals[2][3] = { { 4, 8, 10 }, { 4, 6, 8 } };
    const int th_idx = reduce_mv_pel_precision_highmotion - 1;
    assert(th_idx >= 0 && th_idx < 2);
    if (frame_lowmotion > 0 && frame_lowmotion < 40)
      mv_thresh = 12;
    else
      mv_thresh = (bsize >= BLOCK_32X32)   ? th_vals[th_idx][0]
                  : (bsize >= BLOCK_16X16) ? th_vals[th_idx][1]
                                           : th_vals[th_idx][2];
    if (abs(mv->as_fullmv.row) >= (mv_thresh << 1) ||
        abs(mv->as_fullmv.col) >= (mv_thresh << 1))
      return FULL_PEL;
    else if (abs(mv->as_fullmv.row) >= mv_thresh ||
             abs(mv->as_fullmv.col) >= mv_thresh)
      return HALF_PEL;
  }
  // Reduce MV precision for relatively static (e.g. background), low-complex
  // large areas
  if (cpi->sf.rt_sf.reduce_mv_pel_precision_lowcomplex >= 2) {
    const int qband = x->qindex >> (QINDEX_BITS - 2);
    assert(qband < 4);
    if (x->content_state_sb.source_sad_nonrd <= kVeryLowSad &&
        bsize > BLOCK_16X16 && qband != 0) {
      if (x->source_variance < 500)
        return FULL_PEL;
      else if (x->source_variance < 5000)
        return HALF_PEL;
    }
  } else if (cpi->sf.rt_sf.reduce_mv_pel_precision_lowcomplex >= 1) {
    if (fullpel_performed_well && ref_mv.row == 0 && ref_mv.col == 0 &&
        start_mv.row == 0 && start_mv.col == 0)
      return HALF_PEL;
  }
  return cpi->sf.mv_sf.subpel_force_stop;
}

static INLINE bool use_aggressive_subpel_search_method(
    MACROBLOCK *x, bool use_adaptive_subpel_search,
    bool fullpel_performed_well) {
  if (!use_adaptive_subpel_search) return false;
  const int qband = x->qindex >> (QINDEX_BITS - 2);
  assert(qband < 4);
  if ((qband > 0) && (fullpel_performed_well ||
                      (x->content_state_sb.source_sad_nonrd <= kLowSad) ||
                      (x->source_variance < 100)))
    return true;
  return false;
}

static INLINE void set_force_skip_flag(const AV1_COMP *const cpi,
                                       MACROBLOCK *const x, unsigned int sse,
                                       int *force_skip) {
  if (x->txfm_search_params.tx_mode_search_type == TX_MODE_SELECT &&
      cpi->sf.rt_sf.tx_size_level_based_on_qstep &&
      cpi->sf.rt_sf.tx_size_level_based_on_qstep >= 2) {
    const int qstep = x->plane[AOM_PLANE_Y].dequant_QTX[1] >> (x->e_mbd.bd - 5);
    const unsigned int qstep_sq = qstep * qstep;
    // If the sse is low for low source variance blocks, mark those as
    // transform skip.
    // Note: Though qstep_sq is based on ac qstep, the threshold is kept
    // low so that reliable early estimate of tx skip can be obtained
    // through its comparison with sse.
    if (sse < qstep_sq && x->source_variance < qstep_sq &&
        x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] == 0 &&
        x->color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)] == 0)
      *force_skip = 1;
  }
}

// Adjust the ac_thr according to speed, width, height and normalized sum
static INLINE int ac_thr_factor(int speed, int width, int height,
                                int norm_sum) {
  if (speed >= 8 && norm_sum < 5) {
    if (width <= 640 && height <= 480)
      return 4;
    else
      return 2;
  }
  return 1;
}

// Explicitly enumerate the cases so the compiler can generate SIMD for the
// function. According to the disassembler, gcc generates SSE codes for each of
// the possible block sizes. The hottest case is tx_width 16, which takes up
// about 8% of the self cycle of av1_nonrd_pick_inter_mode_sb. Since
// av1_nonrd_pick_inter_mode_sb takes up about 3% of total encoding time, the
// potential room of improvement for writing AVX2 optimization is only 3% * 8% =
// 0.24% of total encoding time.
static AOM_INLINE void scale_square_buf_vals(int16_t *dst, int tx_width,
                                             const int16_t *src,
                                             int src_stride) {
#define DO_SCALING                                                   \
  do {                                                               \
    for (int idy = 0; idy < tx_width; ++idy) {                       \
      for (int idx = 0; idx < tx_width; ++idx) {                     \
        dst[idy * tx_width + idx] = src[idy * src_stride + idx] * 8; \
      }                                                              \
    }                                                                \
  } while (0)

  if (tx_width == 4) {
    DO_SCALING;
  } else if (tx_width == 8) {
    DO_SCALING;
  } else if (tx_width == 16) {
    DO_SCALING;
  } else {
    assert(0);
  }

#undef DO_SCALING
}

static INLINE void init_mbmi(MB_MODE_INFO *mbmi, PREDICTION_MODE pred_mode,
                             MV_REFERENCE_FRAME ref_frame0,
                             MV_REFERENCE_FRAME ref_frame1,
                             const AV1_COMMON *cm) {
  PALETTE_MODE_INFO *const pmi = &mbmi->palette_mode_info;
  mbmi->ref_mv_idx = 0;
  mbmi->mode = pred_mode;
  mbmi->uv_mode = UV_DC_PRED;
  mbmi->ref_frame[0] = ref_frame0;
  mbmi->ref_frame[1] = ref_frame1;
  pmi->palette_size[PLANE_TYPE_Y] = 0;
  pmi->palette_size[PLANE_TYPE_UV] = 0;
  mbmi->filter_intra_mode_info.use_filter_intra = 0;
  mbmi->mv[0].as_int = mbmi->mv[1].as_int = 0;
  mbmi->motion_mode = SIMPLE_TRANSLATION;
  mbmi->num_proj_ref = 1;
  mbmi->interintra_mode = 0;
  set_default_interp_filters(mbmi, cm->features.interp_filter);
}

static INLINE int get_pred_buffer(PRED_BUFFER *p, int len) {
  for (int buf_idx = 0; buf_idx < len; buf_idx++) {
    if (!p[buf_idx].in_use) {
      p[buf_idx].in_use = 1;
      return buf_idx;
    }
  }
  return -1;
}

static INLINE void free_pred_buffer(PRED_BUFFER *p) {
  if (p != NULL) p->in_use = 0;
}

static INLINE int get_drl_cost(PREDICTION_MODE this_mode, int ref_mv_idx,
                               const MB_MODE_INFO_EXT *mbmi_ext,
                               const int (*const drl_mode_cost0)[2],
                               int8_t ref_frame_type) {
  int cost = 0;
  if (this_mode == NEWMV || this_mode == NEW_NEWMV) {
    for (int idx = 0; idx < 2; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);
        cost += drl_mode_cost0[drl_ctx][ref_mv_idx != idx];
        if (ref_mv_idx == idx) return cost;
      }
    }
    return cost;
  }

  if (have_nearmv_in_inter_mode(this_mode)) {
    for (int idx = 1; idx < 3; ++idx) {
      if (mbmi_ext->ref_mv_count[ref_frame_type] > idx + 1) {
        uint8_t drl_ctx = av1_drl_ctx(mbmi_ext->weight[ref_frame_type], idx);
        cost += drl_mode_cost0[drl_ctx][ref_mv_idx != (idx - 1)];
        if (ref_mv_idx == (idx - 1)) return cost;
      }
    }
    return cost;
  }
  return cost;
}

static INLINE void update_thresh_freq_fact(AV1_COMP *cpi, MACROBLOCK *x,
                                           BLOCK_SIZE bsize,
                                           MV_REFERENCE_FRAME ref_frame,
                                           THR_MODES best_mode_idx,
                                           PREDICTION_MODE mode) {
  const THR_MODES thr_mode_idx = mode_idx[ref_frame][mode_offset(mode)];
  const BLOCK_SIZE min_size = AOMMAX(bsize - 3, BLOCK_4X4);
  const BLOCK_SIZE max_size = AOMMIN(bsize + 6, BLOCK_128X128);
  for (BLOCK_SIZE bs = min_size; bs <= max_size; bs += 3) {
    int *freq_fact = &x->thresh_freq_fact[bs][thr_mode_idx];
    if (thr_mode_idx == best_mode_idx) {
      *freq_fact -= (*freq_fact >> 4);
    } else {
      *freq_fact =
          AOMMIN(*freq_fact + RD_THRESH_INC,
                 cpi->sf.inter_sf.adaptive_rd_thresh * RD_THRESH_MAX_FACT);
    }
  }
}

#if !CONFIG_REALTIME_ONLY
static AOM_INLINE int is_warped_mode_allowed(const AV1_COMP *cpi,
                                             MACROBLOCK *const x,
                                             const MB_MODE_INFO *mbmi) {
  const FeatureFlags *const features = &cpi->common.features;
  const MACROBLOCKD *xd = &x->e_mbd;

  if (cpi->sf.inter_sf.extra_prune_warped) return 0;
  if (has_second_ref(mbmi)) return 0;
  MOTION_MODE last_motion_mode_allowed = SIMPLE_TRANSLATION;

  if (features->switchable_motion_mode) {
    // Determine which motion modes to search if more than SIMPLE_TRANSLATION
    // is allowed.
    last_motion_mode_allowed = motion_mode_allowed(
        xd->global_motion, xd, mbmi, features->allow_warped_motion);
  }

  if (last_motion_mode_allowed == WARPED_CAUSAL) {
    return 1;
  }

  return 0;
}
#endif

static INLINE bool should_prune_intra_modes_using_neighbors(
    const MACROBLOCKD *xd, bool enable_intra_mode_pruning_using_neighbors,
    PREDICTION_MODE this_mode, PREDICTION_MODE above_mode,
    PREDICTION_MODE left_mode) {
  if (!enable_intra_mode_pruning_using_neighbors) return false;

  // Avoid pruning of DC_PRED as it is the most probable mode to win as per the
  // statistics generated for nonrd intra mode evaluations.
  if (this_mode == DC_PRED) return false;

  // Enable the pruning for current mode only if it is not the winner mode of
  // both the neighboring blocks (left/top).
  return xd->up_available && this_mode != above_mode && xd->left_available &&
         this_mode != left_mode;
}

static AOM_INLINE int is_same_gf_and_last_scale(AV1_COMMON *cm) {
  struct scale_factors *const sf_last = get_ref_scale_factors(cm, LAST_FRAME);
  struct scale_factors *const sf_golden =
      get_ref_scale_factors(cm, GOLDEN_FRAME);
  return ((sf_last->x_scale_fp == sf_golden->x_scale_fp) &&
          (sf_last->y_scale_fp == sf_golden->y_scale_fp));
}

// Checks whether Intra mode needs to be pruned based on
// 'intra_y_mode_bsize_mask_nrd' and 'prune_hv_pred_modes_using_blksad'
// speed features.
static INLINE bool is_prune_intra_mode(AV1_COMP *cpi, int mode_index,
                                       int force_intra_check, BLOCK_SIZE bsize,
                                       uint8_t segment_id,
                                       SOURCE_SAD source_sad_nonrd,
                                       uint8_t color_sensitivity[2]) {
  const PREDICTION_MODE this_mode = intra_mode_list[mode_index];
  if (mode_index > 2 || force_intra_check == 0) {
    if (!((1 << this_mode) & cpi->sf.rt_sf.intra_y_mode_bsize_mask_nrd[bsize]))
      return true;

    if (this_mode == DC_PRED) return false;

    if (!cpi->sf.rt_sf.prune_hv_pred_modes_using_src_sad) return false;

    const bool has_color_sensitivity =
        color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_U)] &&
        color_sensitivity[COLOR_SENS_IDX(AOM_PLANE_V)];
    if (has_color_sensitivity &&
        (cpi->rc.frame_source_sad > 1.1 * cpi->rc.avg_source_sad ||
         cyclic_refresh_segment_id_boosted(segment_id) ||
         source_sad_nonrd > kMedSad))
      return false;

    return true;
  }
  return false;
}

static AOM_INLINE int is_filter_search_enabled_blk(
    AV1_COMP *cpi, MACROBLOCK *x, int mi_row, int mi_col, BLOCK_SIZE bsize,
    int segment_id, int cb_pred_filter_search, InterpFilter *filt_select) {
  const AV1_COMMON *const cm = &cpi->common;
  // filt search disabled
  if (!cpi->sf.rt_sf.use_nonrd_filter_search) return 0;
  // filt search purely based on mode properties
  if (!cb_pred_filter_search) return 1;
  MACROBLOCKD *const xd = &x->e_mbd;
  int enable_interp_search = 0;
  if (!(xd->left_mbmi && xd->above_mbmi)) {
    // neighbors info unavailable
    enable_interp_search = 2;
  } else if (!(is_inter_block(xd->left_mbmi) &&
               is_inter_block(xd->above_mbmi))) {
    // neighbor is INTRA
    enable_interp_search = 2;
  } else if (xd->left_mbmi->interp_filters.as_int !=
             xd->above_mbmi->interp_filters.as_int) {
    // filters are different
    enable_interp_search = 2;
  } else if ((cb_pred_filter_search == 1) &&
             (xd->left_mbmi->interp_filters.as_filters.x_filter !=
              EIGHTTAP_REGULAR)) {
    // not regular
    enable_interp_search = 2;
  } else {
    // enable prediction based on chessboard pattern
    if (xd->left_mbmi->interp_filters.as_filters.x_filter == EIGHTTAP_SMOOTH)
      *filt_select = EIGHTTAP_SMOOTH;
    const int bsl = mi_size_wide_log2[bsize];
    enable_interp_search =
        (bool)((((mi_row + mi_col) >> bsl) +
                get_chessboard_index(cm->current_frame.frame_number)) &
               0x1);
    if (cyclic_refresh_segment_id_boosted(segment_id)) enable_interp_search = 1;
  }
  return enable_interp_search;
}

static INLINE void set_compound_mode(MACROBLOCK *x,
                                     MV_REFERENCE_FRAME ref_frame,
                                     MV_REFERENCE_FRAME ref_frame2,
                                     int ref_mv_idx,
                                     int_mv frame_mv[MB_MODE_COUNT][REF_FRAMES],
                                     PREDICTION_MODE this_mode) {
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mi = xd->mi[0];
  mi->ref_frame[0] = ref_frame;
  mi->ref_frame[1] = ref_frame2;
  mi->compound_idx = 1;
  mi->comp_group_idx = 0;
  mi->interinter_comp.type = COMPOUND_AVERAGE;
  MV_REFERENCE_FRAME ref_frame_comp = av1_ref_frame_type(mi->ref_frame);
  if (this_mode == GLOBAL_GLOBALMV) {
    frame_mv[this_mode][ref_frame].as_int = 0;
    frame_mv[this_mode][ref_frame2].as_int = 0;
  } else if (this_mode == NEAREST_NEARESTMV) {
    frame_mv[this_mode][ref_frame].as_int =
        xd->ref_mv_stack[ref_frame_comp][0].this_mv.as_int;
    frame_mv[this_mode][ref_frame2].as_int =
        xd->ref_mv_stack[ref_frame_comp][0].comp_mv.as_int;
  } else if (this_mode == NEAR_NEARMV) {
    frame_mv[this_mode][ref_frame].as_int =
        xd->ref_mv_stack[ref_frame_comp][ref_mv_idx].this_mv.as_int;
    frame_mv[this_mode][ref_frame2].as_int =
        xd->ref_mv_stack[ref_frame_comp][ref_mv_idx].comp_mv.as_int;
  }
}

static AOM_INLINE bool is_globalmv_better(
    PREDICTION_MODE this_mode, MV_REFERENCE_FRAME ref_frame, int rate_mv,
    const ModeCosts *mode_costs,
    const int (*single_inter_mode_costs)[REF_FRAMES],
    const MB_MODE_INFO_EXT *mbmi_ext) {
  const int globalmv_mode_cost =
      single_inter_mode_costs[INTER_OFFSET(GLOBALMV)][ref_frame];
  int this_mode_cost =
      rate_mv + single_inter_mode_costs[INTER_OFFSET(this_mode)][ref_frame];
  if (this_mode == NEWMV || this_mode == NEARMV) {
    const MV_REFERENCE_FRAME rf[2] = { ref_frame, NONE_FRAME };
    this_mode_cost += get_drl_cost(
        NEWMV, 0, mbmi_ext, mode_costs->drl_mode_cost0, av1_ref_frame_type(rf));
  }
  return this_mode_cost > globalmv_mode_cost;
}

static AOM_INLINE bool previous_mode_performed_poorly(
    PREDICTION_MODE mode, MV_REFERENCE_FRAME ref_frame,
    const unsigned int (*vars)[REF_FRAMES],
    const int64_t (*uv_dist)[REF_FRAMES]) {
  unsigned int best_var = UINT_MAX;
  int64_t best_uv_dist = INT64_MAX;
  for (int midx = 0; midx < RTC_INTER_MODES; midx++) {
    best_var = AOMMIN(best_var, vars[midx][ref_frame]);
    best_uv_dist = AOMMIN(best_uv_dist, uv_dist[midx][ref_frame]);
  }
  assert(best_var != UINT_MAX && "Invalid variance data.");
  const float mult = 1.125f;
  bool var_bad = mult * best_var < vars[INTER_OFFSET(mode)][ref_frame];
  if (uv_dist[INTER_OFFSET(mode)][ref_frame] < INT64_MAX &&
      best_uv_dist != uv_dist[INTER_OFFSET(mode)][ref_frame]) {
    // If we have chroma info, then take it into account
    var_bad &= mult * best_uv_dist < uv_dist[INTER_OFFSET(mode)][ref_frame];
  }
  return var_bad;
}

#endif  // AOM_AV1_ENCODER_NONRD_OPT_H_
