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

#ifndef AV1_COMMON_SCAN_H_
#define AV1_COMMON_SCAN_H_

#include "aom/aom_integer.h"
#include "aom_ports/mem.h"

#include "av1/common/enums.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/blockd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_NEIGHBORS 2

extern const SCAN_ORDER av1_default_scan_orders[TX_SIZES];
extern const SCAN_ORDER av1_intra_scan_orders[TX_SIZES_ALL][TX_TYPES];
extern const SCAN_ORDER av1_inter_scan_orders[TX_SIZES_ALL][TX_TYPES];

#if CONFIG_ADAPT_SCAN
#define USE_2X2_PROB 0
#define CACHE_SCAN_PROB 1
#define REDUCED_SET 1
#define SUB_REGION_COUNT 1
#define SUB_FRAME_COUNT 0
#define SIG_REGION 0
#define UNI_RECT 1
#define USE_TOPOLOGICAL_SORT 0
#define USE_LIMIT_SCAN_DISTANCE 0
void av1_update_scan_count_facade(const AV1_COMMON *const cm, MACROBLOCKD *xd,
                                  int mi_row, TX_SIZE tx_size, TX_TYPE tx_type,
                                  const tran_low_t *dqcoeffs, int max_scan);

// embed r + c and coeff_idx info with nonzero probabilities. When sorting the
// nonzero probabilities, if there is a tie, the coefficient with smaller r + c
// will be scanned first
void av1_augment_prob(TX_SIZE tx_size, TX_TYPE tx_type, uint32_t *prob);

#if USE_TOPOLOGICAL_SORT
// apply quick sort on nonzero probabilities to obtain a sort order
void av1_update_sort_order(TX_SIZE tx_size, TX_TYPE tx_type,
                           const uint32_t *non_zero_prob, int16_t *sort_order);

// apply topological sort on the nonzero probabilities sorting order to
// guarantee each to-be-scanned coefficient's upper and left coefficient will be
// scanned before the to-be-scanned coefficient.
void av1_update_scan_order(TX_SIZE tx_size, int16_t *sort_order, int16_t *scan,
                           int16_t *iscan);
#else   // USE_TOPOLOGICAL_SORT
void av1_update_scan_order(TX_SIZE tx_size, TX_TYPE tx_type,
                           uint32_t *non_zero_prob, int16_t *scan,
                           int16_t *iscan);
#endif  // USE_TOPOLOGICAL_SORT

// For each coeff_idx in scan[], update its above and left neighbors in
// neighbors[] accordingly.
void av1_update_neighbors(TX_SIZE tx_size, const int16_t *scan,
                          const int16_t *iscan, int16_t *neighbors);
void av1_init_scan_order(AV1_COMMON *cm);
void av1_adapt_scan_order(AV1_COMMON *cm, FRAME_CONTEXT *ec_ctxs[],
                          int num_tiles);
#if USE_2X2_PROB
void av1_down_sample_scan_count(uint32_t *non_zero_count_ds,
                                const uint32_t *non_zero_count,
                                TX_SIZE tx_size);
#endif  // USE_2X2_PROB
#endif  // CONFIG_ADAPT_SCAN
void av1_deliver_eob_threshold(const AV1_COMMON *cm, MACROBLOCKD *xd);

static INLINE int get_coef_context(const int16_t *neighbors,
                                   const uint8_t *token_cache, int c) {
  return (1 + token_cache[neighbors[MAX_NEIGHBORS * c + 0]] +
          token_cache[neighbors[MAX_NEIGHBORS * c + 1]]) >>
         1;
}

static INLINE const SCAN_ORDER *get_default_scan(TX_SIZE tx_size,
                                                 TX_TYPE tx_type,
                                                 int is_inter) {
  return is_inter ? &av1_inter_scan_orders[tx_size][tx_type]
                  : &av1_intra_scan_orders[tx_size][tx_type];
}

static INLINE int do_adapt_scan(TX_SIZE tx_size, TX_TYPE tx_type) {
  (void)tx_size;
  if (tx_size_2d[tx_size] >= 1024 && tx_type != DCT_DCT) return 0;
  if (tx_size > TX_32X16) return 0;
#if CONFIG_ADAPT_SCAN
#if REDUCED_SET
  const int txw = tx_size_wide[tx_size];
  const int txh = tx_size_high[tx_size];

  if (txw == 16 && txh == 16) return tx_type == DCT_DCT;
  if (txw >= 16 || txh >= 16) return tx_type <= ADST_ADST;
#endif
#endif
  return tx_type < IDTX;
}

#if CONFIG_BLOCK_ADAPT_SCAN
static INLINE const SCAN_ORDER *get_scan(const FRAME_CONTEXT *ec_ctx,
                                         TX_SIZE tx_size, TX_TYPE tx_type,
                                         const MB_MODE_INFO *mbmi) {
  const int is_inter = is_inter_block(mbmi);
  const SCAN_ORDER *scan_order = 0;
  if (ec_ctx) {
    scan_order = is_inter ? ec_ctx->adapt_scan_order_inter[tx_size][tx_type]
                          : ec_ctx->adapt_scan_order_intra[tx_size][tx_type];
  }
  if (!scan_order) {
    scan_order = get_default_scan(tx_size, tx_type, is_inter);
  }
  return scan_order;
}
#else
static INLINE const SCAN_ORDER *get_scan(const AV1_COMMON *cm, TX_SIZE tx_size,
                                         TX_TYPE tx_type,
                                         const MB_MODE_INFO *mbmi) {
  const int is_inter = is_inter_block(mbmi);
#if CONFIG_ADAPT_SCAN
  (void)mbmi;
  (void)is_inter;
  if (!do_adapt_scan(tx_size, tx_type))
    return get_default_scan(tx_size, tx_type, is_inter);
  else
    return &cm->fc->sc[tx_size][tx_type];
#else   // CONFIG_ADAPT_SCAN
  (void)cm;
  return get_default_scan(tx_size, tx_type, is_inter);
#endif  // CONFIG_ADAPT_SCAN
}
#endif

#if CONFIG_BLOCK_ADAPT_SCAN
extern const int8_t av1_adapt_scan_class_inter[TX_SIZES_ALL][TX_TYPES];
extern const int8_t av1_adapt_scan_class_intra[TX_SIZES_ALL][TX_TYPES];
extern const SCAN_ORDER av1_adapt_scan_order[TX_SIZES_ALL][ADAPT_SCAN_ORDERS];

static const int16_t scan1_weight_4x4[7] = {
  0, 8, 17, 26, 34, 40, 43,
};
static const int16_t scan2_weight_4x4[10] = {
  0, 6, 12, 18, 24, 29, 33, 38, 41, 43,
};
static const int16_t scan3_weight_4x4[13] = {
  0, 6, 10, 15, 20, 22, 26, 30, 32, 35, 39, 41, 43,
};
static const int16_t scan1_weight_8x8[15] = {
  0, 8, 17, 26, 37, 48, 60, 72, 84, 94, 101, 108, 112, 116, 118,
};
static const int16_t scan2_weight_8x8[22] = {
  0,  6,  12, 18, 25, 32,  40,  48,  55,  62,  69,
  75, 82, 88, 94, 99, 104, 108, 112, 114, 116, 118,
};
static const int16_t scan3_weight_8x8[29] = {
  0,  6,  10, 15, 21, 26, 32,  39,  44,  49,  54,  59,  63,  68,  72,
  76, 81, 85, 88, 93, 96, 100, 104, 108, 110, 113, 115, 116, 118,
};
static const int16_t scan1_weight_16x16[31] = {
  0,   8,   17,  26,  37,  48,  60,  72,  85,  99,  114,
  129, 144, 160, 176, 193, 209, 224, 237, 249, 260, 270,
  278, 286, 293, 299, 303, 307, 310, 313, 314,
};
static const int16_t scan2_weight_16x16[46] = {
  0,   6,   12,  18,  25,  32,  40,  48,  56,  65,  74,  83,
  93,  103, 113, 123, 134, 143, 153, 162, 171, 180, 189, 197,
  206, 214, 222, 230, 238, 245, 253, 260, 267, 274, 280, 285,
  290, 295, 299, 302, 305, 308, 310, 312, 313, 314,
};
static const int16_t scan3_weight_16x16[61] = {
  0,   6,   10,  15,  21,  26,  32,  39,  45,  51,  59,  66,  73,
  81,  88,  96,  104, 111, 118, 125, 132, 139, 145, 151, 158, 164,
  170, 176, 182, 187, 193, 199, 204, 210, 215, 220, 226, 231, 236,
  241, 247, 252, 257, 262, 267, 272, 277, 281, 286, 290, 293, 297,
  300, 303, 305, 307, 309, 311, 312, 313, 314,
};
static const int16_t scan1_weight_32x32[63] = {
  0,   8,   17,  26,  37,  48,  60,  72,  85,  99,  114, 129, 144,
  160, 176, 193, 210, 228, 246, 264, 283, 302, 322, 341, 361, 382,
  403, 424, 445, 467, 489, 511, 533, 554, 574, 592, 610, 628, 644,
  659, 674, 688, 701, 714, 725, 736, 747, 757, 766, 774, 782, 790,
  796, 802, 808, 813, 817, 821, 824, 826, 828, 830, 831,
};
static const int16_t scan2_weight_32x32[94] = {
  0,   6,   12,  18,  25,  32,  40,  48,  56,  65,  74,  83,  93,  103,
  113, 123, 134, 145, 156, 168, 179, 191, 203, 216, 228, 241, 254, 267,
  280, 294, 307, 321, 335, 348, 361, 374, 387, 399, 412, 424, 436, 448,
  460, 471, 483, 494, 506, 517, 528, 539, 550, 561, 571, 582, 592, 603,
  613, 624, 634, 644, 654, 664, 674, 684, 694, 703, 712, 720, 728, 736,
  744, 751, 758, 764, 770, 776, 782, 787, 792, 797, 801, 805, 809, 812,
  815, 818, 821, 823, 825, 827, 828, 829, 830, 831,
};
static const int16_t scan3_weight_32x32[125] = {
  0,   6,   10,  15,  21,  26,  32,  39,  45,  51,  59,  66,  73,  81,
  88,  96,  105, 113, 122, 130, 139, 148, 158, 167, 176, 186, 196, 206,
  216, 226, 236, 247, 257, 267, 277, 286, 295, 305, 314, 323, 333, 341,
  350, 359, 368, 376, 385, 393, 401, 410, 418, 426, 434, 442, 450, 458,
  465, 473, 481, 489, 496, 504, 511, 518, 526, 533, 541, 548, 555, 562,
  570, 577, 584, 591, 598, 605, 612, 619, 625, 632, 639, 646, 653, 659,
  666, 673, 679, 686, 693, 699, 705, 712, 718, 725, 731, 738, 744, 750,
  755, 761, 766, 771, 776, 780, 785, 789, 793, 797, 800, 803, 807, 810,
  812, 815, 817, 820, 822, 823, 825, 826, 828, 829, 830, 830, 831,
};

static const int16_t *const scan_cost_weight[TX_SIZES_ALL][3] = {
  { scan1_weight_4x4, scan2_weight_4x4, scan3_weight_4x4 },        // TX_4X4
  { scan1_weight_8x8, scan2_weight_8x8, scan3_weight_8x8 },        // TX_8X8
  { scan1_weight_16x16, scan2_weight_16x16, scan3_weight_16x16 },  // TX_16X16
  { scan1_weight_32x32, scan2_weight_32x32, scan3_weight_32x32 },  // TX_32X32
#if CONFIG_TX64X64
  { 0, 0, 0 },  // TX_64X64
#endif
  { 0, 0, 0 },  // TX_4X8
  { 0, 0, 0 },  // TX_8X4
  { 0, 0, 0 },  // TX_8X16
  { 0, 0, 0 },  // TX_16X8
  { 0, 0, 0 },  // TX_16X32
  { 0, 0, 0 },  // TX_32X16
#if CONFIG_TX64X64
  { 0, 0, 0 },  // TX_32X64
  { 0, 0, 0 },  // TX_64X32
#endif
};

static INLINE int av1_do_adapt_scan(const AV1_COMMON *cm, int tx_size,
                                    int tx_type, int is_inter) {
  int scan_class = is_inter ? av1_adapt_scan_class_inter[tx_size][tx_type]
                            : av1_adapt_scan_class_intra[tx_size][tx_type];
  return cm->use_adapt_scan && scan_class >= 0;
}

static INLINE void av1_update_scan_costs(int coeff_idx, int bwl,
                                         int *scan_cost) {
  if (coeff_idx > 0) {
    const int y = coeff_idx >> bwl;
    const int x = coeff_idx - (y << bwl);
    scan_cost[0] = AOMMAX(scan_cost[0], x + 3 * y);
    scan_cost[1] = AOMMAX(scan_cost[1], x + 2 * y);
    scan_cost[2] = AOMMAX(scan_cost[2], x + y);
    scan_cost[3] = AOMMAX(scan_cost[3], 2 * x + y);
    scan_cost[4] = AOMMAX(scan_cost[4], 3 * x + y);
  }
}

static INLINE void av1_accumulate_scan_costs(FRAME_CONTEXT *ec_ctx, int tx_size,
                                             int tx_type, int is_inter,
                                             int *scan_cost) {
  int scan_class = is_inter ? av1_adapt_scan_class_inter[tx_size][tx_type]
                            : av1_adapt_scan_class_intra[tx_size][tx_type];

  const int16_t *const *weight = scan_cost_weight[tx_size];
  // direction dependent weight table to make the costs comparable
  ec_ctx->adapt_scan_cost[scan_class][0] += weight[2][scan_cost[0]];  // 1:3
  ec_ctx->adapt_scan_cost[scan_class][1] += weight[1][scan_cost[1]];  // 1:2
  ec_ctx->adapt_scan_cost[scan_class][2] += weight[0][scan_cost[2]];  // 1:1
  ec_ctx->adapt_scan_cost[scan_class][3] += weight[1][scan_cost[3]];  // 2:1
  ec_ctx->adapt_scan_cost[scan_class][4] += weight[2][scan_cost[4]];  // 3:1
}

SCAN_CONTEXT av1_evaluate_scan_costs(const FRAME_CONTEXT *ec_ctx);
void av1_select_scan_order(FRAME_CONTEXT *ec_ctx, SCAN_CONTEXT left[2],
                           SCAN_CONTEXT above[2]);
void av1_clear_scan_cost(FRAME_CONTEXT *ec_ctx);

#endif  // CONFIG_BLOCK_ADAPT_SCAN

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_SCAN_H_
