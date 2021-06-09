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

#ifndef AOM_AV1_COMMON_ENTROPYMODE_H_
#define AOM_AV1_COMMON_ENTROPYMODE_H_

#include "av1/common/entropy.h"
#include "av1/common/entropymv.h"
#include "av1/common/filter.h"
#include "av1/common/seg_common.h"
#include "aom_dsp/aom_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCK_SIZE_GROUPS 4

#define TX_SIZE_CONTEXTS 3

#if CONFIG_NEW_INTER_MODES
#define INTER_OFFSET(mode) ((mode)-NEARMV)
#define INTER_COMPOUND_OFFSET(mode) (uint8_t)((mode)-NEAR_NEARMV)
#if CONFIG_OPTFLOW_REFINEMENT
#define INTER_OPFL_OFFSET(mode) (uint8_t)((mode)-NEAR_NEARMV_OPTFLOW)
#endif  // CONFIG_OPTFLOW_REFINEMENT
#else
#define INTER_OFFSET(mode) ((mode)-NEARESTMV)
#define INTER_COMPOUND_OFFSET(mode) (uint8_t)((mode)-NEAREST_NEARESTMV)
#endif  // CONFIG_NEW_INTER_MODES
// Number of possible contexts for a color index.
// As can be seen from av1_get_palette_color_index_context(), the possible
// contexts are (2,0,0), (2,2,1), (3,2,0), (4,1,0), (5,0,0). These are mapped to
// a value from 0 to 4 using 'palette_color_index_context_lookup' table.
#define PALETTE_COLOR_INDEX_CONTEXTS 5

// Palette Y mode context for a block is determined by number of neighboring
// blocks (top and/or left) using a palette for Y plane. So, possible Y mode'
// context values are:
// 0 if neither left nor top block uses palette for Y plane,
// 1 if exactly one of left or top block uses palette for Y plane, and
// 2 if both left and top blocks use palette for Y plane.
#define PALETTE_Y_MODE_CONTEXTS 3

// Palette UV mode context for a block is determined by whether this block uses
// palette for the Y plane. So, possible values are:
// 0 if this block doesn't use palette for Y plane.
// 1 if this block uses palette for Y plane (i.e. Y palette size > 0).
#define PALETTE_UV_MODE_CONTEXTS 2

// Map the number of pixels in a block size to a context
//   64(BLOCK_8X8, BLOCK_4x16, BLOCK_16X4)  -> 0
//  128(BLOCK_8X16, BLOCK_16x8)             -> 1
//   ...
// 4096(BLOCK_64X64)                        -> 6
#define PALATTE_BSIZE_CTXS 7

#if CONFIG_DERIVED_INTRA_MODE
#define KF_MODE_CONTEXTS 6
#define KF_IS_DR_MODE_CONTEXTS 3
#else
#define KF_MODE_CONTEXTS 5
#endif

struct AV1Common;

typedef struct {
  const int16_t *scan;
  const int16_t *iscan;
} SCAN_ORDER;

typedef struct frame_contexts {
  aom_cdf_prob txb_skip_cdf[TX_SIZES][TXB_SKIP_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob eob_extra_cdf[TX_SIZES][PLANE_TYPES][EOB_COEF_CONTEXTS]
                            [CDF_SIZE(2)];
  aom_cdf_prob dc_sign_cdf[PLANE_TYPES][DC_SIGN_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob eob_flag_cdf16[PLANE_TYPES][2][CDF_SIZE(5)];
  aom_cdf_prob eob_flag_cdf32[PLANE_TYPES][2][CDF_SIZE(6)];
  aom_cdf_prob eob_flag_cdf64[PLANE_TYPES][2][CDF_SIZE(7)];
  aom_cdf_prob eob_flag_cdf128[PLANE_TYPES][2][CDF_SIZE(8)];
  aom_cdf_prob eob_flag_cdf256[PLANE_TYPES][2][CDF_SIZE(9)];
  aom_cdf_prob eob_flag_cdf512[PLANE_TYPES][2][CDF_SIZE(10)];
  aom_cdf_prob eob_flag_cdf1024[PLANE_TYPES][2][CDF_SIZE(11)];
  aom_cdf_prob coeff_base_eob_cdf[TX_SIZES][PLANE_TYPES][SIG_COEF_CONTEXTS_EOB]
                                 [CDF_SIZE(3)];
  aom_cdf_prob coeff_base_cdf[TX_SIZES][PLANE_TYPES][SIG_COEF_CONTEXTS]
                             [CDF_SIZE(4)];
  aom_cdf_prob coeff_br_cdf[TX_SIZES][PLANE_TYPES][LEVEL_CONTEXTS]
                           [CDF_SIZE(BR_CDF_SIZE)];

#if CONFIG_NEW_INTER_MODES
  aom_cdf_prob inter_single_mode_cdf[INTER_SINGLE_MODE_CONTEXTS]
                                    [CDF_SIZE(INTER_SINGLE_MODES)];
  aom_cdf_prob drl0_cdf[DRL_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob drl1_cdf[DRL_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob drl2_cdf[DRL_MODE_CONTEXTS][CDF_SIZE(2)];
#else
  aom_cdf_prob newmv_cdf[NEWMV_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob zeromv_cdf[GLOBALMV_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob drl_cdf[DRL_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob refmv_cdf[REFMV_MODE_CONTEXTS][CDF_SIZE(2)];
#endif  // CONFIG_NEW_INTER_MODES

#if CONFIG_OPTFLOW_REFINEMENT
  aom_cdf_prob use_optflow_cdf[INTER_COMPOUND_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob inter_compound_mode_cdf[INTER_COMPOUND_MODE_CONTEXTS]
                                      [CDF_SIZE(INTER_COMPOUND_REF_TYPES)];
#else
  aom_cdf_prob inter_compound_mode_cdf[INTER_COMPOUND_MODE_CONTEXTS]
                                      [CDF_SIZE(INTER_COMPOUND_MODES)];
#endif  // CONFIG_OPTFLOW_REFINEMENT
  aom_cdf_prob compound_type_cdf[BLOCK_SIZES_ALL]
                                [CDF_SIZE(MASKED_COMPOUND_TYPES)];
  aom_cdf_prob wedge_idx_cdf[BLOCK_SIZES_ALL][CDF_SIZE(16)];
  aom_cdf_prob interintra_cdf[BLOCK_SIZE_GROUPS][CDF_SIZE(2)];
  aom_cdf_prob wedge_interintra_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob interintra_mode_cdf[BLOCK_SIZE_GROUPS]
                                  [CDF_SIZE(INTERINTRA_MODES)];
  aom_cdf_prob motion_mode_cdf[BLOCK_SIZES_ALL][CDF_SIZE(MOTION_MODES)];
#if CONFIG_EXT_ROTATION
  aom_cdf_prob warp_rotation_flag_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob warp_rotation_degree_cdf[CDF_SIZE(ROTATION_COUNT)];

  aom_cdf_prob globalmv_rotation_flag_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob globalmv_rotation_degree_cdf[CDF_SIZE(ROTATION_COUNT)];

  aom_cdf_prob translation_rotation_flag_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob translation_rotation_degree_cdf[CDF_SIZE(ROTATION_COUNT)];
#endif  // CONFIG_EXT_ROTATION
  aom_cdf_prob obmc_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob palette_y_size_cdf[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
  aom_cdf_prob palette_uv_size_cdf[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
  aom_cdf_prob palette_y_color_index_cdf[PALETTE_SIZES]
                                        [PALETTE_COLOR_INDEX_CONTEXTS]
                                        [CDF_SIZE(PALETTE_COLORS)];
  aom_cdf_prob palette_uv_color_index_cdf[PALETTE_SIZES]
                                         [PALETTE_COLOR_INDEX_CONTEXTS]
                                         [CDF_SIZE(PALETTE_COLORS)];
  aom_cdf_prob palette_y_mode_cdf[PALATTE_BSIZE_CTXS][PALETTE_Y_MODE_CONTEXTS]
                                 [CDF_SIZE(2)];
  aom_cdf_prob palette_uv_mode_cdf[PALETTE_UV_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob comp_inter_cdf[COMP_INTER_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob single_ref_cdf[REF_CONTEXTS][SINGLE_REFS - 1][CDF_SIZE(2)];
  aom_cdf_prob comp_ref_type_cdf[COMP_REF_TYPE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob uni_comp_ref_cdf[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1]
                               [CDF_SIZE(2)];
  aom_cdf_prob comp_ref_cdf[REF_CONTEXTS][FWD_REFS - 1][CDF_SIZE(2)];
  aom_cdf_prob comp_bwdref_cdf[REF_CONTEXTS][BWD_REFS - 1][CDF_SIZE(2)];
#if CONFIG_NEW_TX_PARTITION
  aom_cdf_prob inter_4way_txfm_partition_cdf[2][TXFM_PARTITION_INTER_CONTEXTS]
                                            [CDF_SIZE(4)];
  aom_cdf_prob inter_2way_txfm_partition_cdf[CDF_SIZE(2)];
  aom_cdf_prob inter_2way_rect_txfm_partition_cdf[CDF_SIZE(2)];
#else   // CONFIG_NEW_TX_PARTITION
  aom_cdf_prob txfm_partition_cdf[TXFM_PARTITION_CONTEXTS][CDF_SIZE(2)];
#endif  // CONFIG_NEW_TX_PARTITION
  aom_cdf_prob compound_index_cdf[COMP_INDEX_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob comp_group_idx_cdf[COMP_GROUP_IDX_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob skip_mode_cdfs[SKIP_MODE_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob skip_txfm_cdfs[SKIP_CONTEXTS][CDF_SIZE(2)];
  aom_cdf_prob intra_inter_cdf[INTRA_INTER_CONTEXTS][CDF_SIZE(2)];
  nmv_context nmvc;
  nmv_context ndvc;
  aom_cdf_prob intrabc_cdf[CDF_SIZE(2)];
  struct segmentation_probs seg;
  aom_cdf_prob filter_intra_cdfs[BLOCK_SIZES_ALL][CDF_SIZE(2)];
  aom_cdf_prob filter_intra_mode_cdf[CDF_SIZE(FILTER_INTRA_MODES)];
  aom_cdf_prob switchable_restore_cdf[CDF_SIZE(RESTORE_SWITCHABLE_TYPES)];
  aom_cdf_prob wiener_restore_cdf[CDF_SIZE(2)];
  aom_cdf_prob sgrproj_restore_cdf[CDF_SIZE(2)];
  aom_cdf_prob y_mode_cdf[BLOCK_SIZE_GROUPS][CDF_SIZE(INTRA_MODES)];
  aom_cdf_prob uv_mode_cdf[CFL_ALLOWED_TYPES][INTRA_MODES]
                          [CDF_SIZE(UV_INTRA_MODES)];
  aom_cdf_prob partition_cdf[PARTITION_CONTEXTS][CDF_SIZE(EXT_PARTITION_TYPES)];
#if CONFIG_EXT_RECUR_PARTITIONS
  aom_cdf_prob partition_rec_cdf[PARTITION_CONTEXTS_REC]
                                [CDF_SIZE(PARTITION_TYPES_REC)];
#endif  // CONFIG_EXT_RECUR_PARTITIONS
  aom_cdf_prob switchable_interp_cdf[SWITCHABLE_FILTER_CONTEXTS]
                                    [CDF_SIZE(SWITCHABLE_FILTERS)];
#if CONFIG_DERIVED_INTRA_MODE
  aom_cdf_prob derived_intra_mode_cdf[3][CDF_SIZE(2)];
  aom_cdf_prob uv_derived_intra_mode_cdf[2][CDF_SIZE(2)];
  aom_cdf_prob kf_is_dr_mode_cdf[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS]
                                [CDF_SIZE(2)];
  aom_cdf_prob kf_dr_mode_cdf[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS]
                             [CDF_SIZE(DIRECTIONAL_MODES)];
  aom_cdf_prob kf_none_dr_mode_cdf[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS]
                                  [CDF_SIZE(NONE_DIRECTIONAL_MODES)];
#else
  /* kf_y_cdf is discarded after use, so does not require persistent storage.
       However, we keep it with the other CDFs in this struct since it needs to
       be copied to each tile to support parallelism just like the others.
   */
  aom_cdf_prob kf_y_cdf[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS]
                       [CDF_SIZE(INTRA_MODES)];
#endif  // CONFIG_DERIVED_INTRA_MODE

  aom_cdf_prob angle_delta_cdf[DIRECTIONAL_MODES]
                              [CDF_SIZE(2 * MAX_ANGLE_DELTA + 1)];

#if CONFIG_NEW_TX_PARTITION
  aom_cdf_prob intra_4way_txfm_partition_cdf[2][TX_SIZE_CONTEXTS][CDF_SIZE(4)];
  aom_cdf_prob intra_2way_txfm_partition_cdf[CDF_SIZE(2)];
  aom_cdf_prob intra_2way_rect_txfm_partition_cdf[CDF_SIZE(2)];
#else
  aom_cdf_prob tx_size_cdf[MAX_TX_CATS][TX_SIZE_CONTEXTS]
                          [CDF_SIZE(MAX_TX_DEPTH + 1)];
#endif  // CONFIG_NEW_TX_PARTITION
  aom_cdf_prob delta_q_cdf[CDF_SIZE(DELTA_Q_PROBS + 1)];
  aom_cdf_prob delta_lf_multi_cdf[FRAME_LF_COUNT][CDF_SIZE(DELTA_LF_PROBS + 1)];
  aom_cdf_prob delta_lf_cdf[CDF_SIZE(DELTA_LF_PROBS + 1)];
  aom_cdf_prob intra_ext_tx_cdf[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES]
                               [CDF_SIZE(TX_TYPES)];
  aom_cdf_prob inter_ext_tx_cdf[EXT_TX_SETS_INTER][EXT_TX_SIZES]
                               [CDF_SIZE(TX_TYPES)];
  aom_cdf_prob cfl_sign_cdf[CDF_SIZE(CFL_JOINT_SIGNS)];
  aom_cdf_prob cfl_alpha_cdf[CFL_ALPHA_CONTEXTS][CDF_SIZE(CFL_ALPHABET_SIZE)];
#if CONFIG_FLEX_MVRES
  aom_cdf_prob
      sb_mv_precision_cdf[MV_SUBPEL_PRECISIONS - MV_SUBPEL_HALF_PRECISION]
                         [CDF_SIZE(MV_SUBPEL_PRECISIONS)];
#endif  // CONFIG_FLEX_MVRES

  int initialized;
} FRAME_CONTEXT;

static const int av1_ext_tx_ind[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 3, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 5, 6, 4, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 0 },
  { 3, 4, 5, 8, 6, 7, 9, 10, 11, 0, 1, 2, 0, 0, 0, 0 },
  { 7, 8, 9, 12, 10, 11, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6 },
};

static const int av1_ext_tx_inv[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 9, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 9, 0, 10, 11, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 9, 10, 11, 0, 1, 2, 4, 5, 3, 6, 7, 8, 0, 0, 0, 0 },
  { 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 4, 5, 3, 6, 7, 8 },
};

void av1_set_default_ref_deltas(int8_t *ref_deltas);
void av1_set_default_mode_deltas(int8_t *mode_deltas);
void av1_setup_frame_contexts(struct AV1Common *cm);
void av1_setup_past_independence(struct AV1Common *cm);

// Returns (int)ceil(log2(n)).
// NOTE: This implementation only works for n <= 2^30.
static INLINE int av1_ceil_log2(int n) {
  if (n < 2) return 0;
  int i = 1, p = 2;
  while (p < n) {
    i++;
    p = p << 1;
  }
  return i;
}

#if CONFIG_NEW_INTER_MODES
static INLINE int16_t inter_single_mode_ctx(int16_t mode_ctx) {
  // refmv_ctx values 2 and 4 are mapped to binary 1 while the rest map to 0.
  // This is intended to capture the case of ref_match_count >= 2 in
  // setup_ref_mv_list() function in mvref_common.c as a limited binary
  // context in addition to newmv_ctx and zeromv_ctx.
  // TODO(debargha, elliottk): Measure how much the limited refmv_ctx
  // actually helps
  static const int refmv_ctx_to_isrefmv_ctx[REFMV_MODE_CONTEXTS] = { 0, 0, 1,
                                                                     0, 1, 0 };
  const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
  assert(newmv_ctx < NEWMV_MODE_CONTEXTS);
  const int16_t zeromv_ctx = (mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
  const int16_t refmv_ctx = (mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;
  const int16_t isrefmv_ctx = refmv_ctx_to_isrefmv_ctx[refmv_ctx];
  const int16_t ctx =
      GLOBALMV_MODE_CONTEXTS * ISREFMV_MODE_CONTEXTS * newmv_ctx +
      ISREFMV_MODE_CONTEXTS * zeromv_ctx + isrefmv_ctx;
  assert(ctx < INTER_SINGLE_MODE_CONTEXTS);
  return ctx;
}

#if CONFIG_MVP_INDEPENDENT_PARSING
static INLINE int16_t drl_mode_ctx(int16_t mode_ctx) {
  const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
  assert(newmv_ctx < NEWMV_MODE_CONTEXTS);
  const int16_t ctx = (newmv_ctx >> 1);
  return ctx;
}
#else
static INLINE int16_t drl_mode_ctx(int16_t mode_ctx) {
  const int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
  assert(newmv_ctx < NEWMV_MODE_CONTEXTS);
  const int16_t zeromv_ctx = (mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
  const int16_t ctx = GLOBALMV_MODE_CONTEXTS * newmv_ctx + zeromv_ctx;
  assert(ctx < DRL_MODE_CONTEXTS);
  return ctx;
}
#endif  // CONFIG_MVP_INDEPENDENT_PARSING
#endif  // CONFIG_NEW_INTER_MODES

// Returns the context for palette color index at row 'r' and column 'c',
// along with the 'color_order' of neighbors and the 'color_idx'.
// The 'color_map' is a 2D array with the given 'stride'.
int av1_get_palette_color_index_context(const uint8_t *color_map, int stride,
                                        int r, int c, int palette_size,
                                        uint8_t *color_order, int *color_idx);

// A faster version of av1_get_palette_color_index_context used by the encoder
// exploiting the fact that the encoder does not need to maintain a color order.
int av1_fast_palette_color_index_context(const uint8_t *color_map, int stride,
                                         int r, int c, int *color_idx);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_ENTROPYMODE_H_
