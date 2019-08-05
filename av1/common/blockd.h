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

#ifndef AOM_AV1_COMMON_BLOCKD_H_
#define AOM_AV1_COMMON_BLOCKD_H_

#include "config/aom_config.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_ports/mem.h"
#include "aom_scale/yv12config.h"

#include "av1/common/common_data.h"
#include "av1/common/quant_common.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/mv.h"
#include "av1/common/scale.h"
#include "av1/common/seg_common.h"
#include "av1/common/tile_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define USE_B_QUANT_NO_TRELLIS 1

#define MAX_MB_PLANE 3

#define MAX_DIFFWTD_MASK_BITS 1

#define INTERINTRA_WEDGE_SIGN 0

// DIFFWTD_MASK_TYPES should not surpass 1 << MAX_DIFFWTD_MASK_BITS
enum {
  DIFFWTD_38 = 0,
  DIFFWTD_38_INV,
  DIFFWTD_MASK_TYPES,
} UENUM1BYTE(DIFFWTD_MASK_TYPE);

enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  INTRA_ONLY_FRAME = 2,  // replaces intra-only
  S_FRAME = 3,
  FRAME_TYPES,
} UENUM1BYTE(FRAME_TYPE);

static INLINE int is_comp_ref_allowed(BLOCK_SIZE bsize) {
  return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}

static INLINE int is_inter_mode(PREDICTION_MODE mode) {
  return mode >= INTER_MODE_START && mode < INTER_MODE_END;
}

typedef struct {
  uint8_t *plane[MAX_MB_PLANE];
  int stride[MAX_MB_PLANE];
} BUFFER_SET;

static INLINE int is_inter_singleref_mode(PREDICTION_MODE mode) {
  return mode >= SINGLE_INTER_MODE_START && mode < SINGLE_INTER_MODE_END;
}
static INLINE int is_inter_compound_mode(PREDICTION_MODE mode) {
  return mode >= COMP_INTER_MODE_START && mode < COMP_INTER_MODE_END;
}

static INLINE PREDICTION_MODE compound_ref0_mode(PREDICTION_MODE mode) {
  static PREDICTION_MODE lut[] = {
    MB_MODE_COUNT,  // DC_PRED
    MB_MODE_COUNT,  // V_PRED
    MB_MODE_COUNT,  // H_PRED
    MB_MODE_COUNT,  // D45_PRED
    MB_MODE_COUNT,  // D135_PRED
    MB_MODE_COUNT,  // D113_PRED
    MB_MODE_COUNT,  // D157_PRED
    MB_MODE_COUNT,  // D203_PRED
    MB_MODE_COUNT,  // D67_PRED
    MB_MODE_COUNT,  // SMOOTH_PRED
    MB_MODE_COUNT,  // SMOOTH_V_PRED
    MB_MODE_COUNT,  // SMOOTH_H_PRED
    MB_MODE_COUNT,  // PAETH_PRED
    MB_MODE_COUNT,  // NEARESTMV
    MB_MODE_COUNT,  // NEARMV
    MB_MODE_COUNT,  // GLOBALMV
    MB_MODE_COUNT,  // NEWMV
    NEARESTMV,      // NEAREST_NEARESTMV
    NEARMV,         // NEAR_NEARMV
    NEARESTMV,      // NEAREST_NEWMV
    NEWMV,          // NEW_NEARESTMV
    NEARMV,         // NEAR_NEWMV
    NEWMV,          // NEW_NEARMV
    GLOBALMV,       // GLOBAL_GLOBALMV
    NEWMV,          // NEW_NEWMV
  };
  assert(NELEMENTS(lut) == MB_MODE_COUNT);
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE PREDICTION_MODE compound_ref1_mode(PREDICTION_MODE mode) {
  static PREDICTION_MODE lut[] = {
    MB_MODE_COUNT,  // DC_PRED
    MB_MODE_COUNT,  // V_PRED
    MB_MODE_COUNT,  // H_PRED
    MB_MODE_COUNT,  // D45_PRED
    MB_MODE_COUNT,  // D135_PRED
    MB_MODE_COUNT,  // D113_PRED
    MB_MODE_COUNT,  // D157_PRED
    MB_MODE_COUNT,  // D203_PRED
    MB_MODE_COUNT,  // D67_PRED
    MB_MODE_COUNT,  // SMOOTH_PRED
    MB_MODE_COUNT,  // SMOOTH_V_PRED
    MB_MODE_COUNT,  // SMOOTH_H_PRED
    MB_MODE_COUNT,  // PAETH_PRED
    MB_MODE_COUNT,  // NEARESTMV
    MB_MODE_COUNT,  // NEARMV
    MB_MODE_COUNT,  // GLOBALMV
    MB_MODE_COUNT,  // NEWMV
    NEARESTMV,      // NEAREST_NEARESTMV
    NEARMV,         // NEAR_NEARMV
    NEWMV,          // NEAREST_NEWMV
    NEARESTMV,      // NEW_NEARESTMV
    NEWMV,          // NEAR_NEWMV
    NEARMV,         // NEW_NEARMV
    GLOBALMV,       // GLOBAL_GLOBALMV
    NEWMV,          // NEW_NEWMV
  };
  assert(NELEMENTS(lut) == MB_MODE_COUNT);
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE int have_nearmv_in_inter_mode(PREDICTION_MODE mode) {
  return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV ||
          mode == NEW_NEARMV);
}

static INLINE int have_newmv_in_inter_mode(PREDICTION_MODE mode) {
  return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
          mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE int is_masked_compound_type(COMPOUND_TYPE type) {
  return (type == COMPOUND_WEDGE || type == COMPOUND_DIFFWTD);
}

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

typedef struct {
  // Value of base colors for Y, U, and V
  uint16_t palette_colors[3 * PALETTE_MAX_SIZE];
  // Number of base colors for Y (0) and UV (1)
  uint8_t palette_size[2];
} PALETTE_MODE_INFO;

typedef struct {
  FILTER_INTRA_MODE filter_intra_mode;
  uint8_t use_filter_intra;
} FILTER_INTRA_MODE_INFO;

static const PREDICTION_MODE fimode_to_intradir[FILTER_INTRA_MODES] = {
  DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};

#if CONFIG_ADAPT_FILTER_INTRA
typedef struct {
  uint8_t use_adapt_filter_intra;
  ADAPT_FILTER_INTRA_MODE adapt_filter_intra_mode;
} ADAPT_FILTER_INTRA_MODE_INFO;

static const PREDICTION_MODE afimode_to_intradir[ADAPT_FILTER_INTRA_MODES] = {
  D135_PRED, D203_PRED, D67_PRED, D203_PRED, D67_PRED, D203_PRED, D67_PRED
};
#endif  // CONFIG_ADAPT_FILTER_INTRA

#if CONFIG_RD_DEBUG
#define TXB_COEFF_COST_MAP_SIZE (MAX_MIB_SIZE)
#endif

typedef struct RD_STATS {
  int rate;
  int64_t dist;
  // Please be careful of using rdcost, it's not guaranteed to be set all the
  // time.
  // TODO(angiebird): Create a set of functions to manipulate the RD_STATS. In
  // these functions, make sure rdcost is always up-to-date according to
  // rate/dist.
  int64_t rdcost;
  int64_t sse;
  int skip;  // sse should equal to dist when skip == 1
  int zero_rate;
#if CONFIG_RD_DEBUG
  int txb_coeff_cost[MAX_MB_PLANE];
  int txb_coeff_cost_map[MAX_MB_PLANE][TXB_COEFF_COST_MAP_SIZE]
                        [TXB_COEFF_COST_MAP_SIZE];
#endif  // CONFIG_RD_DEBUG
} RD_STATS;

// This struct is used to group function args that are commonly
// sent together in functions related to interinter compound modes
typedef struct {
  uint8_t *seg_mask;
  int8_t wedge_index;
  int8_t wedge_sign;
  DIFFWTD_MASK_TYPE mask_type;
  COMPOUND_TYPE type;
} INTERINTER_COMPOUND_DATA;

#define INTER_TX_SIZE_BUF_LEN 16
#define TXK_TYPE_BUF_LEN 64
// This structure now relates to 4x4 block regions.
typedef struct MB_MODE_INFO {
  // interinter members
  INTERINTER_COMPOUND_DATA interinter_comp;
  WarpedMotionParams wm_params;
  int_mv mv[2];
  int current_qindex;
  // Only for INTER blocks
  int_interpfilters interp_filters;
  // TODO(debargha): Consolidate these flags
#if CONFIG_RD_DEBUG
  RD_STATS rd_stats;
  int mi_row;
  int mi_col;
#endif
#if CONFIG_INSPECTION
  int16_t tx_skip[TXK_TYPE_BUF_LEN];
#endif
  PALETTE_MODE_INFO palette_mode_info;
  // Common for both INTER and INTRA blocks
  BLOCK_SIZE sb_type;
  PREDICTION_MODE mode;
  // Only for INTRA blocks
  UV_PREDICTION_MODE uv_mode;
  // interintra members
  INTERINTRA_MODE interintra_mode;
#if CONFIG_ADAPT_FILTER_INTRA
  ADAPT_FILTER_INTRA_MODE_INFO adapt_filter_intra_mode_info;
#endif
  MOTION_MODE motion_mode;
  PARTITION_TYPE partition;
  TX_TYPE txk_type[TXK_TYPE_BUF_LEN];
#if CONFIG_VQ4X4
  int gain_sign[MAX_MB_PLANE][TXK_TYPE_BUF_LEN];
  int qgain_idx[MAX_MB_PLANE][TXK_TYPE_BUF_LEN];
  int shape_idx[MAX_MB_PLANE][TXK_TYPE_BUF_LEN];
#endif
  MV_REFERENCE_FRAME ref_frame[2];
  FILTER_INTRA_MODE_INFO filter_intra_mode_info;
  int8_t skip;
  uint8_t inter_tx_size[INTER_TX_SIZE_BUF_LEN];
  TX_SIZE tx_size;
  int8_t delta_lf_from_base;
  int8_t delta_lf[FRAME_LF_COUNT];
  int8_t interintra_wedge_index;
  // The actual prediction angle is the base angle + (angle_delta * step).
  int8_t angle_delta[PLANE_TYPES];
  /* deringing gain *per-superblock* */
  // Joint sign of alpha Cb and alpha Cr
  int8_t cfl_alpha_signs;
  // Index of the alpha Cb and alpha Cr combination
  uint8_t cfl_alpha_idx;
  uint8_t num_proj_ref;
  uint8_t overlappable_neighbors[2];
  // If comp_group_idx=0, indicate if dist_wtd_comp(0) or avg_comp(1) is used.
  uint8_t compound_idx;
  uint8_t use_wedge_interintra : 1;
  uint8_t segment_id : 3;
  uint8_t seg_id_predicted : 1;  // valid only when temporal_update is enabled
  uint8_t skip_mode : 1;
  uint8_t use_intrabc : 1;
  uint8_t ref_mv_idx : 2;
  // Indicate if masked compound is used(1) or not(0).
  uint8_t comp_group_idx : 1;
  int8_t cdef_strength : 4;
#if CONFIG_INTRA_ENTROPY
  uint64_t gradient_hist[8];
  int64_t recon_var;  // Variance of reconstructed pixel values.
#endif                // CONFIG_INTRA_ENTROPY
} MB_MODE_INFO;

static INLINE int is_intrabc_block(const MB_MODE_INFO *mbmi) {
  return mbmi->use_intrabc;
}

static INLINE PREDICTION_MODE get_uv_mode(UV_PREDICTION_MODE mode) {
  assert(mode < UV_INTRA_MODES);
  static const PREDICTION_MODE uv2y[] = {
    DC_PRED,        // UV_DC_PRED
    V_PRED,         // UV_V_PRED
    H_PRED,         // UV_H_PRED
    D45_PRED,       // UV_D45_PRED
    D135_PRED,      // UV_D135_PRED
    D113_PRED,      // UV_D113_PRED
    D157_PRED,      // UV_D157_PRED
    D203_PRED,      // UV_D203_PRED
    D67_PRED,       // UV_D67_PRED
    SMOOTH_PRED,    // UV_SMOOTH_PRED
    SMOOTH_V_PRED,  // UV_SMOOTH_V_PRED
    SMOOTH_H_PRED,  // UV_SMOOTH_H_PRED
    PAETH_PRED,     // UV_PAETH_PRED
    DC_PRED,        // UV_CFL_PRED
    INTRA_INVALID,  // UV_INTRA_MODES
    INTRA_INVALID,  // UV_MODE_INVALID
  };
  return uv2y[mode];
}

static INLINE int is_inter_block(const MB_MODE_INFO *mbmi) {
  return is_intrabc_block(mbmi) || mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int has_second_ref(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[1] > INTRA_FRAME;
}

static INLINE int has_uni_comp_refs(const MB_MODE_INFO *mbmi) {
  return has_second_ref(mbmi) && (!((mbmi->ref_frame[0] >= BWDREF_FRAME) ^
                                    (mbmi->ref_frame[1] >= BWDREF_FRAME)));
}

static INLINE MV_REFERENCE_FRAME comp_ref0(int ref_idx) {
  static const MV_REFERENCE_FRAME lut[] = {
    LAST_FRAME,     // LAST_LAST2_FRAMES,
    LAST_FRAME,     // LAST_LAST3_FRAMES,
    LAST_FRAME,     // LAST_GOLDEN_FRAMES,
    BWDREF_FRAME,   // BWDREF_ALTREF_FRAMES,
    LAST2_FRAME,    // LAST2_LAST3_FRAMES
    LAST2_FRAME,    // LAST2_GOLDEN_FRAMES,
    LAST3_FRAME,    // LAST3_GOLDEN_FRAMES,
    BWDREF_FRAME,   // BWDREF_ALTREF2_FRAMES,
    ALTREF2_FRAME,  // ALTREF2_ALTREF_FRAMES,
  };
  assert(NELEMENTS(lut) == TOTAL_UNIDIR_COMP_REFS);
  return lut[ref_idx];
}

static INLINE MV_REFERENCE_FRAME comp_ref1(int ref_idx) {
  static const MV_REFERENCE_FRAME lut[] = {
    LAST2_FRAME,    // LAST_LAST2_FRAMES,
    LAST3_FRAME,    // LAST_LAST3_FRAMES,
    GOLDEN_FRAME,   // LAST_GOLDEN_FRAMES,
    ALTREF_FRAME,   // BWDREF_ALTREF_FRAMES,
    LAST3_FRAME,    // LAST2_LAST3_FRAMES
    GOLDEN_FRAME,   // LAST2_GOLDEN_FRAMES,
    GOLDEN_FRAME,   // LAST3_GOLDEN_FRAMES,
    ALTREF2_FRAME,  // BWDREF_ALTREF2_FRAMES,
    ALTREF_FRAME,   // ALTREF2_ALTREF_FRAMES,
  };
  assert(NELEMENTS(lut) == TOTAL_UNIDIR_COMP_REFS);
  return lut[ref_idx];
}

PREDICTION_MODE av1_left_block_mode(const MB_MODE_INFO *left_mi);

PREDICTION_MODE av1_above_block_mode(const MB_MODE_INFO *above_mi);

static INLINE int is_global_mv_block(const MB_MODE_INFO *const mbmi,
                                     TransformationType type) {
  const PREDICTION_MODE mode = mbmi->mode;
  const BLOCK_SIZE bsize = mbmi->sb_type;
  const int block_size_allowed =
      AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
  return (mode == GLOBALMV || mode == GLOBAL_GLOBALMV) && type > TRANSLATION &&
         block_size_allowed;
}

#if CONFIG_MISMATCH_DEBUG
static INLINE void mi_to_pixel_loc(int *pixel_c, int *pixel_r, int mi_col,
                                   int mi_row, int tx_blk_col, int tx_blk_row,
                                   int subsampling_x, int subsampling_y) {
  *pixel_c = ((mi_col >> subsampling_x) << MI_SIZE_LOG2) +
             (tx_blk_col << tx_size_wide_log2[0]);
  *pixel_r = ((mi_row >> subsampling_y) << MI_SIZE_LOG2) +
             (tx_blk_row << tx_size_high_log2[0]);
}
#endif

enum { MV_PRECISION_Q3, MV_PRECISION_Q4 } UENUM1BYTE(mv_precision);

struct buf_2d {
  uint8_t *buf;
  uint8_t *buf0;
  int width;
  int height;
  int stride;
};

typedef struct eob_info {
  uint16_t eob;
  uint16_t max_scan_line;
} eob_info;

typedef struct {
  DECLARE_ALIGNED(32, tran_low_t, dqcoeff[MAX_MB_PLANE][MAX_SB_SQUARE]);
  eob_info eob_data[MAX_MB_PLANE]
                   [MAX_SB_SQUARE / (TX_SIZE_W_MIN * TX_SIZE_H_MIN)];
  DECLARE_ALIGNED(16, uint8_t, color_index_map[2][MAX_SB_SQUARE]);
} CB_BUFFER;

typedef struct macroblockd_plane {
  tran_low_t *dqcoeff;
  tran_low_t *dqcoeff_block;
  eob_info *eob_data;
  PLANE_TYPE plane_type;
  int subsampling_x;
  int subsampling_y;
  struct buf_2d dst;
  struct buf_2d pre[2];
  ENTROPY_CONTEXT *above_context;
  ENTROPY_CONTEXT *left_context;

  // The dequantizers below are true dequantizers used only in the
  // dequantization process.  They have the same coefficient
  // shift/scale as TX.
  int16_t seg_dequant_QTX[MAX_SEGMENTS][2];
  uint8_t *color_index_map;

  // block size in pixels
  uint8_t width, height;

  qm_val_t *seg_iqmatrix[MAX_SEGMENTS][TX_SIZES_ALL];
  qm_val_t *seg_qmatrix[MAX_SEGMENTS][TX_SIZES_ALL];
} MACROBLOCKD_PLANE;

#define BLOCK_OFFSET(x, i) \
  ((x) + (i) * (1 << (tx_size_wide_log2[0] + tx_size_high_log2[0])))

typedef struct {
  DECLARE_ALIGNED(16, InterpKernel, vfilter);
  DECLARE_ALIGNED(16, InterpKernel, hfilter);
} WienerInfo;

typedef struct {
  int ep;
  int xqd[2];
} SgrprojInfo;

#if CONFIG_LOOP_RESTORE_CNN
typedef struct {
  FRAME_TYPE frame_type;
  int base_qindex;
} CNNInfo;
#endif  // CONFIG_LOOP_RESTORE_CNN

#if CONFIG_DEBUG
#define CFL_SUB8X8_VAL_MI_SIZE (4)
#define CFL_SUB8X8_VAL_MI_SQUARE \
  (CFL_SUB8X8_VAL_MI_SIZE * CFL_SUB8X8_VAL_MI_SIZE)
#endif  // CONFIG_DEBUG
#define CFL_MAX_BLOCK_SIZE (BLOCK_32X32)
#define CFL_BUF_LINE (32)
#define CFL_BUF_LINE_I128 (CFL_BUF_LINE >> 3)
#define CFL_BUF_LINE_I256 (CFL_BUF_LINE >> 4)
#define CFL_BUF_SQUARE (CFL_BUF_LINE * CFL_BUF_LINE)
typedef struct cfl_ctx {
  // Q3 reconstructed luma pixels (only Q2 is required, but Q3 is used to avoid
  // shifts)
  uint16_t recon_buf_q3[CFL_BUF_SQUARE];
  // Q3 AC contributions (reconstructed luma pixels - tx block avg)
  int16_t ac_buf_q3[CFL_BUF_SQUARE];

  // Cache the DC_PRED when performing RDO, so it does not have to be recomputed
  // for every scaling parameter
  int dc_pred_is_cached[CFL_PRED_PLANES];
  // The DC_PRED cache is disable when decoding
  int use_dc_pred_cache;
  // Only cache the first row of the DC_PRED
  int16_t dc_pred_cache[CFL_PRED_PLANES][CFL_BUF_LINE];

  // Height and width currently used in the CfL prediction buffer.
  int buf_height, buf_width;

  int are_parameters_computed;

  // Chroma subsampling
  int subsampling_x, subsampling_y;

  int mi_row, mi_col;

  // Whether the reconstructed luma pixels need to be stored
  int store_y;

#if CONFIG_DEBUG
  int rate;
#endif  // CONFIG_DEBUG

  int is_chroma_reference;
} CFL_CTX;

typedef struct dist_wtd_comp_params {
  int use_dist_wtd_comp_avg;
  int fwd_offset;
  int bck_offset;
} DIST_WTD_COMP_PARAMS;

struct scale_factors;

// Most/all of the pointers are mere pointers to actual arrays are allocated
// elsewhere. This is mostly for coding convenience.
typedef struct macroblockd {
  struct macroblockd_plane plane[MAX_MB_PLANE];

  TileInfo tile;

  int mi_stride;

  MB_MODE_INFO **mi;
  MB_MODE_INFO *left_mbmi;
  MB_MODE_INFO *above_mbmi;
  MB_MODE_INFO *chroma_left_mbmi;
  MB_MODE_INFO *chroma_above_mbmi;

  int up_available;
  int left_available;
  int chroma_up_available;
  int chroma_left_available;

  /* Distance of MB away from frame edges in subpixels (1/8th pixel)  */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  /* pointers to reference frame scale factors */
  const struct scale_factors *block_ref_scale_factors[2];

  /* pointer to current frame */
  const YV12_BUFFER_CONFIG *cur_buf;

  ENTROPY_CONTEXT *above_context[MAX_MB_PLANE];
  ENTROPY_CONTEXT left_context[MAX_MB_PLANE][MAX_MIB_SIZE];

  PARTITION_CONTEXT *above_seg_context;
  PARTITION_CONTEXT left_seg_context[MAX_MIB_SIZE];

  TXFM_CONTEXT *above_txfm_context;
  TXFM_CONTEXT *left_txfm_context;
  TXFM_CONTEXT left_txfm_context_buffer[MAX_MIB_SIZE];

  WienerInfo wiener_info[MAX_MB_PLANE];
  SgrprojInfo sgrproj_info[MAX_MB_PLANE];

  // block dimension in the unit of mode_info.
  uint8_t n4_w, n4_h;

  uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
  CANDIDATE_MV ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];
  uint16_t weight[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];
  uint8_t is_sec_rect;

  // Counts of each reference frame in the above and left neighboring blocks.
  // NOTE: Take into account both single and comp references.
  uint8_t neighbors_ref_counts[REF_FRAMES];

  FRAME_CONTEXT *tile_ctx;
  /* Bit depth: 8, 10, 12 */
  int bd;

  int qindex[MAX_SEGMENTS];
  int lossless[MAX_SEGMENTS];
  int corrupted;
  int cur_frame_force_integer_mv;
  // same with that in AV1_COMMON
  struct aom_internal_error_info *error_info;
  const WarpedMotionParams *global_motion;
  int delta_qindex;
  int current_qindex;
  // Since actual frame level loop filtering level value is not available
  // at the beginning of the tile (only available during actual filtering)
  // at encoder side.we record the delta_lf (against the frame level loop
  // filtering level) and code the delta between previous superblock's delta
  // lf and current delta lf. It is equivalent to the delta between previous
  // superblock's actual lf and current lf.
  int8_t delta_lf_from_base;
  // For this experiment, we have four frame filter levels for different plane
  // and direction. So, to support the per superblock update, we need to add
  // a few more params as below.
  // 0: delta loop filter level for y plane vertical
  // 1: delta loop filter level for y plane horizontal
  // 2: delta loop filter level for u plane
  // 3: delta loop filter level for v plane
  // To make it consistent with the reference to each filter level in segment,
  // we need to -1, since
  // SEG_LVL_ALT_LF_Y_V = 1;
  // SEG_LVL_ALT_LF_Y_H = 2;
  // SEG_LVL_ALT_LF_U   = 3;
  // SEG_LVL_ALT_LF_V   = 4;
  int8_t delta_lf[FRAME_LF_COUNT];
  int cdef_preset[4];

  DECLARE_ALIGNED(16, uint8_t, seg_mask[2 * MAX_SB_SQUARE]);
  uint8_t *mc_buf[2];
  CFL_CTX cfl;

  DIST_WTD_COMP_PARAMS jcp_param;

  uint16_t cb_offset[MAX_MB_PLANE];
  uint16_t txb_offset[MAX_MB_PLANE];
  uint16_t color_index_map_offset[2];

  CONV_BUF_TYPE *tmp_conv_dst;
  uint8_t *tmp_obmc_bufs[2];
} MACROBLOCKD;

static INLINE int is_cur_buf_hbd(const MACROBLOCKD *xd) {
  return xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH ? 1 : 0;
}

static INLINE uint8_t *get_buf_by_bd(const MACROBLOCKD *xd, uint8_t *buf16) {
  return (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH)
             ? CONVERT_TO_BYTEPTR(buf16)
             : buf16;
}

static INLINE int get_sqr_bsize_idx(BLOCK_SIZE bsize) {
  switch (bsize) {
    case BLOCK_4X4: return 0;
    case BLOCK_8X8: return 1;
    case BLOCK_16X16: return 2;
    case BLOCK_32X32: return 3;
    case BLOCK_64X64: return 4;
    case BLOCK_128X128: return 5;
    default: return SQR_BLOCK_SIZES;
  }
}

// For a square block size 'bsize', returns the size of the sub-blocks used by
// the given partition type. If the partition produces sub-blocks of different
// sizes, then the function returns the largest sub-block size.
// Implements the Partition_Subsize lookup table in the spec (Section 9.3.
// Conversion tables).
// Note: the input block size should be square.
// Otherwise it's considered invalid.
static INLINE BLOCK_SIZE get_partition_subsize(BLOCK_SIZE bsize,
                                               PARTITION_TYPE partition) {
  if (partition == PARTITION_INVALID) {
    return BLOCK_INVALID;
  } else {
    const int sqr_bsize_idx = get_sqr_bsize_idx(bsize);
    return sqr_bsize_idx >= SQR_BLOCK_SIZES
               ? BLOCK_INVALID
               : subsize_lookup[partition][sqr_bsize_idx];
  }
}

static TX_TYPE intra_mode_to_tx_type(const MB_MODE_INFO *mbmi,
                                     PLANE_TYPE plane_type) {
  static const TX_TYPE _intra_mode_to_tx_type[INTRA_MODES] = {
    DCT_DCT,    // DC_PRED
    ADST_DCT,   // V_PRED
    DCT_ADST,   // H_PRED
    DCT_DCT,    // D45_PRED
    ADST_ADST,  // D135_PRED
    ADST_DCT,   // D113_PRED
    DCT_ADST,   // D157_PRED
    DCT_ADST,   // D203_PRED
    ADST_DCT,   // D67_PRED
    ADST_ADST,  // SMOOTH_PRED
    ADST_DCT,   // SMOOTH_V_PRED
    DCT_ADST,   // SMOOTH_H_PRED
    ADST_ADST,  // PAETH_PRED
  };
  const PREDICTION_MODE mode =
      (plane_type == PLANE_TYPE_Y) ? mbmi->mode : get_uv_mode(mbmi->uv_mode);
  assert(mode < INTRA_MODES);
  return _intra_mode_to_tx_type[mode];
}

static INLINE int is_rect_tx(TX_SIZE tx_size) { return tx_size >= TX_SIZES; }

static INLINE int block_signals_txsize(BLOCK_SIZE bsize) {
  return bsize > BLOCK_4X4;
}

// Number of transform types in each set type
static const int av1_num_ext_tx_set[EXT_TX_SET_TYPES] = {
  1, 2, 5, 7, 12, 16,
#if CONFIG_VQ4X4
  0,  // not used
#endif
};

#if CONFIG_MODE_DEP_TX
// av1_num_ext_tx_set is used to indicate the number of symbols in
// inter_ext_tx_cdf, so we use 16 even when MDTXs are used
#if USE_MDTX_INTRA && USE_MDTX_INTER
static const int av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
#if CONFIG_VQ4X4
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
#endif
};
#elif USE_MDTX_INTRA
static const int av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 },
#if CONFIG_VQ4X4
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
#endif
};
#elif USE_MDTX_INTER
static const int av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
#if CONFIG_VQ4X4
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
#endif
};
#endif
#else
static const int av1_ext_tx_used[EXT_TX_SET_TYPES][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
#if CONFIG_VQ4X4
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
#endif
};
#endif

static const uint16_t av1_reduced_intra_tx_used_flag[INTRA_MODES] = {
  0x080F,  // DC_PRED:       0000 1000 0000 1111
  0x040F,  // V_PRED:        0000 0100 0000 1111
  0x080F,  // H_PRED:        0000 1000 0000 1111
  0x020F,  // D45_PRED:      0000 0010 0000 1111
  0x080F,  // D135_PRED:     0000 1000 0000 1111
  0x040F,  // D113_PRED:     0000 0100 0000 1111
  0x080F,  // D157_PRED:     0000 1000 0000 1111
  0x080F,  // D203_PRED:     0000 1000 0000 1111
  0x040F,  // D67_PRED:      0000 0100 0000 1111
  0x080F,  // SMOOTH_PRED:   0000 1000 0000 1111
  0x040F,  // SMOOTH_V_PRED: 0000 0100 0000 1111
  0x080F,  // SMOOTH_H_PRED: 0000 1000 0000 1111
  0x0C0E,  // PAETH_PRED:    0000 1100 0000 1110
};

static const uint16_t av1_ext_tx_used_flag[EXT_TX_SET_TYPES] = {
  0x0001,  // 0000 0000 0000 0001
  0x0201,  // 0000 0010 0000 0001
  0x020F,  // 0000 0010 0000 1111
  0x0E0F,  // 0000 1110 0000 1111
  0x0FFF,  // 0000 1111 1111 1111
  0xFFFF,  // 1111 1111 1111 1111
#if CONFIG_VQ4X4
  0x0000,  // not used
#endif
};

static INLINE TxSetType av1_get_ext_tx_set_type(TX_SIZE tx_size, int is_inter,
                                                int use_reduced_set) {
  const TX_SIZE tx_size_sqr_up = txsize_sqr_up_map[tx_size];
  if (tx_size_sqr_up > TX_32X32) return EXT_TX_SET_DCTONLY;
  if (tx_size_sqr_up == TX_32X32)
    return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DCTONLY;
  if (use_reduced_set)
    return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DTT4_IDTX;
  const TX_SIZE tx_size_sqr = txsize_sqr_map[tx_size];
#if CONFIG_VQ4X4
  if (!is_inter && tx_size_sqr_up == TX_4X4) return EXT_TX_SET_VQ;
#endif
  if (is_inter) {
    return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT9_IDTX_1DDCT
#if CONFIG_MODE_DEP_TX && USE_MDTX_INTER
                                    : EXT_TX_SET_ALL16_MDTX8);
#else
                                    : EXT_TX_SET_ALL16);
#endif
  } else {
    return (tx_size_sqr == TX_16X16 ? EXT_TX_SET_DTT4_IDTX
#if CONFIG_MODE_DEP_TX && USE_MDTX_INTRA
                                    : EXT_TX_SET_DTT4_IDTX_1DDCT_MDTX3);
#else
                                    : EXT_TX_SET_DTT4_IDTX_1DDCT);
#endif
  }
}

// Maps tx set types to the indices.
static const int ext_tx_set_index[2][EXT_TX_SET_TYPES] = {
  {
      // Intra
      0,
      -1,
      2,
      1,
      -1,
      -1,
#if CONFIG_VQ4X4
      3,
#endif
  },
  { // Inter
    0, 3, -1, -1, 2, 1 },
};

static INLINE int get_ext_tx_set(TX_SIZE tx_size, int is_inter,
                                 int use_reduced_set) {
  const TxSetType set_type =
      av1_get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
  return ext_tx_set_index[is_inter][set_type];
}

static INLINE int get_ext_tx_types(TX_SIZE tx_size, int is_inter,
                                   int use_reduced_set) {
  const int set_type =
      av1_get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
  return av1_num_ext_tx_set[set_type];
}

#define TXSIZEMAX(t1, t2) (tx_size_2d[(t1)] >= tx_size_2d[(t2)] ? (t1) : (t2))
#define TXSIZEMIN(t1, t2) (tx_size_2d[(t1)] <= tx_size_2d[(t2)] ? (t1) : (t2))

static INLINE TX_SIZE tx_size_from_tx_mode(BLOCK_SIZE bsize, TX_MODE tx_mode) {
  const TX_SIZE largest_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  const TX_SIZE max_rect_tx_size = max_txsize_rect_lookup[bsize];
  if (bsize == BLOCK_4X4)
    return AOMMIN(max_txsize_lookup[bsize], largest_tx_size);
  if (txsize_sqr_map[max_rect_tx_size] <= largest_tx_size)
    return max_rect_tx_size;
  else
    return largest_tx_size;
}

static const uint8_t mode_to_angle_map[] = {
  0, 90, 180, 45, 135, 113, 157, 203, 67, 0, 0, 0, 0,
};

// Converts block_index for given transform size to index of the block in raster
// order.
static INLINE int av1_block_index_to_raster_order(TX_SIZE tx_size,
                                                  int block_idx) {
  // For transform size 4x8, the possible block_idx values are 0 & 2, because
  // block_idx values are incremented in steps of size 'tx_width_unit x
  // tx_height_unit'. But, for this transform size, block_idx = 2 corresponds to
  // block number 1 in raster order, inside an 8x8 MI block.
  // For any other transform size, the two indices are equivalent.
  return (tx_size == TX_4X8 && block_idx == 2) ? 1 : block_idx;
}

// Inverse of above function.
// Note: only implemented for transform sizes 4x4, 4x8 and 8x4 right now.
static INLINE int av1_raster_order_to_block_index(TX_SIZE tx_size,
                                                  int raster_order) {
  assert(tx_size == TX_4X4 || tx_size == TX_4X8 || tx_size == TX_8X4);
  // We ensure that block indices are 0 & 2 if tx size is 4x8 or 8x4.
  return (tx_size == TX_4X4) ? raster_order : (raster_order > 0) ? 2 : 0;
}

static INLINE TX_TYPE get_default_tx_type(PLANE_TYPE plane_type,
                                          const MACROBLOCKD *xd,
                                          TX_SIZE tx_size,
                                          int is_screen_content_type) {
  const MB_MODE_INFO *const mbmi = xd->mi[0];

  if (is_inter_block(mbmi) || plane_type != PLANE_TYPE_Y ||
      xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32 ||
      is_screen_content_type)
    return DCT_DCT;

  return intra_mode_to_tx_type(mbmi, plane_type);
}

// Implements the get_plane_residual_size() function in the spec (Section
// 5.11.38. Get plane residual size function).
static INLINE BLOCK_SIZE get_plane_block_size(BLOCK_SIZE bsize,
                                              int subsampling_x,
                                              int subsampling_y) {
  if (bsize == BLOCK_INVALID) return BLOCK_INVALID;
  assert(subsampling_x >= 0 && subsampling_x < 2);
  assert(subsampling_y >= 0 && subsampling_y < 2);
  return ss_size_lookup[bsize][subsampling_x][subsampling_y];
}

static INLINE int av1_get_txb_size_index(BLOCK_SIZE bsize, int blk_row,
                                         int blk_col) {
  assert(bsize < BLOCK_SIZES_ALL);
  TX_SIZE txs = max_txsize_rect_lookup[bsize];
  for (int level = 0; level < MAX_VARTX_DEPTH - 1; ++level)
    txs = sub_tx_size_map[txs];
  const int tx_w_log2 = tx_size_wide_log2[txs] - MI_SIZE_LOG2;
  const int tx_h_log2 = tx_size_high_log2[txs] - MI_SIZE_LOG2;
  const int bw_log2 = mi_size_wide_log2[bsize];
  const int stride_log2 = bw_log2 - tx_w_log2;
  const int index =
      ((blk_row >> tx_h_log2) << stride_log2) + (blk_col >> tx_w_log2);
  assert(index < INTER_TX_SIZE_BUF_LEN);
  return index;
}

static INLINE int av1_get_txk_type_index(BLOCK_SIZE bsize, int blk_row,
                                         int blk_col) {
  assert(bsize < BLOCK_SIZES_ALL);
  TX_SIZE txs = max_txsize_rect_lookup[bsize];
  for (int level = 0; level < MAX_VARTX_DEPTH; ++level)
    txs = sub_tx_size_map[txs];
  const int tx_w_log2 = tx_size_wide_log2[txs] - MI_SIZE_LOG2;
  const int tx_h_log2 = tx_size_high_log2[txs] - MI_SIZE_LOG2;
  const int bw_uint_log2 = mi_size_wide_log2[bsize];
  const int stride_log2 = bw_uint_log2 - tx_w_log2;
  const int index =
      ((blk_row >> tx_h_log2) << stride_log2) + (blk_col >> tx_w_log2);
  assert(index < TXK_TYPE_BUF_LEN);
  return index;
}

#if CONFIG_VQ4X4
#define VQ_BLOCK_DEBUG 0
#define VQ_RD_DEBUG 0
#define VQ_GAIN_LEVELS 64
#define VQ_SHAPES 256

// Quantized gain (norm) values of 4x4 blocks
static const int16_t vq_gain_vals[VQ_GAIN_LEVELS] = {
  0,   10,  16,  23,  30,  38,  46,  54,  62,  70,  78,  86,  94,
  102, 110, 117, 125, 133, 141, 149, 157, 165, 173, 181, 189, 197,
  205, 213, 221, 229, 237, 245, 253, 261, 268, 276, 284, 292, 301,
  309, 317, 324, 332, 340, 347, 355, 363, 371, 379, 388, 396, 405,
  415, 423, 431, 439, 447, 453, 467, 485, 512, 551, 609, 729,
};

// Threshold for squared gain values
static const int32_t vq_gain_sq_thresholds[VQ_GAIN_LEVELS] = {
  0,      48,     183,    407,    733,    1179,   1757,   2473,
  3324,   4309,   5425,   6667,   8039,   9530,   11148,  12894,
  14755,  16741,  18852,  21117,  23534,  26055,  28684,  31432,
  34310,  37305,  40429,  43724,  47132,  50611,  54195,  57923,
  61816,  65832,  69967,  74195,  78604,  83205,  87916,  92778,
  97767,  102780, 107822, 112969, 118161, 123521, 129136, 134821,
  140787, 147064, 153602, 160461, 167880, 175225, 182117, 189006,
  196176, 202467, 211645, 226944, 248638, 282665, 337331, 451418,
};

// shape_4x4[i] = The i-th unit-norm codeword * 2^8
static const int32_t shape_4x4[VQ_SHAPES][16] = {
  { 66, 67, 63, 62, 67, 66, 63, 64, 64, 63, 61, 63, 64, 65, 63, 65 },
  { 29, 39, 39, 36, 43, 55, 59, 62, 51, 68, 79, 85, 55, 78, 92, 100 },
  { 28, 68, 77, 68, 31, 76, 84, 75, 26, 68, 76, 69, 26, 63, 72, 68 },
  { 12, 30, 62, 86, 16, 41, 77, 103, 21, 47, 75, 96, 27, 53, 74, 92 },
  { 21, 20, 17, 15, 48, 45, 37, 34, 87, 81, 64, 57, 108, 106, 85, 75 },
  { 4, 9, 24, 44, 12, 19, 48, 86, 17, 30, 74, 122, 19, 39, 93, 143 },
  { 33, 41, 31, 21, 64, 84, 71, 53, 69, 92, 78, 60, 63, 86, 71, 52 },
  { -4, -1, 16, 38, 9, 38, 72, 89, 39, 78, 95, 97, 55, 81, 83, 82 },
  { -10, -12, -12, -10, 2, 2, 2, 2, 69, 86, 82, 78, 86, 106, 104, 103 },
  { -2, 1, 60, 106, -2, -2, 67, 125, -2, -4, 56, 116, -1, -4, 48, 107 },
  { 78, 50, 22, 28, 101, 64, 31, 39, 103, 70, 40, 45, 98, 75, 48, 49 },
  { 74, 90, 83, 76, 79, 87, 78, 76, 47, 51, 48, 50, 25, 30, 33, 37 },
  { 5, 9, 21, 31, 14, 20, 36, 51, 41, 55, 68, 82, 79, 101, 113, 120 },
  { 5, 20, 55, 76, 50, 65, 74, 84, 76, 76, 64, 69, 70, 69, 58, 60 },
  { -3, -3, -5, -8, -9, -14, -17, -14, 1, 12, 34, 58, 61, 112, 141, 154 },
  { -13, 20, 88, 70, -18, 21, 105, 84, -18, 19, 100, 81, -17, 17, 96, 81 },
  { 0, 20, 32, 26, -1, 45, 72, 62, -5, 65, 106, 89, -4, 73, 118, 99 },
  { 3, -3, -1, 37, 8, -7, 20, 109, 6, -11, 43, 154, 1, -10, 55, 151 },
  { -2, -11, -13, -9, 52, 56, 54, 56, 87, 97, 96, 101, 60, 64, 63, 72 },
  { 1, -1, -8, -5, -3, -6, 6, 40, 2, 25, 78, 123, 37, 92, 125, 128 },
  { 21, 24, 23, 26, 74, 96, 95, 93, 64, 86, 92, 95, 12, 16, 23, 31 },
  { -1, -2, -12, -22, -2, -9, -16, 3, -7, -13, 51, 128, -10, 15, 127, 170 },
  { 0, 2, 3, 4, 2, 3, -5, -12, -2, -2, 6, 21, -7, 19, 143, 210 },
  { 44, 75, 85, 80, 30, 67, 94, 106, 3, 27, 65, 94, -14, -5, 27, 62 },
  { 14, 47, 44, 14, 21, 79, 99, 65, 11, 54, 103, 105, 1, 18, 68, 103 },
  { -1, -1, 0, 1, -5, -11, -10, -7, 51, 42, 29, 16, 131, 147, 120, 83 },
  { 17, 17, 10, 6, 62, 62, 23, 5, 103, 112, 36, -2, 110, 132, 51, 2 },
  { 4, -1, -13, -18, 28, 41, 29, 11, 33, 76, 109, 103, 12, 44, 104, 139 },
  { 98, 136, 139, 131, -8, -13, -13, -9, -6, -9, -9, -9, -2, -2, -1, -2 },
  { 49, 93, 36, 8, 56, 113, 43, 6, 52, 110, 44, 9, 53, 111, 54, 19 },
  { 19, 61, 117, 132, 17, 45, 75, 91, 18, 35, 46, 52, 23, 41, 52, 53 },
  { 1, 20, 13, -2, 10, 64, 59, 20, 25, 108, 108, 52, 31, 109, 115, 65 },
  { 56, 74, 69, 64, 72, 104, 113, 115, 17, 29, 40, 48, -14, -19, -19, -15 },
  { 5, 4, 5, 4, 0, -2, 1, 3, -10, -18, -19, -12, 114, 140, 129, 124 },
  { 0, -27, 4, 93, 1, -36, -2, 128, 2, -38, -8, 133, 5, -37, -7, 134 },
  { 0, 22, 63, 90, 22, 68, 115, 133, 33, 64, 76, 70, 14, 20, 14, 5 },
  { 37, 5, -2, 4, 105, 24, -4, 9, 149, 43, -3, 14, 157, 58, 2, 17 },
  { 47, 13, 26, 93, 62, 18, 30, 117, 67, 25, 31, 113, 62, 30, 33, 101 },
  { 63, 60, 7, -28, 72, 87, 54, 13, 51, 73, 85, 79, 40, 58, 78, 93 },
  { 30, 64, 78, 3, 38, 79, 95, 2, 45, 81, 96, 8, 52, 85, 97, 18 },
  { -3, -1, -9, -12, -1, -1, -14, 51, 7, -4, -14, 169, 10, -9, -10, 183 },
  { -3, 1, 2, -3, -2, 4, 0, -7, 5, 3, 10, 89, 12, -7, 44, 235 },
  { 56, 71, 69, 58, 81, 91, 77, 56, 86, 82, 45, 20, 72, 57, 9, -18 },
  { 71, 107, 123, 132, 33, 52, 65, 81, -16, -20, -21, -19, -8, -10, -14, -15 },
  { 2, -12, -15, 4, -18, -17, 30, 75, 3, 56, 104, 104, 77, 123, 95, 57 },
  { 82, 61, -11, -2, 98, 83, -12, -4, 94, 92, -8, -6, 101, 104, 0, -7 },
  { 8, 12, -5, 119, 5, 14, -13, 137, 2, 14, -20, 126, 0, 15, -16, 123 },
  { -4, -8, 4, 0, -12, -5, 57, 52, -27, 5, 128, 106, -30, 11, 136, 106 },
  { 18, 52, 111, 131, 8, 28, 94, 144, 2, 0, 16, 48, 4, 1, -3, 6 },
  { 1, 73, 69, 39, -9, 92, 85, 47, -20, 94, 87, 40, -17, 94, 85, 34 },
  { 7, 44, 153, 199, 0, -3, -10, -4, 4, 2, -5, -11, 2, 3, 4, 7 },
  { 33, 4, -10, 4, 79, 63, 6, -7, 72, 115, 78, 25, 27, 85, 116, 86 },
  { -13, -13, -1, 1, 12, -6, -7, 1, 124, 65, -18, -5, 173, 122, -10, -7 },
  { 122, 10, 5, 10, 139, 3, 8, 12, 130, 0, 12, 10, 115, 4, 14, 11 },
  { -23, -35, -6, 47, -40, -63, -6, 86, -55, -84, -5, 110, -64, -91, 0, 118 },
  { 76, 45, 27, 62, 74, 71, 56, 77, 22, 59, 80, 86, -3, 33, 78, 93 },
  { -12, -6, 3, -2, 54, -14, -1, -1, 168, -18, -10, 1, 182, -18, -15, 3 },
  { 0, 0, 1, 3, 1, -3, -8, -10, -15, -27, -37, -34, 3, 82, 168, 164 },
  { 41, -24, -54, -35, 74, -25, -79, -50, 98, -20, -91, -58, 103, -12, -92,
    -61 },
  { -12, -18, -22, -18, 7, 34, 76, 96, 14, 62, 127, 152, 0, 19, 40, 53 },
  { 33, 53, 58, 50, 18, 32, 36, 34, 6, 26, 43, 48, 39, 93, 136, 135 },
  { -23, -45, -63, -71, -16, -33, -44, -41, 50, 67, 75, 81, 63, 86, 93, 98 },
  { 58, 69, 40, -33, 77, 87, 52, -44, 78, 86, 54, -42, 80, 88, 56, -37 },
  { 1, 1, 3, 6, 1, 4, -3, -21, -1, -7, -32, -34, -9, -36, 33, 246 },
  { 56, 69, 76, 81, 37, 36, 40, 57, 44, 37, 21, 28, 99, 111, 83, 67 },
  { -38, -20, 26, 56, -48, -23, 45, 90, -44, -17, 61, 118, -39, -8, 78, 139 },
  { 107, 130, 101, 67, 85, 94, 68, 40, 10, 6, 1, -3, -4, -8, -9, -8 },
  { 96, 7, -32, -3, 122, 4, -41, -3, 129, 3, -44, -4, 134, 8, -45, -6 },
  { -4, -4, 43, 66, -9, 29, 90, 75, 24, 91, 98, 33, 75, 120, 70, 4 },
  { 102, 119, 124, 122, 32, 35, 38, 42, 15, 18, 21, 23, 22, 29, 33, 35 },
  { -8, -9, -9, -12, -55, -70, -69, -66, -18, -14, -11, -10, 75, 115, 120,
    118 },
  { 40, 18, -18, -31, 86, 46, -28, -48, 119, 67, -34, -57, 128, 74, -34, -58 },
  { -5, -6, -14, -15, -13, 2, 2, -10, -22, 56, 107, 56, -17, 91, 165, 104 },
  { 3, 22, 15, -25, -5, 22, 69, 41, -13, -7, 97, 145, -8, -26, 52, 152 },
  { 0, 2, 10, 19, 89, 96, 89, 83, -2, 2, 13, 26, 89, 96, 89, 85 },
  { 12, -12, 2, 79, 3, -14, 53, 124, -10, 20, 109, 96, -8, 65, 117, 38 },
  { 55, 6, 50, 33, 83, 12, 70, 50, 103, 17, 80, 59, 111, 23, 83, 62 },
  { 0, -5, -60, 55, -2, -9, -87, 86, -4, -9, -102, 102, -7, -12, -106, 107 },
  { 90, 59, -2, -4, 150, 83, -6, -5, 128, 59, -8, -1, 63, 24, -5, -3 },
  { 8, 45, 91, 104, 22, 64, 102, 100, 4, 2, -14, -27, -21, -51, -83, -92 },
  { -20, 7, 126, 179, -18, -22, 39, 112, -6, -18, -29, -19, -4, -9, -14, -17 },
  { -2, -5, -5, 2, 63, 89, 101, 105, 56, 79, 89, 94, -24, -41, -47, -41 },
  { 4, -5, 27, 80, 5, -10, 62, 192, 6, -6, 33, 128, -2, -3, -5, 11 },
  { -9, -13, 96, 7, -9, -19, 130, 4, -8, -20, 137, 2, -7, -14, 138, 3 },
  { -13, 41, 68, -6, -25, 49, 111, 3, -29, 43, 133, 24, -27, 36, 137, 45 },
  { 1, 102, 80, -18, -6, 112, 84, -31, -10, 98, 71, -36, -7, 82, 63, -29 },
  { 1, 0, -2, 0, 1, 22, 13, -3, 17, 116, 76, -8, 32, 174, 119, -4 },
  { -9, 18, 36, -16, -14, 56, 86, -25, -14, 86, 122, -36, -11, 97, 132, -36 },
  { -37, -34, 44, 39, -55, -51, 71, 68, -64, -62, 83, 83, -69, -67, 79, 85 },
  { 68, 64, 7, -6, 54, 110, 54, -4, -17, 77, 107, 33, -38, 4, 103, 90 },
  { 0, 2, 0, -1, -3, -16, -14, -2, -7, 3, 3, -6, 1, 181, 179, 4 },
  { 82, -46, 1, 3, 119, -65, 1, 7, 125, -69, 5, 8, 117, -61, 7, 9 },
  { -6, -8, -9, -7, 8, 6, 6, 12, 109, 139, 135, 124, 2, -1, -3, -4 },
  { 5, 4, 3, 4, 12, 5, -9, -6, 49, 63, 32, 9, 4, 114, 173, 121 },
  { -1, 1, 1, 2, -13, -7, -1, 4, 20, 45, 5, -5, 131, 211, 27, -18 },
  { 86, 116, 111, 94, 27, 26, 15, 7, -50, -77, -78, -67, -20, -31, -28, -20 },
  { 19, 62, -27, -7, 27, 102, -37, -12, 31, 144, -39, -18, 32, 147, -38, -21 },
  { 133, 104, 26, 12, 133, 95, 25, 20, 52, 42, 27, 28, 18, 23, 27, 29 },
  { -70, -83, -65, -48, 19, 30, 36, 36, 65, 85, 85, 83, 52, 68, 71, 73 },
  { -12, 61, 19, -45, -21, 92, 30, -70, -23, 108, 37, -87, -17, 116, 42, -90 },
  { 43, 67, 33, 0, 94, 152, 76, 5, 65, 110, 62, 12, -3, 0, 4, 3 },
  { 3, 87, 115, 31, -2, 99, 149, 44, -5, 42, 79, 29, -6, 4, 25, 11 },
  { 7, 5, -5, -5, -23, -9, -3, -4, -14, -12, -10, -6, 226, 112, -23, -8 },
  { -4, 7, 32, 29, -10, 31, 126, 132, -18, 10, 99, 132, -19, -28, -19, 7 },
  { 62, 76, 79, 75, 18, 18, 21, 31, 79, 95, 98, 95, 38, 43, 46, 51 },
  { -46, -28, 66, 123, -59, -53, 55, 128, -54, -52, 23, 78, -39, -42, 6, 43 },
  { 36, 24, 7, 2, 59, 36, 3, 0, 74, 44, 21, 57, 84, 59, 81, 175 },
  { 36, 26, -41, -94, 64, 62, -3, -67, 79, 84, 49, 12, 84, 96, 73, 54 },
  { 6, -10, 22, 142, 9, -11, 8, 199, 4, -2, -20, 47, 0, 1, -19, -40 },
  { 56, 83, 82, 68, 23, 42, 60, 69, -7, -10, 35, 84, 27, 14, 67, 141 },
  { 47, -30, -37, 56, 70, -38, -56, 79, 84, -32, -71, 87, 92, -23, -74, 84 },
  { 71, 73, 56, 47, 88, 106, 94, 79, -26, -33, -28, -19, -49, -64, -60, -52 },
  { 3, 39, 102, 131, 76, 115, 97, 62, 47, 41, 7, -9, -4, -9, -6, 4 },
  { 75, 101, 13, -53, 83, 115, 8, -75, 52, 77, -6, -73, 25, 41, -17, -65 },
  { 10, -29, 2, 247, -3, -5, -29, -34, -3, -4, -7, -37, -1, -4, -4, -1 },
  { 67, 88, 46, -15, 42, 65, 96, 106, 24, 33, 64, 105, 31, 45, 45, 50 },
  { 79, 98, 99, 85, -3, -9, -13, -10, 77, 96, 95, 83, -7, -15, -22, -20 },
  { -13, -18, -25, -34, -28, -44, -70, -87, -39, -51, -32, -10, -36, -6, 114,
    169 },
  { 84, 63, -36, -75, 26, -9, -67, -84, -50, -66, -71, -79, -50, -60, -70,
    -80 },
  { 54, 89, 112, 115, -50, -82, -102, -96, 1, 1, -3, -4, -2, -2, -3, 0 },
  { 26, 45, -14, 53, 37, 65, -16, 90, 44, 80, -17, 113, 51, 90, -11, 118 },
  { -6, -6, -8, -10, 8, 12, 14, 19, 59, 95, 105, 98, -64, -97, -100, -87 },
  { -56, -60, -52, -57, -73, -82, -72, -75, -49, -59, -52, -49, 41, 71, 83,
    74 },
  { -13, 0, 0, -3, -25, -5, 6, -1, 61, -24, 4, 4, 243, -30, -14, 12 },
  { -34, -52, -35, -15, 2, -12, -19, -20, 92, 114, 52, -6, 101, 136, 79, 4 },
  { 73, -1, -58, -41, 93, -18, -81, -56, 39, -56, -85, -65, -23, -83, -86,
    -75 },
  { -79, -91, -83, -83, -27, -28, -25, -31, 19, 19, 19, 23, 70, 85, 95, 107 },
  { 141, -11, -17, 4, 192, -21, -16, 6, 80, -21, -4, 2, -25, -14, 1, -3 },
  { 160, 188, 61, -4, -9, -9, -3, 4, -14, -15, -4, 1, -1, -2, -4, -4 },
  { 3, -5, -9, 0, 110, 135, 133, 131, 7, 6, 4, 9, 0, 0, -3, -3 },
  { 70, -35, -48, 14, 95, -60, -74, 16, 92, -75, -84, 13, 76, -72, -80, 4 },
  { 37, 96, 98, 53, -25, 20, 104, 126, -41, -59, -6, 69, -7, -34, -53, -29 },
  { 6, 19, 23, 19, -8, 24, 41, 36, 53, 37, 41, 53, 193, 95, 33, 67 },
  { -51, -92, -42, 38, -56, -69, 12, 85, 6, 27, 69, 97, 48, 68, 78, 90 },
  { -35, -53, -58, -45, 62, 82, 83, 77, -47, -69, -74, -63, 53, 64, 70, 68 },
  { 57, 83, 88, 79, -48, -60, -62, -60, -63, -85, -87, -86, 13, 22, 27, 21 },
  { -50, -70, -75, -62, -84, -84, -39, 2, -50, 0, 63, 85, 11, 62, 89, 86 },
  { 24, 66, 70, 54, -15, 43, 89, 85, -65, -40, 57, 94, -61, -72, 25, 91 },
  { 46, 86, 100, 103, -11, -18, -22, -18, -19, -27, -38, -30, 48, 83, 103,
    104 },
  { 16, -18, -43, -25, 63, 41, -53, -79, 74, 112, 15, -89, 39, 109, 85, -31 },
  { 13, -19, -5, 30, 49, -20, -3, 64, 105, -7, 6, 108, 136, 13, 15, 125 },
  { 76, 77, -2, -29, 119, 65, -42, -46, 102, -1, -75, -47, 54, -43, -77, -42 },
  { -27, -48, -40, -21, -8, -58, -87, -69, 69, 38, -34, -60, 94, 124, 90, 43 },
  { 85, 84, 16, -11, 32, 110, 101, 31, -26, 5, 93, 103, -8, -39, -16, 69 },
  { 27, 28, 40, 38, 77, 57, 49, 56, 49, 81, 70, 63, -101, 38, 105, 74 },
  { -25, -50, -64, -43, -54, -79, -44, 11, -75, -73, 38, 100, -63, -35, 78,
    108 },
  { 83, 26, -19, 6, 97, 92, -5, -6, 30, 127, 50, -21, -27, 97, 97, -5 },
  { 18, -46, 45, 58, 28, -64, 63, 80, 35, -75, 71, 93, 35, -74, 73, 100 },
  { -1, -3, -23, -25, 6, -7, -11, 159, 6, -8, -6, 194, 2, -2, -16, 29 },
  { 91, -44, -79, -31, -4, -81, -67, -38, -79, -77, -54, -56, -73, -62, -61,
    -67 },
  { -11, -22, -21, -14, 99, 98, 44, 20, 127, 132, 69, 39, 42, 31, 11, 8 },
  { 108, 163, 45, -27, 76, 120, 32, -27, -16, -16, -10, -14, -23, -27, -18,
    -16 },
  { 15, -10, 30, 221, 12, 17, 23, 109, 10, 24, 22, 12, 6, 11, 22, 28 },
  { -51, -32, 60, -12, -74, -50, 84, -18, -85, -55, 98, -23, -96, -60, 94,
    -31 },
  { 160, 66, -34, -13, 159, 77, -28, -18, 5, 17, -7, -17, -29, -6, -5, -17 },
  { -11, 63, 35, -9, 46, 105, 22, -17, 123, 94, -13, -6, 139, 47, -26, 8 },
  { -8, -4, 1, -8, 25, 78, 54, -11, 62, 160, 85, -31, 47, 110, 34, -51 },
  { 0, 3, 7, 12, 1, 64, 181, 168, -4, -8, -6, -3, 2, -1, -6, -5 },
  { 228, 104, -22, 4, 31, 4, -9, 0, -23, -17, -7, -5, -1, -3, -7, -6 },
  { 14, -5, -6, 61, -12, -49, -6, 104, -33, -36, 49, 114, -9, 52, 125, 109 },
  { -24, -4, 70, 81, 1, 71, 109, 56, 79, 112, 41, -9, 96, 55, -20, -18 },
  { 4, 8, -12, -24, 7, -21, -52, -2, -14, -82, -43, 106, -44, -83, 51, 174 },
  { 36, -17, -84, -64, 55, 64, -34, -102, 12, 71, 90, 12, -9, 17, 87, 113 },
  { -4, -3, 1, 0, -12, -17, -7, 3, 52, 3, -61, -36, 145, 63, -132, -122 },
  { -27, -53, 29, 84, -21, -94, -15, 126, 2, -79, -79, 82, 9, -42, -87, 11 },
  { -16, -25, 20, 59, -22, -24, 63, 137, -27, -34, 22, 57, -39, -78, -106,
    -112 },
  { -32, -32, -24, -18, -65, -82, -65, -45, 14, -32, -75, -73, 126, 116, -10,
    -62 },
  { 2, 1, -3, -7, -13, -17, -6, -1, -35, -50, -32, -18, 131, 189, 78, -30 },
  { 31, 207, 143, -11, 2, 15, 5, 3, -11, -19, -13, 6, 0, 0, -3, 2 },
  { -14, -6, -6, -14, -21, -13, -11, -33, -36, -8, 66, -24, -54, 34, 229, 22 },
  { -2, -3, -4, -6, -28, -41, -44, -37, 55, 94, 134, 152, 20, 27, 48, 58 },
  { -8, 41, 2, 8, -35, 83, 4, 8, -84, 123, 8, 4, -107, 147, 16, -1 },
  { -11, 107, 168, 93, -9, 61, 99, 56, -2, -4, -9, -7, 0, -3, -10, -5 },
  { -23, 12, 100, 133, -37, -52, -51, -47, -29, -54, -83, -98, -30, -43, -51,
    -54 },
  { -2, 1, -1, -41, -2, 10, -18, -54, 3, 11, -56, 65, 15, -1, -68, 220 },
  { 47, -13, 9, -46, 83, -22, 13, -72, 107, -25, 18, -97, 113, -25, 16, -112 },
  { 2, 11, 88, 171, 7, 5, 5, 14, 26, 35, 42, 42, 31, 55, 91, 102 },
  { -13, 56, 96, 5, -10, 8, 99, 98, 43, 9, 42, 125, 61, 55, 33, 79 },
  { -14, -37, -26, 24, -16, 8, 88, 139, 53, 101, 105, 75, 53, 50, 4, -21 },
  { -57, 85, 222, 70, -15, -10, 2, -15, -2, -7, -12, -9, -2, -2, 1, -1 },
  { -26, 117, 3, 21, -28, 136, -2, 25, -18, 130, -4, 30, -16, 107, 0, 29 },
  { -40, -83, 10, 125, -23, -84, -63, 76, -8, -39, -88, -71, -4, -16, -48,
    -85 },
  { 19, 64, 117, 128, -14, -8, 19, 43, -18, -45, -92, -119, 4, 7, -20, -55 },
  { 96, 164, 86, 3, 28, 39, 25, 14, 12, 24, 37, 44, 39, 56, 71, 71 },
  { -34, -16, 0, -1, 103, -6, -14, 3, 228, 3, -25, 6, 8, -22, -9, -1 },
  { -36, 35, 71, -8, -56, -17, 103, 57, -23, -67, 47, 127, 10, -58, -39, 111 },
  { -65, -45, -3, 0, -52, -80, -15, 10, 87, -42, -46, 19, 177, 58, -56, 9 },
  { 52, 67, 78, 87, 58, 71, 75, 79, -42, -62, -68, -51, 30, 46, 61, 68 },
  { -40, -60, -34, -6, 18, -16, -59, -60, 58, 93, 58, -4, -5, 42, 121, 140 },
  { 68, 3, -25, 6, 18, 107, 15, -32, -33, 20, 136, 34, 9, -35, 22, 155 },
  { -13, -13, -3, 5, 6, 8, -5, -2, 180, 179, 23, -11, 10, 11, -5, -5 },
  { 90, 82, -5, -14, 128, 57, -34, -7, 13, -64, -50, 10, -90, -113, -41, 9 },
  { -6, -21, -32, -24, 119, 119, -4, -47, 103, 109, -3, -47, -40, -51, -48,
    -39 },
  { -1, -4, -12, -18, -2, -11, -4, 14, -3, -24, 75, 241, -7, -13, -8, -2 },
  { 45, 47, 22, 39, 95, 63, 18, 61, 128, 28, -19, 70, 110, -19, -40, 75 },
  { -51, -55, -66, -67, -52, -59, -77, -42, -32, -60, -73, 46, -17, -60, -58,
    134 },
  { -25, -34, -42, -49, -39, -34, -45, -79, -24, 56, 30, -90, -11, 133, 93,
    -97 },
  { 243, -50, -33, 15, 15, -32, -1, -4, -34, -9, -1, -10, 0, -4, -6, -8 },
  { 0, 117, 105, -30, -25, 116, 126, -37, -38, 4, 22, -38, -29, -44, -37, -36 },
  { -22, -4, 12, -11, -14, 106, 144, 21, -3, 112, 138, 20, -15, -5, 2, -9 },
  { -111, -83, 103, 179, -22, -30, -23, 3, 10, 4, -22, -34, 2, 3, 4, 3 },
  { -11, -46, 6, 42, 87, -10, -17, 43, 131, 96, 21, 28, 35, 110, 92, 47 },
  { 1, -2, -8, -12, -6, 4, 22, 5, -25, 92, 209, 108, -7, -7, 3, -4 },
  { 20, 51, 99, 122, 1, 25, 59, 87, -41, -16, 15, 35, -115, -94, -26, 8 },
  { -33, -19, 7, 4, 8, -53, -31, 0, 119, 48, -60, -35, 36, 178, 47, -62 },
  { 60, -43, -27, 18, 123, -10, -58, 10, 117, 71, -55, -19, 54, 116, -15, -42 },
  { 6, -16, 12, 188, 13, 1, -9, 81, 1, 15, -17, -109, -6, 6, -6, -101 },
  { -31, -47, -48, -39, 59, 85, 94, 82, -48, -65, -63, -65, -46, -63, -70,
    -80 },
  { -73, 27, 59, -29, -66, 88, 53, -54, 26, 117, 7, -45, 104, 97, -23, -20 },
  { 2, -20, -80, 177, 0, -14, -75, 136, -3, -3, -32, 42, -1, -2, -10, 9 },
  { -1, -2, 9, 5, -3, 4, 10, -56, -3, 7, -43, -113, 13, -32, -109, 186 },
  { -15, -36, -34, 14, -40, -107, -93, 21, -17, -83, -58, 69, 39, 39, 67, 133 },
  { -8, -27, -50, -64, 57, 108, 96, 53, -37, -13, 55, 106, -33, -64, -77, -66 },
  { -3, 0, -3, 1, -48, 2, 4, 0, -94, -40, 8, 3, 206, -93, -38, 15 },
  { 2, -4, -13, -15, 14, 71, 117, 104, 10, 2, -11, -14, 20, 78, 126, 113 },
  { -9, 1, -9, -35, -19, 0, -7, -43, -61, -59, -2, 50, -95, -136, -18, 156 },
  { -14, -13, 61, 112, -32, 1, 119, 123, -28, 32, 82, -6, -10, 33, -1, -95 },
  { 9, -24, 10, 36, -18, -30, 44, 21, -68, 40, 115, -31, -78, 121, 123, -60 },
  { 7, 2, -22, -33, 3, -1, 5, 19, 33, 66, 111, 113, 111, 91, -42, -103 },
  { 5, 14, -2, -12, 84, 214, 108, -18, -5, 8, 3, -8, -1, 1, 2, -1 },
  { -22, 117, 98, -11, 105, 139, 14, -31, 87, 20, -32, -6, -7, -31, -6, 18 },
  { -86, -13, 121, 49, -95, 4, 152, 48, -45, -11, 27, -17, -21, -26, -32, -36 },
  { 19, -39, 43, -15, 41, -60, 68, -18, 67, -82, 106, -13, 76, -88, 122, -17 },
  { -15, 3, -9, -7, 241, 82, -15, 5, 4, 3, -1, -1, -5, -2, 3, 3 },
  { -15, -23, -13, 4, 50, 55, 34, 15, 62, 90, 102, 95, -40, -86, 12, 131 },
  { 83, 106, 70, 32, 50, 41, -18, -52, 38, 17, -54, -89, 99, 103, 35, -14 },
  { 139, 114, -33, -95, -83, -94, -61, -39, -31, -26, -11, -12, 1, 0, -8, -12 },
  { 97, 62, 37, 129, 105, 46, -7, 105, 68, 30, -34, 20, 38, 22, -25, -14 },
  { 113, 133, 21, -46, 68, 17, -83, -96, -25, -50, -41, -15, 0, 19, 49, 61 },
  { -100, -71, -19, -24, -20, -14, -27, -40, 132, 69, -41, -51, 118, 34, -57,
    -42 },
  { -52, -46, 52, 104, -32, -90, -8, 123, 48, 8, -27, 35, 72, 105, 48, 8 },
  { 147, 22, -2, 79, 117, -5, 19, 105, 23, -8, 41, 79, -3, 11, 39, 45 },
  { 75, 43, -26, -3, 108, 13, -49, 35, 74, -39, -22, 102, 29, -43, 51, 137 },
  { -30, 19, 66, -94, -25, 19, 85, -113, 5, 17, 80, 9, 21, 18, 65, 135 },
  { 142, -59, 3, 11, 108, -57, 19, 11, -71, 0, 16, -4, -137, 48, 11, -4 },
  { 31, 89, 4, -47, -27, 49, 140, 22, 1, -27, 41, 167, 7, 2, -28, 11 },
  { 58, -12, -106, -33, 41, -60, -81, 46, -13, -14, 59, 100, -1, 67, 112, 69 },
  { -14, 11, 102, 136, -5, -14, -6, -7, -4, 17, 114, 148, -6, -21, -16, -10 },
  { 107, 114, -4, -101, -4, 28, 102, 137, -8, -12, 1, 23, 4, 6, -2, -8 },
  { 11, 2, -9, -16, 14, 14, 8, -1, -103, -107, -35, 14, 125, 81, -69, -119 },
  { 195, 50, -15, 16, -37, -5, 35, 42, 10, 39, 51, 55, 58, 69, 52, 46 },
  { -10, -8, -2, -1, -2, -3, 61, 160, -5, -4, -9, -17, -6, -11, 56, 180 },
  { -10, -22, -44, -36, -13, -35, 55, 236, -8, -16, -22, -4, -1, -5, -25, -14 },
  { -44, 55, 85, -46, -25, -45, 99, 117, 12, -35, -48, 124, 0, 3, -45, -62 },
  { -44, 4, 110, 138, 120, 92, -17, -53, -63, -47, -8, 8, 7, 6, -1, -7 },
  { 134, 102, 13, -22, -25, -27, -19, -13, 141, 112, 20, -20, -20, -24, -19,
    -14 },
  { -17, -31, 3, 8, 104, -10, -43, 6, 2, 137, 5, -52, -55, -2, 163, 18 },
  { -66, 150, 123, -117, -38, -25, -24, -41, -20, -37, -37, -25, -18, -19, -20,
    -24 },
  { -4, 21, 52, 63, 77, 54, -2, -16, -33, 0, 84, 128, 117, 100, 16, -28 },
  { 55, 63, 45, 28, 34, 74, 104, 97, 3, -4, 18, 49, 122, 77, -48, -54 },
  { 18, 62, 96, 96, -45, -37, 20, 69, -4, -54, -59, -13, 128, 91, -26, -49 },
  { -25, 40, 150, -31, 16, -24, 5, 196, 0, 10, -19, -9, -1, 4, 8, -11 },
  { -10, -19, -3, 33, -5, 27, 83, 99, 108, 86, -31, -86, -77, -87, -63, -41 },
  { 86, 93, 47, 17, -95, -73, 52, 115, 52, 25, -55, -90, -8, 0, 27, 37 },
  { 84, -31, 111, 29, 83, -43, 130, 1, 45, -41, 95, -18, -7, -56, 30, -41 },
  { 32, -105, 112, -46, 62, -81, 87, -18, 35, 5, -20, 57, 7, 69, -74, 72 },
};

static INLINE void update_cw_array(int *sign_arr, int *qgain_idx_arr,
                                   int *shape_idx_arr, BLOCK_SIZE bsize,
                                   int blk_row, int blk_col, int sign,
                                   int gain_idx, int shape_idx) {
  const int blk_idx = av1_get_txk_type_index(bsize, blk_row, blk_col);
  sign_arr[blk_idx] = sign;
  qgain_idx_arr[blk_idx] = gain_idx;
  shape_idx_arr[blk_idx] = shape_idx;
}

static INLINE void vq_qgain_idx_to_symbols(int qgain_idx, int *gain_sym1,
                                           int *gain_sym2) {
  assert(qgain_idx >= 0 && qgain_idx < VQ_GAIN_LEVELS);
  *gain_sym1 = qgain_idx / VQ_GAIN_SYMBOLS_2;
  *gain_sym2 = qgain_idx % VQ_GAIN_SYMBOLS_2;
}

static INLINE void vq_shape_idx_to_symbols(int shape_idx, int *shape_sym1,
                                           int *shape_sym2) {
  assert(shape_idx >= 0 && shape_idx < VQ_SHAPES);
  *shape_sym1 = shape_idx / VQ_SHAPE_SYMBOLS_2;
  *shape_sym2 = shape_idx % VQ_SHAPE_SYMBOLS_2;
}

static INLINE int vq_gain_sqrt_quant(uint64_t gain_sq) {
  // This functions returns an index i that specifies the quantized gain,
  // whose value is vq_gain_vals[i].
  // TODO(kslu): use binary search instead of linear search
  for (int i = 0; i < VQ_GAIN_LEVELS - 1; i++)
    if ((int32_t)gain_sq <= vq_gain_sq_thresholds[i + 1]) return i;
  return VQ_GAIN_LEVELS - 1;
}
#endif  // CONFIG_VQ4X4

static INLINE void update_txk_array(TX_TYPE *txk_type, BLOCK_SIZE bsize,
                                    int blk_row, int blk_col, TX_SIZE tx_size,
                                    TX_TYPE tx_type) {
  const int txk_type_idx = av1_get_txk_type_index(bsize, blk_row, blk_col);
  txk_type[txk_type_idx] = tx_type;

  const int txw = tx_size_wide_unit[tx_size];
  const int txh = tx_size_high_unit[tx_size];
  // The 16x16 unit is due to the constraint from tx_64x64 which sets the
  // maximum tx size for chroma as 32x32. Coupled with 4x1 transform block
  // size, the constraint takes effect in 32x16 / 16x32 size too. To solve
  // the intricacy, cover all the 16x16 units inside a 64 level transform.
  if (txw == tx_size_wide_unit[TX_64X64] ||
      txh == tx_size_high_unit[TX_64X64]) {
    const int tx_unit = tx_size_wide_unit[TX_16X16];
    for (int idy = 0; idy < txh; idy += tx_unit) {
      for (int idx = 0; idx < txw; idx += tx_unit) {
        const int this_index =
            av1_get_txk_type_index(bsize, blk_row + idy, blk_col + idx);
        txk_type[this_index] = tx_type;
      }
    }
  }
}

static INLINE TX_TYPE av1_get_tx_type(PLANE_TYPE plane_type,
                                      const MACROBLOCKD *xd, int blk_row,
                                      int blk_col, TX_SIZE tx_size,
                                      int reduced_tx_set) {
  const MB_MODE_INFO *const mbmi = xd->mi[0];
  const struct macroblockd_plane *const pd = &xd->plane[plane_type];
  const TxSetType tx_set_type =
      av1_get_ext_tx_set_type(tx_size, is_inter_block(mbmi), reduced_tx_set);

  TX_TYPE tx_type;
  if (xd->lossless[mbmi->segment_id] || txsize_sqr_up_map[tx_size] > TX_32X32) {
    tx_type = DCT_DCT;
  } else {
    if (plane_type == PLANE_TYPE_Y) {
      const int txk_type_idx =
          av1_get_txk_type_index(mbmi->sb_type, blk_row, blk_col);
      tx_type = mbmi->txk_type[txk_type_idx];
    } else if (is_inter_block(mbmi)) {
      // scale back to y plane's coordinate
      blk_row <<= pd->subsampling_y;
      blk_col <<= pd->subsampling_x;
      const int txk_type_idx =
          av1_get_txk_type_index(mbmi->sb_type, blk_row, blk_col);
      tx_type = mbmi->txk_type[txk_type_idx];
    } else {
      // In intra mode, uv planes don't share the same prediction mode as y
      // plane, so the tx_type should not be shared
      tx_type = intra_mode_to_tx_type(mbmi, PLANE_TYPE_UV);
    }
  }
  assert(tx_type < TX_TYPES);
  if (!av1_ext_tx_used[tx_set_type][tx_type]) return DCT_DCT;
  return tx_type;
}

void av1_setup_block_planes(MACROBLOCKD *xd, int ss_x, int ss_y,
                            const int num_planes);

static INLINE int bsize_to_max_depth(BLOCK_SIZE bsize) {
  TX_SIZE tx_size = max_txsize_rect_lookup[bsize];
  int depth = 0;
  while (depth < MAX_TX_DEPTH && tx_size != TX_4X4) {
    depth++;
    tx_size = sub_tx_size_map[tx_size];
  }
  return depth;
}

static INLINE int bsize_to_tx_size_cat(BLOCK_SIZE bsize) {
  TX_SIZE tx_size = max_txsize_rect_lookup[bsize];
  assert(tx_size != TX_4X4);
  int depth = 0;
  while (tx_size != TX_4X4) {
    depth++;
    tx_size = sub_tx_size_map[tx_size];
    assert(depth < 10);
  }
  assert(depth <= MAX_TX_CATS);
  return depth - 1;
}

static INLINE TX_SIZE depth_to_tx_size(int depth, BLOCK_SIZE bsize) {
  TX_SIZE max_tx_size = max_txsize_rect_lookup[bsize];
  TX_SIZE tx_size = max_tx_size;
  for (int d = 0; d < depth; ++d) tx_size = sub_tx_size_map[tx_size];
  return tx_size;
}

static INLINE TX_SIZE av1_get_adjusted_tx_size(TX_SIZE tx_size) {
  switch (tx_size) {
    case TX_64X64:
    case TX_64X32:
    case TX_32X64: return TX_32X32;
    case TX_64X16: return TX_32X16;
    case TX_16X64: return TX_16X32;
#if CONFIG_FLEX_PARTITION
    case TX_64X8: return TX_32X8;
    case TX_8X64: return TX_8X32;
    case TX_64X4: return TX_32X4;
    case TX_4X64: return TX_4X32;
#endif  // CONFIG_FLEX_PARTITION
    default: return tx_size;
  }
}

static INLINE TX_SIZE av1_get_max_uv_txsize(BLOCK_SIZE bsize, int subsampling_x,
                                            int subsampling_y) {
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsize, subsampling_x, subsampling_y);
  assert(plane_bsize < BLOCK_SIZES_ALL);
  const TX_SIZE uv_tx = max_txsize_rect_lookup[plane_bsize];
  return av1_get_adjusted_tx_size(uv_tx);
}

static INLINE TX_SIZE av1_get_tx_size(int plane, const MACROBLOCKD *xd) {
  const MB_MODE_INFO *mbmi = xd->mi[0];
  if (xd->lossless[mbmi->segment_id]) return TX_4X4;
  if (plane == 0) return mbmi->tx_size;
  const MACROBLOCKD_PLANE *pd = &xd->plane[plane];
  return av1_get_max_uv_txsize(mbmi->sb_type, pd->subsampling_x,
                               pd->subsampling_y);
}

void av1_reset_skip_context(MACROBLOCKD *xd, int mi_row, int mi_col,
                            BLOCK_SIZE bsize, const int num_planes);

void av1_reset_loop_filter_delta(MACROBLOCKD *xd, int num_planes);

void av1_reset_loop_restoration(MACROBLOCKD *xd, const int num_planes);

typedef void (*foreach_transformed_block_visitor)(int plane, int block,
                                                  int blk_row, int blk_col,
                                                  BLOCK_SIZE plane_bsize,
                                                  TX_SIZE tx_size, void *arg);

void av1_set_contexts(const MACROBLOCKD *xd, struct macroblockd_plane *pd,
                      int plane, BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                      int has_eob, int aoff, int loff);

#define MAX_INTERINTRA_SB_SQUARE 32 * 32
static INLINE int is_interintra_mode(const MB_MODE_INFO *mbmi) {
  return (mbmi->ref_frame[0] > INTRA_FRAME &&
          mbmi->ref_frame[1] == INTRA_FRAME);
}

static INLINE int is_interintra_allowed_bsize(const BLOCK_SIZE bsize) {
  return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}

static INLINE int is_interintra_allowed_mode(const PREDICTION_MODE mode) {
  return (mode >= SINGLE_INTER_MODE_START) && (mode < SINGLE_INTER_MODE_END);
}

static INLINE int is_interintra_allowed_ref(const MV_REFERENCE_FRAME rf[2]) {
  return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}

static INLINE int is_interintra_allowed(const MB_MODE_INFO *mbmi) {
  return is_interintra_allowed_bsize(mbmi->sb_type) &&
         is_interintra_allowed_mode(mbmi->mode) &&
         is_interintra_allowed_ref(mbmi->ref_frame);
}

static INLINE int is_interintra_allowed_bsize_group(int group) {
  int i;
  for (i = 0; i < BLOCK_SIZES_ALL; i++) {
    if (size_group_lookup[i] == group &&
        is_interintra_allowed_bsize((BLOCK_SIZE)i)) {
      return 1;
    }
  }
  return 0;
}

static INLINE int is_interintra_pred(const MB_MODE_INFO *mbmi) {
  return mbmi->ref_frame[0] > INTRA_FRAME &&
         mbmi->ref_frame[1] == INTRA_FRAME && is_interintra_allowed(mbmi);
}

static INLINE int get_vartx_max_txsize(const MACROBLOCKD *xd, BLOCK_SIZE bsize,
                                       int plane) {
  if (xd->lossless[xd->mi[0]->segment_id]) return TX_4X4;
  assert(bsize < BLOCK_SIZES_ALL);
  const TX_SIZE max_txsize = max_txsize_rect_lookup[bsize];
  if (plane == 0) return max_txsize;            // luma
  return av1_get_adjusted_tx_size(max_txsize);  // chroma
}

static INLINE int is_motion_variation_allowed_bsize(BLOCK_SIZE bsize) {
  assert(bsize < BLOCK_SIZES_ALL);
  return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}

static INLINE int is_motion_variation_allowed_compound(
    const MB_MODE_INFO *mbmi) {
  if (!has_second_ref(mbmi))
    return 1;
  else
    return 0;
}

// input: log2 of length, 0(4), 1(8), ...
static const int max_neighbor_obmc[6] = { 0, 1, 2, 3, 4, 4 };

static INLINE int check_num_overlappable_neighbors(const MB_MODE_INFO *mbmi) {
  return !(mbmi->overlappable_neighbors[0] == 0 &&
           mbmi->overlappable_neighbors[1] == 0);
}

static INLINE MOTION_MODE
motion_mode_allowed(const WarpedMotionParams *gm_params, const MACROBLOCKD *xd,
                    const MB_MODE_INFO *mbmi, int allow_warped_motion) {
  if (xd->cur_frame_force_integer_mv == 0) {
    const TransformationType gm_type = gm_params[mbmi->ref_frame[0]].wmtype;
    if (is_global_mv_block(mbmi, gm_type)) return SIMPLE_TRANSLATION;
  }
  if (is_motion_variation_allowed_bsize(mbmi->sb_type) &&
      is_inter_mode(mbmi->mode) && mbmi->ref_frame[1] != INTRA_FRAME &&
      is_motion_variation_allowed_compound(mbmi)) {
    if (!check_num_overlappable_neighbors(mbmi)) return SIMPLE_TRANSLATION;
    assert(!has_second_ref(mbmi));
    if (mbmi->num_proj_ref >= 1 &&
        (allow_warped_motion &&
         !av1_is_scaled(xd->block_ref_scale_factors[0]))) {
      if (xd->cur_frame_force_integer_mv) {
        return OBMC_CAUSAL;
      }
      return WARPED_CAUSAL;
    }
    return OBMC_CAUSAL;
  } else {
    return SIMPLE_TRANSLATION;
  }
}

static INLINE void assert_motion_mode_valid(MOTION_MODE mode,
                                            const WarpedMotionParams *gm_params,
                                            const MACROBLOCKD *xd,
                                            const MB_MODE_INFO *mbmi,
                                            int allow_warped_motion) {
  const MOTION_MODE last_motion_mode_allowed =
      motion_mode_allowed(gm_params, xd, mbmi, allow_warped_motion);

  // Check that the input mode is not illegal
  if (last_motion_mode_allowed < mode)
    assert(0 && "Illegal motion mode selected");
}

static INLINE int is_neighbor_overlappable(const MB_MODE_INFO *mbmi) {
  return (is_inter_block(mbmi));
}

static INLINE int av1_allow_palette(int allow_screen_content_tools,
                                    BLOCK_SIZE sb_type) {
  assert(sb_type < BLOCK_SIZES_ALL);
  return allow_screen_content_tools && block_size_wide[sb_type] <= 64 &&
         block_size_high[sb_type] <= 64 && sb_type >= BLOCK_8X8;
}

// Returns sub-sampled dimensions of the given block.
// The output values for 'rows_within_bounds' and 'cols_within_bounds' will
// differ from 'height' and 'width' when part of the block is outside the
// right
// and/or bottom image boundary.
static INLINE void av1_get_block_dimensions(BLOCK_SIZE bsize, int plane,
                                            const MACROBLOCKD *xd, int *width,
                                            int *height,
                                            int *rows_within_bounds,
                                            int *cols_within_bounds) {
  const int block_height = block_size_high[bsize];
  const int block_width = block_size_wide[bsize];
  const int block_rows = (xd->mb_to_bottom_edge >= 0)
                             ? block_height
                             : (xd->mb_to_bottom_edge >> 3) + block_height;
  const int block_cols = (xd->mb_to_right_edge >= 0)
                             ? block_width
                             : (xd->mb_to_right_edge >> 3) + block_width;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  assert(IMPLIES(plane == PLANE_TYPE_Y, pd->subsampling_x == 0));
  assert(IMPLIES(plane == PLANE_TYPE_Y, pd->subsampling_y == 0));
  assert(block_width >= block_cols);
  assert(block_height >= block_rows);
  const int plane_block_width = block_width >> pd->subsampling_x;
  const int plane_block_height = block_height >> pd->subsampling_y;
  // Special handling for chroma sub8x8.
  const int is_chroma_sub8_x = plane > 0 && plane_block_width < 4;
  const int is_chroma_sub8_y = plane > 0 && plane_block_height < 4;
  if (width) *width = plane_block_width + 2 * is_chroma_sub8_x;
  if (height) *height = plane_block_height + 2 * is_chroma_sub8_y;
  if (rows_within_bounds) {
    *rows_within_bounds =
        (block_rows >> pd->subsampling_y) + 2 * is_chroma_sub8_y;
  }
  if (cols_within_bounds) {
    *cols_within_bounds =
        (block_cols >> pd->subsampling_x) + 2 * is_chroma_sub8_x;
  }
}

/* clang-format off */
typedef aom_cdf_prob (*MapCdf)[PALETTE_COLOR_INDEX_CONTEXTS]
                              [CDF_SIZE(PALETTE_COLORS)];
typedef const int (*ColorCost)[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS]
                              [PALETTE_COLORS];
/* clang-format on */

typedef struct {
  int rows;
  int cols;
  int n_colors;
  int plane_width;
  int plane_height;
  uint8_t *color_map;
  MapCdf map_cdf;
  ColorCost color_cost;
} Av1ColorMapParam;

static INLINE int is_nontrans_global_motion(const MACROBLOCKD *xd,
                                            const MB_MODE_INFO *mbmi) {
  int ref;

  // First check if all modes are GLOBALMV
  if (mbmi->mode != GLOBALMV && mbmi->mode != GLOBAL_GLOBALMV) return 0;

  if (AOMMIN(mi_size_wide[mbmi->sb_type], mi_size_high[mbmi->sb_type]) < 2)
    return 0;

  // Now check if all global motion is non translational
  for (ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
    if (xd->global_motion[mbmi->ref_frame[ref]].wmtype == TRANSLATION) return 0;
  }
  return 1;
}

static INLINE PLANE_TYPE get_plane_type(int plane) {
  return (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
}

static INLINE int av1_get_max_eob(TX_SIZE tx_size) {
  if (tx_size == TX_64X64 || tx_size == TX_64X32 || tx_size == TX_32X64) {
    return 1024;
  }
  if (tx_size == TX_16X64 || tx_size == TX_64X16) {
    return 512;
  }
  return tx_size_2d[tx_size];
}

static INLINE int av1_pixels_to_mi(int pixels) {
  return ALIGN_POWER_OF_TWO(pixels, 3) >> MI_SIZE_LOG2;
}

#if CONFIG_INTRA_ENTROPY
// Calculate histogram of gradient orientations of the reconstructed pixel
// values in current coding block.
void av1_get_gradient_hist(const MACROBLOCKD *const xd,
                           MB_MODE_INFO *const mbmi, BLOCK_SIZE bsize);
// Calculate variance of the reconstructed pixel values in current coding block.
void av1_get_recon_var(const MACROBLOCKD *const xd, MB_MODE_INFO *const mbmi,
                       BLOCK_SIZE bsize);
#endif  // CONFIG_INTRA_ENTROPY

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_BLOCKD_H_
