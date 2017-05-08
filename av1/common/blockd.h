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

#ifndef AV1_COMMON_BLOCKD_H_
#define AV1_COMMON_BLOCKD_H_

#include "./aom_config.h"

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
#if CONFIG_PVQ
#include "av1/common/pvq.h"
#include "av1/common/pvq_state.h"
#include "av1/decoder/decint.h"
#endif
#if CONFIG_CFL
#include "av1/common/cfl.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif

#define SUB8X8_COMP_REF (!(CONFIG_CB4X4 && CONFIG_CHROMA_2X2))

#define MAX_MB_PLANE 3

#if CONFIG_EXT_INTER

#if CONFIG_COMPOUND_SEGMENT
// Set COMPOUND_SEGMENT_TYPE to one of the three
// 0: Uniform
// 1: Difference weighted
#define COMPOUND_SEGMENT_TYPE 1

#if COMPOUND_SEGMENT_TYPE == 0
#define MAX_SEG_MASK_BITS 1
// SEG_MASK_TYPES should not surpass 1 << MAX_SEG_MASK_BITS
typedef enum {
  UNIFORM_45 = 0,
  UNIFORM_45_INV,
  SEG_MASK_TYPES,
} SegMaskType;

#elif COMPOUND_SEGMENT_TYPE == 1
#define MAX_SEG_MASK_BITS 1
// SEG_MASK_TYPES should not surpass 1 << MAX_SEG_MASK_BITS
typedef enum {
  DIFFWTD_42 = 0,
  DIFFWTD_42_INV,
  SEG_MASK_TYPES,
} SegMaskType;

#endif  // COMPOUND_SEGMENT_TYPE
#endif  // CONFIG_COMPOUND_SEGMENT
#endif  // CONFIG_EXT_INTER

typedef enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  FRAME_TYPES,
} FrameType;

static INLINE int is_inter_mode(PredictionMode mode) {
#if CONFIG_EXT_INTER
  return mode >= NEARESTMV && mode <= NEW_NEWMV;
#else
  return mode >= NEARESTMV && mode <= NEWMV;
#endif  // CONFIG_EXT_INTER
}

#if CONFIG_PVQ
typedef struct PvqInfo {
  int theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  OdCoeff y[OD_TXSIZE_MAX * OD_TXSIZE_MAX];
  int nb_bands;
  int off[PVQ_MAX_PARTITIONS];
  int size[PVQ_MAX_PARTITIONS];
  int skip_rest;
  int skip_dir;
  int bs;  // log of the block size minus two,
           // i.e. equivalent to aom's TxSize
  // Block skip info, indicating whether DC/AC, is coded.
  PvqSkipType ac_dc_coded;  // bit0: DC coded, bit1 : AC coded (1 means coded)
  TranLowT dq_dc_residue;
} PvqInfo;

typedef struct PvqQueue {
  PvqInfo *buf;  // buffer for pvq info, stored in encoding order
  int curr_pos;  // curr position to write PvqInfo
  int buf_len;   // allocated buffer length
  int last_pos;  // last written position of PvqInfo in a tile
} PvqQueue;
#endif

typedef struct {
  uint8_t *plane[MAX_MB_PLANE];
  int stride[MAX_MB_PLANE];
} BufferSet;

#if CONFIG_EXT_INTER
static INLINE int is_inter_singleref_mode(PredictionMode mode) {
  return mode >= NEARESTMV && mode <= NEWMV;
}
#if CONFIG_COMPOUND_SINGLEREF
static INLINE int is_inter_singleref_comp_mode(PredictionMode mode) {
  return mode >= SR_NEAREST_NEARMV && mode <= SR_NEW_NEWMV;
}
#endif  // CONFIG_COMPOUND_SINGLEREF
static INLINE int is_inter_compound_mode(PredictionMode mode) {
  return mode >= NEAREST_NEARESTMV && mode <= NEW_NEWMV;
}

static INLINE PredictionMode compound_ref0_mode(PredictionMode mode) {
  static PredictionMode lut[MB_MODE_COUNT] = {
    MB_MODE_COUNT,  // DC_PRED
    MB_MODE_COUNT,  // V_PRED
    MB_MODE_COUNT,  // H_PRED
    MB_MODE_COUNT,  // D45_PRED
    MB_MODE_COUNT,  // D135_PRED
    MB_MODE_COUNT,  // D117_PRED
    MB_MODE_COUNT,  // D153_PRED
    MB_MODE_COUNT,  // D207_PRED
    MB_MODE_COUNT,  // D63_PRED
#if CONFIG_ALT_INTRA
    MB_MODE_COUNT,  // SMOOTH_PRED
#endif              // CONFIG_ALT_INTRA
    MB_MODE_COUNT,  // TM_PRED
    MB_MODE_COUNT,  // NEARESTMV
    MB_MODE_COUNT,  // NEARMV
    MB_MODE_COUNT,  // ZEROMV
    MB_MODE_COUNT,  // NEWMV
#if CONFIG_COMPOUND_SINGLEREF
    NEARESTMV,  // SR_NEAREST_NEARMV
    NEARESTMV,  // SR_NEAREST_NEWMV
    NEARMV,     // SR_NEAR_NEWMV
    ZEROMV,     // SR_ZERO_NEWMV
    NEWMV,      // SR_NEW_NEWMV
#endif          // CONFIG_COMPOUND_SINGLEREF
    NEARESTMV,  // NEAREST_NEARESTMV
    NEARESTMV,  // NEAREST_NEARMV
    NEARMV,     // NEAR_NEARESTMV
    NEARMV,     // NEAR_NEARMV
    NEARESTMV,  // NEAREST_NEWMV
    NEWMV,      // NEW_NEARESTMV
    NEARMV,     // NEAR_NEWMV
    NEWMV,      // NEW_NEARMV
    ZEROMV,     // ZERO_ZEROMV
    NEWMV,      // NEW_NEWMV
  };
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE PredictionMode compound_ref1_mode(PredictionMode mode) {
  static PredictionMode lut[MB_MODE_COUNT] = {
    MB_MODE_COUNT,  // DC_PRED
    MB_MODE_COUNT,  // V_PRED
    MB_MODE_COUNT,  // H_PRED
    MB_MODE_COUNT,  // D45_PRED
    MB_MODE_COUNT,  // D135_PRED
    MB_MODE_COUNT,  // D117_PRED
    MB_MODE_COUNT,  // D153_PRED
    MB_MODE_COUNT,  // D207_PRED
    MB_MODE_COUNT,  // D63_PRED
#if CONFIG_ALT_INTRA
    MB_MODE_COUNT,  // SMOOTH_PRED
#endif              // CONFIG_ALT_INTRA
    MB_MODE_COUNT,  // TM_PRED
    MB_MODE_COUNT,  // NEARESTMV
    MB_MODE_COUNT,  // NEARMV
    MB_MODE_COUNT,  // ZEROMV
    MB_MODE_COUNT,  // NEWMV
#if CONFIG_COMPOUND_SINGLEREF
    NEARMV,     // SR_NEAREST_NEARMV
    NEWMV,      // SR_NEAREST_NEWMV
    NEWMV,      // SR_NEAR_NEWMV
    NEWMV,      // SR_ZERO_NEWMV
    NEWMV,      // SR_NEW_NEWMV
#endif          // CONFIG_COMPOUND_SINGLEREF
    NEARESTMV,  // NEAREST_NEARESTMV
    NEARMV,     // NEAREST_NEARMV
    NEARESTMV,  // NEAR_NEARESTMV
    NEARMV,     // NEAR_NEARMV
    NEWMV,      // NEAREST_NEWMV
    NEARESTMV,  // NEW_NEARESTMV
    NEWMV,      // NEAR_NEWMV
    NEARMV,     // NEW_NEARMV
    ZEROMV,     // ZERO_ZEROMV
    NEWMV,      // NEW_NEWMV
  };
  assert(is_inter_compound_mode(mode));
  return lut[mode];
}

static INLINE int have_nearmv_in_inter_mode(PredictionMode mode) {
  return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAREST_NEARMV ||
          mode == NEAR_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE int have_newmv_in_inter_mode(PredictionMode mode) {
  return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
          mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE int use_masked_motion_search(CompoundType type) {
#if CONFIG_WEDGE
  return (type == COMPOUND_WEDGE);
#else
  (void)type;
  return 0;
#endif
}

static INLINE int is_masked_compound_type(CompoundType type) {
#if CONFIG_COMPOUND_SEGMENT && CONFIG_WEDGE
  return (type == COMPOUND_WEDGE || type == COMPOUND_SEG);
#elif !CONFIG_COMPOUND_SEGMENT && CONFIG_WEDGE
  return (type == COMPOUND_WEDGE);
#elif CONFIG_COMPOUND_SEGMENT && !CONFIG_WEDGE
  return (type == COMPOUND_SEG);
#endif  // CONFIG_COMPOUND_SEGMENT
  (void)type;
  return 0;
}
#else

static INLINE int have_nearmv_in_inter_mode(PredictionMode mode) {
  return (mode == NEARMV);
}

static INLINE int have_newmv_in_inter_mode(PredictionMode mode) {
  return (mode == NEWMV);
}
#endif  // CONFIG_EXT_INTER

/* For keyframes, intra block modes are predicted by the (already decoded)
   modes for the Y blocks to the left and above us; for interframes, there
   is a single probability table. */

typedef struct {
  PredictionMode as_mode;
  IntMv as_mv[2];  // first, second inter predictor motion vectors
  IntMv pred_mv[2];
#if CONFIG_EXT_INTER
  IntMv ref_mv[2];
#endif  // CONFIG_EXT_INTER
} BModeInfo;

typedef int8_t MvReferenceFrame;

#if CONFIG_PALETTE
typedef struct {
  // Number of base colors for Y (0) and UV (1)
  uint8_t palette_size[2];
// Value of base colors for Y, U, and V
#if CONFIG_HIGHBITDEPTH
  uint16_t palette_colors[3 * PALETTE_MAX_SIZE];
#else
  uint8_t palette_colors[3 * PALETTE_MAX_SIZE];
#endif  // CONFIG_HIGHBITDEPTH
  // Only used by encoder to store the color index of the top left pixel.
  // TODO(huisu): move this to encoder
  uint8_t palette_first_color_idx[2];
} PaletteModeInfo;
#endif  // CONFIG_PALETTE

#if CONFIG_FILTER_INTRA
#define USE_3TAP_INTRA_FILTER 1  // 0: 4-tap; 1: 3-tap
typedef struct {
  // 1: an ext intra mode is used; 0: otherwise.
  uint8_t use_filter_intra_mode[PLANE_TYPES];
  FilterIntraMode filter_intra_mode[PLANE_TYPES];
} FilterIntraModeInfo;
#endif  // CONFIG_FILTER_INTRA

#if CONFIG_VAR_TX
#if CONFIG_RD_DEBUG
#define TXB_COEFF_COST_MAP_SIZE (2 * MAX_MIB_SIZE)
#endif
#endif

typedef struct RdStats {
  int rate;
  int64_t dist;
  // Please be careful of using rdcost, it's not guaranteed to be set all the
  // time.
  // TODO(angiebird): Create a set of functions to manipulate the RdStats. In
  // these functions, make sure rdcost is always up-to-date according to
  // rate/dist.
  int64_t rdcost;
  int64_t sse;
  int skip;  // sse should equal to dist when skip == 1
#if CONFIG_RD_DEBUG
  int txb_coeff_cost[MAX_MB_PLANE];
#if CONFIG_VAR_TX
  int txb_coeff_cost_map[MAX_MB_PLANE][TXB_COEFF_COST_MAP_SIZE]
                        [TXB_COEFF_COST_MAP_SIZE];
#endif  // CONFIG_VAR_TX
#endif  // CONFIG_RD_DEBUG
} RdStats;

#if CONFIG_EXT_INTER
// This struct is used to group function args that are commonly
// sent together in functions related to interinter compound modes
typedef struct {
#if CONFIG_WEDGE
  int wedge_index;
  int wedge_sign;
#endif  // CONFIG_WEDGE
#if CONFIG_COMPOUND_SEGMENT
  SegMaskType mask_type;
  uint8_t *seg_mask;
#endif  // CONFIG_COMPOUND_SEGMENT
  CompoundType interinter_compound_type;
} InterinterCompoundData;
#endif  // CONFIG_EXT_INTER

// This structure now relates to 8x8 block regions.
typedef struct MbModeInfo {
  // Common for both INTER and INTRA blocks
  BlockSize sb_type;
  PredictionMode mode;
  TxSize tx_size;
#if CONFIG_VAR_TX
  // TODO(jingning): This effectively assigned a separate entry for each
  // 8x8 block. Apparently it takes much more space than needed.
  TxSize inter_tx_size[MAX_MIB_SIZE][MAX_MIB_SIZE];
  TxSize min_tx_size;
#endif
  int8_t skip;
  int8_t segment_id;
#if CONFIG_SUPERTX
  // Minimum of all segment IDs under the current supertx block.
  int8_t segment_id_supertx;
#endif                      // CONFIG_SUPERTX
  int8_t seg_id_predicted;  // valid only when temporal_update is enabled

  // Only for INTRA blocks
  PredictionMode uv_mode;
#if CONFIG_PALETTE
  PaletteModeInfo palette_mode_info;
#endif  // CONFIG_PALETTE
#if CONFIG_INTRABC
  uint8_t use_intrabc;
#endif  // CONFIG_INTRABC

// Only for INTER blocks
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4];
#else
  InterpFilter interp_filter;
#endif
  MvReferenceFrame ref_frame[2];
  TxType tx_type;
#if CONFIG_TXK_SEL
  TxType txk_type[MAX_SB_SQUARE / (TX_SIZE_W_MIN * TX_SIZE_H_MIN)];
#endif

#if CONFIG_FILTER_INTRA
  FilterIntraModeInfo filter_intra_mode_info;
#endif  // CONFIG_FILTER_INTRA
#if CONFIG_EXT_INTRA
  // The actual prediction angle is the base angle + (angle_delta * step).
  int8_t angle_delta[2];
#if CONFIG_INTRA_INTERP
  // To-Do (huisu): this may be replaced by interp_filter
  IntraFilter intra_filter;
#endif  // CONFIG_INTRA_INTERP
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTER
  // interintra members
  InterintraMode interintra_mode;
  // TODO(debargha): Consolidate these flags
  int use_wedge_interintra;
  int interintra_wedge_index;
  int interintra_wedge_sign;
  // interinter members
  CompoundType interinter_compound_type;
#if CONFIG_WEDGE
  int wedge_index;
  int wedge_sign;
#endif  // CONFIG_WEDGE
#if CONFIG_COMPOUND_SEGMENT
  SegMaskType mask_type;
#endif  // CONFIG_COMPOUND_SEGMENT
#endif  // CONFIG_EXT_INTER
  MotionMode motion_mode;
#if CONFIG_MOTION_VAR
  int overlappable_neighbors[2];
#endif  // CONFIG_MOTION_VAR
  IntMv mv[2];
  IntMv pred_mv[2];
  uint8_t ref_mv_idx;
#if CONFIG_EXT_PARTITION_TYPES
  PartitionType partition;
#endif
#if CONFIG_NEW_QUANT
  int dq_off_index;
  int send_dq_bit;
#endif  // CONFIG_NEW_QUANT
  /* deringing gain *per-superblock* */
  int8_t cdef_strength;
#if CONFIG_DELTA_Q
  int current_q_index;
#if CONFIG_EXT_DELTA_Q
  int current_delta_lf_from_base;
#endif
#endif
#if CONFIG_RD_DEBUG
  RdStats rd_stats;
  int mi_row;
  int mi_col;
#endif
#if CONFIG_WARPED_MOTION
  int num_proj_ref[2];
  WarpedMotionParams wm_params[2];
#endif  // CONFIG_WARPED_MOTION

#if CONFIG_CFL
  // Index of the alpha Cb and alpha Cr combination
  int cfl_alpha_ind;
  // Signs of alpha Cb and alpha Cr
  CflSignType cfl_alpha_signs[CFL_PRED_PLANES];
#endif

  BoundaryType boundary_info;
} MbModeInfo;

typedef struct ModeInfo {
  MbModeInfo mbmi;
  BModeInfo bmi[4];
} ModeInfo;

#if CONFIG_INTRABC
static INLINE int is_intrabc_block(const MbModeInfo *mbmi) {
  return mbmi->use_intrabc;
}
#endif

static INLINE PredictionMode get_y_mode(const ModeInfo *mi, int block) {
#if CONFIG_CB4X4
  (void)block;
  return mi->mbmi.mode;
#else
  return mi->mbmi.sb_type < BLOCK_8X8 ? mi->bmi[block].as_mode : mi->mbmi.mode;
#endif
}

static INLINE int is_inter_block(const MbModeInfo *mbmi) {
#if CONFIG_INTRABC
  if (is_intrabc_block(mbmi)) return 1;
#endif
  return mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int has_second_ref(const MbModeInfo *mbmi) {
  return mbmi->ref_frame[1] > INTRA_FRAME;
}

PredictionMode av1_left_block_mode(const ModeInfo *cur_mi,
                                   const ModeInfo *left_mi, int b);

PredictionMode av1_above_block_mode(const ModeInfo *cur_mi,
                                    const ModeInfo *above_mi, int b);

#if CONFIG_GLOBAL_MOTION
static INLINE int is_global_mv_block(const ModeInfo *mi, int block,
                                     TransformationType type) {
  PredictionMode mode = get_y_mode(mi, block);
#if GLOBAL_SUB8X8_USED
  const int block_size_allowed = 1;
#else
  const BlockSize bsize = mi->mbmi.sb_type;
  const int block_size_allowed = (bsize >= BLOCK_8X8);
#endif  // GLOBAL_SUB8X8_USED
#if CONFIG_EXT_INTER
  return (mode == ZEROMV || mode == ZERO_ZEROMV) && type > TRANSLATION &&
         block_size_allowed;
#else
  return mode == ZEROMV && type > TRANSLATION && block_size_allowed;
#endif  // CONFIG_EXT_INTER
}
#endif  // CONFIG_GLOBAL_MOTION

enum mv_precision { MV_PRECISION_Q3, MV_PRECISION_Q4 };

struct Buf2d {
  uint8_t *buf;
  uint8_t *buf0;
  int width;
  int height;
  int stride;
};

typedef struct MacroblockdPlane {
  TranLowT *dqcoeff;
  PlaneType plane_type;
  int subsampling_x;
  int subsampling_y;
  struct Buf2d dst;
  struct Buf2d pre[2];
  EntropyContext *above_context;
  EntropyContext *left_context;
  int16_t seg_dequant[MAX_SEGMENTS][2];
#if CONFIG_NEW_QUANT
  DequantValTypeNuq seg_dequant_nuq[MAX_SEGMENTS][QUANT_PROFILES][COEF_BANDS];
#endif
#if CONFIG_PALETTE
  uint8_t *color_index_map;
#endif  // CONFIG_PALETTE

  // number of 4x4s in current block
  uint16_t n4_w, n4_h;
  // log2 of n4_w, n4_h
  uint8_t n4_wl, n4_hl;
  // block size in pixels
  uint8_t width, height;

#if CONFIG_AOM_QM
  const QmValT *seg_iqmatrix[MAX_SEGMENTS][2][TX_SIZES];
#endif
  // encoder
  const int16_t *dequant;
#if CONFIG_NEW_QUANT
  const DequantValTypeNuq *dequant_val_nuq[QUANT_PROFILES];
#endif  // CONFIG_NEW_QUANT
#if CONFIG_AOM_QM
  const QmValT *seg_qmatrix[MAX_SEGMENTS][2][TX_SIZES];
#endif

#if CONFIG_PVQ || CONFIG_DAALA_DIST
  DECLARE_ALIGNED(16, int16_t, pred[MAX_SB_SQUARE]);
  // PVQ: forward transformed predicted image, a reference for PVQ.
  TranLowT *pvq_ref_coeff;
#endif
} MacroblockdPlane;

#define BLOCK_OFFSET(x, i) \
  ((x) + (i) * (1 << (tx_size_wide_log2[0] + tx_size_high_log2[0])))

typedef struct RefBuffer {
  // TODO(dkovalev): idx is not really required and should be removed, now it
  // is used in av1_onyxd_if.c
  int idx;
  Yv12BufferConfig *buf;
  struct ScaleFactors sf;
} RefBuffer;

typedef int16_t EobThresholdMD[TX_SIZES_ALL][TX_TYPES];

typedef struct Macroblockd {
  struct MacroblockdPlane plane[MAX_MB_PLANE];
  uint8_t bmode_blocks_wl;
  uint8_t bmode_blocks_hl;

  FrameCounts *counts;
  TileInfo tile;

  int mi_stride;

  ModeInfo **mi;
  ModeInfo *left_mi;
  ModeInfo *above_mi;
  MbModeInfo *left_mbmi;
  MbModeInfo *above_mbmi;

  int up_available;
  int left_available;
#if CONFIG_CHROMA_SUB8X8
  int chroma_up_available;
  int chroma_left_available;
#endif

  const AomProb (*partition_probs)[PARTITION_TYPES - 1];

  /* Distance of MB away from frame edges */
  int mb_to_left_edge;
  int mb_to_right_edge;
  int mb_to_top_edge;
  int mb_to_bottom_edge;

  FrameContext *fc;

  /* pointers to reference frames */
  const RefBuffer *block_refs[2];

  /* pointer to current frame */
  const Yv12BufferConfig *cur_buf;

#if CONFIG_INTRABC
  /* Scale of the current frame with respect to itself */
  struct ScaleFactors sf_identity;
#endif

  EntropyContext *above_context[MAX_MB_PLANE];
  EntropyContext left_context[MAX_MB_PLANE][2 * MAX_MIB_SIZE];

  PartitionContext *above_seg_context;
  PartitionContext left_seg_context[MAX_MIB_SIZE];

#if CONFIG_VAR_TX
  TxfmContext *above_txfm_context;
  TxfmContext *left_txfm_context;
  TxfmContext left_txfm_context_buffer[MAX_MIB_SIZE];

  TxSize max_tx_size;
#if CONFIG_SUPERTX
  TxSize supertx_size;
#endif
#endif

  // block dimension in the unit of mode_info.
  uint8_t n8_w, n8_h;

  uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
  CandidateMv ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];
  uint8_t is_sec_rect;

#if CONFIG_PVQ
  DaalaDecCtx daala_dec;
#endif
#if CONFIG_EC_ADAPT
  FrameContext *tile_ctx;
#endif
  /* Bit depth: 8, 10, 12 */
  int bd;

  int qindex[MAX_SEGMENTS];
  int lossless[MAX_SEGMENTS];
  int corrupted;

  struct AomInternalErrorInfo *error_info;
#if CONFIG_GLOBAL_MOTION
  WarpedMotionParams *global_motion;
#endif  // CONFIG_GLOBAL_MOTION
#if CONFIG_DELTA_Q
  int prev_qindex;
  int delta_qindex;
  int current_qindex;
#if CONFIG_EXT_DELTA_Q
  // Since actual frame level loop filtering level value is not available
  // at the beginning of the tile (only available during actual filtering)
  // at encoder side.we record the delta_lf (against the frame level loop
  // filtering level) and code the delta between previous superblock's delta
  // lf and current delta lf. It is equivalent to the delta between previous
  // superblock's actual lf and current lf.
  int prev_delta_lf_from_base;
  int current_delta_lf_from_base;
#endif
#endif
#if CONFIG_ADAPT_SCAN
  const EobThresholdMD *eob_threshold_md;
#endif

#if CONFIG_EXT_INTER && CONFIG_COMPOUND_SEGMENT
  DECLARE_ALIGNED(16, uint8_t, seg_mask[2 * MAX_SB_SQUARE]);
#endif  // CONFIG_EXT_INTER && CONFIG_COMPOUND_SEGMENT

#if CONFIG_CFL
  CflCtx *cfl;
#endif
} Macroblockd;

static INLINE BlockSize get_subsize(BlockSize bsize, PartitionType partition) {
  if (partition == PARTITION_INVALID)
    return BLOCK_INVALID;
  else
    return subsize_lookup[partition][bsize];
}

static const TxType intra_mode_to_tx_type_context[INTRA_MODES] = {
  DCT_DCT,    // DC
  ADST_DCT,   // V
  DCT_ADST,   // H
  DCT_DCT,    // D45
  ADST_ADST,  // D135
  ADST_DCT,   // D117
  DCT_ADST,   // D153
  DCT_ADST,   // D207
  ADST_DCT,   // D63
#if CONFIG_ALT_INTRA
  ADST_ADST,  // SMOOTH
#endif        // CONFIG_ALT_INTRA
  ADST_ADST,  // TM
};

#if CONFIG_SUPERTX
static INLINE int supertx_enabled(const MbModeInfo *mbmi) {
  TxSize max_tx_size = txsize_sqr_map[mbmi->tx_size];
  return tx_size_wide[max_tx_size] >
         AOMMIN(block_size_wide[mbmi->sb_type], block_size_high[mbmi->sb_type]);
}
#endif  // CONFIG_SUPERTX

#define USE_TXTYPE_SEARCH_FOR_SUB8X8_IN_CB4X4 1

#if CONFIG_RECT_TX
static INLINE int is_rect_tx(TxSize tx_size) { return tx_size >= TX_SIZES; }
#endif  // CONFIG_RECT_TX

#if CONFIG_EXT_TX
#define ALLOW_INTRA_EXT_TX 1

typedef enum {
  // DCT only
  EXT_TX_SET_DCTONLY = 0,
  // DCT + Identity only
  EXT_TX_SET_DCT_IDTX = 1,
  // Discrete Trig transforms w/o flip (4) + Identity (1)
  EXT_TX_SET_DTT4_IDTX = 2,
  // Discrete Trig transforms w/o flip (4) + Identity (1) + 1D Hor/vert DCT (2)
  EXT_TX_SET_DTT4_IDTX_1DDCT = 3,
  // Discrete Trig transforms w/ flip (9) + Identity (1) + 1D Hor/Ver DCT (2)
  EXT_TX_SET_DTT9_IDTX_1DDCT = 4,
  // Discrete Trig transforms w/ flip (9) + Identity (1) + 1D Hor/Ver (6)
  EXT_TX_SET_ALL16 = 5,
  EXT_TX_SET_TYPES
} TxSetType;

// Number of transform types in each set type
static const int num_ext_tx_set[EXT_TX_SET_TYPES] = { 1, 2, 5, 7, 12, 16 };

// Maps intra set index to the set type
static const int ext_tx_set_type_intra[EXT_TX_SETS_INTRA] = {
  EXT_TX_SET_DCTONLY, EXT_TX_SET_DTT4_IDTX_1DDCT, EXT_TX_SET_DTT4_IDTX
};

// Maps inter set index to the set type
static const int ext_tx_set_type_inter[EXT_TX_SETS_INTER] = {
  EXT_TX_SET_DCTONLY, EXT_TX_SET_ALL16, EXT_TX_SET_DTT9_IDTX_1DDCT,
  EXT_TX_SET_DCT_IDTX
};

// Maps set types above to the indices used for intra
static const int ext_tx_set_index_intra[EXT_TX_SET_TYPES] = { 0, -1, 2,
                                                              1, -1, -1 };

// Maps set types above to the indices used for inter
static const int ext_tx_set_index_inter[EXT_TX_SET_TYPES] = {
  0, 3, -1, -1, 2, 1
};

static INLINE TxSetType get_ext_tx_set_type(TxSize tx_size, BlockSize bs,
                                            int is_inter, int use_reduced_set) {
  const TxSize tx_size2 = txsize_sqr_up_map[tx_size];
  tx_size = txsize_sqr_map[tx_size];
#if CONFIG_CB4X4 && USE_TXTYPE_SEARCH_FOR_SUB8X8_IN_CB4X4
  (void)bs;
  if (tx_size > TX_32X32) return EXT_TX_SET_DCTONLY;
#else
  if (tx_size > TX_32X32 || bs < BLOCK_8X8) return EXT_TX_SET_DCTONLY;
#endif
  if (use_reduced_set)
    return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DTT4_IDTX;
  if (tx_size2 == TX_32X32)
    return is_inter ? EXT_TX_SET_DCT_IDTX : EXT_TX_SET_DCTONLY;
  if (is_inter)
    return (tx_size == TX_16X16 ? EXT_TX_SET_DTT9_IDTX_1DDCT
                                : EXT_TX_SET_ALL16);
  else
    return (tx_size == TX_16X16 ? EXT_TX_SET_DTT4_IDTX
                                : EXT_TX_SET_DTT4_IDTX_1DDCT);
}

static INLINE int get_ext_tx_set(TxSize tx_size, BlockSize bs, int is_inter,
                                 int use_reduced_set) {
  const TxSetType set_type =
      get_ext_tx_set_type(tx_size, bs, is_inter, use_reduced_set);
  return is_inter ? ext_tx_set_index_inter[set_type]
                  : ext_tx_set_index_intra[set_type];
}

static const int use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][EXT_TX_SIZES] =
    {
#if CONFIG_CB4X4
      { 1, 1, 1, 1, 1 },  // unused
      { 0, 1, 1, 0, 0 },
      { 0, 0, 0, 1, 0 },
#else
      { 1, 1, 1, 1 },  // unused
      { 1, 1, 0, 0 },
      { 0, 0, 1, 0 },
#endif  // CONFIG_CB4X4
    };

static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][EXT_TX_SIZES] =
    {
#if CONFIG_CB4X4
      { 1, 1, 1, 1, 1 },  // unused
      { 0, 1, 1, 0, 0 },
      { 0, 0, 0, 1, 0 },
      { 0, 0, 0, 0, 1 },
#else
      { 1, 1, 1, 1 },  // unused
      { 1, 1, 0, 0 },
      { 0, 0, 1, 0 },
      { 0, 0, 0, 1 },
#endif  // CONFIG_CB4X4
    };

// Transform types used in each intra set
static const int ext_tx_used_intra[EXT_TX_SETS_INTRA][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
};

// Numbers of transform types used in each intra set
static const int ext_tx_cnt_intra[EXT_TX_SETS_INTRA] = { 1, 7, 5 };

// Transform types used in each inter set
static const int ext_tx_used_inter[EXT_TX_SETS_INTER][TX_TYPES] = {
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
};

// Numbers of transform types used in each inter set
static const int ext_tx_cnt_inter[EXT_TX_SETS_INTER] = { 1, 16, 12, 2 };

// 1D Transforms used in inter set, this needs to be changed if
// ext_tx_used_inter is changed
static const int ext_tx_used_inter_1D[EXT_TX_SETS_INTER][TX_TYPES_1D] = {
  { 1, 0, 0, 0 }, { 1, 1, 1, 1 }, { 1, 1, 1, 1 }, { 1, 0, 0, 1 },
};

static INLINE int get_ext_tx_types(TxSize tx_size, BlockSize bs, int is_inter,
                                   int use_reduced_set) {
  const int set_type =
      get_ext_tx_set_type(tx_size, bs, is_inter, use_reduced_set);
  return num_ext_tx_set[set_type];
}

#if CONFIG_RECT_TX
static INLINE int is_rect_tx_allowed_bsize(BlockSize bsize) {
  static const char LUT[BLOCK_SIZES] = {
#if CONFIG_CB4X4
    0,  // BLOCK_2X2
    0,  // BLOCK_2X4
    0,  // BLOCK_4X2
#endif
    0,  // BLOCK_4X4
    1,  // BLOCK_4X8
    1,  // BLOCK_8X4
    0,  // BLOCK_8X8
    1,  // BLOCK_8X16
    1,  // BLOCK_16X8
    0,  // BLOCK_16X16
    1,  // BLOCK_16X32
    1,  // BLOCK_32X16
    0,  // BLOCK_32X32
    0,  // BLOCK_32X64
    0,  // BLOCK_64X32
    0,  // BLOCK_64X64
#if CONFIG_EXT_PARTITION
    0,  // BLOCK_64X128
    0,  // BLOCK_128X64
    0,  // BLOCK_128X128
#endif  // CONFIG_EXT_PARTITION
  };

  return LUT[bsize];
}

static INLINE int is_rect_tx_allowed(const Macroblockd *xd,
                                     const MbModeInfo *mbmi) {
  return is_rect_tx_allowed_bsize(mbmi->sb_type) &&
         !xd->lossless[mbmi->segment_id];
}
#endif  // CONFIG_RECT_TX
#endif  // CONFIG_EXT_TX

static INLINE TxSize tx_size_from_tx_mode(BlockSize bsize, TxMode tx_mode,
                                          int is_inter) {
  const TxSize largest_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
#if (CONFIG_VAR_TX || CONFIG_EXT_TX) && CONFIG_RECT_TX
  const TxSize max_rect_tx_size = max_txsize_rect_lookup[bsize];
#else
  const TxSize max_tx_size = max_txsize_lookup[bsize];
#endif  // (CONFIG_VAR_TX || CONFIG_EXT_TX) && CONFIG_RECT_TX
  (void)is_inter;
#if CONFIG_VAR_TX && CONFIG_RECT_TX
#if CONFIG_CB4X4
  if (bsize == BLOCK_4X4)
    return AOMMIN(max_txsize_lookup[bsize], largest_tx_size);
#else
  if (bsize < BLOCK_8X8)
    return AOMMIN(max_txsize_lookup[bsize], largest_tx_size);
#endif
  if (txsize_sqr_map[max_rect_tx_size] <= largest_tx_size)
    return max_rect_tx_size;
  else
    return largest_tx_size;
#elif CONFIG_EXT_TX && CONFIG_RECT_TX
  if (txsize_sqr_up_map[max_rect_tx_size] <= largest_tx_size) {
    return max_rect_tx_size;
  } else {
    return largest_tx_size;
  }
#else
  return AOMMIN(max_tx_size, largest_tx_size);
#endif  // CONFIG_VAR_TX && CONFIG_RECT_TX
}

#if CONFIG_EXT_INTRA
#define MAX_ANGLE_DELTA 3
#define ANGLE_STEP 3
extern const int16_t dr_intra_derivative[90];
static const uint8_t mode_to_angle_map[INTRA_MODES] = {
  0, 90, 180, 45, 135, 111, 157, 203, 67, 0,
};
#if CONFIG_INTRA_INTERP
// Returns whether filter selection is needed for a given
// intra prediction angle.
int av1_is_intra_filter_switchable(int angle);
#endif  // CONFIG_INTRA_INTERP
#endif  // CONFIG_EXT_INTRA

#define FIXED_TX_TYPE 0

// Converts block_index for given transform size to index of the block in raster
// order.
static INLINE int av1_block_index_to_raster_order(TxSize tx_size,
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
static INLINE int av1_raster_order_to_block_index(TxSize tx_size,
                                                  int raster_order) {
  assert(tx_size == TX_4X4 || tx_size == TX_4X8 || tx_size == TX_8X4);
  // We ensure that block indices are 0 & 2 if tx size is 4x8 or 8x4.
  return (tx_size == TX_4X4) ? raster_order : (raster_order > 0) ? 2 : 0;
}

static INLINE TxType get_default_tx_type(PlaneType plane_type,
                                         const Macroblockd *xd, int block_idx,
                                         TxSize tx_size) {
  const MbModeInfo *const mbmi = &xd->mi[0]->mbmi;

  if (is_inter_block(mbmi) || plane_type != PLANE_TYPE_Y ||
      xd->lossless[mbmi->segment_id] || tx_size >= TX_32X32)
    return DCT_DCT;

  return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y
                                           ? get_y_mode(xd->mi[0], block_idx)
                                           : mbmi->uv_mode];
}

static INLINE TxType get_tx_type(PlaneType plane_type, const Macroblockd *xd,
                                 int block, TxSize tx_size) {
  const ModeInfo *const mi = xd->mi[0];
  const MbModeInfo *const mbmi = &mi->mbmi;
#if CONFIG_INTRABC
  // TODO(aconverse@google.com): Revisit this decision
  if (is_intrabc_block(mbmi)) return DCT_DCT;
#endif  // CONFIG_INTRABC
#if !CONFIG_TXK_SEL
#if FIXED_TX_TYPE
  const int block_raster_idx = av1_block_index_to_raster_order(tx_size, block);
  return get_default_tx_type(plane_type, xd, block_raster_idx, tx_size);
#elif CONFIG_EXT_TX
#if !CONFIG_CB4X4
  const int block_raster_idx = av1_block_index_to_raster_order(tx_size, block);
#endif  // !CONFIG_CB4X4
  if (xd->lossless[mbmi->segment_id] || txsize_sqr_map[tx_size] > TX_32X32 ||
      (txsize_sqr_map[tx_size] >= TX_32X32 && !is_inter_block(mbmi)))
    return DCT_DCT;
  if (mbmi->sb_type >= BLOCK_8X8 || CONFIG_CB4X4) {
    if (plane_type == PLANE_TYPE_Y) {
#if !ALLOW_INTRA_EXT_TX
      if (is_inter_block(mbmi))
#endif  // ALLOW_INTRA_EXT_TX
        return mbmi->tx_type;
    }

    if (is_inter_block(mbmi)) {
// UV Inter only
#if CONFIG_CB4X4
      if (tx_size < TX_4X4) return DCT_DCT;
#endif
      return (mbmi->tx_type == IDTX && txsize_sqr_map[tx_size] >= TX_32X32)
                 ? DCT_DCT
                 : mbmi->tx_type;
    }
  }

#if CONFIG_CB4X4
  (void)block;
  if (tx_size < TX_4X4)
    return DCT_DCT;
  else
    return intra_mode_to_tx_type_context[mbmi->uv_mode];
#else

  // Sub8x8-Inter/Intra OR UV-Intra
  if (is_inter_block(mbmi))  // Sub8x8-Inter
    return DCT_DCT;
  else  // Sub8x8 Intra OR UV-Intra
    return intra_mode_to_tx_type_context[plane_type == PLANE_TYPE_Y
                                             ? get_y_mode(mi, block_raster_idx)
                                             : mbmi->uv_mode];
#endif  // CONFIG_CB4X4
#else   // CONFIG_EXT_TX
  (void)block;
  if (plane_type != PLANE_TYPE_Y || xd->lossless[mbmi->segment_id] ||
      txsize_sqr_map[tx_size] >= TX_32X32)
    return DCT_DCT;
  return mbmi->tx_type;
#endif  // CONFIG_EXT_TX
#else   // !CONFIG_TXK_SEL
  (void)tx_size;
  TxType tx_type;
  if (plane_type != PLANE_TYPE_Y || xd->lossless[mbmi->segment_id] ||
      mbmi->tx_size >= TX_32X32) {
    tx_type = DCT_DCT;
  } else {
    tx_type = mbmi->txk_type[block];
  }
  assert(tx_type >= DCT_DCT && tx_type < TX_TYPES);
  return tx_type;
#endif  // !CONFIG_TXK_SEL
}

void av1_setup_block_planes(Macroblockd *xd, int ss_x, int ss_y);

static INLINE int tx_size_to_depth(TxSize tx_size) {
  return (int)(tx_size - TX_4X4);
}

static INLINE TxSize depth_to_tx_size(int depth) {
  return (TxSize)(depth + TX_4X4);
}

static INLINE TxSize get_uv_tx_size(const MbModeInfo *mbmi,
                                    const struct MacroblockdPlane *pd) {
  TxSize uv_txsize;
#if CONFIG_CB4X4
  assert(mbmi->tx_size > TX_2X2);
#endif

#if CONFIG_SUPERTX
  if (supertx_enabled(mbmi))
    return uvsupertx_size_lookup[txsize_sqr_map[mbmi->tx_size]]
                                [pd->subsampling_x][pd->subsampling_y];
#endif  // CONFIG_SUPERTX

  uv_txsize = uv_txsize_lookup[mbmi->sb_type][mbmi->tx_size][pd->subsampling_x]
                              [pd->subsampling_y];
#if CONFIG_CB4X4 && !CONFIG_CHROMA_2X2
  uv_txsize = AOMMAX(uv_txsize, TX_4X4);
#endif
  assert(uv_txsize != TX_INVALID);
  return uv_txsize;
}

static INLINE TxSize get_tx_size(int plane, const Macroblockd *xd) {
  const MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  const MacroblockdPlane *pd = &xd->plane[plane];
  const TxSize tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
  return tx_size;
}

static INLINE BlockSize
get_plane_block_size(BlockSize bsize, const struct MacroblockdPlane *pd) {
  return ss_size_lookup[bsize][pd->subsampling_x][pd->subsampling_y];
}

static INLINE void reset_skip_context(Macroblockd *xd, BlockSize bsize) {
  int i;
  for (i = 0; i < MAX_MB_PLANE; i++) {
    struct MacroblockdPlane *const pd = &xd->plane[i];
    const BlockSize plane_bsize = get_plane_block_size(bsize, pd);
    const int txs_wide = block_size_wide[plane_bsize] >> tx_size_wide_log2[0];
    const int txs_high = block_size_high[plane_bsize] >> tx_size_high_log2[0];
    memset(pd->above_context, 0, sizeof(EntropyContext) * txs_wide);
    memset(pd->left_context, 0, sizeof(EntropyContext) * txs_high);
  }
}

typedef void (*ForeachTransformedBlockVisitor)(int plane, int block,
                                               int blk_row, int blk_col,
                                               BlockSize plane_bsize,
                                               TxSize tx_size, void *arg);

void av1_foreach_transformed_block_in_plane(
    const Macroblockd *const xd, BlockSize bsize, int plane,
    ForeachTransformedBlockVisitor visit, void *arg);

#if CONFIG_LV_MAP
void av1_foreach_transformed_block(const Macroblockd *const xd, BlockSize bsize,
                                   int mi_row, int mi_col,
                                   ForeachTransformedBlockVisitor visit,
                                   void *arg);
#endif

#if CONFIG_DAALA_DIST
void av1_foreach_8x8_transformed_block_in_plane(
    const Macroblockd *const xd, BlockSize bsize, int plane,
    ForeachTransformedBlockVisitor visit,
    ForeachTransformedBlockVisitor mi_visit, void *arg);
#endif

#if CONFIG_COEF_INTERLEAVE
static INLINE int get_max_4x4_size(int num_4x4, int mb_to_edge,
                                   int subsampling) {
  return num_4x4 + (mb_to_edge >= 0 ? 0 : mb_to_edge >> (5 + subsampling));
}

void av1_foreach_transformed_block_interleave(
    const Macroblockd *const xd, BlockSize bsize,
    ForeachTransformedBlockVisitor visit, void *arg);
#endif

void av1_set_contexts(const Macroblockd *xd, struct MacroblockdPlane *pd,
                      int plane, TxSize tx_size, int has_eob, int aoff,
                      int loff);

#if CONFIG_EXT_INTER
static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
#if CONFIG_INTERINTRA
  // TODO(debargha): Should this be bsize < BLOCK_LARGEST?
  return (bsize >= BLOCK_8X8) && (bsize < BLOCK_64X64);
#else
  (void)bsize;
  return 0;
#endif  // CONFIG_INTERINTRA
}

static INLINE int is_interintra_allowed_mode(const PredictionMode mode) {
#if CONFIG_INTERINTRA
  return (mode >= NEARESTMV) && (mode <= NEWMV);
#else
  (void)mode;
  return 0;
#endif  // CONFIG_INTERINTRA
}

static INLINE int is_interintra_allowed_ref(const MvReferenceFrame rf[2]) {
#if CONFIG_INTERINTRA
  return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
#else
  (void)rf;
  return 0;
#endif  // CONFIG_INTERINTRA
}

static INLINE int is_interintra_allowed(const MbModeInfo *mbmi) {
  return is_interintra_allowed_bsize(mbmi->sb_type) &&
         is_interintra_allowed_mode(mbmi->mode) &&
         is_interintra_allowed_ref(mbmi->ref_frame);
}

static INLINE int is_interintra_allowed_bsize_group(int group) {
  int i;
  for (i = 0; i < BLOCK_SIZES; i++) {
    if (size_group_lookup[i] == group &&
        is_interintra_allowed_bsize((BlockSize)i)) {
      return 1;
    }
  }
  return 0;
}

static INLINE int is_interintra_pred(const MbModeInfo *mbmi) {
  return (mbmi->ref_frame[1] == INTRA_FRAME) && is_interintra_allowed(mbmi);
}
#endif  // CONFIG_EXT_INTER

#if CONFIG_VAR_TX
static INLINE int get_vartx_max_txsize(const MbModeInfo *const mbmi,
                                       BlockSize bsize) {
#if CONFIG_CB4X4
  (void)mbmi;
  return max_txsize_rect_lookup[bsize];
#endif  // CONFIG_C4X4
  return mbmi->sb_type < BLOCK_8X8 ? max_txsize_rect_lookup[mbmi->sb_type]
                                   : max_txsize_rect_lookup[bsize];
}
#endif  // CONFIG_VAR_TX

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
static INLINE int is_motion_variation_allowed_bsize(BlockSize bsize) {
  return (bsize >= BLOCK_8X8);
}

static INLINE int is_motion_variation_allowed_compound(const MbModeInfo *mbmi) {
  if (!has_second_ref(mbmi))
    return 1;
  else
    return 0;
}

#if CONFIG_MOTION_VAR
// input: log2 of length, 0(4), 1(8), ...
static const int max_neighbor_obmc[6] = { 0, 1, 2, 3, 4, 4 };

static INLINE int check_num_overlappable_neighbors(const MbModeInfo *mbmi) {
  return !(mbmi->overlappable_neighbors[0] == 0 &&
           mbmi->overlappable_neighbors[1] == 0);
}
#endif

static INLINE MotionMode motion_mode_allowed(
#if CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
    int block, const WarpedMotionParams *gm_params,
#endif  // CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
    const ModeInfo *mi) {
  const MbModeInfo *mbmi = &mi->mbmi;
#if CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
  const TransformationType gm_type = gm_params[mbmi->ref_frame[0]].wmtype;
  if (is_global_mv_block(mi, block, gm_type)) return SIMPLE_TRANSLATION;
#endif  // CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
#if CONFIG_EXT_INTER
  if (is_motion_variation_allowed_bsize(mbmi->sb_type) &&
      is_inter_mode(mbmi->mode) && mbmi->ref_frame[1] != INTRA_FRAME &&
      is_motion_variation_allowed_compound(mbmi)) {
#else
  if (is_motion_variation_allowed_bsize(mbmi->sb_type) &&
      is_inter_mode(mbmi->mode) && is_motion_variation_allowed_compound(mbmi)) {
#endif  // CONFIG_EXT_INTER
#if CONFIG_MOTION_VAR
    if (!check_num_overlappable_neighbors(mbmi)) return SIMPLE_TRANSLATION;
#endif
#if CONFIG_WARPED_MOTION
    if (!has_second_ref(mbmi) && mbmi->num_proj_ref[0] >= 1)
      return WARPED_CAUSAL;
    else
#endif  // CONFIG_WARPED_MOTION
#if CONFIG_MOTION_VAR
      return OBMC_CAUSAL;
#else
    return SIMPLE_TRANSLATION;
#endif  // CONFIG_MOTION_VAR
  } else {
    return SIMPLE_TRANSLATION;
  }
}

static INLINE void assert_motion_mode_valid(MotionMode mode,
#if CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
                                            int block,
                                            const WarpedMotionParams *gm_params,
#endif  // CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
                                            const ModeInfo *mi) {
  const MotionMode last_motion_mode_allowed = motion_mode_allowed(
#if CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
      block, gm_params,
#endif  // CONFIG_GLOBAL_MOTION && SEPARATE_GLOBAL_MOTION
      mi);
  // Check that the input mode is not illegal
  if (last_motion_mode_allowed < mode)
    assert(0 && "Illegal motion mode selected");
}

#if CONFIG_MOTION_VAR
static INLINE int is_neighbor_overlappable(const MbModeInfo *mbmi) {
  return (is_inter_block(mbmi));
}
#endif  // CONFIG_MOTION_VAR
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

// Returns sub-sampled dimensions of the given block.
// The output values for 'rows_within_bounds' and 'cols_within_bounds' will
// differ from 'height' and 'width' when part of the block is outside the right
// and/or bottom image boundary.
static INLINE void av1_get_block_dimensions(BlockSize bsize, int plane,
                                            const Macroblockd *xd, int *width,
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
  const struct MacroblockdPlane *const pd = &xd->plane[plane];
  assert(IMPLIES(plane == PLANE_TYPE_Y, pd->subsampling_x == 0));
  assert(IMPLIES(plane == PLANE_TYPE_Y, pd->subsampling_y == 0));
  assert(block_width >= block_cols);
  assert(block_height >= block_rows);
  if (width) *width = block_width >> pd->subsampling_x;
  if (height) *height = block_height >> pd->subsampling_y;
  if (rows_within_bounds) *rows_within_bounds = block_rows >> pd->subsampling_y;
  if (cols_within_bounds) *cols_within_bounds = block_cols >> pd->subsampling_x;
}

#if CONFIG_GLOBAL_MOTION
static INLINE int is_nontrans_global_motion(const Macroblockd *xd) {
  const ModeInfo *mi = xd->mi[0];
  const MbModeInfo *const mbmi = &mi->mbmi;
  int ref;
#if CONFIG_CB4X4
  const int unify_bsize = 1;
#else
  const int unify_bsize = 0;
#endif

  // First check if all modes are ZEROMV
  if (mbmi->sb_type >= BLOCK_8X8 || unify_bsize) {
#if CONFIG_EXT_INTER
    if (mbmi->mode != ZEROMV && mbmi->mode != ZERO_ZEROMV) return 0;
#else
    if (mbmi->mode != ZEROMV) return 0;
#endif  // CONFIG_EXT_INTER
  } else {
#if CONFIG_EXT_INTER
    if (mi->bmi[0].as_mode != ZEROMV || mi->bmi[1].as_mode != ZEROMV ||
        mi->bmi[2].as_mode != ZEROMV || mi->bmi[3].as_mode != ZEROMV ||
        mi->bmi[0].as_mode != ZERO_ZEROMV ||
        mi->bmi[1].as_mode != ZERO_ZEROMV ||
        mi->bmi[2].as_mode != ZERO_ZEROMV || mi->bmi[3].as_mode != ZERO_ZEROMV)
      return 0;
#else
    if (mi->bmi[0].as_mode != ZEROMV || mi->bmi[1].as_mode != ZEROMV ||
        mi->bmi[2].as_mode != ZEROMV || mi->bmi[3].as_mode != ZEROMV)
      return 0;
#endif  // CONFIG_EXT_INTER
  }

#if !GLOBAL_SUB8X8_USED
  if (mbmi->sb_type < BLOCK_8X8) return 0;
#endif

  // Now check if all global motion is non translational
  for (ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
    if (xd->global_motion[mbmi->ref_frame[ref]].wmtype <= TRANSLATION) return 0;
  }
  return 1;
}
#endif  // CONFIG_GLOBAL_MOTION

static INLINE PlaneType get_plane_type(int plane) {
  return (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_BLOCKD_H_
