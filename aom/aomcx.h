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
#ifndef AOM_AOM_AOMCX_H_
#define AOM_AOM_AOMCX_H_

/*!\defgroup aom_encoder AOMedia AOM/AV1 Encoder
 * \ingroup aom
 *
 * @{
 */
#include "aom/aom.h"
#include "aom/aom_encoder.h"
#include "aom/aom_external_partition.h"
#include "aom/aomcx_defines.h"

/*!\file
 * \brief Provides definitions for using AOM or AV1 encoder algorithm within the
 *        aom Codec Interface.
 *
 * Several interfaces are excluded with CONFIG_REALTIME_ONLY build:
 * Global motion
 * Warped motion
 * OBMC
 * TPL model
 * Loop restoration
 *
 * The following features are also disabled with CONFIG_REALTIME_ONLY:
 * CNN
 * 4X rectangular blocks
 * 4X rectangular transform in intra prediction
 */

#ifdef __cplusplus
extern "C" {
#endif

/*!\name Algorithm interface for AV1
 *
 * This interface provides the capability to encode raw AV1 streams.
 *@{
 */

/*!\brief A single instance of the AV1 encoder.
 *\deprecated This access mechanism is provided for backwards compatibility;
 * prefer aom_codec_av1_cx().
 */
extern aom_codec_iface_t aom_codec_av1_cx_algo;

/*!\brief The interface to the AV1 encoder.
 */
extern aom_codec_iface_t *aom_codec_av1_cx(void);
/*!@} - end algorithm interface member group */

/*!\brief  aom region of interest map
 *
 * These defines the data structures for the region of interest map
 *
 * TODO(yaowu): create a unit test for ROI map related APIs
 *
 */
typedef struct aom_roi_map {
  /*! An id between 0 and 7 for each 8x8 region within a frame. */
  unsigned char *roi_map;
  unsigned int rows;              /**< Number of rows. */
  unsigned int cols;              /**< Number of columns. */
  int delta_q[AOM_MAX_SEGMENTS];  /**< Quantizer deltas. */
  int delta_lf[AOM_MAX_SEGMENTS]; /**< Loop filter deltas. */
  /*! Static breakout threshold for each segment. */
  unsigned int static_threshold[AOM_MAX_SEGMENTS];
} aom_roi_map_t;

/*!\brief  aom active region map
 *
 * These defines the data structures for active region map
 *
 */

typedef struct aom_active_map {
  /*!\brief specify an on (1) or off (0) each 16x16 region within a frame */
  unsigned char *active_map;
  unsigned int rows; /**< number of rows */
  unsigned int cols; /**< number of cols */
} aom_active_map_t;

/*!\brief  aom image scaling mode
 *
 * This defines the data structure for image scaling mode
 *
 */
typedef struct aom_scaling_mode {
  AOM_SCALING_MODE h_scaling_mode; /**< horizontal scaling mode */
  AOM_SCALING_MODE v_scaling_mode; /**< vertical scaling mode   */
} aom_scaling_mode_t;

/*!brief Struct for spatial and temporal layer ID */
typedef struct aom_svc_layer_id {
  int spatial_layer_id;  /**< Spatial layer ID */
  int temporal_layer_id; /**< Temporal layer ID */
} aom_svc_layer_id_t;

/*!brief Parameter type for SVC */
typedef struct aom_svc_params {
  int number_spatial_layers;                 /**< Number of spatial layers */
  int number_temporal_layers;                /**< Number of temporal layers */
  int max_quantizers[AOM_MAX_LAYERS];        /**< Max Q for each layer */
  int min_quantizers[AOM_MAX_LAYERS];        /**< Min Q for each layer */
  int scaling_factor_num[AOM_MAX_SS_LAYERS]; /**< Scaling factor-numerator */
  int scaling_factor_den[AOM_MAX_SS_LAYERS]; /**< Scaling factor-denominator */
  /*! Target bitrate for each layer */
  int layer_target_bitrate[AOM_MAX_LAYERS];
  /*! Frame rate factor for each temporal layer */
  int framerate_factor[AOM_MAX_TS_LAYERS];
} aom_svc_params_t;

/*!brief Parameters for setting ref frame config */
typedef struct aom_svc_ref_frame_config {
  // 7 references: LAST_FRAME (0), LAST2_FRAME(1), LAST3_FRAME(2),
  // GOLDEN_FRAME(3), BWDREF_FRAME(4), ALTREF2_FRAME(5), ALTREF_FRAME(6).
  int reference[7]; /**< Reference flag for each of the 7 references. */
  /*! Buffer slot index for each of 7 references. */
  int ref_idx[7];
  int refresh[8]; /**< Refresh flag for each of the 8 slots. */
} aom_svc_ref_frame_config_t;

/*!brief Parameters for setting ref frame compound prediction */
typedef struct aom_svc_ref_frame_comp_pred {
  // Use compound prediction for the ref_frame pairs GOLDEN_LAST (0),
  // LAST2_LAST (1), and ALTREF_LAST (2).
  int use_comp_pred[3]; /**<Compound reference flag. */
} aom_svc_ref_frame_comp_pred_t;

/*!\cond */
/*!\brief Encoder control function parameter type
 *
 * Defines the data types that AOME/AV1E control functions take.
 *
 * \note Additional common controls are defined in aom.h.
 *
 * \note For each control ID "X", a macro-define of
 * AOM_CTRL_X is provided. It is used at compile time to determine
 * if the control ID is supported by the libaom library available,
 * when the libaom version cannot be controlled.
 */
AOM_CTRL_USE_TYPE(AOME_USE_REFERENCE, int)
#define AOM_CTRL_AOME_USE_REFERENCE

AOM_CTRL_USE_TYPE(AOME_SET_ROI_MAP, aom_roi_map_t *)
#define AOM_CTRL_AOME_SET_ROI_MAP

AOM_CTRL_USE_TYPE(AOME_SET_ACTIVEMAP, aom_active_map_t *)
#define AOM_CTRL_AOME_SET_ACTIVEMAP

AOM_CTRL_USE_TYPE(AOME_SET_SCALEMODE, aom_scaling_mode_t *)
#define AOM_CTRL_AOME_SET_SCALEMODE

AOM_CTRL_USE_TYPE(AOME_SET_SPATIAL_LAYER_ID, unsigned int)
#define AOM_CTRL_AOME_SET_SPATIAL_LAYER_ID

AOM_CTRL_USE_TYPE(AOME_SET_CPUUSED, int)
#define AOM_CTRL_AOME_SET_CPUUSED

AOM_CTRL_USE_TYPE(AOME_SET_ENABLEAUTOALTREF, unsigned int)
#define AOM_CTRL_AOME_SET_ENABLEAUTOALTREF

AOM_CTRL_USE_TYPE(AOME_SET_ENABLEAUTOBWDREF, unsigned int)
#define AOM_CTRL_AOME_SET_ENABLEAUTOBWDREF

AOM_CTRL_USE_TYPE(AOME_SET_SHARPNESS, unsigned int)
#define AOM_CTRL_AOME_SET_SHARPNESS

AOM_CTRL_USE_TYPE(AOME_SET_STATIC_THRESHOLD, unsigned int)
#define AOM_CTRL_AOME_SET_STATIC_THRESHOLD

AOM_CTRL_USE_TYPE(AOME_SET_ARNR_MAXFRAMES, unsigned int)
#define AOM_CTRL_AOME_SET_ARNR_MAXFRAMES

AOM_CTRL_USE_TYPE(AOME_SET_ARNR_STRENGTH, unsigned int)
#define AOM_CTRL_AOME_SET_ARNR_STRENGTH

AOM_CTRL_USE_TYPE(AOME_SET_TUNING, int) /* aom_tune_metric */
#define AOM_CTRL_AOME_SET_TUNING

AOM_CTRL_USE_TYPE(AOME_SET_CQ_LEVEL, unsigned int)
#define AOM_CTRL_AOME_SET_CQ_LEVEL

AOM_CTRL_USE_TYPE(AV1E_SET_ROW_MT, unsigned int)
#define AOM_CTRL_AV1E_SET_ROW_MT

AOM_CTRL_USE_TYPE(AV1E_SET_TILE_COLUMNS, unsigned int)
#define AOM_CTRL_AV1E_SET_TILE_COLUMNS

AOM_CTRL_USE_TYPE(AV1E_SET_TILE_ROWS, unsigned int)
#define AOM_CTRL_AV1E_SET_TILE_ROWS

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_TPL_MODEL, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_TPL_MODEL

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_KEYFRAME_FILTERING, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_KEYFRAME_FILTERING

AOM_CTRL_USE_TYPE(AOME_GET_LAST_QUANTIZER, int *)
#define AOM_CTRL_AOME_GET_LAST_QUANTIZER

AOM_CTRL_USE_TYPE(AOME_GET_LAST_QUANTIZER_64, int *)
#define AOM_CTRL_AOME_GET_LAST_QUANTIZER_64

AOM_CTRL_USE_TYPE(AOME_SET_MAX_INTRA_BITRATE_PCT, unsigned int)
#define AOM_CTRL_AOME_SET_MAX_INTRA_BITRATE_PCT

AOM_CTRL_USE_TYPE(AOME_SET_MAX_INTER_BITRATE_PCT, unsigned int)
#define AOM_CTRL_AOME_SET_MAX_INTER_BITRATE_PCT

AOM_CTRL_USE_TYPE(AOME_SET_NUMBER_SPATIAL_LAYERS, int)
#define AOME_CTRL_AOME_SET_NUMBER_SPATIAL_LAYERS

AOM_CTRL_USE_TYPE(AV1E_SET_GF_CBR_BOOST_PCT, unsigned int)
#define AOM_CTRL_AV1E_SET_GF_CBR_BOOST_PCT

AOM_CTRL_USE_TYPE(AV1E_SET_LOSSLESS, unsigned int)
#define AOM_CTRL_AV1E_SET_LOSSLESS

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_CDEF, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_CDEF

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_RESTORATION, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_RESTORATION

AOM_CTRL_USE_TYPE(AV1E_SET_FORCE_VIDEO_MODE, unsigned int)
#define AOM_CTRL_AV1E_SET_FORCE_VIDEO_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_OBMC, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_OBMC

AOM_CTRL_USE_TYPE(AV1E_SET_DISABLE_TRELLIS_QUANT, unsigned int)
#define AOM_CTRL_AV1E_SET_DISABLE_TRELLIS_QUANT

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_QM, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_QM

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DIST_8X8, unsigned int)
#define AOM_CTRL_AV1E_SET_ENABLE_DIST_8X8

AOM_CTRL_USE_TYPE(AV1E_SET_QM_MIN, unsigned int)
#define AOM_CTRL_AV1E_SET_QM_MIN

AOM_CTRL_USE_TYPE(AV1E_SET_QM_MAX, unsigned int)
#define AOM_CTRL_AV1E_SET_QM_MAX

AOM_CTRL_USE_TYPE(AV1E_SET_QM_Y, unsigned int)
#define AOM_CTRL_AV1E_SET_QM_Y

AOM_CTRL_USE_TYPE(AV1E_SET_QM_U, unsigned int)
#define AOM_CTRL_AV1E_SET_QM_U

AOM_CTRL_USE_TYPE(AV1E_SET_QM_V, unsigned int)
#define AOM_CTRL_AV1E_SET_QM_V

AOM_CTRL_USE_TYPE(AV1E_SET_NUM_TG, unsigned int)
#define AOM_CTRL_AV1E_SET_NUM_TG

AOM_CTRL_USE_TYPE(AV1E_SET_MTU, unsigned int)
#define AOM_CTRL_AV1E_SET_MTU

AOM_CTRL_USE_TYPE(AV1E_SET_TIMING_INFO_TYPE, int) /* aom_timing_info_type_t */
#define AOM_CTRL_AV1E_SET_TIMING_INFO_TYPE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_RECT_PARTITIONS, int)
#define AOM_CTRL_AV1E_SET_ENABLE_RECT_PARTITIONS

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_AB_PARTITIONS, int)
#define AOM_CTRL_AV1E_SET_ENABLE_AB_PARTITIONS

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_1TO4_PARTITIONS, int)
#define AOM_CTRL_AV1E_SET_ENABLE_1TO4_PARTITIONS

AOM_CTRL_USE_TYPE(AV1E_SET_MIN_PARTITION_SIZE, int)
#define AOM_CTRL_AV1E_SET_MIN_PARTITION_SIZE

AOM_CTRL_USE_TYPE(AV1E_SET_MAX_PARTITION_SIZE, int)
#define AOM_CTRL_AV1E_SET_MAX_PARTITION_SIZE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_INTRA_EDGE_FILTER, int)
#define AOM_CTRL_AV1E_SET_ENABLE_INTRA_EDGE_FILTER

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_ORDER_HINT, int)
#define AOM_CTRL_AV1E_SET_ENABLE_ORDER_HINT

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_TX64, int)
#define AOM_CTRL_AV1E_SET_ENABLE_TX64

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_FLIP_IDTX, int)
#define AOM_CTRL_AV1E_SET_ENABLE_FLIP_IDTX

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_RECT_TX, int)
#define AOM_CTRL_AV1E_SET_ENABLE_RECT_TX

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DIST_WTD_COMP, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DIST_WTD_COMP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_REF_FRAME_MVS, int)
#define AOM_CTRL_AV1E_SET_ENABLE_REF_FRAME_MVS

AOM_CTRL_USE_TYPE(AV1E_SET_ALLOW_REF_FRAME_MVS, int)
#define AOM_CTRL_AV1E_SET_ALLOW_REF_FRAME_MVS

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DUAL_FILTER, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DUAL_FILTER

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_CHROMA_DELTAQ, int)
#define AOM_CTRL_AV1E_SET_ENABLE_CHROMA_DELTAQ

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_MASKED_COMP, int)
#define AOM_CTRL_AV1E_SET_ENABLE_MASKED_COMP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_ONESIDED_COMP, int)
#define AOM_CTRL_AV1E_SET_ENABLE_ONESIDED_COMP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_INTERINTRA_COMP, int)
#define AOM_CTRL_AV1E_SET_ENABLE_INTERINTRA_COMP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_SMOOTH_INTERINTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_SMOOTH_INTERINTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DIFF_WTD_COMP, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DIFF_WTD_COMP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_INTERINTER_WEDGE, int)
#define AOM_CTRL_AV1E_SET_ENABLE_INTERINTER_WEDGE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_INTERINTRA_WEDGE, int)
#define AOM_CTRL_AV1E_SET_ENABLE_INTERINTRA_WEDGE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_GLOBAL_MOTION, int)
#define AOM_CTRL_AV1E_SET_ENABLE_GLOBAL_MOTION

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_WARPED_MOTION, int)
#define AOM_CTRL_AV1E_SET_ENABLE_WARPED_MOTION

AOM_CTRL_USE_TYPE(AV1E_SET_ALLOW_WARPED_MOTION, int)
#define AOM_CTRL_AV1E_SET_ALLOW_WARPED_MOTION

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_FILTER_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_FILTER_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_SMOOTH_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_SMOOTH_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_PAETH_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_PAETH_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_CFL_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_CFL_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DIAGONAL_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DIAGONAL_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_SUPERRES, int)
#define AOM_CTRL_AV1E_SET_ENABLE_SUPERRES

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_OVERLAY, int)
#define AOM_CTRL_AV1E_SET_ENABLE_OVERLAY

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_PALETTE, int)
#define AOM_CTRL_AV1E_SET_ENABLE_PALETTE

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_INTRABC, int)
#define AOM_CTRL_AV1E_SET_ENABLE_INTRABC

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_ANGLE_DELTA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_ANGLE_DELTA

AOM_CTRL_USE_TYPE(AV1E_SET_FRAME_PARALLEL_DECODING, unsigned int)
#define AOM_CTRL_AV1E_SET_FRAME_PARALLEL_DECODING

AOM_CTRL_USE_TYPE(AV1E_SET_ERROR_RESILIENT_MODE, int)
#define AOM_CTRL_AV1E_SET_ERROR_RESILIENT_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_S_FRAME_MODE, int)
#define AOM_CTRL_AV1E_SET_S_FRAME_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_AQ_MODE, unsigned int)
#define AOM_CTRL_AV1E_SET_AQ_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_DELTAQ_MODE, unsigned int)
#define AOM_CTRL_AV1E_SET_DELTAQ_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_DELTAQ_STRENGTH, unsigned int)
#define AOM_CTRL_AV1E_SET_DELTAQ_STRENGTH

AOM_CTRL_USE_TYPE(AV1E_SET_DELTALF_MODE, unsigned int)
#define AOM_CTRL_AV1E_SET_DELTALF_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_FRAME_PERIODIC_BOOST, unsigned int)
#define AOM_CTRL_AV1E_SET_FRAME_PERIODIC_BOOST

AOM_CTRL_USE_TYPE(AV1E_SET_NOISE_SENSITIVITY, unsigned int)
#define AOM_CTRL_AV1E_SET_NOISE_SENSITIVITY

AOM_CTRL_USE_TYPE(AV1E_SET_TUNE_CONTENT, int) /* aom_tune_content */
#define AOM_CTRL_AV1E_SET_TUNE_CONTENT

AOM_CTRL_USE_TYPE(AV1E_SET_COLOR_PRIMARIES, int)
#define AOM_CTRL_AV1E_SET_COLOR_PRIMARIES

AOM_CTRL_USE_TYPE(AV1E_SET_TRANSFER_CHARACTERISTICS, int)
#define AOM_CTRL_AV1E_SET_TRANSFER_CHARACTERISTICS

AOM_CTRL_USE_TYPE(AV1E_SET_MATRIX_COEFFICIENTS, int)
#define AOM_CTRL_AV1E_SET_MATRIX_COEFFICIENTS

AOM_CTRL_USE_TYPE(AV1E_SET_CHROMA_SAMPLE_POSITION, int)
#define AOM_CTRL_AV1E_SET_CHROMA_SAMPLE_POSITION

AOM_CTRL_USE_TYPE(AV1E_SET_MIN_GF_INTERVAL, unsigned int)
#define AOM_CTRL_AV1E_SET_MIN_GF_INTERVAL

AOM_CTRL_USE_TYPE(AV1E_SET_MAX_GF_INTERVAL, unsigned int)
#define AOM_CTRL_AV1E_SET_MAX_GF_INTERVAL

AOM_CTRL_USE_TYPE(AV1E_GET_ACTIVEMAP, aom_active_map_t *)
#define AOM_CTRL_AV1E_GET_ACTIVEMAP

AOM_CTRL_USE_TYPE(AV1E_SET_COLOR_RANGE, int)
#define AOM_CTRL_AV1E_SET_COLOR_RANGE

#define AOM_CTRL_AV1E_SET_RENDER_SIZE
AOM_CTRL_USE_TYPE(AV1E_SET_RENDER_SIZE, int *)

AOM_CTRL_USE_TYPE(AV1E_SET_SUPERBLOCK_SIZE, unsigned int)
#define AOM_CTRL_AV1E_SET_SUPERBLOCK_SIZE

AOM_CTRL_USE_TYPE(AV1E_GET_SEQ_LEVEL_IDX, int *)
#define AOM_CTRL_AV1E_GET_SEQ_LEVEL_IDX

AOM_CTRL_USE_TYPE(AV1E_GET_BASELINE_GF_INTERVAL, int *)
#define AOM_CTRL_AV1E_GET_BASELINE_GF_INTERVAL

AOM_CTRL_USE_TYPE(AV1E_SET_SINGLE_TILE_DECODING, unsigned int)
#define AOM_CTRL_AV1E_SET_SINGLE_TILE_DECODING

AOM_CTRL_USE_TYPE(AV1E_ENABLE_MOTION_VECTOR_UNIT_TEST, unsigned int)
#define AOM_CTRL_AV1E_ENABLE_MOTION_VECTOR_UNIT_TEST

AOM_CTRL_USE_TYPE(AV1E_ENABLE_EXT_TILE_DEBUG, unsigned int)
#define AOM_CTRL_AV1E_ENABLE_EXT_TILE_DEBUG

AOM_CTRL_USE_TYPE(AV1E_SET_VMAF_MODEL_PATH, const char *)
#define AOM_CTRL_AV1E_SET_VMAF_MODEL_PATH

AOM_CTRL_USE_TYPE(AV1E_SET_FILM_GRAIN_TEST_VECTOR, int)
#define AOM_CTRL_AV1E_SET_FILM_GRAIN_TEST_VECTOR

AOM_CTRL_USE_TYPE(AV1E_SET_FILM_GRAIN_TABLE, const char *)
#define AOM_CTRL_AV1E_SET_FILM_GRAIN_TABLE

AOM_CTRL_USE_TYPE(AV1E_SET_CDF_UPDATE_MODE, unsigned int)
#define AOM_CTRL_AV1E_SET_CDF_UPDATE_MODE

AOM_CTRL_USE_TYPE(AV1E_SET_DENOISE_NOISE_LEVEL, int)
#define AOM_CTRL_AV1E_SET_DENOISE_NOISE_LEVEL

AOM_CTRL_USE_TYPE(AV1E_SET_DENOISE_BLOCK_SIZE, unsigned int)
#define AOM_CTRL_AV1E_SET_DENOISE_BLOCK_SIZE

AOM_CTRL_USE_TYPE(AV1E_SET_CHROMA_SUBSAMPLING_X, unsigned int)
#define AOM_CTRL_AV1E_SET_CHROMA_SUBSAMPLING_X

AOM_CTRL_USE_TYPE(AV1E_SET_CHROMA_SUBSAMPLING_Y, unsigned int)
#define AOM_CTRL_AV1E_SET_CHROMA_SUBSAMPLING_Y

AOM_CTRL_USE_TYPE(AV1E_SET_REDUCED_TX_TYPE_SET, int)
#define AOM_CTRL_AV1E_SET_REDUCED_TX_TYPE_SET

AOM_CTRL_USE_TYPE(AV1E_SET_INTRA_DCT_ONLY, int)
#define AOM_CTRL_AV1E_SET_INTRA_DCT_ONLY

AOM_CTRL_USE_TYPE(AV1E_SET_INTER_DCT_ONLY, int)
#define AOM_CTRL_AV1E_SET_INTER_DCT_ONLY

AOM_CTRL_USE_TYPE(AV1E_SET_INTRA_DEFAULT_TX_ONLY, int)
#define AOM_CTRL_AV1E_SET_INTRA_DEFAULT_TX_ONLY

AOM_CTRL_USE_TYPE(AV1E_SET_QUANT_B_ADAPT, int)
#define AOM_CTRL_AV1E_SET_QUANT_B_ADAPT

AOM_CTRL_USE_TYPE(AV1E_SET_GF_MIN_PYRAMID_HEIGHT, unsigned int)
#define AOM_CTRL_AV1E_SET_GF_MIN_PYRAMID_HEIGHT

AOM_CTRL_USE_TYPE(AV1E_SET_GF_MAX_PYRAMID_HEIGHT, unsigned int)
#define AOM_CTRL_AV1E_SET_GF_MAX_PYRAMID_HEIGHT

AOM_CTRL_USE_TYPE(AV1E_SET_MAX_REFERENCE_FRAMES, int)
#define AOM_CTRL_AV1E_SET_MAX_REFERENCE_FRAMES

AOM_CTRL_USE_TYPE(AV1E_SET_REDUCED_REFERENCE_SET, int)
#define AOM_CTRL_AV1E_SET_REDUCED_REFERENCE_SET

AOM_CTRL_USE_TYPE(AV1E_SET_COEFF_COST_UPD_FREQ, unsigned int)
#define AOM_CTRL_AV1E_SET_COEFF_COST_UPD_FREQ

AOM_CTRL_USE_TYPE(AV1E_SET_MODE_COST_UPD_FREQ, unsigned int)
#define AOM_CTRL_AV1E_SET_MODE_COST_UPD_FREQ

AOM_CTRL_USE_TYPE(AV1E_SET_MV_COST_UPD_FREQ, unsigned int)
#define AOM_CTRL_AV1E_SET_MV_COST_UPD_FREQ

AOM_CTRL_USE_TYPE(AV1E_SET_TARGET_SEQ_LEVEL_IDX, int)
#define AOM_CTRL_AV1E_SET_TARGET_SEQ_LEVEL_IDX

AOM_CTRL_USE_TYPE(AV1E_SET_TIER_MASK, unsigned int)
#define AOM_CTRL_AV1E_SET_TIER_MASK

AOM_CTRL_USE_TYPE(AV1E_SET_MIN_CR, unsigned int)
#define AOM_CTRL_AV1E_SET_MIN_CR

AOM_CTRL_USE_TYPE(AV1E_SET_SVC_LAYER_ID, aom_svc_layer_id_t *)
#define AOME_CTRL_AV1E_SET_SVC_LAYER_ID

AOM_CTRL_USE_TYPE(AV1E_SET_SVC_PARAMS, aom_svc_params_t *)
#define AOME_CTRL_AV1E_SET_SVC_PARAMS

AOM_CTRL_USE_TYPE(AV1E_SET_SVC_REF_FRAME_CONFIG, aom_svc_ref_frame_config_t *)
#define AOME_CTRL_AV1E_SET_SVC_REF_FRAME_CONFIG

AOM_CTRL_USE_TYPE(AV1E_ENABLE_SB_MULTIPASS_UNIT_TEST, unsigned int)
#define AOM_CTRL_AV1E_ENABLE_SB_MULTIPASS_UNIT_TEST

AOM_CTRL_USE_TYPE(AV1E_SET_VBR_CORPUS_COMPLEXITY_LAP, unsigned int)
#define AOM_CTRL_AV1E_SET_VBR_CORPUS_COMPLEXITY_LAP

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DNL_DENOISING, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DNL_DENOISING

AOM_CTRL_USE_TYPE(AV1E_SET_DV_COST_UPD_FREQ, unsigned int)
#define AOM_CTRL_AV1E_SET_DV_COST_UPD_FREQ

AOM_CTRL_USE_TYPE(AV1E_SET_PARTITION_INFO_PATH, const char *)
#define AOM_CTRL_AV1E_SET_PARTITION_INFO_PATH

AOM_CTRL_USE_TYPE(AV1E_SET_EXTERNAL_PARTITION, aom_ext_part_funcs_t *)
#define AOM_CTRL_AV1E_SET_EXTERNAL_PARTITION

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_DIRECTIONAL_INTRA, int)
#define AOM_CTRL_AV1E_SET_ENABLE_DIRECTIONAL_INTRA

AOM_CTRL_USE_TYPE(AV1E_SET_ENABLE_TX_SIZE_SEARCH, int)
#define AOM_CTRL_AV1E_SET_ENABLE_TX_SIZE_SEARCH

AOM_CTRL_USE_TYPE(AV1E_SET_SVC_REF_FRAME_COMP_PRED,
                  aom_svc_ref_frame_comp_pred_t *)
#define AOME_CTRL_AV1E_SET_SVC_REF_FRAME_COMP_PRED

AOM_CTRL_USE_TYPE(AV1E_SET_LOOPFILTER_CONTROL, int)
#define AOM_CTRL_AV1E_SET_LOOPFILTER_CONTROL

AOM_CTRL_USE_TYPE(AOME_GET_LOOPFILTER_LEVEL, int *)
#define AOM_CTRL_AOME_GET_LOOPFILTER_LEVEL

AOM_CTRL_USE_TYPE(AV1E_SET_AUTO_INTRA_TOOLS_OFF, unsigned int)
#define AOM_CTRL_AV1E_SET_AUTO_INTRA_TOOLS_OFF

AOM_CTRL_USE_TYPE(AV1E_SET_RTC_EXTERNAL_RC, int)
#define AOM_CTRL_AV1E_SET_RTC_EXTERNAL_RC

AOM_CTRL_USE_TYPE(AV1E_SET_FP_MT, unsigned int)
#define AOM_CTRL_AV1E_SET_FP_MT

AOM_CTRL_USE_TYPE(AV1E_SET_FP_MT_UNIT_TEST, unsigned int)
#define AOM_CTRL_AV1E_SET_FP_MT_UNIT_TEST

/*!\endcond */
/*! @} - end defgroup aom_encoder */
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AOM_AOMCX_H_
