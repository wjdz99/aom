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

/*!\file
 * \brief Provides the high level interface to wrap decoder algorithms.
 *
 */
#include <stdarg.h>
#include <stdlib.h>

#include "config/aom_config.h"
#include "config/aom_version.h"

#include "aom/aom_integer.h"
#include "aom/aomcx.h"
#include "aom/aomdx.h"
#include "aom/internal/aom_codec_internal.h"

#define SAVE_STATUS(ctx, var) (ctx ? (ctx->err = var) : var)

int aom_codec_version(void) { return VERSION_PACKED; }

const char *aom_codec_version_str(void) { return VERSION_STRING_NOSP; }

const char *aom_codec_version_extra_str(void) { return VERSION_EXTRA; }

const char *aom_codec_iface_name(aom_codec_iface_t *iface) {
  return iface ? iface->name : "<invalid interface>";
}

const char *aom_codec_err_to_string(aom_codec_err_t err) {
  switch (err) {
    case AOM_CODEC_OK: return "Success";
    case AOM_CODEC_ERROR: return "Unspecified internal error";
    case AOM_CODEC_MEM_ERROR: return "Memory allocation error";
    case AOM_CODEC_ABI_MISMATCH: return "ABI version mismatch";
    case AOM_CODEC_INCAPABLE:
      return "Codec does not implement requested capability";
    case AOM_CODEC_UNSUP_BITSTREAM:
      return "Bitstream not supported by this decoder";
    case AOM_CODEC_UNSUP_FEATURE:
      return "Bitstream required feature not supported by this decoder";
    case AOM_CODEC_CORRUPT_FRAME: return "Corrupt frame detected";
    case AOM_CODEC_INVALID_PARAM: return "Invalid parameter";
    case AOM_CODEC_LIST_END: return "End of iterated list";
  }

  return "Unrecognized error code";
}

const char *aom_codec_error(aom_codec_ctx_t *ctx) {
  return (ctx) ? aom_codec_err_to_string(ctx->err)
               : aom_codec_err_to_string(AOM_CODEC_INVALID_PARAM);
}

const char *aom_codec_error_detail(aom_codec_ctx_t *ctx) {
  if (ctx && ctx->err)
    return ctx->priv ? ctx->priv->err_detail : ctx->err_detail;

  return NULL;
}

aom_codec_err_t aom_codec_destroy(aom_codec_ctx_t *ctx) {
  aom_codec_err_t res;

  if (!ctx)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv)
    res = AOM_CODEC_ERROR;
  else {
    ctx->iface->destroy((aom_codec_alg_priv_t *)ctx->priv);

    ctx->iface = NULL;
    ctx->name = NULL;
    ctx->priv = NULL;
    res = AOM_CODEC_OK;
  }

  return SAVE_STATUS(ctx, res);
}

aom_codec_caps_t aom_codec_get_caps(aom_codec_iface_t *iface) {
  return (iface) ? iface->caps : 0;
}

aom_codec_err_t aom_codec_control_set_int(aom_codec_ctx_t *ctx, int ctrl_id,
                                          int val) {
  switch (ctrl_id) {
    case AOME_SET_CPUUSED:
    case AOME_SET_NUMBER_SPATIAL_LAYERS:
    case AOME_SET_TUNING:
    case AOM_SET_DBG_COLOR_B_MODES:
    case AOM_SET_DBG_COLOR_MB_MODES:
    case AOM_SET_DBG_COLOR_REF_FRAME:
    case AOM_SET_DBG_DISPLAY_MV:
    case AV1D_SET_OPERATING_POINT:
    case AV1D_SET_OUTPUT_ALL_LAYERS:
    case AV1D_SET_SKIP_FILM_GRAIN:
    case AV1E_SET_ALLOW_REF_FRAME_MVS:
    case AV1E_SET_ALLOW_WARPED_MOTION:
    case AV1E_SET_CHROMA_SAMPLE_POSITION:
    case AV1E_SET_COLOR_PRIMARIES:
    case AV1E_SET_COLOR_RANGE:
    case AV1E_SET_DENOISE_NOISE_LEVEL:
    case AV1E_SET_ENABLE_1TO4_PARTITIONS:
    case AV1E_SET_ENABLE_AB_PARTITIONS:
    case AV1E_SET_ENABLE_ANGLE_DELTA:
    case AV1E_SET_ENABLE_CFL_INTRA:
    case AV1E_SET_ENABLE_CHROMA_DELTAQ:
    case AV1E_SET_ENABLE_DIFF_WTD_COMP:
    case AV1E_SET_ENABLE_DIST_WTD_COMP:
    case AV1E_SET_ENABLE_DUAL_FILTER:
    case AV1E_SET_ENABLE_FILTER_INTRA:
    case AV1E_SET_ENABLE_FLIP_IDTX:
    case AV1E_SET_ENABLE_GLOBAL_MOTION:
    case AV1E_SET_ENABLE_INTERINTER_WEDGE:
    case AV1E_SET_ENABLE_INTERINTRA_COMP:
    case AV1E_SET_ENABLE_INTERINTRA_WEDGE:
    case AV1E_SET_ENABLE_INTRABC:
    case AV1E_SET_ENABLE_INTRA_EDGE_FILTER:
    case AV1E_SET_ENABLE_MASKED_COMP:
    case AV1E_SET_ENABLE_ONESIDED_COMP:
    case AV1E_SET_ENABLE_ORDER_HINT:
    case AV1E_SET_ENABLE_OVERLAY:
    case AV1E_SET_ENABLE_PAETH_INTRA:
    case AV1E_SET_ENABLE_PALETTE:
    case AV1E_SET_ENABLE_RECT_PARTITIONS:
    case AV1E_SET_ENABLE_REF_FRAME_MVS:
    case AV1E_SET_ENABLE_SMOOTH_INTERINTRA:
    case AV1E_SET_ENABLE_SMOOTH_INTRA:
    case AV1E_SET_ENABLE_SUPERRES:
    case AV1E_SET_ENABLE_TX64:
    case AV1E_SET_ENABLE_WARPED_MOTION:
    case AV1E_SET_ERROR_RESILIENT_MODE:
    case AV1E_SET_FILM_GRAIN_TEST_VECTOR:
    case AV1E_SET_INTER_DCT_ONLY:
    case AV1E_SET_INTRA_DCT_ONLY:
    case AV1E_SET_INTRA_DEFAULT_TX_ONLY:
    case AV1E_SET_MATRIX_COEFFICIENTS:
    case AV1E_SET_MAX_PARTITION_SIZE:
    case AV1E_SET_MAX_REFERENCE_FRAMES:
    case AV1E_SET_MIN_PARTITION_SIZE:
    case AV1E_SET_QUANT_B_ADAPT:
    case AV1E_SET_REDUCED_REFERENCE_SET:
    case AV1E_SET_REDUCED_TX_TYPE_SET:
    case AV1E_SET_S_FRAME_MODE:
    case AV1E_SET_TARGET_SEQ_LEVEL_IDX:
    case AV1E_SET_TIMING_INFO_TYPE:
    case AV1E_SET_TRANSFER_CHARACTERISTICS:
    case AV1E_SET_TUNE_CONTENT:
    case AV1_SET_DECODE_TILE_COL:
    case AV1_SET_DECODE_TILE_ROW: return aom_codec_control(ctx, ctrl_id, val);
    default:
      if (ctx == NULL) {
        return AOM_CODEC_INVALID_PARAM;
      }
      ctx->err = AOM_CODEC_INVALID_PARAM;
      return ctx->err;
  }
}

aom_codec_err_t aom_codec_control_get_int(aom_codec_ctx_t *ctx, int ctrl_id,
                                          int *val) {
  switch (ctrl_id) {
    case AOMD_GET_FRAME_CORRUPTED:
    case AOMD_GET_LAST_REF_UPDATES:
    case AOMD_GET_LAST_REF_USED:
    case AOMD_GET_LAST_QUANTIZER:
    case AOME_GET_LAST_QUANTIZER:
    case AOME_GET_LAST_QUANTIZER_64:
    case AV1D_GET_DISPLAY_SIZE:
    case AV1D_GET_FRAME_SIZE:
    case AV1E_GET_SEQ_LEVEL_IDX: return aom_codec_control(ctx, ctrl_id, val);
    default:
      if (ctx == NULL) {
        return AOM_CODEC_INVALID_PARAM;
      }
      ctx->err = AOM_CODEC_INVALID_PARAM;
      return ctx->err;
  }
}

aom_codec_err_t aom_codec_control_set_uint(aom_codec_ctx_t *ctx, int ctrl_id,
                                           unsigned int val) {
  switch (ctrl_id) {
    case AOME_SET_ARNR_MAXFRAMES:
    case AOME_SET_ARNR_STRENGTH:
    case AOME_SET_CQ_LEVEL:
    case AOME_SET_ENABLEAUTOALTREF:
    case AOME_SET_ENABLEAUTOBWDREF:
    case AOME_SET_MAX_INTRA_BITRATE_PCT:
    case AOME_SET_SHARPNESS:
    case AOME_SET_SPATIAL_LAYER_ID:
    case AOME_SET_STATIC_THRESHOLD:
    case AV1D_SET_IS_ANNEXB:
    case AV1D_SET_ROW_MT:
    case AV1E_SET_AQ_MODE:
    case AV1E_SET_CDF_UPDATE_MODE:
    case AV1E_SET_CHROMA_SUBSAMPLING_X:
    case AV1E_SET_CHROMA_SUBSAMPLING_Y:
    case AV1E_SET_COEFF_COST_UPD_FREQ:
    case AV1E_SET_DELTALF_MODE:
    case AV1E_SET_DELTAQ_MODE:
    case AV1E_SET_DENOISE_BLOCK_SIZE:
    case AV1E_SET_DISABLE_TRELLIS_QUANT:
    case AV1E_SET_ENABLE_CDEF:
    case AV1E_SET_ENABLE_DIST_8X8:
    case AV1E_SET_ENABLE_KEYFRAME_FILTERING:
    case AV1E_SET_ENABLE_OBMC:
    case AV1E_SET_ENABLE_QM:
    case AV1E_SET_ENABLE_RESTORATION:
    case AV1E_SET_ENABLE_TPL_MODEL:
    case AV1E_SET_FORCE_VIDEO_MODE:
    case AV1E_SET_FRAME_PARALLEL_DECODING:
    case AV1E_SET_FRAME_PERIODIC_BOOST:
    case AV1E_SET_GF_CBR_BOOST_PCT:
    case AV1E_SET_GF_MAX_PYRAMID_HEIGHT:
    case AV1E_SET_GF_MIN_PYRAMID_HEIGHT:
    case AV1E_SET_LOSSLESS:
    case AV1E_SET_MAX_GF_INTERVAL:
    case AV1E_SET_MAX_INTER_BITRATE_PCT:
    case AV1E_SET_MIN_CR:
    case AV1E_SET_MIN_GF_INTERVAL:
    case AV1E_SET_MODE_COST_UPD_FREQ:
    case AV1E_SET_MTU:
    case AV1E_SET_MV_COST_UPD_FREQ:
    case AV1E_SET_NOISE_SENSITIVITY:
    case AV1E_SET_NUM_TG:
    case AV1E_SET_QM_MAX:
    case AV1E_SET_QM_MIN:
    case AV1E_SET_QM_U:
    case AV1E_SET_QM_V:
    case AV1E_SET_QM_Y:
    case AV1E_SET_ROW_MT:
    case AV1E_SET_SINGLE_TILE_DECODING:
    case AV1E_SET_SUPERBLOCK_SIZE:
    case AV1E_SET_TIER_MASK:
    case AV1E_SET_TILE_COLUMNS:
    case AV1E_SET_TILE_ROWS:
    case AV1E_SET_VBR_CORPUS_COMPLEXITY_LAP:
    case AV1_SET_TILE_MODE: return aom_codec_control(ctx, ctrl_id, val);
    default:
      if (ctx == NULL) {
        return AOM_CODEC_INVALID_PARAM;
      }
      ctx->err = AOM_CODEC_INVALID_PARAM;
      return ctx->err;
  }
}

aom_codec_err_t aom_codec_control_get_uint(aom_codec_ctx_t *ctx, int ctrl_id,
                                           unsigned int *val) {
  switch (ctrl_id) {
    case AV1D_GET_BIT_DEPTH:
    case AV1D_GET_TILE_COUNT:
    case AV1D_GET_TILE_SIZE: return aom_codec_control(ctx, ctrl_id, val);
    default:
      if (ctx == NULL) {
        return AOM_CODEC_INVALID_PARAM;
      }
      ctx->err = AOM_CODEC_INVALID_PARAM;
      return ctx->err;
  }
}

aom_codec_err_t aom_codec_control_set_ptr(aom_codec_ctx_t *ctx, int ctrl_id,
                                          void *val, size_t size) {
  size_t expected_size;
  switch (ctrl_id) {
    case AOME_SET_ACTIVEMAP: expected_size = sizeof(aom_active_map_t); break;
    case AOME_SET_ROI_MAP: expected_size = sizeof(aom_roi_map_t); break;
    case AOME_SET_SCALEMODE: expected_size = sizeof(aom_scaling_mode_t); break;
    case AV1D_SET_EXT_REF_PTR:
      expected_size = sizeof(av1_ext_ref_frame_t);
      break;
    case AV1E_SET_RENDER_SIZE: expected_size = sizeof(int32_t[2]); break;
    case AV1E_SET_SVC_LAYER_ID:
      expected_size = sizeof(aom_svc_layer_id_t);
      break;
    case AV1E_SET_SVC_PARAMS: expected_size = sizeof(aom_svc_params_t); break;
    case AV1E_SET_SVC_REF_FRAME_CONFIG:
      expected_size = sizeof(aom_svc_ref_frame_config_t);
      break;
    case AV1_SET_INSPECTION_CALLBACK:
      expected_size = sizeof(aom_inspect_init);
      break;
    case AV1_SET_REFERENCE: expected_size = sizeof(av1_ref_frame_t); break;
    default: expected_size = 0;
  }
  if (expected_size != 0 && expected_size == size) {
    return aom_codec_control(ctx, ctrl_id, val);
  }
  if (ctx == NULL) {
    return AOM_CODEC_INVALID_PARAM;
  }
  ctx->err = AOM_CODEC_INVALID_PARAM;
  return ctx->err;
}

aom_codec_err_t aom_codec_control_get_ptr(aom_codec_ctx_t *ctx, int ctrl_id,
                                          void *val, size_t size) {
  size_t expected_size;
  switch (ctrl_id) {
    case AV1D_GET_FRAME_HEADER_INFO:
      expected_size = sizeof(aom_tile_data);
      break;
    case AV1D_GET_IMG_FORMAT: expected_size = sizeof(aom_img_fmt_t); break;
    case AV1D_GET_TILE_DATA: expected_size = sizeof(aom_tile_data); break;
    case AV1_GET_ACCOUNTING: expected_size = sizeof(Accounting *); break;
    case AV1_GET_NEW_FRAME_IMAGE: expected_size = sizeof(aom_image_t); break;
    case AV1_GET_REFERENCE: expected_size = sizeof(av1_ref_frame_t); break;
    default: expected_size = 0;
  }
  if (expected_size != 0 && expected_size == size) {
    return aom_codec_control(ctx, ctrl_id, val);
  }
  if (ctx == NULL) {
    return AOM_CODEC_INVALID_PARAM;
  }
  ctx->err = AOM_CODEC_INVALID_PARAM;
  return ctx->err;
}

aom_codec_err_t aom_codec_control(aom_codec_ctx_t *ctx, int ctrl_id, ...) {
  aom_codec_err_t res;

  if (!ctx || !ctrl_id)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv || !ctx->iface->ctrl_maps)
    res = AOM_CODEC_ERROR;
  else {
    aom_codec_ctrl_fn_map_t *entry;

    res = AOM_CODEC_ERROR;

    for (entry = ctx->iface->ctrl_maps; entry && entry->fn; entry++) {
      if (!entry->ctrl_id || entry->ctrl_id == ctrl_id) {
        va_list ap;

        va_start(ap, ctrl_id);
        res = entry->fn((aom_codec_alg_priv_t *)ctx->priv, ap);
        va_end(ap);
        break;
      }
    }
  }

  return SAVE_STATUS(ctx, res);
}

void aom_internal_error(struct aom_internal_error_info *info,
                        aom_codec_err_t error, const char *fmt, ...) {
  va_list ap;

  info->error_code = error;
  info->has_detail = 0;

  if (fmt) {
    size_t sz = sizeof(info->detail);

    info->has_detail = 1;
    va_start(ap, fmt);
    vsnprintf(info->detail, sz - 1, fmt, ap);
    va_end(ap);
    info->detail[sz - 1] = '\0';
  }

  if (info->setjmp) longjmp(info->jmp, info->error_code);
}

void aom_merge_corrupted_flag(int *corrupted, int value) {
  *corrupted |= value;
}

const char *aom_obu_type_to_string(OBU_TYPE type) {
  switch (type) {
    case OBU_SEQUENCE_HEADER: return "OBU_SEQUENCE_HEADER";
    case OBU_TEMPORAL_DELIMITER: return "OBU_TEMPORAL_DELIMITER";
    case OBU_FRAME_HEADER: return "OBU_FRAME_HEADER";
    case OBU_REDUNDANT_FRAME_HEADER: return "OBU_REDUNDANT_FRAME_HEADER";
    case OBU_FRAME: return "OBU_FRAME";
    case OBU_TILE_GROUP: return "OBU_TILE_GROUP";
    case OBU_METADATA: return "OBU_METADATA";
    case OBU_TILE_LIST: return "OBU_TILE_LIST";
    case OBU_PADDING: return "OBU_PADDING";
    default: break;
  }
  return "<Invalid OBU Type>";
}
