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
 * \brief Provides the high level interface to wrap encoder algorithms.
 *
 */
#include "./aom_config.h"

#if HAVE_FEXCEPT
#define _GNU_SOURCE
#include <fenv.h>
#endif

#include <limits.h>
#include <string.h>
#include "aom/internal/aom_codec_internal.h"

#define SAVE_STATUS(ctx, var) (ctx ? (ctx->err = var) : var)

static AomCodecAlgPrivT *get_alg_priv(AomCodecCtxT *ctx) {
  return (AomCodecAlgPrivT *)ctx->priv;
}

AomCodecErrT aom_codec_enc_init_ver(AomCodecCtxT *ctx, AomCodecIfaceT *iface,
                                    const AomCodecEncCfgT *cfg,
                                    AomCodecFlagsT flags, int ver) {
  AomCodecErrT res;

  if (ver != AOM_ENCODER_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if (!ctx || !iface || !cfg)
    res = AOM_CODEC_INVALID_PARAM;
  else if (iface->abi_version != AOM_CODEC_INTERNAL_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if (!(iface->caps & AOM_CODEC_CAP_ENCODER))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_PSNR) && !(iface->caps & AOM_CODEC_CAP_PSNR))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_OUTPUT_PARTITION) &&
           !(iface->caps & AOM_CODEC_CAP_OUTPUT_PARTITION))
    res = AOM_CODEC_INCAPABLE;
  else {
    ctx->iface = iface;
    ctx->name = iface->name;
    ctx->priv = NULL;
    ctx->init_flags = flags;
    ctx->config.enc = cfg;
    res = ctx->iface->init(ctx, NULL);

    if (res) {
      ctx->err_detail = ctx->priv ? ctx->priv->err_detail : NULL;
      aom_codec_destroy(ctx);
    }
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_enc_init_multi_ver(AomCodecCtxT *ctx,
                                          AomCodecIfaceT *iface,
                                          AomCodecEncCfgT *cfg, int num_enc,
                                          AomCodecFlagsT flags,
                                          AomRationalT *dsf, int ver) {
  AomCodecErrT res = AOM_CODEC_OK;

  if (ver != AOM_ENCODER_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if (!ctx || !iface || !cfg || (num_enc > 16 || num_enc < 1))
    res = AOM_CODEC_INVALID_PARAM;
  else if (iface->abi_version != AOM_CODEC_INTERNAL_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if (!(iface->caps & AOM_CODEC_CAP_ENCODER))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_PSNR) && !(iface->caps & AOM_CODEC_CAP_PSNR))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_OUTPUT_PARTITION) &&
           !(iface->caps & AOM_CODEC_CAP_OUTPUT_PARTITION))
    res = AOM_CODEC_INCAPABLE;
  else {
    int i;
    void *mem_loc = NULL;

    if (!(res = iface->enc.mr_get_mem_loc(cfg, &mem_loc))) {
      for (i = 0; i < num_enc; i++) {
        AomCodecPrivEncMrCfgT mr_cfg;

        /* Validate down-sampling factor. */
        if (dsf->num < 1 || dsf->num > 4096 || dsf->den < 1 ||
            dsf->den > dsf->num) {
          res = AOM_CODEC_INVALID_PARAM;
          break;
        }

        mr_cfg.mr_low_res_mode_info = mem_loc;
        mr_cfg.mr_total_resolutions = num_enc;
        mr_cfg.mr_encoder_id = num_enc - 1 - i;
        mr_cfg.mr_down_sampling_factor.num = dsf->num;
        mr_cfg.mr_down_sampling_factor.den = dsf->den;

        /* Force Key-frame synchronization. Namely, encoder at higher
         * resolution always use the same frame_type chosen by the
         * lowest-resolution encoder.
         */
        if (mr_cfg.mr_encoder_id) cfg->kf_mode = AOM_KF_DISABLED;

        ctx->iface = iface;
        ctx->name = iface->name;
        ctx->priv = NULL;
        ctx->init_flags = flags;
        ctx->config.enc = cfg;
        res = ctx->iface->init(ctx, &mr_cfg);

        if (res) {
          const char *error_detail = ctx->priv ? ctx->priv->err_detail : NULL;
          /* Destroy current ctx */
          ctx->err_detail = error_detail;
          aom_codec_destroy(ctx);

          /* Destroy already allocated high-level ctx */
          while (i) {
            ctx--;
            ctx->err_detail = error_detail;
            aom_codec_destroy(ctx);
            i--;
          }
        }

        if (res) break;

        ctx++;
        cfg++;
        dsf++;
      }
      ctx--;
    }
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_enc_config_default(AomCodecIfaceT *iface,
                                          AomCodecEncCfgT *cfg,
                                          unsigned int usage) {
  AomCodecErrT res;
  AomCodecEncCfgMapT *map;
  int i;

  if (!iface || !cfg || usage > INT_MAX)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!(iface->caps & AOM_CODEC_CAP_ENCODER))
    res = AOM_CODEC_INCAPABLE;
  else {
    res = AOM_CODEC_INVALID_PARAM;

    for (i = 0; i < iface->enc.cfg_map_count; ++i) {
      map = iface->enc.cfg_maps + i;
      if (map->usage == (int)usage) {
        *cfg = map->cfg;
        cfg->g_usage = usage;
        res = AOM_CODEC_OK;
        break;
      }
    }
  }

  return res;
}

/* clang-format off */
#define FLOATING_POINT_BEGIN_SCOPE do {
#define FLOATING_POINT_END_SCOPE } while (0);
/* clang-format on */

#if ARCH_X86 || ARCH_X86_64
/* On X86, disable the x87 unit's internal 80 bit precision for better
 * consistency with the SSE unit's 64 bit precision.
 */
#include "aom_ports/x86.h"
#define FLOATING_POINT_SET_PRECISION \
  unsigned short x87_orig_mode = x87_set_double_precision();
#define FLOATING_POINT_RESTORE_PRECISION x87_set_control_word(x87_orig_mode);
#else
#define FLOATING_POINT_SET_PRECISION
#define FLOATING_POINT_RESTORE_PRECISION
#endif  // ARCH_X86 || ARCH_X86_64

#if HAVE_FEXCEPT && CONFIG_DEBUG
#define FLOATING_POINT_SET_EXCEPTIONS \
  const int float_excepts = feenableexcept(FE_DIVBYZERO);
#define FLOATING_POINT_RESTORE_EXCEPTIONS feenableexcept(float_excepts);
#else
#define FLOATING_POINT_SET_EXCEPTIONS
#define FLOATING_POINT_RESTORE_EXCEPTIONS
#endif  // HAVE_FEXCEPT && CONFIG_DEBUG

#define FLOATING_POINT_INIT    \
  FLOATING_POINT_BEGIN_SCOPE   \
  FLOATING_POINT_SET_PRECISION \
  FLOATING_POINT_SET_EXCEPTIONS

#define FLOATING_POINT_RESTORE      \
  FLOATING_POINT_RESTORE_EXCEPTIONS \
  FLOATING_POINT_RESTORE_PRECISION  \
  FLOATING_POINT_END_SCOPE

AomCodecErrT aom_codec_encode(AomCodecCtxT *ctx, const AomImageT *img,
                              AomCodecPtsT pts, unsigned long duration,
                              AomEncFrameFlagsT flags, unsigned long deadline) {
  AomCodecErrT res = AOM_CODEC_OK;

  if (!ctx || (img && !duration))
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv)
    res = AOM_CODEC_ERROR;
  else if (!(ctx->iface->caps & AOM_CODEC_CAP_ENCODER))
    res = AOM_CODEC_INCAPABLE;
  else {
    unsigned int num_enc = ctx->priv->enc.total_encoders;

    /* Execute in a normalized floating point environment, if the platform
     * requires it.
     */
    FLOATING_POINT_INIT

    if (num_enc == 1)
      res = ctx->iface->enc.encode(get_alg_priv(ctx), img, pts, duration, flags,
                                   deadline);
    else {
      /* Multi-resolution encoding:
       * Encode multi-levels in reverse order. For example,
       * if mr_total_resolutions = 3, first encode level 2,
       * then encode level 1, and finally encode level 0.
       */
      int i;

      ctx += num_enc - 1;
      if (img) img += num_enc - 1;

      for (i = num_enc - 1; i >= 0; i--) {
        if ((res = ctx->iface->enc.encode(get_alg_priv(ctx), img, pts, duration,
                                          flags, deadline)))
          break;

        ctx--;
        if (img) img--;
      }
      ctx++;
    }

    FLOATING_POINT_RESTORE
  }

  return SAVE_STATUS(ctx, res);
}

const AomCodecCxPktT *aom_codec_get_cx_data(AomCodecCtxT *ctx,
                                            AomCodecIterT *iter) {
  const AomCodecCxPktT *pkt = NULL;

  if (ctx) {
    if (!iter)
      ctx->err = AOM_CODEC_INVALID_PARAM;
    else if (!ctx->iface || !ctx->priv)
      ctx->err = AOM_CODEC_ERROR;
    else if (!(ctx->iface->caps & AOM_CODEC_CAP_ENCODER))
      ctx->err = AOM_CODEC_INCAPABLE;
    else
      pkt = ctx->iface->enc.get_cx_data(get_alg_priv(ctx), iter);
  }

  if (pkt && pkt->kind == AOM_CODEC_CX_FRAME_PKT) {
    // If the application has specified a destination area for the
    // compressed data, and the codec has not placed the data there,
    // and it fits, copy it.
    AomCodecPrivT *const priv = ctx->priv;
    char *const dst_buf = (char *)priv->enc.cx_data_dst_buf.buf;

    if (dst_buf && pkt->data.raw.buf != dst_buf &&
        pkt->data.raw.sz + priv->enc.cx_data_pad_before +
                priv->enc.cx_data_pad_after <=
            priv->enc.cx_data_dst_buf.sz) {
      AomCodecCxPktT *modified_pkt = &priv->enc.cx_data_pkt;

      memcpy(dst_buf + priv->enc.cx_data_pad_before, pkt->data.raw.buf,
             pkt->data.raw.sz);
      *modified_pkt = *pkt;
      modified_pkt->data.raw.buf = dst_buf;
      modified_pkt->data.raw.sz +=
          priv->enc.cx_data_pad_before + priv->enc.cx_data_pad_after;
      pkt = modified_pkt;
    }

    if (dst_buf == pkt->data.raw.buf) {
      priv->enc.cx_data_dst_buf.buf = dst_buf + pkt->data.raw.sz;
      priv->enc.cx_data_dst_buf.sz -= pkt->data.raw.sz;
    }
  }

  return pkt;
}

AomCodecErrT aom_codec_set_cx_data_buf(AomCodecCtxT *ctx,
                                       const AomFixedBufT *buf,
                                       unsigned int pad_before,
                                       unsigned int pad_after) {
  if (!ctx || !ctx->priv) return AOM_CODEC_INVALID_PARAM;

  if (buf) {
    ctx->priv->enc.cx_data_dst_buf = *buf;
    ctx->priv->enc.cx_data_pad_before = pad_before;
    ctx->priv->enc.cx_data_pad_after = pad_after;
  } else {
    ctx->priv->enc.cx_data_dst_buf.buf = NULL;
    ctx->priv->enc.cx_data_dst_buf.sz = 0;
    ctx->priv->enc.cx_data_pad_before = 0;
    ctx->priv->enc.cx_data_pad_after = 0;
  }

  return AOM_CODEC_OK;
}

const AomImageT *aom_codec_get_preview_frame(AomCodecCtxT *ctx) {
  AomImageT *img = NULL;

  if (ctx) {
    if (!ctx->iface || !ctx->priv)
      ctx->err = AOM_CODEC_ERROR;
    else if (!(ctx->iface->caps & AOM_CODEC_CAP_ENCODER))
      ctx->err = AOM_CODEC_INCAPABLE;
    else if (!ctx->iface->enc.get_preview)
      ctx->err = AOM_CODEC_INCAPABLE;
    else
      img = ctx->iface->enc.get_preview(get_alg_priv(ctx));
  }

  return img;
}

AomFixedBufT *aom_codec_get_global_headers(AomCodecCtxT *ctx) {
  AomFixedBufT *buf = NULL;

  if (ctx) {
    if (!ctx->iface || !ctx->priv)
      ctx->err = AOM_CODEC_ERROR;
    else if (!(ctx->iface->caps & AOM_CODEC_CAP_ENCODER))
      ctx->err = AOM_CODEC_INCAPABLE;
    else if (!ctx->iface->enc.get_glob_hdrs)
      ctx->err = AOM_CODEC_INCAPABLE;
    else
      buf = ctx->iface->enc.get_glob_hdrs(get_alg_priv(ctx));
  }

  return buf;
}

AomCodecErrT aom_codec_enc_config_set(AomCodecCtxT *ctx,
                                      const AomCodecEncCfgT *cfg) {
  AomCodecErrT res;

  if (!ctx || !ctx->iface || !ctx->priv || !cfg)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!(ctx->iface->caps & AOM_CODEC_CAP_ENCODER))
    res = AOM_CODEC_INCAPABLE;
  else
    res = ctx->iface->enc.cfg_set(get_alg_priv(ctx), cfg);

  return SAVE_STATUS(ctx, res);
}

int aom_codec_pkt_list_add(struct AomCodecPktList *list,
                           const struct AomCodecCxPkt *pkt) {
  if (list->cnt < list->max) {
    list->pkts[list->cnt++] = *pkt;
    return 0;
  }

  return 1;
}

const AomCodecCxPktT *aom_codec_pkt_list_get(struct AomCodecPktList *list,
                                             AomCodecIterT *iter) {
  const AomCodecCxPktT *pkt;

  if (!(*iter)) {
    *iter = list->pkts;
  }

  pkt = (const AomCodecCxPktT *)*iter;

  if ((size_t)(pkt - list->pkts) < list->cnt)
    *iter = pkt + 1;
  else
    pkt = NULL;

  return pkt;
}
