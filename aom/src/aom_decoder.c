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
#include <string.h>
#include "aom/internal/aom_codec_internal.h"

#define SAVE_STATUS(ctx, var) (ctx ? (ctx->err = var) : var)

static AomCodecAlgPrivT *get_alg_priv(AomCodecCtxT *ctx) {
  return (AomCodecAlgPrivT *)ctx->priv;
}

AomCodecErrT aom_codec_dec_init_ver(AomCodecCtxT *ctx, AomCodecIfaceT *iface,
                                    const AomCodecDecCfgT *cfg,
                                    AomCodecFlagsT flags, int ver) {
  AomCodecErrT res;

  if (ver != AOM_DECODER_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if (!ctx || !iface)
    res = AOM_CODEC_INVALID_PARAM;
  else if (iface->abi_version != AOM_CODEC_INTERNAL_ABI_VERSION)
    res = AOM_CODEC_ABI_MISMATCH;
  else if ((flags & AOM_CODEC_USE_POSTPROC) &&
           !(iface->caps & AOM_CODEC_CAP_POSTPROC))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_ERROR_CONCEALMENT) &&
           !(iface->caps & AOM_CODEC_CAP_ERROR_CONCEALMENT))
    res = AOM_CODEC_INCAPABLE;
  else if ((flags & AOM_CODEC_USE_INPUT_FRAGMENTS) &&
           !(iface->caps & AOM_CODEC_CAP_INPUT_FRAGMENTS))
    res = AOM_CODEC_INCAPABLE;
  else if (!(iface->caps & AOM_CODEC_CAP_DECODER))
    res = AOM_CODEC_INCAPABLE;
  else {
    memset(ctx, 0, sizeof(*ctx));
    ctx->iface = iface;
    ctx->name = iface->name;
    ctx->priv = NULL;
    ctx->init_flags = flags;
    ctx->config.dec = cfg;

    res = ctx->iface->init(ctx, NULL);
    if (res) {
      ctx->err_detail = ctx->priv ? ctx->priv->err_detail : NULL;
      aom_codec_destroy(ctx);
    }
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_peek_stream_info(AomCodecIfaceT *iface,
                                        const uint8_t *data,
                                        unsigned int data_sz,
                                        AomCodecStreamInfoT *si) {
  AomCodecErrT res;

  if (!iface || !data || !data_sz || !si) {
    res = AOM_CODEC_INVALID_PARAM;
  } else {
    /* Set default/unknown values */
    si->w = 0;
    si->h = 0;

    res = iface->dec.peek_si(data, data_sz, si);
  }

  return res;
}

AomCodecErrT aom_codec_get_stream_info(AomCodecCtxT *ctx,
                                       AomCodecStreamInfoT *si) {
  AomCodecErrT res;

  if (!ctx || !si) {
    res = AOM_CODEC_INVALID_PARAM;
  } else if (!ctx->iface || !ctx->priv) {
    res = AOM_CODEC_ERROR;
  } else {
    /* Set default/unknown values */
    si->w = 0;
    si->h = 0;

    res = ctx->iface->dec.get_si(get_alg_priv(ctx), si);
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_decode(AomCodecCtxT *ctx, const uint8_t *data,
                              unsigned int data_sz, void *user_priv,
                              long deadline) {
  AomCodecErrT res;

  /* Sanity checks */
  /* NULL data ptr allowed if data_sz is 0 too */
  if (!ctx || (!data && data_sz) || (data && !data_sz))
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv)
    res = AOM_CODEC_ERROR;
  else {
    res = ctx->iface->dec.decode(get_alg_priv(ctx), data, data_sz, user_priv,
                                 deadline);
  }

  return SAVE_STATUS(ctx, res);
}

AomImageT *aom_codec_get_frame(AomCodecCtxT *ctx, AomCodecIterT *iter) {
  AomImageT *img;

  if (!ctx || !iter || !ctx->iface || !ctx->priv)
    img = NULL;
  else
    img = ctx->iface->dec.get_frame(get_alg_priv(ctx), iter);

  return img;
}

AomCodecErrT aom_codec_register_put_frame_cb(AomCodecCtxT *ctx,
                                             AomCodecPutFrameCbFnT cb,
                                             void *user_priv) {
  AomCodecErrT res;

  if (!ctx || !cb)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv ||
           !(ctx->iface->caps & AOM_CODEC_CAP_PUT_FRAME))
    res = AOM_CODEC_ERROR;
  else {
    ctx->priv->dec.put_frame_cb.u.put_frame = cb;
    ctx->priv->dec.put_frame_cb.user_priv = user_priv;
    res = AOM_CODEC_OK;
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_register_put_slice_cb(AomCodecCtxT *ctx,
                                             AomCodecPutSliceCbFnT cb,
                                             void *user_priv) {
  AomCodecErrT res;

  if (!ctx || !cb)
    res = AOM_CODEC_INVALID_PARAM;
  else if (!ctx->iface || !ctx->priv ||
           !(ctx->iface->caps & AOM_CODEC_CAP_PUT_SLICE))
    res = AOM_CODEC_ERROR;
  else {
    ctx->priv->dec.put_slice_cb.u.put_slice = cb;
    ctx->priv->dec.put_slice_cb.user_priv = user_priv;
    res = AOM_CODEC_OK;
  }

  return SAVE_STATUS(ctx, res);
}

AomCodecErrT aom_codec_set_frame_buffer_functions(
    AomCodecCtxT *ctx, AomGetFrameBufferCbFnT cb_get,
    AomReleaseFrameBufferCbFnT cb_release, void *cb_priv) {
  AomCodecErrT res;

  if (!ctx || !cb_get || !cb_release) {
    res = AOM_CODEC_INVALID_PARAM;
  } else if (!ctx->iface || !ctx->priv ||
             !(ctx->iface->caps & AOM_CODEC_CAP_EXTERNAL_FRAME_BUFFER)) {
    res = AOM_CODEC_ERROR;
  } else {
    res = ctx->iface->dec.set_fb_fn(get_alg_priv(ctx), cb_get, cb_release,
                                    cb_priv);
  }

  return SAVE_STATUS(ctx, res);
}
