/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef _CX_IFACE_HELPER_H_
#define _CX_IFACE_HELPER_H_
#include "aom/internal/aom_codec_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

AV1_COMP *getCPI(const aom_codec_alg_priv_t *priv);

#ifdef __cplusplus
}
#endif

#endif
