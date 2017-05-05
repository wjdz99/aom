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

#ifndef AV1_DECODER_DETOKENIZE_H_
#define AV1_DECODER_DETOKENIZE_H_

#include "./aom_config.h"
#if !CONFIG_PVQ || CONFIG_VAR_TX
#include "av1/decoder/decoder.h"
#include "av1/common/scan.h"
#endif  // !CONFIG_PVQ

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_PALETTE
void av1_decode_palette_tokens(Macroblockd *const xd, int plane, AomReader *r);
#endif  // CONFIG_PALETTE

#if !CONFIG_PVQ || CONFIG_VAR_TX
int av1_decode_block_tokens(Av1Common *cm, Macroblockd *const xd, int plane,
                            const ScanOrder *sc, int x, int y, TxSize tx_size,
                            TxType tx_type, int16_t *max_scan_line,
                            AomReader *r, int seg_id);
#endif  // !CONFIG_PVQ
#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AV1_DECODER_DETOKENIZE_H_
