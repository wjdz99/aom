/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_TRANSFORM_SEARCH_H_
#define AOM_AV1_ENCODER_TRANSFORM_SEARCH_H_

#include "aom/aom_integer.h"
#include "av1/encoder/block.h"
#include "av1/encoder/encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

void model_rd_from_sse(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                       BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                       int num_samples, int *rate, int64_t *dist);

void model_rd_for_sb(const AV1_COMP *const cpi, BLOCK_SIZE bsize, MACROBLOCK *x,
                     MACROBLOCKD *xd, int plane_from, int plane_to,
                     int *out_rate_sum, int64_t *out_dist_sum,
                     int *skip_txfm_sb, int64_t *skip_sse_sb, int *plane_rate,
                     int64_t *plane_sse, int64_t *plane_dist);

void model_rd_with_curvfit(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                           BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                           int num_samples, int *rate, int64_t *dist);

void model_rd_for_sb_with_curvfit(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                  MACROBLOCK *x, MACROBLOCKD *xd,
                                  int plane_from, int plane_to,
                                  int *out_rate_sum, int64_t *out_dist_sum,
                                  int *skip_txfm_sb, int64_t *skip_sse_sb,
                                  int *plane_rate, int64_t *plane_sse,
                                  int64_t *plane_dist);

void model_rd_with_surffit(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                           BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                           int num_samples, int *rate, int64_t *dist);

void model_rd_for_sb_with_surffit(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                                  MACROBLOCK *x, MACROBLOCKD *xd,
                                  int plane_from, int plane_to,
                                  int *out_rate_sum, int64_t *out_dist_sum,
                                  int *skip_txfm_sb, int64_t *skip_sse_sb,
                                  int *plane_rate, int64_t *plane_sse,
                                  int64_t *plane_dist);

void model_rd_with_dnn(const AV1_COMP *const cpi, const MACROBLOCK *const x,
                       BLOCK_SIZE plane_bsize, int plane, int64_t sse,
                       int num_samples, int *rate, int64_t *dist);

void model_rd_for_sb_with_dnn(const AV1_COMP *const cpi, BLOCK_SIZE bsize,
                              MACROBLOCK *x, MACROBLOCKD *xd, int plane_from,
                              int plane_to, int *out_rate_sum,
                              int64_t *out_dist_sum, int *skip_txfm_sb,
                              int64_t *skip_sse_sb, int *plane_rate,
                              int64_t *plane_sse, int64_t *plane_dist);

void get_2x2_normalized_sses_and_sads(
    const AV1_COMP *const cpi, BLOCK_SIZE tx_bsize, const uint8_t *const src,
    int src_stride, const uint8_t *const dst, int dst_stride,
    const int16_t *const src_diff, int diff_stride, double *const sse_norm_arr,
    double *const sad_norm_arr);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_TRANSFORM_SEARCH_H_
