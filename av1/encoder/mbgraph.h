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

#ifndef AV1_ENCODER_MBGRAPH_H_
#define AV1_ENCODER_MBGRAPH_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  struct {
    int err;
    union {
      IntMv mv;
      PredictionMode mode;
    } m;
  } ref[TOTAL_REFS_PER_FRAME];
} MbgraphMbStats;

typedef struct { MbgraphMbStats *mb_stats; } MbgraphFrameStats;

struct Av1Comp;

void av1_update_mbgraph_stats(struct Av1Comp *cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_MBGRAPH_H_
