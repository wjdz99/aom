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

#ifndef RATE_HIST_H_
#define RATE_HIST_H_

#include "aom/aom_encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

struct rate_hist;

struct rate_hist *init_rate_histogram(const AomCodecEncCfgT *cfg,
                                      const AomRationalT *fps);

void destroy_rate_histogram(struct rate_hist *hist);

void update_rate_histogram(struct rate_hist *hist, const AomCodecEncCfgT *cfg,
                           const AomCodecCxPktT *pkt);

void show_q_histogram(const int counts[64], int max_buckets);

void show_rate_histogram(struct rate_hist *hist, const AomCodecEncCfgT *cfg,
                         int max_buckets);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // RATE_HIST_H_
