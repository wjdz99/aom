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

#ifndef VP10_ENCODER_SUBEXP_H_
#define VP10_ENCODER_SUBEXP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "aom_dsp/prob.h"

struct vpx_writer;

void av1_write_prob_diff_update(struct vpx_writer *w, vpx_prob newp,
                                 vpx_prob oldp);

void av1_cond_prob_diff_update(struct vpx_writer *w, vpx_prob *oldp,
                                const unsigned int ct[2]);

int av1_prob_diff_update_savings_search(const unsigned int *ct, vpx_prob oldp,
                                         vpx_prob *bestp, vpx_prob upd);

int av1_prob_diff_update_savings_search_model(const unsigned int *ct,
                                               const vpx_prob *oldp,
                                               vpx_prob *bestp, vpx_prob upd,
                                               int stepsize);

int av1_cond_prob_diff_update_savings(vpx_prob *oldp,
                                       const unsigned int ct[2]);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_SUBEXP_H_
