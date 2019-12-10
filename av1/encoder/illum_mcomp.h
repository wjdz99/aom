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

#ifndef AOM_AV1_ENCODER_ILLUM_MCOMP_H_
#define AOM_AV1_ENCODER_ILLUM_MCOMP_H_

#include "av1/encoder/mcomp.h"

#ifdef __cplusplus
extern "C" {
#endif

int av1_illum_full_pixel_search(const struct AV1_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE bsize, MV *mvp_full, int step_param,
                                int method, int run_mesh_search,
                                int error_per_bit, int *cost_list,
                                const MV *ref_mv, int var_max, int rd,
                                int x_pos, int y_pos, int intra,
                                const search_site_config *cfg);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_MCOMP_H_
