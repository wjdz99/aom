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

#ifndef AV1_MV_CDF_TABLES_H
#define AV1_MV_CDF_TABLES_H

#include "av1/common/entropy.h"

#if CONFIG_EC_MULTISYMBOL
static const aom_cdf_prob default_mv_joint_cdf[CDF_SIZE(MV_JOINTS)] = {
  4096, 11264, 19328, 32768, 0
};

static const aom_cdf_prob default_mv_class_cdf[CDF_SIZE(MV_CLASSES)] = {
  28672, 30976, 31858, 32320, 32551, 32656, 32740, 32757, 32762, 32767, 32768, 0
};

static const aom_cdf_prob
    default_mv_class0_fp_cdf[CLASS0_SIZE][CDF_SIZE(MV_FP_SIZE)] = {
      { 16384, 24576, 26624, 32768, 0 }, { 12288, 21248, 24128, 32768, 0 },
    };

static const aom_cdf_prob default_mv_fp_cdf[CDF_SIZE(MV_FP_SIZE)] = {
  8192, 17408, 21248, 32768, 0
};

#endif  // CONFIG_EC_MULTISYMBOL

#endif  // AV1_MV_CDF_TABLES_H
