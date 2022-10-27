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
#include <assert.h>

#include "av1/encoder/cost.h"
#include "av1/common/entropy.h"

// round(-log2(i/256. + 1/512.) * (1 << AV1_PROB_COST_SHIFT)); i = 128~255.
const uint16_t av1_prob_cost[128] = {
  509, 503, 498, 492, 486, 481, 475, 470, 465, 459, 454, 448, 443, 438, 433,
  428, 422, 417, 412, 407, 402, 397, 392, 387, 383, 378, 373, 368, 364, 359,
  354, 349, 345, 340, 336, 331, 327, 322, 318, 313, 309, 305, 300, 296, 292,
  287, 283, 279, 275, 271, 266, 262, 258, 254, 250, 246, 242, 238, 234, 230,
  226, 222, 218, 214, 211, 207, 203, 199, 195, 192, 188, 184, 181, 177, 173,
  170, 166, 162, 159, 155, 152, 148, 145, 141, 138, 134, 131, 127, 124, 120,
  117, 114, 110, 107, 104, 100, 97,  94,  90,  87,  84,  81,  78,  74,  71,
  68,  65,  62,  59,  55,  52,  49,  46,  43,  40,  37,  34,  31,  28,  25,
  22,  19,  16,  13,  10,  7,   4,   1,
};

void av1_cost_tokens_from_cdf(int *costs, const aom_cdf_prob *cdf,
                              const int *inv_map) {
  int i;
  aom_cdf_prob prev_cdf = 0;
  for (i = 0;; ++i) {
    aom_cdf_prob p15 = AOM_ICDF(cdf[i]) - prev_cdf;
    p15 = (p15 < EC_MIN_PROB) ? EC_MIN_PROB : p15;
    prev_cdf = AOM_ICDF(cdf[i]);

    if (inv_map)
      costs[inv_map[i]] = av1_cost_symbol(p15);
    else
      costs[i] = av1_cost_symbol(p15);

    // Stop once we reach the end of the CDF
    if (cdf[i] == AOM_ICDF(CDF_PROB_TOP)) break;
  }
}
