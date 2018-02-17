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

#ifndef AOM_DSP_BITWRITER_H_
#define AOM_DSP_BITWRITER_H_

#include <assert.h>
#include "./aom_config.h"

#include "aom_dsp/daalaboolwriter.h"
#include "aom_dsp/prob.h"

#if CONFIG_RD_DEBUG
#include "av1/common/blockd.h"
#include "av1/encoder/cost.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct daala_writer aom_writer;

typedef struct TOKEN_STATS {
  int cost;
#if CONFIG_RD_DEBUG
  int txb_coeff_cost_map[TXB_COEFF_COST_MAP_SIZE][TXB_COEFF_COST_MAP_SIZE];
#endif
} TOKEN_STATS;

static INLINE void init_token_stats(TOKEN_STATS *token_stats) {
#if CONFIG_RD_DEBUG
  int r, c;
  for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r) {
    for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c) {
      token_stats->txb_coeff_cost_map[r][c] = 0;
    }
  }
#endif
  token_stats->cost = 0;
}

static INLINE void aom_start_encode(aom_writer *bc, uint8_t *buffer) {
  aom_daala_start_encode(bc, buffer);
}

static INLINE void aom_stop_encode(aom_writer *bc) {
  aom_daala_stop_encode(bc);
}

static INLINE void aom_write(aom_writer *br, int bit, int probability) {
  aom_daala_write(br, bit, probability);
}

static INLINE void aom_write_bit(aom_writer *w, int bit) {
  aom_write(w, bit, 128);  // aom_prob_half
}

static INLINE void aom_write_literal(aom_writer *w, int data, int bits) {
  int bit;

  for (bit = bits - 1; bit >= 0; bit--) aom_write_bit(w, 1 & (data >> bit));
}

static INLINE void aom_write_cdf(aom_writer *w, int symb,
                                 const aom_cdf_prob *cdf, int nsymbs) {
  daala_write_symbol(w, symb, cdf, nsymbs);
}

static INLINE void aom_write_symbol(aom_writer *w, int symb, aom_cdf_prob *cdf,
                                    int nsymbs) {
  aom_write_cdf(w, symb, cdf, nsymbs);
  if (w->allow_update_cdf) update_cdf(cdf, symb, nsymbs);
}

static INLINE void aom_write_multibit(aom_writer *w, int nbits, int data) {
  static const aom_cdf_prob multibit_cdf[4][17] = {
    { AOM_CDF2(0x4000), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { AOM_CDF4(0x2000, 0x4000, 0x6000), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { AOM_CDF8(0x1000, 0x2000, 0x3000, 0x4000, 0x5000, 0x6000, 0x7000), 0, 0, 0,
      0, 0, 0, 0, 0 },
    { AOM_CDF16(0x800, 0x1000, 0x1800, 0x2000, 0x2800, 0x3000, 0x3800, 0x4000,
                0x4800, 0x5000, 0x5800, 0x6000, 0x6800, 0x7000, 0x7800) }
  };
  assert(nbits <= 4);
  aom_write_cdf(w, data, multibit_cdf[nbits - 1], 1 << nbits);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITWRITER_H_
