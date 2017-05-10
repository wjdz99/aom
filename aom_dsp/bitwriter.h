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
#if CONFIG_EC_ADAPT && !CONFIG_EC_MULTISYMBOL
#error "CONFIG_EC_ADAPT is enabled without enabling CONFIG_EC_MULTISYMBOL"
#endif

#if CONFIG_ANS
#include "aom_dsp/buf_ans.h"
#elif CONFIG_DAALA_EC
#include "aom_dsp/daalaboolwriter.h"
#else
#include "aom_dsp/dkboolwriter.h"
#endif
#include "aom_dsp/prob.h"

#if CONFIG_RD_DEBUG
#include "av1/common/blockd.h"
#include "av1/encoder/cost.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_ANS
typedef struct BufAnsCoder AomWriter;
#elif CONFIG_DAALA_EC
typedef struct DaalaWriter AomWriter;
#else
typedef struct AomDkWriter AomWriter;
#endif

typedef struct TokenStats {
  int cost;
#if CONFIG_VAR_TX
#if CONFIG_RD_DEBUG
  int txb_coeff_cost_map[TXB_COEFF_COST_MAP_SIZE][TXB_COEFF_COST_MAP_SIZE];
#endif
#endif
} TokenStats;

static INLINE void init_token_stats(TokenStats *token_stats) {
#if CONFIG_VAR_TX
#if CONFIG_RD_DEBUG
  int r, c;
  for (r = 0; r < TXB_COEFF_COST_MAP_SIZE; ++r) {
    for (c = 0; c < TXB_COEFF_COST_MAP_SIZE; ++c) {
      token_stats->txb_coeff_cost_map[r][c] = 0;
    }
  }
#endif
#endif
  token_stats->cost = 0;
}

static INLINE void aom_start_encode(AomWriter *bc, uint8_t *buffer) {
#if CONFIG_ANS
  (void)bc;
  (void)buffer;
  assert(0 && "buf_ans requires a more complicated startup procedure");
#elif CONFIG_DAALA_EC
  aom_daala_start_encode(bc, buffer);
#else
  aom_dk_start_encode(bc, buffer);
#endif
}

static INLINE void aom_stop_encode(AomWriter *bc) {
#if CONFIG_ANS
  (void)bc;
  assert(0 && "buf_ans requires a more complicated shutdown procedure");
#elif CONFIG_DAALA_EC
  aom_daala_stop_encode(bc);
#else
  aom_dk_stop_encode(bc);
#endif
}

static INLINE void aom_write(AomWriter *br, int bit, int probability) {
#if CONFIG_ANS
  buf_rabs_write(br, bit, probability);
#elif CONFIG_DAALA_EC
  aom_daala_write(br, bit, probability);
#else
  aom_dk_write(br, bit, probability);
#endif
}

static INLINE void aom_write_record(AomWriter *br, int bit, int probability,
                                    TokenStats *token_stats) {
  aom_write(br, bit, probability);
#if CONFIG_RD_DEBUG
  token_stats->cost += av1_cost_bit(probability, bit);
#else
  (void)token_stats;
#endif
}

static INLINE void aom_write_bit(AomWriter *w, int bit) {
#if CONFIG_ANS
  buf_rabs_write_bit(w, bit);
#elif CONFIG_DAALA_EC && CONFIG_RAWBITS
  // Note this uses raw bits and is not the same as aom_daala_write(r, 128);
  aom_daala_write_bit(w, bit);
#else
  aom_write(w, bit, 128);  // aom_prob_half
#endif
}

static INLINE void aom_write_bit_record(AomWriter *w, int bit,
                                        TokenStats *token_stats) {
  aom_write_bit(w, bit);
#if CONFIG_RD_DEBUG
  token_stats->cost += av1_cost_bit(128, bit);  // aom_prob_half
#else
  (void)token_stats;
#endif
}

static INLINE void aom_write_literal(AomWriter *w, int data, int bits) {
  int bit;

  for (bit = bits - 1; bit >= 0; bit--) aom_write_bit(w, 1 & (data >> bit));
}

static INLINE void aom_write_tree_as_bits(AomWriter *w, const AomTreeIndex *tr,
                                          const AomProb *probs, int bits,
                                          int len, AomTreeIndex i) {
  do {
    const int bit = (bits >> --len) & 1;
    aom_write(w, bit, probs[i >> 1]);
    i = tr[i + bit];
  } while (len);
}

static INLINE void aom_write_tree_as_bits_record(AomWriter *w,
                                                 const AomTreeIndex *tr,
                                                 const AomProb *probs, int bits,
                                                 int len, AomTreeIndex i,
                                                 TokenStats *token_stats) {
  do {
    const int bit = (bits >> --len) & 1;
    aom_write_record(w, bit, probs[i >> 1], token_stats);
    i = tr[i + bit];
  } while (len);
}

#if CONFIG_EC_MULTISYMBOL
static INLINE void aom_write_cdf(AomWriter *w, int symb, const AomCdfProb *cdf,
                                 int nsymbs) {
#if CONFIG_ANS
  (void)nsymbs;
  assert(cdf);
  const AomCdfProb cum_prob = symb > 0 ? cdf[symb - 1] : 0;
  const AomCdfProb prob = cdf[symb] - cum_prob;
  buf_rans_write(w, cum_prob, prob);
#elif CONFIG_DAALA_EC
  daala_write_symbol(w, symb, cdf, nsymbs);
#else
#error \
    "CONFIG_EC_MULTISYMBOL is selected without a valid backing entropy " \
  "coder. Enable daala_ec or ans for a valid configuration."
#endif
}

static INLINE void aom_write_symbol(AomWriter *w, int symb, AomCdfProb *cdf,
                                    int nsymbs) {
  aom_write_cdf(w, symb, cdf, nsymbs);
#if CONFIG_EC_ADAPT
  update_cdf(cdf, symb, nsymbs);
#endif
}

static INLINE void aom_write_tree_as_cdf(AomWriter *w, const AomTreeIndex *tree,
                                         const AomProb *probs, int bits,
                                         int len, AomTreeIndex i) {
  AomTreeIndex root;
  root = i;
  do {
    AomCdfProb cdf[16];
    AomTreeIndex index[16];
    int path[16];
    int dist[16];
    int nsymbs;
    int symb;
    int j;
    /* Compute the CDF of the binary tree using the given probabilities. */
    nsymbs = tree_to_cdf(tree, probs, root, cdf, index, path, dist);
    /* Find the symbol to code. */
    symb = -1;
    for (j = 0; j < nsymbs; j++) {
      /* If this symbol codes a leaf node,  */
      if (index[j] <= 0) {
        if (len == dist[j] && path[j] == bits) {
          symb = j;
          break;
        }
      } else {
        if (len > dist[j] && path[j] == bits >> (len - dist[j])) {
          symb = j;
          break;
        }
      }
    }
    OD_ASSERT(symb != -1);
    aom_write_cdf(w, symb, cdf, nsymbs);
    bits &= (1 << (len - dist[symb])) - 1;
    len -= dist[symb];
  } while (len);
}

#endif  // CONFIG_EC_MULTISYMBOL

static INLINE void aom_write_tree(AomWriter *w, const AomTreeIndex *tree,
                                  const AomProb *probs, int bits, int len,
                                  AomTreeIndex i) {
#if CONFIG_EC_MULTISYMBOL
  aom_write_tree_as_cdf(w, tree, probs, bits, len, i);
#else
  aom_write_tree_as_bits(w, tree, probs, bits, len, i);
#endif
}

static INLINE void aom_write_tree_record(AomWriter *w, const AomTreeIndex *tree,
                                         const AomProb *probs, int bits,
                                         int len, AomTreeIndex i,
                                         TokenStats *token_stats) {
#if CONFIG_EC_MULTISYMBOL
  (void)token_stats;
  aom_write_tree_as_cdf(w, tree, probs, bits, len, i);
#else
  aom_write_tree_as_bits_record(w, tree, probs, bits, len, i, token_stats);
#endif
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITWRITER_H_
