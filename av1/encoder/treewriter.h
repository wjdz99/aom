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

#ifndef AV1_ENCODER_TREEWRITER_H_
#define AV1_ENCODER_TREEWRITER_H_

#include "aom_dsp/bitwriter.h"
#if CONFIG_DAALA_EC
#include "aom_dsp/prob.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void av1_tree_probs_from_distribution(aom_tree tree,
                                      unsigned int branch_ct[/* n - 1 */][2],
                                      const unsigned int num_events[/* n */]);

struct av1_token {
  int value;
  int len;
};

void av1_tokens_from_tree(struct av1_token *, const aom_tree_index *);

#if CONFIG_DAALA_EC
static INLINE void daala_write_tree(aom_writer *w, const aom_tree_index *tree,
                                    const aom_prob *probs, int bits, int len,
                                    aom_tree_index i) {
  aom_tree_index root;
  root = i;
  do {
    uint16_t cdf[16];
    aom_tree_index index[16];
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
      }
      else {
        if (len > dist[j] && path[j] == bits >> (len - dist[j])) {
          symb = j;
          break;
        }
      }
    }
    OD_ASSERT(symb != -1);
    od_ec_encode_cdf_q15(&w->ec, symb, cdf, nsymbs);
    bits &= (1 << (len - dist[symb])) - 1;
    len -= dist[symb];
  }
  while (len);
}
#endif

static INLINE void av1_write_tree(aom_writer *w, const aom_tree_index *tree,
                                  const aom_prob *probs, int bits, int len,
                                  aom_tree_index i) {
  do {
    const int bit = (bits >> --len) & 1;
    aom_write(w, bit, probs[i >> 1]);
    i = tree[i + bit];
  } while (len);
}

static INLINE void av1_write_token(aom_writer *w, const aom_tree_index *tree,
                                   const aom_prob *probs,
                                   const struct av1_token *token) {
#if CONFIG_DAALA_EC
  daala_write_tree(w, tree, probs, token->value, token->len, 0);
#else
  av1_write_tree(w, tree, probs, token->value, token->len, 0);
#endif
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_TREEWRITER_H_
