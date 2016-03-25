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

#ifdef __cplusplus
extern "C" {
#endif

void av1_tree_probs_from_distribution(vpx_tree tree,
                                       unsigned int branch_ct[/* n - 1 */][2],
                                       const unsigned int num_events[/* n */]);

struct av1_token {
  int value;
  int len;
};

void av1_tokens_from_tree(struct av1_token *, const vpx_tree_index *);

static INLINE void av1_write_tree(vpx_writer *w, const vpx_tree_index *tree,
                                   const vpx_prob *probs, int bits, int len,
                                   vpx_tree_index i) {
  do {
    const int bit = (bits >> --len) & 1;
    vpx_write(w, bit, probs[i >> 1]);
    i = tree[i + bit];
  } while (len);
}

static INLINE void av1_write_token(vpx_writer *w, const vpx_tree_index *tree,
                                    const vpx_prob *probs,
                                    const struct av1_token *token) {
  av1_write_tree(w, tree, probs, token->value, token->len, 0);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_TREEWRITER_H_
