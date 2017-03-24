/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom_dsp/bitwriter.h"

#include "av1/common/common.h"

void aom_write_primitive_symmetric(aom_writer *w, int word,
                                   unsigned int abs_bits) {
  if (word == 0) {
    aom_write_bit(w, 0);
  } else {
    const int x = abs(word);
    const int s = word < 0;
    aom_write_bit(w, 1);
    aom_write_bit(w, s);
    aom_write_literal(w, x - 1, abs_bits);
  }
}

// Encodes a value v in [0, n-1] quasi-uniformly
static void write_primitive_quniform(aom_writer *w, uint16_t n, uint16_t v) {
  if (n <= 1) return;
  const int l = get_msb(n - 1) + 1;
  const int m = (1 << l) - n;
  if (v < m) {
    aom_write_literal(w, v, l - 1);
  } else {
    aom_write_literal(w, m + ((v - m) >> 1), l - 1);
    aom_write_bit(w, (v - m) & 1);
  }
}

// Encodes a value v in [0, n-1] based on a reference ref also in [0, n-1]
// The closest p values of v from ref are coded using a p-ary quasi-unoform
// short code while the remaining n-p values are coded with a longer code.
void aom_write_primitive_refbilevel(aom_writer *w, uint16_t n, uint16_t p,
                                    uint16_t ref, uint16_t v) {
  if (n <= 1) return;
  assert(p > 0 && p < n);
  assert(ref < n);
  int lolimit = ref - p / 2;
  int hilimit = lolimit + p - 1;
  if (lolimit < 0) {
    lolimit = 0;
    hilimit = p - 1;
  } else if (hilimit >= n) {
    hilimit = n - 1;
    lolimit = n - p;
  }
  if (v >= lolimit && v <= hilimit) {
    aom_write_bit(w, 1);
    v = v - lolimit;
    write_primitive_quniform(w, p, v);
  } else {
    aom_write_bit(w, 0);
    if (v > hilimit) v -= p;
    write_primitive_quniform(w, n - p, v);
  }
}
