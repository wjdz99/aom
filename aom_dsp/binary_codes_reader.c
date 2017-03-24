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

#include "aom_dsp/bitreader.h"

#include "av1/common/common.h"

int aom_read_primitive_symmetric(aom_reader *r, unsigned int mag_bits) {
  if (aom_read_bit(r, ACCT_STR_NAME)) {
    int s = aom_read_bit(r, ACCT_STR_NAME);
    int x = aom_read_literal(r, mag_bits, ACCT_STR_NAME) + 1;
    return (s > 0 ? -x : x);
  } else {
    return 0;
  }
}

static uint16_t read_primitive_quniform(aom_reader *r, uint16_t n) {
  if (n <= 1) return 0;
  const int l = get_msb(n - 1) + 1;
  const int m = (1 << l) - n;
  const int v = aom_read_literal(r, l - 1, NULL);
  return v < m ? v : (v << 1) - m + aom_read_bit(r, NULL);
}

uint16_t aom_read_primitive_refbilevel(aom_reader *r, uint16_t n, uint16_t p,
                                       uint16_t ref) {
  if (n <= 1) return 0;
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
  int v;
  if (aom_read_bit(r, NULL)) {
    v = read_primitive_quniform(r, p) + lolimit;
  } else {
    v = read_primitive_quniform(r, n - p);
    if (v >= lolimit) v += p;
  }
  return v;
}
