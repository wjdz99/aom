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

#ifndef AOM_DSP_ANS_H_
#define AOM_DSP_ANS_H_
// An implementation of Asymmetric Numeral Systems
// http://arxiv.org/abs/1311.2540v2

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/prob.h"
#include "aom_ports/mem_ops.h"

#define ANS_DIV(dividend, divisor) ((dividend) / (divisor))

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

struct AnsCoder {
  uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

struct AnsDecoder {
  const uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

typedef uint8_t AnsP8;
#define ans_p8_precision 256u
#define ans_p8_shift 8

#define l_base (ans_p8_precision * 4)  // l_base % precision must be 0
#define io_base 256
// Range I = { l_base, l_base + 1, ..., l_base * io_base - 1 }

static INLINE void ans_write_init(struct AnsCoder *const ans,
                                  uint8_t *const buf) {
  ans->buf = buf;
  ans->buf_offset = 0;
  ans->state = l_base;
}

static INLINE int ans_write_end(struct AnsCoder *const ans) {
  uint32_t state;
  assert(ans->state >= l_base);
  assert(ans->state < l_base * io_base);
  state = ans->state - l_base;
  if (state < (1 << 6)) {
    ans->buf[ans->buf_offset] = (0x00 << 6) + state;
    return ans->buf_offset + 1;
  } else if (state < (1 << 14)) {
    mem_put_le16(ans->buf + ans->buf_offset, (0x01 << 14) + state);
    return ans->buf_offset + 2;
  } else if (state < (1 << 22)) {
    mem_put_le24(ans->buf + ans->buf_offset, (0x02 << 22) + state);
    return ans->buf_offset + 3;
  } else {
    assert(0 && "State is too large to be serialized");
    return ans->buf_offset;
  }
}

// uABS with normalization
static INLINE void uabs_write(struct AnsCoder *ans, int val, AnsP8 p0) {
  AnsP8 p = ans_p8_precision - p0;
  const unsigned l_s = val ? p : p0;
  while (ans->state >= l_base / ans_p8_precision * io_base * l_s) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  if (!val)
    ans->state = ANS_DIV(ans->state * ans_p8_precision, p0);
  else
    ans->state = ANS_DIV((ans->state + 1) * ans_p8_precision + p - 1, p) - 1;
}

static INLINE int uabs_read(struct AnsDecoder *ans, AnsP8 p0) {
  AnsP8 p = ans_p8_precision - p0;
  int s;
  unsigned xp, sp;
  unsigned state = ans->state;
  while (state < l_base && ans->buf_offset > 0) {
    state = state * io_base + ans->buf[--ans->buf_offset];
  }
  sp = state * p;
  xp = sp / ans_p8_precision;
  s = (sp & 0xFF) >= p0;
  if (s)
    ans->state = xp;
  else
    ans->state = state - xp;
  return s;
}

static INLINE int uabs_read_bit(struct AnsDecoder *ans) {
  int s;
  unsigned state = ans->state;
  while (state < l_base && ans->buf_offset > 0) {
    state = state * io_base + ans->buf[--ans->buf_offset];
  }
  s = (int)(state & 1);
  ans->state = state >> 1;
  return s;
}

static INLINE int uabs_read_literal(struct AnsDecoder *ans, int bits) {
  int literal = 0, bit;
  assert(bits < 31);

  // TODO(aconverse): Investigate ways to read/write literals faster,
  // e.g. 8-bit chunks.
  for (bit = bits - 1; bit >= 0; bit--) literal |= uabs_read_bit(ans) << bit;

  return literal;
}

// TODO(aconverse): Replace trees with tokensets.
static INLINE int uabs_read_tree(struct AnsDecoder *ans,
                                 const aom_tree_index *tree,
                                 const AnsP8 *probs) {
  aom_tree_index i = 0;

  while ((i = tree[i + uabs_read(ans, probs[i >> 1])]) > 0) continue;

  return -i;
}

static INLINE int ans_read_init(struct AnsDecoder *const ans,
                                const uint8_t *const buf, int offset) {
  unsigned x;
  if (offset < 1) return 1;
  ans->buf = buf;
  x = buf[offset - 1] >> 6;
  if (x == 0) {
    ans->buf_offset = offset - 1;
    ans->state = buf[offset - 1] & 0x3F;
  } else if (x == 1) {
    if (offset < 2) return 1;
    ans->buf_offset = offset - 2;
    ans->state = mem_get_le16(buf + offset - 2) & 0x3FFF;
  } else if (x == 2) {
    if (offset < 3) return 1;
    ans->buf_offset = offset - 3;
    ans->state = mem_get_le24(buf + offset - 3) & 0x3FFFFF;
  } else {
    // x == 3 implies this byte is a superframe marker
    return 1;
  }
  ans->state += l_base;
  if (ans->state >= l_base * io_base) return 1;
  return 0;
}

static INLINE int ans_read_end(struct AnsDecoder *const ans) {
  return ans->state == l_base;
}

static INLINE int ans_reader_has_error(const struct AnsDecoder *const ans) {
  return ans->state < l_base && ans->buf_offset == 0;
}
#undef ANS_DIV
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // AOM_DSP_ANS_H_
