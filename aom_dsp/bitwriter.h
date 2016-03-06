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

#include "aom_ports/mem.h"

#include "aom_dsp/prob.h"
#include "aom_dsp/aom_dsp_common.h"
#if CONFIG_DAALA_EC
#include "aom_dsp/entenc.h"
#endif

#define AOM_MEASURE_EC_OVERHEAD 0

#ifdef __cplusplus
extern "C" {
#endif

typedef struct aom_writer {
  unsigned int lowvalue;
  unsigned int range;
  int count;
  unsigned int pos;
  uint8_t *buffer;
#if CONFIG_DAALA_EC
  od_ec_enc ec;
#endif
#if AOM_MEASURE_EC_OVERHEAD
  double entropy;
  int nb_symbols;
#endif
} aom_writer;

void aom_start_encode(aom_writer *bc, uint8_t *buffer);
void aom_stop_encode(aom_writer *bc);

static INLINE void aom_write(aom_writer *br, int bit, int probability) {
#if CONFIG_DAALA_EC
  if (probability == 128) {
    od_ec_enc_bits(&br->ec, bit, 1);
  }
  else {
    int p = (32768*probability + (256 - probability)) >> 8;
    od_ec_encode_bool_q15(&br->ec, bit, p);
  }
#else
  unsigned int split;
  int count = br->count;
  unsigned int range = br->range;
  unsigned int lowvalue = br->lowvalue;
  register int shift;

  split = 1 + (((range - 1) * probability) >> 8);

  range = split;

  if (bit) {
    lowvalue += split;
    range = br->range - split;
  }

  shift = aom_norm[range];

  range <<= shift;
  count += shift;

  if (count >= 0) {
    int offset = shift - count;

    if ((lowvalue << (offset - 1)) & 0x80000000) {
      int x = br->pos - 1;

      while (x >= 0 && br->buffer[x] == 0xff) {
        br->buffer[x] = 0;
        x--;
      }

      br->buffer[x] += 1;
    }

    br->buffer[br->pos++] = (lowvalue >> (24 - offset));
    lowvalue <<= offset;
    shift = count;
    lowvalue &= 0xffffff;
    count -= 8;
  }

  lowvalue <<= shift;
  br->count = count;
  br->lowvalue = lowvalue;
  br->range = range;
#endif
#if AOM_MEASURE_EC_OVERHEAD
  assert(probability > 0 && probability < 256);
  br->entropy -= AOMLOG2((double)(bit ? 256 - probability : probability)/256.0);
  br->nb_symbols++;
#endif
}

static INLINE void aom_write_bit(aom_writer *w, int bit) {
  aom_write(w, bit, 128);  // aom_prob_half
}

static INLINE void aom_write_literal(aom_writer *w, int data, int bits) {
  int bit;

  for (bit = bits - 1; bit >= 0; bit--) aom_write_bit(w, 1 & (data >> bit));
}

#define aom_write_prob(w, v) aom_write_literal((w), (v), 8)

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITWRITER_H_
