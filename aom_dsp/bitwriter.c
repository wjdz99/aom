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
#include <stdio.h>

#include "./bitwriter.h"

void aom_start_encode(aom_writer *br, uint8_t *source) {
  br->lowvalue = 0;
  br->range = 255;
  br->count = -24;
  br->buffer = source;
  br->pos = 0;
  aom_write_bit(br, 0);
#if AOM_MEASURE_EC_OVERHEAD
  br->entropy = 0;
  br->nb_symbols = 0;
#endif
}

void aom_stop_encode(aom_writer *br) {
  int i;

  for (i = 0; i < 32; i++) aom_write_bit(br, 0);

  // Ensure there's no ambigous collision with any index marker bytes
  if ((br->buffer[br->pos - 1] & 0xe0) == 0xc0) br->buffer[br->pos++] = 0;

#if AOM_MEASURE_EC_OVERHEAD
  {
    uint32_t tell;
    tell = br->pos << 3;
    fprintf(stderr, "overhead: %f%%\n", 100*(tell - br->entropy)/br->entropy);
    fprintf(stderr, "efficiency: %f bits/symbol\n",
     (double)tell/br->nb_symbols);
  }
#endif
}
