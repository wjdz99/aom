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
#ifndef AV1_COMMON_CLPF_H_
#define AV1_COMMON_CLPF_H_

#include "av1/common/reconinter.h"

#define clip(n, l, h) \
  ((h) < ((n) > (l) ? (n) : (l)) ? (h) : ((n) > (l) ? (n) : (l)))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#if __GNUC__
static unsigned int log2i(uint32_t x) { return 31 - __builtin_clz(x); }
#else
static unsigned int log2i(uint32_t n) {
  int c = 0;
  assert(n > 0);
  while (n >>= 1) c++;
  return c;
}
#endif

#define MAX_FB_SIZE 128

int clpf_maxbits(const AV1_COMMON *cm);
int clpf_sample(int X, int A, int B, int C, int D, int E, int F, int b);
void clpf_block(const uint8_t *src, uint8_t *dst, int stride, int x0, int y0,
                int sizex, int sizey, int width, int height,
                unsigned int strength);
int av1_clpf_frame(const YV12_BUFFER_CONFIG *dst, const YV12_BUFFER_CONFIG *rec,
                   const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                   int enable_fb_flag, unsigned int strength,
                   unsigned int fb_size_log2, char *blocks,
                   int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                                   const YV12_BUFFER_CONFIG *,
                                   const AV1_COMMON *cm, int, int, int,
                                   unsigned int, unsigned int, char *));

#endif
