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

#include <math.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "./config.h"
#endif

#include "./aom_dsp_rtcd.h"
#include "./av1_rtcd.h"
#include "./cdef.h"

/* Generated from gen_filter_tables.c. */
const int cdef_directions[8][3] = {
    {-1 * OD_FILT_BSTRIDE + 1, -2 * OD_FILT_BSTRIDE + 2,
     -3 * OD_FILT_BSTRIDE + 3},
    {0 * OD_FILT_BSTRIDE + 1, -1 * OD_FILT_BSTRIDE + 2,
     -1 * OD_FILT_BSTRIDE + 3},
    {0 * OD_FILT_BSTRIDE + 1, 0 * OD_FILT_BSTRIDE + 2, 0 * OD_FILT_BSTRIDE + 3},
    {0 * OD_FILT_BSTRIDE + 1, 1 * OD_FILT_BSTRIDE + 2, 1 * OD_FILT_BSTRIDE + 3},
    {1 * OD_FILT_BSTRIDE + 1, 2 * OD_FILT_BSTRIDE + 2, 3 * OD_FILT_BSTRIDE + 3},
    {1 * OD_FILT_BSTRIDE + 0, 2 * OD_FILT_BSTRIDE + 1, 3 * OD_FILT_BSTRIDE + 1},
    {1 * OD_FILT_BSTRIDE + 0, 2 * OD_FILT_BSTRIDE + 0, 3 * OD_FILT_BSTRIDE + 0},
    {1 * OD_FILT_BSTRIDE + 0, 2 * OD_FILT_BSTRIDE - 1, 3 * OD_FILT_BSTRIDE - 1},
};

/* Detect direction. 0 means 45-degree up-right, 2 is horizontal, and so on.
   The search minimizes the weighted variance along all the lines in a
   particular direction, i.e. the squared error between the input and a
   "predicted" block where each pixel is replaced by the average along a line
   in a particular direction. Since each direction have the same sum(x^2) term,
   that term is never computed. See Section 2, step 2, of:
   http://jmvalin.ca/notes/intra_paint.pdf */
int cdef_find_dir_c(const uint16_t *img, int stride, int32_t *var,
                    int coeff_shift) {
  int i;
  int32_t cost[8] = { 0 };
  int partial[8][15] = { { 0 } };
  int32_t best_cost = 0;
  int best_dir = 0;
  /* Instead of dividing by n between 2 and 8, we multiply by 3*5*7*8/n.
     The output is then 840 times larger, but we don't care for finding
     the max. */
  static const int div_table[] = { 0, 840, 420, 280, 210, 168, 140, 120, 105 };
  for (i = 0; i < 8; i++) {
    int j;
    for (j = 0; j < 8; j++) {
      int x;
      /* We subtract 128 here to reduce the maximum range of the squared
         partial sums. */
      x = (img[i * stride + j] >> coeff_shift) - 128;
      partial[0][i + j] += x;
      partial[1][i + j / 2] += x;
      partial[2][i] += x;
      partial[3][3 + i - j / 2] += x;
      partial[4][7 + i - j] += x;
      partial[5][3 - i / 2 + j] += x;
      partial[6][j] += x;
      partial[7][i / 2 + j] += x;
    }
  }
  for (i = 0; i < 8; i++) {
    cost[2] += partial[2][i] * partial[2][i];
    cost[6] += partial[6][i] * partial[6][i];
  }
  cost[2] *= div_table[8];
  cost[6] *= div_table[8];
  for (i = 0; i < 7; i++) {
    cost[0] += (partial[0][i] * partial[0][i] +
                partial[0][14 - i] * partial[0][14 - i]) *
               div_table[i + 1];
    cost[4] += (partial[4][i] * partial[4][i] +
                partial[4][14 - i] * partial[4][14 - i]) *
               div_table[i + 1];
  }
  cost[0] += partial[0][7] * partial[0][7] * div_table[8];
  cost[4] += partial[4][7] * partial[4][7] * div_table[8];
  for (i = 1; i < 8; i += 2) {
    int j;
    for (j = 0; j < 4 + 1; j++) {
      cost[i] += partial[i][3 + j] * partial[i][3 + j];
    }
    cost[i] *= div_table[8];
    for (j = 0; j < 4 - 1; j++) {
      cost[i] += (partial[i][j] * partial[i][j] +
                  partial[i][10 - j] * partial[i][10 - j]) *
                 div_table[2 * j + 2];
    }
  }
  for (i = 0; i < 8; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
  /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
  *var = best_cost - cost[(best_dir + 4) & 7];
  /* We'd normally divide by 840, but dividing by 1024 is close enough
     for what we're going to do with this. */
  *var >>= 10;
  return best_dir;
}

/* Smooth in the direction detected. */
void cdef_filter_block_c(uint16_t *y, int ystride, const uint16_t *in,
                         int pri_strength, int sec_strength, int dir,
                         int pri_damping, int sec_damping, int bsize, int max) {
  int i, j, k;
  static const int pri_taps[3] = { 4, 3 };
  static const int sec_taps[3] = { 2, 1 };
  const int s = OD_FILT_BSTRIDE;
  for (i = 0; i < 4 << (bsize == BLOCK_8X8); i++) {
    for (j = 0; j < 4 << (bsize == BLOCK_8X8); j++) {
      int16_t sum = 0;
      uint16_t yy;
      uint16_t xx = in[i * s + j];
      for (k = 0; k < 2; k++) {
        uint16_t p0 = in[i * s + j + cdef_directions[dir][k]];
        uint16_t p1 = in[i * s + j - cdef_directions[dir][k]];
        uint16_t s0 = in[i * s + j + cdef_directions[(dir + 2) & 7][k]];
        uint16_t s1 = in[i * s + j - cdef_directions[(dir + 2) & 7][k]];
        uint16_t s2 = in[i * s + j + cdef_directions[(dir + 6) & 7][k]];
        uint16_t s3 = in[i * s + j - cdef_directions[(dir + 6) & 7][k]];
        sum += pri_taps[k] * constrain(p0 - xx, pri_strength, pri_damping);
        sum += pri_taps[k] * constrain(p1 - xx, pri_strength, pri_damping);
        sum += sec_taps[k] * constrain(s0 - xx, sec_strength, sec_damping);
        sum += sec_taps[k] * constrain(s1 - xx, sec_strength, sec_damping);
        sum += sec_taps[k] * constrain(s2 - xx, sec_strength, sec_damping);
        sum += sec_taps[k] * constrain(s3 - xx, sec_strength, sec_damping);
      }
      yy = clamp(xx + ((8 + sum - (sum < 0)) >> 4), 0, max);
      y[i * ystride + j] = yy;
    }
  }
}

/* Compute the primary filter strength for an 8x8 block based on the
   directional variance difference. A high variance difference means
   that we have a highly directional pattern (e.g. a high contrast
   edge), so we can apply more deringing. A low variance means that we
   either have a low contrast edge, or a non-directional texture, so
   we want to be careful not to blur. */
static INLINE int od_adjust_strength(int strength, int32_t var) {
  const int i = var >> 6 ? AOMMIN(get_msb(var >> 6), 12) : 0;
  /* We use the variance of 8x8 blocks to adjust the strength. */
  return var ? (strength * (4 + i) + 8) >> 4 : 0;
}

void copy_8x8_16bit_to_16bit_c(uint16_t *dst, int dstride, const uint16_t *src,
                               int sstride) {
  int i, j;
  for (i = 0; i < 8; i++)
    for (j = 0; j < 8; j++) dst[i * dstride + j] = src[i * sstride + j];
}

void copy_4x4_16bit_to_16bit_c(uint16_t *dst, int dstride, const uint16_t *src,
                               int sstride) {
  int i, j;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) dst[i * dstride + j] = src[i * sstride + j];
}

void copy_dering_16bit_to_16bit(uint16_t *dst, int dstride, uint16_t *src,
                                dering_list *dlist, int dering_count,
                                int bsize) {
  int bi, bx, by;

  if (bsize == BLOCK_8X8) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_8x8_16bit_to_16bit(&dst[(by << 3) * dstride + (bx << 3)], dstride,
                              &src[bi << (3 + 3)], 8);
    }
  } else if (bsize == BLOCK_4X8) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_16bit(&dst[(by << 3) * dstride + (bx << 2)], dstride,
                              &src[bi << (3 + 2)], 4);
      copy_4x4_16bit_to_16bit(&dst[((by << 3) + 4) * dstride + (bx << 2)],
                              dstride, &src[(bi << (3 + 2)) + 4 * 4], 4);
    }
  } else if (bsize == BLOCK_8X4) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_16bit(&dst[(by << 2) * dstride + (bx << 3)], dstride,
                              &src[bi << (2 + 3)], 8);
      copy_4x4_16bit_to_16bit(&dst[(by << 2) * dstride + (bx << 3) + 4],
                              dstride, &src[(bi << (2 + 3)) + 4], 8);
    }
  } else {
    assert(bsize == BLOCK_4X4);
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_16bit(&dst[(by << 2) * dstride + (bx << 2)], dstride,
                              &src[bi << (2 + 2)], 4);
    }
  }
}

void copy_8x8_16bit_to_8bit_c(uint8_t *dst, int dstride, const uint16_t *src,
                              int sstride) {
  int i, j;
  for (i = 0; i < 8; i++)
    for (j = 0; j < 8; j++)
      dst[i * dstride + j] = (uint8_t)src[i * sstride + j];
}

void copy_4x4_16bit_to_8bit_c(uint8_t *dst, int dstride, const uint16_t *src,
                              int sstride) {
  int i, j;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      dst[i * dstride + j] = (uint8_t)src[i * sstride + j];
}

static void copy_dering_16bit_to_8bit(uint8_t *dst, int dstride,
                                      const uint16_t *src, dering_list *dlist,
                                      int dering_count, int bsize) {
  int bi, bx, by;
  if (bsize == BLOCK_8X8) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_8x8_16bit_to_8bit(&dst[(by << 3) * dstride + (bx << 3)], dstride,
                             &src[bi << (3 + 3)], 8);
    }
  } else if (bsize == BLOCK_4X8) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_8bit(&dst[(by << 3) * dstride + (bx << 2)], dstride,
                             &src[bi << (3 + 2)], 4);
      copy_4x4_16bit_to_8bit(&dst[((by << 3) + 4) * dstride + (bx << 2)],
                             dstride, &src[(bi << (3 + 2)) + 4 * 4], 4);
    }
  } else if (bsize == BLOCK_8X4) {
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_8bit(&dst[(by << 2) * dstride + (bx << 3)], dstride,
                             &src[bi << (2 + 3)], 8);
      copy_4x4_16bit_to_8bit(&dst[(by << 2) * dstride + (bx << 3) + 4], dstride,
                             &src[(bi << (2 + 3)) + 4], 8);
    }
  } else {
    assert(bsize == BLOCK_4X4);
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      copy_4x4_16bit_to_8bit(&dst[(by << 2) * dstride + (bx << 2)], dstride,
                             &src[bi << (2 * 2)], 4);
    }
  }
}

int get_filter_skip(int level) {
  int filter_skip = level & 1;
  if (level == 1) filter_skip = 0;
  return filter_skip;
}

void od_dering(uint8_t *dst, int dstride, uint16_t *y, uint16_t *in, int xdec,
               int ydec, int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
               int *dirinit, int var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
               int pli, dering_list *dlist, int dering_count, int level,
               int sec_strength, int pri_damping, int sec_damping,
               int coeff_shift, int hbd) {
  int bi;
  int bx;
  int by;
  int bsize, bsizex, bsizey;

  int pri_strength = (level >> 1) << coeff_shift;
  int filter_skip = get_filter_skip(level);
  if (level == 1)
    pri_strength = 24 << coeff_shift;

  sec_damping += coeff_shift - (pli != AOM_PLANE_Y);
  pri_damping += coeff_shift - (pli != AOM_PLANE_Y);
  bsize =
      ydec ? (xdec ? BLOCK_4X4 : BLOCK_8X4) : (xdec ? BLOCK_4X8 : BLOCK_8X8);
  bsizex = 3 - xdec;
  bsizey = 3 - ydec;

  if (dirinit && pri_strength == 0 && sec_strength == 0) {
    // If we're here, both primary and secondary strengths are 0, and
    // we still haven't written anything to y[] yet, so we just copy
    // the input to y[]. This is necessary only for av1_cdef_search()
    // and only av1_cdef_search() sets dirinit.
    for (bi = 0; bi < dering_count; bi++) {
      by = dlist[bi].by;
      bx = dlist[bi].bx;
      int iy, ix;
      // TODO(stemidts/jmvalin): SIMD optimisations
      for (iy = 0; iy < 1 << bsizey; iy++)
        for (ix = 0; ix < 1 << bsizex; ix++)
          y[(bi << (bsizex + bsizey)) + (iy << bsizex) + ix] =
              in[((by << bsizey) + iy) * OD_FILT_BSTRIDE + (bx << bsizex) + ix];
    }
    return;
  }

  if (pli == 0) {
    if (!dirinit || !*dirinit) {
      for (bi = 0; bi < dering_count; bi++) {
        by = dlist[bi].by;
        bx = dlist[bi].bx;
        dir[by][bx] = cdef_find_dir(&in[8 * by * OD_FILT_BSTRIDE + 8 * bx],
                                    OD_FILT_BSTRIDE, &var[by][bx], coeff_shift);
      }
      if (dirinit) *dirinit = 1;
    }
  }

  assert(bsize == BLOCK_8X8 || bsize == BLOCK_4X4);
  for (bi = 0; bi < dering_count; bi++) {
    int t = !filter_skip && dlist[bi].skip ? 0 : pri_strength;
    int s = !filter_skip && dlist[bi].skip ? 0 : sec_strength;
    by = dlist[bi].by;
    bx = dlist[bi].bx;
    cdef_filter_block_c(&y[bi << (bsizex + bsizey)], 1 << bsizex,
                        &in[(by * OD_FILT_BSTRIDE << bsizey) + (bx << bsizex)],
                        (pli ? t : od_adjust_strength(t, var[by][bx])), s,
                        t ? dir[by][bx] : 0, pri_damping, sec_damping, bsize,
                        (256 << coeff_shift) - 1);
  }
  if (dst) {
    if (hbd) {
      copy_dering_16bit_to_16bit((uint16_t *)dst, dstride, y, dlist,
                                 dering_count, bsize);
    } else {
      copy_dering_16bit_to_8bit(dst, dstride, y, dlist, dering_count, bsize);
    }
  }
}
