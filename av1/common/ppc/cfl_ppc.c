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

#include <altivec.h>

#include "./av1_rtcd.h"

#include "av1/common/cfl.h"

typedef vector signed short int16x8_t;
typedef vector unsigned short uint16x8_t;
typedef vector signed int int32x4_t;
typedef vector unsigned int uint32x4_t;
typedef vector unsigned long uint64x2_t;

static int16x8_t vec_ld_64(const int16_t *a) {
  vector unsigned long vec = { 0, 0 };
  vec[0] = *((long *)a);
  return (int16x8_t)vec;
}

static void vec_st_64(const int16x8_t a, int16_t *b) {
  *((long *)b) = ((uint64x2_t)a)[0];
}

static INLINE void subtract_average_vsx(int16_t *pred_buf, int width,
                                        int height, int round_offset,
                                        int num_pel_log2) {
  const int16_t *end = pred_buf + height * CFL_BUF_LINE;
  const int16_t *sum_buf = pred_buf;

  int32x4_t sum_32x4 = { 0, 0, 0, 0 };
  int32x4_t sum = { 0, 0, 0, 0 };
  int16x8_t row;
  do {
    if (width == 4) {
      row = vec_ld_64(sum_buf);
    } else {
      row = vec_vsx_ld(0, sum_buf);
      if (width >= 16) {
        int16x8_t row2 = vec_vsx_ld(0, sum_buf + 8);
        sum_32x4 = vec_sum4s(row2, sum_32x4);
      }
      if (width == 32) {
        int16x8_t row3 = vec_vsx_ld(0, sum_buf + 16);
        int16x8_t row4 = vec_vsx_ld(0, sum_buf + 24);
        sum_32x4 = vec_sum4s(row3, sum_32x4);
        sum_32x4 = vec_sum4s(row4, sum_32x4);
      }
      //    printf("%d %d %d %d %d %d %d %d\n", row[0], row[1], row[2], row[3],
      //           row[4], row[5], row[6], row[7]);
    }
    sum_32x4 = vec_sum4s(row, sum_32x4);
    //   printf("%d %d %d %d\n", sum_32x4[0], sum_32x4[1], sum_32x4[2],
    //   sum_32x4[3]);
  } while ((sum_buf += CFL_BUF_LINE) < end);
  sum = vec_sums(sum_32x4, sum);
  //  printf("%d %d %d %d\n", sum_32x4[0], sum_32x4[1], sum_32x4[2],
  //  sum_32x4[3]);

  const int16_t avg = (sum[3] + round_offset) >> num_pel_log2;
  const int16x8_t vec_avg = { avg, avg, avg, avg, avg, avg, avg, avg };
  do {
    if (width == 4) {
      vec_st_64(vec_sub(vec_vsx_ld(0, pred_buf), vec_avg), pred_buf);
    } else {
      vec_vsx_st(vec_sub(vec_vsx_ld(0, pred_buf), vec_avg), 0, pred_buf);
      if (width >= 16) {
        vec_vsx_st(vec_sub(vec_vsx_ld(0, pred_buf + 8), vec_avg), 0,
                   pred_buf + 8);
      }
      if (width == 32) {
        vec_vsx_st(vec_sub(vec_vsx_ld(0, pred_buf + 16), vec_avg), 0,
                   pred_buf + 16);
        vec_vsx_st(vec_sub(vec_vsx_ld(0, pred_buf + 24), vec_avg), 0,
                   pred_buf + 24);
      }
    }
  } while ((pred_buf += CFL_BUF_LINE) < end);
}

CFL_SUB_AVG_FN(vsx)
