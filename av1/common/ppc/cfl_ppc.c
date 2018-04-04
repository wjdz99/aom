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

#define OFF_0 0
#define OFF_1 16
#define OFF_2 32
#define OFF_3 48
#define CFL_BUF_LINE_BYTES 64

typedef vector int16_t int16x8_t;
typedef vector uint16_t uint16x8_t;
typedef vector int32_t int32x4_t;
typedef vector uint32_t uint32x4_t;
typedef vector uint64_t uint64x2_t;

static int16x8_t vec_ld_64(const int16_t *a) {
  uint64x2_t vec = { 0, 0 };
  vec[0] = *((int64_t *)a);
  return (int16x8_t)vec;
}

static void vec_st_64(const int16x8_t a, int16_t *b) {
  *((int64_t *)b) = ((uint64x2_t)a)[0];
}

static INLINE void subtract_average_vsx(int16_t *pred_buf, int width,
                                        int height, int round_offset,
                                        int num_pel_log2) {
  const int16_t *end = pred_buf + height * CFL_BUF_LINE;
  const int16_t *sum_buf = pred_buf;
  const int32x4_t zeros = vec_splats(0);
  const int32x4_t offset = vec_splats(round_offset);
  const uint32x4_t div_shift = vec_splats((uint32_t)num_pel_log2);

  int32x4_t sum_32x4 = { 0, 0, 0, 0 };
  do {
    int16x8_t v0;
    int16x8_t v1;
    if (width == 4) {
      v0 = vec_ld_64(sum_buf);
      v1 = vec_ld_64(sum_buf + CFL_BUF_LINE);
    } else {
      v0 = vec_vsx_ld(OFF_0, sum_buf);
      v1 = vec_vsx_ld(OFF_0 + CFL_BUF_LINE_BYTES, sum_buf);
      if (width >= 16) {
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_1, sum_buf), sum_32x4);
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_1 + CFL_BUF_LINE_BYTES, sum_buf),
                             sum_32x4);
      }
      if (width == 32) {
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_2, sum_buf), sum_32x4);
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_2 + CFL_BUF_LINE_BYTES, sum_buf),
                             sum_32x4);
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_3, sum_buf), sum_32x4);
        sum_32x4 = vec_sum4s(vec_vsx_ld(OFF_3 + CFL_BUF_LINE_BYTES, sum_buf),
                             sum_32x4);
      }
    }
    sum_32x4 = vec_sum4s(v0, sum_32x4);
    sum_32x4 = vec_sum4s(v1, sum_32x4);
  } while ((sum_buf += (CFL_BUF_LINE * 2)) < end);

  const int32x4_t sum = vec_splat(vec_sums(sum_32x4, zeros), 3);
  const int32x4_t avg = vec_sr(vec_add(sum, offset), div_shift);
  const int16x8_t vec_avg = vec_pack(avg, avg);
  do {
    if (width == 4) {
      vec_st_64(vec_sub(vec_vsx_ld(OFF_0, pred_buf), vec_avg), pred_buf);
    } else {
      vec_vsx_st(vec_sub(vec_vsx_ld(OFF_0, pred_buf), vec_avg), 0, pred_buf);
      if (width >= 16) {
        vec_vsx_st(vec_sub(vec_vsx_ld(OFF_1, pred_buf), vec_avg), OFF_1,
                   pred_buf);
      }
      if (width == 32) {
        vec_vsx_st(vec_sub(vec_vsx_ld(OFF_2, pred_buf), vec_avg), OFF_2,
                   pred_buf);
        vec_vsx_st(vec_sub(vec_vsx_ld(OFF_3, pred_buf), vec_avg), OFF_3,
                   pred_buf);
      }
    }
  } while ((pred_buf += CFL_BUF_LINE) < end);
}

CFL_SUB_AVG_FN(vsx)
