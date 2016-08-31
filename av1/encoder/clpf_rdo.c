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

#include "av1/common/clpf.h"
#include "aom/aom_integer.h"
#include "av1/common/quant_common.h"

// Calculate the error of a filtered and unfiltered block
static void detect_clpf(const uint8_t *rec, const uint8_t *org, int x0, int y0,
                        int width, int height, int so, int stride, int *sum0,
                        int *sum1, unsigned int strength, int size) {
  int x, y;
  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * so + x];
      int X = rec[y * stride + x];
      int A = rec[AOMMAX(0, y - 1) * stride + x];
      int B = rec[y * stride + AOMMAX(0, x - 2)];
      int C = rec[y * stride + AOMMAX(0, x - 1)];
      int D = rec[y * stride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * stride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * stride + x];
      int delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      int Y = X + delta;
      *sum0 += (O - X) * (O - X);
      *sum1 += (O - Y) * (O - Y);
    }
  }
}

static void detect_multi_clpf(const uint8_t *rec, const uint8_t *org, int x0,
                              int y0, int width, int height, int so, int stride,
                              int *sum, int size) {
  int x, y;

  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * so + x];
      int X = rec[y * stride + x];
      int A = rec[AOMMAX(0, y - 1) * stride + x];
      int B = rec[y * stride + AOMMAX(0, x - 2)];
      int C = rec[y * stride + AOMMAX(0, x - 1)];
      int D = rec[y * stride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * stride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * stride + x];
      int delta1 = av1_clpf_sample(X, A, B, C, D, E, F, 1);
      int delta2 = av1_clpf_sample(X, A, B, C, D, E, F, 2);
      int delta3 = av1_clpf_sample(X, A, B, C, D, E, F, 4);
      int F1 = X + delta1;
      int F2 = X + delta2;
      int F3 = X + delta3;
      sum[0] += (O - X) * (O - X);
      sum[1] += (O - F1) * (O - F1);
      sum[2] += (O - F2) * (O - F2);
      sum[3] += (O - F3) * (O - F3);
    }
  }
}

int av1_clpf_decision(int k, int l, const YV12_BUFFER_CONFIG *rec,
                      const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                      int block_size, int w, int h, unsigned int strength,
                      unsigned int fb_size_log2, uint8_t *res, int comp) {
  int m, n, sum0 = 0, sum1 = 0;
  int size = 8 >> (comp && (rec->subsampling_x || rec->subsampling_x));
  uint8_t *rec_buffer =
      comp ? (comp == 1 ? rec->u_buffer : rec->v_buffer) : rec->y_buffer;
  uint8_t *org_buffer =
      comp ? (comp == 1 ? org->u_buffer : org->v_buffer) : org->y_buffer;
  int rec_width = comp ? rec->uv_crop_width : rec->y_crop_width;
  int rec_height = comp ? rec->uv_crop_height : rec->y_crop_height;
  int rec_stride = comp ? rec->uv_stride : rec->y_stride;
  int org_stride = comp ? org->uv_stride : org->y_stride;

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = (l << fb_size_log2) + n * block_size;
      int ypos = (k << fb_size_log2) + m * block_size;
      if (!cm->mi_grid_visible[ypos / MAX_MIB_SIZE * cm->mi_stride +
                               xpos / MAX_MIB_SIZE]
               ->mbmi.skip)
        detect_clpf(rec_buffer, org_buffer, xpos, ypos, rec_width, rec_height,
                    org_stride, rec_stride, &sum0, &sum1, strength, size);
    }
  }
  *res = sum1 < sum0;
  return *res;
}

// Calculate the square error of all filter settings.  Result:
// res[0][0]   : unfiltered
// res[0][1-3] : strength=1,2,4, no signals
// (Only for luma:)
// res[1][0]   : (bit count, fb size = 128)
// res[1][1-3] : strength=1,2,4, fb size = 128
// res[2][0]   : (bit count, fb size = 64)
// res[2][1-3] : strength=1,2,4, fb size = 64
// res[3][0]   : (bit count, fb size = 32)
// res[3][1-3] : strength=1,2,4, fb size = 32
static int clpf_rdo(int y, int x, const YV12_BUFFER_CONFIG *rec,
                    const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                    unsigned int block_size, unsigned int fb_size_log2, int w,
                    int h, int64_t res[4][4], int comp) {
  int c, m, n, filtered = 0;
  int size = 8 >> (comp && (rec->subsampling_x || rec->subsampling_x));
  int sum[4];
  int bslog = get_msb(block_size);
  uint8_t *rec_buffer =
      comp ? (comp == 1 ? rec->u_buffer : rec->v_buffer) : rec->y_buffer;
  uint8_t *org_buffer =
      comp ? (comp == 1 ? org->u_buffer : org->v_buffer) : org->y_buffer;
  int rec_width = comp ? rec->uv_crop_width : rec->y_crop_width;
  int rec_height = comp ? rec->uv_crop_height : rec->y_crop_height;
  int rec_stride = comp ? rec->uv_stride : rec->y_stride;
  int org_stride = comp ? org->uv_stride : org->y_stride;
  sum[0] = sum[1] = sum[2] = sum[3] = 0;
  if (!comp && fb_size_log2 > (unsigned int)get_msb(MAX_FB_SIZE) - 3) {
    int w1, h1, w2, h2, i, sum1, sum2, sum3, oldfiltered;

    fb_size_log2--;
    w1 = AOMMIN(1 << (fb_size_log2 - bslog), w);
    h1 = AOMMIN(1 << (fb_size_log2 - bslog), h);
    w2 = AOMMIN(w - (1 << (fb_size_log2 - bslog)), w >> 1);
    h2 = AOMMIN(h - (1 << (fb_size_log2 - bslog)), h >> 1);
    i = get_msb(MAX_FB_SIZE) - fb_size_log2;
    sum1 = res[i][1];
    sum2 = res[i][2];
    sum3 = res[i][3];
    oldfiltered = res[i][0];
    res[i][0] = 0;

    filtered = clpf_rdo(y, x, rec, org, cm, block_size, fb_size_log2, w1, h1,
                        res, comp);
    if (1 << (fb_size_log2 - bslog) < w)
      filtered |= clpf_rdo(y, x + (1 << fb_size_log2), rec, org, cm, block_size,
                           fb_size_log2, w2, h1, res, comp);
    if (1 << (fb_size_log2 - bslog) < h) {
      filtered |= clpf_rdo(y + (1 << fb_size_log2), x, rec, org, cm, block_size,
                           fb_size_log2, w1, h2, res, comp);
      filtered |=
          clpf_rdo(y + (1 << fb_size_log2), x + (1 << fb_size_log2), rec, org,
                   cm, block_size, fb_size_log2, w2, h2, res, comp);
    }

    res[i][1] = AOMMIN(sum1 + res[i][0], res[i][1]);
    res[i][2] = AOMMIN(sum2 + res[i][0], res[i][2]);
    res[i][3] = AOMMIN(sum3 + res[i][0], res[i][3]);
    res[i][0] = oldfiltered + filtered;  // Number of signal bits
    return filtered;
  }

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = x + n * block_size;
      int ypos = y + m * block_size;
      if (!cm->mi_grid_visible[ypos / MAX_MIB_SIZE * cm->mi_stride +
                               xpos / MAX_MIB_SIZE]
               ->mbmi.skip) {
        detect_multi_clpf(rec_buffer, org_buffer, xpos, ypos, rec_width,
                          rec_height, org_stride, rec_stride, sum, size);
        filtered = 1;
      }
    }
  }

  for (c = 0; c < (comp ? 1 : 4); c++) {
    res[c][0] += sum[0];
    res[c][1] += sum[1];
    res[c][2] += sum[2];
    res[c][3] += sum[3];
  }
  return filtered;
}

void av1_clpf_test_frame(const YV12_BUFFER_CONFIG *rec,
                         const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                         int *best_strength, int *best_bs, int comp) {
  int i, j, k, l;
  int64_t best, sums[4][4];
  int width = comp ? rec->uv_crop_width : rec->y_crop_width;
  int height = comp ? rec->uv_crop_height : rec->y_crop_height;
  const int bs = 8 >> (comp && (rec->subsampling_x || rec->subsampling_x));
  const int bslog = get_msb(bs);
  int fb_size_log2 = get_msb(MAX_FB_SIZE);
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;

  memset(sums, 0, sizeof(sums));

  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      // Calculate the block size after frame border clipping
      int h =
          AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      int w =
          AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      h += !h << fb_size_log2;
      w += !w << fb_size_log2;
      clpf_rdo(k << fb_size_log2, l << fb_size_log2, rec, org, cm, bs,
               fb_size_log2, w >> bslog, h >> bslog, sums, comp);
    }
  }

  for (j = 0; j < 4; j++) {
    static double lambda_square[] = {
      // exp((x / 8.8))
      1.0000, 1.1203, 1.2552, 1.4062, 1.5755, 1.7651, 1.9775, 2.2154, 2.4821,
      2.7808, 3.1154, 3.4903, 3.9104, 4.3810, 4.9082, 5.4989, 6.1606, 6.9020,
      7.7327, 8.6633, 9.7058, 10.8739, 12.1825, 13.6486, 15.2911, 17.1313,
      19.1930, 21.5028, 24.0905, 26.9897, 30.2377, 33.8767, 37.9536, 42.5211,
      47.6383, 53.3713, 59.7942, 66.9901, 75.0520, 84.0842, 94.2032, 105.5401,
      118.2413, 132.4710, 148.4132, 166.2739, 186.2841, 208.7024, 233.8186,
      261.9574, 293.4825, 328.8016, 368.3711, 412.7026, 462.3692, 518.0128,
      580.3529, 650.1953, 728.4428, 816.1070, 914.3211, 1024.3548, 1147.6304,
      1285.7415
    };

    // Estimate the bit costs and adjust the square errors
    double lambda =
        lambda_square[av1_get_qindex(&cm->seg, 0, cm->base_qindex) >> 2];
    int c, cost = (int)((lambda * (sums[j][0] + 2 + 2 * (j > 0)) + 0.5));
    for (c = 0; c < 4; c++)
      sums[j][c] = ((sums[j][c] + (j && c) * cost) << 4) + j * 4 + c;
  }

  best = (int64_t)1 << 62;
  for (i = 0; i < (comp ? 1 : 4); i++)
    for (j = 0; j < 4; j++)
      if ((!i || j) && sums[i][j] < best) best = sums[i][j];
  best &= 15;
  *best_bs = (best > 3) * (5 + (best < 12) + (best < 8));
  *best_strength = best ? 1 << ((best - 1) & 3) : 0;
}
