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
#include "./aom_dsp_rtcd.h"
#include "aom/aom_image.h"
#include "aom/aom_integer.h"
#include "av1/common/quant_common.h"

// Calculate the error of a filtered and unfiltered block
void aom_clpf_detect_c(const uint8_t *rec, const uint8_t *org, int rstride,
                       int ostride, int x0, int y0, int width, int height,
                       int *sum0, int *sum1, unsigned int strength, int size) {
  int x, y;
  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x];
      int X = rec[y * rstride + x];
      int A = rec[AOMMAX(0, y - 1) * rstride + x];
      int B = rec[y * rstride + AOMMAX(0, x - 2)];
      int C = rec[y * rstride + AOMMAX(0, x - 1)];
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x];
      int delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      int Y = X + delta;
      *sum0 += (O - X) * (O - X);
      *sum1 += (O - Y) * (O - Y);
    }
  }
}

void aom_clpf_detect_multi_c(const uint8_t *rec, const uint8_t *org,
                             int rstride, int ostride, int x0, int y0,
                             int width, int height, int *sum, int size) {
  int x, y;

  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x];
      int X = rec[y * rstride + x];
      int A = rec[AOMMAX(0, y - 1) * rstride + x];
      int B = rec[y * rstride + AOMMAX(0, x - 2)];
      int C = rec[y * rstride + AOMMAX(0, x - 1)];
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)];
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)];
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x];
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

#if CONFIG_AOM_HIGHBITDEPTH
// Identical to aom_clpf_detect_c() apart from "rec" and "org".
void aom_clpf_detect_hbd_c(const uint16_t *rec, const uint16_t *org,
                           int rstride, int ostride, int x0, int y0, int width,
                           int height, int *sum0, int *sum1,
                           unsigned int strength, int shift, int size) {
  int x, y;
  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x] >> shift;
      int X = rec[y * rstride + x] >> shift;
      int A = rec[AOMMAX(0, y - 1) * rstride + x] >> shift;
      int B = rec[y * rstride + AOMMAX(0, x - 2)] >> shift;
      int C = rec[y * rstride + AOMMAX(0, x - 1)] >> shift;
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)] >> shift;
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)] >> shift;
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x] >> shift;
      int delta = av1_clpf_sample(X, A, B, C, D, E, F, strength >> shift);
      int Y = X + delta;
      *sum0 += (O - X) * (O - X);
      *sum1 += (O - Y) * (O - Y);
    }
  }
}

// aom_clpf_detect_multi_c() apart from "rec" and "org".
void aom_clpf_detect_multi_hbd_c(const uint16_t *rec, const uint16_t *org,
                                 int rstride, int ostride, int x0, int y0,
                                 int width, int height, int *sum, int shift,
                                 int size) {
  int x, y;

  for (y = y0; y < y0 + size; y++) {
    for (x = x0; x < x0 + size; x++) {
      int O = org[y * ostride + x] >> shift;
      int X = rec[y * rstride + x] >> shift;
      int A = rec[AOMMAX(0, y - 1) * rstride + x] >> shift;
      int B = rec[y * rstride + AOMMAX(0, x - 2)] >> shift;
      int C = rec[y * rstride + AOMMAX(0, x - 1)] >> shift;
      int D = rec[y * rstride + AOMMIN(width - 1, x + 1)] >> shift;
      int E = rec[y * rstride + AOMMIN(width - 1, x + 2)] >> shift;
      int F = rec[AOMMIN(height - 1, y + 1) * rstride + x] >> shift;
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
#endif

int av1_clpf_decision(int k, int l, const YV12_BUFFER_CONFIG *rec,
                      const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                      int block_size, int w, int h, unsigned int strength,
                      unsigned int fb_size_log2, uint8_t *res) {
  int m, n, sum0 = 0, sum1 = 0;

  for (m = 0; m < h; m++) {
    for (n = 0; n < w; n++) {
      int xpos = (l << fb_size_log2) + n * block_size;
      int ypos = (k << fb_size_log2) + m * block_size;
      if (!cm->mi_grid_visible[ypos / MI_SIZE * cm->mi_stride + xpos / MI_SIZE]
          ->mbmi.skip) {
#if CONFIG_AOM_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
          aom_clpf_detect_hbd(CONVERT_TO_SHORTPTR(rec->y_buffer),
                              CONVERT_TO_SHORTPTR(org->y_buffer), rec->y_stride,
                              org->y_stride, xpos, ypos, rec->y_crop_width,
                              rec->y_crop_height, &sum0, &sum1, strength,
                              cm->bit_depth - 8, block_size);
        } else {
          aom_clpf_detect(rec->y_buffer, org->y_buffer, rec->y_stride,
                          org->y_stride, xpos, ypos, rec->y_crop_width,
                          rec->y_crop_height, &sum0, &sum1, strength,
                          block_size);
        }
#else
        aom_clpf_detect(rec->y_buffer, org->y_buffer, rec->y_stride,
                        org->y_stride, xpos, ypos, rec->y_crop_width,
                        rec->y_crop_height, &sum0, &sum1, strength, block_size);
#endif
      }
    }
  }
  *res = sum1 < sum0 + (sum0 >> 8);
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
                    int h, int64_t res[4][4], int plane) {
  int c, m, n, filtered = 0;
  int sum[4];
  const int subx = plane != AOM_PLANE_Y && rec->subsampling_x;
  const int suby = plane != AOM_PLANE_Y && rec->subsampling_y;
  int bslog = get_msb(block_size);
  uint8_t *rec_buffer =
      plane != AOM_PLANE_Y
          ? (plane == AOM_PLANE_U ? rec->u_buffer : rec->v_buffer)
          : rec->y_buffer;
  uint8_t *org_buffer =
      plane != AOM_PLANE_Y
          ? (plane == AOM_PLANE_U ? org->u_buffer : org->v_buffer)
          : org->y_buffer;
  int rec_width = plane != AOM_PLANE_Y ? rec->uv_crop_width : rec->y_crop_width;
  int rec_height =
      plane != AOM_PLANE_Y ? rec->uv_crop_height : rec->y_crop_height;
  int rec_stride = plane != AOM_PLANE_Y ? rec->uv_stride : rec->y_stride;
  int org_stride = plane != AOM_PLANE_Y ? org->uv_stride : org->y_stride;
  sum[0] = sum[1] = sum[2] = sum[3] = 0;
  if (plane == AOM_PLANE_Y &&
      fb_size_log2 > (unsigned int)get_msb(MAX_FB_SIZE) - 3) {
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
                        res, plane);
    if (1 << (fb_size_log2 - bslog) < w)
      filtered |= clpf_rdo(y, x + (1 << fb_size_log2), rec, org, cm, block_size,
                           fb_size_log2, w2, h1, res, plane);
    if (1 << (fb_size_log2 - bslog) < h) {
      filtered |= clpf_rdo(y + (1 << fb_size_log2), x, rec, org, cm, block_size,
                           fb_size_log2, w1, h2, res, plane);
      filtered |=
          clpf_rdo(y + (1 << fb_size_log2), x + (1 << fb_size_log2), rec, org,
                   cm, block_size, fb_size_log2, w2, h2, res, plane);
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
      if (!cm->mi_grid_visible[(ypos << suby) / MI_SIZE * cm->mi_stride +
                               (xpos << subx) / MI_SIZE]
          ->mbmi.skip) {
#if CONFIG_AOM_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
          aom_clpf_detect_multi_hbd(
              CONVERT_TO_SHORTPTR(rec_buffer), CONVERT_TO_SHORTPTR(org_buffer),
              rec_stride, org_stride, xpos, ypos, rec_width, rec_height, sum,
              cm->bit_depth - 8, block_size);
        } else {
          aom_clpf_detect_multi(rec_buffer, org_buffer, rec_stride, org_stride,
                                xpos, ypos, rec_width, rec_height, sum,
                                block_size);
        }
#else
        aom_clpf_detect_multi(rec_buffer, org_buffer, rec_stride, org_stride,
                              xpos, ypos, rec_width, rec_height, sum,
                              block_size);
#endif
        filtered = 1;
      }
    }
  }

  for (c = 0; c < (plane == AOM_PLANE_Y ? 4 : 1); c++) {
    res[c][0] += sum[0];
    res[c][1] += sum[1];
    res[c][2] += sum[2];
    res[c][3] += sum[3];
  }
  return filtered;
}

void av1_clpf_test_frame(const YV12_BUFFER_CONFIG *rec,
                         const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                         int *best_strength, int *best_bs, int plane) {
  int c, j, k, l;
  int64_t best, sums[4][4];
  int width = plane != AOM_PLANE_Y ? rec->uv_crop_width : rec->y_crop_width;
  int height = plane != AOM_PLANE_Y ? rec->uv_crop_height : rec->y_crop_height;
  const int bs = MI_SIZE;
  const int bslog = get_msb(bs);
  int fb_size_log2 = get_msb(MAX_FB_SIZE);
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;

  memset(sums, 0, sizeof(sums));

  if (plane != AOM_PLANE_Y)
    // Use a block size of MI_SIZE regardless of the subsampling.  This
    // This is accurate enough to determine the best strength and
    // we don't need to add SIMD optimisations for 4x4 blocks.
    clpf_rdo(0, 0, rec, org, cm, bs, fb_size_log2, width >> bslog,
             height >> bslog, sums, plane);
  else
    for (k = 0; k < num_fb_ver; k++) {
      for (l = 0; l < num_fb_hor; l++) {
        // Calculate the block size after frame border clipping
        int h =
            AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
        int w =
            AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
        h += !h << fb_size_log2;
        w += !w << fb_size_log2;
        clpf_rdo(k << fb_size_log2, l << fb_size_log2, rec, org, cm, MI_SIZE,
                 fb_size_log2, w >> bslog, h >> bslog, sums, plane);
      }
    }

  if (plane != AOM_PLANE_Y)  // Slightly favour unfiltered chroma
    sums[0][0] -= sums[0][0] >> 7;

  for (j = 0; j < 4; j++) {
    static const double lambda_square[] = {
      // exp(x + 3 / 34)
      1.0922, 1.1248, 1.1584, 1.1930, 1.2286, 1.2653, 1.3030, 1.3419, 1.3820,
      1.4232, 1.4657, 1.5095, 1.5545, 1.6009, 1.6487, 1.6979, 1.7486, 1.8008,
      1.8546, 1.9099, 1.9669, 2.0256, 2.0861, 2.1484, 2.2125, 2.2785, 2.3465,
      2.4166, 2.4887, 2.5630, 2.6395, 2.7183, 2.7994, 2.8830, 2.9690, 3.0577,
      3.1489, 3.2429, 3.3397, 3.4394, 3.5421, 3.6478, 3.7567, 3.8688, 3.9843,
      4.1032, 4.2257, 4.3518, 4.4817, 4.6155, 4.7532, 4.8951, 5.0412, 5.1917,
      5.3467, 5.5062, 5.6706, 5.8399, 6.0142, 6.1937, 6.3786, 6.5689, 6.7650,
      6.9669, 7.1749, 7.3891, 7.6096, 7.8367, 8.0707, 8.3116, 8.5596, 8.8151,
      9.0783, 9.3492, 9.6283, 9.9157, 10.212, 10.516, 10.830, 11.154, 11.487,
      11.829, 12.182, 12.546, 12.921, 13.306, 13.703, 14.112, 14.534, 14.968,
      15.414, 15.874, 16.348, 16.836, 17.339, 17.856, 18.389, 18.938, 19.503,
      20.086, 20.685, 21.302, 21.938, 22.593, 23.268, 23.962, 24.677, 25.414,
      26.172, 26.954, 27.758, 28.587, 29.440, 30.319, 31.224, 32.156, 33.115,
      34.104, 35.122, 36.170, 37.250, 38.362, 39.507, 40.686, 41.900, 43.151,
      44.439, 45.765, 47.131, 48.538, 49.987, 51.479, 53.016, 54.598, 56.228,
      57.906, 59.635, 61.415, 63.248, 65.136, 67.080, 69.082, 71.144, 73.268,
      75.454, 77.707, 80.026, 82.415, 84.875, 87.408, 90.017, 92.704, 95.471,
      98.321, 101.26, 104.28, 107.39, 110.60, 113.90, 117.30, 120.80, 124.40,
      128.12, 131.94, 135.88, 139.93, 144.11, 148.41, 152.84, 157.41, 162.10,
      166.94, 171.93, 177.06, 182.34, 187.78, 193.39, 199.16, 205.11, 211.23,
      217.53, 224.03, 230.71, 237.60, 244.69, 252.00, 259.52, 267.26, 275.24,
      283.46, 291.92, 300.63, 309.60, 318.85, 328.36, 338.16, 348.26, 358.65,
      369.36, 380.38, 391.74, 403.43, 415.47, 427.87, 440.64, 453.8,  467.34,
      481.29, 495.66, 510.45, 525.69, 541.38, 557.54, 574.18, 591.32, 608.97,
      627.14, 645.86, 665.14, 685.00, 705.44, 726.5,  748.18, 770.51, 793.51,
      817.2,  841.59, 866.71, 892.58, 919.22, 946.66, 974.92, 1004.0, 1034.0,
      1064.8, 1096.6, 1129.4, 1163.1, 1197.8, 1233.5, 1270.4, 1308.3, 1347.3,
      1387.5, 1429.0, 1471.6, 1515.5, 1560.8, 1607.4, 1655.3, 1704.8, 1755.6,
      1808.0, 1862.0, 1917.6, 1974.8
    };

    // Estimate the bit costs and adjust the square errors
    double lambda =
        lambda_square[av1_get_qindex(&cm->seg, 0, cm->base_qindex)];
    int i, cost = (int)(lambda * (sums[j][0] + 6 + 2 * (j > 0)));
    for (i = 0; i < 4; i++)
      sums[j][i] = ((sums[j][i] + (i && j) * cost) << 4) + j * 4 + i;
  }

  best = (int64_t)1 << 62;
  for (c = 0; c < (plane == AOM_PLANE_Y ? 4 : 1); c++)
    for (j = 0; j < 4; j++)
      if ((!c || j) && sums[c][j] < best) best = sums[c][j];
  best &= 15;
  if (best_bs)
    *best_bs = (best > 3) * (5 + (best < 12) + (best < 8));
  *best_strength = best ? 1 << ((best - 1) & 3) : 0;
}
