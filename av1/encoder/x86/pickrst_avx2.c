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

#include <immintrin.h>  // avx2

#include "./av1_rtcd.h"
#include "./aom_dsp_rtcd.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "aom_ports/system_state.h"

static double find_average_highbd(const uint16_t *src, int h_start, int h_end,
                                  int v_start, int v_end, int stride) {
  uint64_t sum = 0;
  double avg = 0;
  int i, j;
  aom_clear_system_state();
  for (i = v_start; i < v_end; i++)
    for (j = h_start; j < h_end; j++) sum += src[i * stride + j];
  avg = (double)sum / ((v_end - v_start) * (h_end - h_start));
  return avg;
}

void compute_stats_highbd_avx2(int wiener_win, const uint8_t *dgd8,
                               const uint8_t *src8, int h_start, int h_end,
                               int v_start, int v_end, int dgd_stride,
                               int src_stride, double *M, double *H) {
  const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
  const uint16_t *dgd = CONVERT_TO_SHORTPTR(dgd8);
  const double avg =
      find_average_highbd(dgd, h_start, h_end, v_start, v_end, dgd_stride);
  const int wiener_win2 = wiener_win * wiener_win;
  const int wiener_halfwin = (wiener_win >> 1);
  double Y[WIENER_WIN2];
  assert(wiener_win == 7 || wiener_win == 5);
  memset(M, 0, sizeof(*M) * wiener_win2);
  memset(H, 0, sizeof(*H) * wiener_win2 * wiener_win2);
  for (int i = v_start; i < v_end; i++) {
    for (int j = h_start; j < h_end; j++) {
      const uint16_t *dgd2 =
          &dgd[(i - wiener_halfwin) * dgd_stride + (j - wiener_halfwin)];
      if (wiener_win == 7) {
#define KLOOP(K)                                 \
  Y[0 * 7 + K] = dgd2[K * dgd_stride + 0] - avg; \
  Y[1 * 7 + K] = dgd2[K * dgd_stride + 1] - avg; \
  Y[2 * 7 + K] = dgd2[K * dgd_stride + 2] - avg; \
  Y[3 * 7 + K] = dgd2[K * dgd_stride + 3] - avg; \
  Y[4 * 7 + K] = dgd2[K * dgd_stride + 4] - avg; \
  Y[5 * 7 + K] = dgd2[K * dgd_stride + 5] - avg; \
  Y[6 * 7 + K] = dgd2[K * dgd_stride + 6] - avg
        KLOOP(0);
        KLOOP(1);
        KLOOP(2);
        KLOOP(3);
        KLOOP(4);
        KLOOP(5);
        KLOOP(6);
#undef KLOOP
      } else if (wiener_win == 5) {
#define KLOOP(K)                                 \
  Y[0 * 5 + K] = dgd2[K * dgd_stride + 0] - avg; \
  Y[1 * 5 + K] = dgd2[K * dgd_stride + 1] - avg; \
  Y[2 * 5 + K] = dgd2[K * dgd_stride + 2] - avg; \
  Y[3 * 5 + K] = dgd2[K * dgd_stride + 3] - avg; \
  Y[4 * 5 + K] = dgd2[K * dgd_stride + 4] - avg;
        KLOOP(0);
        KLOOP(1);
        KLOOP(2);
        KLOOP(3);
        KLOOP(4);
#undef KLOOP
      } else {
        assert(0);
      }
      const double X = (double)src[i * src_stride + j] - avg;
      const __m128d vX = _mm_set_sd(X);
      for (int k = 0; k < wiener_win2; ++k) {
        // H is a symmetric matrix, so we only need to fill out the upper
        // triangle here. We can copy it down to the lower triangle outside
        // the (i, j) loops.
        const double *Y2 = &Y[k];
        const __m128d xYk = _mm_load_sd(Y2);
        const int s = wiener_win2 - k;
        const size_t n16 = s / 16;
        const __m256d vYk = _mm256_broadcastsd_pd(xYk);
        _mm_store_sd(&M[k], _mm_fmadd_sd(xYk, vX, _mm_load_sd(&M[k])));
        double *H2 = &H[k * wiener_win2 + k];
        if (n16) {
          __m256d vh0 = _mm256_loadu_pd(&H2[0]);
          __m256d vh1 = _mm256_loadu_pd(&H2[4]);
          __m256d vh2 = _mm256_loadu_pd(&H2[8]);
          __m256d vh3 = _mm256_loadu_pd(&H2[12]);
          for (size_t m = 0; m < n16 - 1; ++m, H2 += 16, Y2 += 16) {
            __m256d nvh0 = _mm256_loadu_pd(&H2[16]);
            __m256d nvh1 = _mm256_loadu_pd(&H2[20]);
            __m256d nvh2 = _mm256_loadu_pd(&H2[24]);
            __m256d nvh3 = _mm256_loadu_pd(&H2[28]);
            _mm256_storeu_pd(
                &H2[0], _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[0]), vh0));
            _mm256_storeu_pd(
                &H2[4], _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[4]), vh1));
            _mm256_storeu_pd(
                &H2[8], _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[8]), vh2));
            _mm256_storeu_pd(
                &H2[12], _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[12]), vh3));
            vh0 = nvh0;
            vh1 = nvh1;
            vh2 = nvh2;
            vh3 = nvh3;
          }
          _mm256_storeu_pd(&H2[0],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[0]), vh0));
          _mm256_storeu_pd(&H2[4],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[4]), vh1));
          _mm256_storeu_pd(&H2[8],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[8]), vh2));
          _mm256_storeu_pd(&H2[12],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[12]), vh3));
          H2 += 16;
          Y2 += 16;
        }
        if (s & 8) {
          __m256d vh0 = _mm256_loadu_pd(&H2[0]);
          __m256d vh1 = _mm256_loadu_pd(&H2[4]);
          _mm256_storeu_pd(&H2[0],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[0]), vh0));
          _mm256_storeu_pd(&H2[4],
                           _mm256_fmadd_pd(vYk, _mm256_loadu_pd(&Y2[4]), vh1));
          H2 += 8;
          Y2 += 8;
        }
        if (s & 4) {
          __m256d vh0 = _mm256_loadu_pd(H2);
          _mm256_storeu_pd(H2, _mm256_fmadd_pd(vYk, _mm256_loadu_pd(Y2), vh0));
          H2 += 4;
          Y2 += 4;
        }
        if (s & 2) {
          __m128d vh0 = _mm_loadu_pd(H2);
          _mm_storeu_pd(H2, _mm_fmadd_pd(_mm256_castpd256_pd128(vYk),
                                         _mm_loadu_pd(Y2), vh0));
          H2 += 2;
          Y2 += 2;
        }
        if (s & 1) {
          __m128d vh0 = _mm_load_sd(H2);
          _mm_store_sd(H2, _mm_fmadd_sd(_mm256_castpd256_pd128(vYk),
                                        _mm_load_sd(Y2), vh0));
        }
      }
    }
  }
  for (int k = 0; k < wiener_win2; ++k) {
    for (int l = k + 1; l < wiener_win2; ++l) {
      H[l * wiener_win2 + k] = H[k * wiener_win2 + l];
    }
  }
}
