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
#include <stdint.h>
#include <string.h>

#include "av1/encoder/palette.h"
#include "av1/encoder/random.h"

#define AV1_K_MEANS_RENAME(func, dim) func##_dim##dim##_c
#define AV1_K_MEANS_RENAME2(func, dim) func##_dim##dim

#define VAR_DIM(AV1_K_MEANS_DIM)                                               \
  static int AV1_K_MEANS_RENAME(calc_dist, AV1_K_MEANS_DIM)(const int *p1,     \
                                                            const int *p2) {   \
    int dist = 0;                                                              \
    for (int i = 0; i < AV1_K_MEANS_DIM; ++i) {                                \
      const int diff = p1[i] - p2[i];                                          \
      dist += diff * diff;                                                     \
    }                                                                          \
    return dist;                                                               \
  }                                                                            \
                                                                               \
  void AV1_K_MEANS_RENAME(av1_calc_indices, AV1_K_MEANS_DIM)(                  \
      const int *data, const int *centroids, uint8_t *indices, int n, int k) { \
    for (int i = 0; i < n; ++i) {                                              \
      int min_dist = AV1_K_MEANS_RENAME(calc_dist, AV1_K_MEANS_DIM)(           \
          data + i * AV1_K_MEANS_DIM, centroids);                              \
      indices[i] = 0;                                                          \
      for (int j = 1; j < k; ++j) {                                            \
        const int this_dist = AV1_K_MEANS_RENAME(calc_dist, AV1_K_MEANS_DIM)(  \
            data + i * AV1_K_MEANS_DIM, centroids + j * AV1_K_MEANS_DIM);      \
        if (this_dist < min_dist) {                                            \
          min_dist = this_dist;                                                \
          indices[i] = j;                                                      \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }                                                                            \
                                                                               \
  static void AV1_K_MEANS_RENAME(calc_centroids, AV1_K_MEANS_DIM)(             \
      const int *data, int *centroids, const uint8_t *indices, int n, int k) { \
    int i, j;                                                                  \
    int count[PALETTE_MAX_SIZE] = { 0 };                                       \
    unsigned int rand_state = (unsigned int)data[0];                           \
    assert(n <= 32768);                                                        \
    memset(centroids, 0, sizeof(centroids[0]) * k * AV1_K_MEANS_DIM);          \
                                                                               \
    for (i = 0; i < n; ++i) {                                                  \
      const int index = indices[i];                                            \
      assert(index < k);                                                       \
      ++count[index];                                                          \
      for (j = 0; j < AV1_K_MEANS_DIM; ++j) {                                  \
        centroids[index * AV1_K_MEANS_DIM + j] +=                              \
            data[i * AV1_K_MEANS_DIM + j];                                     \
      }                                                                        \
    }                                                                          \
                                                                               \
    for (i = 0; i < k; ++i) {                                                  \
      if (count[i] == 0) {                                                     \
        memcpy(centroids + i * AV1_K_MEANS_DIM,                                \
               data + (lcg_rand16(&rand_state) % n) * AV1_K_MEANS_DIM,         \
               sizeof(centroids[0]) * AV1_K_MEANS_DIM);                        \
      } else {                                                                 \
        for (j = 0; j < AV1_K_MEANS_DIM; ++j) {                                \
          centroids[i * AV1_K_MEANS_DIM + j] =                                 \
              DIVIDE_AND_ROUND(centroids[i * AV1_K_MEANS_DIM + j], count[i]);  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }                                                                            \
                                                                               \
  static int64_t AV1_K_MEANS_RENAME(calc_total_dist, AV1_K_MEANS_DIM)(         \
      const int *data, const int *centroids, const uint8_t *indices, int n,    \
      int k) {                                                                 \
    int64_t dist = 0;                                                          \
    (void)k;                                                                   \
    for (int i = 0; i < n; ++i) {                                              \
      dist += AV1_K_MEANS_RENAME(calc_dist, AV1_K_MEANS_DIM)(                  \
          data + i * AV1_K_MEANS_DIM,                                          \
          centroids + indices[i] * AV1_K_MEANS_DIM);                           \
    }                                                                          \
    return dist;                                                               \
  }                                                                            \
                                                                               \
  void AV1_K_MEANS_RENAME(av1_k_means, AV1_K_MEANS_DIM)(                       \
      const int *data, int *centroids, uint8_t *indices, int n, int k,         \
      int max_itr) {                                                           \
    int pre_centroids[2 * PALETTE_MAX_SIZE];                                   \
    uint8_t pre_indices[MAX_SB_SQUARE];                                        \
                                                                               \
    AV1_K_MEANS_RENAME(av1_calc_indices, AV1_K_MEANS_DIM)                      \
    (data, centroids, indices, n, k);                                          \
    int64_t this_dist = AV1_K_MEANS_RENAME(calc_total_dist, AV1_K_MEANS_DIM)(  \
        data, centroids, indices, n, k);                                       \
                                                                               \
    for (int i = 0; i < max_itr; ++i) {                                        \
      const int64_t pre_dist = this_dist;                                      \
      memcpy(pre_centroids, centroids,                                         \
             sizeof(pre_centroids[0]) * k * AV1_K_MEANS_DIM);                  \
      memcpy(pre_indices, indices, sizeof(pre_indices[0]) * n);                \
                                                                               \
      AV1_K_MEANS_RENAME(calc_centroids, AV1_K_MEANS_DIM)                      \
      (data, centroids, indices, n, k);                                        \
      AV1_K_MEANS_RENAME(av1_calc_indices, AV1_K_MEANS_DIM)                    \
      (data, centroids, indices, n, k);                                        \
      this_dist = AV1_K_MEANS_RENAME(calc_total_dist, AV1_K_MEANS_DIM)(        \
          data, centroids, indices, n, k);                                     \
                                                                               \
      if (this_dist > pre_dist) {                                              \
        memcpy(centroids, pre_centroids,                                       \
               sizeof(pre_centroids[0]) * k * AV1_K_MEANS_DIM);                \
        memcpy(indices, pre_indices, sizeof(pre_indices[0]) * n);              \
        break;                                                                 \
      }                                                                        \
      if (!memcmp(centroids, pre_centroids,                                    \
                  sizeof(pre_centroids[0]) * k * AV1_K_MEANS_DIM))             \
        break;                                                                 \
    }                                                                          \
  }

VAR_DIM(1);
VAR_DIM(2);
#undef VAR_DIM
