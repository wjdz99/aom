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

#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "noise_util.h"

#include "aom_dsp/fwd_txfm.h"
#include "aom_dsp/inv_txfm.h"

#include "aom_dsp/noise_util.h"
#include "aom_mem/aom_mem.h"

#if CONFIG_DENOISE
#define NOISE_USE_FFTW 1

#ifdef NOISE_USE_FFTW
#include <fftw3.h>
#endif
#endif

#define AOM_NOISE_TX_SCALE 255.0

// Return normally distrbuted values with standard deviation of sigma.
double aom_randn(double sigma) {
  while (1) {
    const double u = 2.0 * ((double)rand()) / RAND_MAX - 1.0;
    const double v = 2.0 * ((double)rand()) / RAND_MAX - 1.0;
    const double s = u * u + v * v;
    if (s > 0 && s < 1) {
      return sigma * (u * sqrt(-2.0 * log(s) / s));
    }
  }
  return 0;
}

#if CONFIG_DENOISE
double aom_noise_psd_get_default_value(int block_size, double factor) {
#ifdef NOISE_USE_FFTW
  fprintf(stderr, "Level: %f\n",
          (factor * factor / 10000) * 0.25 * block_size * block_size);
  return (factor * factor / 10000) * 0.25 * block_size * block_size /
         2;  /// AOM_NOISE_TX_SCALE;
  // (double)block_size / AOM_NOISE_TX_SCALE;
#else
  return (1 << 12) * 2.5 * (double)block_size / 255.;
#endif
}

#ifdef NOISE_USE_FFTW
#warning NOISE_USE_FFTW
struct aom_noise_tx_t {
  double *block;
  fftw_complex *tx_block;
  fftw_plan for_plan;
  fftw_plan inv_plan;
  int block_size;
};

struct aom_noise_tx_t *aom_noise_tx_malloc(int block_size) {
  struct aom_noise_tx_t *noise_tx =
      (struct aom_noise_tx_t *)aom_malloc(sizeof(struct aom_noise_tx_t));
  noise_tx->block_size = block_size;
  noise_tx->block =
      (double *)aom_malloc(sizeof(double) * block_size * block_size);
  noise_tx->tx_block = (fftw_complex *)aom_malloc(sizeof(fftw_complex) *
                                                  block_size * block_size);
  noise_tx->for_plan = fftw_plan_dft_r2c_2d(
      block_size, block_size, noise_tx->block, noise_tx->tx_block, 0);
  noise_tx->inv_plan = fftw_plan_dft_c2r_2d(
      block_size, block_size, noise_tx->tx_block, noise_tx->block, 0);
  if (!noise_tx->block || !noise_tx->tx_block || !noise_tx->for_plan ||
      !noise_tx->inv_plan) {
    aom_noise_tx_free(noise_tx);
    aom_free(noise_tx);
    return 0;
  }
  return noise_tx;
}

void aom_noise_tx_forward(struct aom_noise_tx_t *noise_tx, double *data) {
  const int block_size = noise_tx->block_size;
  double *block = noise_tx->block;
  memcpy(block, data, sizeof(double) * block_size * block_size);
  memset(noise_tx->tx_block, 0, sizeof(fftw_complex) * block_size * block_size);
  fftw_execute(noise_tx->for_plan);
}

void aom_noise_tx_filter(struct aom_noise_tx_t *noise_tx, double *psd) {
  const int block_size = noise_tx->block_size;
  const double kBeta = 1.1;
  const double kEps = 1e-8;
  int i = 0;
  for (i = 0; i < block_size * block_size; ++i) {
    fftw_complex *c = noise_tx->tx_block + i;
    const double p = (*c)[0] * (*c)[0] + (*c)[1] * (*c)[1];
    if (p > kBeta * psd[i] && p > 1e-6) {
      noise_tx->tx_block[i][0] *= (p - psd[i]) / AOMMAX(p, kEps);
      noise_tx->tx_block[i][1] *= (p - psd[i]) / AOMMAX(p, kEps);
    } else {
      noise_tx->tx_block[i][0] *= (kBeta - 1.0) / kBeta;
      noise_tx->tx_block[i][1] *= (kBeta - 1.0) / kBeta;
    }
  }
}

void aom_noise_tx_inverse(struct aom_noise_tx_t *noise_tx, double *data) {
  const int n = noise_tx->block_size * noise_tx->block_size;
  int i = 0;
  memset(noise_tx->block, 0, sizeof(double) * n);
  fftw_execute(noise_tx->inv_plan);
  // fftw returns values scaled by n; undo this scaling.
  for (i = 0; i < n; ++i) {
    data[i] = noise_tx->block[i] / n;
  }
}

void aom_noise_tx_add_energy(const struct aom_noise_tx_t *noise_tx,
                             double *psd) {
  const int block_size = noise_tx->block_size;
  for (int yb = 0; yb < block_size; ++yb) {
    for (int xb = 0; xb < block_size; ++xb) {
      fftw_complex *c = noise_tx->tx_block + yb * block_size + xb;
      psd[yb * block_size + xb] += (*c)[0] * (*c)[0] + (*c)[1] * (*c)[1];
    }
  }
}

void aom_noise_tx_free(struct aom_noise_tx_t *noise_tx) {
  if (!noise_tx) return;
  fftw_destroy_plan(noise_tx->for_plan);
  fftw_destroy_plan(noise_tx->inv_plan);
  aom_free(noise_tx->tx_block);
  aom_free(noise_tx->block);
  memset(noise_tx, 0, sizeof(struct aom_noise_tx_t));
}

#else

struct aom_noise_tx_t {
  int16_t *block;
  tran_low_t *tx_block;
  int block_size;
};

// Do an idct but keep the results in a 16-bit integer.
void aom_idct32x32_uint16(const tran_low_t *input, tran_low_t *dest,
                          int stride) {
  tran_low_t out[32 * 32];
  tran_low_t *outptr = out;
  int i, j;
  tran_low_t temp_in[32], temp_out[32];

  // Rows
  for (i = 0; i < 32; ++i) {
    int16_t zero_coeff[16];
    for (j = 0; j < 16; ++j) zero_coeff[j] = input[2 * j] | input[2 * j + 1];
    for (j = 0; j < 8; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 4; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];
    for (j = 0; j < 2; ++j)
      zero_coeff[j] = zero_coeff[2 * j] | zero_coeff[2 * j + 1];

    if (zero_coeff[0] | zero_coeff[1])
      aom_idct32_c(input, outptr);
    else
      memset(outptr, 0, sizeof(tran_low_t) * 32);
    input += 32;
    outptr += 32;
  }

  // Columns
  for (i = 0; i < 32; ++i) {
    for (j = 0; j < 32; ++j) temp_in[j] = out[j * 32 + i];
    aom_idct32_c(temp_in, temp_out);
    for (j = 0; j < 32; ++j) {
      dest[j * stride + i] = temp_out[j];
    }
  }
}

struct aom_noise_tx_t *aom_noise_tx_malloc(int block_size) {
  struct aom_noise_tx_t *noise_tx =
      (struct aom_noise_tx_t *)aom_malloc(sizeof(struct aom_noise_tx_t));
  noise_tx->block_size = block_size;
  noise_tx->block =
      (int16_t *)aom_malloc(sizeof(int16_t) * block_size * block_size);
  noise_tx->tx_block =
      (tran_low_t *)aom_malloc(sizeof(tran_low_t) * block_size * block_size);
  if (!noise_tx->block || !noise_tx->tx_block) {
    aom_noise_tx_free(noise_tx);
    aom_free(noise_tx);
    return 0;
  }
  return noise_tx;
}

void aom_noise_tx_forward(struct aom_noise_tx_t *noise_tx, double *data) {
  const int block_size = noise_tx->block_size;
  for (int row = 0; row < block_size; ++row) {
    for (int col = 0; col < block_size; ++col) {
      noise_tx->block[row * block_size + col] =
          (int16_t)(255 * 4 * data[row * block_size + col]);
    }
  }
  memset(noise_tx->tx_block, 0, sizeof(tran_low_t) * block_size * block_size);
  aom_fdct32x32_c(noise_tx->block, noise_tx->tx_block, block_size);
}

void aom_noise_tx_filter(struct aom_noise_tx_t *noise_tx, double *psd) {
  const int block_size = noise_tx->block_size;
  const double kBeta = 1.1;
  const double kEps = 1e-8;
  for (int i = 0; i < block_size * block_size; ++i) {
    tran_low_t *c = noise_tx->tx_block + i;
    const double p = c[0] * c[0];
    if (psd[i] == 0) continue;

    if (p > kBeta * psd[i] * 2 && p > 1e-6) {
      noise_tx->tx_block[i] *= (p - psd[i] * 2) / AOMMAX(p, kEps);
    } else {
      noise_tx->tx_block[i] *= (kBeta - 1.0) / kBeta;
    }
  }
}

void aom_noise_tx_inverse(struct aom_noise_tx_t *noise_tx, double *data) {
  const int block_size = noise_tx->block_size;
  // for (int row = 0; row < block_size; ++row) {
  // for (int col = 0; col < block_size; ++col) {
  //    noise_tx->buffer2[row * block_size + col] = 128;
  //  }
  //}
  // memset(noise_tx->buffer2, 0, sizeof(uint8_t) * block_size * block_size);
  // aom_idct32x32_1024_add_c(noise_tx->tx_block, noise_tx->buffer2,
  // block_size);

  tran_low_t out[32 * 32];
  aom_idct32x32_uint16(noise_tx->tx_block, out, block_size);
  for (int row = 0; row < block_size; ++row) {
    for (int col = 0; col < block_size; ++col) {
      data[row * block_size + col] =
          ((double)out[row * block_size + col]) / (255.0 * (1 << 8));
    }
  }
}

void aom_noise_tx_add_energy(const struct aom_noise_tx_t *noise_tx,
                             double *psd) {
  const int block_size = noise_tx->block_size;
  for (int yb = 0; yb < block_size; ++yb) {
    for (int xb = 0; xb < block_size; ++xb) {
      tran_low_t *c = noise_tx->tx_block + yb * block_size + xb;
      psd[yb * block_size + xb] += *c * *c;
    }
  }
}

void aom_noise_tx_free(struct aom_noise_tx_t *noise_tx) {
  if (!noise_tx) return;
  aom_free(noise_tx->tx_block);
  aom_free(noise_tx->block);
  memset(noise_tx, 0, sizeof(struct aom_noise_tx_t));
}
#endif
#endif  // CONFIG_DENOISE

double aom_normalized_cross_correlation(const double *a, const double *b,
                                        int n) {
  double c = 0;
  double a_len = 0;
  double b_len = 0;
  for (int i = 0; i < n; ++i) {
    a_len += a[i] * a[i];
    b_len += b[i] * b[i];
    c += a[i] * b[i];
  }
  return c / (sqrt(a_len) * sqrt(b_len));
}

void aom_noise_synth(int lag, int n, const int (*coords)[2],
                     const double *coeffs, double *data, int w, int h) {
  const int pad_size = 3 * lag;
  const int padded_w = w + pad_size;
  const int padded_h = h + pad_size;
  int x = 0, y = 0;
  double *padded = (double *)aom_malloc(padded_w * padded_h * sizeof(*padded));

  for (y = 0; y < padded_h; ++y) {
    for (x = 0; x < padded_w; ++x) {
      padded[y * padded_w + x] = aom_randn(1.0);
    }
  }
  for (y = lag; y < padded_h; ++y) {
    for (x = lag; x < padded_w; ++x) {
      double sum = 0;
      int i = 0;
      for (i = 0; i < n; ++i) {
        const int dx = coords[i][0];
        const int dy = coords[i][1];
        sum += padded[(y + dy) * padded_w + (x + dx)] * coeffs[i];
      }
      padded[y * padded_w + x] += sum;
    }
  }
  // Copy over the padded rows to the output
  for (y = 0; y < h; ++y) {
    memcpy(data + y * w, padded + y * padded_w, sizeof(*data) * w);
  }
  aom_free(padded);
}

#if CONFIG_DENOISE
int aom_noise_psd_get(const double *data, int w, int h, double *psd,
                      int block_size) {
  const int n = block_size * block_size;
  int num_blocks = 0;
  int y, x, row, col, i;
  struct aom_noise_tx_t *noise_tx = aom_noise_tx_malloc(block_size);
  const double scale = 1.;  // TODO(birkbeck): fix this with DCT
  double *input = 0;
  if (!noise_tx) {
    fprintf(stderr, "Unable to obtain noise_tx for block_size=%d\n",
            block_size);
    return 0;
  }
  input = aom_malloc(sizeof(double) * n);
  if (!input) {
    fprintf(stderr, "Unable to allcoate memory for input block\n");
    aom_noise_tx_free(noise_tx);
    aom_free(noise_tx);
    return 0;
  }
  memset(input, 0, sizeof(tran_high_t) * n);
  memset(psd, 0, sizeof(double) * block_size * block_size);

  for (y = 0; y < h - block_size; y += block_size / 2) {
    for (x = 0; x < w - block_size; x += block_size / 2) {
      for (row = 0; row < block_size; ++row) {
        for (col = 0; col < block_size; ++col) {
          input[row * block_size + col] =
              (scale * data[(y + row) * w + (x + col)]);
        }
      }
      aom_noise_tx_forward(noise_tx, input);
      aom_noise_tx_add_energy(noise_tx, psd);
      num_blocks++;
    }
  }
  for (i = 0; i < block_size * block_size; ++i) {
    psd[i] /= num_blocks;
  }

  aom_noise_tx_free(noise_tx);
  aom_free(noise_tx);
  aom_free(input);
  return 1;
}
#endif

double aom_noise_psd_compare(const double *a, const double *b, int block_size) {
  return aom_normalized_cross_correlation(a, b, block_size * block_size);
}

int aom_noise_data_validate(const double *data, int w, int h,
                            const double *ref_psd, int block_size) {
  const double kVarianceThreshold = 2;
  const double kMeanThreshold = 2;
  const double kCorrelationThreshold = 0.98;

  int x = 0, y = 0;
  int ret_value = 1;
  double var = 0, mean = 0, correlation = 0;
  double *mean_x, *mean_y, *var_x, *var_y;
  double *data_psd =
      (double *)aom_malloc(sizeof(double) * block_size * block_size);
#if CONFIG_DENOISE
  aom_noise_psd_get(data, w, h, data_psd, block_size);
  correlation = aom_noise_psd_compare(data_psd, ref_psd, block_size);
#endif
  // Check that noise variance is not increasing in x or y
  // and that the data is zero mean.
  mean_x = (double *)aom_malloc(sizeof(*mean_x) * w);
  var_x = (double *)aom_malloc(sizeof(*var_x) * w);
  mean_y = (double *)aom_malloc(sizeof(*mean_x) * h);
  var_y = (double *)aom_malloc(sizeof(*var_y) * h);

  memset(mean_x, 0, sizeof(*mean_x) * w);
  memset(var_x, 0, sizeof(*var_x) * w);
  memset(mean_y, 0, sizeof(*mean_y) * h);
  memset(var_y, 0, sizeof(*var_y) * h);

  for (y = 0; y < h; ++y) {
    for (x = 0; x < w; ++x) {
      const double d = data[y * w + x];
      var_x[x] += d * d;
      var_y[y] += d * d;
      mean_x[x] += d;
      mean_y[y] += d;
      var += d * d;
      mean += d;
    }
  }
  mean /= (w * h);
  var = var / (w * h) - mean * mean;

  for (y = 0; y < h; ++y) {
    mean_y[y] /= h;
    var_y[y] = var_y[y] / h - mean_y[y] * mean_y[y];
    if (fabs(var_y[y] - var) >= kVarianceThreshold) {
      fprintf(stderr, "Variance distance too large %f %f\n", var_y[y], var);
      ret_value = 0;
      break;
    }
    if (fabs(mean_y[y] - mean) >= kMeanThreshold) {
      fprintf(stderr, "Mean distance too large %f %f\n", mean_y[y], mean);
      ret_value = 0;
      break;
    }
  }

  for (x = 0; x < w; ++x) {
    mean_x[x] /= w;
    var_x[x] = var_x[x] / w - mean_x[x] * mean_x[x];
    if (fabs(var_x[x] - var) >= kVarianceThreshold) {
      fprintf(stderr, "Variance distance too large %f %f\n", var_x[x], var);
      ret_value = 0;
      break;
    }
    if (fabs(mean_x[x] - mean) >= kMeanThreshold) {
      fprintf(stderr, "Mean distance too large %f %f\n", mean_x[x], mean);
      ret_value = 0;
      break;
    }
  }
#ifdef CONFIG_DENOISE
  if (correlation < kCorrelationThreshold) {
    fprintf(stderr, "Correlation = %lf is too low (should be >= %lf)\n",
            correlation, kCorrelationThreshold);
    ret_value = 0;
  }
#endif

  aom_free(mean_x);
  aom_free(mean_y);
  aom_free(var_x);
  aom_free(var_y);
  aom_free(data_psd);

  return ret_value;
}
