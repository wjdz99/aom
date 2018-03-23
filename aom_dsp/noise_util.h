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

#ifndef AOM_DSP_NOISE_UTIL_H_
#define AOM_DSP_NOISE_UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

// aom_noise_tx_t is an abstraction of a transform that is used for denoising.
// It is meant to be lightweight and does hold the transformed data (as
// the user should not be manipulating the transformed data directly).
struct aom_noise_tx_t;

// Allocates and returns a aom_noise_tx_t useful for denoising the given
// block_size. The resulting aom_noise_tx_t should be free'd with
// aom_noise_tx_free.
struct aom_noise_tx_t *aom_noise_tx_malloc(int block_size);
void aom_noise_tx_free(struct aom_noise_tx_t *aom_noise_tx);

// Transforms the internal data and holds it the aom_noise_tx's internal buffer.
void aom_noise_tx_forward(struct aom_noise_tx_t *aom_noise_tx, double *data);

// Filters aom_noise_tx's internal data using the provided noise power spectral
// density. The PSD must be at least block_size * block_size and should be
// populated with a constant or via estimates taken from
// aom_noise_tx_add_energy.
void aom_noise_tx_filter(struct aom_noise_tx_t *aom_noise_tx, double *psd);

// Performs an inverse transform using the internal transform data.
void aom_noise_tx_inverse(struct aom_noise_tx_t *aom_noise_tx, double *data);

// Aggregates the power of the buffered transform data into the psd buffer.
void aom_noise_tx_add_energy(const struct aom_noise_tx_t *aom_noise_tx,
                             double *psd);

// Returns a default value suitable for denosing a transform of the given
// block_size.
double aom_noise_psd_get_default_value(int block_size, double factor);

// Computes normalized cross correlation of two vectors a and b of length n.
double aom_normalized_cross_correlation(const double *a, const double *b,
                                        int n);

// Synthesizes noise using the auto-regressive filter of the given lag,
// with the provided n coefficients sampled at the given coords.
void aom_noise_synth(int lag, int n, const int (*coords)[2],
                     const double *coeffs, double *data, int w, int h);

// Computes a noise_power spectral density (PSD) of the given data buffer of
// size w x h using a transform of the given block_size.
int aom_noise_psd_get(const double *data, int w, int h, double *psd,
                      int block_size);

// Compares two power spectral densities using normalized cross correlation
// and returns the block size.
double aom_noise_psd_compare(const double *a, const double *b, int block_size);

// Validates the correlated noise in the data buffer of size (w, h) against a
// reference noise PSD using the given block_size.
int aom_noise_data_validate(const double *data, int w, int h,
                            const double *ref_psd, int block_size);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#endif  // AOM_DSP_NOISE_UTIL_H_
