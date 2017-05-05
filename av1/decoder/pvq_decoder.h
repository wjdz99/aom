/*
 * Copyright (c) 2001-2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/* clang-format off */

#if !defined(_pvq_decoder_H)
# define _pvq_decoder_H (1)
# include "aom_dsp/bitreader.h"
# include "aom_dsp/entdec.h"
# include "av1/common/pvq.h"
# include "av1/decoder/decint.h"

#define aom_read_symbol_pvq(r, cdf, nsymbs, ACCT_STR_NAME) \
  aom_read_symbol_pvq_(r, cdf, nsymbs ACCT_STR_ARG(ACCT_STR_NAME))

int aom_read_symbol_pvq_(AomReader *r, AomCdfProb *cdf, int nsymbs
  ACCT_STR_PARAM);

void aom_decode_band_pvq_splits(AomReader *r, OdPvqCodewordCtx *adapt,
 OdCoeff *y, int n, int k, int level);

#define aom_laplace_decode_special(r, decay, ACCT_STR_NAME) \
  aom_laplace_decode_special_(r, decay ACCT_STR_ARG(ACCT_STR_NAME))

int aom_laplace_decode_special_(AomReader *r, unsigned decay ACCT_STR_PARAM);

void od_pvq_decode(DaalaDecCtx *dec, OdCoeff *ref, OdCoeff *out, int q0,
    int pli, int bs, const OdVal16 *beta, int is_keyframe,
    unsigned int *flags, PvqSkipType ac_dc_coded, const int16_t *qm,
    const int16_t *qm_inv);

#endif
