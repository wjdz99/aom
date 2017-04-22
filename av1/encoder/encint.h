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

#if !defined(_encint_H)
# define _encint_H (1)

typedef struct DaalaEncCtx od_enc_ctx;
typedef struct od_params_ctx od_params_ctx;
typedef struct OdRollbackBuffer od_rollback_buffer;

# include "aom_dsp/entenc.h"
# include "av1/common/odintrin.h"
# include "av1/common/pvq_state.h"

struct DaalaEncCtx{
  /* Stores context-adaptive CDFs for PVQ. */
  OdState state;
  /* AOM entropy encoder. */
  aom_writer w;
  int use_activity_masking;
  /* Mode of quantization matrice : FLAT (0) or HVS (1) */
  int qm;
  /*Normalized PVQ lambda for use where we've already performed
     quantization.*/
  double pvq_norm_lambda;
  double pvq_norm_lambda_dc;
};

// from daalaenc.h
/**The encoder context.*/
typedef struct DaalaEncCtx daala_enc_ctx;

/** Holds important encoder information so we can roll back decisions */
struct OdRollbackBuffer {
  od_ec_enc ec;
  OdAdaptCtx adapt;
};

void od_encode_checkpoint(const DaalaEncCtx *enc, OdRollbackBuffer *rbuf);
void od_encode_rollback(DaalaEncCtx *enc, const OdRollbackBuffer *rbuf);

#endif
