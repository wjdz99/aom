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

#if !defined(_pvq_H)
# define _pvq_H (1)
# include "generic_code.h"
# include "odintrin.h"

extern const uint16_t EXP_CDF_TABLE[][16];
extern const uint16_t LAPLACE_OFFSET[];

#if CONFIG_DAALA_DIST
#define AV1_PVQ_ENABLE_ACTIVITY_MASKING (1)
#else
#define AV1_PVQ_ENABLE_ACTIVITY_MASKING (0)
#endif

# define PVQ_MAX_PARTITIONS (1 + 3*(OD_TXSIZES-1))

# define OD_NOREF_ADAPT_SPEED (4)
/* Normalized lambda for PVQ quantizer. Since we normalize the gain by q, the
   distortion is normalized by q^2 and lambda does not need the q^2 factor.
   At high rate, this would be log(2)/6, but we're using a slightly more
   aggressive value, closer to:
   Li, Xiang, et al. "Laplace distribution based Lagrangian rate distortion
   optimization for hybrid video coding." Circuits and Systems for Video
   Technology, IEEE Transactions on 19.2 (2009): 193-205.
   */
# define OD_PVQ_LAMBDA (.1146)

#define OD_PVQ_SKIP_ZERO 1
#define OD_PVQ_SKIP_COPY 2

/* Maximum size for coding a PVQ band. */
#define OD_MAX_PVQ_SIZE (1024)

#if defined(OD_FLOAT_PVQ)
#define OD_QM_SHIFT (15)
#else
#define OD_QM_SHIFT (11)
#endif
#define OD_QM_SCALE (1 << OD_QM_SHIFT)
#if defined(OD_FLOAT_PVQ)
#define OD_QM_SCALE_1 (1./OD_QM_SCALE)
#endif
#define OD_QM_SCALE_MAX 32767
#define OD_QM_INV_SHIFT (12)
#define OD_QM_INV_SCALE (1 << OD_QM_INV_SHIFT)
#if defined(OD_FLOAT_PVQ)
#define OD_QM_INV_SCALE_1 (1./OD_QM_INV_SCALE)
#endif
#define OD_QM_OFFSET(bs) ((((1 << 2*bs) - 1) << 2*OD_LOG_BSIZE0)/3)
#define OD_QM_STRIDE (OD_QM_OFFSET(OD_TXSIZES))
#define OD_QM_BUFFER_SIZE (2*OD_QM_STRIDE)

#if !defined(OD_FLOAT_PVQ)
#define OD_THETA_SHIFT (15)
#define OD_THETA_SCALE ((1 << OD_THETA_SHIFT)*2./M_PI)
#define OD_MAX_THETA_SCALE (1 << OD_THETA_SHIFT)
#define OD_TRIG_SCALE (32768)
#define OD_BETA_SHIFT (12)
#define OD_BETA_SCALE_1 (1./(1 << OD_BETA_SHIFT))
/*Multiplies 16-bit a by 32-bit b and keeps bits [16:64-OD_BETA_SHIFT-1].*/
#define OD_MULT16_32_QBETA(a, b) \
 ((int16_t)(a)*(int64_t)(int32_t)(b) >> OD_BETA_SHIFT)
# define OD_MULT16_16_QBETA(a, b) \
  ((((int16_t)(a))*((int32_t)(int16_t)(b))) >> OD_BETA_SHIFT)
#define OD_CGAIN_SHIFT (8)
#define OD_CGAIN_SCALE (1 << OD_CGAIN_SHIFT)
#else
#define OD_BETA_SCALE_1 (1.)
#define OD_THETA_SCALE (1)
#define OD_TRIG_SCALE (1)
#define OD_CGAIN_SCALE (1)
#endif
#define OD_THETA_SCALE_1 (1./OD_THETA_SCALE)
#define OD_TRIG_SCALE_1 (1./OD_TRIG_SCALE)
#define OD_CGAIN_SCALE_1 (1./OD_CGAIN_SCALE)
#define OD_CGAIN_SCALE_2 (OD_CGAIN_SCALE_1*OD_CGAIN_SCALE_1)

/* Largest PVQ partition is half the coefficients of largest block size. */
#define MAXN (OD_TXSIZE_MAX*OD_TXSIZE_MAX/2)

#define OD_COMPAND_SHIFT (8 + OD_COEFF_SHIFT)
#define OD_COMPAND_SCALE (1 << OD_COMPAND_SHIFT)
#define OD_COMPAND_SCALE_1 (1./OD_COMPAND_SCALE)

#define OD_QM_SIZE (OD_TXSIZES*(OD_TXSIZES + 1))

#define OD_FLAT_QM 0
#define OD_HVS_QM  1

# define OD_NSB_ADAPT_CTXS (4)

# define OD_ADAPT_K_Q8        0
# define OD_ADAPT_SUM_EX_Q8   1
# define OD_ADAPT_COUNT_Q8    2
# define OD_ADAPT_COUNT_EX_Q8 3

# define OD_ADAPT_NO_VALUE (-2147483647-1)

typedef enum {
  PVQ_SKIP = 0x0,
  DC_CODED = 0x1,
  AC_CODED = 0x2,
  AC_DC_CODED = 0x3,
} PvqSkipType;

typedef struct OdPvqAdaptCtx  OdPvqAdaptCtx;
typedef struct OdPvqCodewordCtx OdPvqCodewordCtx;

struct OdPvqCodewordCtx {
  int                 pvq_adapt[2*OD_TXSIZES*OD_NSB_ADAPT_CTXS];
  /* CDFs are size 16 despite the fact that we're using less than that. */
  uint16_t            pvq_k1_cdf[12][CDF_SIZE(16)];
  uint16_t            pvq_split_cdf[22*7][CDF_SIZE(8)];
};

struct OdPvqAdaptCtx {
  OdPvqCodewordCtx pvq_codeword_ctx;
  GenericEncoder     pvq_param_model[3];
  int                 pvq_ext[OD_TXSIZES*PVQ_MAX_PARTITIONS];
  int                 pvq_exg[OD_NPLANES_MAX][OD_TXSIZES][PVQ_MAX_PARTITIONS];
  uint16_t pvq_gaintheta_cdf[2*OD_TXSIZES*PVQ_MAX_PARTITIONS][CDF_SIZE(16)];
  uint16_t pvq_skip_dir_cdf[2*(OD_TXSIZES-1)][CDF_SIZE(7)];
};

typedef struct OdQmEntry {
  int interp_q;
  int scale_q8;
  const unsigned char *qm_q4;
} OdQmEntry;

extern const OdQmEntry OD_DEFAULT_QMS[2][2][OD_NPLANES_MAX];

void od_adapt_pvq_ctx_reset(OdPvqAdaptCtx *state, int is_keyframe);
int od_pvq_size_ctx(int n);
int od_pvq_k1_ctx(int n, int orig_size);

OdVal16 od_pvq_sin(OdVal32 x);
OdVal16 od_pvq_cos(OdVal32 x);
#if !defined(OD_FLOAT_PVQ)
int od_vector_log_mag(const OdCoeff *x, int n);
#endif

void od_interp_qm(unsigned char *out, int q, const OdQmEntry *entry1,
                  const OdQmEntry *entry2);

int od_qm_get_index(int bs, int band);

extern const OdVal16 *const OD_PVQ_BETA[2][OD_NPLANES_MAX][OD_TXSIZES + 1];

void od_init_qm(int16_t *x, int16_t *x_inv, const int *qm);
int od_compute_householder(OdVal16 *r, int n, OdVal32 gr, int *sign,
 int shift);
void od_apply_householder(OdVal16 *out, const OdVal16 *x, const OdVal16 *r,
 int n);
void od_pvq_synthesis_partial(OdCoeff *xcoeff, const OdCoeff *ypulse,
                                  const OdVal16 *r, int n,
                                  int noref, OdVal32 g,
                                  OdVal32 theta, int m, int s,
                                  const int16_t *qm_inv);
OdVal32 od_gain_expand(OdVal32 cg, int q0, OdVal16 beta);
OdVal32 od_pvq_compute_gain(const OdVal16 *x, int n, int q0, OdVal32 *g,
 OdVal16 beta, int bshift);
int od_pvq_compute_max_theta(OdVal32 qcg, OdVal16 beta);
OdVal32 od_pvq_compute_theta(int t, int max_theta);
int od_pvq_compute_k(OdVal32 qcg, int itheta, int noref, int n, OdVal16 beta);

int od_vector_is_null(const OdCoeff *x, int len);
int od_qm_offset(int bs, int xydec);

#endif
