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
#include <math.h>
#include <stdio.h>

#include "config/av1_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/bitops.h"
#include "aom_ports/mem.h"
#include "aom_ports/system_state.h"

#include "av1/common/common.h"
#include "av1/common/entropy.h"
#include "av1/common/entropymode.h"
#include "av1/common/mvref_common.h"
#include "av1/common/pred_common.h"
#include "av1/common/quant_common.h"
#include "av1/common/reconinter.h"
#include "av1/common/reconintra.h"
#include "av1/common/seg_common.h"

#include "av1/encoder/av1_quantize.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/encodemb.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/tokenize.h"

#define RD_THRESH_POW 1.25

// The baseline rd thresholds for breaking out of the rd loop for
// certain modes are assumed to be based on 8x8 blocks.
// This table is used to correct for block size.
// The factors here are << 2 (2 = x0.5, 32 = x8 etc).
static const uint8_t rd_thresh_block_size_factor[BLOCK_SIZES_ALL] = {
  2, 3, 3, 4, 6, 6, 8, 12, 12, 16, 24, 24, 32, 48, 48, 64, 4, 4, 8, 8, 16, 16
};

static const int use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][EXT_TX_SIZES] =
    {
      { 1, 1, 1, 1 },  // unused
      { 1, 1, 0, 0 },
      { 0, 0, 1, 0 },
    };

static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][EXT_TX_SIZES] =
    {
      { 1, 1, 1, 1 },  // unused
      { 1, 1, 0, 0 },
      { 0, 0, 1, 0 },
      { 0, 0, 0, 1 },
    };

static const int av1_ext_tx_set_idx_to_type[2][AOMMAX(EXT_TX_SETS_INTRA,
                                                      EXT_TX_SETS_INTER)] = {
  {
      // Intra
      EXT_TX_SET_DCTONLY,
      EXT_TX_SET_DTT4_IDTX_1DDCT,
      EXT_TX_SET_DTT4_IDTX,
  },
  {
      // Inter
      EXT_TX_SET_DCTONLY,
      EXT_TX_SET_ALL16,
      EXT_TX_SET_DTT9_IDTX_1DDCT,
      EXT_TX_SET_DCT_IDTX,
  },
};

void av1_fill_mode_rates(AV1_COMMON *const cm, MACROBLOCK *x,
                         FRAME_CONTEXT *fc) {
  int i, j;

  for (i = 0; i < PARTITION_CONTEXTS; ++i)
    av1_cost_tokens_from_cdf(x->partition_cost[i], fc->partition_cdf[i], NULL);

  if (cm->current_frame.skip_mode_info.skip_mode_flag) {
    for (i = 0; i < SKIP_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->skip_mode_cost[i], fc->skip_mode_cdfs[i],
                               NULL);
    }
  }

  for (i = 0; i < SKIP_CONTEXTS; ++i) {
    av1_cost_tokens_from_cdf(x->skip_cost[i], fc->skip_cdfs[i], NULL);
  }

  for (i = 0; i < KF_MODE_CONTEXTS; ++i)
    for (j = 0; j < KF_MODE_CONTEXTS; ++j)
      av1_cost_tokens_from_cdf(x->y_mode_costs[i][j], fc->kf_y_cdf[i][j], NULL);

  for (i = 0; i < BLOCK_SIZE_GROUPS; ++i)
    av1_cost_tokens_from_cdf(x->mbmode_cost[i], fc->y_mode_cdf[i], NULL);
  for (i = 0; i < CFL_ALLOWED_TYPES; ++i)
    for (j = 0; j < INTRA_MODES; ++j)
      av1_cost_tokens_from_cdf(x->intra_uv_mode_cost[i][j],
                               fc->uv_mode_cdf[i][j], NULL);

  av1_cost_tokens_from_cdf(x->filter_intra_mode_cost, fc->filter_intra_mode_cdf,
                           NULL);
  for (i = 0; i < BLOCK_SIZES_ALL; ++i) {
    if (av1_filter_intra_allowed_bsize(cm, i))
      av1_cost_tokens_from_cdf(x->filter_intra_cost[i],
                               fc->filter_intra_cdfs[i], NULL);
  }

  for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
    av1_cost_tokens_from_cdf(x->switchable_interp_costs[i],
                             fc->switchable_interp_cdf[i], NULL);

  for (i = 0; i < PALATTE_BSIZE_CTXS; ++i) {
    av1_cost_tokens_from_cdf(x->palette_y_size_cost[i],
                             fc->palette_y_size_cdf[i], NULL);
    av1_cost_tokens_from_cdf(x->palette_uv_size_cost[i],
                             fc->palette_uv_size_cdf[i], NULL);
    for (j = 0; j < PALETTE_Y_MODE_CONTEXTS; ++j) {
      av1_cost_tokens_from_cdf(x->palette_y_mode_cost[i][j],
                               fc->palette_y_mode_cdf[i][j], NULL);
    }
  }

  for (i = 0; i < PALETTE_UV_MODE_CONTEXTS; ++i) {
    av1_cost_tokens_from_cdf(x->palette_uv_mode_cost[i],
                             fc->palette_uv_mode_cdf[i], NULL);
  }

  for (i = 0; i < PALETTE_SIZES; ++i) {
    for (j = 0; j < PALETTE_COLOR_INDEX_CONTEXTS; ++j) {
      av1_cost_tokens_from_cdf(x->palette_y_color_cost[i][j],
                               fc->palette_y_color_index_cdf[i][j], NULL);
      av1_cost_tokens_from_cdf(x->palette_uv_color_cost[i][j],
                               fc->palette_uv_color_index_cdf[i][j], NULL);
    }
  }

  int sign_cost[CFL_JOINT_SIGNS];
  av1_cost_tokens_from_cdf(sign_cost, fc->cfl_sign_cdf, NULL);
  for (int joint_sign = 0; joint_sign < CFL_JOINT_SIGNS; joint_sign++) {
    int *cost_u = x->cfl_cost[joint_sign][CFL_PRED_U];
    int *cost_v = x->cfl_cost[joint_sign][CFL_PRED_V];
    if (CFL_SIGN_U(joint_sign) == CFL_SIGN_ZERO) {
      memset(cost_u, 0, CFL_ALPHABET_SIZE * sizeof(*cost_u));
    } else {
      const aom_cdf_prob *cdf_u = fc->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
      av1_cost_tokens_from_cdf(cost_u, cdf_u, NULL);
    }
    if (CFL_SIGN_V(joint_sign) == CFL_SIGN_ZERO) {
      memset(cost_v, 0, CFL_ALPHABET_SIZE * sizeof(*cost_v));
    } else {
      const aom_cdf_prob *cdf_v = fc->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];
      av1_cost_tokens_from_cdf(cost_v, cdf_v, NULL);
    }
    for (int u = 0; u < CFL_ALPHABET_SIZE; u++)
      cost_u[u] += sign_cost[joint_sign];
  }

  for (i = 0; i < MAX_TX_CATS; ++i)
    for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
      av1_cost_tokens_from_cdf(x->tx_size_cost[i][j], fc->tx_size_cdf[i][j],
                               NULL);

  for (i = 0; i < TXFM_PARTITION_CONTEXTS; ++i) {
    av1_cost_tokens_from_cdf(x->txfm_partition_cost[i],
                             fc->txfm_partition_cdf[i], NULL);
  }

  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    int s;
    for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
      if (use_inter_ext_tx_for_txsize[s][i]) {
        av1_cost_tokens_from_cdf(
            x->inter_tx_type_costs[s][i], fc->inter_ext_tx_cdf[s][i],
            av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[1][s]]);
      }
    }
    for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
      if (use_intra_ext_tx_for_txsize[s][i]) {
        for (j = 0; j < INTRA_MODES; ++j) {
          av1_cost_tokens_from_cdf(
              x->intra_tx_type_costs[s][i][j], fc->intra_ext_tx_cdf[s][i][j],
              av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[0][s]]);
        }
      }
    }
  }
  for (i = 0; i < DIRECTIONAL_MODES; ++i) {
    av1_cost_tokens_from_cdf(x->angle_delta_cost[i], fc->angle_delta_cdf[i],
                             NULL);
  }
  av1_cost_tokens_from_cdf(x->switchable_restore_cost,
                           fc->switchable_restore_cdf, NULL);
  av1_cost_tokens_from_cdf(x->wiener_restore_cost, fc->wiener_restore_cdf,
                           NULL);
  av1_cost_tokens_from_cdf(x->sgrproj_restore_cost, fc->sgrproj_restore_cdf,
                           NULL);
  av1_cost_tokens_from_cdf(x->intrabc_cost, fc->intrabc_cdf, NULL);

  if (!frame_is_intra_only(cm)) {
    for (i = 0; i < COMP_INTER_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->comp_inter_cost[i], fc->comp_inter_cdf[i],
                               NULL);
    }

    for (i = 0; i < REF_CONTEXTS; ++i) {
      for (j = 0; j < SINGLE_REFS - 1; ++j) {
        av1_cost_tokens_from_cdf(x->single_ref_cost[i][j],
                                 fc->single_ref_cdf[i][j], NULL);
      }
    }

    for (i = 0; i < COMP_REF_TYPE_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->comp_ref_type_cost[i],
                               fc->comp_ref_type_cdf[i], NULL);
    }

    for (i = 0; i < UNI_COMP_REF_CONTEXTS; ++i) {
      for (j = 0; j < UNIDIR_COMP_REFS - 1; ++j) {
        av1_cost_tokens_from_cdf(x->uni_comp_ref_cost[i][j],
                                 fc->uni_comp_ref_cdf[i][j], NULL);
      }
    }

    for (i = 0; i < REF_CONTEXTS; ++i) {
      for (j = 0; j < FWD_REFS - 1; ++j) {
        av1_cost_tokens_from_cdf(x->comp_ref_cost[i][j], fc->comp_ref_cdf[i][j],
                                 NULL);
      }
    }

    for (i = 0; i < REF_CONTEXTS; ++i) {
      for (j = 0; j < BWD_REFS - 1; ++j) {
        av1_cost_tokens_from_cdf(x->comp_bwdref_cost[i][j],
                                 fc->comp_bwdref_cdf[i][j], NULL);
      }
    }

    for (i = 0; i < INTRA_INTER_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->intra_inter_cost[i], fc->intra_inter_cdf[i],
                               NULL);
    }

    for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->newmv_mode_cost[i], fc->newmv_cdf[i], NULL);
    }

    for (i = 0; i < GLOBALMV_MODE_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->zeromv_mode_cost[i], fc->zeromv_cdf[i], NULL);
    }

    for (i = 0; i < REFMV_MODE_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->refmv_mode_cost[i], fc->refmv_cdf[i], NULL);
    }

    for (i = 0; i < DRL_MODE_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->drl_mode_cost0[i], fc->drl_cdf[i], NULL);
    }
    for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
      av1_cost_tokens_from_cdf(x->inter_compound_mode_cost[i],
                               fc->inter_compound_mode_cdf[i], NULL);
    for (i = 0; i < BLOCK_SIZES_ALL; ++i)
      av1_cost_tokens_from_cdf(x->compound_type_cost[i],
                               fc->compound_type_cdf[i], NULL);
    for (i = 0; i < BLOCK_SIZES_ALL; ++i) {
      if (get_interinter_wedge_bits(i)) {
        av1_cost_tokens_from_cdf(x->wedge_idx_cost[i], fc->wedge_idx_cdf[i],
                                 NULL);
      }
    }
    for (i = 0; i < BLOCK_SIZE_GROUPS; ++i) {
      av1_cost_tokens_from_cdf(x->interintra_cost[i], fc->interintra_cdf[i],
                               NULL);
      av1_cost_tokens_from_cdf(x->interintra_mode_cost[i],
                               fc->interintra_mode_cdf[i], NULL);
    }
    for (i = 0; i < BLOCK_SIZES_ALL; ++i) {
      av1_cost_tokens_from_cdf(x->wedge_interintra_cost[i],
                               fc->wedge_interintra_cdf[i], NULL);
    }
    for (i = BLOCK_8X8; i < BLOCK_SIZES_ALL; i++) {
      av1_cost_tokens_from_cdf(x->motion_mode_cost[i], fc->motion_mode_cdf[i],
                               NULL);
    }
    for (i = BLOCK_8X8; i < BLOCK_SIZES_ALL; i++) {
      av1_cost_tokens_from_cdf(x->motion_mode_cost1[i], fc->obmc_cdf[i], NULL);
    }
    for (i = 0; i < COMP_INDEX_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->comp_idx_cost[i], fc->compound_index_cdf[i],
                               NULL);
    }
    for (i = 0; i < COMP_GROUP_IDX_CONTEXTS; ++i) {
      av1_cost_tokens_from_cdf(x->comp_group_idx_cost[i],
                               fc->comp_group_idx_cdf[i], NULL);
    }
  }
}

// Values are now correlated to quantizer.
static int sad_per_bit16lut_8[QINDEX_RANGE];
static int sad_per_bit4lut_8[QINDEX_RANGE];
static int sad_per_bit16lut_10[QINDEX_RANGE];
static int sad_per_bit4lut_10[QINDEX_RANGE];
static int sad_per_bit16lut_12[QINDEX_RANGE];
static int sad_per_bit4lut_12[QINDEX_RANGE];

static void init_me_luts_bd(int *bit16lut, int *bit4lut, int range,
                            aom_bit_depth_t bit_depth) {
  int i;
  // Initialize the sad lut tables using a formulaic calculation for now.
  // This is to make it easier to resolve the impact of experimental changes
  // to the quantizer tables.
  for (i = 0; i < range; i++) {
    const double q = av1_convert_qindex_to_q(i, bit_depth);
    bit16lut[i] = (int)(0.0418 * q + 2.4107);
    bit4lut[i] = (int)(0.063 * q + 2.742);
  }
}

void av1_init_me_luts(void) {
  init_me_luts_bd(sad_per_bit16lut_8, sad_per_bit4lut_8, QINDEX_RANGE,
                  AOM_BITS_8);
  init_me_luts_bd(sad_per_bit16lut_10, sad_per_bit4lut_10, QINDEX_RANGE,
                  AOM_BITS_10);
  init_me_luts_bd(sad_per_bit16lut_12, sad_per_bit4lut_12, QINDEX_RANGE,
                  AOM_BITS_12);
}

static const int rd_boost_factor[16] = { 64, 32, 32, 32, 24, 16, 12, 12,
                                         8,  8,  4,  4,  2,  2,  1,  0 };
static const int rd_frame_type_factor[FRAME_UPDATE_TYPES] = {
  128, 144, 128, 128, 144,
  // TODO(zoeliu): To adjust further following factor values.
  128, 128, 128,
  // TODO(weitinglin): We should investigate if the values should be the same
  //                   as the value used by OVERLAY frame
  144,  // INTNL_OVERLAY_UPDATE
  128   // INTNL_ARF_UPDATE
};

int av1_compute_rd_mult_based_on_qindex(const AV1_COMP *cpi, int qindex) {
  const int q = av1_dc_quant_Q3(qindex, 0, cpi->common.seq_params.bit_depth);
  int rdmult = q * q;
  rdmult = rdmult * 3 + (rdmult * 2 / 3);
  switch (cpi->common.seq_params.bit_depth) {
    case AOM_BITS_8: break;
    case AOM_BITS_10: rdmult = ROUND_POWER_OF_TWO(rdmult, 4); break;
    case AOM_BITS_12: rdmult = ROUND_POWER_OF_TWO(rdmult, 8); break;
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
  return rdmult > 0 ? rdmult : 1;
}

int av1_compute_rd_mult(const AV1_COMP *cpi, int qindex) {
  int64_t rdmult = av1_compute_rd_mult_based_on_qindex(cpi, qindex);
  if (cpi->oxcf.pass == 2 &&
      (cpi->common.current_frame.frame_type != KEY_FRAME)) {
    const GF_GROUP *const gf_group = &cpi->twopass.gf_group;
    const FRAME_UPDATE_TYPE frame_type = gf_group->update_type[gf_group->index];
    const int boost_index = AOMMIN(15, (cpi->rc.gfu_boost / 100));

    rdmult = (rdmult * rd_frame_type_factor[frame_type]) >> 7;
    rdmult += ((rdmult * rd_boost_factor[boost_index]) >> 7);
  }
  return (int)rdmult;
}

int av1_get_adaptive_rdmult(const AV1_COMP *cpi, double beta) {
  const AV1_COMMON *cm = &cpi->common;
  int64_t q =
      av1_dc_quant_Q3(cm->base_qindex, 0, cpi->common.seq_params.bit_depth);
  int64_t rdmult = 0;

  switch (cpi->common.seq_params.bit_depth) {
    case AOM_BITS_8: rdmult = (int)((88 * q * q / beta) / 24); break;
    case AOM_BITS_10:
      rdmult = ROUND_POWER_OF_TWO((int)((88 * q * q / beta) / 24), 4);
      break;
    default:
      assert(cpi->common.seq_params.bit_depth == AOM_BITS_12);
      rdmult = ROUND_POWER_OF_TWO((int)((88 * q * q / beta) / 24), 8);
      break;
  }

  if (cpi->oxcf.pass == 2 &&
      (cpi->common.current_frame.frame_type != KEY_FRAME)) {
    const GF_GROUP *const gf_group = &cpi->twopass.gf_group;
    const FRAME_UPDATE_TYPE frame_type = gf_group->update_type[gf_group->index];
    const int boost_index = AOMMIN(15, (cpi->rc.gfu_boost / 100));

    rdmult = (rdmult * rd_frame_type_factor[frame_type]) >> 7;
    rdmult += ((rdmult * rd_boost_factor[boost_index]) >> 7);
  }
  if (rdmult < 1) rdmult = 1;
  return (int)rdmult;
}

static int compute_rd_thresh_factor(int qindex, aom_bit_depth_t bit_depth) {
  double q;
  switch (bit_depth) {
    case AOM_BITS_8: q = av1_dc_quant_Q3(qindex, 0, AOM_BITS_8) / 4.0; break;
    case AOM_BITS_10: q = av1_dc_quant_Q3(qindex, 0, AOM_BITS_10) / 16.0; break;
    case AOM_BITS_12: q = av1_dc_quant_Q3(qindex, 0, AOM_BITS_12) / 64.0; break;
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
  // TODO(debargha): Adjust the function below.
  return AOMMAX((int)(pow(q, RD_THRESH_POW) * 5.12), 8);
}

void av1_initialize_me_consts(const AV1_COMP *cpi, MACROBLOCK *x, int qindex) {
  switch (cpi->common.seq_params.bit_depth) {
    case AOM_BITS_8:
      x->sadperbit16 = sad_per_bit16lut_8[qindex];
      x->sadperbit4 = sad_per_bit4lut_8[qindex];
      break;
    case AOM_BITS_10:
      x->sadperbit16 = sad_per_bit16lut_10[qindex];
      x->sadperbit4 = sad_per_bit4lut_10[qindex];
      break;
    case AOM_BITS_12:
      x->sadperbit16 = sad_per_bit16lut_12[qindex];
      x->sadperbit4 = sad_per_bit4lut_12[qindex];
      break;
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
  }
}

static void set_block_thresholds(const AV1_COMMON *cm, RD_OPT *rd) {
  int i, bsize, segment_id;

  for (segment_id = 0; segment_id < MAX_SEGMENTS; ++segment_id) {
    const int qindex =
        clamp(av1_get_qindex(&cm->seg, segment_id, cm->base_qindex) +
                  cm->y_dc_delta_q,
              0, MAXQ);
    const int q = compute_rd_thresh_factor(qindex, cm->seq_params.bit_depth);

    for (bsize = 0; bsize < BLOCK_SIZES_ALL; ++bsize) {
      // Threshold here seems unnecessarily harsh but fine given actual
      // range of values used for cpi->sf.thresh_mult[].
      const int t = q * rd_thresh_block_size_factor[bsize];
      const int thresh_max = INT_MAX / t;

      for (i = 0; i < MAX_MODES; ++i)
        rd->threshes[segment_id][bsize][i] = rd->thresh_mult[i] < thresh_max
                                                 ? rd->thresh_mult[i] * t / 4
                                                 : INT_MAX;
    }
  }
}

void av1_fill_coeff_costs(MACROBLOCK *x, FRAME_CONTEXT *fc,
                          const int num_planes) {
  const int nplanes = AOMMIN(num_planes, PLANE_TYPES);
  for (int eob_multi_size = 0; eob_multi_size < 7; ++eob_multi_size) {
    for (int plane = 0; plane < nplanes; ++plane) {
      LV_MAP_EOB_COST *pcost = &x->eob_costs[eob_multi_size][plane];

      for (int ctx = 0; ctx < 2; ++ctx) {
        aom_cdf_prob *pcdf;
        switch (eob_multi_size) {
          case 0: pcdf = fc->eob_flag_cdf16[plane][ctx]; break;
          case 1: pcdf = fc->eob_flag_cdf32[plane][ctx]; break;
          case 2: pcdf = fc->eob_flag_cdf64[plane][ctx]; break;
          case 3: pcdf = fc->eob_flag_cdf128[plane][ctx]; break;
          case 4: pcdf = fc->eob_flag_cdf256[plane][ctx]; break;
          case 5: pcdf = fc->eob_flag_cdf512[plane][ctx]; break;
          case 6:
          default: pcdf = fc->eob_flag_cdf1024[plane][ctx]; break;
        }
        av1_cost_tokens_from_cdf(pcost->eob_cost[ctx], pcdf, NULL);
      }
    }
  }
  for (int tx_size = 0; tx_size < TX_SIZES; ++tx_size) {
    for (int plane = 0; plane < nplanes; ++plane) {
      LV_MAP_COEFF_COST *pcost = &x->coeff_costs[tx_size][plane];

      for (int ctx = 0; ctx < TXB_SKIP_CONTEXTS; ++ctx)
        av1_cost_tokens_from_cdf(pcost->txb_skip_cost[ctx],
                                 fc->txb_skip_cdf[tx_size][ctx], NULL);

      for (int ctx = 0; ctx < SIG_COEF_CONTEXTS_EOB; ++ctx)
        av1_cost_tokens_from_cdf(pcost->base_eob_cost[ctx],
                                 fc->coeff_base_eob_cdf[tx_size][plane][ctx],
                                 NULL);
      for (int ctx = 0; ctx < SIG_COEF_CONTEXTS; ++ctx)
        av1_cost_tokens_from_cdf(pcost->base_cost[ctx],
                                 fc->coeff_base_cdf[tx_size][plane][ctx], NULL);

      for (int ctx = 0; ctx < EOB_COEF_CONTEXTS; ++ctx)
        av1_cost_tokens_from_cdf(pcost->eob_extra_cost[ctx],
                                 fc->eob_extra_cdf[tx_size][plane][ctx], NULL);

      for (int ctx = 0; ctx < DC_SIGN_CONTEXTS; ++ctx)
        av1_cost_tokens_from_cdf(pcost->dc_sign_cost[ctx],
                                 fc->dc_sign_cdf[plane][ctx], NULL);

      for (int ctx = 0; ctx < LEVEL_CONTEXTS; ++ctx) {
        int br_rate[BR_CDF_SIZE];
        int prev_cost = 0;
        int i, j;
        av1_cost_tokens_from_cdf(br_rate, fc->coeff_br_cdf[tx_size][plane][ctx],
                                 NULL);
        // printf("br_rate: ");
        // for(j = 0; j < BR_CDF_SIZE; j++)
        //  printf("%4d ", br_rate[j]);
        // printf("\n");
        for (i = 0; i < COEFF_BASE_RANGE; i += BR_CDF_SIZE - 1) {
          for (j = 0; j < BR_CDF_SIZE - 1; j++) {
            pcost->lps_cost[ctx][i + j] = prev_cost + br_rate[j];
          }
          prev_cost += br_rate[j];
        }
        pcost->lps_cost[ctx][i] = prev_cost;
        // printf("lps_cost: %d %d %2d : ", tx_size, plane, ctx);
        // for (i = 0; i <= COEFF_BASE_RANGE; i++)
        //  printf("%5d ", pcost->lps_cost[ctx][i]);
        // printf("\n");
      }
    }
  }
}

void av1_initialize_cost_tables(const AV1_COMMON *const cm, MACROBLOCK *x) {
  if (cm->cur_frame_force_integer_mv) {
    av1_build_nmv_cost_table(x->nmv_vec_cost, x->nmvcost, &cm->fc->nmvc,
                             MV_SUBPEL_NONE);
  } else {
    av1_build_nmv_cost_table(
        x->nmv_vec_cost,
        cm->allow_high_precision_mv ? x->nmvcost_hp : x->nmvcost, &cm->fc->nmvc,
        cm->allow_high_precision_mv);
  }
}

void av1_initialize_rd_consts(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->td.mb;
  RD_OPT *const rd = &cpi->rd;

  aom_clear_system_state();

  rd->RDMULT = av1_compute_rd_mult(cpi, cm->base_qindex + cm->y_dc_delta_q);

  set_error_per_bit(x, rd->RDMULT);

  set_block_thresholds(cm, rd);

  av1_initialize_cost_tables(cm, x);

  if (frame_is_intra_only(cm) && cm->allow_screen_content_tools &&
      cpi->oxcf.pass != 1) {
    int *dvcost[2] = { &cpi->dv_cost[0][MV_MAX], &cpi->dv_cost[1][MV_MAX] };
    av1_build_nmv_cost_table(cpi->dv_joint_cost, dvcost, &cm->fc->ndvc,
                             MV_SUBPEL_NONE);
  }

  if (cpi->oxcf.pass != 1) {
    for (int i = 0; i < TRANS_TYPES; ++i)
      // IDENTITY: 1 bit
      // TRANSLATION: 3 bits
      // ROTZOOM: 2 bits
      // AFFINE: 3 bits
      cpi->gmtype_cost[i] = (1 + (i > 0 ? (i == ROTZOOM ? 1 : 2) : 0))
                            << AV1_PROB_COST_SHIFT;
  }
}

static void model_rd_norm(int xsq_q10, int *r_q10, int *d_q10) {
  // NOTE: The tables below must be of the same size.

  // The functions described below are sampled at the four most significant
  // bits of x^2 + 8 / 256.

  // Normalized rate:
  // This table models the rate for a Laplacian source with given variance
  // when quantized with a uniform quantizer with given stepsize. The
  // closed form expression is:
  // Rn(x) = H(sqrt(r)) + sqrt(r)*[1 + H(r)/(1 - r)],
  // where r = exp(-sqrt(2) * x) and x = qpstep / sqrt(variance),
  // and H(x) is the binary entropy function.
  static const int rate_tab_q10[] = {
    65536, 6086, 5574, 5275, 5063, 4899, 4764, 4651, 4553, 4389, 4255, 4142,
    4044,  3958, 3881, 3811, 3748, 3635, 3538, 3453, 3376, 3307, 3244, 3186,
    3133,  3037, 2952, 2877, 2809, 2747, 2690, 2638, 2589, 2501, 2423, 2353,
    2290,  2232, 2179, 2130, 2084, 2001, 1928, 1862, 1802, 1748, 1698, 1651,
    1608,  1530, 1460, 1398, 1342, 1290, 1243, 1199, 1159, 1086, 1021, 963,
    911,   864,  821,  781,  745,  680,  623,  574,  530,  490,  455,  424,
    395,   345,  304,  269,  239,  213,  190,  171,  154,  126,  104,  87,
    73,    61,   52,   44,   38,   28,   21,   16,   12,   10,   8,    6,
    5,     3,    2,    1,    1,    1,    0,    0,
  };
  // Normalized distortion:
  // This table models the normalized distortion for a Laplacian source
  // with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expression is:
  // Dn(x) = 1 - 1/sqrt(2) * x / sinh(x/sqrt(2))
  // where x = qpstep / sqrt(variance).
  // Note the actual distortion is Dn * variance.
  static const int dist_tab_q10[] = {
    0,    0,    1,    1,    1,    2,    2,    2,    3,    3,    4,    5,
    5,    6,    7,    7,    8,    9,    11,   12,   13,   15,   16,   17,
    18,   21,   24,   26,   29,   31,   34,   36,   39,   44,   49,   54,
    59,   64,   69,   73,   78,   88,   97,   106,  115,  124,  133,  142,
    151,  167,  184,  200,  215,  231,  245,  260,  274,  301,  327,  351,
    375,  397,  418,  439,  458,  495,  528,  559,  587,  613,  637,  659,
    680,  717,  749,  777,  801,  823,  842,  859,  874,  899,  919,  936,
    949,  960,  969,  977,  983,  994,  1001, 1006, 1010, 1013, 1015, 1017,
    1018, 1020, 1022, 1022, 1023, 1023, 1023, 1024,
  };
  static const int xsq_iq_q10[] = {
    0,      4,      8,      12,     16,     20,     24,     28,     32,
    40,     48,     56,     64,     72,     80,     88,     96,     112,
    128,    144,    160,    176,    192,    208,    224,    256,    288,
    320,    352,    384,    416,    448,    480,    544,    608,    672,
    736,    800,    864,    928,    992,    1120,   1248,   1376,   1504,
    1632,   1760,   1888,   2016,   2272,   2528,   2784,   3040,   3296,
    3552,   3808,   4064,   4576,   5088,   5600,   6112,   6624,   7136,
    7648,   8160,   9184,   10208,  11232,  12256,  13280,  14304,  15328,
    16352,  18400,  20448,  22496,  24544,  26592,  28640,  30688,  32736,
    36832,  40928,  45024,  49120,  53216,  57312,  61408,  65504,  73696,
    81888,  90080,  98272,  106464, 114656, 122848, 131040, 147424, 163808,
    180192, 196576, 212960, 229344, 245728,
  };
  const int tmp = (xsq_q10 >> 2) + 8;
  const int k = get_msb(tmp) - 3;
  const int xq = (k << 3) + ((tmp >> k) & 0x7);
  const int one_q10 = 1 << 10;
  const int a_q10 = ((xsq_q10 - xsq_iq_q10[xq]) << 10) >> (2 + k);
  const int b_q10 = one_q10 - a_q10;
  *r_q10 = (rate_tab_q10[xq] * b_q10 + rate_tab_q10[xq + 1] * a_q10) >> 10;
  *d_q10 = (dist_tab_q10[xq] * b_q10 + dist_tab_q10[xq + 1] * a_q10) >> 10;
}

void av1_model_rd_from_var_lapndz(int64_t var, unsigned int n_log2,
                                  unsigned int qstep, int *rate,
                                  int64_t *dist) {
  // This function models the rate and distortion for a Laplacian
  // source with given variance when quantized with a uniform quantizer
  // with given stepsize. The closed form expressions are in:
  // Hang and Chen, "Source Model for transform video coder and its
  // application - Part I: Fundamental Theory", IEEE Trans. Circ.
  // Sys. for Video Tech., April 1997.
  if (var == 0) {
    *rate = 0;
    *dist = 0;
  } else {
    int d_q10, r_q10;
    static const uint32_t MAX_XSQ_Q10 = 245727;
    const uint64_t xsq_q10_64 =
        (((uint64_t)qstep * qstep << (n_log2 + 10)) + (var >> 1)) / var;
    const int xsq_q10 = (int)AOMMIN(xsq_q10_64, MAX_XSQ_Q10);
    model_rd_norm(xsq_q10, &r_q10, &d_q10);
    *rate = ROUND_POWER_OF_TWO(r_q10 << n_log2, 10 - AV1_PROB_COST_SHIFT);
    *dist = (var * (int64_t)d_q10 + 512) >> 10;
  }
}

static double interp_cubic(const double *p, double x) {
  return p[1] + 0.5 * x *
                    (p[2] - p[0] +
                     x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                          x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

static double interp_bicubic(const double *p, int p_stride, double x,
                             double y) {
  double q[4];
  q[0] = interp_cubic(p, x);
  q[1] = interp_cubic(p + p_stride, x);
  q[2] = interp_cubic(p + 2 * p_stride, x);
  q[3] = interp_cubic(p + 3 * p_stride, x);
  return interp_cubic(q, y);
}

static const uint8_t bsize_model_cat_lookup[BLOCK_SIZES_ALL] = {
  0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 0, 0, 1, 1, 2, 2
};

static const double interp_rgrid_surf[4][33 * 18] = {
  {
      22.327766,   25.923067,   22.621252,   24.194063,   33.896771,
      36.319449,   36.366516,   36.979779,   39.464116,   41.447575,
      46.889195,   48.242725,   48.263586,   48.247683,   47.671045,
      43.549538,   34.704173,   27.228206,   26.770352,   31.083589,
      27.132589,   29.003093,   40.611265,   43.538468,   43.604120,
      44.349390,   47.587462,   50.774400,   56.488014,   57.845609,
      57.850690,   57.274052,   53.144594,   43.870329,   39.652468,
      32.632756,   26.841850,   31.332654,   27.873044,   28.764785,
      38.785286,   43.170358,   43.830838,   44.969559,   48.914529,
      55.248689,   57.715924,   57.999751,   57.440703,   53.297908,
      44.019158,   39.898232,   39.511037,   33.825638,   26.852364,
      32.010800,   30.562706,   27.509767,   31.043052,   41.244780,
      44.290230,   46.940215,   49.690755,   56.343742,   57.979010,
      57.726394,   54.443900,   44.304913,   39.902663,   39.634229,
      41.771555,   38.753461,   26.854985,   32.180082,   31.245253,
      27.244444,   29.124383,   40.764816,   44.408952,   47.699710,
      50.900230,   56.625450,   57.723177,   55.589891,   48.900788,
      41.048654,   39.651996,   41.946139,   50.014160,   53.210943,
      26.855026,   32.194419,   32.011426,   30.277402,   29.877105,
      40.762578,   44.408839,   48.787255,   55.239820,   57.449038,
      55.607836,   49.190974,   42.194715,   39.937752,   41.950569,
      50.377665,   76.657003,   105.411775,  26.880862,   32.257610,
      34.445056,   40.463171,   34.053336,   40.413470,   43.244692,
      48.884515,   56.080787,   55.866764,   50.336939,   42.484901,
      39.959949,   41.955000,   50.377188,   77.031907,   131.719110,
      129.653489,  28.545717,   33.899273,   33.497577,   37.244333,
      40.609556,   40.872757,   40.799357,   50.888230,   54.551623,
      51.368178,   47.108051,   41.112398,   41.972866,   50.377188,
      76.996760,   129.827763,  144.790408,  121.709632,  35.139643,
      41.089108,   38.857853,   45.773848,   48.204882,   50.281504,
      49.790855,   61.085675,   52.973518,   48.627377,   44.175499,
      42.734512,   50.401937,   76.996975,   129.686810,  135.933948,
      100.408000,  76.304392,   37.278986,   46.411705,   51.480947,
      54.188399,   54.637112,   58.639801,   57.338896,   61.288437,
      60.309154,   54.247756,   52.843246,   55.810227,   77.781493,
      129.775529,  136.065397,  98.664431,   81.823654,   70.731319,
      43.343464,   56.934639,   63.041391,   65.880816,   67.307934,
      69.485021,   72.504386,   73.994384,   77.404783,   76.762250,
      72.919071,   87.604487,   124.655467,  141.139733,  109.506655,
      101.573425,  110.558962,  110.954184,  56.275468,   75.586969,
      79.134744,   85.760823,   91.432396,   96.692219,   101.593076,
      108.688433,  115.626170,  123.146328,  126.933619,  122.721080,
      120.884622,  141.208229,  150.944485,  190.008523,  215.816865,
      187.903946,  71.745710,   97.211914,   112.384864,  132.155864,
      150.936701,  165.780213,  174.134719,  180.493771,  189.913560,
      202.662326,  210.344006,  211.196275,  207.805823,  208.800962,
      226.054324,  239.983737,  245.673992,  205.828309,  124.479827,
      159.126024,  207.167074,  247.957043,  280.027846,  300.565341,
      309.125014,  308.844803,  308.242035,  320.398283,  326.903952,
      331.848231,  318.933315,  276.715654,  252.364942,  248.615769,
      248.016097,  206.870967,  231.034935,  314.086494,  381.796325,
      426.224473,  458.884493,  482.152708,  488.527863,  483.299350,
      476.274792,  476.459988,  478.876631,  469.963234,  429.787653,
      394.277143,  324.432401,  307.597458,  306.548302,  255.676476,
      413.968464,  554.882940,  605.627570,  637.923272,  657.421397,
      679.157516,  687.181830,  679.896017,  669.870485,  658.258846,
      646.474815,  625.548342,  565.207304,  545.127806,  561.608002,
      566.401720,  565.022957,  471.257149,  582.721167,  767.921972,
      824.949807,  849.221107,  863.486811,  878.828835,  890.366118,
      887.907263,  879.589391,  871.849383,  856.178689,  810.744328,
      743.088736,  742.075795,  754.486011,  757.009409,  755.131391,
      629.817369,  703.678130,  915.595218,  1015.693601, 1059.130673,
      1068.368861, 1064.177407, 1078.635964, 1091.553657, 1093.451802,
      1095.298009, 1080.827027, 1019.437255, 937.748773,  932.787441,
      946.364480,  956.007929,  955.406277,  796.879732,  823.683255,
      1039.054698, 1140.890662, 1242.449688, 1273.871014, 1255.349022,
      1259.115468, 1281.155098, 1294.597851, 1308.719628, 1308.385046,
      1248.994570, 1149.675281, 1114.567257, 1152.087224, 1190.457022,
      1195.269769, 997.101000,  942.355457,  1199.851925, 1281.731783,
      1384.658504, 1474.599214, 1479.461120, 1448.850602, 1463.795485,
      1482.123236, 1500.951075, 1522.013080, 1486.301163, 1375.733353,
      1326.076757, 1329.957795, 1346.628852, 1373.858211, 1151.886846,
      1036.706394, 1323.757789, 1414.185291, 1515.619007, 1617.286638,
      1650.141206, 1668.093106, 1653.480989, 1653.698132, 1692.389087,
      1716.607222, 1697.110553, 1607.303836, 1543.859439, 1533.502435,
      1565.481866, 1689.137857, 1435.374865, 1123.038914, 1421.142142,
      1532.687860, 1635.400944, 1693.543551, 1792.484818, 1818.433605,
      1832.818987, 1881.506293, 1886.105218, 1901.578698, 1879.745193,
      1797.544294, 1746.289310, 1755.344122, 1779.963154, 1855.519162,
      1564.239506, 1199.400912, 1524.660922, 1646.416078, 1734.737227,
      1808.709132, 1846.503686, 1909.385219, 2066.433073, 2134.005568,
      2107.387402, 2091.425244, 2053.032173, 1935.803667, 1822.906046,
      1882.648932, 1915.949512, 1928.985812, 1612.128526, 1284.024933,
      1620.824219, 1731.806554, 1820.161672, 1855.255123, 1916.872117,
      2118.561996, 2313.234943, 2398.766982, 2400.007640, 2302.576938,
      2236.162502, 2108.546124, 1966.471643, 2012.027085, 2065.374469,
      2078.087903, 1735.287563, 1365.312552, 1714.219671, 1816.795361,
      1855.737333, 1918.451238, 2120.106526, 2324.598160, 2473.530807,
      2638.825497, 2557.986746, 2509.463813, 2454.667835, 2316.706304,
      2140.184009, 2128.835809, 2143.879492, 2142.854643, 1787.750158,
      1437.419836, 1779.555341, 1847.489680, 1918.325091, 2120.113693,
      2324.617694, 2473.846782, 2649.468965, 2620.410085, 2648.675582,
      2699.137109, 2675.278121, 2523.476319, 2297.248430, 2253.688417,
      2254.913300, 2249.223343, 1875.971541, 1493.896696, 1834.866071,
      1916.269602, 2119.813096, 2323.539728, 2473.577989, 2649.467145,
      2620.622628, 2653.420533, 2728.513109, 2781.855529, 2834.544067,
      2736.568264, 2371.324644, 2282.442721, 2281.314152, 2275.463779,
      1897.849564, 1539.549127, 1911.564001, 2119.764580, 2322.461269,
      2469.254923, 2648.389179, 2620.605916, 2653.421457, 2728.592640,
      2783.540151, 2851.509097, 2878.698355, 2782.517962, 2382.893912,
      2283.266840, 2281.725921, 2275.870459, 1898.188755, 1606.706318,
      2115.038111, 2322.457102, 2468.986130, 2647.311213, 2620.337122,
      2653.417289, 2728.592640, 2783.540431, 2851.529588, 2880.473679,
      2882.960994, 2783.531114, 2384.093347, 2283.529255, 2281.729952,
      2275.870459, 1898.188755, 1803.536504, 2317.627627, 2468.986065,
      2647.307046, 2620.320410, 2653.413122, 2728.592575, 2783.540431,
      2851.529588, 2880.473679, 2884.504202, 2883.328227, 2800.050192,
      2450.337864, 2300.047472, 2281.986048, 2275.870459, 1898.188755,
      1966.177469, 2463.537233, 2647.307046, 2620.320410, 2653.413122,
      2728.592575, 2783.540431, 2851.529588, 2880.473679, 2884.504202,
      2884.615338, 2884.356134, 2866.294709, 2716.004317, 2366.291989,
      2283.013093, 2275.870459, 1898.188755, 2070.297745, 2634.518769,
      2613.296256, 2646.290355, 2721.246305, 2776.138593, 2844.126888,
      2873.070979, 2877.101502, 2877.212638, 2877.213500, 2877.209529,
      2875.414196, 2775.106201, 2376.694617, 2277.409697, 2270.029964,
      1893.317493, 1882.250271, 2170.027657, 2193.231808, 2261.827187,
      2302.303482, 2366.707610, 2395.596133, 2399.626655, 2399.737792,
      2399.738653, 2399.738653, 2399.738653, 2398.451542, 2315.432712,
      1982.494696, 1899.475867, 1893.317493, 1579.120622,
  },
  {
      6.262675,    7.439091,    7.698993,    7.778980,    7.875715,
      8.108790,    8.147575,    8.082632,    8.066267,    8.066015,
      8.066246,    8.074624,    8.134507,    8.261873,    8.294361,
      8.316923,    8.486638,    7.554566,    7.508544,    8.905528,
      9.176430,    9.350165,    9.591205,    9.747009,    9.719927,
      9.678525,    9.671020,    9.671137,    9.679515,    9.739514,
      9.874331,    9.936702,    9.966750,    10.171737,   11.338029,
      11.467385,   7.527119,    8.873752,    8.979175,    9.328706,
      9.652844,    9.732245,    9.546719,    9.653964,    9.695250,
      9.703934,    9.764280,    9.899213,    9.961701,    9.992213,
      10.197316,   11.407407,   14.315840,   13.772040,   7.534136,
      8.888990,    8.931233,    9.315212,    9.653045,    9.720054,
      9.497227,    9.641624,    9.695638,    9.734397,    9.891761,
      9.961586,    9.992213,    10.197316,   11.407407,   14.376951,
      17.713669,   18.671825,   7.563129,    9.005062,    8.959378,
      9.315468,    9.653044,    9.719863,    9.496461,    9.641664,
      9.703666,    9.771917,    9.931703,    9.991749,    10.197316,
      11.407407,   14.376951,   17.794155,   23.863171,   24.804787,
      7.575519,    9.039883,    8.971080,    9.320311,    9.658319,
      9.723488,    9.502759,    9.651691,    9.763223,    9.901695,
      9.984296,    10.197200,   11.407407,   14.376951,   17.794155,
      24.004415,   33.915049,   41.788173,   7.911491,    9.418600,
      9.243672,    9.552386,    9.983055,    9.961797,    9.904549,
      9.824931,    9.881775,    9.979604,    10.197659,   11.407416,
      14.376951,   17.794155,   24.004415,   34.071492,   51.878735,
      49.657853,   9.447816,    11.112602,   10.543783,   10.725970,
      11.517217,   11.215939,   11.680662,   10.843212,   10.313657,
      10.421884,   11.458569,   14.378861,   17.794173,   24.004480,
      34.075659,   52.056293,   60.036524,   51.837690,   10.927916,
      12.683731,   12.328553,   12.655561,   13.307360,   13.322625,
      13.537403,   13.711979,   13.305593,   13.204373,   15.006696,
      17.897589,   24.010120,   34.079895,   52.325100,   61.292102,
      63.559020,   56.464217,   13.381673,   17.290200,   16.284541,
      16.077053,   16.277702,   17.126057,   18.184339,   19.539514,
      20.935187,   21.151970,   21.564043,   25.769198,   34.394665,
      52.350783,   62.391113,   68.167594,   69.387417,   58.740769,
      17.569465,   29.121669,   22.254776,   20.293023,   21.582658,
      24.040717,   27.336878,   31.707605,   37.148429,   40.852285,
      40.180300,   43.475897,   55.261472,   62.907110,   69.908542,
      77.453519,   78.876581,   65.820970,   20.688111,   35.464135,
      28.590042,   29.859684,   37.106004,   44.830617,   51.510355,
      61.454756,   72.885129,   82.427833,   85.145105,   85.301898,
      84.775608,   83.937031,   90.524072,   114.021545,  119.511209,
      99.753183,   30.938439,   47.989519,   54.896116,   72.132343,
      93.728428,   111.210180,  119.201301,  126.976908,  141.443652,
      153.999978,  162.436182,  164.552551,  155.760626,  145.606839,
      146.674969,  159.408450,  162.282038,  135.153611,  66.474280,
      87.571887,   138.334800,  186.410316,  219.663563,  240.006854,
      241.793309,  235.259026,  238.160006,  251.168220,  257.448854,
      249.457596,  215.322054,  188.238161,  183.699402,  188.854501,
      196.691767,  148.687823,  133.133673,  199.547259,  286.243461,
      347.003147,  385.366447,  406.557124,  398.244822,  376.377281,
      368.917319,  373.523444,  379.457093,  360.137476,  305.217051,
      276.392793,  272.716800,  287.679613,  319.261623,  206.033537,
      346.181038,  460.387661,  509.533364,  533.866444,  560.116896,
      583.225963,  573.575178,  541.342602,  528.416398,  530.900691,
      527.802607,  498.829026,  430.766443,  415.867812,  417.852673,
      428.528796,  438.923560,  356.633328,  527.219175,  691.563636,
      753.092236,  744.517322,  738.701940,  756.684249,  752.102200,
      723.385479,  710.808139,  715.868264,  703.995966,  655.340555,
      567.538643,  537.616534,  535.209423,  537.063320,  536.720334,
      448.894212,  655.965217,  846.209812,  948.937510,  978.175024,
      930.312738,  914.627044,  918.590256,  903.713896,  901.062558,
      915.566498,  906.979046,  838.104351,  730.315122,  708.505523,
      708.794808,  708.833589,  707.028700,  589.719945,  782.378840,
      1017.745398, 1120.303734, 1184.005590, 1155.480331, 1098.853178,
      1081.056778, 1077.002270, 1083.125009, 1111.096783, 1121.676326,
      1053.013849, 926.412140,  897.880661,  897.578485,  897.689993,
      895.624943,  747.042658,  885.551729,  1156.520033, 1232.700221,
      1310.599983, 1364.300269, 1321.579159, 1274.030969, 1256.716378,
      1268.341472, 1291.131471, 1325.417154, 1281.277473, 1128.270819,
      1034.205749, 1019.431365, 1026.301994, 1039.242827, 869.853512,
      993.073680,  1262.530901, 1324.307314, 1409.447946, 1479.963846,
      1520.862525, 1513.541034, 1452.262977, 1454.056131, 1475.674490,
      1505.504388, 1484.838325, 1346.369054, 1291.400107, 1289.232319,
      1320.381767, 1390.386950, 1174.251661, 1074.639996, 1347.053773,
      1418.178262, 1480.851536, 1569.820521, 1687.670043, 1689.515476,
      1648.902063, 1653.757784, 1672.791230, 1689.212374, 1654.078557,
      1530.523723, 1477.003827, 1471.696416, 1490.511936, 1548.859722,
      1304.586825, 1136.193865, 1416.762638, 1476.736548, 1545.957291,
      1679.037014, 1730.853678, 1766.921891, 1849.258564, 1878.365328,
      1854.066857, 1880.388655, 1833.035858, 1690.940430, 1597.114513,
      1596.717317, 1609.516368, 1618.881033, 1352.718943, 1189.974265,
      1460.621832, 1541.313734, 1663.771250, 1690.450024, 1763.090588,
      1896.056429, 2005.497847, 1962.300376, 1973.891926, 2098.063701,
      2057.407525, 1873.225589, 1702.411123, 1734.875462, 1774.876151,
      1776.803671, 1482.058363, 1224.117022, 1534.445462, 1662.921738,
      1687.914916, 1753.320046, 1895.613175, 2017.923275, 1973.823791,
      1968.099428, 2269.375830, 2388.604831, 2304.438904, 2087.104224,
      1861.935304, 1843.461426, 1827.084051, 1816.350155, 1514.845817,
      1280.363797, 1655.992005, 1687.853454, 1753.282051, 1895.539646,
      2022.923236, 1994.041530, 1970.321478, 2285.332836, 2541.497572,
      2552.864045, 2494.010394, 2312.711822, 2064.081024, 2015.226628,
      1899.704219, 1864.017207, 1554.281285, 1357.391307, 1671.427169,
      1753.087029, 1895.539646, 2022.942670, 1995.294998, 1975.347503,
      2286.620309, 2547.557386, 2605.473170, 2673.059800, 2681.030506,
      2564.159534, 2186.218206, 2094.120966, 2063.864325, 2050.875846,
      1710.431918, 1390.553019, 1745.790833, 1895.491017, 2022.942670,
      1995.295300, 1975.366936, 2286.698245, 2547.578195, 2605.653144,
      2679.462677, 2717.948438, 2721.774852, 2619.188434, 2208.598571,
      2105.257885, 2102.961379, 2097.441280, 1749.368423, 1471.767838,
      1891.196999, 2022.941916, 1995.295300, 1975.366936, 2286.698245,
      2547.578195, 2605.653144, 2679.464068, 2718.040940, 2723.901908,
      2722.741253, 2637.030354, 2277.099964, 2122.418316, 2103.829164,
      2098.163195, 1749.972089, 1603.711220, 2018.484778, 1995.295300,
      1975.366936, 2286.698245, 2547.578195, 2605.653144, 2679.464068,
      2718.040940, 2723.901908, 2724.065973, 2723.799587, 2705.210435,
      2550.528787, 2190.598397, 2104.886218, 2098.163195, 1749.972089,
      1704.935985, 1990.434085, 1975.366936, 2286.698245, 2547.578195,
      2605.653144, 2679.464068, 2718.040940, 2723.901908, 2724.065973,
      2724.067252, 2724.063166, 2722.215375, 2618.972446, 2208.656304,
      2105.413375, 2098.167281, 1749.972089, 1644.139660, 1963.419473,
      2280.060667, 2540.670333, 2598.663735, 2672.473379, 2711.050251,
      2716.911219, 2717.075284, 2717.076564, 2717.076564, 2717.076564,
      2715.751844, 2630.043752, 2270.376941, 2116.752347, 2093.042323,
      1745.481191, 1297.803162, 1845.422979, 2112.545912, 2153.105955,
      2221.655846, 2260.150200, 2266.011168, 2266.175233, 2266.176512,
      2266.176512, 2266.176512, 2266.176512, 2265.908847, 2247.587360,
      2110.170145, 1819.476889, 1746.538245, 1455.817819,
  },
  {
      1.613984,    1.935090,    1.938435,    1.938782,    1.966672,
      1.982255,    2.015457,    2.022916,    1.979904,    1.808882,
      1.766270,    1.766077,    1.791534,    1.893256,    1.948358,
      2.165754,    2.715814,    2.665682,    1.936673,    2.320401,
      2.317574,    2.298246,    2.351479,    2.376785,    2.416396,
      2.416937,    2.339842,    2.160318,    2.117728,    2.122793,
      2.168957,    2.304755,    2.547391,    3.139157,    4.127166,
      4.807955,    2.041962,    2.351356,    2.322225,    2.297399,
      2.359420,    2.397272,    2.425675,    2.389126,    2.209170,
      2.131947,    2.128295,    2.169387,    2.290159,    2.548211,
      3.145001,    4.156470,    6.297385,    6.995383,    2.444633,
      2.453221,    2.324141,    2.297521,    2.373793,    2.455017,
      2.439811,    2.380317,    2.175087,    2.129088,    2.169388,
      2.290081,    2.547899,    3.144923,    4.156469,    6.342137,
      9.881919,    13.075095,   2.569957,    2.578604,    2.349816,
      2.298176,    2.377466,    2.469247,    2.434899,    2.346294,
      2.171702,    2.170049,    2.290081,    2.547899,    3.144923,
      4.156469,    6.342137,    9.944302,    17.098819,   19.046234,
      2.670532,    2.984006,    2.468490,    2.314588,    2.384346,
      2.468030,    2.403018,    2.216778,    2.179424,    2.290593,
      2.549042,    3.145208,    4.156473,    6.342137,    9.944302,
      17.210443,   26.245986,   32.942300,   2.609600,    3.025575,
      2.607787,    2.451351,    2.502184,    2.559416,    2.520158,
      2.316537,    2.333242,    2.573958,    3.218940,    4.174834,
      6.342422,    9.944302,    17.210443,   26.369867,   40.932662,
      39.306613,   2.323681,    2.815857,    2.915570,    2.898623,
      3.013214,    3.137214,    3.259217,    3.097019,    2.942892,
      3.358477,    4.475945,    6.416996,    9.945458,    17.210443,
      26.369867,   41.056733,   47.309232,   40.195621,   2.645339,
      3.408980,    3.642128,    3.604415,    3.771696,    4.212864,
      4.794327,    5.212609,    5.086837,    5.212025,    6.833239,
      10.036015,   17.213840,   26.369901,   41.056733,   47.433304,
      48.198240,   40.220136,   3.446585,    4.969326,    4.713654,
      4.685640,    4.963227,    5.978995,    7.550985,    9.367771,
      10.908408,   11.665776,   12.789162,   18.243908,   26.531265,
      41.109704,   47.494165,   48.385813,   48.316701,   40.371652,
      3.943882,    6.034464,    5.889472,    6.511478,    8.258414,
      10.939460,   14.668837,   19.425375,   25.122294,   29.744153,
      30.257586,   32.420852,   42.854069,   50.921156,   52.457060,
      52.651998,   54.491207,   50.154942,   4.677004,    7.475672,
      9.339550,    13.670512,   21.624305,   29.993455,   38.218236,
      48.705632,   58.485547,   66.829161,   68.644891,   65.435354,
      63.345149,   74.861396,   78.016941,   78.708481,   86.423562,
      90.956742,   8.358829,    13.593058,   28.694070,   49.715243,
      73.129397,   91.059591,   99.427656,   106.033834,  116.932863,
      125.909619,  127.521092,  123.755674,  118.591601,  131.294759,
      134.956752,  135.235986,  137.110126,  119.114603,  34.423222,
      51.827467,   100.886563,  154.837905,  191.573446,  212.546870,
      208.408590,  195.047709,  194.729770,  204.466358,  206.391060,
      201.457385,  192.776688,  204.238285,  207.844866,  211.889264,
      227.311130,  192.976194,  152.530527,  198.621256,  240.702674,
      299.221661,  340.783457,  363.243261,  344.006632,  308.941723,
      296.581930,  302.646592,  309.171713,  296.898229,  272.669291,
      264.795441,  264.339193,  280.316000,  343.502316,  299.790306,
      318.584501,  438.526552,  454.052901,  463.644859,  488.971114,
      517.737455,  496.088354,  448.188527,  428.548670,  428.809858,
      432.465332,  414.957102,  364.217175,  367.050948,  370.940201,
      374.988997,  389.961698,  328.561859,  476.626713,  643.786670,
      688.485027,  666.142956,  636.563325,  658.830965,  647.939850,
      602.716530,  581.553952,  583.966331,  580.631307,  551.434654,
      472.543100,  447.089231,  445.224840,  445.291548,  444.396349,
      370.699966,  599.819910,  789.624848,  886.777881,  896.633336,
      816.615786,  786.364001,  788.117725,  761.855575,  743.659120,
      753.405047,  754.404395,  706.050427,  616.949550,  589.495474,
      594.893030,  596.873475,  595.305782,  496.297817,  700.277519,
      944.489816,  1058.699894, 1119.768103, 1055.452015, 972.587595,
      938.636840,  917.891447,  907.316821,  926.846721,  936.676467,
      877.375543,  772.656560,  742.876277,  771.247054,  779.285108,
      773.353888,  631.095162,  819.341176,  1085.985219, 1174.247742,
      1230.190785, 1251.005393, 1209.813272, 1126.913150, 1086.784971,
      1078.922485, 1095.538470, 1115.182761, 1063.999674, 947.196025,
      924.920170,  933.489013,  939.765148,  938.071469,  730.098872,
      928.338108,  1175.811526, 1244.373712, 1303.863309, 1392.260318,
      1419.063789, 1361.830192, 1280.983792, 1261.447278, 1271.249636,
      1280.938671, 1242.741711, 1106.747098, 1052.249862, 1046.607578,
      1068.568361, 1149.777779, 963.356904,  993.790193,  1246.393259,
      1326.237603, 1405.484043, 1475.560884, 1554.122846, 1575.403227,
      1531.489225, 1496.096686, 1481.123286, 1467.647936, 1401.019577,
      1266.139501, 1203.276643, 1194.879698, 1216.428241, 1298.879373,
      1100.912917, 1051.668893, 1327.490245, 1412.938704, 1482.720032,
      1564.735644, 1643.243081, 1713.805000, 1735.203701, 1768.893909,
      1723.894482, 1700.136055, 1594.509861, 1438.390775, 1363.960477,
      1357.361978, 1379.089851, 1397.097230, 1168.877421, 1123.985063,
      1421.579651, 1486.097340, 1566.198606, 1643.590059, 1720.037210,
      1796.360637, 1964.066319, 1975.726755, 1974.558540, 1948.645959,
      1873.432254, 1669.073524, 1536.603783, 1535.978858, 1607.415138,
      1621.100690, 1352.360079, 1202.245810, 1485.634132, 1566.938715,
      1643.571242, 1716.697428, 1783.151348, 1969.343107, 2038.819154,
      2137.625726, 2138.166832, 2163.338914, 2138.134532, 1904.821576,
      1705.863604, 1662.457761, 1676.837189, 1676.864315, 1398.644613,
      1248.692000, 1563.460473, 1643.584711, 1716.658988, 1782.212368,
      1966.018111, 2038.420777, 2144.567626, 2174.128911, 2239.527266,
      2336.322349, 2321.786112, 2153.538702, 1928.371089, 1762.802130,
      1731.371882, 1726.515947, 1440.000745, 1320.496114, 1640.418472,
      1716.643278, 1780.558394, 1959.428411, 2040.298220, 2158.900759,
      2177.829635, 2241.771056, 2367.204227, 2476.770268, 2489.425171,
      2385.814787, 2049.280607, 1936.352503, 1927.037919, 1921.967879,
      1603.016463, 1380.620258, 1712.765176, 1780.458044, 1952.796134,
      2013.922986, 2166.638751, 2235.357216, 2256.141651, 2367.527429,
      2482.141240, 2522.163364, 2525.441439, 2434.356482, 2069.785382,
      1977.547028, 1975.779927, 1970.703987, 1643.664793, 1441.019348,
      1776.247259, 1952.770533, 2012.269213, 2160.062029, 2237.286701,
      2270.486281, 2371.110687, 2482.198011, 2522.244370, 2527.330342,
      2526.291935, 2450.190500, 2130.613969, 1993.274567, 1976.769552,
      1971.459584, 1644.294999, 1489.931103, 1948.331529, 2012.268815,
      2160.036389, 2237.184737, 2270.516195, 2371.333084, 2482.253565,
      2522.245231, 2527.330342, 2527.468132, 2527.231537, 2510.726521,
      2373.387119, 2053.810589, 1977.708094, 1971.459584, 1644.294999,
      1660.688668, 2007.809242, 2160.036389, 2237.184737, 2270.516195,
      2371.333084, 2482.253565, 2522.245231, 2527.330342, 2527.468132,
      2527.469193, 2527.465564, 2525.824940, 2434.157168, 2069.843921,
      1978.176148, 1971.463213, 1644.294999, 1671.699856, 2149.137365,
      2231.105063, 2264.100466, 2364.847980, 2475.767400, 2515.759066,
      2520.844177, 2520.981967, 2520.983028, 2520.983028, 2520.983028,
      2519.806830, 2443.707964, 2124.365460, 1987.964600, 1966.634317,
      1640.075296, 1514.149810, 1854.381627, 1871.960909, 1951.032837,
      2057.477572, 2097.400823, 2102.485934, 2102.623724, 2102.624785,
      2102.624785, 2102.624785, 2102.624785, 2102.387129, 2086.119769,
      1964.109185, 1706.007217, 1641.013838, 1367.904079,
  },
  {
      0.189299,    0.229228,    0.238048,    0.268887,    0.281142,
      0.300935,    0.309796,    0.310918,    0.310937,    0.311146,
      0.312713,    0.316330,    0.328262,    0.409035,    0.618740,
      0.949941,    1.653170,    1.784429,    0.218491,    0.272725,
      0.285358,    0.321015,    0.331577,    0.359439,    0.371414,
      0.372781,    0.372819,    0.373806,    0.378200,    0.391098,
      0.472222,    0.691840,    1.063346,    1.823813,    2.823813,
      3.615378,    0.216849,    0.272844,    0.285990,    0.316324,
      0.310304,    0.354847,    0.372284,    0.373746,    0.373989,
      0.376045,    0.391284,    0.473184,    0.692820,    1.064479,
      1.825572,    2.844380,    4.787408,    5.485676,    0.216260,
      0.270611,    0.285414,    0.314938,    0.304784,    0.353470,
      0.372266,    0.373954,    0.375850,    0.390606,    0.473384,
      0.692915,    1.064481,    1.825572,    2.844380,    4.823394,
      7.806808,    10.489172,   0.214035,    0.261700,    0.283234,
      0.314927,    0.304827,    0.353579,    0.372421,    0.375084,
      0.390508,    0.479709,    0.718336,    1.070820,    1.825670,
      2.844380,    4.823394,    7.856992,    13.726022,   15.320282,
      0.213703,    0.260423,    0.285850,    0.318083,    0.313461,
      0.362455,    0.382302,    0.389398,    0.479353,    0.744049,
      1.173988,    1.851395,    2.844778,    4.823394,    7.856992,
      13.815504,   21.091861,   26.419131,   0.227492,    0.278836,
      0.311462,    0.348960,    0.367754,    0.427027,    0.461678,
      0.494050,    0.745559,    1.200067,    1.954674,    2.870530,
      4.823794,    7.856992,    13.815504,   21.190814,   32.801624,
      31.418061,   0.290771,    0.353917,    0.389185,    0.451632,
      0.495950,    0.621270,    0.749745,    0.914889,    1.324104,
      2.064863,    3.188577,    4.902607,    7.858213,    13.815515,
      21.191535,   32.903044,   37.773998,   31.958897,   0.350426,
      0.434934,    0.512010,    0.597215,    0.711539,    1.024520,
      1.585378,    2.436157,    3.408353,    3.996570,    5.333157,
      7.946050,    13.818735,   21.190997,   32.949186,   38.061432,
      38.526094,   32.605583,   0.437661,    0.517781,    0.704742,
      0.855991,    1.122258,    1.995200,    3.979176,    6.542714,
      9.494608,    10.456686,   10.736623,   14.888952,   21.334225,
      32.871045,   38.230390,   39.391925,   39.379174,   32.974928,
      0.667323,    0.636232,    1.268594,    1.746815,    2.937984,
      5.779973,    10.942591,   15.872349,   22.428396,   25.699069,
      25.684790,   27.119510,   33.924927,   38.054713,   39.703480,
      40.888018,   41.028675,   34.225022,   1.530388,    1.313254,
      3.527730,    6.089944,    11.805656,   20.408891,   30.911018,
      39.783315,   48.102426,   53.401737,   55.090997,   52.915124,
      48.954513,   48.732610,   49.927551,   53.481688,   54.187545,
      45.199686,   3.087916,    4.497456,    15.982772,   29.034626,
      48.215660,   66.344749,   75.409963,   82.193674,   89.037887,
      93.713379,   96.467289,   94.900069,   91.861402,   88.584111,
      87.067073,   87.664751,   87.533539,   72.983182,   2.872201,
      18.556404,   64.090577,   102.775677,  137.958868,  158.042073,
      149.718018,  140.527287,  139.918791,  146.041339,  150.632261,
      149.457331,  134.255846,  118.984381,  118.784779,  121.985029,
      122.311521,  102.015710,  28.248313,   82.866817,   156.997146,
      210.276224,  244.575300,  258.524081,  232.667923,  208.658317,
      201.365405,  205.904842,  213.280236,  209.802029,  189.592969,
      182.386667,  210.577536,  228.169569,  230.336431,  192.146182,
      146.572031,  274.657225,  316.176128,  321.869827,  323.467609,
      330.950486,  308.073551,  283.468692,  268.826362,  267.880448,
      274.132681,  266.193452,  247.853062,  248.027919,  253.809445,
      257.546064,  257.561163,  214.827324,  309.277171,  480.944806,
      524.048101,  464.259230,  388.281957,  375.099540,  371.054799,
      356.566115,  341.321587,  334.827793,  337.655829,  325.593464,
      302.627911,  292.344854,  279.803428,  276.864255,  276.118017,
      230.296243,  432.046593,  592.338624,  585.768072,  548.784487,
      459.941646,  413.142484,  416.871830,  424.795092,  418.953520,
      407.383657,  403.849069,  388.320341,  366.732043,  363.308134,
      360.767312,  359.978099,  359.016964,  299.394311,  521.324465,
      626.755756,  614.329239,  621.435226,  569.128464,  488.002514,
      473.099634,  481.085277,  489.656962,  483.918994,  477.495606,
      454.853845,  432.036758,  430.932862,  430.966415,  427.792608,
      425.131640,  351.787315,  518.127342,  628.818563,  686.411730,
      687.346046,  631.427027,  572.715518,  540.215653,  544.293907,
      557.302840,  561.205715,  562.357584,  535.124712,  509.451692,
      498.328186,  493.473768,  482.205757,  480.687889,  390.959754,
      529.127869,  751.253047,  881.133714,  701.807620,  642.631866,
      667.156870,  596.489950,  615.101817,  628.508517,  640.930921,
      640.551065,  618.437963,  581.162398,  570.006951,  568.475806,
      571.668393,  594.144327,  498.057496,  652.162883,  924.023974,
      855.250607,  665.669781,  669.076324,  612.926757,  679.040560,
      718.892048,  702.360171,  691.865382,  709.136901,  702.130142,
      659.699006,  651.662187,  652.186398,  653.836117,  658.623386,
      550.627649,  802.998685,  860.698507,  678.293710,  607.121030,
      594.591940,  686.681634,  786.710365,  816.566385,  872.590613,
      781.402241,  791.410541,  807.089387,  778.069658,  799.987135,
      809.079659,  815.183031,  820.047297,  685.188807,  711.459934,
      679.282151,  603.313719,  577.397989,  682.387869,  789.593080,
      844.025092,  1016.505215, 1135.590748, 949.725207,  935.359433,
      959.265515,  931.846999,  934.428712,  947.128474,  972.392443,
      997.418130,  836.741232,  585.699187,  615.165071,  577.551008,
      682.119545,  789.526840,  844.117007,  1020.532519, 1189.979456,
      1138.517203, 1091.417593, 1150.972523, 1119.151224, 1025.484936,
      1014.404351, 1039.592996, 1051.302268, 1055.550286, 881.589383,
      576.316845,  596.076655,  682.543736,  789.959965,  844.225008,
      1020.534959, 1190.124723, 1145.130845, 1132.232575, 1290.602055,
      1481.914983, 1335.022841, 1156.455266, 1117.559441, 1121.429355,
      1122.935564, 1120.182318, 934.306300,  747.712630,  749.295605,
      798.128192,  872.163909,  1027.501009, 1190.232723, 1145.132268,
      1132.368163, 1294.536558, 1542.091529, 1600.666915, 1508.348816,
      1262.168946, 1153.030876, 1138.920317, 1138.735763, 1135.813805,
      947.325005,  669.156553,  851.621796,  909.983290,  1139.686823,
      1218.169312, 1145.565392, 1132.368163, 1294.537224, 1542.263000,
      1610.114508, 1606.261139, 1589.809935, 1474.120526, 1204.263721,
      1139.972407, 1138.978724, 1136.055799, 947.526839,  760.844456,
      1063.059416, 1185.238619, 1246.704066, 1152.531442, 1132.476164,
      1294.537224, 1542.263000, 1610.116443, 1606.586635, 1604.444633,
      1602.988695, 1525.357641, 1217.010752, 1140.170005, 1138.978724,
      1136.055799, 947.526839,  928.443510,  1221.585752, 1256.432493,
      1153.113721, 1132.584165, 1294.538899, 1542.263000, 1610.116443,
      1606.586635, 1604.447698, 1604.380734, 1603.188910, 1526.348163,
      1218.001274, 1140.370763, 1138.981789, 1136.055799, 947.526839,
      1042.376643, 1254.249913, 1153.262875, 1132.586477, 1294.538899,
      1542.263000, 1610.116443, 1606.586635, 1604.447698, 1604.380734,
      1604.380190, 1603.386604, 1539.099410, 1269.138768, 1153.122010,
      1139.179483, 1136.055799, 947.526839,  1046.078187, 1146.822724,
      1128.497689, 1290.389157, 1538.145186, 1605.999173, 1602.469365,
      1600.330428, 1600.263463, 1600.262920, 1600.262920, 1600.062161,
      1586.122699, 1470.303634, 1201.135821, 1137.046321, 1133.140375,
      945.095231,  782.287926,  915.627388,  1026.661938, 1270.486433,
      1340.399801, 1336.905029, 1334.766092, 1334.699127, 1334.698584,
      1334.698584, 1334.698584, 1334.695519, 1333.507304, 1268.825734,
      1012.409168, 948.517361,  945.095231,  788.256261,
  },
};

static const double interp_dgrid_surf[33 * 18] = {
  10.732740, 12.852623, 12.871726, 12.831927, 12.830211, 12.861430, 12.867990,
  12.867782, 12.867835, 12.871858, 12.902024, 12.965509, 12.968408, 12.881384,
  12.675279, 12.474284, 12.072515, 9.889579,  12.929247, 15.425118, 15.433135,
  15.390981, 15.406660, 15.426355, 15.428419, 15.428079, 15.428421, 15.447409,
  15.525862, 15.547404, 15.464003, 15.247174, 15.002634, 14.598756, 14.118559,
  11.268177, 12.978390, 15.468997, 15.473333, 15.454267, 15.541287, 15.489730,
  15.468482, 15.467889, 15.472135, 15.505528, 15.572188, 15.503756, 15.287214,
  15.042508, 14.637955, 14.189846, 13.405506, 10.701531, 12.982633, 15.485081,
  15.477421, 15.460239, 15.564988, 15.495640, 15.468632, 15.471851, 15.501737,
  15.557136, 15.499691, 15.287073, 15.042507, 14.637955, 14.189846, 13.433340,
  12.496832, 9.340326,  12.998657, 15.549387, 15.493687, 15.460800, 15.565345,
  15.495777, 15.468954, 15.486779, 15.553263, 15.494038, 15.265509, 15.037143,
  14.637872, 14.189846, 13.433340, 12.521317, 10.919583, 8.161390,  13.003041,
  15.569001, 15.513672, 15.475296, 15.564665, 15.498759, 15.474738, 15.493960,
  15.478845, 15.243554, 14.950019, 14.616148, 14.189509, 13.433340, 12.521317,
  10.935392, 9.181084,  5.696463,  13.025275, 15.607663, 15.597233, 15.542220,
  15.563239, 15.514427, 15.497799, 15.490465, 15.246155, 14.930528, 14.538075,
  14.170042, 13.433039, 12.521317, 10.935392, 9.194731,  6.576716,  4.574338,
  13.113925, 15.713455, 15.693931, 15.596521, 15.565416, 15.507961, 15.468824,
  15.340517, 14.936491, 14.475185, 13.941369, 13.376111, 12.520434, 10.935388,
  9.194521,  6.589278,  5.438675,  4.363233,  13.137197, 15.737061, 15.707853,
  15.616963, 15.562968, 15.400673, 15.227099, 14.857736, 14.271325, 13.781258,
  13.226211, 12.485885, 10.934308, 9.194299,  6.575550,  5.396014,  5.158836,
  4.126388,  13.104712, 15.634910, 15.679537, 15.625973, 15.556098, 15.283626,
  14.694503, 13.972933, 13.092234, 12.473881, 11.999564, 10.706415, 9.155210,
  6.573035,  5.329847,  4.895610,  4.799558,  3.959161,  12.958260, 15.206588,
  15.513489, 15.505583, 15.269032, 14.750269, 13.877349, 12.965270, 11.607588,
  10.603269, 10.036216, 8.703476,  6.416604,  5.268632,  4.822348,  4.524073,
  4.451947,  3.711706,  12.807672, 14.961824, 15.187909, 14.953962, 14.257449,
  13.329376, 12.073392, 10.697850, 9.364881,  8.346431,  7.905747,  7.072510,
  5.280775,  4.629840,  4.381936,  4.289014,  4.265648,  3.557610,  12.361927,
  14.462730, 14.064972, 13.188092, 11.798536, 10.490230, 9.102954,  7.699975,
  6.613675,  5.889225,  5.627229,  5.268535,  4.695451,  4.394366,  3.869183,
  3.740339,  3.726922,  3.108088,  11.201816, 13.059030, 11.064286, 9.438105,
  7.955206,  6.858176,  5.846328,  4.871307,  4.234502,  3.869846,  3.714787,
  3.583539,  3.381026,  3.276404,  3.140935,  3.082152,  2.966113,  2.451496,
  9.462887,  10.165517, 7.319754,  5.729278,  4.721657,  4.033048,  3.468795,
  2.936523,  2.592624,  2.426184,  2.335125,  2.274122,  2.185747,  2.044925,
  2.009657,  1.900832,  1.464277,  1.131493,  5.769254,  5.782902,  4.275958,
  3.318562,  2.748567,  2.319186,  1.998873,  1.729693,  1.539165,  1.438280,
  1.394711,  1.343040,  1.238374,  1.120804,  1.131274,  1.113517,  1.003169,
  0.814304,  3.405195,  3.407397,  2.558915,  1.968935,  1.599929,  1.347929,
  1.145740,  0.997683,  0.896194,  0.830237,  0.802812,  0.778310,  0.711895,
  0.691497,  0.834230,  0.872105,  0.868793,  0.724269,  2.217809,  2.136880,
  1.546442,  1.187157,  0.965022,  0.798981,  0.668734,  0.570746,  0.512161,
  0.474407,  0.450901,  0.434308,  0.414115,  0.391743,  0.422912,  0.432405,
  0.431444,  0.359846,  1.417932,  1.351022,  1.004332,  0.743108,  0.588869,
  0.480837,  0.397131,  0.333965,  0.291170,  0.267766,  0.252477,  0.239263,
  0.227263,  0.201775,  0.199561,  0.200498,  0.199938,  0.166745,  0.901960,
  0.868050,  0.686228,  0.503925,  0.372170,  0.295718,  0.241543,  0.201026,
  0.170280,  0.150486,  0.140379,  0.132175,  0.125480,  0.119661,  0.131991,
  0.134575,  0.130280,  0.107826,  0.641236,  0.620162,  0.482665,  0.354209,
  0.255574,  0.196477,  0.153945,  0.123824,  0.102894,  0.087692,  0.077917,
  0.072577,  0.069383,  0.069821,  0.073402,  0.069871,  0.051978,  0.039665,
  0.480116,  0.455004,  0.329705,  0.246672,  0.187615,  0.135982,  0.104112,
  0.080926,  0.063595,  0.053878,  0.046111,  0.040719,  0.037623,  0.035449,
  0.035022,  0.032254,  0.021085,  0.015280,  0.353002,  0.322302,  0.238023,
  0.174887,  0.130831,  0.094217,  0.073249,  0.055308,  0.040213,  0.033475,
  0.028578,  0.024379,  0.021134,  0.017955,  0.017207,  0.016567,  0.014650,
  0.011836,  0.245838,  0.230739,  0.172581,  0.127909,  0.092769,  0.070133,
  0.053009,  0.035489,  0.027619,  0.022279,  0.017952,  0.016305,  0.014345,
  0.012433,  0.010805,  0.009845,  0.009623,  0.008019,  0.178516,  0.170268,
  0.127522,  0.092700,  0.068902,  0.048110,  0.034054,  0.026432,  0.020014,
  0.014301,  0.012149,  0.011950,  0.012031,  0.013364,  0.009456,  0.008477,
  0.008477,  0.007070,  0.132411,  0.126829,  0.092693,  0.068895,  0.047694,
  0.032383,  0.026011,  0.019928,  0.013674,  0.010013,  0.008868,  0.009095,
  0.008532,  0.007607,  0.006655,  0.007667,  0.007967,  0.006649,  0.097664,
  0.091906,  0.068890,  0.047686,  0.031919,  0.024152,  0.019464,  0.013665,
  0.009934,  0.008107,  0.006770,  0.006583,  0.006336,  0.005713,  0.005629,
  0.005948,  0.006013,  0.005016,  0.070957,  0.068482,  0.047685,  0.031917,
  0.024037,  0.019005,  0.013551,  0.009933,  0.008104,  0.006666,  0.006121,
  0.006052,  0.005963,  0.005620,  0.005535,  0.005539,  0.005526,  0.004609,
  0.053665,  0.047454,  0.031917,  0.024037,  0.019003,  0.013544,  0.009931,
  0.008104,  0.006666,  0.006119,  0.006046,  0.006042,  0.005958,  0.005620,
  0.005534,  0.005533,  0.005519,  0.004603,  0.036016,  0.031740,  0.024037,
  0.019003,  0.013544,  0.009931,  0.008104,  0.006666,  0.006119,  0.006046,
  0.006044,  0.006043,  0.005972,  0.005676,  0.005548,  0.005533,  0.005519,
  0.004603,  0.024022,  0.023917,  0.019003,  0.013544,  0.009931,  0.008104,
  0.006666,  0.006119,  0.006046,  0.006044,  0.006044,  0.006043,  0.006028,
  0.005901,  0.005604,  0.005534,  0.005519,  0.004603,  0.018890,  0.018901,
  0.013523,  0.009914,  0.008089,  0.006650,  0.006104,  0.006030,  0.006028,
  0.006028,  0.006028,  0.006028,  0.006027,  0.005942,  0.005604,  0.005520,
  0.005504,  0.004591,  0.013391,  0.011887,  0.008567,  0.007022,  0.005649,
  0.005103,  0.005030,  0.005028,  0.005028,  0.005028,  0.005028,  0.005028,
  0.005027,  0.004956,  0.004674,  0.004604,  0.004591,  0.003829,
};

void av1_model_rd_surffit(BLOCK_SIZE bsize, double xm, double yl,
                          double *rate_f, double *dist_f) {
  const double x_start = -0.5;
  const double x_end = 16.5;
  const double x_step = 1.0;
  const double y_start = -15.5;
  const double y_end = 16.5;
  const double y_step = 1.0;
  const double epsilon = 1e-6;
  const int stride = (int)rint((x_end - x_start) / x_step) + 1;
  const int cat = bsize_model_cat_lookup[bsize];
  (void)y_end;

  xm = AOMMAX(xm, x_start + x_step + epsilon);
  xm = AOMMIN(xm, x_end - x_step - epsilon);
  yl = AOMMAX(yl, y_start + y_step + epsilon);
  yl = AOMMIN(yl, y_end - y_step - epsilon);

  const double y = (yl - y_start) / y_step;
  const double x = (xm - x_start) / x_step;

  const int yi = (int)floor(y);
  const int xi = (int)floor(x);
  assert(xi > 0);
  assert(yi > 0);

  const double yo = y - yi;
  const double xo = x - xi;
  const double *prate = &interp_rgrid_surf[cat][(yi - 1) * stride + (xi - 1)];
  const double *pdist = &interp_dgrid_surf[(yi - 1) * stride + (xi - 1)];
  *rate_f = interp_bicubic(prate, stride, xo, yo);
  *dist_f = interp_bicubic(pdist, stride, xo, yo);
}

static const double interp_rgrid_curv[4][65] = {
  {
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      41.525408,   51.597692,   49.566271,   54.632979,   60.321507,
      67.730678,   75.766165,   85.324032,   96.600012,   120.839562,
      173.917577,  255.974908,  354.107573,  458.063476,  562.345966,
      668.568424,  772.072881,  878.598490,  982.202274,  1082.708946,
      1188.037853, 1287.702240, 1395.588773, 1490.825830, 1584.231230,
      1691.386090, 1766.822555, 1869.630904, 1926.743565, 2002.949495,
      2047.431137, 2138.486068, 2154.743767, 2209.242472, 2278.252010,
      2298.028834, 2302.326180, 2293.979995, 2275.826226, 2250.700821,
      2221.439725, 2190.878887, 2161.854252, 2137.201768, 2119.757381,
      2112.357039, 2117.836689, 2139.032277, 2178.779750, 2239.915056,
  },
  {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      11.561347,    12.578139,    14.205101,    16.770584,    19.094853,
      21.330863,    23.298907,    26.901921,    34.501017,    57.891733,
      112.234763,   194.853189,   288.302032,   380.499422,   472.625309,
      560.226809,   647.928463,   734.155122,   817.489721,   906.265783,
      999.260562,   1094.489206,  1197.062998,  1293.296825,  1378.926484,
      1472.760990,  1552.663779,  1635.196884,  1692.451951,  1759.741063,
      1822.162720,  1916.515921,  1966.686071,  2031.647506,  2031.381029,
      2067.971335,  2203.662704,  2500.257936,  3019.559830,  3823.371186,
      4973.494802,  6531.733478,  8559.890013,  11119.767206, 14273.167855,
      18081.894761, 22607.750723, 27912.538538, 34058.061008, 41106.120930,
  },
  {
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      0.000000,    0.000000,    0.000000,    0.000000,    0.000000,
      3.281800,    3.765589,    4.342578,    5.145582,    5.611038,
      6.642238,    7.945977,    11.800522,   17.346624,   37.501413,
      87.216800,   165.860942,  253.865564,  332.039345,  408.518863,
      478.120452,  547.268590,  616.067676,  680.022540,  753.863541,
      834.529973,  919.489191,  1008.264989, 1092.230318, 1173.971886,
      1249.514122, 1330.510941, 1399.523249, 1466.923387, 1530.533471,
      1586.515722, 1695.197774, 1746.648696, 1837.136959, 1909.056910,
      1974.948082, 2063.374132, 2178.496387, 2324.476176, 2505.474827,
      2725.653666, 2989.174023, 3300.197225, 3662.884600, 4081.397476,
      4559.897180, 5102.545042, 5713.502387, 6396.930546, 7156.990844,
  },
  {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.614483,     0.842937,     1.050824,     1.326663,     1.717750,
      2.530591,     3.582302,     6.995373,     9.973335,     24.042464,
      56.598240,    113.680735,   180.018689,   231.050567,   266.101082,
      294.957934,   323.326511,   349.434429,   380.443211,   408.171987,
      441.214916,   475.716772,   512.900000,   551.186939,   592.364455,
      624.527378,   661.940693,   679.185473,   724.800679,   764.781792,
      873.050019,   950.299001,   939.292954,   1052.406153,  1030.816617,
      1086.316710,  1275.467594,  1671.923018,  2349.336727,  3381.362469,
      4841.653990,  6803.865037,  9341.649358,  12528.660698, 16438.552805,
      21144.979426, 26721.594308, 33242.051197, 40780.003840, 49409.105984,
  },
};

static const double interp_dgrid_curv[65] = {
  14.604855, 14.604855, 14.604855, 14.604855, 14.604855, 14.604855, 14.604855,
  14.604855, 14.604855, 14.604855, 14.604855, 14.604855, 14.555776, 14.533692,
  14.439920, 14.257791, 13.977230, 13.623229, 13.064884, 12.355411, 11.560773,
  10.728960, 9.861975,  8.643612,  6.916021,  5.154769,  3.734940,  2.680051,
  1.925506,  1.408410,  1.042223,  0.767641,  0.565392,  0.420116,  0.310427,
  0.231711,  0.172999,  0.128293,  0.094992,  0.072171,  0.052972,  0.039354,
  0.029555,  0.022857,  0.016832,  0.013297,  0.000000,  0.000000,  0.000000,
  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
  0.000000,  0.000000,
};

void av1_model_rd_curvfit(BLOCK_SIZE bsize, double xqr, double *rate_f,
                          double *distbysse_f) {
  const double x_start = -15.5;
  const double x_end = 16.5;
  const double x_step = 0.5;
  const double epsilon = 1e-6;
  const int cat = bsize_model_cat_lookup[bsize];
  (void)x_end;

  xqr = AOMMAX(xqr, x_start + x_step + epsilon);
  xqr = AOMMIN(xqr, x_end - x_step - epsilon);
  const double x = (xqr - x_start) / x_step;
  const int xi = (int)floor(x);
  const double xo = x - xi;

  assert(xi > 0);

  const double *prate = &interp_rgrid_curv[cat][(xi - 1)];
  const double *pdist = &interp_dgrid_curv[(xi - 1)];
  *rate_f = interp_cubic(prate, xo);
  *distbysse_f = interp_cubic(pdist, xo);
}

static void get_entropy_contexts_plane(BLOCK_SIZE plane_bsize,
                                       const struct macroblockd_plane *pd,
                                       ENTROPY_CONTEXT t_above[MAX_MIB_SIZE],
                                       ENTROPY_CONTEXT t_left[MAX_MIB_SIZE]) {
  const int num_4x4_w = block_size_wide[plane_bsize] >> tx_size_wide_log2[0];
  const int num_4x4_h = block_size_high[plane_bsize] >> tx_size_high_log2[0];
  const ENTROPY_CONTEXT *const above = pd->above_context;
  const ENTROPY_CONTEXT *const left = pd->left_context;

  memcpy(t_above, above, sizeof(ENTROPY_CONTEXT) * num_4x4_w);
  memcpy(t_left, left, sizeof(ENTROPY_CONTEXT) * num_4x4_h);
}

void av1_get_entropy_contexts(BLOCK_SIZE bsize,
                              const struct macroblockd_plane *pd,
                              ENTROPY_CONTEXT t_above[MAX_MIB_SIZE],
                              ENTROPY_CONTEXT t_left[MAX_MIB_SIZE]) {
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(bsize, pd->subsampling_x, pd->subsampling_y);
  get_entropy_contexts_plane(plane_bsize, pd, t_above, t_left);
}

void av1_mv_pred(const AV1_COMP *cpi, MACROBLOCK *x, uint8_t *ref_y_buffer,
                 int ref_y_stride, int ref_frame, BLOCK_SIZE block_size) {
  int i;
  int zero_seen = 0;
  int best_sad = INT_MAX;
  int this_sad = INT_MAX;
  int max_mv = 0;
  uint8_t *src_y_ptr = x->plane[0].src.buf;
  uint8_t *ref_y_ptr;
  MV pred_mv[MAX_MV_REF_CANDIDATES + 1];
  int num_mv_refs = 0;
  const MV_REFERENCE_FRAME ref_frames[2] = { ref_frame, NONE_FRAME };
  const int_mv ref_mv =
      av1_get_ref_mv_from_stack(0, ref_frames, 0, x->mbmi_ext);
  const int_mv ref_mv1 =
      av1_get_ref_mv_from_stack(0, ref_frames, 1, x->mbmi_ext);

  pred_mv[num_mv_refs++] = ref_mv.as_mv;
  if (ref_mv.as_int != ref_mv1.as_int) {
    pred_mv[num_mv_refs++] = ref_mv1.as_mv;
  }
  if (cpi->sf.adaptive_motion_search && block_size < x->max_partition_size)
    pred_mv[num_mv_refs++] = x->pred_mv[ref_frame];

  assert(num_mv_refs <= (int)(sizeof(pred_mv) / sizeof(pred_mv[0])));

  // Get the sad for each candidate reference mv.
  for (i = 0; i < num_mv_refs; ++i) {
    const MV *this_mv = &pred_mv[i];
    int fp_row, fp_col;
    fp_row = (this_mv->row + 3 + (this_mv->row >= 0)) >> 3;
    fp_col = (this_mv->col + 3 + (this_mv->col >= 0)) >> 3;
    max_mv = AOMMAX(max_mv, AOMMAX(abs(this_mv->row), abs(this_mv->col)) >> 3);

    if (fp_row == 0 && fp_col == 0 && zero_seen) continue;
    zero_seen |= (fp_row == 0 && fp_col == 0);

    ref_y_ptr = &ref_y_buffer[ref_y_stride * fp_row + fp_col];
    // Find sad for current vector.
    this_sad = cpi->fn_ptr[block_size].sdf(src_y_ptr, x->plane[0].src.stride,
                                           ref_y_ptr, ref_y_stride);
    // Note if it is the best so far.
    if (this_sad < best_sad) {
      best_sad = this_sad;
    }
  }

  // Note the index of the mv that worked best in the reference list.
  x->max_mv_context[ref_frame] = max_mv;
  x->pred_mv_sad[ref_frame] = best_sad;
}

void av1_setup_pred_block(const MACROBLOCKD *xd,
                          struct buf_2d dst[MAX_MB_PLANE],
                          const YV12_BUFFER_CONFIG *src, int mi_row, int mi_col,
                          const struct scale_factors *scale,
                          const struct scale_factors *scale_uv,
                          const int num_planes) {
  int i;

  dst[0].buf = src->y_buffer;
  dst[0].stride = src->y_stride;
  dst[1].buf = src->u_buffer;
  dst[2].buf = src->v_buffer;
  dst[1].stride = dst[2].stride = src->uv_stride;

  for (i = 0; i < num_planes; ++i) {
    setup_pred_plane(dst + i, xd->mi[0]->sb_type, dst[i].buf,
                     i ? src->uv_crop_width : src->y_crop_width,
                     i ? src->uv_crop_height : src->y_crop_height,
                     dst[i].stride, mi_row, mi_col, i ? scale_uv : scale,
                     xd->plane[i].subsampling_x, xd->plane[i].subsampling_y);
  }
}

int av1_raster_block_offset(BLOCK_SIZE plane_bsize, int raster_block,
                            int stride) {
  const int bw = mi_size_wide_log2[plane_bsize];
  const int y = 4 * (raster_block >> bw);
  const int x = 4 * (raster_block & ((1 << bw) - 1));
  return y * stride + x;
}

int16_t *av1_raster_block_offset_int16(BLOCK_SIZE plane_bsize, int raster_block,
                                       int16_t *base) {
  const int stride = block_size_wide[plane_bsize];
  return base + av1_raster_block_offset(plane_bsize, raster_block, stride);
}

YV12_BUFFER_CONFIG *av1_get_scaled_ref_frame(const AV1_COMP *cpi,
                                             int ref_frame) {
  assert(ref_frame >= LAST_FRAME && ref_frame <= ALTREF_FRAME);
  RefCntBuffer *const scaled_buf = cpi->scaled_ref_buf[ref_frame - 1];
  const RefCntBuffer *const ref_buf =
      get_ref_frame_buf(&cpi->common, ref_frame);
  return (scaled_buf != ref_buf && scaled_buf != NULL) ? &scaled_buf->buf
                                                       : NULL;
}

int av1_get_switchable_rate(const AV1_COMMON *const cm, MACROBLOCK *x,
                            const MACROBLOCKD *xd) {
  if (cm->interp_filter == SWITCHABLE) {
    const MB_MODE_INFO *const mbmi = xd->mi[0];
    int inter_filter_cost = 0;
    int dir;

    for (dir = 0; dir < 2; ++dir) {
      const int ctx = av1_get_pred_context_switchable_interp(xd, dir);
      const InterpFilter filter =
          av1_extract_interp_filter(mbmi->interp_filters, dir);
      inter_filter_cost += x->switchable_interp_costs[ctx][filter];
    }
    return SWITCHABLE_INTERP_RATE_FACTOR * inter_filter_cost;
  } else {
    return 0;
  }
}

void av1_set_rd_speed_thresholds(AV1_COMP *cpi) {
  int i;
  RD_OPT *const rd = &cpi->rd;
  SPEED_FEATURES *const sf = &cpi->sf;

  // Set baseline threshold values.
  for (i = 0; i < MAX_MODES; ++i) rd->thresh_mult[i] = cpi->oxcf.mode == 0;

  if (sf->adaptive_rd_thresh) {
    rd->thresh_mult[THR_NEARESTMV] = 300;
    rd->thresh_mult[THR_NEARESTL2] = 300;
    rd->thresh_mult[THR_NEARESTL3] = 300;
    rd->thresh_mult[THR_NEARESTB] = 300;
    rd->thresh_mult[THR_NEARESTA2] = 300;
    rd->thresh_mult[THR_NEARESTA] = 300;
    rd->thresh_mult[THR_NEARESTG] = 300;
  } else {
    rd->thresh_mult[THR_NEARESTMV] = 0;
    rd->thresh_mult[THR_NEARESTL2] = 0;
    rd->thresh_mult[THR_NEARESTL3] = 0;
    rd->thresh_mult[THR_NEARESTB] = 0;
    rd->thresh_mult[THR_NEARESTA2] = 0;
    rd->thresh_mult[THR_NEARESTA] = 0;
    rd->thresh_mult[THR_NEARESTG] = 0;
  }

  rd->thresh_mult[THR_NEWMV] += 1000;
  rd->thresh_mult[THR_NEWL2] += 1000;
  rd->thresh_mult[THR_NEWL3] += 1000;
  rd->thresh_mult[THR_NEWB] += 1000;
  rd->thresh_mult[THR_NEWA2] = 1000;
  rd->thresh_mult[THR_NEWA] += 1000;
  rd->thresh_mult[THR_NEWG] += 1000;

  rd->thresh_mult[THR_NEARMV] += 1000;
  rd->thresh_mult[THR_NEARL2] += 1000;
  rd->thresh_mult[THR_NEARL3] += 1000;
  rd->thresh_mult[THR_NEARB] += 1000;
  rd->thresh_mult[THR_NEARA2] = 1000;
  rd->thresh_mult[THR_NEARA] += 1000;
  rd->thresh_mult[THR_NEARG] += 1000;

  rd->thresh_mult[THR_GLOBALMV] += 2000;
  rd->thresh_mult[THR_GLOBALL2] += 2000;
  rd->thresh_mult[THR_GLOBALL3] += 2000;
  rd->thresh_mult[THR_GLOBALB] += 2000;
  rd->thresh_mult[THR_GLOBALA2] = 2000;
  rd->thresh_mult[THR_GLOBALG] += 2000;
  rd->thresh_mult[THR_GLOBALA] += 2000;

  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLA] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL2A] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL3A] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTGA] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLB] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL2B] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL3B] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTGB] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLA2] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL2A2] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTL3A2] += 1000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTGA2] += 1000;

  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLL2] += 2000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLL3] += 2000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTLG] += 2000;
  rd->thresh_mult[THR_COMP_NEAREST_NEARESTBA] += 2000;

  rd->thresh_mult[THR_COMP_NEAR_NEARLA] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLA] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLA] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWLA] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARLA] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWLA] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLA] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL2A] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL2A] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL2A] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL2A] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL2A] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL2A] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL2A] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL3A] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL3A] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL3A] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL3A] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL3A] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL3A] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL3A] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARGA] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWGA] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTGA] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWGA] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARGA] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWGA] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALGA] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARLB] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLB] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLB] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWLB] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARLB] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWLB] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLB] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL2B] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL2B] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL2B] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL2B] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL2B] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL2B] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL2B] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL3B] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL3B] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL3B] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL3B] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL3B] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL3B] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL3B] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARGB] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWGB] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTGB] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWGB] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARGB] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWGB] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALGB] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARLA2] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLA2] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLA2] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWLA2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARLA2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWLA2] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLA2] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL2A2] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL2A2] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL2A2] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL2A2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL2A2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL2A2] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL2A2] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARL3A2] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWL3A2] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTL3A2] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWL3A2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARL3A2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWL3A2] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALL3A2] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARGA2] += 1200;
  rd->thresh_mult[THR_COMP_NEAREST_NEWGA2] += 1500;
  rd->thresh_mult[THR_COMP_NEW_NEARESTGA2] += 1500;
  rd->thresh_mult[THR_COMP_NEAR_NEWGA2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEARGA2] += 1700;
  rd->thresh_mult[THR_COMP_NEW_NEWGA2] += 2000;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALGA2] += 2500;

  rd->thresh_mult[THR_COMP_NEAR_NEARLL2] += 1600;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLL2] += 2000;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLL2] += 2000;
  rd->thresh_mult[THR_COMP_NEAR_NEWLL2] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEARLL2] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEWLL2] += 2400;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLL2] += 3200;

  rd->thresh_mult[THR_COMP_NEAR_NEARLL3] += 1600;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLL3] += 2000;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLL3] += 2000;
  rd->thresh_mult[THR_COMP_NEAR_NEWLL3] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEARLL3] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEWLL3] += 2400;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLL3] += 3200;

  rd->thresh_mult[THR_COMP_NEAR_NEARLG] += 1600;
  rd->thresh_mult[THR_COMP_NEAREST_NEWLG] += 2000;
  rd->thresh_mult[THR_COMP_NEW_NEARESTLG] += 2000;
  rd->thresh_mult[THR_COMP_NEAR_NEWLG] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEARLG] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEWLG] += 2400;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALLG] += 3200;

  rd->thresh_mult[THR_COMP_NEAR_NEARBA] += 1600;
  rd->thresh_mult[THR_COMP_NEAREST_NEWBA] += 2000;
  rd->thresh_mult[THR_COMP_NEW_NEARESTBA] += 2000;
  rd->thresh_mult[THR_COMP_NEAR_NEWBA] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEARBA] += 2200;
  rd->thresh_mult[THR_COMP_NEW_NEWBA] += 2400;
  rd->thresh_mult[THR_COMP_GLOBAL_GLOBALBA] += 3200;

  rd->thresh_mult[THR_DC] += 1000;
  rd->thresh_mult[THR_PAETH] += 1000;
  rd->thresh_mult[THR_SMOOTH] += 2000;
  rd->thresh_mult[THR_SMOOTH_V] += 2000;
  rd->thresh_mult[THR_SMOOTH_H] += 2000;
  rd->thresh_mult[THR_H_PRED] += 2000;
  rd->thresh_mult[THR_V_PRED] += 2000;
  rd->thresh_mult[THR_D135_PRED] += 2500;
  rd->thresh_mult[THR_D203_PRED] += 2500;
  rd->thresh_mult[THR_D157_PRED] += 2500;
  rd->thresh_mult[THR_D67_PRED] += 2500;
  rd->thresh_mult[THR_D113_PRED] += 2500;
  rd->thresh_mult[THR_D45_PRED] += 2500;
}

void av1_set_rd_speed_thresholds_sub8x8(AV1_COMP *cpi) {
  static const int thresh_mult[MAX_REFS] = { 2500, 2500, 2500, 2500, 2500,
                                             2500, 2500, 4500, 4500, 4500,
                                             4500, 4500, 4500, 4500, 4500,
                                             4500, 4500, 4500, 4500, 2500 };
  RD_OPT *const rd = &cpi->rd;
  memcpy(rd->thresh_mult_sub8x8, thresh_mult, sizeof(thresh_mult));
}

void av1_update_rd_thresh_fact(const AV1_COMMON *const cm,
                               int (*factor_buf)[MAX_MODES], int rd_thresh,
                               int bsize, int best_mode_index) {
  if (rd_thresh > 0) {
    const int top_mode = MAX_MODES;
    int mode;
    for (mode = 0; mode < top_mode; ++mode) {
      const BLOCK_SIZE min_size = AOMMAX(bsize - 1, BLOCK_4X4);
      const BLOCK_SIZE max_size =
          AOMMIN(bsize + 2, (int)cm->seq_params.sb_size);
      BLOCK_SIZE bs;
      for (bs = min_size; bs <= max_size; ++bs) {
        int *const fact = &factor_buf[bs][mode];
        if (mode == best_mode_index) {
          *fact -= (*fact >> 4);
        } else {
          *fact = AOMMIN(*fact + RD_THRESH_INC, rd_thresh * RD_THRESH_MAX_FACT);
        }
      }
    }
  }
}

int av1_get_intra_cost_penalty(int qindex, int qdelta,
                               aom_bit_depth_t bit_depth) {
  const int q = av1_dc_quant_Q3(qindex, qdelta, bit_depth);
  switch (bit_depth) {
    case AOM_BITS_8: return 20 * q;
    case AOM_BITS_10: return 5 * q;
    case AOM_BITS_12: return ROUND_POWER_OF_TWO(5 * q, 2);
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
}
