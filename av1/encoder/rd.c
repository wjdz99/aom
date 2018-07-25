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

  if (cm->skip_mode_flag) {
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

int av1_compute_rd_mult(const AV1_COMP *cpi, int qindex) {
  const int64_t q =
      av1_dc_quant_Q3(qindex, 0, cpi->common.seq_params.bit_depth);
  int64_t rdmult = 0;
  switch (cpi->common.seq_params.bit_depth) {
    case AOM_BITS_8: rdmult = 88 * q * q / 24; break;
    case AOM_BITS_10: rdmult = ROUND_POWER_OF_TWO(88 * q * q / 24, 4); break;
    case AOM_BITS_12: rdmult = ROUND_POWER_OF_TWO(88 * q * q / 24, 8); break;
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
  if (cpi->oxcf.pass == 2 && (cpi->common.frame_type != KEY_FRAME)) {
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

void av1_set_mvcost(MACROBLOCK *x, int ref, int ref_mv_idx) {
  (void)ref;
  (void)ref_mv_idx;
  x->mvcost = x->mv_cost_stack;
  x->nmvjointcost = x->nmv_vec_cost;
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

void av1_initialize_rd_consts(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &cpi->td.mb;
  RD_OPT *const rd = &cpi->rd;

  aom_clear_system_state();

  rd->RDMULT = av1_compute_rd_mult(cpi, cm->base_qindex + cm->y_dc_delta_q);

  set_error_per_bit(x, rd->RDMULT);

  set_block_thresholds(cm, rd);

  if (cm->cur_frame_force_integer_mv) {
    av1_build_nmv_cost_table(x->nmv_vec_cost, x->nmvcost, &cm->fc->nmvc,
                             MV_SUBPEL_NONE);
  } else {
    av1_build_nmv_cost_table(
        x->nmv_vec_cost,
        cm->allow_high_precision_mv ? x->nmvcost_hp : x->nmvcost, &cm->fc->nmvc,
        cm->allow_high_precision_mv);
  }

  x->mvcost = x->mv_cost_stack;
  x->nmvjointcost = x->nmv_vec_cost;

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

void av1_model_rd_surffit(double xm, double yl, double *rate_f,
                          double *dist_f) {
  static const double interp_rgrid[33 * 35] = {
    0.038171,    0.049072,    0.052072,    0.084358,    0.162366,
    0.187087,    0.170802,    0.162876,    0.162198,    0.162187,
    0.162187,    0.162187,    0.162187,    0.162187,    0.162187,
    0.162227,    0.162800,    0.164810,    0.175632,    0.217940,
    0.285146,    0.332011,    0.379473,    0.478469,    0.588887,
    0.755301,    0.963999,    1.214563,    1.518387,    1.939313,
    2.445811,    3.180269,    4.077490,    5.275091,    4.948256,
    0.049072,    0.063084,    0.066942,    0.108447,    0.208881,
    0.242669,    0.224823,    0.211543,    0.208665,    0.208501,
    0.208501,    0.208501,    0.208501,    0.208541,    0.209114,
    0.211124,    0.220551,    0.244179,    0.283204,    0.367077,
    0.453293,    0.562537,    0.682172,    0.870435,    1.100440,
    1.392839,    1.750432,    2.229129,    2.805020,    3.638528,
    4.658296,    6.155703,    8.083745,    10.339025,   9.538375,
    0.049829,    0.064058,    0.067931,    0.109481,    0.210916,
    0.251101,    0.241335,    0.220497,    0.212408,    0.211730,
    0.211719,    0.211759,    0.212332,    0.214342,    0.223769,
    0.247397,    0.285848,    0.362042,    0.436477,    0.558116,
    0.686825,    0.876276,    1.106922,    1.400839,    1.759913,
    2.241516,    2.821143,    3.658665,    4.683255,    6.187545,
    8.124102,    10.665043,   13.896879,   17.700488,   16.195219,
    0.049087,    0.063316,    0.067084,    0.100610,    0.189528,
    0.244625,    0.248640,    0.228941,    0.215127,    0.211948,
    0.212332,    0.214342,    0.223769,    0.247397,    0.285848,
    0.362042,    0.436437,    0.557542,    0.685430,    0.875703,
    1.106882,    1.400839,    1.759913,    2.241516,    2.821143,
    3.658665,    4.683255,    6.187545,    8.124102,    10.665043,
    13.896879,   18.162710,   23.413544,   29.545185,   26.811907,
    0.038397,    0.052568,    0.062605,    0.083901,    0.148597,
    0.229011,    0.259640,    0.257320,    0.230413,    0.220965,
    0.234724,    0.263523,    0.291909,    0.362456,    0.436437,
    0.557542,    0.685430,    0.875703,    1.106882,    1.400839,
    1.759913,    2.241516,    2.821143,    3.658665,    4.683255,
    6.187545,    8.124102,    10.665043,   13.896879,   18.162710,
    23.413544,   30.298334,   38.534352,   47.685325,   42.323100,
    0.012332,    0.025733,    0.056178,    0.099613,    0.173321,
    0.247341,    0.300856,    0.342398,    0.323078,    0.347921,
    0.494074,    0.683892,    0.598268,    0.602648,    0.695324,
    0.876314,    1.106882,    1.400839,    1.759913,    2.241516,
    2.821143,    3.658665,    4.683255,    6.187545,    8.124102,
    10.665043,   13.896879,   18.162710,   23.413544,   30.298334,
    38.534352,   48.846370,   60.363101,   73.254956,   64.164358,
    0.002246,    0.015029,    0.068712,    0.171481,    0.271946,
    0.339499,    0.479767,    0.644396,    0.658373,    0.836457,
    1.581683,    2.452201,    2.086447,    1.694026,    1.447640,
    1.497531,    1.778926,    2.242625,    2.821143,    3.658665,
    4.683255,    6.187545,    8.124102,    10.665043,   13.896879,
    18.162710,   23.413544,   30.298334,   38.534352,   48.846370,
    60.363101,   74.939096,   90.327939,   106.977142,  92.977795,
    0.012335,    0.033611,    0.135349,    0.342848,    0.486845,
    0.596723,    0.915178,    1.247348,    1.327046,    1.754073,
    3.154289,    4.872001,    5.541837,    6.325433,    5.082333,
    3.649384,    3.166897,    3.708390,    4.691986,    6.188092,
    8.124102,    10.665043,   13.896879,   18.162710,   23.413544,
    30.298334,   38.534352,   48.846370,   60.363101,   74.939096,
    90.327939,   109.313474,  129.251357,  149.313642,  128.293539,
    0.040354,    0.084205,    0.251028,    0.554419,    0.776738,
    1.026705,    1.414103,    1.832330,    2.314690,    2.997173,
    3.919792,    5.240364,    8.188776,    12.722185,   11.911753,
    8.779421,    6.679524,    6.971329,    8.444145,    10.784412,
    13.926373,   18.164552,   23.413544,   30.298334,   38.534352,
    48.846370,   60.363101,   74.939096,   90.327939,   109.313474,
    129.251357,  152.275141,  174.105269,  190.049693,  156.946750,
    0.107447,    0.167574,    0.330945,    0.676585,    1.079902,
    1.531304,    2.003952,    2.834885,    3.827170,    4.956427,
    6.494385,    8.170701,    10.482754,   14.104276,   16.214422,
    16.660132,   15.987997,   16.368765,   17.657407,   20.357016,
    24.264229,   30.520948,   38.575758,   48.848745,   60.363101,
    74.939096,   90.327939,   109.313474,  129.251357,  152.275141,
    174.105269,  193.418962,  208.870054,  216.230173,  171.862790,
    0.341098,    0.415157,    0.481479,    0.878803,    1.415967,
    2.131000,    3.306241,    5.365666,    7.779505,    10.934053,
    14.804469,   18.404667,   21.899245,   26.170660,   30.929966,
    36.084083,   39.216197,   40.621465,   41.374039,   43.796512,
    46.918406,   52.654644,   61.527294,   75.248979,   90.395325,
    109.317589,  129.251357,  152.272646,  174.066860,  193.295688,
    208.746781,  219.639472,  224.931174,  223.149723,  173.888032,
    0.540513,    0.697446,    0.989957,    1.965982,    3.148747,
    5.328609,    9.041807,    15.368045,   23.561619,   32.406497,
    40.853709,   49.165669,   56.869934,   63.790466,   71.372008,
    79.883657,   85.696937,   88.690570,   89.473356,   92.502817,
    95.110386,   96.563245,   101.094315,  114.324137,  130.911331,
    152.662099,  174.072034,  193.229561,  208.191492,  217.865325,
    223.157027,  225.964523,  225.685821,  219.481112,  169.808050,
    0.498783,    1.148068,    3.124141,    7.808820,    15.281308,
    27.140175,   40.473653,   58.070288,   77.253440,   95.763726,
    111.722796,  124.075039,  133.559202,  141.363819,  148.659031,
    154.316640,  159.181013,  163.623815,  164.900391,  167.100743,
    168.831105,  167.898083,  168.104387,  177.525209,  188.104312,
    199.897282,  209.763779,  218.264847,  222.025274,  221.664651,
    221.370341,  221.468162,  221.001602,  216.415871,  167.949760,
    1.747112,    4.532788,    14.562938,   39.338122,   76.652684,
    115.565748,  145.146701,  174.608599,  197.986056,  217.437797,
    234.289462,  245.196131,  250.996376,  253.886997,  254.707487,
    255.072534,  256.956564,  259.249479,  259.510152,  257.365560,
    255.545813,  254.214135,  254.004717,  259.521894,  261.791511,
    257.612065,  246.633008,  236.947743,  227.439013,  223.014990,
    221.848504,  220.470540,  219.428533,  215.755413,  167.790572,
    7.815572,    16.246922,   49.547590,   129.851278,  230.959693,
    305.159949,  344.005687,  367.313988,  375.701807,  382.346383,
    387.136116,  387.091152,  384.434876,  381.107310,  377.546049,
    375.653809,  375.515615,  376.890815,  376.118853,  368.495019,
    362.277366,  352.176659,  348.494588,  347.954779,  345.995677,
    337.811084,  323.415069,  301.955510,  278.456251,  270.141129,
    263.013756,  243.761907,  229.286920,  222.168742,  172.644105,
    21.934681,   50.004338,   132.763281,  279.366550,  425.227348,
    524.648519,  577.438640,  591.179884,  581.671429,  573.611764,
    564.917647,  552.397941,  536.875359,  523.145442,  515.010717,
    513.537022,  514.293613,  517.389016,  522.072517,  519.593681,
    503.256918,  475.115969,  462.228083,  451.558872,  442.993486,
    435.206120,  432.689614,  426.213530,  418.651638,  419.814900,
    407.200911,  361.283520,  326.129371,  312.335075,  242.531670,
    75.986081,   175.921721,  322.356766,  500.634776,  643.831711,
    733.799267,  783.805467,  801.014867,  797.696092,  785.559278,
    767.348736,  744.402107,  717.370246,  692.430609,  677.016262,
    670.921939,  670.792096,  677.646659,  688.543430,  693.717124,
    685.557417,  657.393128,  621.492018,  585.066441,  565.858916,
    559.094010,  557.009881,  556.519086,  560.678682,  573.698147,
    577.177972,  560.502666,  546.184575,  534.249671,  415.403124,
    246.355125,  460.803043,  608.765500,  756.413120,  882.442614,
    952.454434,  987.232377,  1000.615592, 1002.686101, 997.308471,
    983.729618,  957.783123,  922.524263,  887.394326,  864.137545,
    852.409291,  850.024753,  856.014441,  866.303600,  876.931269,
    882.140915,  864.608827,  824.832334,  771.420997,  737.166989,
    722.792348,  694.342826,  681.880598,  682.980272,  688.585636,
    692.308127,  691.315450,  688.577595,  677.131804,  526.669777,
    483.653682,  742.358204,  877.049422,  991.487111,  1084.394490,
    1144.539097, 1182.968962, 1204.239613, 1211.200946, 1206.247331,
    1193.072795, 1170.920373, 1140.261588, 1106.305912, 1076.101332,
    1054.510698, 1048.154993, 1053.396518, 1063.935591, 1074.795605,
    1078.891689, 1067.592355, 1037.881748, 981.108723,  935.094900,
    920.917961,  893.443584,  880.511252,  879.721308,  880.113007,
    879.669006,  869.105529,  843.231159,  819.948140,  637.239665,
    697.772566,  977.970594,  1083.496381, 1161.652937, 1230.167830,
    1288.203239, 1336.632247, 1378.589994, 1408.103456, 1417.056088,
    1411.039575, 1387.585710, 1354.452436, 1322.009452, 1293.616109,
    1275.253155, 1265.943075, 1262.224098, 1272.188377, 1285.066609,
    1290.675713, 1280.336299, 1247.453278, 1187.067117, 1132.124350,
    1114.630020, 1106.169913, 1102.638436, 1102.385639, 1102.429208,
    1101.258090, 1076.904342, 1014.582685, 973.746014,  756.051025,
    851.306188,  1146.317556, 1230.198236, 1290.864542, 1346.447838,
    1398.504631, 1445.133157, 1497.754966, 1559.775380, 1603.069267,
    1620.291950, 1605.403547, 1578.286990, 1542.354326, 1508.786274,
    1492.063928, 1487.939788, 1489.682895, 1493.072159, 1494.613859,
    1507.823509, 1505.584143, 1470.309223, 1406.085261, 1336.644514,
    1302.469253, 1308.813296, 1314.926052, 1315.392774, 1316.019817,
    1324.304921, 1335.616811, 1318.762177, 1288.868777, 1001.999901,
    959.016679,  1273.168942, 1346.455320, 1398.631130, 1445.650788,
    1498.478317, 1562.909296, 1623.142072, 1680.447605, 1721.810303,
    1755.204397, 1769.211406, 1778.104297, 1766.642422, 1738.741197,
    1717.857937, 1710.236505, 1714.258279, 1717.474634, 1724.019280,
    1732.787334, 1720.169420, 1685.854917, 1626.922978, 1545.843160,
    1492.430146, 1505.559346, 1517.923903, 1519.054991, 1523.198435,
    1551.466726, 1606.740002, 1627.074821, 1603.118717, 1246.981214,
    1044.992484, 1378.986148, 1445.650788, 1498.478317, 1562.905876,
    1623.101542, 1680.597755, 1724.535050, 1766.846125, 1799.291326,
    1840.556551, 1876.141487, 1911.036306, 1947.526342, 1951.018348,
    1950.537521, 1959.568139, 1958.965788, 1953.568405, 1965.656530,
    1972.250001, 1961.400641, 1926.088507, 1857.530470, 1772.756141,
    1709.660824, 1716.571645, 1729.396458, 1733.667828, 1772.030326,
    1872.710999, 1932.343099, 1943.985658, 1915.055316, 1489.667918,
    1118.572430, 1477.173763, 1562.905876, 1623.101542, 1680.597755,
    1724.535050, 1766.853555, 1799.420926, 1841.258263, 1879.045394,
    1920.850254, 1970.362634, 2004.976790, 2062.550568, 2107.204760,
    2153.774211, 2193.917449, 2211.121171, 2228.916129, 2253.909478,
    2254.248218, 2252.916520, 2236.315719, 2187.488868, 2106.225619,
    2037.198160, 2034.871532, 2062.795678, 2079.788012, 2172.205201,
    2395.779419, 2488.960979, 2495.956619, 2458.061921, 1912.057554,
    1214.838274, 1600.442149, 1680.597755, 1724.535050, 1766.853555,
    1799.420926, 1841.258263, 1879.045394, 1920.855141, 1970.452959,
    2005.402346, 2063.330126, 2113.842853, 2179.458446, 2248.611542,
    2327.366181, 2388.724005, 2440.677892, 2524.118750, 2576.307950,
    2576.719056, 2596.020300, 2634.629756, 2640.274818, 2577.285869,
    2509.150877, 2504.030416, 2537.466010, 2553.961840, 2592.600622,
    2684.256946, 2721.872180, 2723.543288, 2679.893111, 2083.891722,
    1302.801328, 1699.942153, 1766.853555, 1799.420926, 1841.258263,
    1879.045394, 1920.855141, 1970.452959, 2005.402346, 2063.330126,
    2113.837696, 2179.413334, 2248.890836, 2329.448358, 2399.641650,
    2478.850560, 2595.235911, 2679.568364, 2740.312802, 2818.355464,
    2863.160355, 2906.002388, 2986.561374, 3024.215688, 2987.515868,
    2907.385333, 2864.343117, 2813.558337, 2792.180701, 2793.278262,
    2799.646860, 2801.317969, 2787.869302, 2713.064868, 2100.024732,
    1365.500937, 1773.320173, 1841.258263, 1879.045394, 1920.855141,
    1970.452959, 2005.402346, 2063.330126, 2113.837696, 2179.413334,
    2248.890836, 2329.448358, 2399.642580, 2478.899019, 2596.027077,
    2684.895418, 2755.948695, 2838.280005, 2882.113399, 2953.215625,
    3073.258508, 3159.401336, 3227.635233, 3263.499499, 3241.985591,
    3120.058224, 3023.869972, 2973.105928, 2956.222602, 2954.098829,
    2954.013289, 2951.709538, 2918.554148, 2795.309988, 2149.114462,
    1422.985339, 1851.860349, 1920.855141, 1970.452959, 2005.402346,
    2063.330126, 2113.837696, 2179.413334, 2248.890836, 2329.448358,
    2399.642580, 2478.899019, 2596.027077, 2684.895418, 2755.963471,
    2838.521352, 2882.838281, 2951.345920, 3064.244306, 3147.849282,
    3230.126596, 3314.414892, 3354.723170, 3367.419188, 3358.623829,
    3276.391240, 3152.166151, 3087.991449, 3040.880858, 3021.866291,
    3017.640683, 3015.462497, 3001.746368, 2923.691083, 2263.864968,
    1481.766942, 1941.894582, 2005.402346, 2063.330126, 2113.837696,
    2179.413334, 2248.890836, 2329.448358, 2399.642580, 2478.899019,
    2596.027077, 2684.895418, 2755.963471, 2838.521352, 2882.838281,
    2951.345920, 3064.227502, 3147.592864, 3229.030079, 3309.705711,
    3348.982095, 3377.057864, 3388.202081, 3389.624541, 3388.020247,
    3367.629100, 3320.328561, 3259.428423, 3152.136201, 3091.236063,
    3043.913532, 3023.135527, 3019.726256, 2971.489970, 2310.716755,
    1540.831361, 2033.066978, 2113.837696, 2179.413334, 2248.890836,
    2329.448358, 2399.642580, 2478.899019, 2596.027077, 2684.895418,
    2755.963471, 2838.521352, 2882.838281, 2951.345920, 3064.227502,
    3147.592864, 3229.030079, 3309.705711, 3348.959495, 3376.754599,
    3387.727422, 3390.154702, 3390.861864, 3390.911000, 3390.825460,
    3389.594375, 3385.368767, 3366.420007, 3320.243021, 3258.197338,
    3134.418526, 3048.139141, 3026.195858, 2976.057446, 2314.053146,
    1628.183900, 2147.759256, 2248.890836, 2329.448358, 2399.642580,
    2478.899019, 2596.027077, 2684.895418, 2755.963471, 2838.521352,
    2882.838281, 2951.345920, 3064.227502, 3147.592864, 3229.030079,
    3309.705711, 3348.959495, 3376.754599, 3387.727422, 3390.154702,
    3390.861864, 3390.911000, 3390.911000, 3390.911000, 3390.911000,
    3390.911000, 3390.825460, 3389.594375, 3385.368767, 3363.425484,
    3277.146098, 3153.367286, 3091.236063, 2998.000729, 2316.515316,
    1706.975639, 2255.427545, 2357.373384, 2435.278083, 2552.894270,
    2640.996371, 2709.013904, 2790.252600, 2833.774121, 2900.561641,
    3012.736061, 3096.052287, 3177.489503, 3258.165135, 3297.418918,
    3325.214022, 3336.186845, 3338.614125, 3339.321287, 3339.370423,
    3339.370423, 3339.370423, 3339.370423, 3339.370423, 3339.370423,
    3339.370423, 3339.370423, 3339.370423, 3339.284883, 3336.822713,
    3316.110515, 3268.787984, 3206.656761, 3037.663003, 2300.739448,
    1369.205637, 1813.683691, 1944.559541, 2013.207515, 2088.250159,
    2158.461204, 2158.079893, 2205.881822, 2306.608780, 2365.169206,
    2436.429007, 2516.397477, 2555.651261, 2583.446364, 2594.419187,
    2596.846468, 2597.553630, 2597.602766, 2597.602766, 2597.602766,
    2597.602766, 2597.602766, 2597.602766, 2597.602766, 2597.602766,
    2597.602766, 2597.602766, 2597.602766, 2597.602766, 2597.517226,
    2596.286140, 2592.060532, 2570.117250, 2448.666423, 1819.934272,
  };
  static const double interp_dgrid[33 * 35] = {
    4.111763,  5.285911,  5.421834,  6.203874,  8.137816,  9.376572,  10.549511,
    11.028043, 11.068966, 11.069607, 11.069607, 11.069607, 11.069607, 11.069607,
    11.069607, 11.069580, 11.069187, 11.067859, 11.061131, 11.034544, 10.999584,
    10.995443, 11.007412, 11.029740, 11.031992, 11.014306, 10.988094, 10.920677,
    10.875280, 10.796029, 10.694138, 10.554224, 10.405829, 10.014510, 7.624060,
    5.285911,  6.795349,  6.970085,  7.975443,  10.452589, 11.923896, 13.245225,
    14.046959, 14.220754, 14.230628, 14.230628, 14.230628, 14.230628, 14.230601,
    14.230208, 14.228880, 14.223106, 14.209299, 14.193618, 14.164159, 14.149616,
    14.165608, 14.171265, 14.162844, 14.140845, 14.074547, 14.021543, 13.923206,
    13.815704, 13.656777, 13.478183, 13.203391, 12.873335, 12.324199, 9.378643,
    5.367495,  6.900230,  7.076585,  8.083020,  10.553868, 11.770764, 12.658513,
    13.920297, 14.408703, 14.449625, 14.450267, 14.450239, 14.449847, 14.448519,
    14.442744, 14.428937, 14.413649, 14.389443, 14.382960, 14.390473, 14.389923,
    14.380735, 14.358973, 14.293318, 14.240607, 14.142348, 14.034318, 13.874064,
    13.695080, 13.418967, 13.086812, 12.711691, 12.290180, 11.632403, 8.777221,
    5.321351,  6.854085,  7.047162,  7.865552,  10.017932, 11.412335, 12.190427,
    13.471605, 14.269633, 14.440576, 14.449847, 14.448519, 14.442744, 14.428937,
    14.413649, 14.389443, 14.382988, 14.390865, 14.390877, 14.381127, 14.359000,
    14.293318, 14.240607, 14.142348, 14.034318, 13.874064, 13.695080, 13.418967,
    13.086812, 12.711691, 12.290180, 11.793243, 11.242660, 10.449251, 7.759961,
    4.656612,  6.180904,  6.788188,  7.401810,  8.936924,  10.936831, 11.860304,
    12.782907, 13.980826, 14.409784, 14.431769, 14.415547, 14.408765, 14.389112,
    14.382988, 14.390865, 14.390877, 14.381127, 14.359000, 14.293318, 14.240607,
    14.142348, 14.034318, 13.874064, 13.695080, 13.418967, 13.086812, 12.711691,
    12.290180, 11.793243, 11.242660, 10.587401, 9.877580,  9.027090,  6.668149,
    3.035871,  4.439617,  5.986609,  7.318365,  9.023158,  10.933023, 11.873896,
    12.719181, 13.825048, 14.243951, 14.224353, 14.145742, 14.273208, 14.369806,
    14.387640, 14.380945, 14.359000, 14.293318, 14.240607, 14.142348, 14.034318,
    13.874064, 13.695080, 13.418967, 13.086812, 12.711691, 12.290180, 11.793243,
    11.242660, 10.587401, 9.877580,  9.139919,  8.395541,  7.440598,  5.378937,
    2.399050,  3.520786,  5.106620,  7.264754,  9.371436,  11.003332, 11.993763,
    12.744840, 13.490990, 13.789592, 13.534294, 13.177114, 13.619523, 14.064851,
    14.256814, 14.270717, 14.237335, 14.142181, 14.034318, 13.874064, 13.695080,
    13.418967, 13.086812, 12.711691, 12.290180, 11.793243, 11.242660, 10.587401,
    9.877580,  9.139919,  8.395541,  7.526351,  6.687473,  5.676038,  3.942031,
    2.408970,  3.439918,  4.908654,  7.295600,  9.493708,  11.047767, 11.989073,
    12.626884, 13.235572, 13.464634, 12.932467, 12.230970, 12.436694, 12.698073,
    13.310701, 13.831051, 13.979672, 13.871570, 13.695592, 13.419008, 13.086812,
    12.711691, 12.290180, 11.793243, 11.242660, 10.587401, 9.877580,  9.139919,
    8.395541,  7.526351,  6.687473,  5.736069,  4.856927,  3.990997,  2.721462,
    2.350457,  3.378298,  4.892127,  7.324143,  9.519243,  11.054103, 12.029321,
    12.733103, 13.247763, 13.439132, 13.279063, 12.982729, 12.404525, 11.297051,
    11.789822, 12.902933, 13.456238, 13.406731, 13.092388, 12.710679, 12.289666,
    11.793210, 11.242660, 10.587401, 9.877580,  9.139919,  8.395541,  7.526351,
    6.687473,  5.736069,  4.856927,  4.034050,  3.378202,  2.834415,  1.969789,
    2.396832,  3.423453,  4.922082,  7.393900,  9.599680,  11.158470, 12.245940,
    12.934559, 13.324933, 13.509052, 13.497501, 13.356077, 13.015855, 12.269659,
    12.052673, 12.205783, 12.507289, 12.598712, 12.264833, 11.766014, 11.229728,
    10.585106, 9.877662,  9.139936,  8.395541,  7.526351,  6.687473,  5.736069,
    4.856927,  4.034050,  3.378202,  2.868595,  2.494118,  2.240219,  1.663727,
    2.343165,  3.383423,  4.954288,  7.468584,  9.718328,  11.331646, 12.411226,
    13.000297, 13.276785, 13.310868, 13.139850, 12.928071, 12.696874, 12.221448,
    11.690038, 11.127237, 10.945462, 11.055842, 10.886695, 10.418341, 9.760101,
    9.104807,  8.394342,  7.527638,  6.688210,  5.736124,  4.856927,  4.033962,
    3.376849,  2.864250,  2.489774,  2.269461,  2.134461,  2.010454,  1.539838,
    2.149440,  3.202890,  4.912555,  7.499374,  9.789360,  11.355115, 12.279560,
    12.627570, 12.613710, 12.381159, 12.000994, 11.541156, 11.063136, 10.495633,
    9.876339,  9.106561,  8.554012,  8.429198,  8.500408,  8.406931,  7.887574,
    7.339313,  6.656567,  5.757767,  4.872202,  4.033921,  3.373438,  2.861690,
    2.470205,  2.206939,  2.071938,  2.018345,  1.959884,  1.821455,  1.382649,
    1.712478,  2.887924,  4.817530,  7.388206,  9.504125,  10.759708, 11.324473,
    11.209867, 10.721015, 10.101397, 9.460077,  8.803309,  8.197636,  7.682507,
    7.158105,  6.551958,  6.046410,  5.769328,  5.718300,  5.743889,  5.570121,
    5.279368,  4.775826,  4.153023,  3.462642,  2.847365,  2.418215,  2.183597,
    2.022818,  1.866233,  1.807803,  1.800176,  1.781400,  1.710956,  1.317016,
    1.514737,  2.850764,  4.677608,  6.843667,  8.357078,  8.941536,  8.973117,
    8.402416,  7.661619,  6.986591,  6.327657,  5.724844,  5.232076,  4.841967,
    4.512318,  4.218085,  3.968328,  3.781122,  3.652350,  3.592157,  3.485109,
    3.388087,  3.242951,  2.987120,  2.529736,  2.084546,  1.849609,  1.783321,
    1.765497,  1.724375,  1.708476,  1.713392,  1.715595,  1.687100,  1.311407,
    1.474474,  2.541837,  4.070860,  5.423715,  5.977272,  5.994783,  5.799174,
    5.250738,  4.723447,  4.300767,  3.923553,  3.562443,  3.221338,  2.938382,
    2.733331,  2.584102,  2.473217,  2.384763,  2.276152,  2.165501,  2.113679,
    2.062540,  1.996723,  1.877760,  1.664670,  1.476382,  1.435178,  1.456002,
    1.493136,  1.507210,  1.541560,  1.631957,  1.681868,  1.664835,  1.295367,
    1.415633,  2.191884,  3.189803,  3.722179,  3.721101,  3.507909,  3.269287,
    3.008462,  2.780443,  2.537032,  2.350791,  2.175407,  1.999644,  1.832543,
    1.681843,  1.567108,  1.491283,  1.437136,  1.385755,  1.328447,  1.277267,
    1.205540,  1.169675,  1.111642,  1.038942,  0.985774,  0.943086,  0.935795,
    0.948222,  0.952437,  1.029996,  1.247877,  1.369139,  1.369078,  1.065790,
    1.482158,  1.987374,  2.372558,  2.413073,  2.240173,  2.062185,  1.894405,
    1.737246,  1.599939,  1.481328,  1.389379,  1.281469,  1.195932,  1.115465,
    1.035558,  0.962719,  0.897782,  0.848492,  0.819022,  0.796491,  0.766239,
    0.726512,  0.697314,  0.660666,  0.644613,  0.634202,  0.560377,  0.525828,
    0.518590,  0.502100,  0.516968,  0.601150,  0.650649,  0.649290,  0.505403,
    1.618021,  1.849147,  1.765419,  1.575349,  1.360237,  1.242913,  1.152876,
    1.045707,  0.946684,  0.869254,  0.804963,  0.747487,  0.703272,  0.652710,
    0.613435,  0.577369,  0.541679,  0.508778,  0.479818,  0.458128,  0.444150,
    0.430717,  0.411765,  0.391253,  0.396149,  0.413876,  0.361522,  0.332692,
    0.328187,  0.320404,  0.314518,  0.318412,  0.322341,  0.318292,  0.247628,
    1.333242,  1.475332,  1.265179,  1.051776,  0.883980,  0.779891,  0.709055,
    0.644228,  0.584546,  0.528786,  0.478774,  0.439313,  0.407237,  0.379431,
    0.358196,  0.334206,  0.316247,  0.300246,  0.282611,  0.266699,  0.253154,
    0.242524,  0.233198,  0.223931,  0.236887,  0.265337,  0.222827,  0.197187,
    0.195122,  0.194569,  0.194251,  0.197825,  0.206959,  0.207512,  0.161621,
    0.920728,  1.028569,  0.876210,  0.733074,  0.620275,  0.533376,  0.469060,
    0.414339,  0.367299,  0.328733,  0.297061,  0.268816,  0.243679,  0.223612,
    0.207085,  0.192853,  0.182532,  0.171389,  0.162139,  0.154293,  0.146011,
    0.138206,  0.130593,  0.125114,  0.142748,  0.174177,  0.132188,  0.106477,
    0.104571,  0.104567,  0.105149,  0.114165,  0.136357,  0.143278,  0.111946,
    0.632495,  0.717677,  0.619933,  0.528229,  0.454409,  0.393433,  0.345540,
    0.295280,  0.244137,  0.209002,  0.187819,  0.171325,  0.153958,  0.137587,
    0.123852,  0.113304,  0.104915,  0.097878,  0.092612,  0.087242,  0.082770,
    0.078765,  0.074396,  0.071723,  0.088001,  0.114869,  0.079149,  0.057311,
    0.055692,  0.055638,  0.055131,  0.057020,  0.065397,  0.068049,  0.053137,
    0.448957,  0.517272,  0.454397,  0.393236,  0.344856,  0.295183,  0.245383,
    0.204067,  0.169322,  0.144875,  0.127127,  0.114301,  0.099387,  0.087020,
    0.077906,  0.069657,  0.062609,  0.057371,  0.053216,  0.049741,  0.047056,
    0.044351,  0.042074,  0.040953,  0.051115,  0.067684,  0.045860,  0.032492,
    0.031500,  0.031359,  0.029479,  0.025179,  0.023943,  0.023708,  0.018456,
    0.330989,  0.385313,  0.344856,  0.295183,  0.245389,  0.204158,  0.169494,
    0.144019,  0.123509,  0.107783,  0.090185,  0.076478,  0.065002,  0.056152,
    0.050418,  0.044764,  0.039350,  0.035224,  0.031879,  0.029241,  0.027115,
    0.025335,  0.023945,  0.022904,  0.025821,  0.030782,  0.022585,  0.017775,
    0.017411,  0.017193,  0.016023,  0.013988,  0.013207,  0.012954,  0.010076,
    0.253451,  0.289310,  0.245389,  0.204158,  0.169494,  0.144019,  0.123505,
    0.107729,  0.089963,  0.075795,  0.063263,  0.053069,  0.045471,  0.039093,
    0.034681,  0.030165,  0.025758,  0.022778,  0.020286,  0.018049,  0.016250,
    0.014890,  0.013780,  0.012902,  0.012800,  0.012986,  0.011189,  0.010393,
    0.010358,  0.009957,  0.008927,  0.008396,  0.008315,  0.008185,  0.006367,
    0.175383,  0.199661,  0.169494,  0.144019,  0.123505,  0.107729,  0.089963,
    0.075795,  0.063262,  0.053045,  0.045365,  0.038742,  0.033758,  0.028787,
    0.024120,  0.020381,  0.017402,  0.015468,  0.013450,  0.011657,  0.010278,
    0.009131,  0.008185,  0.007510,  0.007019,  0.006621,  0.006436,  0.006945,
    0.007201,  0.007055,  0.006654,  0.006490,  0.006489,  0.006416,  0.004999,
    0.121132,  0.140913,  0.123505,  0.107729,  0.089963,  0.075795,  0.063262,
    0.053045,  0.045365,  0.038742,  0.033757,  0.028774,  0.024065,  0.020223,
    0.016995,  0.014691,  0.012372,  0.010658,  0.009267,  0.007893,  0.006702,
    0.005832,  0.005173,  0.004681,  0.004316,  0.004069,  0.004263,  0.005233,
    0.005658,  0.005677,  0.005649,  0.005648,  0.005804,  0.006089,  0.004856,
    0.089342,  0.105529,  0.089963,  0.075795,  0.063262,  0.053045,  0.045365,
    0.038742,  0.033757,  0.028774,  0.024065,  0.020223,  0.016994,  0.014684,
    0.012344,  0.010568,  0.009054,  0.007571,  0.006284,  0.005398,  0.004547,
    0.003971,  0.003589,  0.003302,  0.003158,  0.003171,  0.003333,  0.003718,
    0.003880,  0.003893,  0.003893,  0.003919,  0.004301,  0.005143,  0.004291,
    0.063790,  0.074102,  0.063262,  0.053045,  0.045365,  0.038742,  0.033757,
    0.028774,  0.024065,  0.020223,  0.016994,  0.014684,  0.012344,  0.010568,
    0.009054,  0.007568,  0.006269,  0.005342,  0.004433,  0.003912,  0.003521,
    0.003136,  0.002943,  0.002859,  0.002836,  0.002901,  0.003019,  0.003100,
    0.003153,  0.003171,  0.003175,  0.003187,  0.003344,  0.003666,  0.002971,
    0.045190,  0.051881,  0.045365,  0.038742,  0.033757,  0.028774,  0.024065,
    0.020223,  0.016994,  0.014684,  0.012344,  0.010568,  0.009054,  0.007568,
    0.006269,  0.005342,  0.004433,  0.003909,  0.003515,  0.003150,  0.002975,
    0.002851,  0.002801,  0.002793,  0.002793,  0.002811,  0.002854,  0.002909,
    0.003007,  0.003062,  0.003105,  0.003125,  0.003138,  0.003116,  0.002432,
    0.033019,  0.037950,  0.033757,  0.028774,  0.024065,  0.020223,  0.016994,
    0.014684,  0.012344,  0.010568,  0.009054,  0.007568,  0.006269,  0.005342,
    0.004433,  0.003909,  0.003515,  0.003150,  0.002975,  0.002852,  0.002804,
    0.002793,  0.002790,  0.002790,  0.002790,  0.002791,  0.002795,  0.002812,
    0.002854,  0.002910,  0.003023,  0.003101,  0.003121,  0.003077,  0.002395,
    0.024815,  0.028201,  0.024065,  0.020223,  0.016994,  0.014684,  0.012344,
    0.010568,  0.009054,  0.007568,  0.006269,  0.005342,  0.004433,  0.003909,
    0.003515,  0.003150,  0.002975,  0.002852,  0.002804,  0.002793,  0.002790,
    0.002790,  0.002790,  0.002790,  0.002790,  0.002790,  0.002790,  0.002791,
    0.002795,  0.002815,  0.002893,  0.003006,  0.003062,  0.003057,  0.002392,
    0.017158,  0.019639,  0.016856,  0.014570,  0.012255,  0.010491,  0.008991,
    0.007511,  0.006216,  0.005296,  0.004391,  0.003867,  0.003473,  0.003107,
    0.002933,  0.002810,  0.002762,  0.002751,  0.002748,  0.002747,  0.002747,
    0.002747,  0.002747,  0.002747,  0.002747,  0.002747,  0.002747,  0.002747,
    0.002748,  0.002750,  0.002769,  0.002812,  0.002868,  0.002934,  0.002338,
    0.010306,  0.012127,  0.010260,  0.008850,  0.007701,  0.006408,  0.005314,
    0.004477,  0.003622,  0.003209,  0.002860,  0.002497,  0.002323,  0.002200,
    0.002151,  0.002140,  0.002137,  0.002137,  0.002137,  0.002137,  0.002137,
    0.002137,  0.002137,  0.002137,  0.002137,  0.002137,  0.002137,  0.002137,
    0.002137,  0.002137,  0.002138,  0.002142,  0.002162,  0.002204,  0.001792,
  };
  const double x_start = -0.5;
  const double x_end = 16.5;
  const double x_step = 0.5;
  const double y_start = -15.5;
  const double y_end = 16.5;
  const double y_step = 1.0;
  const int stride = (int)rint((x_end - x_start) / x_step) + 1;
  (void)y_end;

  const double y = (yl - y_start) / y_step;
  const double x = (xm - x_start) / x_step;

  const int yi = (int)floor(y);
  const int xi = (int)floor(x);
  const double yo = y - yi;
  const double xo = x - xi;
  const double *prate = &interp_rgrid[(yi - 1) * stride + (xi - 1)];
  const double *pdist = &interp_dgrid[(yi - 1) * stride + (xi - 1)];
  *rate_f = interp_bicubic(prate, stride, xo, yo);
  *dist_f = interp_bicubic(pdist, stride, xo, yo);
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
  const AV1_COMMON *const cm = &cpi->common;
  const int scaled_idx = cpi->scaled_ref_idx[ref_frame - 1];
  const int ref_idx = get_ref_frame_buf_idx(cpi, ref_frame);
  return (scaled_idx != ref_idx && scaled_idx != INVALID_IDX)
             ? &cm->buffer_pool->frame_bufs[scaled_idx].buf
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

  rd->thresh_mult[THR_DC] += 1000;

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

  rd->thresh_mult[THR_PAETH] += 1000;

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

  rd->thresh_mult[THR_H_PRED] += 2000;
  rd->thresh_mult[THR_V_PRED] += 2000;
  rd->thresh_mult[THR_D135_PRED] += 2500;
  rd->thresh_mult[THR_D203_PRED] += 2500;
  rd->thresh_mult[THR_D157_PRED] += 2500;
  rd->thresh_mult[THR_D67_PRED] += 2500;
  rd->thresh_mult[THR_D113_PRED] += 2500;
  rd->thresh_mult[THR_D45_PRED] += 2500;

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
