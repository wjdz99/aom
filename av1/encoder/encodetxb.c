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

#include "av1/common/scan.h"
#include "av1/common/blockd.h"
#include "av1/common/idct.h"
#include "av1/common/pred_common.h"
#include "av1/encoder/bitstream.h"
#include "av1/encoder/encodeframe.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/subexp.h"
#include "av1/encoder/tokenize.h"

void av1_alloc_txb_buf(Av1Comp *cpi) {
#if 0
  Av1Common *cm = &cpi->common;
  int mi_block_size = 1 << MI_SIZE_LOG2;
  // TODO(angiebird): Make sure cm->subsampling_x/y is set correctly, and then
  // use precise buffer size according to cm->subsampling_x/y
  int pixel_stride = mi_block_size * cm->mi_cols;
  int pixel_height = mi_block_size * cm->mi_rows;
  int i;
  for (i = 0; i < MAX_MB_PLANE; ++i) {
    CHECK_MEM_ERROR(
        cm, cpi->tcoeff_buf[i],
        aom_malloc(sizeof(*cpi->tcoeff_buf[i]) * pixel_stride * pixel_height));
  }
#else
  (void)cpi;
#endif
}

void av1_free_txb_buf(Av1Comp *cpi) {
#if 0
  int i;
  for (i = 0; i < MAX_MB_PLANE; ++i) {
    aom_free(cpi->tcoeff_buf[i]);
  }
#else
  (void)cpi;
#endif
}

static void write_golomb(AomWriter *w, int level) {
  int x = level + 1;
  int i = x;
  int length = 0;

  while (i) {
    i >>= 1;
    ++length;
  }
  assert(length > 0);

  for (i = 0; i < length - 1; ++i) aom_write_bit(w, 0);

  for (i = length - 1; i >= 0; --i) aom_write_bit(w, (x >> i) & 0x01);
}

void av1_write_coeffs_txb(const Av1Common *const cm, Macroblockd *xd,
                          AomWriter *w, int block, int plane,
                          const TranLowT *tcoeff, uint16_t eob,
                          TXB_CTX *TxbCtx) {
  AomProb *nz_map;
  AomProb *eob_flag;
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  const PlaneType plane_type = get_plane_type(plane);
  const TxSize tx_size = get_tx_size(plane, xd);
  const TxType tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const ScanOrder *const scan_order =
      get_scan(cm, tx_size, tx_type, is_inter_block(mbmi));
  const int16_t *scan = scan_order->scan;
  int c;
  int is_nz;
  const int bwl = b_width_log2_lookup[txsize_to_bsize[tx_size]] + 2;
  const int seg_eob = tx_size_2d[tx_size];
  uint8_t txb_mask[32 * 32] = { 0 };
  uint16_t update_eob = 0;

  aom_write(w, eob == 0, cm->fc->txb_skip[tx_size][TxbCtx->txb_skip_ctx]);

  if (eob == 0) return;
#if CONFIG_TXK_SEL
  av1_write_tx_type(cm, xd, block, plane, w);
#endif

  nz_map = cm->fc->nz_map[tx_size][plane_type];
  eob_flag = cm->fc->eob_flag[tx_size][plane_type];

  for (c = 0; c < eob; ++c) {
    int coeff_ctx = get_nz_map_ctx(tcoeff, txb_mask, scan[c], bwl);
    int eob_ctx = get_eob_ctx(tcoeff, scan[c], bwl);

    TranLowT v = tcoeff[scan[c]];
    is_nz = (v != 0);

    if (c == seg_eob - 1) break;

    aom_write(w, is_nz, nz_map[coeff_ctx]);

    if (is_nz) {
      aom_write(w, c == (eob - 1), eob_flag[eob_ctx]);
    }
    txb_mask[scan[c]] = 1;
  }

  int i;
  for (i = 0; i < NUM_BASE_LEVELS; ++i) {
    AomProb *coeff_base = cm->fc->coeff_base[tx_size][plane_type][i];

    update_eob = 0;
    for (c = eob - 1; c >= 0; --c) {
      TranLowT v = tcoeff[scan[c]];
      TranLowT level = abs(v);
      int sign = (v < 0) ? 1 : 0;
      int ctx;

      if (level <= i) continue;

      ctx = get_base_ctx(tcoeff, scan[c], bwl, i + 1);

      if (level == i + 1) {
        aom_write(w, 1, coeff_base[ctx]);
        if (c == 0) {
          aom_write(w, sign, cm->fc->dc_sign[plane_type][TxbCtx->dc_sign_ctx]);
        } else {
          aom_write_bit(w, sign);
        }
        continue;
      }
      aom_write(w, 0, coeff_base[ctx]);
      update_eob = AOMMAX(update_eob, c);
    }
  }

  for (c = update_eob; c >= 0; --c) {
    TranLowT v = tcoeff[scan[c]];
    TranLowT level = abs(v);
    int sign = (v < 0) ? 1 : 0;
    int idx;
    int ctx;

    if (level <= NUM_BASE_LEVELS) continue;

    if (c == 0) {
      aom_write(w, sign, cm->fc->dc_sign[plane_type][TxbCtx->dc_sign_ctx]);
    } else {
      aom_write_bit(w, sign);
    }

    // level is above 1.
    ctx = get_level_ctx(tcoeff, scan[c], bwl);
    for (idx = 0; idx < COEFF_BASE_RANGE; ++idx) {
      if (level == (idx + 1 + NUM_BASE_LEVELS)) {
        aom_write(w, 1, cm->fc->coeff_lps[tx_size][plane_type][ctx]);
        break;
      }
      aom_write(w, 0, cm->fc->coeff_lps[tx_size][plane_type][ctx]);
    }
    if (idx < COEFF_BASE_RANGE) continue;

    // use 0-th order Golomb code to handle the residual level.
    write_golomb(w, level - COEFF_BASE_RANGE - 1 - NUM_BASE_LEVELS);
  }
}

void av1_write_coeffs_mb(const Av1Common *const cm, Macroblock *x, AomWriter *w,
                         int plane) {
  Macroblockd *xd = &x->e_mbd;
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  BlockSize bsize = mbmi->sb_type;
  struct MacroblockdPlane *pd = &xd->plane[plane];

#if CONFIG_CB4X4
  const BlockSize plane_bsize = get_plane_block_size(bsize, pd);
#else
  const BlockSize plane_bsize =
      get_plane_block_size(AOMMAX(bsize, BLOCK_8X8), pd);
#endif
  const int max_blocks_wide = max_block_wide(xd, plane_bsize, plane);
  const int max_blocks_high = max_block_high(xd, plane_bsize, plane);
  TxSize tx_size = get_tx_size(plane, xd);
  const int bkw = tx_size_wide_unit[tx_size];
  const int bkh = tx_size_high_unit[tx_size];
  const int step = tx_size_wide_unit[tx_size] * tx_size_high_unit[tx_size];
  int row, col;
  int block = 0;
  for (row = 0; row < max_blocks_high; row += bkh) {
    for (col = 0; col < max_blocks_wide; col += bkw) {
      TranLowT *tcoeff = BLOCK_OFFSET(x->mbmi_ext->tcoeff[plane], block);
      uint16_t eob = x->mbmi_ext->eobs[plane][block];
      TXB_CTX TxbCtx = { x->mbmi_ext->txb_skip_ctx[plane][block],
                         x->mbmi_ext->dc_sign_ctx[plane][block] };
      av1_write_coeffs_txb(cm, xd, w, block, plane, tcoeff, eob, &TxbCtx);
      block += step;
    }
  }
}

static INLINE void get_base_ctx_set(const TranLowT *tcoeffs,
                                    int c,  // raster order
                                    const int bwl,
                                    int ctx_set[NUM_BASE_LEVELS]) {
  const int row = c >> bwl;
  const int col = c - (row << bwl);
  const int stride = 1 << bwl;
  int mag[NUM_BASE_LEVELS] = { 0 };
  int idx;
  TranLowT abs_coeff;
  int i;

  for (idx = 0; idx < BASE_CONTEXT_POSITION_NUM; ++idx) {
    int ref_row = row + base_ref_offset[idx][0];
    int ref_col = col + base_ref_offset[idx][1];
    int pos = (ref_row << bwl) + ref_col;

    if (ref_row < 0 || ref_col < 0 || ref_row >= stride || ref_col >= stride)
      continue;

    abs_coeff = abs(tcoeffs[pos]);

    for (i = 0; i < NUM_BASE_LEVELS; ++i) {
      ctx_set[i] += abs_coeff > i;
      if (base_ref_offset[idx][0] >= 0 && base_ref_offset[idx][1] >= 0)
        mag[i] |= abs_coeff > (i + 1);
    }
  }

  for (i = 0; i < NUM_BASE_LEVELS; ++i) {
    ctx_set[i] = (ctx_set[i] + 1) >> 1;

    if (row == 0 && col == 0)
      ctx_set[i] = (ctx_set[i] << 1) + mag[i];
    else if (row == 0)
      ctx_set[i] = 8 + (ctx_set[i] << 1) + mag[i];
    else if (col == 0)
      ctx_set[i] = 18 + (ctx_set[i] << 1) + mag[i];
    else
      ctx_set[i] = 28 + (ctx_set[i] << 1) + mag[i];
  }
  return;
}

int av1_cost_coeffs_txb(const Av1Comp *const cpi, Macroblock *x, int plane,
                        int block, TXB_CTX *TxbCtx) {
  const Av1Common *const cm = &cpi->common;
  Macroblockd *const xd = &x->e_mbd;
  const TxSize tx_size = get_tx_size(plane, xd);
  const PlaneType plane_type = get_plane_type(plane);
  const TxType tx_type = get_tx_type(plane_type, xd, block, tx_size);
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  const struct MacroblockPlane *p = &x->plane[plane];
  const int eob = p->eobs[block];
  const TranLowT *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  int c, cost;
  const int seg_eob = AOMMIN(eob, tx_size_2d[tx_size] - 1);
  int txb_skip_ctx = TxbCtx->txb_skip_ctx;
  AomProb *nz_map = xd->fc->nz_map[tx_size][plane_type];

  const int bwl = b_width_log2_lookup[txsize_to_bsize[tx_size]] + 2;
  // txb_mask is only initialized for once here. After that, it will be set when
  // coding zero map and then reset when coding level 1 info.
  uint8_t txb_mask[32 * 32] = { 0 };
  AomProb(*coeff_base)[COEFF_BASE_CONTEXTS] =
      xd->fc->coeff_base[tx_size][plane_type];

  const ScanOrder *const scan_order =
      get_scan(cm, tx_size, tx_type, is_inter_block(mbmi));
  const int16_t *scan = scan_order->scan;

  cost = 0;

  if (eob == 0) {
    cost = av1_cost_bit(xd->fc->txb_skip[tx_size][txb_skip_ctx], 1);
    return cost;
  }

  cost = av1_cost_bit(xd->fc->txb_skip[tx_size][txb_skip_ctx], 0);

#if CONFIG_TXK_SEL
  cost += av1_tx_type_cost(cpi, xd, mbmi->sb_type, plane, tx_size, tx_type);
#endif

  for (c = 0; c < eob; ++c) {
    TranLowT v = qcoeff[scan[c]];
    int is_nz = (v != 0);
    int level = abs(v);

    if (c < seg_eob) {
      int coeff_ctx = get_nz_map_ctx(qcoeff, txb_mask, scan[c], bwl);
      cost += av1_cost_bit(nz_map[coeff_ctx], is_nz);
    }

    if (is_nz) {
      int ctx_ls[NUM_BASE_LEVELS] = { 0 };
      int sign = (v < 0) ? 1 : 0;

      // sign bit cost
      if (c == 0) {
        int dc_sign_ctx = TxbCtx->dc_sign_ctx;

        cost += av1_cost_bit(xd->fc->dc_sign[plane_type][dc_sign_ctx], sign);
      } else {
        cost += av1_cost_bit(128, sign);
      }

      get_base_ctx_set(qcoeff, scan[c], bwl, ctx_ls);

      int i;
      for (i = 0; i < NUM_BASE_LEVELS; ++i) {
        if (level <= i) continue;

        if (level == i + 1) {
          cost += av1_cost_bit(coeff_base[i][ctx_ls[i]], 1);
          continue;
        }
        cost += av1_cost_bit(coeff_base[i][ctx_ls[i]], 0);
      }

      if (level > NUM_BASE_LEVELS) {
        int idx;
        int ctx;

        ctx = get_level_ctx(qcoeff, scan[c], bwl);

        for (idx = 0; idx < COEFF_BASE_RANGE; ++idx) {
          if (level == (idx + 1 + NUM_BASE_LEVELS)) {
            cost +=
                av1_cost_bit(xd->fc->coeff_lps[tx_size][plane_type][ctx], 1);
            break;
          }
          cost += av1_cost_bit(xd->fc->coeff_lps[tx_size][plane_type][ctx], 0);
        }

        if (idx >= COEFF_BASE_RANGE) {
          // residual cost
          int r = level - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
          int ri = r;
          int length = 0;

          while (ri) {
            ri >>= 1;
            ++length;
          }

          for (ri = 0; ri < length - 1; ++ri) cost += av1_cost_bit(128, 0);

          for (ri = length - 1; ri >= 0; --ri)
            cost += av1_cost_bit(128, (r >> ri) & 0x01);
        }
      }

      if (c < seg_eob) {
        int eob_ctx = get_eob_ctx(qcoeff, scan[c], bwl);
        cost += av1_cost_bit(xd->fc->eob_flag[tx_size][plane_type][eob_ctx],
                             c == (eob - 1));
      }
    }

    txb_mask[scan[c]] = 1;
  }

  return cost;
}

typedef struct TxbParams {
  const Av1Comp *cpi;
  ThreadData *td;
  int rate;
} TxbParams;

int av1_get_txb_entropy_context(const TranLowT *qcoeff,
                                const ScanOrder *scan_order, int eob) {
  const int16_t *scan = scan_order->scan;
  int cul_level = 0;
  int c;
  for (c = 0; c < eob; ++c) {
    cul_level += abs(qcoeff[scan[c]]);
  }

  cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);
  set_dc_sign(&cul_level, qcoeff[0]);

  return cul_level;
}

static void update_txb_context(int plane, int block, int blk_row, int blk_col,
                               BlockSize plane_bsize, TxSize tx_size,
                               void *arg) {
  TxbParams *const args = arg;
  const Av1Comp *cpi = args->cpi;
  const Av1Common *cm = &cpi->common;
  ThreadData *const td = args->td;
  Macroblock *const x = &td->mb;
  Macroblockd *const xd = &x->e_mbd;
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  struct MacroblockPlane *p = &x->plane[plane];
  struct MacroblockdPlane *pd = &xd->plane[plane];
  const uint16_t eob = p->eobs[block];
  const TranLowT *qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  const PlaneType plane_type = pd->plane_type;
  const TxType tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const ScanOrder *const scan_order =
      get_scan(cm, tx_size, tx_type, is_inter_block(mbmi));
  (void)plane_bsize;

  int cul_level = av1_get_txb_entropy_context(qcoeff, scan_order, eob);
  av1_set_contexts(xd, pd, plane, tx_size, cul_level, blk_col, blk_row);
}

static void update_and_record_txb_context(int plane, int block, int blk_row,
                                          int blk_col, BlockSize plane_bsize,
                                          TxSize tx_size, void *arg) {
  TxbParams *const args = arg;
  const Av1Comp *cpi = args->cpi;
  const Av1Common *cm = &cpi->common;
  ThreadData *const td = args->td;
  Macroblock *const x = &td->mb;
  Macroblockd *const xd = &x->e_mbd;
  struct MacroblockPlane *p = &x->plane[plane];
  struct MacroblockdPlane *pd = &xd->plane[plane];
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  int eob = p->eobs[block], update_eob = 0;
  const PlaneType plane_type = pd->plane_type;
  const TranLowT *qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  TranLowT *tcoeff = BLOCK_OFFSET(x->mbmi_ext->tcoeff[plane], block);
  const int segment_id = mbmi->segment_id;
  const TxType tx_type = get_tx_type(plane_type, xd, block, tx_size);
  const ScanOrder *const scan_order =
      get_scan(cm, tx_size, tx_type, is_inter_block(mbmi));
  const int16_t *scan = scan_order->scan;
  const int seg_eob = get_tx_eob(&cpi->common.seg, segment_id, tx_size);
  int c, i;
  TXB_CTX TxbCtx;
  get_txb_ctx(plane_bsize, tx_size, plane, pd->above_context + blk_col,
              pd->left_context + blk_row, &TxbCtx);
  const int bwl = b_width_log2_lookup[txsize_to_bsize[tx_size]] + 2;
  int cul_level = 0;
  unsigned int(*nz_map_count)[SIG_COEF_CONTEXTS][2];
  uint8_t txb_mask[32 * 32] = { 0 };

  nz_map_count = &td->counts->nz_map[tx_size][plane_type];

  memcpy(tcoeff, qcoeff, sizeof(*tcoeff) * seg_eob);

  ++td->counts->txb_skip[tx_size][TxbCtx.txb_skip_ctx][eob == 0];
  x->mbmi_ext->txb_skip_ctx[plane][block] = TxbCtx.txb_skip_ctx;

  x->mbmi_ext->eobs[plane][block] = eob;

  if (eob == 0) {
    av1_set_contexts(xd, pd, plane, tx_size, 0, blk_col, blk_row);
    return;
  }

#if CONFIG_TXK_SEL
  av1_update_tx_type_count(cm, xd, block, plane, mbmi->sb_type, tx_size,
                           td->counts);
#endif

  for (c = 0; c < eob; ++c) {
    TranLowT v = qcoeff[scan[c]];
    int is_nz = (v != 0);
    int coeff_ctx = get_nz_map_ctx(tcoeff, txb_mask, scan[c], bwl);
    int eob_ctx = get_eob_ctx(tcoeff, scan[c], bwl);

    if (c == seg_eob - 1) break;

    ++(*nz_map_count)[coeff_ctx][is_nz];

    if (is_nz) {
      ++td->counts->eob_flag[tx_size][plane_type][eob_ctx][c == (eob - 1)];
    }
    txb_mask[scan[c]] = 1;
  }

  // Reverse process order to handle coefficient level and sign.
  for (i = 0; i < NUM_BASE_LEVELS; ++i) {
    update_eob = 0;
    for (c = eob - 1; c >= 0; --c) {
      TranLowT v = qcoeff[scan[c]];
      TranLowT level = abs(v);
      int ctx;

      if (level <= i) continue;

      ctx = get_base_ctx(tcoeff, scan[c], bwl, i + 1);

      if (level == i + 1) {
        ++td->counts->coeff_base[tx_size][plane_type][i][ctx][1];
        if (c == 0) {
          int dc_sign_ctx = TxbCtx.dc_sign_ctx;

          ++td->counts->dc_sign[plane_type][dc_sign_ctx][v < 0];
          x->mbmi_ext->dc_sign_ctx[plane][block] = dc_sign_ctx;
        }
        cul_level += level;
        continue;
      }
      ++td->counts->coeff_base[tx_size][plane_type][i][ctx][0];
      update_eob = AOMMAX(update_eob, c);
    }
  }

  for (c = update_eob; c >= 0; --c) {
    TranLowT v = qcoeff[scan[c]];
    TranLowT level = abs(v);
    int idx;
    int ctx;

    if (level <= NUM_BASE_LEVELS) continue;

    cul_level += level;
    if (c == 0) {
      int dc_sign_ctx = TxbCtx.dc_sign_ctx;

      ++td->counts->dc_sign[plane_type][dc_sign_ctx][v < 0];
      x->mbmi_ext->dc_sign_ctx[plane][block] = dc_sign_ctx;
    }

    // level is above 1.
    ctx = get_level_ctx(tcoeff, scan[c], bwl);
    for (idx = 0; idx < COEFF_BASE_RANGE; ++idx) {
      if (level == (idx + 1 + NUM_BASE_LEVELS)) {
        ++td->counts->coeff_lps[tx_size][plane_type][ctx][1];
        break;
      }
      ++td->counts->coeff_lps[tx_size][plane_type][ctx][0];
    }
    if (idx < COEFF_BASE_RANGE) continue;

    // use 0-th order Golomb code to handle the residual level.
  }

  cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);

  // DC value
  set_dc_sign(&cul_level, tcoeff[0]);
  av1_set_contexts(xd, pd, plane, tx_size, cul_level, blk_col, blk_row);

#if CONFIG_ADAPT_SCAN
  // Since dqcoeff is not available here, we pass qcoeff into
  // av1_update_scan_count_facade(). The update behavior should be the same
  // because av1_update_scan_count_facade() only cares if coefficients are zero
  // or not.
  av1_update_scan_count_facade((Av1Common *)cm, td->counts, tx_size, tx_type,
                               qcoeff, eob);
#endif
}

void av1_update_txb_context(const Av1Comp *cpi, ThreadData *td, RunType dry_run,
                            BlockSize bsize, int *rate, int mi_row,
                            int mi_col) {
  const Av1Common *const cm = &cpi->common;
  Macroblock *const x = &td->mb;
  Macroblockd *const xd = &x->e_mbd;
  MbModeInfo *const mbmi = &xd->mi[0]->mbmi;
  const int ctx = av1_get_skip_context(xd);
  const int skip_inc =
      !segfeature_active(&cm->seg, mbmi->segment_id, SEG_LVL_SKIP);
  struct TxbParams arg = { cpi, td, 0 };
  (void)rate;
  (void)mi_row;
  (void)mi_col;
  if (mbmi->skip) {
    if (!dry_run) td->counts->skip[ctx][1] += skip_inc;
    reset_skip_context(xd, bsize);
    return;
  }

  if (!dry_run) {
    td->counts->skip[ctx][0] += skip_inc;
    av1_foreach_transformed_block(xd, bsize, mi_row, mi_col,
                                  update_and_record_txb_context, &arg);
  } else if (dry_run == DRY_RUN_NORMAL) {
    av1_foreach_transformed_block(xd, bsize, mi_row, mi_col, update_txb_context,
                                  &arg);
  } else {
    printf("DRY_RUN_COSTCOEFFS is not supported yet\n");
    assert(0);
  }
}

static void find_new_prob(unsigned int *branch_cnt, AomProb *oldp, int *savings,
                          int *update, AomWriter *const bc) {
  const AomProb upd = DIFF_UPDATE_PROB;
  int u = 0;
  AomProb newp = get_binary_prob(branch_cnt[0], branch_cnt[1]);
  int s = av1_prob_diff_update_savings_search(branch_cnt, *oldp, &newp, upd, 1);

  if (s > 0 && newp != *oldp) u = 1;

  if (u)
    *savings += s - (int)(av1_cost_zero(upd));  // TODO(jingning): 1?
  else
    *savings -= (int)(av1_cost_zero(upd));

  if (update) {
    ++update[u];
    return;
  }

  aom_write(bc, u, upd);
  if (u) {
    /* send/use new probability */
    av1_write_prob_diff_update(bc, newp, *oldp);
    *oldp = newp;
  }
}

static void write_txb_probs(AomWriter *const bc, Av1Comp *cpi, TxSize tx_size) {
  FrameContext *fc = cpi->common.fc;
  FrameCounts *counts = cpi->td.counts;
  int savings = 0;
  int update[2] = { 0, 0 };
  int plane, ctx, level;

  for (ctx = 0; ctx < TXB_SKIP_CONTEXTS; ++ctx) {
    find_new_prob(counts->txb_skip[tx_size][ctx], &fc->txb_skip[tx_size][ctx],
                  &savings, update, bc);
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < SIG_COEF_CONTEXTS; ++ctx) {
      find_new_prob(counts->nz_map[tx_size][plane][ctx],
                    &fc->nz_map[tx_size][plane][ctx], &savings, update, bc);
    }
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < EOB_COEF_CONTEXTS; ++ctx) {
      find_new_prob(counts->eob_flag[tx_size][plane][ctx],
                    &fc->eob_flag[tx_size][plane][ctx], &savings, update, bc);
    }
  }

  for (level = 0; level < NUM_BASE_LEVELS; ++level) {
    for (plane = 0; plane < PLANE_TYPES; ++plane) {
      for (ctx = 0; ctx < COEFF_BASE_CONTEXTS; ++ctx) {
        find_new_prob(counts->coeff_base[tx_size][plane][level][ctx],
                      &fc->coeff_base[tx_size][plane][level][ctx], &savings,
                      update, bc);
      }
    }
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < LEVEL_CONTEXTS; ++ctx) {
      find_new_prob(counts->coeff_lps[tx_size][plane][ctx],
                    &fc->coeff_lps[tx_size][plane][ctx], &savings, update, bc);
    }
  }

  // Decide if to update the model for this tx_size
  if (update[1] == 0 || savings < 0) {
    aom_write_bit(bc, 0);
    return;
  }
  aom_write_bit(bc, 1);

  for (ctx = 0; ctx < TXB_SKIP_CONTEXTS; ++ctx) {
    find_new_prob(counts->txb_skip[tx_size][ctx], &fc->txb_skip[tx_size][ctx],
                  &savings, NULL, bc);
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < SIG_COEF_CONTEXTS; ++ctx) {
      find_new_prob(counts->nz_map[tx_size][plane][ctx],
                    &fc->nz_map[tx_size][plane][ctx], &savings, NULL, bc);
    }
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < EOB_COEF_CONTEXTS; ++ctx) {
      find_new_prob(counts->eob_flag[tx_size][plane][ctx],
                    &fc->eob_flag[tx_size][plane][ctx], &savings, NULL, bc);
    }
  }

  for (level = 0; level < NUM_BASE_LEVELS; ++level) {
    for (plane = 0; plane < PLANE_TYPES; ++plane) {
      for (ctx = 0; ctx < COEFF_BASE_CONTEXTS; ++ctx) {
        find_new_prob(counts->coeff_base[tx_size][plane][level][ctx],
                      &fc->coeff_base[tx_size][plane][level][ctx], &savings,
                      NULL, bc);
      }
    }
  }

  for (plane = 0; plane < PLANE_TYPES; ++plane) {
    for (ctx = 0; ctx < LEVEL_CONTEXTS; ++ctx) {
      find_new_prob(counts->coeff_lps[tx_size][plane][ctx],
                    &fc->coeff_lps[tx_size][plane][ctx], &savings, NULL, bc);
    }
  }
}

void av1_write_txb_probs(Av1Comp *cpi, AomWriter *w) {
  const TxMode tx_mode = cpi->common.tx_mode;
  const TxSize max_tx_size = tx_mode_to_biggest_tx_size[tx_mode];
  TxSize tx_size;
  int ctx, plane;

  for (plane = 0; plane < PLANE_TYPES; ++plane)
    for (ctx = 0; ctx < DC_SIGN_CONTEXTS; ++ctx)
      av1_cond_prob_diff_update(w, &cpi->common.fc->dc_sign[plane][ctx],
                                cpi->td.counts->dc_sign[plane][ctx], 1);

  for (tx_size = TX_4X4; tx_size <= max_tx_size; ++tx_size)
    write_txb_probs(w, cpi, tx_size);
}

#if CONFIG_TXK_SEL
int64_t av1_search_txk_type(const Av1Comp *cpi, Macroblock *x, int plane,
                            int block, int blk_row, int blk_col,
                            BlockSize plane_bsize, TxSize tx_size,
                            const EntropyContext *a, const EntropyContext *l,
                            int use_fast_coef_costing, RdStats *rd_stats) {
  const Av1Common *cm = &cpi->common;
  Macroblockd *xd = &x->e_mbd;
  MbModeInfo *mbmi = &xd->mi[0]->mbmi;
  TxType txk_start = DCT_DCT;
  TxType txk_end = TX_TYPES - 1;
  TxType best_tx_type = txk_start;
  int64_t best_rd = INT64_MAX;
  const int coeff_ctx = combine_entropy_contexts(*a, *l);
  TxType tx_type;
  for (tx_type = txk_start; tx_type <= txk_end; ++tx_type) {
    if (plane == 0) mbmi->txk_type[block] = tx_type;
    TxType ref_tx_type = get_tx_type(get_plane_type(plane), xd, block, tx_size);
    if (tx_type != ref_tx_type) {
      // use get_tx_type() to check if the tx_type is valid for the current mode
      // if it's not, we skip it here.
      continue;
    }
    RdStats this_rd_stats;
    av1_invalid_rd_stats(&this_rd_stats);
    av1_xform_quant(cm, x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                    coeff_ctx, AV1_XFORM_QUANT_FP);
    av1_optimize_b(cm, x, plane, block, tx_size, coeff_ctx);
    av1_dist_block(cpi, x, plane, plane_bsize, block, blk_row, blk_col, tx_size,
                   &this_rd_stats.dist, &this_rd_stats.sse,
                   OUTPUT_HAS_PREDICTED_PIXELS);
    const ScanOrder *scan_order =
        get_scan(cm, tx_size, tx_type, is_inter_block(mbmi));
    this_rd_stats.rate = av1_cost_coeffs(
        cpi, x, plane, block, tx_size, scan_order, a, l, use_fast_coef_costing);
    int rd =
        RDCOST(x->rdmult, x->rddiv, this_rd_stats.rate, this_rd_stats.dist);
    if (rd < best_rd) {
      best_rd = rd;
      *rd_stats = this_rd_stats;
      best_tx_type = tx_type;
    }
  }
  if (plane == 0) mbmi->txk_type[block] = best_tx_type;
  // TODO(angiebird): Instead of re-call av1_xform_quant and av1_optimize_b,
  // copy the best result in the above tx_type search for loop
  av1_xform_quant(cm, x, plane, block, blk_row, blk_col, plane_bsize, tx_size,
                  coeff_ctx, AV1_XFORM_QUANT_FP);
  av1_optimize_b(cm, x, plane, block, tx_size, coeff_ctx);
  if (!is_inter_block(mbmi)) {
    // intra mode needs decoded result such that the next transform block
    // can use it for prediction.
    av1_inverse_transform_block_facade(xd, plane, block, blk_row, blk_col,
                                       x->plane[plane].eobs[block]);
  }
  return best_rd;
}
#endif  // CONFIG_TXK_SEL
