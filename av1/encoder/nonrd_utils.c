/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "config/aom_dsp_rtcd.h"

#include "av1/common/reconinter.h"

#include "av1/encoder/encodemv.h"
#include "av1/encoder/nonrd_opt.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/nonrd_utils.h"

static const SCAN_ORDER av1_fast_idtx_scan_order_16x16 = {
  av1_fast_idtx_scan_16x16, av1_fast_idtx_iscan_16x16
};

#define DECLARE_BLOCK_YRD_BUFFERS()                      \
  DECLARE_ALIGNED(64, tran_low_t, dqcoeff_buf[16 * 16]); \
  DECLARE_ALIGNED(64, tran_low_t, qcoeff_buf[16 * 16]);  \
  DECLARE_ALIGNED(64, tran_low_t, coeff_buf[16 * 16]);   \
  uint16_t eob[1];

#define DECLARE_BLOCK_YRD_VARS()                                          \
  /* When is_tx_8x8_dual_applicable is true, we compute the txfm for the  \
   * entire bsize and write macroblock_plane::coeff. So low_coeff is kept \
   * as a non-const so we can reassign it to macroblock_plane::coeff. */  \
  int16_t *low_coeff = (int16_t *)coeff_buf;                              \
  int16_t *const low_qcoeff = (int16_t *)qcoeff_buf;                      \
  int16_t *const low_dqcoeff = (int16_t *)dqcoeff_buf;                    \
  const int diff_stride = bw;

#define DECLARE_LOOP_VARS_BLOCK_YRD() \
  const int16_t *src_diff = &p->src_diff[(r * diff_stride + c) << 2];

static AOM_FORCE_INLINE void update_yrd_loop_vars(
    MACROBLOCK *x, int *skippable, int step, int ncoeffs,
    int16_t *const low_coeff, int16_t *const low_qcoeff,
    int16_t *const low_dqcoeff, RD_STATS *this_rdc, int *eob_cost,
    int tx_blk_id) {
  const int is_txfm_skip = (ncoeffs == 0);
  *skippable &= is_txfm_skip;
  x->txfm_search_info.blk_skip[tx_blk_id] = is_txfm_skip;
  *eob_cost += get_msb(ncoeffs + 1);
  if (ncoeffs == 1)
    this_rdc->rate += (int)abs(low_qcoeff[0]);
  else if (ncoeffs > 1)
    this_rdc->rate += aom_satd_lp(low_qcoeff, step << 4);

  this_rdc->dist += av1_block_error_lp(low_coeff, low_dqcoeff, step << 4) >> 2;
}

static INLINE void aom_process_hadamard_lp_8x16(MACROBLOCK *x,
                                                int max_blocks_high,
                                                int max_blocks_wide,
                                                int num_4x4_w, int step,
                                                int block_step) {
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  const int bw = 4 * num_4x4_w;
  const int num_4x4 = AOMMIN(num_4x4_w, max_blocks_wide);
  int block = 0;

  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0; c < num_4x4; c += 2 * block_step) {
      const int16_t *src_diff = &p->src_diff[(r * bw + c) << 2];
      int16_t *low_coeff = (int16_t *)p->coeff + BLOCK_OFFSET(block);
      aom_hadamard_lp_8x8_dual(src_diff, (ptrdiff_t)bw, low_coeff);
      block += 2 * step;
    }
  }
}

#if CONFIG_AV1_HIGHBITDEPTH
#define DECLARE_BLOCK_YRD_HBD_VARS()     \
  tran_low_t *const coeff = coeff_buf;   \
  tran_low_t *const qcoeff = qcoeff_buf; \
  tran_low_t *const dqcoeff = dqcoeff_buf;

static AOM_FORCE_INLINE void update_yrd_loop_vars_hbd(
    MACROBLOCK *x, int *skippable, int step, int ncoeffs,
    tran_low_t *const coeff, tran_low_t *const qcoeff,
    tran_low_t *const dqcoeff, RD_STATS *this_rdc, int *eob_cost,
    int tx_blk_id) {
  const int is_txfm_skip = (ncoeffs == 0);
  *skippable &= is_txfm_skip;
  x->txfm_search_info.blk_skip[tx_blk_id] = is_txfm_skip;
  *eob_cost += get_msb(ncoeffs + 1);

  int64_t dummy;
  if (ncoeffs == 1)
    this_rdc->rate += (int)abs(qcoeff[0]);
  else if (ncoeffs > 1)
    this_rdc->rate += aom_satd(qcoeff, step << 4);

  this_rdc->dist += av1_block_error(coeff, dqcoeff, step << 4, &dummy) >> 2;
}
#endif

/*!\brief Calculates RD Cost using Hadamard transform.
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Calculates RD Cost using Hadamard transform. For low bit depth this function
 * uses low-precision set of functions (16-bit) and 32 bit for high bit depth
 * \param[in]    x              Pointer to structure holding all the data for
                                the current macroblock
 * \param[in]    this_rdc       Pointer to calculated RD Cost
 * \param[in]    skippable      Pointer to a flag indicating possible tx skip
 * \param[in]    bsize          Current block size
 * \param[in]    tx_size        Transform size
 * \param[in]    is_inter_mode  Flag to indicate inter mode
 *
 * \remark Nothing is returned. Instead, calculated RD cost is placed to
 * \c this_rdc. \c skippable flag is set if there is no non-zero quantized
 * coefficients for Hadamard transform
 */
void av1_block_yrd(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                   BLOCK_SIZE bsize, TX_SIZE tx_size, int is_inter_mode) {
  MACROBLOCKD *xd = &x->e_mbd;
  const struct macroblockd_plane *pd = &xd->plane[AOM_PLANE_Y];
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  assert(bsize < BLOCK_SIZES_ALL);
  const int num_4x4_w = mi_size_wide[bsize];
  const int num_4x4_h = mi_size_high[bsize];
  const int step = 1 << (tx_size << 1);
  const int block_step = (1 << tx_size);
  const int row_step = step * num_4x4_w >> tx_size;
  int block = 0;
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >> 5);
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 : xd->mb_to_bottom_edge >> 5);
  int eob_cost = 0;
  const int bw = 4 * num_4x4_w;
  const int bh = 4 * num_4x4_h;
  const int use_hbd = is_cur_buf_hbd(xd);
  int num_blk_skip_w = num_4x4_w;
  int sh_blk_skip = 0;
  if (is_inter_mode) {
    num_blk_skip_w = num_4x4_w >> 1;
    sh_blk_skip = 1;
  }

#if CONFIG_AV1_HIGHBITDEPTH
  if (use_hbd) {
    aom_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride);
  } else {
    aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                       pd->dst.buf, pd->dst.stride);
  }
#else
  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
#endif

  // Keep the intermediate value on the stack here. Writing directly to
  // skippable causes speed regression due to load-and-store issues in
  // update_yrd_loop_vars.
  int temp_skippable = 1;
  this_rdc->dist = 0;
  this_rdc->rate = 0;
  // For block sizes 8x16 or above, Hadamard txfm of two adjacent 8x8 blocks
  // can be done per function call. Hence the call of Hadamard txfm is
  // abstracted here for the specified cases.
  int is_tx_8x8_dual_applicable =
      (tx_size == TX_8X8 && block_size_wide[bsize] >= 16 &&
       block_size_high[bsize] >= 8);

#if CONFIG_AV1_HIGHBITDEPTH
  // As of now, dual implementation of hadamard txfm is available for low
  // bitdepth.
  if (use_hbd) is_tx_8x8_dual_applicable = 0;
#endif

  if (is_tx_8x8_dual_applicable) {
    aom_process_hadamard_lp_8x16(x, max_blocks_high, max_blocks_wide, num_4x4_w,
                                 step, block_step);
  }

  const SCAN_ORDER *const scan_order = &av1_scan_orders[tx_size][DCT_DCT];
  DECLARE_BLOCK_YRD_BUFFERS()
  DECLARE_BLOCK_YRD_VARS()
#if CONFIG_AV1_HIGHBITDEPTH
  DECLARE_BLOCK_YRD_HBD_VARS()
#else
  (void)use_hbd;
#endif

  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0, s = 0; c < max_blocks_wide; c += block_step, s += step) {
      DECLARE_LOOP_VARS_BLOCK_YRD()

      switch (tx_size) {
#if CONFIG_AV1_HIGHBITDEPTH
        case TX_16X16:
          if (use_hbd) {
            aom_hadamard_16x16(src_diff, diff_stride, coeff);
            av1_quantize_fp(coeff, 16 * 16, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob,
                            // default_scan_fp_16x16_transpose and
                            // av1_default_iscan_fp_16x16_transpose have to be
                            // used together.
                            default_scan_fp_16x16_transpose,
                            av1_default_iscan_fp_16x16_transpose);
          } else {
            aom_hadamard_lp_16x16(src_diff, diff_stride, low_coeff);
            av1_quantize_lp(low_coeff, 16 * 16, p->round_fp_QTX,
                            p->quant_fp_QTX, low_qcoeff, low_dqcoeff,
                            p->dequant_QTX, eob,
                            // default_scan_lp_16x16_transpose and
                            // av1_default_iscan_lp_16x16_transpose have to be
                            // used together.
                            default_scan_lp_16x16_transpose,
                            av1_default_iscan_lp_16x16_transpose);
          }
          break;
        case TX_8X8:
          if (use_hbd) {
            aom_hadamard_8x8(src_diff, diff_stride, coeff);
            av1_quantize_fp(
                coeff, 8 * 8, p->zbin_QTX, p->round_fp_QTX, p->quant_fp_QTX,
                p->quant_shift_QTX, qcoeff, dqcoeff, p->dequant_QTX, eob,
                default_scan_8x8_transpose, av1_default_iscan_8x8_transpose);
          } else {
            if (is_tx_8x8_dual_applicable) {
              // The coeffs are pre-computed for the whole block, so re-assign
              // low_coeff to the appropriate location.
              const int block_offset = BLOCK_OFFSET(block + s);
              low_coeff = (int16_t *)p->coeff + block_offset;
            } else {
              aom_hadamard_lp_8x8(src_diff, diff_stride, low_coeff);
            }
            av1_quantize_lp(
                low_coeff, 8 * 8, p->round_fp_QTX, p->quant_fp_QTX, low_qcoeff,
                low_dqcoeff, p->dequant_QTX, eob,
                // default_scan_8x8_transpose and
                // av1_default_iscan_8x8_transpose have to be used together.
                default_scan_8x8_transpose, av1_default_iscan_8x8_transpose);
          }
          break;
        default:
          assert(tx_size == TX_4X4);
          // In tx_size=4x4 case, aom_fdct4x4 and aom_fdct4x4_lp generate
          // normal coefficients order, so we don't need to change the scan
          // order here.
          if (use_hbd) {
            aom_fdct4x4(src_diff, coeff, diff_stride);
            av1_quantize_fp(coeff, 4 * 4, p->zbin_QTX, p->round_fp_QTX,
                            p->quant_fp_QTX, p->quant_shift_QTX, qcoeff,
                            dqcoeff, p->dequant_QTX, eob, scan_order->scan,
                            scan_order->iscan);
          } else {
            aom_fdct4x4_lp(src_diff, low_coeff, diff_stride);
            av1_quantize_lp(low_coeff, 4 * 4, p->round_fp_QTX, p->quant_fp_QTX,
                            low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                            scan_order->scan, scan_order->iscan);
          }
          break;
#else
        case TX_16X16:
          aom_hadamard_lp_16x16(src_diff, diff_stride, low_coeff);
          av1_quantize_lp(low_coeff, 16 * 16, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          default_scan_lp_16x16_transpose,
                          av1_default_iscan_lp_16x16_transpose);
          break;
        case TX_8X8:
          if (is_tx_8x8_dual_applicable) {
            // The coeffs are pre-computed for the whole block, so re-assign
            // low_coeff to the appropriate location.
            const int block_offset = BLOCK_OFFSET(block + s);
            low_coeff = (int16_t *)p->coeff + block_offset;
          } else {
            aom_hadamard_lp_8x8(src_diff, diff_stride, low_coeff);
          }
          av1_quantize_lp(low_coeff, 8 * 8, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          default_scan_8x8_transpose,
                          av1_default_iscan_8x8_transpose);
          break;
        default:
          aom_fdct4x4_lp(src_diff, low_coeff, diff_stride);
          av1_quantize_lp(low_coeff, 4 * 4, p->round_fp_QTX, p->quant_fp_QTX,
                          low_qcoeff, low_dqcoeff, p->dequant_QTX, eob,
                          scan_order->scan, scan_order->iscan);
          break;
#endif
      }
      assert(*eob <= 1024);
#if CONFIG_AV1_HIGHBITDEPTH
      if (use_hbd)
        update_yrd_loop_vars_hbd(x, &temp_skippable, step, *eob, coeff, qcoeff,
                                 dqcoeff, this_rdc, &eob_cost,
                                 (r * num_blk_skip_w + c) >> sh_blk_skip);
      else
#endif
        update_yrd_loop_vars(x, &temp_skippable, step, *eob, low_coeff,
                             low_qcoeff, low_dqcoeff, this_rdc, &eob_cost,
                             (r * num_blk_skip_w + c) >> sh_blk_skip);
    }
    block += row_step;
  }

  this_rdc->skip_txfm = *skippable = temp_skippable;
  if (this_rdc->sse < INT64_MAX) {
    this_rdc->sse = (this_rdc->sse << 6) >> 2;
    if (temp_skippable) {
      this_rdc->dist = 0;
      this_rdc->dist = this_rdc->sse;
      return;
    }
  }

  // If skippable is set, rate gets clobbered later.
  this_rdc->rate <<= (2 + AV1_PROB_COST_SHIFT);
  this_rdc->rate += (eob_cost << AV1_PROB_COST_SHIFT);
}

// Explicitly enumerate the cases so the compiler can generate SIMD for the
// function. According to the disassembler, gcc generates SSE codes for each of
// the possible block sizes. The hottest case is tx_width 16, which takes up
// about 8% of the self cycle of av1_nonrd_pick_inter_mode_sb. Since
// av1_nonrd_pick_inter_mode_sb takes up about 3% of total encoding time, the
// potential room of improvement for writing AVX2 optimization is only 3% * 8% =
// 0.24% of total encoding time.
static AOM_INLINE void scale_square_buf_vals(int16_t *dst, int tx_width,
                                             const int16_t *src,
                                             int src_stride) {
#define DO_SCALING                                                   \
  do {                                                               \
    for (int idy = 0; idy < tx_width; ++idy) {                       \
      for (int idx = 0; idx < tx_width; ++idx) {                     \
        dst[idy * tx_width + idx] = src[idy * src_stride + idx] * 8; \
      }                                                              \
    }                                                                \
  } while (0)

  if (tx_width == 4) {
    DO_SCALING;
  } else if (tx_width == 8) {
    DO_SCALING;
  } else if (tx_width == 16) {
    DO_SCALING;
  } else {
    assert(0);
  }

#undef DO_SCALING
}

/*!\brief Calculates RD Cost when the block uses Identity transform.
 * Note that thie function is only for low bit depth encoding, since it
 * is called in real-time mode for now, which sets high bit depth to 0:
 * -DCONFIG_AV1_HIGHBITDEPTH=0
 *
 * \ingroup nonrd_mode_search
 * \callgraph
 * \callergraph
 * Calculates RD Cost. For low bit depth this function
 * uses low-precision set of functions (16-bit) and 32 bit for high bit depth
 * \param[in]    x              Pointer to structure holding all the data for
                                the current macroblock
 * \param[in]    this_rdc       Pointer to calculated RD Cost
 * \param[in]    skippable      Pointer to a flag indicating possible tx skip
 * \param[in]    bsize          Current block size
 * \param[in]    tx_size        Transform size
 *
 * \remark Nothing is returned. Instead, calculated RD cost is placed to
 * \c this_rdc. \c skippable flag is set if all coefficients are zero.
 */
void av1_block_yrd_idtx(MACROBLOCK *x, RD_STATS *this_rdc, int *skippable,
                        BLOCK_SIZE bsize, TX_SIZE tx_size) {
  MACROBLOCKD *xd = &x->e_mbd;
  const struct macroblockd_plane *pd = &xd->plane[AOM_PLANE_Y];
  struct macroblock_plane *const p = &x->plane[AOM_PLANE_Y];
  assert(bsize < BLOCK_SIZES_ALL);
  const int num_4x4_w = mi_size_wide[bsize];
  const int num_4x4_h = mi_size_high[bsize];
  const int step = 1 << (tx_size << 1);
  const int block_step = (1 << tx_size);
  const int max_blocks_wide =
      num_4x4_w + (xd->mb_to_right_edge >= 0 ? 0 : xd->mb_to_right_edge >> 5);
  const int max_blocks_high =
      num_4x4_h + (xd->mb_to_bottom_edge >= 0 ? 0 : xd->mb_to_bottom_edge >> 5);
  int eob_cost = 0;
  const int bw = 4 * num_4x4_w;
  const int bh = 4 * num_4x4_h;
  const int num_blk_skip_w = num_4x4_w >> 1;
  const int sh_blk_skip = 1;
  // Keep the intermediate value on the stack here. Writing directly to
  // skippable causes speed regression due to load-and-store issues in
  // update_yrd_loop_vars.
  int temp_skippable = 1;
  int tx_wd = 0;
  const SCAN_ORDER *scan_order = NULL;
  switch (tx_size) {
    case TX_64X64:
      assert(0);  // Not implemented
      break;
    case TX_32X32:
      assert(0);  // Not used
      break;
    case TX_16X16:
      scan_order = &av1_fast_idtx_scan_order_16x16;
      tx_wd = 16;
      break;
    case TX_8X8:
      scan_order = &av1_fast_idtx_scan_order_8x8;
      tx_wd = 8;
      break;
    default:
      assert(tx_size == TX_4X4);
      scan_order = &av1_fast_idtx_scan_order_4x4;
      tx_wd = 4;
      break;
  }
  assert(scan_order != NULL);

  this_rdc->dist = 0;
  this_rdc->rate = 0;
  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
  // Keep track of the row and column of the blocks we use so that we know
  // if we are in the unrestricted motion border.
  DECLARE_BLOCK_YRD_BUFFERS()
  DECLARE_BLOCK_YRD_VARS()
  for (int r = 0; r < max_blocks_high; r += block_step) {
    for (int c = 0, s = 0; c < max_blocks_wide; c += block_step, s += step) {
      DECLARE_LOOP_VARS_BLOCK_YRD()
      scale_square_buf_vals(low_coeff, tx_wd, src_diff, diff_stride);
      av1_quantize_lp(low_coeff, tx_wd * tx_wd, p->round_fp_QTX,
                      p->quant_fp_QTX, low_qcoeff, low_dqcoeff, p->dequant_QTX,
                      eob, scan_order->scan, scan_order->iscan);
      assert(*eob <= 1024);
      update_yrd_loop_vars(x, &temp_skippable, step, *eob, low_coeff,
                           low_qcoeff, low_dqcoeff, this_rdc, &eob_cost,
                           (r * num_blk_skip_w + c) >> sh_blk_skip);
    }
  }
  this_rdc->skip_txfm = *skippable = temp_skippable;
  if (this_rdc->sse < INT64_MAX) {
    this_rdc->sse = (this_rdc->sse << 6) >> 2;
    if (temp_skippable) {
      this_rdc->dist = 0;
      this_rdc->dist = this_rdc->sse;
      return;
    }
  }
  // If skippable is set, rate gets clobbered later.
  this_rdc->rate <<= (2 + AV1_PROB_COST_SHIFT);
  this_rdc->rate += (eob_cost << AV1_PROB_COST_SHIFT);
}

#if CONFIG_INTERNAL_STATS
void av1_store_coding_context_nonrd(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx,
                                    int mode_index) {
#else
void av1_store_coding_context_nonrd(MACROBLOCK *x, PICK_MODE_CONTEXT *ctx) {
#endif  // CONFIG_INTERNAL_STATS
  MACROBLOCKD *const xd = &x->e_mbd;
  TxfmSearchInfo *txfm_info = &x->txfm_search_info;

  // Take a snapshot of the coding context so it can be
  // restored if we decide to encode this way
  ctx->rd_stats.skip_txfm = txfm_info->skip_txfm;

  ctx->skippable = txfm_info->skip_txfm;
#if CONFIG_INTERNAL_STATS
  ctx->best_mode_index = mode_index;
#endif  // CONFIG_INTERNAL_STATS
  ctx->mic = *xd->mi[0];
  ctx->skippable = txfm_info->skip_txfm;
  av1_copy_mbmi_ext_to_mbmi_ext_frame(&ctx->mbmi_ext_best, &x->mbmi_ext,
                                      av1_ref_frame_type(xd->mi[0]->ref_frame));
}

int64_t av1_model_rd_for_sb_uv(AV1_COMP *cpi, BLOCK_SIZE plane_bsize,
                               MACROBLOCK *x, MACROBLOCKD *xd,
                               RD_STATS *this_rdc, int start_plane,
                               int stop_plane) {
  // Note our transform coeffs are 8 times an orthogonal transform.
  // Hence quantizer step is also 8 times. To get effective quantizer
  // we need to divide by 8 before sending to modeling function.
  unsigned int sse;
  int rate;
  int64_t dist;
  int plane;
  int64_t tot_sse = 0;

  this_rdc->rate = 0;
  this_rdc->dist = 0;
  this_rdc->skip_txfm = 0;

  for (plane = start_plane; plane <= stop_plane; ++plane) {
    struct macroblock_plane *const p = &x->plane[plane];
    struct macroblockd_plane *const pd = &xd->plane[plane];
    const uint32_t dc_quant = p->dequant_QTX[0];
    const uint32_t ac_quant = p->dequant_QTX[1];
    const BLOCK_SIZE bs = plane_bsize;
    unsigned int var;
    if (!x->color_sensitivity[COLOR_SENS_IDX(plane)]) continue;

    var = cpi->ppi->fn_ptr[bs].vf(p->src.buf, p->src.stride, pd->dst.buf,
                                  pd->dst.stride, &sse);
    assert(sse >= var);
    tot_sse += sse;

    av1_model_rd_from_var_lapndz(sse - var, num_pels_log2_lookup[bs],
                                 dc_quant >> 3, &rate, &dist);

    this_rdc->rate += rate >> 1;
    this_rdc->dist += dist << 3;

    av1_model_rd_from_var_lapndz(var, num_pels_log2_lookup[bs], ac_quant >> 3,
                                 &rate, &dist);

    this_rdc->rate += rate;
    this_rdc->dist += dist << 4;
  }

  if (this_rdc->rate == 0) {
    this_rdc->skip_txfm = 1;
  }

  if (RDCOST(x->rdmult, this_rdc->rate, this_rdc->dist) >=
      RDCOST(x->rdmult, 0, tot_sse << 4)) {
    this_rdc->rate = 0;
    this_rdc->dist = tot_sse << 4;
    this_rdc->skip_txfm = 1;
  }

  return tot_sse;
}
