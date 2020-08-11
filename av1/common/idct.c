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

#include <math.h>

#include "config/aom_dsp_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom_ports/mem.h"
#include "av1/common/av1_inv_txfm1d_cfg.h"
#include "av1/common/av1_txfm.h"
#include "av1/common/blockd.h"
#include "av1/common/enums.h"
#include "av1/common/idct.h"
#if CONFIG_DSPL_RESIDUAL
#include "av1/common/resize.h"
#include "av1/common/scan.h"
#endif

int av1_get_tx_scale(const TX_SIZE tx_size) {
  const int pels = tx_size_2d[tx_size];
  // Largest possible pels is 4096 (64x64).
  return (pels > 256) + (pels > 1024);
}

// NOTE: The implementation of all inverses need to be aware of the fact
// that input and output could be the same buffer.

// idct
void av1_highbd_iwht4x4_add(const tran_low_t *input, uint8_t *dest, int stride,
                            int eob, int bd) {
  if (eob > 1)
    av1_highbd_iwht4x4_16_add(input, dest, stride, bd);
  else
    av1_highbd_iwht4x4_1_add(input, dest, stride, bd);
}

void av1_highbd_inv_txfm_add_4x4_c(const tran_low_t *input, uint8_t *dest,
                                   int stride, const TxfmParam *txfm_param) {
  assert(av1_ext_tx_used[txfm_param->tx_set_type][txfm_param->tx_type]);
  int eob = txfm_param->eob;
  int bd = txfm_param->bd;
  int lossless = txfm_param->lossless;
  const int32_t *src = cast_to_int32(input);
  const TX_TYPE tx_type = txfm_param->tx_type;
  if (lossless) {
    assert(tx_type == DCT_DCT);
    av1_highbd_iwht4x4_add(input, dest, stride, eob, bd);
    return;
  }

  av1_inv_txfm2d_add_4x4_c(src, CONVERT_TO_SHORTPTR(dest), stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_4x8_c(const tran_low_t *input, uint8_t *dest,
                                   int stride, const TxfmParam *txfm_param) {
  assert(av1_ext_tx_used[txfm_param->tx_set_type][txfm_param->tx_type]);
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_4x8_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                           txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_8x4_c(const tran_low_t *input, uint8_t *dest,
                                   int stride, const TxfmParam *txfm_param) {
  assert(av1_ext_tx_used[txfm_param->tx_set_type][txfm_param->tx_type]);
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_8x4_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                           txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_16x32_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_16x32_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_32x16_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_32x16_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_16x4_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_16x4_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_4x16_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_4x16_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_32x8_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_32x8_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_8x32_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_8x32_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_32x64_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_32x64_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_64x32_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_64x32_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_16x64_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_16x64_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_64x16_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_64x16_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                             txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_8x8_c(const tran_low_t *input, uint8_t *dest,
                                   int stride, const TxfmParam *txfm_param) {
  int bd = txfm_param->bd;
  const TX_TYPE tx_type = txfm_param->tx_type;
  const int32_t *src = cast_to_int32(input);

  av1_inv_txfm2d_add_8x8_c(src, CONVERT_TO_SHORTPTR(dest), stride, tx_type, bd);
}

void av1_highbd_inv_txfm_add_16x16_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  int bd = txfm_param->bd;
  const TX_TYPE tx_type = txfm_param->tx_type;
  const int32_t *src = cast_to_int32(input);

  av1_inv_txfm2d_add_16x16_c(src, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
}

void av1_highbd_inv_txfm_add_8x16_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_8x16_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_16x8_c(const tran_low_t *input, uint8_t *dest,
                                    int stride, const TxfmParam *txfm_param) {
  const int32_t *src = cast_to_int32(input);
  av1_inv_txfm2d_add_16x8_c(src, CONVERT_TO_SHORTPTR(dest), stride,
                            txfm_param->tx_type, txfm_param->bd);
}

void av1_highbd_inv_txfm_add_32x32_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int bd = txfm_param->bd;
  const TX_TYPE tx_type = txfm_param->tx_type;
  const int32_t *src = cast_to_int32(input);

  av1_inv_txfm2d_add_32x32_c(src, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
}

void av1_highbd_inv_txfm_add_64x64_c(const tran_low_t *input, uint8_t *dest,
                                     int stride, const TxfmParam *txfm_param) {
  const int bd = txfm_param->bd;
  const TX_TYPE tx_type = txfm_param->tx_type;
  const int32_t *src = cast_to_int32(input);
  assert(tx_type == DCT_DCT);
  av1_inv_txfm2d_add_64x64_c(src, CONVERT_TO_SHORTPTR(dest), stride, tx_type,
                             bd);
}

static void init_txfm_param(const MACROBLOCKD *xd, int plane, TX_SIZE tx_size,
                            TX_TYPE tx_type, int eob, int reduced_tx_set,
                            TxfmParam *txfm_param) {
  (void)plane;
  txfm_param->tx_type = tx_type;
  txfm_param->tx_size = tx_size;
  txfm_param->eob = eob;
  txfm_param->lossless = xd->lossless[xd->mi[0]->segment_id];
  txfm_param->bd = xd->bd;
  txfm_param->is_hbd = is_cur_buf_hbd(xd);
  txfm_param->tx_set_type = av1_get_ext_tx_set_type(
      txfm_param->tx_size, is_inter_block(xd->mi[0]), reduced_tx_set);
}

void av1_highbd_inv_txfm_add_c(const tran_low_t *input, uint8_t *dest,
                               int stride, const TxfmParam *txfm_param) {
  assert(av1_ext_tx_used[txfm_param->tx_set_type][txfm_param->tx_type]);
  const TX_SIZE tx_size = txfm_param->tx_size;
  switch (tx_size) {
    case TX_32X32:
      av1_highbd_inv_txfm_add_32x32_c(input, dest, stride, txfm_param);
      break;
    case TX_16X16:
      av1_highbd_inv_txfm_add_16x16_c(input, dest, stride, txfm_param);
      break;
    case TX_8X8:
      av1_highbd_inv_txfm_add_8x8_c(input, dest, stride, txfm_param);
      break;
    case TX_4X8:
      av1_highbd_inv_txfm_add_4x8_c(input, dest, stride, txfm_param);
      break;
    case TX_8X4:
      av1_highbd_inv_txfm_add_8x4_c(input, dest, stride, txfm_param);
      break;
    case TX_8X16:
      av1_highbd_inv_txfm_add_8x16_c(input, dest, stride, txfm_param);
      break;
    case TX_16X8:
      av1_highbd_inv_txfm_add_16x8_c(input, dest, stride, txfm_param);
      break;
    case TX_16X32:
      av1_highbd_inv_txfm_add_16x32_c(input, dest, stride, txfm_param);
      break;
    case TX_32X16:
      av1_highbd_inv_txfm_add_32x16_c(input, dest, stride, txfm_param);
      break;
    case TX_64X64:
      av1_highbd_inv_txfm_add_64x64_c(input, dest, stride, txfm_param);
      break;
    case TX_32X64:
      av1_highbd_inv_txfm_add_32x64_c(input, dest, stride, txfm_param);
      break;
    case TX_64X32:
      av1_highbd_inv_txfm_add_64x32_c(input, dest, stride, txfm_param);
      break;
    case TX_16X64:
      av1_highbd_inv_txfm_add_16x64_c(input, dest, stride, txfm_param);
      break;
    case TX_64X16:
      av1_highbd_inv_txfm_add_64x16_c(input, dest, stride, txfm_param);
      break;
    case TX_4X4:
      // this is like av1_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      av1_highbd_inv_txfm_add_4x4_c(input, dest, stride, txfm_param);
      break;
    case TX_16X4:
      av1_highbd_inv_txfm_add_16x4_c(input, dest, stride, txfm_param);
      break;
    case TX_4X16:
      av1_highbd_inv_txfm_add_4x16_c(input, dest, stride, txfm_param);
      break;
    case TX_8X32:
      av1_highbd_inv_txfm_add_8x32_c(input, dest, stride, txfm_param);
      break;
    case TX_32X8:
      av1_highbd_inv_txfm_add_32x8_c(input, dest, stride, txfm_param);
      break;
    default: assert(0 && "Invalid transform size"); break;
  }
}

void av1_inv_txfm_add_c(const tran_low_t *dqcoeff, uint8_t *dst, int stride,
                        const TxfmParam *txfm_param) {
  const TX_SIZE tx_size = txfm_param->tx_size;
  DECLARE_ALIGNED(32, uint16_t, tmp[MAX_TX_SQUARE]);
  int tmp_stride = MAX_TX_SIZE;
  int w = tx_size_wide[tx_size];
  int h = tx_size_high[tx_size];
  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      tmp[r * tmp_stride + c] = dst[r * stride + c];
    }
  }

  av1_highbd_inv_txfm_add(dqcoeff, CONVERT_TO_BYTEPTR(tmp), tmp_stride,
                          txfm_param);

  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      dst[r * stride + c] = (uint8_t)tmp[r * tmp_stride + c];
    }
  }
}

void av1_inverse_transform_block(const MACROBLOCKD *xd,
                                 const tran_low_t *dqcoeff, int plane,
                                 TX_TYPE tx_type, TX_SIZE tx_size, uint8_t *dst,
                                 int stride, int eob, int reduced_tx_set) {
  if (!eob) return;

  assert(eob <= av1_get_max_eob(tx_size));

  TxfmParam txfm_param;

#if CONFIG_DSPL_RESIDUAL
  DSPL_TYPE dspl_type = xd->mi[0]->dspl_type;
  if (plane == 0 && dspl_type == DSPL_XY) {
    const TX_SIZE new_tx_size = half_tx_size_map[tx_size];
    const uint8_t txw = tx_size_wide[tx_size], txh = tx_size_high[tx_size];
    const uint8_t dspl_txw = tx_size_wide[new_tx_size],
                  dspl_txh = tx_size_high[new_tx_size];

    int fn_stride = MAX_TX_SIZE;
    assert(fn_stride >= dspl_txw);

    // Buffers
    DECLARE_ALIGNED(32, tran_low_t, scan_buf[MAX_TX_SQUARE]);
    DECLARE_ALIGNED(32, tran_low_t, dqcoeff_buf[MAX_TX_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, buf_diff[MAX_TX_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, buf_up_diff[MAX_TX_SQUARE]);

    int dspl_eob = eob;
    init_txfm_param(xd, plane, tx_size, tx_type, dspl_eob, reduced_tx_set,
                    &txfm_param);

    // Unpack coefficients
    const SCAN_ORDER *scan_order = get_scan(tx_size, tx_type),
                     *dspl_scan_order = get_scan(new_tx_size, tx_type);
    int size = av1_get_max_eob(tx_size);
    memset(scan_buf, 0, size * sizeof(tran_low_t));
    scan_array(dqcoeff, scan_buf, eob, scan_order);
    memset(dqcoeff_buf, 0, txw * txh * sizeof(tran_low_t));
    iscan_array(scan_buf, dqcoeff_buf, dspl_eob, dspl_scan_order);

    // Invert and compute difference
    av1_inverse_transform_block_diff(xd, dqcoeff_buf, plane, tx_type,
                                     new_tx_size, buf_diff, fn_stride, dspl_eob,
                                     reduced_tx_set);

    // Upsample
    av1_signed_up2(buf_diff, dspl_txh, dspl_txw, fn_stride, buf_up_diff,
                   fn_stride, 1, 1, txfm_param.bd);

    // Add to output
    for (int r = 0; r < txh; ++r) {
      for (int c = 0; c < txw; ++c) {
        const tran_low_t residue = (tran_low_t)buf_up_diff[r * fn_stride + c];
        const uint16_t out =
            highbd_clip_pixel_add(dst[r * stride + c], residue, txfm_param.bd);
        dst[r * stride + c] = (uint8_t)out;
      }
    }

  } else {
    init_txfm_param(xd, plane, tx_size, tx_type, eob, reduced_tx_set,
                    &txfm_param);
    assert(av1_ext_tx_used[txfm_param.tx_set_type][txfm_param.tx_type]);

    if (txfm_param.is_hbd) {
      av1_highbd_inv_txfm_add(dqcoeff, dst, stride, &txfm_param);
    } else {
      av1_inv_txfm_add(dqcoeff, dst, stride, &txfm_param);
    }
  }
#else
  init_txfm_param(xd, plane, tx_size, tx_type, eob, reduced_tx_set,
                  &txfm_param);
  assert(av1_ext_tx_used[txfm_param.tx_set_type][txfm_param.tx_type]);

  if (txfm_param.is_hbd) {
    av1_highbd_inv_txfm_add(dqcoeff, dst, stride, &txfm_param);
  } else {
    av1_inv_txfm_add(dqcoeff, dst, stride, &txfm_param);
  }
#endif
}

#if CONFIG_DSPL_RESIDUAL
void av1_inverse_transform_block_diff(const MACROBLOCKD *xd,
                                      const tran_low_t *dqcoeff, int plane,
                                      TX_TYPE tx_type, TX_SIZE tx_size,
                                      int16_t *dst, int stride, int eob,
                                      int reduced_tx_set) {
  if (!eob) return;

  // TODO(singhprakhar): extend this method to higher bitdepths
  assert(xd->bd == 8);

  const uint8_t txw = tx_size_wide[tx_size], txh = tx_size_high[tx_size];

  DECLARE_ALIGNED(32, uint16_t, recon_diffu[MAX_SB_SQUARE]);
  uint8_t *recon = CONVERT_TO_BYTEPTR(recon_diffu);

  for (int r = 0; r < txh; ++r)
    for (int c = 0; c < txw; ++c) recon_diffu[r * MAX_TX_SIZE + c] = 1023;

  TxfmParam txfm_param;
  init_txfm_param(xd, plane, tx_size, tx_type, eob, reduced_tx_set,
                  &txfm_param);
  txfm_param.bd = 12;
  txfm_param.is_hbd = 1;

  av1_highbd_inv_txfm_add(dqcoeff, recon, MAX_TX_SIZE, &txfm_param);

  for (int r = 0; r < txh; ++r)
    for (int c = 0; c < txw; ++c)
      dst[r * stride + c] = (int16_t)recon_diffu[r * MAX_TX_SIZE + c] - 1023;
}
#endif
