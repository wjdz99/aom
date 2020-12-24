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

#include "av1/common/cfl.h"
#include "av1/common/common_data.h"
#include "av1/common/onyxc_int.h"

#include "config/av1_rtcd.h"

void cfl_init(CFL_CTX *cfl, const SequenceHeader *seq_params) {
  assert(block_size_wide[CFL_MAX_BLOCK_SIZE] == CFL_BUF_LINE);
  assert(block_size_high[CFL_MAX_BLOCK_SIZE] == CFL_BUF_LINE);

  memset(&cfl->recon_buf_q3, 0, sizeof(cfl->recon_buf_q3));
  memset(&cfl->ac_buf_q3, 0, sizeof(cfl->ac_buf_q3));
#if CONFIG_CFL_SEARCH_VERSION_1
  memset(&cfl->recon_Neighbor, 0, sizeof(NeiborPix) * (CFL_BUF_LINE << 1));
#endif
  cfl->subsampling_x = seq_params->subsampling_x;
  cfl->subsampling_y = seq_params->subsampling_y;
  cfl->are_parameters_computed = 0;
  cfl->store_y = 0;
  // The DC_PRED cache is disabled by default and is only enabled in
  // cfl_rd_pick_alpha
#if CONFIG_CFL_SEARCH_VERSION_1
  cfl->use_search_res_cache = 0;
  cfl->search_res_is_cached[CFL_PRED_U] = 0;
  cfl->search_res_is_cached[CFL_PRED_V] = 0;
  cfl->use_dc = 0;
  cfl->use_left = 0;
  cfl->use_up = 0;
#else
  cfl->use_dc_pred_cache = 0;
  cfl->dc_pred_is_cached[CFL_PRED_U] = 0;
  cfl->dc_pred_is_cached[CFL_PRED_V] = 0;
#endif
}
#if CONFIG_CFL_SEARCH_VERSION_1
void cfl_get_dc(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                CFL_PRED_TYPE pred_plane, TX_SIZE txsize) {
  int sum = 0;
  uint8_t *pdst = dst;
  const int width = tx_size_wide[txsize];
  const int height = tx_size_high[txsize];

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      sum += pdst[j];
    }
    pdst += dst_stride;
  }

  uint8_t shift = tx_size_wide_log2[txsize] + tx_size_high_log2[txsize];
  uint8_t offset = 1 << (shift - 1);
  uint8_t dc = (sum + offset) >> shift;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      dst[i] = dc;
    }
    dst += dst_stride;
  }
}

#endif
#if CONFIG_CFL_SEARCH_VERSION_1
void cfl_search(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                CFL_PRED_TYPE pred_plane,
                TX_SIZE txsize) {  // plane为0是U分量，为1是V分量

  int width = tx_size_wide[txsize];
  int height = tx_size_high[txsize];
  assert(pred_plane < CFL_PRED_PLANES);
  assert(width <= CFL_BUF_LINE);
  assert(height <= CFL_BUF_LINE);

  CFL_CTX *const cfl = &xd->cfl;

  uint16_t *recon = cfl->recon_buf_q3;  //亮度重建像素
  NeiborPix *rec_nb = cfl->recon_Neighbor;  //亮度重建上一行和左一列像素

  int flag = 0;
  int DValue;
  // 搜索最佳值
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      DValue = 1024;
      if (cfl->use_up) {
        for (int k = 0; k < width; k++)  //遍历上一行像素
        {
          if (abs(recon[i] - rec_nb[k].LumaPixel) < DValue) {
            DValue = abs(recon[i] - rec_nb[k].LumaPixel);
            flag = k;
          } else if (abs(recon[i] - rec_nb[k].LumaPixel) == DValue) {
            if ((((i + 1) - rec_nb[k].positionx) *
                     ((i + 1) - rec_nb[k].positionx) +
                 ((j + 1) - rec_nb[k].positiony) *
                     ((j + 1) - rec_nb[k].positiony)) <
                (((i + 1) - rec_nb[flag].positionx) *
                     ((i + 1) - rec_nb[flag].positionx) +
                 ((j + 1) - rec_nb[flag].positiony) *
                     ((j + 1) - rec_nb[flag].positiony))) {
              flag = k;
            }
          }
        }
      }

      if (cfl->use_left) {
        for (int k = 0; k < height; k++)  //遍历左一列像素
        {
          if (abs(recon[i] - rec_nb[CFL_BUF_LINE + k].LumaPixel) < DValue) {
            DValue = abs(recon[i] - rec_nb[CFL_BUF_LINE + k].LumaPixel);
            flag = k + CFL_BUF_LINE;
          } else if (abs(recon[i] - rec_nb[CFL_BUF_LINE + k].LumaPixel) ==
                     DValue) {
            if ((((i + 1) - rec_nb[CFL_BUF_LINE + k].positionx) *
                     ((i + 1) - rec_nb[CFL_BUF_LINE + k].positionx) +
                 ((j + 1) - rec_nb[CFL_BUF_LINE + k].positiony) *
                     ((j + 1) - rec_nb[CFL_BUF_LINE + k].positiony)) <
                (((i + 1) - rec_nb[flag].positionx) *
                     ((i + 1) - rec_nb[flag].positionx) +
                 ((j + 1) - rec_nb[flag].positiony) *
                     ((j + 1) - rec_nb[flag].positiony))) {
              flag = k + CFL_BUF_LINE;
            }
          }
        }
      }

      if (pred_plane == CFL_PRED_U) {
        dst[i] = rec_nb[flag].ChromCbPixel;
      } else  // 0是U分量
      {
        dst[i] = rec_nb[flag].ChromCrPixel;
      }
    }
    recon += CFL_BUF_LINE;
    dst += dst_stride;
  }
}

void cfl_store_search_res(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                          CFL_PRED_TYPE pred_plane, TX_SIZE txsize) {
  int width = tx_size_wide[txsize];
  int height = tx_size_high[txsize];
  CFL_CTX *const cfl = &xd->cfl;

  int16_t *res_cache = cfl->search_res_cache[pred_plane];
  for (int i = 0; i < height; i++) {
    memcpy(res_cache, dst, width);
    res_cache += CFL_BUF_LINE;
    dst += dst_stride;
  }
}

static void cfl_load_search_res_lbd(const int16_t *search_res_cache,
                                    uint8_t *dst, int dst_stride, int width,
                                    int height) {
  for (int j = 0; j < height; j++) {
    memcpy(dst, search_res_cache, width);
    dst += dst_stride;
    search_res_cache += CFL_BUF_LINE;
  }
}

static void cfl_load_search_res_hbd(const int16_t *search_res_cache,
                                    uint16_t *dst, int dst_stride, int width,
                                    int height) {
  const size_t num_bytes = width << 1;
  for (int j = 0; j < height; j++) {
    memcpy(dst, search_res_cache, num_bytes);
    dst += dst_stride;
    search_res_cache += CFL_BUF_LINE;
  }
}
void cfl_load_search_res(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                         TX_SIZE tx_size, CFL_PRED_TYPE pred_plane) {
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  assert(pred_plane < CFL_PRED_PLANES);
  assert(width <= CFL_BUF_LINE);
  assert(height <= CFL_BUF_LINE);

  if (is_cur_buf_hbd(xd)) {
    uint16_t *dst_16 = CONVERT_TO_SHORTPTR(dst);
    cfl_load_search_res_hbd(xd->cfl.search_res_cache[pred_plane], dst_16,
                            dst_stride, width, height);
    return;
  }
  cfl_load_search_res_lbd(xd->cfl.search_res_cache[pred_plane], dst, dst_stride,
                          width, height);
}
#else
void cfl_store_dc_pred(MACROBLOCKD *const xd, const uint8_t *input,
                       CFL_PRED_TYPE pred_plane, int width) {
  assert(pred_plane < CFL_PRED_PLANES);
  assert(width <= CFL_BUF_LINE);

  if (is_cur_buf_hbd(xd)) {
    uint16_t *const input_16 = CONVERT_TO_SHORTPTR(input);
    memcpy(xd->cfl.dc_pred_cache[pred_plane], input_16, width << 1);
    return;
  }

  memcpy(xd->cfl.dc_pred_cache[pred_plane], input, width);
}

static void cfl_load_dc_pred_lbd(const int16_t *dc_pred_cache, uint8_t *dst,
                                 int dst_stride, int width, int height) {
  for (int j = 0; j < height; j++) {
    memcpy(dst, dc_pred_cache, width);
    dst += dst_stride;
  }
}

static void cfl_load_dc_pred_hbd(const int16_t *dc_pred_cache, uint16_t *dst,
                                 int dst_stride, int width, int height) {
  const size_t num_bytes = width << 1;
  for (int j = 0; j < height; j++) {
    memcpy(dst, dc_pred_cache, num_bytes);
    dst += dst_stride;
  }
}
void cfl_load_dc_pred(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                      TX_SIZE tx_size, CFL_PRED_TYPE pred_plane) {
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  assert(pred_plane < CFL_PRED_PLANES);
  assert(width <= CFL_BUF_LINE);
  assert(height <= CFL_BUF_LINE);
  if (is_cur_buf_hbd(xd)) {
    uint16_t *dst_16 = CONVERT_TO_SHORTPTR(dst);
    cfl_load_dc_pred_hbd(xd->cfl.dc_pred_cache[pred_plane], dst_16, dst_stride,
                         width, height);
    return;
  }
  cfl_load_dc_pred_lbd(xd->cfl.dc_pred_cache[pred_plane], dst, dst_stride,
                       width, height);
}
#endif
// Due to frame boundary issues, it is possible that the total area covered by
// chroma exceeds that of luma. When this happens, we fill the missing pixels by
// repeating the last columns and/or rows.
static INLINE void cfl_pad(CFL_CTX *cfl, int width, int height) {
  const int diff_width = width - cfl->buf_width;
  const int diff_height = height - cfl->buf_height;

  if (diff_width > 0) {
    const int min_height = height - diff_height;
    uint16_t *recon_buf_q3 = cfl->recon_buf_q3 + (width - diff_width);
    for (int j = 0; j < min_height; j++) {
      const uint16_t last_pixel = recon_buf_q3[-1];
      assert(recon_buf_q3 + diff_width <= cfl->recon_buf_q3 + CFL_BUF_SQUARE);
      for (int i = 0; i < diff_width; i++) {
        recon_buf_q3[i] = last_pixel;
      }
      recon_buf_q3 += CFL_BUF_LINE;
    }
    cfl->buf_width = width;
  }
  if (diff_height > 0) {
    uint16_t *recon_buf_q3 =
        cfl->recon_buf_q3 + ((height - diff_height) * CFL_BUF_LINE);
    for (int j = 0; j < diff_height; j++) {
      const uint16_t *last_row_q3 = recon_buf_q3 - CFL_BUF_LINE;
      assert(recon_buf_q3 + width <= cfl->recon_buf_q3 + CFL_BUF_SQUARE);
      for (int i = 0; i < width; i++) {
        recon_buf_q3[i] = last_row_q3[i];
      }
      recon_buf_q3 += CFL_BUF_LINE;
    }
    cfl->buf_height = height;
  }
}

static void subtract_average_c(const uint16_t *src, int16_t *dst, int width,
                               int height, int round_offset, int num_pel_log2) {
  int sum = round_offset;
  const uint16_t *recon = src;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      sum += recon[i];
    }
    recon += CFL_BUF_LINE;
  }
  const int avg = sum >> num_pel_log2;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      dst[i] = src[i] - avg;
    }
    src += CFL_BUF_LINE;
    dst += CFL_BUF_LINE;
  }
}

CFL_SUB_AVG_FN(c)

static INLINE int cfl_idx_to_alpha(uint8_t alpha_idx, int8_t joint_sign,
                                   CFL_PRED_TYPE pred_type) {
  const int alpha_sign = (pred_type == CFL_PRED_U) ? CFL_SIGN_U(joint_sign)
                                                   : CFL_SIGN_V(joint_sign);
  if (alpha_sign == CFL_SIGN_ZERO) return 0;
  const int abs_alpha_q3 =
      (pred_type == CFL_PRED_U) ? CFL_IDX_U(alpha_idx) : CFL_IDX_V(alpha_idx);
  return (alpha_sign == CFL_SIGN_POS) ? abs_alpha_q3 + 1 : -abs_alpha_q3 - 1;
}

static INLINE void cfl_predict_lbd_c(const int16_t *ac_buf_q3, uint8_t *dst,
                                     int dst_stride, int alpha_q3, int width,
                                     int height) {
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      dst[i] = clip_pixel(get_scaled_luma_q0(alpha_q3, ac_buf_q3[i]) + dst[i]);
    }
    dst += dst_stride;
    ac_buf_q3 += CFL_BUF_LINE;
  }
}

CFL_PREDICT_FN(c, lbd)

void cfl_predict_hbd_c(const int16_t *ac_buf_q3, uint16_t *dst, int dst_stride,
                       int alpha_q3, int bit_depth, int width, int height) {
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      dst[i] = clip_pixel_highbd(
          get_scaled_luma_q0(alpha_q3, ac_buf_q3[i]) + dst[i], bit_depth);
    }
    dst += dst_stride;
    ac_buf_q3 += CFL_BUF_LINE;
  }
}

CFL_PREDICT_FN(c, hbd)

static void cfl_compute_parameters(MACROBLOCKD *const xd, TX_SIZE tx_size) {
  CFL_CTX *const cfl = &xd->cfl;
  // Do not call cfl_compute_parameters multiple time on the same values.
  assert(cfl->are_parameters_computed == 0);

  cfl_pad(cfl, tx_size_wide[tx_size], tx_size_high[tx_size]);
  cfl_get_subtract_average_fn(tx_size)(cfl->recon_buf_q3, cfl->ac_buf_q3);
  cfl->are_parameters_computed = 1;
}

void cfl_predict_block(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                       TX_SIZE tx_size, int plane) {
  CFL_CTX *const cfl = &xd->cfl;
  MB_MODE_INFO *mbmi = xd->mi[0];
  assert(is_cfl_allowed(xd));

  if (!cfl->are_parameters_computed) cfl_compute_parameters(xd, tx_size);

  const int alpha_q3 =
      cfl_idx_to_alpha(mbmi->cfl_alpha_idx, mbmi->cfl_alpha_signs, plane - 1);
  assert((tx_size_high[tx_size] - 1) * CFL_BUF_LINE + tx_size_wide[tx_size] <=
         CFL_BUF_SQUARE);
  if (is_cur_buf_hbd(xd)) {
    uint16_t *dst_16 = CONVERT_TO_SHORTPTR(dst);
    cfl_get_predict_hbd_fn(tx_size)(cfl->ac_buf_q3, dst_16, dst_stride,
                                    alpha_q3, xd->bd);
    return;
  }
  cfl_get_predict_lbd_fn(tx_size)(cfl->ac_buf_q3, dst, dst_stride, alpha_q3);
}

static void cfl_luma_subsampling_420_lbd_c(const uint8_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  for (int j = 0; j < height; j += 2) {
    for (int i = 0; i < width; i += 2) {
      const int bot = i + input_stride;
      output_q3[i >> 1] =
          (input[i] + input[i + 1] + input[bot] + input[bot + 1]) << 1;
    }
    input += input_stride << 1;
    output_q3 += CFL_BUF_LINE;
  }
}

static void cfl_luma_subsampling_422_lbd_c(const uint8_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i += 2) {
      output_q3[i >> 1] = (input[i] + input[i + 1]) << 2;
    }
    input += input_stride;
    output_q3 += CFL_BUF_LINE;
  }
}

static void cfl_luma_subsampling_444_lbd_c(const uint8_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      output_q3[i] = input[i] << 3;
    }
    input += input_stride;
    output_q3 += CFL_BUF_LINE;
  }
}

static void cfl_luma_subsampling_420_hbd_c(const uint16_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  for (int j = 0; j < height; j += 2) {
    for (int i = 0; i < width; i += 2) {
      const int bot = i + input_stride;
      output_q3[i >> 1] =
          (input[i] + input[i + 1] + input[bot] + input[bot + 1]) << 1;
    }
    input += input_stride << 1;
    output_q3 += CFL_BUF_LINE;
  }
}

static void cfl_luma_subsampling_422_hbd_c(const uint16_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i += 2) {
      output_q3[i >> 1] = (input[i] + input[i + 1]) << 2;
    }
    input += input_stride;
    output_q3 += CFL_BUF_LINE;
  }
}

static void cfl_luma_subsampling_444_hbd_c(const uint16_t *input,
                                           int input_stride,
                                           uint16_t *output_q3, int width,
                                           int height) {
  assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      output_q3[i] = input[i] << 3;
    }
    input += input_stride;
    output_q3 += CFL_BUF_LINE;
  }
}

CFL_GET_SUBSAMPLE_FUNCTION(c)

static INLINE cfl_subsample_hbd_fn cfl_subsampling_hbd(TX_SIZE tx_size,
                                                       int sub_x, int sub_y) {
  if (sub_x == 1) {
    if (sub_y == 1) {
      return cfl_get_luma_subsampling_420_hbd(tx_size);
    }
    return cfl_get_luma_subsampling_422_hbd(tx_size);
  }
  return cfl_get_luma_subsampling_444_hbd(tx_size);
}

static INLINE cfl_subsample_lbd_fn cfl_subsampling_lbd(TX_SIZE tx_size,
                                                       int sub_x, int sub_y) {
  if (sub_x == 1) {
    if (sub_y == 1) {
      return cfl_get_luma_subsampling_420_lbd(tx_size);
    }
    return cfl_get_luma_subsampling_422_lbd(tx_size);
  }
  return cfl_get_luma_subsampling_444_lbd(tx_size);
}

static void cfl_store(CFL_CTX *cfl, const uint8_t *input, int input_stride,
                      int row, int col, TX_SIZE tx_size, int use_hbd) {
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const int tx_off_log2 = tx_size_wide_log2[0];
  const int sub_x = cfl->subsampling_x;
  const int sub_y = cfl->subsampling_y;
  const int store_row = row << (tx_off_log2 - sub_y);
  const int store_col = col << (tx_off_log2 - sub_x);
  const int store_height = height >> sub_y;
  const int store_width = width >> sub_x;

  // Invalidate current parameters
  cfl->are_parameters_computed = 0;

  // Store the surface of the pixel buffer that was written to, this way we
  // can manage chroma overrun (e.g. when the chroma surfaces goes beyond the
  // frame boundary)
  if (col == 0 && row == 0) {
    cfl->buf_width = store_width;
    cfl->buf_height = store_height;
  } else {
    cfl->buf_width = OD_MAXI(store_col + store_width, cfl->buf_width);
    cfl->buf_height = OD_MAXI(store_row + store_height, cfl->buf_height);
  }

  // Check that we will remain inside the pixel buffer.
  assert(store_row + store_height <= CFL_BUF_LINE);
  assert(store_col + store_width <= CFL_BUF_LINE);

  // Store the input into the CfL pixel buffer
  uint16_t *recon_buf_q3 =
      cfl->recon_buf_q3 + (store_row * CFL_BUF_LINE + store_col);

  if (use_hbd) {
    cfl_subsampling_hbd(tx_size, sub_x, sub_y)(CONVERT_TO_SHORTPTR(input),
                                               input_stride, recon_buf_q3);
  } else {
    cfl_subsampling_lbd(tx_size, sub_x, sub_y)(input, input_stride,
                                               recon_buf_q3);
  }
}
#if CONFIG_CFL_SEARCH_VERSION_1
void cfl_store_nb(MACROBLOCKD *const xd, BLOCK_SIZE bsize) {
  CFL_CTX *const cfl = &xd->cfl;
  cfl->use_dc = 0;
  cfl->use_left = 0;
  cfl->use_up = 0;

  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  struct macroblockd_plane *const pdU = &xd->plane[AOM_PLANE_U];
  struct macroblockd_plane *const pdV = &xd->plane[AOM_PLANE_V];

  const int subWidth = pdU->width;
  const int subHeight = pdU->height;

  NeiborPix *recNb = cfl->recon_Neighbor;

  uint8_t *dst = &pd->dst.buf[0];
  const int dst_stride = pd->dst.stride;
  uint8_t *above_ref = dst - dst_stride;
  uint8_t *left_ref = dst - 1;
  const int ref_stride = dst_stride;

  uint8_t *pChromaCbSrc = &pdU->dst.buf[0];
  const int pChromaCbStride = pdU->dst.stride;

  uint8_t *pChromaCrSrc = &pdV->dst.buf[0];
  const int pChromaCrStride = pdV->dst.stride;
  uint8_t *pCurrCb = NULL;
  uint8_t *pCurrCr = NULL;

  if (xd->chroma_up_available) {
    cfl->use_up = 1;
    pCurrCb = pChromaCbSrc - pChromaCbStride;
    pCurrCr = pChromaCrSrc - pChromaCrStride;
    if (block_size_wide[bsize] == 4) {
      for (int i = 0; i < subWidth; i++) {
        recNb[i].LumaPixel = above_ref[i] << 3;
        recNb[i].positionx = i;
        recNb[i].positiony = -1;
      }
    } else {
      for (int i = 0; i < subWidth; i++) {
        recNb[i].LumaPixel = (above_ref[2 * i] + above_ref[2 * i + 1] +
                              (above_ref - ref_stride)[2 * i] +
                              (above_ref - ref_stride)[2 * i + 1])
                             << 1;
        recNb[i].positionx = i;
        recNb[i].positiony = -1;
      }
    }

    for (int i = 0; i < subWidth; i++) {
      recNb[i].ChromCbPixel = pCurrCb[i];
      recNb[i].ChromCrPixel = pCurrCr[i];
    }
  }

  if (xd->chroma_left_available) {
    cfl->use_left = 1;
    pCurrCb = pChromaCbSrc - 1;
    pCurrCr = pChromaCrSrc - 1;
    if (block_size_high[bsize] == 4) {
      for (int i = 0; i < subHeight; i++) {
        recNb[CFL_BUF_LINE + i].LumaPixel = left_ref[0] << 3;
        recNb[CFL_BUF_LINE + i].positionx = -1;
        recNb[CFL_BUF_LINE + i].positiony = i;
        left_ref += ref_stride;
      }
    } else {
      for (int i = 0; i < subHeight; i++) {
        recNb[CFL_BUF_LINE + i].LumaPixel =
            (left_ref[0] + (left_ref + ref_stride)[0] + (left_ref - 1)[0] +
             (left_ref + ref_stride - 1)[0])
            << 1;
        recNb[CFL_BUF_LINE + i].positionx = -1;
        recNb[CFL_BUF_LINE + i].positiony = i;
        left_ref += ref_stride << 1;
      }
    }

    for (int i = 0; i < subHeight; i++) {
      recNb[CFL_BUF_LINE + i].ChromCbPixel = pCurrCb[0];
      recNb[CFL_BUF_LINE + i].ChromCrPixel = pCurrCr[0];
      pCurrCb += pChromaCbStride;
      pCurrCr += pChromaCrStride;
    }
  }

  if (!cfl->use_left && !cfl->use_up) {
    cfl->use_dc = 1;
  }
}
#endif
void cfl_store_tx(MACROBLOCKD *const xd, int row, int col, TX_SIZE tx_size) {
  CFL_CTX *const cfl = &xd->cfl;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  uint8_t *dst =
      &pd->dst.buf[(row * pd->dst.stride + col) << tx_size_wide_log2[0]];
  const int mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
  const int mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
  const int row_offset = mi_row - xd->mi[0]->chroma_ref_info.mi_row_chroma_base;
  const int col_offset = mi_col - xd->mi[0]->chroma_ref_info.mi_col_chroma_base;

  cfl_store(cfl, dst, pd->dst.stride, row + row_offset, col + col_offset,
            tx_size, is_cur_buf_hbd(xd));
}

void cfl_store_block(MACROBLOCKD *const xd, BLOCK_SIZE bsize, TX_SIZE tx_size) {
  CFL_CTX *const cfl = &xd->cfl;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  const int width = max_intra_block_width(xd, bsize, AOM_PLANE_Y, tx_size);
  const int height = max_intra_block_height(xd, bsize, AOM_PLANE_Y, tx_size);
  const int mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
  const int mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
  const int row_offset = mi_row - xd->mi[0]->chroma_ref_info.mi_row_chroma_base;
  const int col_offset = mi_col - xd->mi[0]->chroma_ref_info.mi_col_chroma_base;

  tx_size = get_tx_size(width, height);
  cfl_store(cfl, pd->dst.buf, pd->dst.stride, row_offset, col_offset, tx_size,
            is_cur_buf_hbd(xd));
}
