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

#include "aom/internal/aom_codec_internal.h"

void cfl_init(CFL_CTX *cfl, AV1_COMMON *cm) {
  if (!((cm->subsampling_x == 0 && cm->subsampling_y == 0) ||
        (cm->subsampling_x == 1 && cm->subsampling_y == 1))) {
    aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
                       "Only 4:4:4 and 4:2:0 are currently supported by CfL");
  }
  memset(&cfl->ac_con_q3, 0, sizeof(cfl->ac_con_q3));
  cfl->subsampling_x = cm->subsampling_x;
  cfl->subsampling_y = cm->subsampling_y;
  cfl->are_parameters_computed = 0;
  cfl->store_y = 0;
#if CONFIG_CHROMA_SUB8X8 && CONFIG_DEBUG
  cfl_clear_sub8x8_val(cfl);
#endif  // CONFIG_CHROMA_SUB8X8 && CONFIG_DEBUG
}

// TODO(ltrudeau) Add high bit depth version of this function
static int cfl_luma_subsampling_420(const uint8_t *input, int16_t *output,
                                    int input_stride, int width, int height) {
  int sum = 0;
  for (int j = 0; j < height; j += 2) {
    int c = 0;
    for (int i = 0; i < width; i += 2) {
      // TODO(ltrudeau) Use sum instead of average when subsampling
      int sub = ((input[i] + input[i + 1] + input[i + input_stride] +
                  input[i + input_stride + 1] + 2) >>
                 2)
                << 3;

      // Worst case, subsampling sum fits in 13bits (LBD)
      assert(sub < (1 << 13));

      sum += sub;
      output[c++] = sub;
    }
    input += input_stride << 1;
    output += MAX_SB_SIZE;
  }
  // Worst case, subsampling sum fits in 23bits (LBD)
  assert(sum < (1 << 23));
  return sum;
}

// Due to frame boundary issues, it is possible that the total area of
// covered by chroma exceeds that of luma. When this happens, we fill the
// missing by repeating the last columns and/or rows.
static void cfl_pad_ac_con(CFL_CTX *cfl) {
  const int width = cfl->uv_width;
  const int height = cfl->uv_height;

  assert(cfl->ac_width > 0);
  assert(cfl->ac_height > 0);

  const int diff_width = width - cfl->ac_width;
  const int diff_height = height - cfl->ac_height;

  if (diff_width > 0) {
    int16_t *output = cfl->ac_con_q3 + (width - diff_width);

    for (int j = 0; j < height; j++) {
      int last_pixel = output[-1];
      for (int i = 0; i < diff_width; i++) {
        output[i] = last_pixel;
      }
      output += MAX_SB_SIZE;
    }
  }

  if (diff_height > 0) {
    int16_t *output = cfl->ac_con_q3 + (height - diff_height) * MAX_SB_SIZE;

    for (int j = 0; j < diff_height; j++) {
      int16_t *last_row = output - MAX_SB_SIZE;
      for (int i = 0; i < width; i++) {
        output[i] = last_row[i];
      }
      output += MAX_SB_SIZE;
    }
  }
  cfl->ac_width = width;
  cfl->ac_height = height;
}

// CfL computes its own block-level DC_PRED. This is required to compute both
// alpha_cb and alpha_cr before the prediction are computed.
static void cfl_dc_pred(MACROBLOCKD *xd, BLOCK_SIZE plane_bsize) {
  const struct macroblockd_plane *const pd_u = &xd->plane[AOM_PLANE_U];
  const struct macroblockd_plane *const pd_v = &xd->plane[AOM_PLANE_V];

  const uint8_t *const dst_u = pd_u->dst.buf;
  const uint8_t *const dst_v = pd_v->dst.buf;

  const int dst_u_stride = pd_u->dst.stride;
  const int dst_v_stride = pd_v->dst.stride;

  CFL_CTX *const cfl = xd->cfl;

  // Compute DC_PRED until block boundary. We can't assume the neighbor will use
  // the same transform size.
  const int width = max_block_wide(xd, plane_bsize, AOM_PLANE_U)
                    << tx_size_wide_log2[0];
  const int height = max_block_high(xd, plane_bsize, AOM_PLANE_U)
                     << tx_size_high_log2[0];
  // Number of pixel on the top and left borders.
  const int num_pel = width + height;

  int sum_u = 0;
  int sum_v = 0;

// Match behavior of build_intra_predictors (reconintra.c) at superblock
// boundaries:
//
// 127 127 127 .. 127 127 127 127 127 127
// 129  A   B  ..  Y   Z
// 129  C   D  ..  W   X
// 129  E   F  ..  U   V
// 129  G   H  ..  S   T   T   T   T   T
// ..

#if CONFIG_CHROMA_SUB8X8
  if (xd->chroma_up_available && xd->mb_to_right_edge >= 0) {
#else
  if (xd->up_available && xd->mb_to_right_edge >= 0) {
#endif
    // TODO(ltrudeau) replace this with DC_PRED assembly
    for (int i = 0; i < width; i++) {
      sum_u += dst_u[-dst_u_stride + i];
      sum_v += dst_v[-dst_v_stride + i];
    }
  } else {
    sum_u = width * 127;
    sum_v = width * 127;
  }

#if CONFIG_CHROMA_SUB8X8
  if (xd->chroma_left_available && xd->mb_to_bottom_edge >= 0) {
#else
  if (xd->left_available && xd->mb_to_bottom_edge >= 0) {
#endif
    for (int i = 0; i < height; i++) {
      sum_u += dst_u[i * dst_u_stride - 1];
      sum_v += dst_v[i * dst_v_stride - 1];
    }
  } else {
    sum_u += height * 129;
    sum_v += height * 129;
  }

  // TODO(ltrudeau) Because of max_block_wide and max_block_high, num_pel will
  // not be a power of two. So these divisions will have to use a lookup table.
  cfl->dc_pred[CFL_PRED_U] = (sum_u + (num_pel >> 1)) / num_pel;
  cfl->dc_pred[CFL_PRED_V] = (sum_v + (num_pel >> 1)) / num_pel;
}

static INLINE int cfl_idx_to_alpha(int alpha_idx, int joint_sign,
                                   CFL_PRED_TYPE pred_type) {
  const int alpha_sign = (pred_type == CFL_PRED_U) ? CFL_SIGN_U(joint_sign)
                                                   : CFL_SIGN_V(joint_sign);
  if (alpha_sign == CFL_SIGN_ZERO) return 0;
  const int abs_alpha_q3 =
      (pred_type == CFL_PRED_U) ? CFL_IDX_U(alpha_idx) : CFL_IDX_V(alpha_idx);
  return (alpha_sign == CFL_SIGN_POS) ? abs_alpha_q3 + 1 : -abs_alpha_q3 - 1;
}

// Predict the current transform block using CfL.
void cfl_predict_block(MACROBLOCKD *const xd, uint8_t *dst, int dst_stride,
                       int row, int col, TX_SIZE tx_size, int plane) {
  CFL_CTX *const cfl = xd->cfl;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  col = col << tx_size_wide_log2[0];
  row = row << tx_size_high_log2[0];

  // CfL parameters must be computed before prediction can be done.
  assert(cfl->are_parameters_computed == 1);

  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];
  const int16_t *ac_con_q3 = cfl->ac_con_q3 + (row * MAX_SB_SIZE + col);

  const int dc_pred = cfl->dc_pred[plane - 1];
  const int alpha_q3 =
      cfl_idx_to_alpha(mbmi->cfl_alpha_idx, mbmi->cfl_alpha_signs, plane - 1);

  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      // TODO(ltrudeau) add support for HBD.
      int scaled_luma_q6 = alpha_q3 * ac_con_q3[i];
      int cfl_pred = ROUND_POWER_OF_TWO_SIGNED(scaled_luma_q6, 6) + dc_pred;
      dst[i] = clip_pixel(cfl_pred);
    }
    dst += dst_stride;
    ac_con_q3 += MAX_SB_SIZE;
  }
}

static INLINE void cfl_store(CFL_CTX *cfl, const uint8_t *input,
                             int input_stride, int row, int col, int width,
                             int height) {
  const int sub_x = cfl->subsampling_x;
  const int sub_y = cfl->subsampling_y;
  const int store_width = width >> sub_x;
  const int store_height = height >> sub_y;
  const int store_col = col << (tx_size_wide_log2[0] - sub_x);
  const int store_row = row << (tx_size_high_log2[0] - sub_y);

  // Check that we will remain inside the pixel buffer.
  assert(store_row + store_height <= (MAX_SB_SIZE >> sub_y));
  assert(store_col + store_width <= (MAX_SB_SIZE >> sub_x));

#if CONFIG_DEBUG
  // Invalidate current parameters
  cfl->are_parameters_computed = 0;
#endif

  // Track the width and height of the AC contributions currently in the buffer.
  // This way, we can manage chroma overrun when chroma surface exceeds the luma
  // surface. (which can occur at frame boundaries)
  if (col == 0 && row == 0) {
    cfl->ac_width = store_width;
    cfl->ac_height = store_height;
  } else {
    cfl->ac_width = OD_MAXI(store_col + store_width, cfl->ac_width);
    cfl->ac_height = OD_MAXI(store_row + store_height, cfl->ac_height);
  }

  int16_t *ac_con_q3 = cfl->ac_con_q3 + (store_row * MAX_SB_SIZE + store_col);

  int sum = 0;
  if (sub_y == 0 && sub_x == 0) {  // In 4:4:4, pixels match 1 to 1
    for (int j = 0; j < store_height; j++) {
      for (int i = 0; i < store_width; i++) {
        ac_con_q3[i] = input[i] << 3;
        sum += input[i];
      }
      ac_con_q3 += MAX_SB_SIZE;
      input += input_stride;
    }
  } else if (sub_y == 1 && sub_x == 1) {
    sum =
        cfl_luma_subsampling_420(input, ac_con_q3, input_stride, width, height);
  } else {
    assert(0);
  }

  // TODO(ltrudeau) replace division with shifts
  int num_pel = store_height * store_width;
  int avg_q3 = (sum + (num_pel >> 1)) / num_pel;
  for (int j = 0; j < store_height; j++) {
    for (int i = 0; i < store_width; i++) {
      ac_con_q3[i] = ac_con_q3[i] - avg_q3;
    }
    ac_con_q3 += MAX_SB_SIZE;
  }
}

#if CONFIG_CHROMA_SUB8X8
// Adjust the row and column of blocks smaller than 8X8, as chroma-referenced
// and non-chroma-referenced blocks are stored together.
static INLINE void sub8x8_adjust_offset(const CFL_CTX *cfl, BLOCK_SIZE bsize,
                                        int *row_out, int *col_out) {
  assert(bsize < BLOCK_8X8);
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];

  // The following code is adapted from the is_chroma_reference() function.
  if ((cfl->mi_row &
       0x01)        // Increment the row index for odd indexed 4X4 blocks
      && (bh == 4)  // But not for 4X8 blocks
      && cfl->subsampling_y) {  // And only when chroma is subsampled
    assert(*row_out == 0);
    (*row_out)++;
  }

  if ((cfl->mi_col &
       0x01)        // Increment the col index for odd indexed 4X4 blocks
      && (bw == 4)  // But not for 8X4 blocks
      && cfl->subsampling_x) {  // And only when chroma is subsampled
    assert(*col_out == 0);
    (*col_out)++;
  }
}
#if CONFIG_DEBUG
static INLINE void sub8x8_set_val(CFL_CTX *cfl, int row, int col, int val_high,
                                  int val_wide) {
  for (int val_r = 0; val_r < val_high; val_r++) {
    assert(row + val_r < 2);
    int row_off = (row + val_r) * 2;
    for (int val_c = 0; val_c < val_wide; val_c++) {
      assert(col + val_c < 2);
      cfl->sub8x8_val[row_off + col + val_c] = 1;
    }
  }
}
#endif  // CONFIG_DEBUG
#endif  // CONFIG_CHROMA_SUB8X8

void cfl_store_tx(MACROBLOCKD *const xd, int row, int col, TX_SIZE tx_size,
                  BLOCK_SIZE bsize) {
  CFL_CTX *const cfl = xd->cfl;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  uint8_t *dst =
      &pd->dst.buf[(row * pd->dst.stride + col) << tx_size_wide_log2[0]];
  (void)bsize;
#if CONFIG_CHROMA_SUB8X8
  if (bsize < BLOCK_8X8) {
    sub8x8_adjust_offset(cfl, bsize, &row, &col);
#if CONFIG_DEBUG
    sub8x8_set_val(cfl, row, col, tx_size_high_unit[tx_size],
                   tx_size_wide_unit[tx_size]);
#endif  // CONFIG_DEBUG
  }
#endif
  cfl_store(cfl, dst, pd->dst.stride, row, col, tx_size_wide[tx_size],
            tx_size_high[tx_size]);
}

void cfl_store_block(MACROBLOCKD *const xd, BLOCK_SIZE bsize, TX_SIZE tx_size) {
  CFL_CTX *const cfl = xd->cfl;
  struct macroblockd_plane *const pd = &xd->plane[AOM_PLANE_Y];
  int row = 0;
  int col = 0;
#if CONFIG_CHROMA_SUB8X8
  bsize = AOMMAX(BLOCK_4X4, bsize);
  if (bsize < BLOCK_8X8) {
    sub8x8_adjust_offset(cfl, bsize, &row, &col);
#if CONFIG_DEBUG
    const int val_high = block_size_high[bsize] / block_size_high[BLOCK_4X4];
    const int val_wide = block_size_wide[bsize] / block_size_wide[BLOCK_4X4];
    sub8x8_set_val(cfl, row, col, val_high, val_wide);
#endif  // CONFIG_DEBUG
  }
#endif  // CONFIG_CHROMA_SUB8X8
  const int width = max_intra_block_width(xd, bsize, AOM_PLANE_Y, tx_size);
  const int height = max_intra_block_height(xd, bsize, AOM_PLANE_Y, tx_size);
  cfl_store(cfl, pd->dst.buf, pd->dst.stride, row, col, width, height);
}

void cfl_compute_parameters(MACROBLOCKD *const xd, TX_SIZE tx_size) {
  CFL_CTX *const cfl = xd->cfl;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;

  // Do not call cfl_compute_parameters multiple time on the same values.
  assert(cfl->are_parameters_computed == 0);

#if CONFIG_CHROMA_SUB8X8
  const BLOCK_SIZE plane_bsize = AOMMAX(
      BLOCK_4X4, get_plane_block_size(mbmi->sb_type, &xd->plane[AOM_PLANE_U]));
#if CONFIG_DEBUG
  if (mbmi->sb_type < BLOCK_8X8) {
    const int val_high =
        block_size_high[BLOCK_8X8] / block_size_high[BLOCK_4X4];
    const int val_wide =
        block_size_wide[BLOCK_8X8] / block_size_wide[BLOCK_4X4];
    for (int val_r = 0; val_r < val_high; val_r++) {
      for (int val_c = 0; val_c < val_wide; val_c++) {
        assert(cfl->sub8x8_val[(val_r * val_wide) + val_c] == 1);
      }
    }
    cfl_clear_sub8x8_val(cfl);
  }
#endif  // CONFIG_DEBUG
#else
  const BLOCK_SIZE plane_bsize =
      get_plane_block_size(mbmi->sb_type, &xd->plane[AOM_PLANE_U]);
#endif
  // AOM_PLANE_U is used, but both planes will have the same sizes.
  cfl->uv_width = max_intra_block_width(xd, plane_bsize, AOM_PLANE_U, tx_size);
  cfl->uv_height =
      max_intra_block_height(xd, plane_bsize, AOM_PLANE_U, tx_size);

#if CONFIG_DEBUG
  if (mbmi->sb_type >= BLOCK_8X8) {
    assert(cfl->ac_width <= cfl->uv_width << cfl->subsampling_x);
    assert(cfl->ac_height <= cfl->uv_height << cfl->subsampling_y);
  }
#endif  // CONFIG_DEBUG

  cfl_dc_pred(xd, plane_bsize);
  cfl_pad_ac_con(cfl);
  cfl->are_parameters_computed = 1;
}
