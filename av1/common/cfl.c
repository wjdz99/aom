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
  memset(&cfl->y_pix, 0, sizeof(uint8_t) * MAX_SB_SQUARE);
  cfl->subsampling_x = cm->subsampling_x;
  cfl->subsampling_y = cm->subsampling_y;
}

// CfL computes its own block-level DC_PRED. This is required to compute both
// alpha_cb and alpha_cr before the prediction are computed.
void cfl_dc_pred(MACROBLOCKD *xd, BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  const struct macroblockd_plane *const pd_u = &xd->plane[AOM_PLANE_U];
  const struct macroblockd_plane *const pd_v = &xd->plane[AOM_PLANE_V];

  const uint8_t *const dst_u = pd_u->dst.buf;
  const uint8_t *const dst_v = pd_v->dst.buf;

  const int dst_u_stride = pd_u->dst.stride;
  const int dst_v_stride = pd_v->dst.stride;

  const int block_width = (plane_bsize != BLOCK_INVALID)
                              ? block_size_wide[plane_bsize]
                              : tx_size_wide[tx_size];
  const int block_height = (plane_bsize != BLOCK_INVALID)
                               ? block_size_high[plane_bsize]
                               : tx_size_high[tx_size];

  // Number of pixel on the top and left borders.
  const double num_pel = block_width + block_height;

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
    for (int i = 0; i < block_width; i++) {
      sum_u += dst_u[-dst_u_stride + i];
      sum_v += dst_v[-dst_v_stride + i];
    }
  } else {
    sum_u = block_width * 127;
    sum_v = block_width * 127;
  }

#if CONFIG_CHROMA_SUB8X8
  if (xd->chroma_left_available && xd->mb_to_bottom_edge >= 0) {
#else
  if (xd->left_available && xd->mb_to_bottom_edge >= 0) {
#endif
    for (int i = 0; i < block_height; i++) {
      sum_u += dst_u[i * dst_u_stride - 1];
      sum_v += dst_v[i * dst_v_stride - 1];
    }
  } else {
    sum_u += block_height * 129;
    sum_v += block_height * 129;
  }

  xd->cfl->dc_pred[CFL_PRED_U] = sum_u / num_pel;
  xd->cfl->dc_pred[CFL_PRED_V] = sum_v / num_pel;
}

// Predict the current transform block using CfL.
void cfl_predict_block(const CFL_CTX *cfl, uint8_t *dst, int dst_stride,
                       int row, int col, TX_SIZE tx_size, double dc_pred,
                       double alpha) {
  const int width = tx_size_wide[tx_size];
  const int height = tx_size_high[tx_size];

  const double y_avg = cfl_load(cfl, dst, dst_stride, row, col, width, height);

  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      dst[i] = (uint8_t)(alpha * (dst[i] - y_avg) + dc_pred + 0.5);
    }
    dst += dst_stride;
  }
}

static inline void copy_block(uint8_t *dst, int dst_stride, const uint8_t *src,
                              int src_stride, int width, int height) {
  for (int j = 0; j < height; j++) {
    memcpy(dst, src, width);
    dst += dst_stride;
    src += src_stride;
  }
}

static inline void copy_block_420(uint8_t *dst, int dst_stride,
                                  const uint8_t *src, int src_stride, int width,
                                  int height) {
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      int t = i << 1;
      int b = t + src_stride;
      dst[i] = OD_SHR_ROUND(src[t] + src[t + 1]         // Top row
                                + src[b] + src[b + 1],  // Bottom row
                            2);
    }
    dst += dst_stride;
    src += src_stride << 1;
  }
}

static inline void pad_block_col(uint8_t *block, int block_stride, int width,
                                 int height, int diff_width) {
  block += width - diff_width;
  for (int j = 0; j < height; j++) {
    memset(block, block[-1], diff_width);
    block += block_stride;
  }
}

static inline void pad_block_row(uint8_t *block, int block_stride, int width,
                                 int height, int diff_height) {
  block += (height - diff_height) * block_stride;
  for (int j = 0; j < diff_height; j++) {
    memcpy(block, block - block_stride, width);
    block += block_stride;
  }
}

static inline double block_average(uint8_t *block, int block_stride, int width,
                                   int height) {
  int sum = 0;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      sum += block[i];
    }
    block += block_stride;
  }
  return sum / (double)(width * height);
}

void cfl_store(CFL_CTX *cfl, const uint8_t *input, int input_stride, int row,
               int col, TX_SIZE tx_size) {
  const int tx_width = tx_size_wide[tx_size];
  const int tx_height = tx_size_high[tx_size];
  const int tx_off_log2 = tx_size_wide_log2[0];

  // Store the input into the CfL pixel buffer
  uint8_t *y_pix = &cfl->y_pix[(row * MAX_SB_SIZE + col) << tx_off_log2];

  // Check that we remain inside the pixel buffer.
  assert(MAX_SB_SIZE * (row + tx_height - 1) + col + tx_width - 1 <
         MAX_SB_SQUARE);

  copy_block(y_pix, MAX_SB_SIZE, input, input_stride, tx_width, tx_height);

  // Store the surface of the pixel buffer that was written to, this way we
  // can manage chroma overrun (e.g. when the chroma surfaces goes beyond the
  // frame boundary)
  if (col == 0 && row == 0) {
    cfl->y_width = tx_width;
    cfl->y_height = tx_height;
  } else {
    cfl->y_width = OD_MAXI((col << tx_off_log2) + tx_width, cfl->y_width);
    cfl->y_height = OD_MAXI((row << tx_off_log2) + tx_height, cfl->y_height);
  }
}

// Load from the CfL pixel buffer into output
double cfl_load(const CFL_CTX *cfl, uint8_t *output, int output_stride, int row,
                int col, int width, int height) {
  const int sub_x = cfl->subsampling_x;
  const int sub_y = cfl->subsampling_y;
  const int off_log2 = tx_size_wide_log2[0];

  const uint8_t *y_pix;

  // TODO(ltrudeau) add support for 4:2:2
  if (sub_y == 0 && sub_x == 0) {
    y_pix = &cfl->y_pix[(row * MAX_SB_SIZE + col) << off_log2];
    // For 4:4:4, match pixels 1 to 1
    copy_block(output, output_stride, y_pix, MAX_SB_SIZE, width, height);
  } else if (sub_y == 1 && sub_x == 1) {
    y_pix = &cfl->y_pix[(row * MAX_SB_SIZE + col) << (off_log2 + sub_y)];
    copy_block_420(output, output_stride, y_pix, MAX_SB_SIZE, width, height);
  } else {
    assert(0);  // Unsupported chroma subsampling
  }
  // Due to frame boundary issues, it is possible that the total area of
  // covered by Chroma exceeds that of Luma. When this happens, we write over
  // the broken data by repeating the last columns and/or rows.
  //
  // Note that in order to manage the case where both rows and columns
  // overrun,
  // we apply rows first. This way, when the rows overrun the bottom of the
  // frame, the columns will be copied over them.
  const int uv_width = (col << off_log2) + width;
  const int uv_height = (row << off_log2) + height;

  const int diff_width = uv_width - (cfl->y_width >> sub_x);
  const int diff_height = uv_height - (cfl->y_height >> sub_y);

  if (diff_width > 0) {
    pad_block_col(output, output_stride, width, height, diff_width);
  }

  if (diff_height > 0) {
    pad_block_row(output, output_stride, width, height, diff_height);
  }

  return block_average(output, output_stride, width, height);
}
