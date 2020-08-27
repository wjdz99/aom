/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/common/mfqe.h"
#include "av1/common/resize.h"

#include "av1/encoder/rdopt.h"

#define MSE_MAX_BITS 16
#define MSE_MAX (1 << 16)

// Dynamically allocate memory for the buffer to store the blurred version of
// the given buffer representing a single plane, then return its pointer.
static Y_BUFFER_CONFIG get_gaussian_blur(uint8_t *src, int stride, int h, int w,
                                         int high_bd, int bd) {
  Y_BUFFER_CONFIG dst;
  dst.stride = stride;
  dst.height = h;
  dst.width = w;
  dst.buffer = (uint8_t *)aom_memalign(32, sizeof(uint8_t) * stride * w);

  int w_size;
  int h_size;
  uint8_t *src_ptr;
  uint8_t *dst_ptr;
  for (int row = 0; row < w; row += 32) {
    for (int col = 0; col < h; col += 32) {
      w_size = AOMMIN(32, w - row);
      h_size = AOMMIN(32, h - col);
      src_ptr = src + col * stride + row;
      dst_ptr = dst.buffer + col * stride + row;

      av1_gaussian_blur(src_ptr, stride, w_size, h_size, dst_ptr, high_bd, bd);
    }
  }
  return dst;
}

// Dynamically allocate memory for the buffer to store the resized version of
// the given buffer representing a single plane, then return its pointer.
static Y_BUFFER_CONFIG get_resize_plane(const uint8_t *src, int stride, int h,
                                        int w, int scale) {
  Y_BUFFER_CONFIG dst;
  dst.stride = stride * scale;
  dst.height = h * scale;
  dst.width = w * scale;

  int num_buf_bytes = dst.stride * dst.width;
  dst.buffer = (uint8_t *)aom_memalign(32, sizeof(uint8_t) * num_buf_bytes);
  av1_resize_plane(src, h, w, stride, dst.buffer, dst.height, dst.width,
                   dst.stride);
  return dst;
}

// Returns the mean squared error between the given blocks in two buffers. If
// the row and column parameters are not valid indices, return MSE_MAX.
static double get_mse_block(Y_BUFFER_CONFIG buf1, Y_BUFFER_CONFIG buf2,
                            int16_t mb_row_1, int16_t mb_col_1,
                            int16_t mb_row_2, int16_t mb_col_2,
                            int bsize) {
  // Check if rows and columns are valid, return MSE_MAX if not.
  if ((mb_row_1 < 0) || (mb_col_1 < 0) || (mb_row_2 < 0) || (mb_col_2 < 0) ||
      (mb_row_1 >= buf1.height - bsize) ||
      (mb_row_2 >= buf2.height - bsize) ||
      (mb_col_1 >= buf1.width - bsize) ||
      (mb_col_2 >= buf2.width - bsize))
    return MSE_MAX;

  double mse = 0.0;
  int start_id_1 = buf1.stride * mb_col_1 + mb_row_1;
  int start_id_2 = buf2.stride * mb_col_2 + mb_row_2;

  uint8_t p1;
  uint8_t p2;
  double v1;
  double v2;
  int divisor = bsize * bsize;

  for (int row = 0; row < bsize; ++row) {
    for (int col = 0; col < bsize; ++col) {
      p1 = buf1.buffer[start_id_1 + buf1.stride * col + row];
      p2 = buf2.buffer[start_id_2 + buf2.stride * col + row];
      v1 = (double)p1;
      v2 = (double)p2;
      mse += (v1 - v2) * (v1 - v2) / divisor;
    }
  }
  return mse;
}

// Perform initial grid search to obtain the full-pixel motion vector for every
// block in the current frame, going through diamond points around the block.
static void full_pixel_grid_search(MV_MFQE *mvr, int16_t mb_row, int16_t mb_col,
                                   Y_BUFFER_CONFIG cur, Y_BUFFER_CONFIG refs[],
                                   int bsize) {
  MV_MFQE mvr_best = *mvr;
  double this_mse;
  double best_mse = MSE_MAX;
  int16_t dr;
  int16_t dc;
  int16_t mb_row_ref;
  int16_t mb_col_ref;

  for (int ref_index = 0; ref_index < 3; ++ref_index) {
    for (int point = 0; point < 13; ++point) {
      dr = grid_search_rows[point];
      dc = grid_search_cols[point];
      mb_row_ref = mb_row + mvr->mv.row + dr * bsize;
      mb_col_ref = mb_col + mvr->mv.col + dc * bsize;

      this_mse = get_mse_block(cur, refs[ref_index], mb_row, mb_col, mb_row_ref,
                               mb_col_ref, bsize);

      // Store the motion vector with lowest mean squared error.
      if (this_mse < best_mse) {
        best_mse = this_mse;
        mvr_best.mv.row = mvr->mv.row + dr * bsize;
        mvr_best.mv.col = mvr->mv.col + dc * bsize;
        mvr_best.ref_index = ref_index;
      }
    }
  }

  *mvr = mvr_best;
}

// Perform full pixel motion vector search in the low frequency version of the
// current frame and reference frames.
static void full_pixel_search(MV_MFQE *mvr, int16_t mb_row, int16_t mb_col,
                              Y_BUFFER_CONFIG cur, Y_BUFFER_CONFIG refs[],
                              int bsize) {
  MV_MFQE mvr_init = *mvr;
  full_pixel_grid_search(&mvr_init, mb_row, mb_col, cur, refs, bsize);

  mvr->valid = 0;
  mvr->ref_index = mvr_init.ref_index;

  int mb_row_ref;
  int mb_col_ref;
  int search_size = bsize / 2;
  double this_mse;
  double best_mse = MSE_MAX;
  const double threshold = 0.8;

  for (int dr = -search_size; dr <= search_size; ++dr) {
    for (int dc = -search_size; dc <= search_size; ++dc) {
      mb_row_ref = mb_row + mvr_init.mv.row + dr;
      mb_col_ref = mb_col + mvr_init.mv.col + dc;

      this_mse = get_mse_block(cur, refs[mvr_init.ref_index], mb_row, mb_col,
                               mb_row_ref, mb_col_ref, bsize);

      if (this_mse < best_mse && this_mse < threshold) {
        best_mse = this_mse;
        mvr->valid = 1;
        mvr->mv.row = mvr_init.mv.row + dr;
        mvr->mv.col = mvr_init.mv.col + dc;
      }
    }
  }
}

// Perform finer-grained motion vector search at subpel level, then save the
// updated motion vector in MV_MFQE.
static void sub_pixel_search(MV_MFQE *mvr, int16_t mb_row, int16_t mb_col,
                             Y_BUFFER_CONFIG cur, Y_BUFFER_CONFIG refs[],
                             int bsize, int scale) {
  mb_row *= scale;
  mb_col *= scale;
  bsize *= scale;

  int search_size = scale / 2;
  int mb_row_ref;
  int mb_col_ref;

  double this_mse;
  double best_mse = MSE_MAX;
  MV_MFQE mv_best = *mvr;

  // Search for the nearby search_size pixels in subpel accuracy, which is half
  // the size of the scale. For example, if the image is scaled by a factor of
  // 8, the algorithm will search the adjacent 4 pixels in the resized image.
  for (int dr = -search_size; dr <= search_size; ++dr) {
    for (int dc = -search_size; dc <= search_size; ++dc) {
      mb_row_ref = mb_row + mvr->mv.row * scale + dr;
      mb_col_ref = mb_col + mvr->mv.col * scale + dc;

      this_mse = get_mse_block(cur, refs[mvr->ref_index], mb_row, mb_col,
                               mb_row_ref, mb_col_ref, bsize);

      if (this_mse < best_mse) {
        best_mse = this_mse;
        mv_best.subpel_x_qn = dr;
        mv_best.subpel_y_qn = dc;
      }
    }
  }

  mvr->subpel_x_qn = mv_best.subpel_x_qn;
  mvr->subpel_y_qn = mv_best.subpel_y_qn;
}

// Replace the block in current frame using the block from the reference frame,
// using subpel motion vectors and interpolating the reference block.
static void replace_block_subpel(Y_BUFFER_CONFIG tmp,
                                 Y_BUFFER_CONFIG refs_sub[], MV_MFQE *mvr,
                                 int16_t mb_row, int16_t mb_col, int bsize,
                                 ConvolveParams *conv_params,
                                 int_interpfilters *interp_filters) {
  int16_t mb_row_ref = mb_row + mvr->mv.row;
  int16_t mb_col_ref = mb_col + mvr->mv.col;
  Y_BUFFER_CONFIG ref = refs_sub[mvr->ref_index];

  uint8_t *src = ref.buffer + mb_row_ref * ref.stride + mb_col_ref;
  uint8_t *dst = tmp.buffer + mb_row * tmp.stride + mb_col;
  int src_stride = ref.stride;
  int dst_stride = tmp.stride;

  av1_convolve_2d_facade(src, src_stride, dst, dst_stride, bsize, bsize, bsize,
                         bsize, *interp_filters, mvr->subpel_x_qn, 0,
                         mvr->subpel_y_qn, 0, 0, conv_params, 0);
}

void av1_apply_loop_mfqe(Y_BUFFER_CONFIG *tmp, RefCntBuffer *ref_frames[],
                         int bsize, int scale, int high_bd, int bd) {
  Y_BUFFER_CONFIG tmp_low = get_gaussian_blur(
      tmp->buffer, tmp->stride, tmp->height, tmp->width, high_bd, bd);
  Y_BUFFER_CONFIG tmp_sub = get_resize_plane(tmp_low.buffer, tmp->stride,
                                             tmp->height, tmp->width, scale);

  Y_BUFFER_CONFIG refs_low[3];  // Contains blurred versions of ref frames.
  Y_BUFFER_CONFIG refs_sub[3];  // Contains resized versions of ref frames.
  YV12_BUFFER_CONFIG *ref;
  for (int i = 0; i < 3; i++) {
    ref = &ref_frames[i]->buf;
    refs_low[i] = get_gaussian_blur(ref->y_buffer, ref->y_stride, ref->y_height,
                                    ref->y_width, high_bd, bd);
    refs_sub[i] = get_resize_plane(ref->y_buffer, ref->y_stride, ref->y_height,
                                   ref->y_width, scale);
  }

  // Set up parameters for av1_convolve_2d_facade.
  ConvolveParams conv_params = get_conv_params(0, AOM_PLANE_Y, bd);
  int_interpfilters interp_filters;
  interp_filters.as_filters.x_filter = EIGHTTAP_REGULAR;
  interp_filters.as_filters.y_filter = EIGHTTAP_REGULAR;
  MV_MFQE mvr;

  int16_t num_rows = tmp->height / bsize;
  int16_t num_cols = tmp->width / bsize;

  for (int16_t mb_row = 0; mb_row < num_rows; ++mb_row) {
    for (int16_t mb_col = 0; mb_col < num_cols; ++mb_col) {
      mvr = kZeroMvMFQE;
      full_pixel_search(&mvr, mb_row, mb_col, tmp_low, refs_low, bsize);

      if (!mvr.valid) continue;  // Pass if mse is larger than threshold.
      sub_pixel_search(&mvr, mb_row, mb_col, tmp_sub, refs_sub, bsize,
                       scale);
      replace_block_subpel(*tmp, refs_sub, &mvr, mb_row, mb_col, bsize,
                           &conv_params, &interp_filters);
    }
  }

  // Free dynamically allocated buffers.
  aom_free(tmp_low.buffer);
  aom_free(tmp_sub.buffer);
  for (int i = 0; i < 3; i++) {
    aom_free(refs_low[i].buffer);
    aom_free(refs_sub[i].buffer);
  }
}

void av1_decode_restore_mfqe(AV1_COMMON *cm, int scale, int bsize) {
  assert(cm->use_mfqe);

  YV12_BUFFER_CONFIG *cur = &cm->cur_frame->buf;
  Y_BUFFER_CONFIG cur_frame = { cur->y_buffer, cur->y_stride, cur->y_height,
                                cur->y_width };

  RefCntBuffer *ref_frames[7];
  int num_ref_frames = 0;
  MV_REFERENCE_FRAME ref_frame;
  for (ref_frame = LAST_FRAME; ref_frame < ALTREF_FRAME; ++ref_frame) {
    RefCntBuffer *ref = get_ref_frame_buf(cm, ref_frame);
    if (ref) ref_frames[num_ref_frames++] = ref;
  }
  assert(num_ref_frames >= 3);

  // Assert that pointers to RefCntBuffer are valid, then sort the reference
  // frames based on their base_qindex, from lowest to highest.
  for (int i = 0; i < num_ref_frames; i++) assert(ref_frames[i] != NULL);
  qsort(ref_frames, num_ref_frames, sizeof(RefCntBuffer *), cmpref);

  // Perform In-Loop Multi-Frame Quality Enhancement on tmp.
  av1_apply_loop_mfqe(&cur_frame, ref_frames, bsize, scale,
                      cm->seq_params.use_highbitdepth,
                      cm->seq_params.bit_depth);
}
