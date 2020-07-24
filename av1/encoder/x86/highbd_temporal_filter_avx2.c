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

#include <assert.h>
#include <immintrin.h>

#include "config/av1_rtcd.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/temporal_filter.h"

#define SSE_STRIDE (BW + 2)

DECLARE_ALIGNED(32, static const uint32_t, sse_bytemask[4][8]) = {
  { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0, 0 },
  { 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0 },
  { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0 },
  { 0, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF }
};

static AOM_FORCE_INLINE void get_squared_error_16x16_avx2(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height,
    uint32_t *frame_sse, const unsigned int sse_stride) {
  (void)block_width;
  const uint16_t *src1 = CONVERT_TO_SHORTPTR(frame1);
  const uint16_t *src2 = CONVERT_TO_SHORTPTR(frame2);
  uint32_t *dst = frame_sse;
  for (int i = 0; i < block_height; i++) {
    __m256i v_src1 = _mm256_loadu_si256((__m256i *)src1);
    __m256i v_src2 = _mm256_loadu_si256((__m256i *)src2);
    __m256i v_diff = _mm256_sub_epi16(v_src1, v_src2);
    __m256i v_mullo = _mm256_mullo_epi16(v_diff, v_diff);
    __m256i v_mulhi = _mm256_mulhi_epi16(v_diff, v_diff);

    __m256i v_lo = _mm256_unpacklo_epi16(v_mullo, v_mulhi);
    __m256i v_hi = _mm256_unpackhi_epi16(v_mullo, v_mulhi);
    __m256i diff_lo =
        _mm256_inserti128_si256(v_lo, _mm256_extracti128_si256(v_hi, 0), 1);
    __m256i diff_hi =
        _mm256_inserti128_si256(v_hi, _mm256_extracti128_si256(v_lo, 1), 0);

    _mm256_storeu_si256((__m256i *)dst, diff_lo);
    dst += 8;
    _mm256_storeu_si256((__m256i *)dst, diff_hi);

    src1 += stride, src2 += stride2;
    dst += sse_stride - 8;
  }
}

static AOM_FORCE_INLINE void get_squared_error_32x32_avx2(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height,
    uint32_t *frame_sse, const unsigned int sse_stride) {
  (void)block_width;
  const uint16_t *src1 = CONVERT_TO_SHORTPTR(frame1);
  const uint16_t *src2 = CONVERT_TO_SHORTPTR(frame2);
  uint32_t *dst = frame_sse;
  for (int i = 0; i < block_height; i++) {
    __m256i v_src1 = _mm256_loadu_si256((__m256i *)src1);
    __m256i v_src2 = _mm256_loadu_si256((__m256i *)src2);
    __m256i v_diff = _mm256_sub_epi16(v_src1, v_src2);
    __m256i v_mullo = _mm256_mullo_epi16(v_diff, v_diff);
    __m256i v_mulhi = _mm256_mulhi_epi16(v_diff, v_diff);

    __m256i v_lo = _mm256_unpacklo_epi16(v_mullo, v_mulhi);
    __m256i v_hi = _mm256_unpackhi_epi16(v_mullo, v_mulhi);
    __m256i diff_lo =
        _mm256_inserti128_si256(v_lo, _mm256_extracti128_si256(v_hi, 0), 1);
    __m256i diff_hi =
        _mm256_inserti128_si256(v_hi, _mm256_extracti128_si256(v_lo, 1), 0);

    _mm256_storeu_si256((__m256i *)dst, diff_lo);
    dst += 8;
    _mm256_storeu_si256((__m256i *)dst, diff_hi);

    src1 += 16;
    src2 += 16;
    v_src1 = _mm256_loadu_si256((__m256i *)src1);
    v_src2 = _mm256_loadu_si256((__m256i *)src2);
    v_diff = _mm256_sub_epi16(v_src1, v_src2);
    v_mullo = _mm256_mullo_epi16(v_diff, v_diff);
    v_mulhi = _mm256_mulhi_epi16(v_diff, v_diff);

    v_lo = _mm256_unpacklo_epi16(v_mullo, v_mulhi);
    v_hi = _mm256_unpackhi_epi16(v_mullo, v_mulhi);
    diff_lo =
        _mm256_inserti128_si256(v_lo, _mm256_extracti128_si256(v_hi, 0), 1);
    diff_hi =
        _mm256_inserti128_si256(v_hi, _mm256_extracti128_si256(v_lo, 1), 0);

    dst += 8;
    _mm256_storeu_si256((__m256i *)dst, diff_lo);
    dst += 8;
    _mm256_storeu_si256((__m256i *)dst, diff_hi);

    src1 += stride - 16;
    src2 += stride2 - 16;
    dst += sse_stride - 24;
  }
}

static AOM_FORCE_INLINE __m256i xx_load_and_pad(uint32_t *src, int col,
                                                int block_width) {
  __m256i v256tmp = _mm256_loadu_si256((__m256i *)(src));
  if (col == 0) {
    // For the first column, replicate the first element twice to the left
    v256tmp = _mm256_shuffle_epi32(v256tmp, 0x40);
    v256tmp = _mm256_inserti128_si256(v256tmp,
                                      _mm_loadu_si128((__m128i *)(src + 2)), 1);
  }
  if (col == block_width - 4) {
    // For the last column, replicate the last element twice to the right
    v256tmp = _mm256_shuffle_epi32(v256tmp, 0x54);
    v256tmp =
        _mm256_inserti128_si256(v256tmp, _mm_loadu_si128((__m128i *)(src)), 0);
  }
  return v256tmp;
}

static AOM_FORCE_INLINE int32_t xx_mask_and_hadd(__m256i vsum, int i) {
  // Mask the required 5 values inside the vector
  __m256i vtmp = _mm256_and_si256(vsum, *(__m256i *)sse_bytemask[i]);
  __m128i v128a, v128b;
  // Extract 256b as two 128b registers A and B
  v128a = _mm256_castsi256_si128(vtmp);
  v128b = _mm256_extracti128_si256(vtmp, 1);
  // A = [A0+B0, A1+B1, A2+B2, A3+B3]
  v128a = _mm_add_epi32(v128a, v128b);
  // B = [A2+B2, A3+B3, 0, 0]
  v128b = _mm_srli_si128(v128a, 8);
  // A = [A0+B0+A2+B2, A1+B1+A3+B3, X, X]
  v128a = _mm_add_epi32(v128a, v128b);
  // B = [A1+B1+A3+B3, 0, 0, 0]
  v128b = _mm_srli_si128(v128a, 4);
  // A = [A0+B0+A2+B2+A1+B1+A3+B3, X, X, X]
  v128a = _mm_add_epi32(v128a, v128b);
  return _mm_extract_epi32(v128a, 0);
}

static void highbd_apply_temporal_filter(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height,
    const int min_frame_size, const double sigma, const MV *subblock_mvs,
    const int *subblock_mses, const int q_factor, const int filter_strength,
    unsigned int *accumulator, uint16_t *count, uint32_t *luma_sq_error,
    uint32_t *chroma_sq_error, int plane, int ss_x_shift, int ss_y_shift,
    int bd) {
  assert(((block_width == 32) && (block_height == 32)) ||
         ((block_width == 16) && (block_height == 16)));
  if (plane > PLANE_TYPE_Y) assert(chroma_sq_error != NULL);

  uint32_t acc_5x5_sse[BH][BW];
  uint32_t *frame_sse =
      (plane == PLANE_TYPE_Y) ? luma_sq_error : chroma_sq_error;

  if (block_width == 32) {
    get_squared_error_32x32_avx2(frame1, stride, frame2, stride2, block_width,
                                 block_height, frame_sse, SSE_STRIDE);
  } else {
    get_squared_error_16x16_avx2(frame1, stride, frame2, stride2, block_width,
                                 block_height, frame_sse, SSE_STRIDE);
  }

  __m256i vsrc[5];

  const double n_decay = 0.5 + log(2 * sigma + 5.0);
  const double q_decay =
      CLIP(pow((double)q_factor / TF_Q_DECAY_THRESHOLD, 2), 1e-5, 1);
  const double s_decay =
      CLIP(pow((double)filter_strength / TF_STRENGTH_THRESHOLD, 2), 1e-5, 1);

  // Traverse 4 columns at a time
  // First and last columns will require padding
  for (int col = 0; col < block_width; col += 4) {
    uint32_t *src = (col) ? frame_sse + col - 2 : frame_sse;

    // Load and pad(for first and last col) 3 rows from the top
    for (int i = 2; i < 5; i++) {
      vsrc[i] = xx_load_and_pad(src, col, block_width);
      src += SSE_STRIDE;
    }

    // Copy first row to first 2 vectors
    vsrc[0] = vsrc[2];
    vsrc[1] = vsrc[2];

    for (int row = 0; row < block_height; row++) {
      __m256i vsum = _mm256_setzero_si256();

      // Add 5 consecutive rows
      for (int i = 0; i < 5; i++) {
        vsum = _mm256_add_epi32(vsum, vsrc[i]);
      }

      // Push all elements by one element to the top
      for (int i = 0; i < 4; i++) {
        vsrc[i] = vsrc[i + 1];
      }

      // Load next row to the last element
      if (row <= block_width - 4) {
        vsrc[4] = xx_load_and_pad(src, col, block_width);
        src += SSE_STRIDE;
      } else {
        vsrc[4] = vsrc[3];
      }

      // Accumulate the sum horizontally
      for (int i = 0; i < 4; i++) {
        acc_5x5_sse[row][col + i] = xx_mask_and_hadd(vsum, i);
      }
    }
  }

  uint16_t *frame2s = CONVERT_TO_SHORTPTR(frame2);

  for (int i = 0, k = 0; i < block_height; i++) {
    for (int j = 0; j < block_width; j++, k++) {
      const int pixel_value = frame2s[i * stride2 + j];

      int diff_sse = acc_5x5_sse[i][j];
      int num_ref_pixels = TF_WINDOW_LENGTH * TF_WINDOW_LENGTH;

      // Filter U-plane and V-plane using Y-plane. This is because motion
      // search is only done on Y-plane, so the information from Y-plane will
      // be more accurate.
      if (plane != PLANE_TYPE_Y) {
        for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
          for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
            const int yy = (i << ss_y_shift) + ii;  // Y-coord on Y-plane.
            const int xx = (j << ss_x_shift) + jj;  // X-coord on Y-plane.
            diff_sse += luma_sq_error[yy * SSE_STRIDE + xx];
            ++num_ref_pixels;
          }
        }
      }

      // Scale down the difference for high bit depth input.
      diff_sse >>= (bd - 8) * (bd - 8);

      const double window_error = (double)(diff_sse) / num_ref_pixels;
      const int subblock_idx =
          (i >= block_height / 2) * 2 + (j >= block_width / 2);
      const double block_error = (double)subblock_mses[subblock_idx];
      const double combined_error =
          (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) /
          (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / TF_SEARCH_ERROR_NORM_WEIGHT;

      const MV mv = subblock_mvs[subblock_idx];
      const double distance = sqrt(pow(mv.row, 2) + pow(mv.col, 2));
      const double distance_threshold =
          (double)AOMMAX(min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
      const double d_factor = AOMMAX(distance / distance_threshold, 1);

      const double scaled_error =
          AOMMIN(combined_error * d_factor / n_decay / q_decay / s_decay, 7);
      const int weight = (int)(exp(-scaled_error) * TF_WEIGHT_SCALE);

      count[k] += weight;
      accumulator[k] += weight * pixel_value;
    }
  }
}

void av1_highbd_apply_temporal_filter_avx2(
    const YV12_BUFFER_CONFIG *frame_to_filter, const MACROBLOCKD *mbd,
    const BLOCK_SIZE block_size, const int mb_row, const int mb_col,
    const int num_planes, const double *noise_levels, const MV *subblock_mvs,
    const int *subblock_mses, const int q_factor, const int filter_strength,
    const uint8_t *pred, uint32_t *accum, uint16_t *count) {
  const int is_high_bitdepth = frame_to_filter->flags & YV12_FLAG_HIGHBITDEPTH;
  assert(block_size == BLOCK_32X32 && "Only support 32x32 block with avx2!");
  assert(TF_WINDOW_LENGTH == 5 && "Only support window length 5 with avx2!");
  assert(num_planes >= 1 && num_planes <= MAX_MB_PLANE);
  (void)is_high_bitdepth;

  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  const int min_frame_size = AOMMIN(frame_height, frame_width);
  uint32_t luma_sq_error[SSE_STRIDE * BH];
  uint32_t *chroma_sq_error =
      (num_planes > 0)
          ? (uint32_t *)aom_malloc(SSE_STRIDE * BH * sizeof(uint32_t))
          : NULL;

  for (int plane = 0; plane < num_planes; ++plane) {
    const uint32_t plane_h = mb_height >> mbd->plane[plane].subsampling_y;
    const uint32_t plane_w = mb_width >> mbd->plane[plane].subsampling_x;
    const uint32_t frame_stride = frame_to_filter->strides[plane == 0 ? 0 : 1];
    const int frame_offset = mb_row * plane_h * frame_stride + mb_col * plane_w;

    const uint8_t *ref = frame_to_filter->buffers[plane] + frame_offset;
    const int ss_x_shift =
        mbd->plane[plane].subsampling_x - mbd->plane[0].subsampling_x;
    const int ss_y_shift =
        mbd->plane[plane].subsampling_y - mbd->plane[0].subsampling_y;

    highbd_apply_temporal_filter(
        ref, frame_stride, pred + mb_pels * plane, plane_w, plane_w, plane_h,
        min_frame_size, noise_levels[plane], subblock_mvs, subblock_mses,
        q_factor, filter_strength, accum + mb_pels * plane,
        count + mb_pels * plane, luma_sq_error, chroma_sq_error, plane,
        ss_x_shift, ss_y_shift, mbd->bd);
  }
  if (chroma_sq_error != NULL) aom_free(chroma_sq_error);
}
