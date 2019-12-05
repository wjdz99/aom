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
#include <immintrin.h>  // AVX2

#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"
#include "av1/encoder/cost.h"
#include "av1/encoder/mcomp.h"

// The function calculates 8 mvsad_err_cost instances if at least 1 of 8 sads in
// sad_array is  < bestsad
// then if at least 1 of resuting errors is < bestsad it adds
// the resulting error to sad_array. The return value shows the necessity of min
// value in sad array calculation: 1 - we need, 0 - no we needn't
int mvsad_err_cost_avx2(MACROBLOCK *x, const search_site *ss, const MV *this_mv,
                        int sad_per_bit, unsigned int bestsad,
                        unsigned int *sad_array) {
  __m256i sad_256 = _mm256_load_si256(
      (__m256i *)sad_array);  // ss[i].mv.row, ss[i].mv.col, offset etc
  const __m256i bestsad256 = _mm256_set1_epi32(bestsad);
  __m256i mask = _mm256_cmpgt_epi32(bestsad256, sad_256);
  int needsad_mask = _mm256_movemask_epi8(mask);
  if (needsad_mask) {
    __m256i ss_mv_256 = _mm256_loadu_si256(
        (__m256i *)ss);  // ss[i].mv.row, ss[i].mv.col, offset etc
    __m256i ss_mv_256_1 = _mm256_loadu_si256(
        (__m256i *)(ss + 4));  // ss[i+4].mv.row, ss[i+4].mv.col,
                               // offset etc best_mv row,col, etc

    // this_mv.row, this_mv.col,  etc
    const __m256i this_mv_256 = _mm256_set1_epi32(*(int *)this_mv);

    ss_mv_256 = _mm256_shuffle_epi32(ss_mv_256, _MM_SHUFFLE(3, 1, 2, 0));
    ss_mv_256_1 = _mm256_shuffle_epi32(ss_mv_256_1, _MM_SHUFFLE(3, 1, 2, 0));
    ss_mv_256 = _mm256_unpacklo_epi64(
        ss_mv_256, ss_mv_256_1);  // ss[i].mv.row, ss[i].mv.col, etc
    ss_mv_256 = _mm256_permute4x64_epi64(ss_mv_256, _MM_SHUFFLE(3, 1, 2, 0));
    ss_mv_256 = _mm256_sub_epi16(ss_mv_256, this_mv_256);

    __m256i diff_256 = _mm256_slli_epi16(ss_mv_256, 3);  // row,col, etc
    DECLARE_ALIGNED(32, int16_t, diff_array[16]);
    _mm256_store_si256((__m256i *)diff_array, diff_256);

    // (!!mv->col) | ((!!mv->row) << 1)
    __m256i zero = _mm256_setzero_si256();
    __m256i get_mv_joint_256 =
        _mm256_srli_epi16(_mm256_xor_si256(_mm256_cmpeq_epi16(diff_256, zero),
                                           _mm256_cmpeq_epi16(zero, zero)),
                          15);  // (!!mv->row), (!!mv->col) etc
    get_mv_joint_256 = _mm256_slli_epi16(
        get_mv_joint_256, 1);  // (!!mv->row) << 1), (!!mv->col) << 1) etc
    get_mv_joint_256 = _mm256_or_si256(
        get_mv_joint_256,
        _mm256_srli_epi32(get_mv_joint_256,
                          17));  // ((!!mv->row) << 1) moved to col positions
    DECLARE_ALIGNED(32, uint16_t, tmp[16]) = { 0xffff, 0, 0xffff, 0, 0xffff, 0,
                                               0xffff, 0, 0xffff, 0, 0xffff, 0,
                                               0xffff, 0, 0xffff, 0 };
    get_mv_joint_256 = _mm256_and_si256(
        get_mv_joint_256,
        *(__m256i *)tmp);  // 4 16 bit results  for each 128 bit
    __m256i nmv_vec_cost256 = _mm256_broadcastsi128_si256(
        _mm_loadu_si128((__m128i *)x->nmv_vec_cost));
    nmv_vec_cost256 =
        _mm256_permutevar8x32_epi32(nmv_vec_cost256, get_mv_joint_256);

    __m256i mv_row_256, mv_col_256;
    mv_row_256 = _mm256_setr_epi32(
        x->mv_cost_stack[0][diff_array[0]], x->mv_cost_stack[0][diff_array[2]],
        x->mv_cost_stack[0][diff_array[4]], x->mv_cost_stack[0][diff_array[6]],
        x->mv_cost_stack[0][diff_array[8]], x->mv_cost_stack[0][diff_array[10]],
        x->mv_cost_stack[0][diff_array[12]],
        x->mv_cost_stack[0][diff_array[14]]);
    mv_col_256 = _mm256_setr_epi32(
        x->mv_cost_stack[1][diff_array[1]], x->mv_cost_stack[1][diff_array[3]],
        x->mv_cost_stack[1][diff_array[5]], x->mv_cost_stack[1][diff_array[7]],
        x->mv_cost_stack[1][diff_array[9]], x->mv_cost_stack[1][diff_array[11]],
        x->mv_cost_stack[1][diff_array[13]],
        x->mv_cost_stack[1][diff_array[15]]);

    __m256i mv_cost_256 = _mm256_add_epi32(
        _mm256_add_epi32(nmv_vec_cost256, mv_row_256), mv_col_256);
    // multiplication by sad_per_bit, faster tthis way
    __m256i mv_cost_spb_256 = mv_cost_256;
    for (int s = 0; s < sad_per_bit - 1; s++) {
      mv_cost_spb_256 = _mm256_add_epi32(mv_cost_256, mv_cost_spb_256);
    }
    const __m256i mv_shift_const_256 =
        _mm256_set1_epi32(1 << (AV1_PROB_COST_SHIFT - 1));
    __m256i mv_err_cost_256 =
        _mm256_srli_epi32(_mm256_add_epi32(mv_cost_spb_256, mv_shift_const_256),
                          AV1_PROB_COST_SHIFT);
    mask = _mm256_cmpgt_epi32(bestsad256, mv_err_cost_256);
    if (_mm256_movemask_epi8(mask)) {
      mv_err_cost_256 = _mm256_add_epi32(mv_err_cost_256, sad_256);
      _mm256_store_si256((__m256i *)sad_array, mv_err_cost_256);
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

void diamond_search_all_in_bestsad_avx2(int i, MACROBLOCK *x,
                                        int searches_per_step,
                                        const search_site *ss, const MV this_mv,
                                        int sad_per_bit,
                                        unsigned int *sad_array,
                                        unsigned int *bestsad, int *best_site) {
  if (mvsad_err_cost_avx2(x, (ss + i), &this_mv, sad_per_bit, *bestsad,
                          sad_array)) {
    for (int t = 0; t < searches_per_step; t++) {
      if (sad_array[t] < *bestsad) {
        *bestsad = sad_array[t];
        *best_site = i + t;
      }
    }
  }
}

void aom_sadMxNx4d_avx2(int M, int N, const uint8_t *src, int src_stride,
                        const uint8_t *const ref[4], int ref_stride,
                        uint32_t res[4]) {
  __m256i src_reg, ref0_reg, ref1_reg, ref2_reg, ref3_reg;
  __m256i sum_ref0, sum_ref1, sum_ref2, sum_ref3;
  int i, j;
  const uint8_t *ref0, *ref1, *ref2, *ref3;
  static int nn = 0;
  // printf("sadMxNx4d_avx2 %d \n", nn);
  nn++;
  ref0 = ref[0];
  ref1 = ref[1];
  ref2 = ref[2];
  ref3 = ref[3];
  sum_ref0 = _mm256_setzero_si256();
  sum_ref2 = _mm256_setzero_si256();
  sum_ref1 = _mm256_setzero_si256();
  sum_ref3 = _mm256_setzero_si256();

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j += 32) {
      // load src and all refs
      src_reg = _mm256_loadu_si256((const __m256i *)(src + j));
      ref0_reg = _mm256_loadu_si256((const __m256i *)(ref0 + j));
      ref1_reg = _mm256_loadu_si256((const __m256i *)(ref1 + j));
      ref2_reg = _mm256_loadu_si256((const __m256i *)(ref2 + j));
      ref3_reg = _mm256_loadu_si256((const __m256i *)(ref3 + j));

      // sum of the absolute differences between every ref-i to src
      ref0_reg = _mm256_sad_epu8(ref0_reg, src_reg);
      ref1_reg = _mm256_sad_epu8(ref1_reg, src_reg);
      ref2_reg = _mm256_sad_epu8(ref2_reg, src_reg);
      ref3_reg = _mm256_sad_epu8(ref3_reg, src_reg);
      // sum every ref-i
      sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
      sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
      sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
      sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);
    }
    src += src_stride;
    ref0 += ref_stride;
    ref1 += ref_stride;
    ref2 += ref_stride;
    ref3 += ref_stride;
  }
  {
    __m128i sum;
    __m256i sum_mlow, sum_mhigh;
    // in sum_ref-i the result is saved in the first 4 bytes
    // the other 4 bytes are zeroed.
    // sum_ref1 and sum_ref3 are shifted left by 4 bytes
    sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
    sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

    // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
    sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
    sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

    // merge every 64 bit from each sum_ref-i
    sum_mlow = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
    sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

    // add the low 64 bit to the high 64 bit
    sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

    // add the low 128 bit to the high 128 bit
    sum = _mm_add_epi32(_mm256_castsi256_si128(sum_mlow),
                        _mm256_extractf128_si256(sum_mlow, 1));

    _mm_storeu_si128((__m128i *)(res), sum);
  }
}

#define sadMxN_avx2(m, n)                                                      \
  void aom_sad##m##x##n##x4d_avx2(const uint8_t *src, int src_stride,          \
                                  const uint8_t *const ref[4], int ref_stride, \
                                  uint32_t res[4]) {                           \
    aom_sadMxNx4d_avx2(m, n, src, src_stride, ref, ref_stride, res);           \
  }

sadMxN_avx2(32, 8);
sadMxN_avx2(32, 16);
sadMxN_avx2(32, 32);
sadMxN_avx2(32, 64);

sadMxN_avx2(64, 16);
sadMxN_avx2(64, 32);
sadMxN_avx2(64, 64);
sadMxN_avx2(64, 128);

sadMxN_avx2(128, 64);
sadMxN_avx2(128, 128);
