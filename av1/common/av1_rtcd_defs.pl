##
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
sub av1_common_forward_decls() {
print <<EOF
/*
 * AV1
 */

#include "aom/aom_integer.h"
#include "aom_dsp/txfm_common.h"
#include "av1/common/common.h"
#include "av1/common/enums.h"
#include "av1/common/quant_common.h"
#include "av1/common/filter.h"
#include "av1/common/convolve.h"
#include "av1/common/av1_txfm.h"
#include "av1/common/odintrin.h"
#include "av1/common/restoration.h"
#include "av1/common/entropymode.h"

struct macroblockd;

/* Encoder forward decls */
struct macroblock;
struct txfm_param;
struct aom_variance_vtable;
struct search_site_config;
struct yv12_buffer_config;
struct NN_CONFIG;
typedef struct NN_CONFIG NN_CONFIG;

struct CNN_CONFIG;
typedef struct CNN_CONFIG CNN_CONFIG;
struct CNN_LAYER_CONFIG;
typedef struct CNN_LAYER_CONFIG CNN_LAYER_CONFIG;
struct CNN_THREAD_DATA;
typedef struct CNN_THREAD_DATA CNN_THREAD_DATA;
struct CNN_BRANCH_CONFIG;
typedef struct CNN_BRANCH_CONFIG CNN_BRANCH_CONFIG;
struct CNN_MULTI_OUT;
typedef struct CNN_MULTI_OUT CNN_MULTI_OUT;

enum { NONE, RELU, SOFTSIGN, SIGMOID } UENUM1BYTE(ACTIVATION);

/* Function pointers return by CfL functions */
typedef void (*cfl_subsample_lbd_fn)(const uint8_t *input, int input_stride,
                                     uint16_t *output_q3);

typedef void (*cfl_subsample_hbd_fn)(const uint16_t *input, int input_stride,
                                     uint16_t *output_q3);

typedef void (*cfl_subtract_average_fn)(const uint16_t *src, int16_t *dst);

typedef void (*cfl_predict_lbd_fn)(const int16_t *src, uint8_t *dst,
                                   int dst_stride, int alpha_q3);

typedef void (*cfl_predict_hbd_fn)(const int16_t *src, uint16_t *dst,
                                   int dst_stride, int alpha_q3, int bd);
EOF
}
forward_decls qw/av1_common_forward_decls/;

# functions that are 64 bit only.
$mmx_x86_64 = $sse2_x86_64 = $ssse3_x86_64 = $avx_x86_64 = $avx2_x86_64 = '';
if ($opts{arch} eq "x86_64") {
  $mmx_x86_64 = 'mmx';
  $sse2_x86_64 = 'sse2';
  $ssse3_x86_64 = 'ssse3';
  $avx_x86_64 = 'avx';
  $avx2_x86_64 = 'avx2';
}

add_proto qw/void av1_convolve_horiz_rs/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const int16_t *x_filters, int x0_qn, int x_step_qn";
specialize qw/av1_convolve_horiz_rs sse4_1/;

add_proto qw/void av1_highbd_convolve_horiz_rs/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const int16_t *x_filters, int x0_qn, int x_step_qn, int bd";
specialize qw/av1_highbd_convolve_horiz_rs sse4_1/;

add_proto qw/void av1_wiener_convolve_add_src/,       "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, const ConvolveParams *conv_params";

add_proto qw/void av1_highbd_wiener_convolve_add_src/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, const ConvolveParams *conv_params, int bd";

specialize qw/av1_wiener_convolve_add_src sse2 avx2 neon/;
specialize qw/av1_highbd_wiener_convolve_add_src avx2 ssse3/;

add_proto qw/void av1_wiener_hp_convolve_add_src/,       "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, const ConvolveParams *conv_params";

add_proto qw/void av1_highbd_wiener_hp_convolve_add_src/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, const ConvolveParams *conv_params, int bd";

# directional intra predictor functions
add_proto qw/void av1_dr_prediction_z1/, "uint8_t *dst, ptrdiff_t stride, int bw, int bh, const uint8_t *above, const uint8_t *left, int upsample_above, int dx, int dy";
specialize qw/av1_dr_prediction_z1 avx2/;
add_proto qw/void av1_dr_prediction_z2/, "uint8_t *dst, ptrdiff_t stride, int bw, int bh, const uint8_t *above, const uint8_t *left, int upsample_above, int upsample_left, int dx, int dy";
specialize qw/av1_dr_prediction_z2 avx2/;
add_proto qw/void av1_dr_prediction_z3/, "uint8_t *dst, ptrdiff_t stride, int bw, int bh, const uint8_t *above, const uint8_t *left, int upsample_left, int dx, int dy";
specialize qw/av1_dr_prediction_z3 avx2/;

# FILTER_INTRA predictor functions
add_proto qw/void av1_filter_intra_predictor/, "uint8_t *dst, ptrdiff_t stride, TX_SIZE tx_size, const uint8_t *above, const uint8_t *left, int mode";
specialize qw/av1_filter_intra_predictor sse4_1/;

# High bitdepth functions

#
# Sub Pixel Filters
#
add_proto qw/void av1_highbd_convolve_copy/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";

add_proto qw/void av1_highbd_convolve_avg/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";

add_proto qw/void av1_highbd_convolve8/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
specialize qw/av1_highbd_convolve8/, "$sse2_x86_64";

add_proto qw/void av1_highbd_convolve8_horiz/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
specialize qw/av1_highbd_convolve8_horiz/, "$sse2_x86_64";

add_proto qw/void av1_highbd_convolve8_vert/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
specialize qw/av1_highbd_convolve8_vert/, "$sse2_x86_64";

#inv txfm
add_proto qw/void av1_inv_txfm_add/, "const tran_low_t *dqcoeff, uint8_t *dst, int stride, const TxfmParam *txfm_param";
specialize qw/av1_inv_txfm_add ssse3 avx2 neon/;

add_proto qw/void av1_highbd_inv_txfm_add/, "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add sse4_1 avx2/;

add_proto qw/void av1_highbd_inv_txfm_add_4x4/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_4x4 sse4_1/;
add_proto qw/void av1_highbd_inv_txfm_add_8x8/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_8x8 sse4_1/;
add_proto qw/void av1_highbd_inv_txfm_add_4x8/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_4x8 sse4_1/;
add_proto qw/void av1_highbd_inv_txfm_add_8x4/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_8x4 sse4_1/;
add_proto qw/void av1_highbd_inv_txfm_add_4x16/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_4x16 sse4_1/;
add_proto qw/void av1_highbd_inv_txfm_add_16x4/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
specialize qw/av1_highbd_inv_txfm_add_16x4 sse4_1/;
if (aom_config("CONFIG_FLEX_PARTITION") eq "yes") {
  add_proto qw/void av1_highbd_inv_txfm_add_4x32/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
  specialize qw/av1_highbd_inv_txfm_add_4x32 sse4_1/;
  add_proto qw/void av1_highbd_inv_txfm_add_32x4/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
  specialize qw/av1_highbd_inv_txfm_add_32x4 sse4_1/;
  add_proto qw/void av1_highbd_inv_txfm_add_4x64/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
  specialize qw/av1_highbd_inv_txfm_add_4x64 sse4_1/;
  add_proto qw/void av1_highbd_inv_txfm_add_64x4/,  "const tran_low_t *input, uint8_t *dest, int stride, const TxfmParam *txfm_param";
  specialize qw/av1_highbd_inv_txfm_add_64x4 sse4_1/;
}

add_proto qw/void av1_highbd_iwht4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
add_proto qw/void av1_highbd_iwht4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";

  add_proto qw/void av1_inv_txfm2d_add_4x8/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_8x4/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_8x16/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_16x8/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_16x32/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_32x16/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_4x4/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_inv_txfm2d_add_4x4 sse4_1/;
  add_proto qw/void av1_inv_txfm2d_add_8x8/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_inv_txfm2d_add_8x8 sse4_1/;
  add_proto qw/void av1_inv_txfm2d_add_16x16/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_32x32/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";

  add_proto qw/void av1_inv_txfm2d_add_64x64/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_32x64/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_64x32/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_16x64/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_64x16/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";

  add_proto qw/void av1_inv_txfm2d_add_4x16/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_16x4/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_8x32/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_inv_txfm2d_add_32x8/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  if (aom_config("CONFIG_FLEX_PARTITION") eq "yes") {
    add_proto qw/void av1_inv_txfm2d_add_4x32/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_inv_txfm2d_add_32x4/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_inv_txfm2d_add_8x64/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_inv_txfm2d_add_64x8/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_inv_txfm2d_add_4x64/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_inv_txfm2d_add_64x4/, "const int32_t *input, uint16_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  }

# directional intra predictor functions
add_proto qw/void av1_highbd_dr_prediction_z1/, "uint16_t *dst, ptrdiff_t stride, int bw, int bh, const uint16_t *above, const uint16_t *left, int upsample_above, int dx, int dy, int bd";
specialize qw/av1_highbd_dr_prediction_z1 avx2/;
add_proto qw/void av1_highbd_dr_prediction_z2/, "uint16_t *dst, ptrdiff_t stride, int bw, int bh, const uint16_t *above, const uint16_t *left, int upsample_above, int upsample_left, int dx, int dy, int bd";
specialize qw/av1_highbd_dr_prediction_z2 avx2/;
add_proto qw/void av1_highbd_dr_prediction_z3/, "uint16_t *dst, ptrdiff_t stride, int bw, int bh, const uint16_t *above, const uint16_t *left, int upsample_left, int dx, int dy, int bd";
specialize qw/av1_highbd_dr_prediction_z3 avx2/;

# build compound seg mask functions
add_proto qw/void av1_build_compound_diffwtd_mask/, "uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const uint8_t *src0, int src0_stride, const uint8_t *src1, int src1_stride, int h, int w";
specialize qw/av1_build_compound_diffwtd_mask sse4_1 avx2/;

add_proto qw/void av1_build_compound_diffwtd_mask_highbd/, "uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const uint8_t *src0, int src0_stride, const uint8_t *src1, int src1_stride, int h, int w, int bd";
specialize qw/av1_build_compound_diffwtd_mask_highbd ssse3 avx2/;

add_proto qw/void av1_build_compound_diffwtd_mask_d16/, "uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0, int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w, ConvolveParams *conv_params, int bd";
specialize qw/av1_build_compound_diffwtd_mask_d16 sse4_1 avx2 neon/;

# Helper functions.
add_proto qw/void av1_round_shift_array/, "int32_t *arr, int size, int bit";
specialize "av1_round_shift_array", qw/sse4_1 neon/;

#
# Encoder functions below this point.
#
if (aom_config("CONFIG_AV1_ENCODER") eq "yes") {

  # ENCODEMB INVOKE

  # the transform coefficients are held in 32-bit
  # values, so the assembler code for  av1_block_error can no longer be used.
  add_proto qw/int64_t av1_block_error/, "const tran_low_t *coeff, const tran_low_t *dqcoeff, intptr_t block_size, int64_t *ssz";
  specialize qw/av1_block_error sse2 avx2 neon/;

  if (aom_config("CONFIG_EXTQUANT") eq "yes") {
    add_proto qw/void av1_quantize_fp/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int32_t *zbin_ptr, const int32_t *round_ptr, const int32_t *quant_ptr, const int32_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int32_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    add_proto qw/void av1_quantize_fp_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int32_t *zbin_ptr, const int32_t *round_ptr, const int32_t *quant_ptr, const int32_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int32_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    add_proto qw/void av1_quantize_fp_64x64/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int32_t *zbin_ptr, const int32_t *round_ptr, const int32_t *quant_ptr, const int32_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int32_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  } else {
    add_proto qw/void av1_quantize_fp/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    specialize qw/av1_quantize_fp sse2 avx2 neon/;

    add_proto qw/void av1_quantize_fp_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    specialize qw/av1_quantize_fp_32x32 avx2/;

    add_proto qw/void av1_quantize_fp_64x64/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    specialize qw/av1_quantize_fp_64x64 avx2/;
  }

  # fdct functions

  add_proto qw/void av1_fwht4x4/, "const int16_t *input, tran_low_t *output, int stride";

  #fwd txfm
  add_proto qw/void av1_lowbd_fwd_txfm/, "const int16_t *src_diff, tran_low_t *coeff, int diff_stride, TxfmParam *txfm_param";
  specialize qw/av1_lowbd_fwd_txfm sse2 sse4_1 avx2/;

  add_proto qw/void av1_fwd_txfm2d_4x8/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_4x8 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_8x4/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_8x4 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_8x16/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_8x16 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_16x8/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_16x8 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_16x32/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_16x32 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_32x16/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_32x16 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_4x16/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_4x16 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_16x4/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_16x4 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_8x32/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_8x32 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_32x8/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_32x8 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_4x4/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_4x4 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_8x8/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_8x8 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_16x16/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_16x16 sse4_1/;
  add_proto qw/void av1_fwd_txfm2d_32x32/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  specialize qw/av1_fwd_txfm2d_32x32 sse4_1 avx2/;

  add_proto qw/void av1_fwd_txfm2d_64x64/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_fwd_txfm2d_32x64/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_fwd_txfm2d_64x32/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_fwd_txfm2d_16x64/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  add_proto qw/void av1_fwd_txfm2d_64x16/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";

  if (aom_config("CONFIG_NEW_TX64X64") eq "") {
    specialize qw/av1_fwd_txfm2d_64x64 sse4_1 avx2/;
    specialize qw/av1_fwd_txfm2d_32x64 sse4_1/;
    specialize qw/av1_fwd_txfm2d_64x32 sse4_1/;
    specialize qw/av1_fwd_txfm2d_16x64 sse4_1/;
    specialize qw/av1_fwd_txfm2d_64x16 sse4_1/;
  }
  if (aom_config("CONFIG_FLEX_PARTITION") eq "yes") {
    add_proto qw/void av1_fwd_txfm2d_4x32/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_fwd_txfm2d_32x4/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_fwd_txfm2d_8x64/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_fwd_txfm2d_64x8/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_fwd_txfm2d_4x64/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
    add_proto qw/void av1_fwd_txfm2d_64x4/, "const int16_t *input, int32_t *output, int stride, TX_TYPE tx_type, PREDICTION_MODE mode, int bd";
  }

  #
  # Motion search
  #
  add_proto qw/int av1_diamond_search_sad/, "const struct AV1Common *const cm, struct macroblock *x, const struct search_site_config *cfg,  MV *ref_mv, MV *best_mv, int search_param, int sad_per_bit, int *num00, const struct aom_variance_vtable *fn_ptr, const MV *center_mv";

  add_proto qw/int av1_full_range_search/, "const struct macroblock *x, const struct search_site_config *cfg, MV *ref_mv, MV *best_mv, int search_param, int sad_per_bit, int *num00, const struct aom_variance_vtable *fn_ptr, const MV *center_mv";

  if (aom_config("CONFIG_REALTIME_ONLY") ne "yes") {
    add_proto qw/void av1_apply_temporal_filter/, "const uint8_t *y_frame1, int y_stride, const uint8_t *y_pred, int y_buf_stride, const uint8_t *u_frame1, const uint8_t *v_frame1, int uv_stride, const uint8_t *u_pred, const uint8_t *v_pred, int uv_buf_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, int strength, const int *blk_fw, int use_32x32, uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator, uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count";
    specialize qw/av1_apply_temporal_filter sse4_1/;
  }

  if (aom_config("CONFIG_EXTQUANT") eq "yes") {
    add_proto qw/void av1_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int32_t *zbin_ptr, const int32_t *round_ptr, const int32_t *quant_ptr, const int32_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int32_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan, const qm_val_t * qm_ptr, const qm_val_t * iqm_ptr, int log_scale";
  } else {
    add_proto qw/void av1_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan, const qm_val_t * qm_ptr, const qm_val_t * iqm_ptr, int log_scale";
  }

  # ENCODEMB INVOKE

  add_proto qw/int64_t av1_highbd_block_error/, "const tran_low_t *coeff, const tran_low_t *dqcoeff, intptr_t block_size, int64_t *ssz, int bd";
  specialize qw/av1_highbd_block_error sse2 avx2/;

  if (aom_config("CONFIG_REALTIME_ONLY") ne "yes") {
    add_proto qw/void av1_highbd_apply_temporal_filter/, "const uint8_t *yf, int y_stride, const uint8_t *yp, int y_buf_stride, const uint8_t *uf, const uint8_t *vf, int uv_stride, const uint8_t *up, const uint8_t *vp, int uv_buf_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, int strength, const int *blk_fw, int use_32x32, uint32_t *y_accumulator, uint16_t *y_count, uint32_t *u_accumulator, uint16_t *u_count, uint32_t *v_accumulator, uint16_t *v_count";
    specialize qw/av1_highbd_apply_temporal_filter sse4_1/;
  }

  if (aom_config("CONFIG_EXTQUANT") eq "yes") {
    add_proto qw/void av1_highbd_quantize_fp/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int32_t *zbin_ptr, const int32_t *round_ptr, const int32_t *quant_ptr, const int32_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int32_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan, int log_scale";
  } else {
    add_proto qw/void av1_highbd_quantize_fp/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan, int log_scale";
    specialize qw/av1_highbd_quantize_fp sse4_1 avx2/;
  }

  add_proto qw/void av1_highbd_fwht4x4/, "const int16_t *input, tran_low_t *output, int stride";

  # End av1_high encoder functions

  # txb
  add_proto qw/void av1_get_nz_map_contexts/, "const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TX_SIZE tx_size, const TX_CLASS tx_class, int8_t *const coeff_contexts";
  specialize qw/av1_get_nz_map_contexts sse2/;
  add_proto qw/void av1_txb_init_levels/, "const tran_low_t *const coeff, const int width, const int height, uint8_t *const levels";
  specialize qw/av1_txb_init_levels sse4_1 avx2/;

  add_proto qw/uint64_t av1_wedge_sse_from_residuals/, "const int16_t *r1, const int16_t *d, const uint8_t *m, int N";
  specialize qw/av1_wedge_sse_from_residuals sse2 avx2/;
  add_proto qw/int8_t av1_wedge_sign_from_residuals/, "const int16_t *ds, const uint8_t *m, int N, int64_t limit";
  specialize qw/av1_wedge_sign_from_residuals sse2 avx2/;
  add_proto qw/void av1_wedge_compute_delta_squares/, "int16_t *d, const int16_t *a, const int16_t *b, int N";
  specialize qw/av1_wedge_compute_delta_squares sse2 avx2/;

  # hash
  add_proto qw/uint32_t av1_get_crc32c_value/, "void *crc_calculator, uint8_t *p, size_t length";
  specialize qw/av1_get_crc32c_value sse4_2/;

  add_proto qw/void av1_compute_stats/,  "int wiener_win, const uint8_t *dgd8, const uint8_t *src8, int h_start, int h_end, int v_start, int v_end, int dgd_stride, int src_stride, int64_t *M, int64_t *H";
  specialize qw/av1_compute_stats sse4_1 avx2/;

  add_proto qw/void av1_compute_stats_highbd/,  "int wiener_win, const uint8_t *dgd8, const uint8_t *src8, int h_start, int h_end, int v_start, int v_end, int dgd_stride, int src_stride, int64_t *M, int64_t *H, aom_bit_depth_t bit_depth";
  specialize qw/av1_compute_stats_highbd sse4_1 avx2/;

  add_proto qw/int64_t av1_lowbd_pixel_proj_error/, " const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int xq[2], const sgr_params_type *params";
  specialize qw/av1_lowbd_pixel_proj_error sse4_1 avx2/;

  add_proto qw/int64_t av1_highbd_pixel_proj_error/, " const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int xq[2], const sgr_params_type *params";
  specialize qw/av1_highbd_pixel_proj_error sse4_1 avx2/;

  add_proto qw/void av1_get_horver_correlation_full/, " const int16_t *diff, int stride, int w, int h, float *hcorr, float *vcorr";
  specialize qw/av1_get_horver_correlation_full sse4_1 avx2/;

  add_proto qw/void av1_nn_predict/, " const float *input_nodes, const NN_CONFIG *const nn_config, int reduce_prec, float *const output";
  specialize qw/av1_nn_predict sse3/;
}
# end encoder functions

# CNN functions

add_proto qw/void av1_cnn_activate/, " float **input, int channels, int width, int height, int stride, ACTIVATION layer_activation";
add_proto qw/void av1_cnn_add/, " float **input, int channels, int width, int height, int stride, const float **add";
add_proto qw/void av1_cnn_predict/, " const float **input, int in_width, int in_height, int in_stride, const CNN_CONFIG *cnn_config, const CNN_THREAD_DATA *thread_data, CNN_MULTI_OUT *output_struct";
add_proto qw/void av1_cnn_convolve/, " const float **input, int in_width, int in_height, int in_stride, const CNN_LAYER_CONFIG *layer_config, float **output, int out_stride, int start_idx, int step";
add_proto qw/void av1_cnn_deconvolve/, " const float **input, int in_width, int in_height, int in_stride, const CNN_LAYER_CONFIG *layer_config, float **output, int out_stride";
add_proto qw/void av1_cnn_batchnorm/, "float **image, int channels, int width, int height, int stride, const float *gamma, const float *beta, const float *mean, const float *std";

# Deringing Functions

add_proto qw/int cdef_find_dir/, "const uint16_t *img, int stride, int32_t *var, int coeff_shift";
add_proto qw/void cdef_filter_block/, "uint8_t *dst8, uint16_t *dst16, int dstride, const uint16_t *in, int pri_strength, int sec_strength, int dir, int pri_damping, int sec_damping, int bsize, int coeff_shift";

add_proto qw/void cdef_copy_rect8_8bit_to_16bit/, "uint16_t *dst, int dstride, const uint8_t *src, int sstride, int v, int h";
add_proto qw/void cdef_copy_rect8_16bit_to_16bit/, "uint16_t *dst, int dstride, const uint16_t *src, int sstride, int v, int h";

# VS compiling for 32 bit targets does not support vector types in
# structs as arguments, which makes the v256 type of the intrinsics
# hard to support, so optimizations for this target are disabled.
if ($opts{config} !~ /libs-x86-win32-vs.*/) {
  specialize qw/cdef_find_dir sse2 ssse3 sse4_1 avx2 neon/;
  specialize qw/cdef_filter_block sse2 ssse3 sse4_1 avx2 neon/;
  specialize qw/cdef_copy_rect8_8bit_to_16bit sse2 ssse3 sse4_1 avx2 neon/;
  specialize qw/cdef_copy_rect8_16bit_to_16bit sse2 ssse3 sse4_1 avx2 neon/;
}

# WARPED_MOTION / GLOBAL_MOTION functions

add_proto qw/void av1_warp_affine/, "const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta";
specialize qw/av1_warp_affine sse4_1 avx2 neon/;

if (aom_config("CONFIG_EXT_WARP") eq "yes") {
  add_proto qw/void av1_ext_warp_affine/, "const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, ConvolveParams *conv_params";
}

add_proto qw/void av1_highbd_warp_affine/, "const int32_t *mat, const uint16_t *ref, int width, int height, int stride, uint16_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, int bd, ConvolveParams *conv_params, int16_t alpha, int16_t beta, int16_t gamma, int16_t delta";
specialize qw/av1_highbd_warp_affine sse4_1/;

if (aom_config("CONFIG_EXT_WARP") eq "yes") {
  add_proto qw/void av1_ext_highbd_warp_affine/, "const int32_t *mat, const uint16_t *ref, int width, int height, int stride, uint16_t *pred, int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x, int subsampling_y, int bd, ConvolveParams *conv_params";
}

add_proto qw/int64_t av1_calc_frame_error/, "const uint8_t *const ref, int stride, const uint8_t *const dst, int p_width, int p_height, int p_stride";
specialize qw/av1_calc_frame_error sse2 avx2/;

if (aom_config("CONFIG_AV1_ENCODER") eq "yes") {
  add_proto qw/double av1_compute_cross_correlation/, "unsigned char *im1, int stride1, int x1, int y1, unsigned char *im2, int stride2, int x2, int y2";
  specialize qw/av1_compute_cross_correlation sse4_1 avx2/;
}

# LOOP_RESTORATION functions

add_proto qw/void av1_apply_selfguided_restoration/, "const uint8_t *dat, int width, int height, int stride, int eps, const int *xqd, uint8_t *dst, int dst_stride, int32_t *tmpbuf, int bit_depth, int highbd";
specialize qw/av1_apply_selfguided_restoration sse4_1 avx2 neon/;

add_proto qw/int av1_selfguided_restoration/, "const uint8_t *dgd8, int width, int height,
                                 int dgd_stride, int32_t *flt0, int32_t *flt1, int flt_stride,
                                 int sgr_params_idx, int bit_depth, int highbd";
specialize qw/av1_selfguided_restoration sse4_1 avx2 neon/;

# CONVOLVE_ROUND/COMPOUND_ROUND functions

add_proto qw/void av1_convolve_2d_sr/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int subpel_y_qn, ConvolveParams *conv_params";
add_proto qw/void av1_convolve_x_sr/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn, ConvolveParams *conv_params";
add_proto qw/void av1_convolve_y_sr/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn";
add_proto qw/void av1_dist_wtd_convolve_2d/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int subpel_y_qn, ConvolveParams *conv_params";
add_proto qw/void av1_dist_wtd_convolve_2d_copy/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, ConvolveParams *conv_params";
add_proto qw/void av1_dist_wtd_convolve_x/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x,const int subpel_x_qn, ConvolveParams *conv_params";
add_proto qw/void av1_dist_wtd_convolve_y/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn, ConvolveParams *conv_params";
add_proto qw/void av1_highbd_convolve_2d_copy_sr/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h";
add_proto qw/void av1_highbd_convolve_2d_sr/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int subpel_y_qn, ConvolveParams *conv_params, int bd";
add_proto qw/void av1_highbd_convolve_x_sr/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn, ConvolveParams *conv_params, int bd";
add_proto qw/void av1_highbd_convolve_y_sr/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn, int bd";
add_proto qw/void av1_highbd_dist_wtd_convolve_2d/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int subpel_y_qn, ConvolveParams *conv_params, int bd";
add_proto qw/void av1_highbd_dist_wtd_convolve_x/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn, ConvolveParams *conv_params, int bd";
add_proto qw/void av1_highbd_dist_wtd_convolve_y/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn, ConvolveParams *conv_params, int bd";
add_proto qw/void av1_highbd_dist_wtd_convolve_2d_copy/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, ConvolveParams *conv_params, int bd";

add_proto qw/void av1_convolve_2d_scale/, "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int x_step_qn, const int subpel_y_qn, const int y_step_qn, ConvolveParams *conv_params";
add_proto qw/void av1_highbd_convolve_2d_scale/, "const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int subpel_x_qn, const int x_step_qn, const int subpel_y_qn, const int y_step_qn, ConvolveParams *conv_params, int bd";

specialize qw/av1_convolve_2d_sr sse2 avx2 neon/;
specialize qw/av1_convolve_x_sr sse2 avx2 neon/;
specialize qw/av1_convolve_y_sr sse2 avx2 neon/;
specialize qw/av1_convolve_2d_scale sse4_1/;
specialize qw/av1_dist_wtd_convolve_2d sse2 ssse3 avx2 neon/;
specialize qw/av1_dist_wtd_convolve_2d_copy sse2 avx2/;
specialize qw/av1_dist_wtd_convolve_x sse2 avx2/;
specialize qw/av1_dist_wtd_convolve_y sse2 avx2/;
specialize qw/av1_highbd_convolve_2d_copy_sr sse2 avx2/;
if (aom_config("CONFIG_EXT_PARTITIONS") eq "") {
  specialize qw/av1_dist_wtd_convolve_2d_copy neon/;
  specialize qw/av1_dist_wtd_convolve_x neon/;
  specialize qw/av1_dist_wtd_convolve_y neon/;
} # CONFIG_EXT_PARTITIONS
specialize qw/av1_highbd_convolve_2d_sr ssse3 avx2/;
specialize qw/av1_highbd_convolve_x_sr ssse3 avx2/;
specialize qw/av1_highbd_convolve_y_sr ssse3 avx2/;
specialize qw/av1_highbd_convolve_2d_scale sse4_1/;
specialize qw/av1_highbd_dist_wtd_convolve_2d sse4_1 avx2/;
specialize qw/av1_highbd_dist_wtd_convolve_x sse4_1 avx2/;
specialize qw/av1_highbd_dist_wtd_convolve_y sse4_1 avx2/;
specialize qw/av1_highbd_dist_wtd_convolve_2d_copy sse4_1 avx2/;

# INTRA_EDGE functions
add_proto qw/void av1_filter_intra_edge/, "uint8_t *p, int sz, int strength";
specialize qw/av1_filter_intra_edge sse4_1/;
add_proto qw/void av1_upsample_intra_edge/, "uint8_t *p, int sz";
specialize qw/av1_upsample_intra_edge sse4_1/;

add_proto qw/void av1_filter_intra_edge_high/, "uint16_t *p, int sz, int strength";
specialize qw/av1_filter_intra_edge_high sse4_1/;
add_proto qw/void av1_upsample_intra_edge_high/, "uint16_t *p, int sz, int bd";
specialize qw/av1_upsample_intra_edge_high sse4_1/;

# CFL
add_proto qw/cfl_subtract_average_fn cfl_get_subtract_average_fn/, "TX_SIZE tx_size";
specialize qw/cfl_get_subtract_average_fn sse2 avx2 neon vsx/;

add_proto qw/cfl_subsample_lbd_fn cfl_get_luma_subsampling_420_lbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_420_lbd ssse3 avx2 neon/;

add_proto qw/cfl_subsample_lbd_fn cfl_get_luma_subsampling_422_lbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_422_lbd ssse3 avx2 neon/;

add_proto qw/cfl_subsample_lbd_fn cfl_get_luma_subsampling_444_lbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_444_lbd ssse3 avx2 neon/;

add_proto qw/cfl_subsample_hbd_fn cfl_get_luma_subsampling_420_hbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_420_hbd ssse3 avx2 neon/;

add_proto qw/cfl_subsample_hbd_fn cfl_get_luma_subsampling_422_hbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_422_hbd ssse3 avx2 neon/;

add_proto qw/cfl_subsample_hbd_fn cfl_get_luma_subsampling_444_hbd/, "TX_SIZE tx_size";
specialize qw/cfl_get_luma_subsampling_444_hbd ssse3 avx2 neon/;

add_proto qw/cfl_predict_lbd_fn cfl_get_predict_lbd_fn/, "TX_SIZE tx_size";
specialize qw/cfl_get_predict_lbd_fn ssse3 avx2 neon/;

add_proto qw/cfl_predict_hbd_fn cfl_get_predict_hbd_fn/, "TX_SIZE tx_size";
specialize qw/cfl_get_predict_hbd_fn ssse3 avx2 neon/;

if (aom_config("CONFIG_INTRA_ENTROPY") eq "yes") {
add_proto qw/void av1_nn_fc_forward/, "FC_LAYER_EM *layer, const float *input, float *output";
specialize qw/av1_nn_fc_forward sse4_1/;

add_proto qw/void av1_nn_softmax_em/, "const float *input, float *output, int n";
specialize qw/av1_nn_softmax_em sse4_1/;

if (aom_config("CONFIG_USE_SMALL_MODEL") ne "yes") {
add_proto qw/void av1_get_gradient_hist_lbd/, "const uint8_t *dst, int stride, int rows, int cols, uint64_t *hist";
specialize qw/av1_get_gradient_hist_lbd sse4_1/;
}
}

1;
