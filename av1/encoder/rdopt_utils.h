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

#ifndef AOM_AV1_ENCODER_RDOPT_UTILS_H_
#define AOM_AV1_ENCODER_RDOPT_UTILS_H_

#include "aom/aom_integer.h"
#include "av1/encoder/block.h"
#include "av1/common/blockd.h"
#include "av1/encoder/encoder.h"
#include "config/aom_dsp_rtcd.h"

#ifdef __cplusplus
extern "C" {
#endif

static AOM_INLINE double get_sad_norm(const int16_t *diff, int stride, int w, int h) {
    double sum = 0.0;
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            sum += abs(diff[j * stride + i]);
        }
    }
    assert(w > 0 && h > 0);
    return sum / (w * h);
}

static AOM_INLINE double get_sse_norm(const int16_t *diff, int stride, int w, int h) {
    double sum = 0.0;
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            const int err = diff[j * stride + i];
            sum += err * err;
        }
    }
    assert(w > 0 && h > 0);
    return sum / (w * h);
}

static AOM_INLINE void get_2x2_normalized_sses_and_sads(
    const AV1_COMP *const cpi, BLOCK_SIZE tx_bsize, const uint8_t *const src,
    int src_stride, const uint8_t *const dst, int dst_stride,
    const int16_t *const src_diff, int diff_stride, double *const sse_norm_arr,
    double *const sad_norm_arr) {
    const BLOCK_SIZE tx_bsize_half =
        get_partition_subsize(tx_bsize, PARTITION_SPLIT);
    if (tx_bsize_half == BLOCK_INVALID) {  // manually calculate stats
        const int half_width = block_size_wide[tx_bsize] / 2;
        const int half_height = block_size_high[tx_bsize] / 2;
        for (int row = 0; row < 2; ++row) {
            for (int col = 0; col < 2; ++col) {
                const int16_t *const this_src_diff =
                    src_diff + row * half_height * diff_stride + col * half_width;
                if (sse_norm_arr) {
                    sse_norm_arr[row * 2 + col] =
                        get_sse_norm(this_src_diff, diff_stride, half_width, half_height);
                }
                if (sad_norm_arr) {
                    sad_norm_arr[row * 2 + col] =
                        get_sad_norm(this_src_diff, diff_stride, half_width, half_height);
                }
            }
        }
    }
    else {  // use function pointers to calculate stats
        const int half_width = block_size_wide[tx_bsize_half];
        const int half_height = block_size_high[tx_bsize_half];
        const int num_samples_half = half_width * half_height;
        for (int row = 0; row < 2; ++row) {
            for (int col = 0; col < 2; ++col) {
                const uint8_t *const this_src =
                    src + row * half_height * src_stride + col * half_width;
                const uint8_t *const this_dst =
                    dst + row * half_height * dst_stride + col * half_width;

                if (sse_norm_arr) {
                    unsigned int this_sse;
                    cpi->fn_ptr[tx_bsize_half].vf(this_src, src_stride, this_dst,
                        dst_stride, &this_sse);
                    sse_norm_arr[row * 2 + col] = (double)this_sse / num_samples_half;
                }

                if (sad_norm_arr) {
                    const unsigned int this_sad = cpi->fn_ptr[tx_bsize_half].sdf(
                        this_src, src_stride, this_dst, dst_stride);
                    sad_norm_arr[row * 2 + col] = (double)this_sad / num_samples_half;
                }
            }
        }
    }
}

static int64_t get_sse(const AV1_COMP *cpi, const MACROBLOCK *x) {
    const AV1_COMMON *cm = &cpi->common;
    const int num_planes = av1_num_planes(cm);
    const MACROBLOCKD *xd = &x->e_mbd;
    const MB_MODE_INFO *mbmi = xd->mi[0];
    int64_t total_sse = 0;
    for (int plane = 0; plane < num_planes; ++plane) {
        const struct macroblock_plane *const p = &x->plane[plane];
        const struct macroblockd_plane *const pd = &xd->plane[plane];
        const BLOCK_SIZE bs = get_plane_block_size(mbmi->sb_type, pd->subsampling_x,
            pd->subsampling_y);
        unsigned int sse;

        if (x->skip_chroma_rd && plane) continue;

        cpi->fn_ptr[bs].vf(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
            &sse);
        total_sse += sse;
    }
    total_sse <<= 4;
    return total_sse;
}

static AOM_INLINE int64_t av1_block_error_c(const tran_low_t *coeff, const tran_low_t *dqcoeff,
    intptr_t block_size, int64_t *ssz) {
    int i;
    int64_t error = 0, sqcoeff = 0;

    for (i = 0; i < block_size; i++) {
        const int diff = coeff[i] - dqcoeff[i];
        error += diff * diff;
        sqcoeff += coeff[i] * coeff[i];
    }

    *ssz = sqcoeff;
    return error;
}

#if CONFIG_AV1_HIGHBITDEPTH
static int64_t av1_highbd_block_error_c(const tran_low_t *coeff,
    const tran_low_t *dqcoeff, intptr_t block_size,
    int64_t *ssz, int bd) {
    int i;
    int64_t error = 0, sqcoeff = 0;
    int shift = 2 * (bd - 8);
    int rounding = shift > 0 ? 1 << (shift - 1) : 0;

    for (i = 0; i < block_size; i++) {
        const int64_t diff = coeff[i] - dqcoeff[i];
        error += diff * diff;
        sqcoeff += (int64_t)coeff[i] * (int64_t)coeff[i];
    }
    assert(error >= 0 && sqcoeff >= 0);
    error = (error + rounding) >> shift;
    sqcoeff = (sqcoeff + rounding) >> shift;

    *ssz = sqcoeff;
    return error;
}
#endif

static AOM_INLINE double get_highbd_diff_mean(const uint8_t *src8, int src_stride,
    const uint8_t *dst8, int dst_stride, int w,
    int h) {
    const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    const uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
    double sum = 0.0;
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            const int diff = src[j * src_stride + i] - dst[j * dst_stride + i];
            sum += diff;
        }
    }
    assert(w > 0 && h > 0);
    return sum / (w * h);
}

static AOM_INLINE double get_diff_mean(const uint8_t *src, int src_stride,
    const uint8_t *dst, int dst_stride, int w, int h) {
    double sum = 0.0;
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            const int diff = src[j * src_stride + i] - dst[j * dst_stride + i];
            sum += diff;
        }
    }
    assert(w > 0 && h > 0);
    return sum / (w * h);
}

static int64_t calculate_sse(MACROBLOCKD *const xd,
    const struct macroblock_plane *p,
    struct macroblockd_plane *pd, const int bw,
    const int bh) {
    int64_t sse = 0;
    const int shift = xd->bd - 8;
#if CONFIG_AV1_HIGHBITDEPTH
    if (is_cur_buf_hbd(xd)) {
        sse = aom_highbd_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride,
            bw, bh);
    }
    else {
        sse =
            aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
    }
#else
    sse = aom_sse(p->src.buf, p->src.stride, pd->dst.buf, pd->dst.stride, bw, bh);
#endif
    sse = ROUND_POWER_OF_TWO(sse, shift * 2);
    return sse;
}


#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_RDOPT_UTILS_H_
