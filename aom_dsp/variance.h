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

#ifndef AOM_DSP_VARIANCE_H_
#define AOM_DSP_VARIANCE_H_

#include "./aom_config.h"

#include "aom/aom_integer.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FILTER_BITS 7
#define FILTER_WEIGHT 128

typedef unsigned int (*AomSadFnT)(const uint8_t *a, int a_stride,
                                  const uint8_t *b, int b_stride);

typedef unsigned int (*AomSadAvgFnT)(const uint8_t *a, int a_stride,
                                     const uint8_t *b, int b_stride,
                                     const uint8_t *second_pred);

typedef void (*AomCopy32xnFnT)(const uint8_t *a, int a_stride, uint8_t *b,
                               int b_stride, int n);

typedef void (*AomSadMultiFnT)(const uint8_t *a, int a_stride, const uint8_t *b,
                               int b_stride, unsigned int *sad_array);

typedef void (*AomSadMultiDFnT)(const uint8_t *a, int a_stride,
                                const uint8_t *const b_array[], int b_stride,
                                unsigned int *sad_array);

typedef unsigned int (*AomVarianceFnT)(const uint8_t *a, int a_stride,
                                       const uint8_t *b, int b_stride,
                                       unsigned int *sse);

typedef unsigned int (*AomSubpixvarianceFnT)(const uint8_t *a, int a_stride,
                                             int xoffset, int yoffset,
                                             const uint8_t *b, int b_stride,
                                             unsigned int *sse);

typedef unsigned int (*AomSubpAvgVarianceFnT)(const uint8_t *a, int a_stride,
                                              int xoffset, int yoffset,
                                              const uint8_t *b, int b_stride,
                                              unsigned int *sse,
                                              const uint8_t *second_pred);

#if CONFIG_AV1 && CONFIG_EXT_INTER
typedef unsigned int (*AomMaskedSadFnT)(const uint8_t *src, int src_stride,
                                        const uint8_t *ref, int ref_stride,
                                        const uint8_t *msk_ptr, int msk_stride);
typedef unsigned int (*AomMaskedVarianceFnT)(const uint8_t *src, int src_stride,
                                             const uint8_t *ref, int ref_stride,
                                             const uint8_t *msk, int msk_stride,
                                             unsigned int *sse);
typedef unsigned int (*AomMaskedSubpixvarianceFnT)(
    const uint8_t *src, int src_stride, int xoffset, int yoffset,
    const uint8_t *ref, int ref_stride, const uint8_t *msk, int msk_stride,
    unsigned int *sse);
#endif  // CONFIG_AV1 && CONFIG_EXT_INTER

#if CONFIG_AV1 && CONFIG_MOTION_VAR
typedef unsigned int (*AomObmcSadFnT)(const uint8_t *pred, int pred_stride,
                                      const int32_t *wsrc, const int32_t *msk);
typedef unsigned int (*AomObmcVarianceFnT)(const uint8_t *pred, int pred_stride,
                                           const int32_t *wsrc,
                                           const int32_t *msk,
                                           unsigned int *sse);
typedef unsigned int (*AomObmcSubpixvarianceFnT)(
    const uint8_t *pred, int pred_stride, int xoffset, int yoffset,
    const int32_t *wsrc, const int32_t *msk, unsigned int *sse);
#endif  // CONFIG_AV1 && CONFIG_MOTION_VAR

#if CONFIG_AV1
typedef struct AomVarianceVtable {
  AomSadFnT sdf;
  AomSadAvgFnT sdaf;
  AomVarianceFnT vf;
  AomSubpixvarianceFnT svf;
  AomSubpAvgVarianceFnT svaf;
  AomSadMultiFnT sdx3f;
  AomSadMultiFnT sdx8f;
  AomSadMultiDFnT sdx4df;
#if CONFIG_EXT_INTER
  AomMaskedSadFnT msdf;
  AomMaskedVarianceFnT mvf;
  AomMaskedSubpixvarianceFnT msvf;
#endif  // CONFIG_EXT_INTER
#if CONFIG_MOTION_VAR
  AomObmcSadFnT osdf;
  AomObmcVarianceFnT ovf;
  AomObmcSubpixvarianceFnT osvf;
#endif  // CONFIG_MOTION_VAR
} AomVarianceFnPtrT;
#endif  // CONFIG_AV1

void aom_highbd_var_filter_block2d_bil_first_pass(
    const uint8_t *src_ptr8, uint16_t *output_ptr,
    unsigned int src_pixels_per_line, int pixel_step,
    unsigned int output_height, unsigned int output_width,
    const uint8_t *filter);

void aom_highbd_var_filter_block2d_bil_second_pass(
    const uint16_t *src_ptr, uint16_t *output_ptr,
    unsigned int src_pixels_per_line, unsigned int pixel_step,
    unsigned int output_height, unsigned int output_width,
    const uint8_t *filter);

uint32_t aom_sse_odd_size(const uint8_t *a, int a_stride, const uint8_t *b,
                          int b_stride, int w, int h);

#if CONFIG_HIGHBITDEPTH
uint64_t aom_highbd_sse_odd_size(const uint8_t *a, int a_stride,
                                 const uint8_t *b, int b_stride, int w, int h);
#endif  // CONFIG_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_VARIANCE_H_
