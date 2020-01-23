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

#ifndef AOM_AV1_ENCODER_TEMPORAL_FILTER_H_
#define AOM_AV1_ENCODER_TEMPORAL_FILTER_H_

#ifdef __cplusplus
extern "C" {
#endif

#define ARNR_FILT_QINDEX 128

// Block size used in temporal filtering
#define TF_BLOCK BLOCK_32X32
#define BH 32
#define BW 32

#define NUM_KEY_FRAME_DENOISING 7

// Window size for temporal filtering on YUV planes.
// This is particually used for function `av1_apply_temporal_filter_yuv()`.
#define YUV_FILTER_WINDOW_LENGTH 3

// Window size for temporal filtering on Y planes.
// This is particually used for function `av1_apply_temporal_filter_yonly()`.
#define YONLY_FILTER_WINDOW_LENGTH 3

#define ENABLE_PLANEWISE_STRATEGY 1
// Window size for plane-wise temporal filtering.
// This is particually used for function `av1_apply_temporal_filter_planewise()`
#define PLANEWISE_FILTER_WINDOW_LENGTH 5
// A scale factor used in plane-wise temporal filtering to raise the filter
// weight from `double` with range [0, 1] to `int` with range [0, 1000].
#define PLANEWISE_FILTER_WEIGHT_SCALE 1000

#define NOISE_ESTIMATION_EDGE_THRESHOLD 50
// Estimates noise level from a given frame (ONLY using Y-plane).
// This is an adaptation of the mehtod in the following paper:
// Shen-Chuan Tai, Shih-Ming Yang, "A fast method for image noise
// estimation using Laplacian operator and adaptive edge detection",
// Proc. 3rd International Symposium on Communications, Control and
// Signal Processing, 2008, St Julians, Malta.
// Inputs:
//   frame: Pointer to the frame to estimate noise level from.
//   bit_depth: Actual bit-depth of the frame.
// Returns:
//   The estimated noise, or -1.0 if there are too few smooth pixels.
double av1_estimate_noise(const YV12_BUFFER_CONFIG *frame, const int bit_depth);

int av1_temporal_filter(AV1_COMP *cpi, int distance,
                        int *show_existing_alt_ref);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_TEMPORAL_FILTER_H_
