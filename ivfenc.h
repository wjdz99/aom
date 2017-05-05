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
#ifndef IVFENC_H_
#define IVFENC_H_

#include "./tools_common.h"

struct AomCodecEncCfg;
struct AomCodecCxPkt;

#ifdef __cplusplus
extern "C" {
#endif

void ivf_write_file_header(FILE *outfile, const struct AomCodecEncCfg *cfg,
                           uint32_t fourcc, int frame_cnt);

void ivf_write_frame_header(FILE *outfile, int64_t pts, size_t frame_size);

void ivf_write_frame_size(FILE *outfile, size_t frame_size);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // IVFENC_H_
