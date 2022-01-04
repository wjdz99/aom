/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_COMMON_BLOCKD_DEFINES_H_
#define AOM_AV1_COMMON_BLOCKD_DEFINES_H_

#include "aom_ports/mem_bytes.h"

#ifdef __cplusplus
extern "C" {
#endif

#define USE_B_QUANT_NO_TRELLIS 1

#define MAX_MB_PLANE 3

#define MAX_DIFFWTD_MASK_BITS 1

#define INTERINTRA_WEDGE_SIGN 0

#define DEFAULT_INTER_TX_TYPE DCT_DCT

#define MAX_PALETTE_BLOCK_WIDTH 64

#define MAX_PALETTE_BLOCK_HEIGHT 64

/*!\cond */

// DIFFWTD_MASK_TYPES should not surpass 1 << MAX_DIFFWTD_MASK_BITS
enum {
  DIFFWTD_38 = 0,
  DIFFWTD_38_INV,
  DIFFWTD_MASK_TYPES,
} UENUM1BYTE(DIFFWTD_MASK_TYPE);

enum {
  KEY_FRAME = 0,
  INTER_FRAME = 1,
  INTRA_ONLY_FRAME = 2,  // replaces intra-only
  S_FRAME = 3,
  FRAME_TYPES,
} UENUM1BYTE(FRAME_TYPE);

#define INTER_TX_SIZE_BUF_LEN 16
#define TXK_TYPE_BUF_LEN 64
/*!\endcond */

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_BLOCKD_DEFINES_H_
