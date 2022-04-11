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

#include <arm_neon.h>
#include <assert.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom/aom_integer.h"

// -----------------------------------------------------------------------------
// PAETH

static INLINE void highbd_paeth_4or8_x_h_neon(uint16_t *dest, ptrdiff_t stride,
                                              const uint16_t *const top_row,
                                              const uint16_t *const left_column,
                                              int width, int height) {
  const uint16x8_t top_left = vdupq_n_u16(top_row[-1]);
  const uint16x8_t top_left_x2 = vdupq_n_u16(top_row[-1] + top_row[-1]);
  uint16x8_t top;
  if (width == 4) {
    top = vcombine_u16(vld1_u16(top_row), vdup_n_u16(0));
  } else {  // width == 8
    top = vld1q_u16(top_row);
  }

  for (int y = 0; y < height; ++y) {
    const uint16x8_t left = vdupq_n_u16(left_column[y]);

    const uint16x8_t left_dist = vabdq_u16(top, top_left);
    const uint16x8_t top_dist = vabdq_u16(left, top_left);
    const uint16x8_t top_left_dist =
        vabdq_u16(vaddq_u16(top, left), top_left_x2);

    const uint16x8_t left_le_top = vcleq_u16(left_dist, top_dist);
    const uint16x8_t left_le_top_left = vcleq_u16(left_dist, top_left_dist);
    const uint16x8_t top_le_top_left = vcleq_u16(top_dist, top_left_dist);

    // if (left_dist <= top_dist && left_dist <= top_left_dist)
    const uint16x8_t left_mask = vandq_u16(left_le_top, left_le_top_left);
    //   dest[x] = left_column[y];
    // Fill all the unused spaces with 'top'. They will be overwritten when
    // the positions for top_left are known.
    uint16x8_t result = vbslq_u16(left_mask, left, top);
    // else if (top_dist <= top_left_dist)
    //   dest[x] = top_row[x];
    // Add these values to the mask. They were already set.
    const uint16x8_t left_or_top_mask = vorrq_u16(left_mask, top_le_top_left);
    // else
    //   dest[x] = top_left;
    result = vbslq_u16(left_or_top_mask, result, top_left);

    if (width == 4) {
      vst1_u16(dest, vget_low_u16(result));
    } else {  // width == 8
      vst1q_u16(dest, result);
    }
    dest += stride;
  }
}

#define HIGHBD_PAETH_NXM(W, H)                                  \
  void aom_highbd_paeth_predictor_##W##x##H##_neon(             \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,   \
      const uint16_t *left, int bd) {                           \
    (void)bd;                                                   \
    highbd_paeth_4or8_x_h_neon(dst, stride, above, left, W, H); \
  }

HIGHBD_PAETH_NXM(4, 4)
HIGHBD_PAETH_NXM(4, 8)
HIGHBD_PAETH_NXM(8, 4)
HIGHBD_PAETH_NXM(8, 8)
HIGHBD_PAETH_NXM(8, 16)

#if !CONFIG_REALTIME_ONLY
HIGHBD_PAETH_NXM(4, 16)
HIGHBD_PAETH_NXM(8, 32)
#endif

static INLINE uint16x8_t x_le_top_left(const uint16x8_t x_dist,
                                       const uint16x8_t top_left_dist_low) {
  return vcleq_u16(x_dist, top_left_dist_low);
}

// Select the closest values and collect them.
static INLINE uint16x8_t select_paeth(const uint16x8_t top,
                                      const uint16x8_t left,
                                      const uint16x8_t top_left,
                                      const uint16x8_t left_le_top,
                                      const uint16x8_t left_le_top_left,
                                      const uint16x8_t top_le_top_left) {
  // if (left_dist <= top_dist && left_dist <= top_left_dist)
  const uint16x8_t left_mask = vandq_u16(left_le_top, left_le_top_left);
  //   dest[x] = left_column[y];
  // Fill all the unused spaces with 'top'. They will be overwritten when
  // the positions for top_left are known.
  uint16x8_t result = vbslq_u16(left_mask, left, top);
  // else if (top_dist <= top_left_dist)
  //   dest[x] = top_row[x];
  // Add these values to the mask. They were already set.
  const uint16x8_t left_or_top_mask = vorrq_u16(left_mask, top_le_top_left);
  // else
  //   dest[x] = top_left;
  return vbslq_u16(left_or_top_mask, result, top_left);
}

// Generate numbered and high/low versions of top_left_dist.
#define TOP_LEFT_DIST(num)                     \
  const uint16x8_t top_left_##num##_dist_low = \
      vabdq_u16(vaddq_u16((top[num]), (left)), top_left_x2);

// Generate numbered versions of XLeTopLeft with x = left.
#define LEFT_LE_TOP_LEFT(num)               \
  const uint16x8_t left_le_top_left_##num = \
      x_le_top_left(left_##num##_dist, top_left_##num##_dist_low)

// Generate numbered versions of XLeTopLeft with x = top.
#define TOP_LE_TOP_LEFT(num)               \
  const uint16x8_t top_le_top_left_##num = \
      x_le_top_left(top_dist, top_left_##num##_dist_low)

#define PAETH_PREDICTOR(num)                                          \
  const uint16x8_t left_##num##_dist = vabdq_u16(top[num], top_left); \
  TOP_LEFT_DIST(num);                                                 \
  const uint16x8_t left_##num##_le_top =                              \
      vcleq_u16(left_##num##_dist, top_dist);                         \
  LEFT_LE_TOP_LEFT(num);                                              \
  TOP_LE_TOP_LEFT(num);                                               \
  const uint16x8_t result_##num =                                     \
      select_paeth(top[num], left, top_left, left_##num##_le_top,     \
                   left_le_top_left_##num, top_le_top_left_##num);    \
  vst1q_u16(dest + (num * 8), result_##num)

static INLINE void highbd_paeth16_plus_x_h_neon(
    uint16_t *dest, ptrdiff_t stride, const uint16_t *const top_row,
    const uint16_t *const left_column, int width, int height) {
  const uint16x8_t top_left = vdupq_n_u16(top_row[-1]);
  const uint16x8_t top_left_x2 = vdupq_n_u16(top_row[-1] + top_row[-1]);
  uint16x8_t top[8];
  top[0] = vld1q_u16(top_row);
  top[1] = vld1q_u16(top_row + 8);
  if (width > 16) {
    top[2] = vld1q_u16(top_row + 16);
    top[3] = vld1q_u16(top_row + 24);
    if (width == 64) {
      top[4] = vld1q_u16(top_row + 32);
      top[5] = vld1q_u16(top_row + 40);
      top[6] = vld1q_u16(top_row + 48);
      top[7] = vld1q_u16(top_row + 56);
    }
  }

  for (int y = 0; y < height; ++y) {
    const uint16x8_t left = vdupq_n_u16(left_column[y]);
    const uint16x8_t top_dist = vabdq_u16(left, top_left);
    PAETH_PREDICTOR(0);
    PAETH_PREDICTOR(1);
    if (width > 16) {
      PAETH_PREDICTOR(2);
      PAETH_PREDICTOR(3);
      if (width == 64) {
        PAETH_PREDICTOR(4);
        PAETH_PREDICTOR(5);
        PAETH_PREDICTOR(6);
        PAETH_PREDICTOR(7);
      }
    }
    dest += stride;
  }
}

#define HIGHBD_PAETH_NXM_WIDE(W, H)                               \
  void aom_highbd_paeth_predictor_##W##x##H##_neon(               \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,     \
      const uint16_t *left, int bd) {                             \
    (void)bd;                                                     \
    highbd_paeth16_plus_x_h_neon(dst, stride, above, left, W, H); \
  }

HIGHBD_PAETH_NXM_WIDE(16, 8)
HIGHBD_PAETH_NXM_WIDE(16, 16)
HIGHBD_PAETH_NXM_WIDE(16, 32)
HIGHBD_PAETH_NXM_WIDE(32, 16)
HIGHBD_PAETH_NXM_WIDE(32, 32)
HIGHBD_PAETH_NXM_WIDE(32, 64)
HIGHBD_PAETH_NXM_WIDE(64, 32)
HIGHBD_PAETH_NXM_WIDE(64, 64)

#if !CONFIG_REALTIME_ONLY
HIGHBD_PAETH_NXM_WIDE(16, 4)
HIGHBD_PAETH_NXM_WIDE(16, 64)
HIGHBD_PAETH_NXM_WIDE(32, 8)
HIGHBD_PAETH_NXM_WIDE(64, 16)
#endif
