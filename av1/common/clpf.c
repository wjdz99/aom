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
#include "av1/common/clpf.h"
#include "aom_dsp/aom_dsp_common.h"

int av1_clpf_maxbits(const AV1_COMMON *cm) {
  return get_msb(
             ALIGN_POWER_OF_TWO(cm->mi_cols * MAX_MIB_SIZE, cm->clpf_size + 4) *
                 ALIGN_POWER_OF_TWO(cm->mi_rows * MAX_MIB_SIZE,
                                    cm->clpf_size + 4) >>
             (cm->clpf_size * 2 + 8)) +
         1;
}

int av1_clpf_sample(int X, int A, int B, int C, int D, int E, int F, int b) {
  int delta = 4 * clamp(A - X, -b, b) + clamp(B - X, -b, b) +
              3 * clamp(C - X, -b, b) + 3 * clamp(D - X, -b, b) +
              clamp(E - X, -b, b) + 4 * clamp(F - X, -b, b);
  return (8 + delta - (delta < 0)) >> 4;
}

static void clpf_block(const uint8_t *src, uint8_t *dst, int stride, int x0,
                       int y0, int sizex, int sizey, int width, int height,
                       unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + sizey; y++) {
    for (x = x0; x < x0 + sizex; x++) {
      int X = src[y * stride + x];
      int A = src[AOMMAX(0, y - 1) * stride + x];
      int B = src[y * stride + AOMMAX(0, x - 2)];
      int C = src[y * stride + AOMMAX(0, x - 1)];
      int D = src[y * stride + AOMMIN(width - 1, x + 1)];
      int E = src[y * stride + AOMMIN(width - 1, x + 2)];
      int F = src[AOMMIN(height - 1, y + 1) * stride + x];
      int delta;
      delta = av1_clpf_sample(X, A, B, C, D, E, F, strength);
      dst[y * stride + x] = X + delta;
    }
  }
}

// Return number of filtered blocks
int av1_clpf_frame(
    const YV12_BUFFER_CONFIG *dst, const YV12_BUFFER_CONFIG *rec,
    const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm, int enable_fb_flag,
    unsigned int strength, unsigned int fb_size_log2, uint8_t *blocks, int comp,
    int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                    const YV12_BUFFER_CONFIG *, const AV1_COMMON *cm, int, int,
                    int, unsigned int, unsigned int, uint8_t *, int)) {
  /* Constrained low-pass filter (CLPF) */
  int c, k, l, m, n;
  int width = comp ? rec->uv_crop_width : rec->y_crop_width;
  int height = comp ? rec->uv_crop_height : rec->y_crop_height;
  int xpos, ypos;
  int stride = comp ? rec->uv_stride : rec->y_stride;
  const int bs = comp && (rec->subsampling_x || rec->subsampling_y) ? 4 : 8;
  const int bslog = get_msb(bs);
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int block_index = 0;
  uint8_t *rec_buffer =
      comp ? (comp == 1 ? rec->u_buffer : rec->v_buffer) : rec->y_buffer;
  uint8_t *dst_buffer =
      comp ? (comp == 1 ? dst->u_buffer : dst->v_buffer) : dst->y_buffer;

  // Iterate over all filter blocks
  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      int h, w;
      int allskip = 1;
      for (m = 0; allskip && m < (1 << fb_size_log2) / MAX_MIB_SIZE; m++) {
        for (n = 0; allskip && n < (1 << fb_size_log2) / MAX_MIB_SIZE; n++) {
          xpos = (l << fb_size_log2) + n * MAX_MIB_SIZE;
          ypos = (k << fb_size_log2) + m * MAX_MIB_SIZE;
          if (xpos < width && ypos < height) {
            allskip &= cm->mi_grid_visible[ypos / MAX_MIB_SIZE * cm->mi_stride +
                                           xpos / MAX_MIB_SIZE]
                           ->mbmi.skip;
          }
        }
      }

      // Calculate the actual filter block size near frame edges
      h = AOMMIN(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      w = AOMMIN(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      h += !h << fb_size_log2;
      w += !w << fb_size_log2;
      if (!allskip &&  // Do not filter the block if all is skip encoded
          (!enable_fb_flag ||
           decision(k, l, rec, org, cm, bs, w >> bslog, h >> bslog, strength,
                    fb_size_log2, blocks + block_index, comp))) {
        // Iterate over all smaller blocks inside the filter block
        for (m = 0; m<(h + bs - 1)>> bslog; m++) {
          for (n = 0; n<(w + bs - 1)>> bslog; n++) {
            xpos = (l << fb_size_log2) + n * bs;
            ypos = (k << fb_size_log2) + m * bs;
            if (!cm->mi_grid_visible[ypos / MAX_MIB_SIZE * cm->mi_stride +
                                     xpos / MAX_MIB_SIZE]
                     ->mbmi.skip) {
              // Not skip block, apply the filter
              clpf_block(rec_buffer, dst_buffer, stride, xpos, ypos, bs, bs,
                         width, height, strength);
            } else {  // Skip block, copy instead
              for (c = 0; c < bs; c++)
                if (bs == 4)
                  *(uint32_t *)(dst_buffer + (ypos + c) * stride + xpos) =
                      *(uint32_t *)(rec_buffer + (ypos + c) * stride + xpos);
                else
                  *(uint64_t *)(dst_buffer + (ypos + c) * stride + xpos) =
                      *(uint64_t *)(rec_buffer + (ypos + c) * stride + xpos);
            }
          }
        }
      } else {  // Entire filter block is skip, copy
        for (m = 0; m < h; m++)
          memcpy(dst_buffer + ((k << fb_size_log2) + m) * stride +
                     (l << fb_size_log2),
                 rec_buffer + ((k << fb_size_log2) + m) * stride +
                     (l << fb_size_log2),
                 w);
      }
      block_index += !allskip;  // Count number of blocks filtered
    }
  }

  return block_index;
}
