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

int clpf_maxbits(const AV1_COMMON *cm) {
  int bsm1 = (1 << (cm->clpf_size + 4)) - 1;
  return log2i(((cm->mi_cols * MI_BLOCK_SIZE + bsm1) & ~bsm1) *
                   ((cm->mi_rows * MI_BLOCK_SIZE + bsm1) & ~bsm1) >>
               (cm->clpf_size * 2 + 8)) +
         1;
}

int clpf_sample(int X, int A, int B, int C, int D, int E, int F, int b) {
  int delta = 4 * clip(A - X, -b, b) + clip(B - X, -b, b) +
              3 * clip(C - X, -b, b) + 3 * clip(D - X, -b, b) +
              clip(E - X, -b, b) + 4 * clip(F - X, -b, b);
  return (8 + delta - (delta < 0)) >> 4;
}

void clpf_block(const uint8_t *src, uint8_t *dst, int stride, int x0, int y0,
                int sizex, int sizey, int width, int height,
                unsigned int strength) {
  int x, y;
  for (y = y0; y < y0 + sizey; y++) {
    for (x = x0; x < x0 + sizex; x++) {
      int X = src[y * stride + x];
      int A = src[max(0, y - 1) * stride + x];
      int B = src[y * stride + max(0, x - 2)];
      int C = src[y * stride + max(0, x - 1)];
      int D = src[y * stride + min(width - 1, x + 1)];
      int E = src[y * stride + min(width - 1, x + 2)];
      int F = src[min(height - 1, y + 1) * stride + x];
      int delta;
      delta = clpf_sample(X, A, B, C, D, E, F, strength);
      dst[y * stride + x] = X + delta;
    }
  }
}

// Return number of filtered blocks
int av1_clpf_frame(const YV12_BUFFER_CONFIG *dst, const YV12_BUFFER_CONFIG *rec,
                   const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                   int enable_fb_flag, unsigned int strength,
                   unsigned int fb_size_log2, char *blocks,
                   int (*decision)(int, int, const YV12_BUFFER_CONFIG *,
                                   const YV12_BUFFER_CONFIG *,
                                   const AV1_COMMON *cm, int, int, int,
                                   unsigned int, unsigned int, char *)) {
  /* Constrained low-pass filter (CLPF) */
  int c, k, l, m, n;
  int width = rec->y_crop_width;
  int height = rec->y_crop_height;
  int xpos, ypos;
  int stride_y = rec->y_stride;
  int stride_c = rec->uv_stride;
  int subx = rec->subsampling_x;
  int suby = rec->subsampling_y;
  const int bs = 8;  // MI_BLOCK_SIZE;
  int num_fb_hor = (width + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int num_fb_ver = (height + (1 << fb_size_log2) - bs) >> fb_size_log2;
  int block_index = 0;
  const YV12_BUFFER_CONFIG *old_dst = 0;
  YV12_BUFFER_CONFIG tmp;

  // In-place filtering?
  if (rec == dst) {
    tmp = *dst;
    tmp.y_buffer = aom_malloc(dst->y_height * dst->y_stride);
    tmp.u_buffer = aom_malloc(dst->uv_height * dst->uv_stride);
    tmp.v_buffer = aom_malloc(dst->uv_height * dst->uv_stride);
    old_dst = dst;
    dst = &tmp;
  }

  // Iterate over all filter blocks
  for (k = 0; k < num_fb_ver; k++) {
    for (l = 0; l < num_fb_hor; l++) {
      int h, w;
      int allskip = 1;
      for (m = 0; allskip && m < (1 << fb_size_log2) / bs; m++) {
        for (n = 0; allskip && n < (1 << fb_size_log2) / bs; n++) {
          xpos = (l << fb_size_log2) + n * bs;
          ypos = (k << fb_size_log2) + m * bs;
          if (xpos < width && ypos < height) {
            allskip &=
                cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
                    ->mbmi.skip;
          }
        }
      }

      // Calculate the actual filter block size near frame edges
      h = min(height, (k + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      w = min(width, (l + 1) << fb_size_log2) & ((1 << fb_size_log2) - 1);
      h += !h << fb_size_log2;
      w += !w << fb_size_log2;
      if (!allskip &&  // Do not filter the block if all is skip encoded
          (!enable_fb_flag ||
           decision(k, l, rec, org, cm, bs, w / bs, h / bs, strength,
                    fb_size_log2, blocks + block_index))) {
        // Iterate over all smaller blocks inside the filter block
        for (m = 0; m < (h + bs - 1) / bs; m++) {
          for (n = 0; n < (w + bs - 1) / bs; n++) {
            xpos = (l << fb_size_log2) + n * bs;
            ypos = (k << fb_size_log2) + m * bs;
            if (!cm->mi_grid_visible[ypos / bs * cm->mi_stride + xpos / bs]
                     ->mbmi.skip) {
              // Not skip block, apply the filter
              clpf_block(rec->y_buffer, dst->y_buffer, stride_y, xpos, ypos, bs,
                         bs, width, height, strength);
              clpf_block(rec->u_buffer, dst->u_buffer, stride_c, xpos >> subx,
                         ypos >> suby, bs >> subx, bs >> suby, width >> subx,
                         height >> suby, strength);
              clpf_block(rec->v_buffer, dst->v_buffer, stride_c, xpos >> subx,
                         ypos >> suby, bs >> subx, bs >> suby, width >> subx,
                         height >> suby, strength);
            } else {  // Skip block, copy instead
              for (c = 0; c < bs; c++)
                *(uint64_t *)(dst->y_buffer + (ypos + c) * stride_y + xpos) =
                    *(uint64_t *)(rec->y_buffer + (ypos + c) * stride_y + xpos);
              if (subx) {
                for (c = 0; c < (bs >> suby); c++) {
                  *(uint32_t *)(dst->u_buffer +
                                ((ypos >> suby) + c) * stride_c +
                                (xpos >> subx)) =
                      *(uint32_t *)(rec->u_buffer +
                                    ((ypos >> suby) + c) * stride_c +
                                    (xpos >> subx));
                  *(uint32_t *)(dst->v_buffer +
                                ((ypos >> suby) + c) * stride_c +
                                (xpos >> subx)) =
                      *(uint32_t *)(rec->v_buffer +
                                    ((ypos >> suby) + c) * stride_c +
                                    (xpos >> subx));
                }
              } else {
                for (c = 0; c < (bs >> suby); c++) {
                  *(uint64_t *)(dst->u_buffer +
                                ((ypos >> suby) + c) * stride_c +
                                (xpos >> subx)) =
                      *(uint64_t *)(rec->u_buffer +
                                    ((ypos >> suby) + c) * stride_c +
                                    (xpos >> subx));
                  *(uint64_t *)(dst->v_buffer +
                                ((ypos >> suby) + c) * stride_c +
                                (xpos >> subx)) =
                      *(uint64_t *)(rec->v_buffer +
                                    ((ypos >> suby) + c) * stride_c +
                                    (xpos >> subx));
                }
              }
            }
          }
        }
      } else {  // Entire filter block is skip, copy
        for (m = 0; m < h; m++)
          memcpy(dst->y_buffer + ((k << fb_size_log2) + m) * stride_y +
                     (l << fb_size_log2),
                 rec->y_buffer + ((k << fb_size_log2) + m) * stride_y +
                     (l << fb_size_log2),
                 w);
        for (m = 0; m < (h >> suby); m++) {
          memcpy(
              dst->u_buffer + (((k << fb_size_log2) >> suby) + m) * stride_c +
                  ((l << fb_size_log2) >> subx),
              rec->u_buffer + (((k << fb_size_log2) >> suby) + m) * stride_c +
                  ((l << fb_size_log2) >> subx),
              w >> subx);
          memcpy(
              dst->v_buffer + (((k << fb_size_log2) >> suby) + m) * stride_c +
                  ((l << fb_size_log2) >> subx),
              rec->v_buffer + (((k << fb_size_log2) >> suby) + m) * stride_c +
                  ((l << fb_size_log2) >> subx),
              w >> subx);
        }
      }
      block_index += !allskip;  // Count number of blocks filtered
    }
  }

  // Needed for in-place filtering.  The frame copy can be avoided if the caller
  // provides a fresh destination frame and then update the original's pointers.
  if (old_dst) {
    memcpy(old_dst->y_buffer, tmp.y_buffer, dst->y_height * dst->y_stride);
    memcpy(old_dst->u_buffer, tmp.u_buffer, dst->uv_height * dst->uv_stride);
    memcpy(old_dst->v_buffer, tmp.v_buffer, dst->uv_height * dst->uv_stride);
    aom_free(tmp.y_buffer);
    aom_free(tmp.u_buffer);
    aom_free(tmp.v_buffer);
  }
  return block_index;
}
