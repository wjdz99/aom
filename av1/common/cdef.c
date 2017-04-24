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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "./aom_scale_rtcd.h"
#include "aom/aom_integer.h"
#include "av1/common/cdef.h"
#include "av1/common/cdef_block.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/reconinter.h"

int sb_all_skip(const AV1_COMMON *const cm, int mi_row, int mi_col) {
  int r, c;
  int maxc, maxr;
  int skip = 1;
  maxc = cm->mi_cols - mi_col;
  maxr = cm->mi_rows - mi_row;
#if CONFIG_EXT_PARTITION
  if (maxr > cm->mib_size_log2) maxr = cm->mib_size_log2;
  if (maxc > cm->mib_size_log2) maxc = cm->mib_size_log2;
#else
  if (maxr > MAX_MIB_SIZE) maxr = MAX_MIB_SIZE;
  if (maxc > MAX_MIB_SIZE) maxc = MAX_MIB_SIZE;
#endif

  for (r = 0; r < maxr; r++) {
    for (c = 0; c < maxc; c++) {
      skip = skip &&
             cm->mi_grid_visible[(mi_row + r) * cm->mi_stride + mi_col + c]
                 ->mbmi.skip;
    }
  }
  return skip;
}

static int is_8x8_block_skip(MODE_INFO **grid, int mi_row, int mi_col,
                             int mi_stride) {
  int is_skip = 1;
  for (int r = 0; r < mi_size_high[BLOCK_8X8]; ++r)
    for (int c = 0; c < mi_size_wide[BLOCK_8X8]; ++c)
      is_skip &= grid[(mi_row + r) * mi_stride + (mi_col + c)]->mbmi.skip;

  return is_skip;
}

int sb_compute_cdef_list(const AV1_COMMON *const cm, int mi_row, int mi_col,
                         cdef_list *dlist, int filter_skip) {
  int r, c;
  int maxc, maxr;
  MODE_INFO **grid;
  int count = 0;
  grid = cm->mi_grid_visible;
  maxc = cm->mi_cols - mi_col;
  maxr = cm->mi_rows - mi_row;
#if CONFIG_EXT_PARTITION
  if (maxr > cm->mib_size_log2) maxr = cm->mib_size_log2;
  if (maxc > cm->mib_size_log2) maxc = cm->mib_size_log2;
#else
  if (maxr > MAX_MIB_SIZE) maxr = MAX_MIB_SIZE;
  if (maxc > MAX_MIB_SIZE) maxc = MAX_MIB_SIZE;
#endif

  const int r_step = mi_size_high[BLOCK_8X8];
  const int c_step = mi_size_wide[BLOCK_8X8];
  const int r_shift = (r_step == 2);
  const int c_shift = (c_step == 2);

  assert(r_step == 1 || r_step == 2);
  assert(c_step == 1 || c_step == 2);

  if (filter_skip) {
    for (r = 0; r < maxr; r += r_step) {
      for (c = 0; c < maxc; c += c_step) {
        dlist[count].by = r >> r_shift;
        dlist[count].bx = c >> c_shift;
        dlist[count].skip =
            is_8x8_block_skip(grid, mi_row + r, mi_col + c, cm->mi_stride);
        count++;
      }
    }
  } else {
    for (r = 0; r < maxr; r += r_step) {
      for (c = 0; c < maxc; c += c_step) {
        if (!is_8x8_block_skip(grid, mi_row + r, mi_col + c, cm->mi_stride)) {
          dlist[count].by = r >> r_shift;
          dlist[count].bx = c >> c_shift;
          dlist[count].skip = 0;
          count++;
        }
      }
    }
  }
  return count;
}

void copy_rect8_8bit_to_16bit_c(uint16_t *dst, int dstride, const uint8_t *src,
                                int sstride, int v, int h) {
  int i, j;
  for (i = 0; i < v; i++) {
    for (j = 0; j < h; j++) {
      dst[i * dstride + j] = src[i * sstride + j];
    }
  }
}

void copy_rect8_16bit_to_16bit_c(uint16_t *dst, int dstride,
                                 const uint16_t *src, int sstride, int v,
                                 int h) {
  int i, j;
  for (i = 0; i < v; i++) {
    for (j = 0; j < h; j++) {
      dst[i * dstride + j] = src[i * sstride + j];
    }
  }
}

void copy_sb8_16(UNUSED AV1_COMMON *cm, uint16_t *dst, int dstride,
                 const uint8_t *src, int src_voffset, int src_hoffset,
                 int sstride, int vsize, int hsize) {
#if CONFIG_HIGHBITDEPTH
  if (cm->use_highbitdepth) {
    const uint16_t *base =
        &CONVERT_TO_SHORTPTR(src)[src_voffset * sstride + src_hoffset];
    copy_rect8_16bit_to_16bit(dst, dstride, base, sstride, vsize, hsize);
  } else {
#endif
    const uint8_t *base = &src[src_voffset * sstride + src_hoffset];
    copy_rect8_8bit_to_16bit(dst, dstride, base, sstride, vsize, hsize);
#if CONFIG_HIGHBITDEPTH
  }
#endif
}

static INLINE void fill_rect(uint16_t *dst, int dstride, int v, int h,
                             uint16_t x) {
  int i, j;
  for (i = 0; i < v; i++) {
    for (j = 0; j < h; j++) {
      dst[i * dstride + j] = x;
    }
  }
}

static INLINE void copy_rect(uint16_t *dst, int dstride, const uint16_t *src,
                             int sstride, int v, int h) {
  int i, j;
  for (i = 0; i < v; i++) {
    for (j = 0; j < h; j++) {
      dst[i * dstride + j] = src[i * sstride + j];
    }
  }
}

void av1_cdef_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm,
                    MACROBLOCKD *xd) {
  int sbr, sbc;
  int nhsb, nvsb;
  uint16_t src[CDEF_INBUF_SIZE];
  uint16_t *linebuf[3];
  uint16_t *colbuf[3];
  cdef_list dlist[MAX_MIB_SIZE * MAX_MIB_SIZE];
  unsigned char *row_cdef, *prev_row_cdef, *curr_row_cdef;
  int cdef_count;
  int dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int var[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int stride;
  int mi_wide_l2[3];
  int mi_high_l2[3];
  int xdec[3];
  int ydec[3];
  int pli;
  int cdef_left;
  int coeff_shift = AOMMAX(cm->bit_depth - 8, 0);
  int nplanes = 3;
  int chroma_cdef = xd->plane[1].subsampling_x == xd->plane[1].subsampling_y &&
                    xd->plane[2].subsampling_x == xd->plane[2].subsampling_y;
  nvsb = (cm->mi_rows + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  nhsb = (cm->mi_cols + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  av1_setup_dst_planes(xd->plane, frame, 0, 0);
  row_cdef = aom_malloc(sizeof(*row_cdef) * (nhsb + 2) * 2);
  memset(row_cdef, 1, sizeof(*row_cdef) * (nhsb + 2) * 2);
  prev_row_cdef = row_cdef + 1;
  curr_row_cdef = prev_row_cdef + nhsb + 2;
  for (pli = 0; pli < nplanes; pli++) {
    xdec[pli] = xd->plane[pli].subsampling_x;
    ydec[pli] = xd->plane[pli].subsampling_y;
    mi_wide_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_x;
    mi_high_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_y;
  }
  stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
  for (pli = 0; pli < nplanes; pli++) {
    linebuf[pli] = aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
    colbuf[pli] = aom_malloc(
        sizeof(*colbuf) *
        ((MAX_SB_SIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) * CDEF_HBORDER);
  }
  for (sbr = 0; sbr < nvsb; sbr++) {
    for (pli = 0; pli < nplanes; pli++) {
      const int block_height =
          (MAX_MIB_SIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
      fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER,
                CDEF_VERY_LARGE);
    }
    cdef_left = 1;
    for (sbc = 0; sbc < nhsb; sbc++) {
      int level, sec_strength;
      int uv_level, uv_sec_strength;
      int nhb, nvb;
      int cstart = 0;
      curr_row_cdef[sbc] = 0;
      if (cm->mi_grid_visible[MAX_MIB_SIZE * sbr * cm->mi_stride +
                              MAX_MIB_SIZE * sbc] == NULL ||
          cm->mi_grid_visible[MAX_MIB_SIZE * sbr * cm->mi_stride +
                              MAX_MIB_SIZE * sbc]
                  ->mbmi.cdef_strength == -1) {
        cdef_left = 0;
        continue;
      }
      if (!cdef_left) cstart = -CDEF_HBORDER;
      nhb = AOMMIN(MAX_MIB_SIZE, cm->mi_cols - MAX_MIB_SIZE * sbc);
      nvb = AOMMIN(MAX_MIB_SIZE, cm->mi_rows - MAX_MIB_SIZE * sbr);
      int tile_top, tile_left, tile_bottom, tile_right;
      int mi_idx = MAX_MIB_SIZE * sbr * cm->mi_stride + MAX_MIB_SIZE * sbc;
      BOUNDARY_TYPE boundary_tl =
          cm->mi_grid_visible[MAX_MIB_SIZE * sbr * cm->mi_stride +
                              MAX_MIB_SIZE * sbc]
              ->mbmi.boundary_info;
      tile_top = boundary_tl & TILE_ABOVE_BOUNDARY;
      tile_left = boundary_tl & TILE_LEFT_BOUNDARY;
      /* Right and bottom information appear unreliable, so we use the top
         and left flags for the next superblocks. */
      if (sbr != nvsb - 1 &&
          cm->mi_grid_visible[mi_idx + MAX_MIB_SIZE * cm->mi_stride])
        tile_bottom = cm->mi_grid_visible[mi_idx + MAX_MIB_SIZE * cm->mi_stride]
                          ->mbmi.boundary_info &
                      TILE_ABOVE_BOUNDARY;
      else
        tile_bottom = 1;
      if (sbc != nhsb - 1 && cm->mi_grid_visible[mi_idx + MAX_MIB_SIZE])
        tile_right =
            cm->mi_grid_visible[mi_idx + MAX_MIB_SIZE]->mbmi.boundary_info &
            TILE_LEFT_BOUNDARY;
      else
        tile_right = 1;
      const int mbmi_cdef_strength =
          cm->mi_grid_visible[MAX_MIB_SIZE * sbr * cm->mi_stride +
                              MAX_MIB_SIZE * sbc]
              ->mbmi.cdef_strength;
      level = cm->cdef_strengths[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
      sec_strength =
          cm->cdef_strengths[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
      sec_strength += sec_strength == 3;
      uv_level = cm->cdef_uv_strengths[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
      uv_sec_strength =
          cm->cdef_uv_strengths[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
      uv_sec_strength += uv_sec_strength == 3;
      if ((level == 0 && sec_strength == 0 && uv_level == 0 &&
           uv_sec_strength == 0) ||
          (cdef_count = sb_compute_cdef_list(
               cm, sbr * MAX_MIB_SIZE, sbc * MAX_MIB_SIZE, dlist,
               get_filter_skip(level) || get_filter_skip(uv_level))) == 0) {
        cdef_left = 0;
        continue;
      }

      curr_row_cdef[sbc] = 1;
      for (pli = 0; pli < nplanes; pli++) {
        int coffset;
        int rend, cend;
        int pri_damping = cm->cdef_pri_damping;
        int sec_damping = cm->cdef_sec_damping;
        int hsize = nhb << mi_wide_l2[pli];
        int vsize = nvb << mi_high_l2[pli];

        if (pli) {
          if (chroma_cdef)
            level = uv_level;
          else
            level = 0;
          sec_strength = uv_sec_strength;
        }

        if (sbc == nhsb - 1)
          cend = hsize;
        else
          cend = hsize + CDEF_HBORDER;

        if (sbr == nvsb - 1)
          rend = vsize;
        else
          rend = vsize + CDEF_VBORDER;

        coffset = sbc * MAX_MIB_SIZE << mi_wide_l2[pli];
        if (sbc == nhsb - 1) {
          /* On the last superblock column, fill in the right border with
             CDEF_VERY_LARGE to avoid filtering with the outside. */
          fill_rect(&src[cend + CDEF_HBORDER], CDEF_BSTRIDE,
                    rend + CDEF_VBORDER, hsize + CDEF_HBORDER - cend,
                    CDEF_VERY_LARGE);
        }
        if (sbr == nvsb - 1) {
          /* On the last superblock row, fill in the bottom border with
             CDEF_VERY_LARGE to avoid filtering with the outside. */
          fill_rect(&src[(rend + CDEF_VBORDER) * CDEF_BSTRIDE], CDEF_BSTRIDE,
                    CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
        }
        /* Copy in the pixels we need from the current superblock for
           deringing.*/
        copy_sb8_16(cm,
                    &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                    CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                    (MAX_MIB_SIZE << mi_high_l2[pli]) * sbr, coffset + cstart,
                    xd->plane[pli].dst.stride, rend, cend - cstart);
        if (!prev_row_cdef[sbc]) {
          copy_sb8_16(cm, &src[CDEF_HBORDER], CDEF_BSTRIDE,
                      xd->plane[pli].dst.buf,
                      (MAX_MIB_SIZE << mi_high_l2[pli]) * sbr - CDEF_VBORDER,
                      coffset, xd->plane[pli].dst.stride, CDEF_VBORDER, hsize);
        } else if (sbr > 0) {
          copy_rect(&src[CDEF_HBORDER], CDEF_BSTRIDE, &linebuf[pli][coffset],
                    stride, CDEF_VBORDER, hsize);
        } else {
          fill_rect(&src[CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER, hsize,
                    CDEF_VERY_LARGE);
        }
        if (!prev_row_cdef[sbc - 1]) {
          copy_sb8_16(cm, src, CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                      (MAX_MIB_SIZE << mi_high_l2[pli]) * sbr - CDEF_VBORDER,
                      coffset - CDEF_HBORDER, xd->plane[pli].dst.stride,
                      CDEF_VBORDER, CDEF_HBORDER);
        } else if (sbr > 0 && sbc > 0) {
          copy_rect(src, CDEF_BSTRIDE, &linebuf[pli][coffset - CDEF_HBORDER],
                    stride, CDEF_VBORDER, CDEF_HBORDER);
        } else {
          fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, CDEF_HBORDER,
                    CDEF_VERY_LARGE);
        }
        if (!prev_row_cdef[sbc + 1]) {
          copy_sb8_16(cm, &src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                      CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                      (MAX_MIB_SIZE << mi_high_l2[pli]) * sbr - CDEF_VBORDER,
                      coffset + hsize, xd->plane[pli].dst.stride, CDEF_VBORDER,
                      CDEF_HBORDER);
        } else if (sbr > 0 && sbc < nhsb - 1) {
          copy_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE,
                    &linebuf[pli][coffset + hsize], stride, CDEF_VBORDER,
                    CDEF_HBORDER);
        } else {
          fill_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER,
                    CDEF_HBORDER, CDEF_VERY_LARGE);
        }
        if (cdef_left) {
          /* If we deringed the superblock on the left then we need to copy in
             saved pixels. */
          copy_rect(src, CDEF_BSTRIDE, colbuf[pli], CDEF_HBORDER,
                    rend + CDEF_VBORDER, CDEF_HBORDER);
        }
        /* Saving pixels in case we need to dering the superblock on the
            right. */
        copy_rect(colbuf[pli], CDEF_HBORDER, src + hsize, CDEF_BSTRIDE,
                  rend + CDEF_VBORDER, CDEF_HBORDER);
        copy_sb8_16(
            cm, &linebuf[pli][coffset], stride, xd->plane[pli].dst.buf,
            (MAX_MIB_SIZE << mi_high_l2[pli]) * (sbr + 1) - CDEF_VBORDER,
            coffset, xd->plane[pli].dst.stride, CDEF_VBORDER, hsize);

        if (tile_top) {
          fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, hsize + 2 * CDEF_HBORDER,
                    CDEF_VERY_LARGE);
        }
        if (tile_left) {
          fill_rect(src, CDEF_BSTRIDE, vsize + 2 * CDEF_VBORDER, CDEF_HBORDER,
                    CDEF_VERY_LARGE);
        }
        if (tile_bottom) {
          fill_rect(&src[(vsize + CDEF_VBORDER) * CDEF_BSTRIDE], CDEF_BSTRIDE,
                    CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
        }
        if (tile_right) {
          fill_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE,
                    vsize + 2 * CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
        }
#if CONFIG_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
          cdef_filter_sb(
              NULL,
              &CONVERT_TO_SHORTPTR(
                  xd->plane[pli]
                      .dst.buf)[xd->plane[pli].dst.stride *
                                    (MAX_MIB_SIZE * sbr << mi_high_l2[pli]) +
                                (sbc * MAX_MIB_SIZE << mi_wide_l2[pli])],
              xd->plane[pli].dst.stride,
              &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER], xdec[pli],
              ydec[pli], dir, NULL, var, pli, dlist, cdef_count, level,
              sec_strength, pri_damping, sec_damping, coeff_shift);
        } else {
#endif
          cdef_filter_sb(
              &xd->plane[pli]
                   .dst.buf[xd->plane[pli].dst.stride *
                                (MAX_MIB_SIZE * sbr << mi_high_l2[pli]) +
                            (sbc * MAX_MIB_SIZE << mi_wide_l2[pli])],
              NULL, xd->plane[pli].dst.stride,
              &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER], xdec[pli],
              ydec[pli], dir, NULL, var, pli, dlist, cdef_count, level,
              sec_strength, pri_damping, sec_damping, coeff_shift);

#if CONFIG_HIGHBITDEPTH
        }
#endif
      }
      cdef_left = 1;
    }
    {
      unsigned char *tmp;
      tmp = prev_row_cdef;
      prev_row_cdef = curr_row_cdef;
      curr_row_cdef = tmp;
    }
  }
  aom_free(row_cdef);
  for (pli = 0; pli < nplanes; pli++) {
    aom_free(linebuf[pli]);
    aom_free(colbuf[pli]);
  }
}
