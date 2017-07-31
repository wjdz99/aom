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
#if !CONFIG_CDEF_SINGLEPASS
#include "av1/common/od_dering.h"
#else
#include "av1/common/cdef_block.h"
#endif
#include "av1/common/onyxc_int.h"
#include "av1/common/reconinter.h"

int sb_all_skip(const AV1_COMMON *const cm, int mi_row, int mi_col) {
  int r, c;
  int maxc, maxr;
  int skip = 1;
  maxc = cm->mi_cols - mi_col;
  maxr = cm->mi_rows - mi_row;

  maxr = AOMMIN(maxr, MI_SIZE_64X64);
  maxc = AOMMIN(maxc, MI_SIZE_64X64);

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

#if !CONFIG_CDEF_SINGLEPASS
int sb_compute_dering_list(const AV1_COMMON *const cm, int mi_row, int mi_col,
                           dering_list *dlist, int filter_skip) {
#else
int sb_compute_cdef_list(const AV1_COMMON *const cm, int mi_row, int mi_col,
                         cdef_list *dlist, int filter_skip) {
#endif
  int r, c;
  int maxc, maxr;
  MODE_INFO **grid;
  int count = 0;
  grid = cm->mi_grid_visible;
  maxc = cm->mi_cols - mi_col;
  maxr = cm->mi_rows - mi_row;

  maxr = AOMMIN(maxr, MI_SIZE_64X64);
  maxc = AOMMIN(maxc, MI_SIZE_64X64);

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

static void copy_sb8_16(UNUSED AV1_COMMON *cm, uint16_t *dst, int dstride,
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
#if !CONFIG_CDEF_SINGLEPASS
  int sbr, sbc;
  int nhsb, nvsb;
  uint16_t src[OD_DERING_INBUF_SIZE];
#else
  int fbr, fbc;
  int nhfb, nvfb;
  uint16_t src[CDEF_INBUF_SIZE];
#endif
  uint16_t *linebuf[3];
  uint16_t *colbuf[3];
#if !CONFIG_CDEF_SINGLEPASS
  dering_list dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
  unsigned char *row_dering, *prev_row_dering, *curr_row_dering;
  int dering_count;
  int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS] = { { 0 } };
  int var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS] = { { 0 } };
#else
  cdef_list dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
  unsigned char *row_cdef, *prev_row_cdef, *curr_row_cdef;
  int cdef_count;
  int dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int var[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
#endif
  int stride;
  int mi_wide_l2[3];
  int mi_high_l2[3];
  int xdec[3];
  int ydec[3];
  int pli;
#if !CONFIG_CDEF_SINGLEPASS
  int dering_left;
#else
  int cdef_left;
#endif
  int coeff_shift = AOMMAX(cm->bit_depth - 8, 0);
#if !CONFIG_CDEF_SINGLEPASS
  int nplanes = 3;
  int chroma_dering =
      xd->plane[1].subsampling_x == xd->plane[1].subsampling_y &&
      xd->plane[2].subsampling_x == xd->plane[2].subsampling_y;
  nvsb = (cm->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  nhsb = (cm->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
#else
  int nplanes = MAX_MB_PLANE;
  int chroma_cdef = xd->plane[1].subsampling_x == xd->plane[1].subsampling_y &&
                    xd->plane[2].subsampling_x == xd->plane[2].subsampling_y;
  nvfb = (cm->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  nhfb = (cm->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
#endif
  av1_setup_dst_planes(xd->plane, cm->sb_size, frame, 0, 0);
#if !CONFIG_CDEF_SINGLEPASS
  row_dering = aom_malloc(sizeof(*row_dering) * (nhsb + 2) * 2);
  memset(row_dering, 1, sizeof(*row_dering) * (nhsb + 2) * 2);
  prev_row_dering = row_dering + 1;
  curr_row_dering = prev_row_dering + nhsb + 2;
#else
  row_cdef = aom_malloc(sizeof(*row_cdef) * (nhfb + 2) * 2);
  memset(row_cdef, 1, sizeof(*row_cdef) * (nhfb + 2) * 2);
  prev_row_cdef = row_cdef + 1;
  curr_row_cdef = prev_row_cdef + nhfb + 2;
#endif
  for (pli = 0; pli < nplanes; pli++) {
    xdec[pli] = xd->plane[pli].subsampling_x;
    ydec[pli] = xd->plane[pli].subsampling_y;
    mi_wide_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_x;
    mi_high_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_y;
  }
#if !CONFIG_CDEF_SINGLEPASS
  stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * OD_FILT_HBORDER;
#else
  stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
#endif
  for (pli = 0; pli < nplanes; pli++) {
#if !CONFIG_CDEF_SINGLEPASS
    linebuf[pli] = aom_malloc(sizeof(*linebuf) * OD_FILT_VBORDER * stride);
#else
    linebuf[pli] = aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
#endif
    colbuf[pli] =
        aom_malloc(sizeof(*colbuf) *
#if !CONFIG_CDEF_SINGLEPASS
                   ((MAX_SB_SIZE << mi_high_l2[pli]) + 2 * OD_FILT_VBORDER) *
                   OD_FILT_HBORDER);
#else
                   ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
                   CDEF_HBORDER);
#endif
  }
#if !CONFIG_CDEF_SINGLEPASS
  for (sbr = 0; sbr < nvsb; sbr++) {
#else
  for (fbr = 0; fbr < nvfb; fbr++) {
#endif
    for (pli = 0; pli < nplanes; pli++) {
      const int block_height =
#if !CONFIG_CDEF_SINGLEPASS
          (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * OD_FILT_VBORDER;
      fill_rect(colbuf[pli], OD_FILT_HBORDER, block_height, OD_FILT_HBORDER,
                OD_DERING_VERY_LARGE);
#else
          (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
      fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER,
                CDEF_VERY_LARGE);
#endif
    }
#if !CONFIG_CDEF_SINGLEPASS
    dering_left = 1;
    for (sbc = 0; sbc < nhsb; sbc++) {
      int level, clpf_strength;
      int uv_level, uv_clpf_strength;
#else
    cdef_left = 1;
    for (fbc = 0; fbc < nhfb; fbc++) {
      int level, sec_strength;
      int uv_level, uv_sec_strength;
#endif
      int nhb, nvb;
      int cstart = 0;
#if !CONFIG_CDEF_SINGLEPASS
      curr_row_dering[sbc] = 0;
      if (cm->mi_grid_visible[MI_SIZE_64X64 * sbr * cm->mi_stride +
                              MI_SIZE_64X64 * sbc] == NULL ||
          cm->mi_grid_visible[MI_SIZE_64X64 * sbr * cm->mi_stride +
                              MI_SIZE_64X64 * sbc]
#else
      curr_row_cdef[fbc] = 0;
      if (cm->mi_grid_visible[MI_SIZE_64X64 * fbr * cm->mi_stride +
                              MI_SIZE_64X64 * fbc] == NULL ||
          cm->mi_grid_visible[MI_SIZE_64X64 * fbr * cm->mi_stride +
                              MI_SIZE_64X64 * fbc]
#endif
                  ->mbmi.cdef_strength == -1) {
#if !CONFIG_CDEF_SINGLEPASS
        dering_left = 0;
#else
        cdef_left = 0;
#endif
        continue;
      }
#if !CONFIG_CDEF_SINGLEPASS
      if (!dering_left) cstart = -OD_FILT_HBORDER;
      nhb = AOMMIN(MI_SIZE_64X64, cm->mi_cols - MI_SIZE_64X64 * sbc);
      nvb = AOMMIN(MI_SIZE_64X64, cm->mi_rows - MI_SIZE_64X64 * sbr);
#else
      if (!cdef_left) cstart = -CDEF_HBORDER;
      nhb = AOMMIN(MI_SIZE_64X64, cm->mi_cols - MI_SIZE_64X64 * fbc);
      nvb = AOMMIN(MI_SIZE_64X64, cm->mi_rows - MI_SIZE_64X64 * fbr);
#endif
      int tile_top, tile_left, tile_bottom, tile_right;
#if !CONFIG_CDEF_SINGLEPASS
      int mi_idx = MI_SIZE_64X64 * sbr * cm->mi_stride + MI_SIZE_64X64 * sbc;
#else
      int mi_idx = MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc;
#endif
      MODE_INFO *const mi_tl = cm->mi + mi_idx;
      BOUNDARY_TYPE boundary_tl = mi_tl->mbmi.boundary_info;
      tile_top = boundary_tl & TILE_ABOVE_BOUNDARY;
      tile_left = boundary_tl & TILE_LEFT_BOUNDARY;

#if !CONFIG_CDEF_SINGLEPASS
      if (sbr != nvsb - 1 &&
#else
      if (fbr != nvfb - 1 &&
#endif
          (&cm->mi[mi_idx + (MI_SIZE_64X64 - 1) * cm->mi_stride]))
        tile_bottom = cm->mi[mi_idx + (MI_SIZE_64X64 - 1) * cm->mi_stride]
                          .mbmi.boundary_info &
                      TILE_BOTTOM_BOUNDARY;
      else
        tile_bottom = 1;

#if !CONFIG_CDEF_SINGLEPASS
      if (sbc != nhsb - 1 && (&cm->mi[mi_idx + MI_SIZE_64X64 - 1]))
#else
      if (fbc != nhfb - 1 && (&cm->mi[mi_idx + MI_SIZE_64X64 - 1]))
#endif
        tile_right = cm->mi[mi_idx + MI_SIZE_64X64 - 1].mbmi.boundary_info &
                     TILE_RIGHT_BOUNDARY;
      else
        tile_right = 1;

      const int mbmi_cdef_strength =
#if !CONFIG_CDEF_SINGLEPASS
          cm->mi_grid_visible[MI_SIZE_64X64 * sbr * cm->mi_stride +
                              MI_SIZE_64X64 * sbc]
#else
          cm->mi_grid_visible[MI_SIZE_64X64 * fbr * cm->mi_stride +
                              MI_SIZE_64X64 * fbc]
#endif
              ->mbmi.cdef_strength;
#if !CONFIG_CDEF_SINGLEPASS
      level = cm->cdef_strengths[mbmi_cdef_strength] / CLPF_STRENGTHS;
      clpf_strength = cm->cdef_strengths[mbmi_cdef_strength] % CLPF_STRENGTHS;
      clpf_strength += clpf_strength == 3;
      uv_level = cm->cdef_uv_strengths[mbmi_cdef_strength] / CLPF_STRENGTHS;
      uv_clpf_strength =
          cm->cdef_uv_strengths[mbmi_cdef_strength] % CLPF_STRENGTHS;
      uv_clpf_strength += uv_clpf_strength == 3;
      if ((level == 0 && clpf_strength == 0 && uv_level == 0 &&
           uv_clpf_strength == 0) ||
          (dering_count = sb_compute_dering_list(
               cm, sbr * MI_SIZE_64X64, sbc * MI_SIZE_64X64, dlist,
               get_filter_skip(level) || get_filter_skip(uv_level))) == 0) {
        dering_left = 0;
#else
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
               cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist,
               (level & 1) || (uv_level & 1))) == 0) {
        cdef_left = 0;
#endif
        continue;
      }

#if !CONFIG_CDEF_SINGLEPASS
      curr_row_dering[sbc] = 1;
#else
      curr_row_cdef[fbc] = 1;
#endif
      for (pli = 0; pli < nplanes; pli++) {
#if !CONFIG_CDEF_SINGLEPASS
        uint16_t dst[MAX_SB_SIZE * MAX_SB_SIZE];
#endif
        int coffset;
        int rend, cend;
#if !CONFIG_CDEF_SINGLEPASS
        int clpf_damping = cm->cdef_clpf_damping;
        int dering_damping = cm->cdef_dering_damping;
#else
        int pri_damping = cm->cdef_pri_damping;
        int sec_damping = cm->cdef_sec_damping;
#endif
        int hsize = nhb << mi_wide_l2[pli];
        int vsize = nvb << mi_high_l2[pli];

        if (pli) {
#if !CONFIG_CDEF_SINGLEPASS
          if (chroma_dering)
#else
          if (chroma_cdef)
#endif
            level = uv_level;
          else
            level = 0;
#if !CONFIG_CDEF_SINGLEPASS
          clpf_strength = uv_clpf_strength;
#else
          sec_strength = uv_sec_strength;
#endif
        }

#if !CONFIG_CDEF_SINGLEPASS
        if (sbc == nhsb - 1)
#else
        if (fbc == nhfb - 1)
#endif
          cend = hsize;
        else
#if !CONFIG_CDEF_SINGLEPASS
          cend = hsize + OD_FILT_HBORDER;
#else
          cend = hsize + CDEF_HBORDER;
#endif

#if !CONFIG_CDEF_SINGLEPASS
        if (sbr == nvsb - 1)
#else
        if (fbr == nvfb - 1)
#endif
          rend = vsize;
        else
#if !CONFIG_CDEF_SINGLEPASS
          rend = vsize + OD_FILT_VBORDER;
#else
          rend = vsize + CDEF_VBORDER;
#endif

#if !CONFIG_CDEF_SINGLEPASS
        coffset = sbc * MI_SIZE_64X64 << mi_wide_l2[pli];
        if (sbc == nhsb - 1) {
#else
        coffset = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
        if (fbc == nhfb - 1) {
#endif
/* On the last superblock column, fill in the right border with
   OD_DERING_VERY_LARGE to avoid filtering with the outside. */
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[cend + OD_FILT_HBORDER], OD_FILT_BSTRIDE,
                    rend + OD_FILT_VBORDER, hsize + OD_FILT_HBORDER - cend,
                    OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[cend + CDEF_HBORDER], CDEF_BSTRIDE,
                    rend + CDEF_VBORDER, hsize + CDEF_HBORDER - cend,
                    CDEF_VERY_LARGE);
#endif
        }
#if !CONFIG_CDEF_SINGLEPASS
        if (sbr == nvsb - 1) {
#else
        if (fbr == nvfb - 1) {
#endif
/* On the last superblock row, fill in the bottom border with
   OD_DERING_VERY_LARGE to avoid filtering with the outside. */
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[(rend + OD_FILT_VBORDER) * OD_FILT_BSTRIDE],
                    OD_FILT_BSTRIDE, OD_FILT_VBORDER,
                    hsize + 2 * OD_FILT_HBORDER, OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[(rend + CDEF_VBORDER) * CDEF_BSTRIDE], CDEF_BSTRIDE,
                    CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
#endif
        }
/* Copy in the pixels we need from the current superblock for
   deringing.*/
#if !CONFIG_CDEF_SINGLEPASS
        copy_sb8_16(
            cm,
            &src[OD_FILT_VBORDER * OD_FILT_BSTRIDE + OD_FILT_HBORDER + cstart],
            OD_FILT_BSTRIDE, xd->plane[pli].dst.buf,
            (MI_SIZE_64X64 << mi_high_l2[pli]) * sbr, coffset + cstart,
            xd->plane[pli].dst.stride, rend, cend - cstart);
        if (!prev_row_dering[sbc]) {
          copy_sb8_16(
              cm, &src[OD_FILT_HBORDER], OD_FILT_BSTRIDE,
              xd->plane[pli].dst.buf,
              (MI_SIZE_64X64 << mi_high_l2[pli]) * sbr - OD_FILT_VBORDER,
              coffset, xd->plane[pli].dst.stride, OD_FILT_VBORDER, hsize);
        } else if (sbr > 0) {
          copy_rect(&src[OD_FILT_HBORDER], OD_FILT_BSTRIDE,
                    &linebuf[pli][coffset], stride, OD_FILT_VBORDER, hsize);
#else
        copy_sb8_16(cm,
                    &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                    CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr, coffset + cstart,
                    xd->plane[pli].dst.stride, rend, cend - cstart);
        if (!prev_row_cdef[fbc]) {
          copy_sb8_16(cm, &src[CDEF_HBORDER], CDEF_BSTRIDE,
                      xd->plane[pli].dst.buf,
                      (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                      coffset, xd->plane[pli].dst.stride, CDEF_VBORDER, hsize);
        } else if (fbr > 0) {
          copy_rect(&src[CDEF_HBORDER], CDEF_BSTRIDE, &linebuf[pli][coffset],
                    stride, CDEF_VBORDER, hsize);
#endif
        } else {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[OD_FILT_HBORDER], OD_FILT_BSTRIDE, OD_FILT_VBORDER,
                    hsize, OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER, hsize,
                    CDEF_VERY_LARGE);
#endif
        }
#if !CONFIG_CDEF_SINGLEPASS
        if (!prev_row_dering[sbc - 1]) {
          copy_sb8_16(
              cm, src, OD_FILT_BSTRIDE, xd->plane[pli].dst.buf,
              (MI_SIZE_64X64 << mi_high_l2[pli]) * sbr - OD_FILT_VBORDER,
              coffset - OD_FILT_HBORDER, xd->plane[pli].dst.stride,
              OD_FILT_VBORDER, OD_FILT_HBORDER);
        } else if (sbr > 0 && sbc > 0) {
          copy_rect(src, OD_FILT_BSTRIDE,
                    &linebuf[pli][coffset - OD_FILT_HBORDER], stride,
                    OD_FILT_VBORDER, OD_FILT_HBORDER);
#else
        if (!prev_row_cdef[fbc - 1]) {
          copy_sb8_16(cm, src, CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                      (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                      coffset - CDEF_HBORDER, xd->plane[pli].dst.stride,
                      CDEF_VBORDER, CDEF_HBORDER);
        } else if (fbr > 0 && fbc > 0) {
          copy_rect(src, CDEF_BSTRIDE, &linebuf[pli][coffset - CDEF_HBORDER],
                    stride, CDEF_VBORDER, CDEF_HBORDER);
#endif
        } else {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(src, OD_FILT_BSTRIDE, OD_FILT_VBORDER, OD_FILT_HBORDER,
                    OD_DERING_VERY_LARGE);
#else
          fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, CDEF_HBORDER,
                    CDEF_VERY_LARGE);
#endif
        }
#if !CONFIG_CDEF_SINGLEPASS
        if (!prev_row_dering[sbc + 1]) {
          copy_sb8_16(
              cm, &src[OD_FILT_HBORDER + (nhb << mi_wide_l2[pli])],
              OD_FILT_BSTRIDE, xd->plane[pli].dst.buf,
              (MI_SIZE_64X64 << mi_high_l2[pli]) * sbr - OD_FILT_VBORDER,
              coffset + hsize, xd->plane[pli].dst.stride, OD_FILT_VBORDER,
              OD_FILT_HBORDER);
        } else if (sbr > 0 && sbc < nhsb - 1) {
          copy_rect(&src[hsize + OD_FILT_HBORDER], OD_FILT_BSTRIDE,
                    &linebuf[pli][coffset + hsize], stride, OD_FILT_VBORDER,
                    OD_FILT_HBORDER);
#else
        if (!prev_row_cdef[fbc + 1]) {
          copy_sb8_16(cm, &src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                      CDEF_BSTRIDE, xd->plane[pli].dst.buf,
                      (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                      coffset + hsize, xd->plane[pli].dst.stride, CDEF_VBORDER,
                      CDEF_HBORDER);
        } else if (fbr > 0 && fbc < nhfb - 1) {
          copy_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE,
                    &linebuf[pli][coffset + hsize], stride, CDEF_VBORDER,
                    CDEF_HBORDER);
#endif
        } else {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[hsize + OD_FILT_HBORDER], OD_FILT_BSTRIDE,
                    OD_FILT_VBORDER, OD_FILT_HBORDER, OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER,
                    CDEF_HBORDER, CDEF_VERY_LARGE);
#endif
        }
#if !CONFIG_CDEF_SINGLEPASS
        if (dering_left) {
#else
        if (cdef_left) {
#endif
/* If we deringed the superblock on the left then we need to copy in
   saved pixels. */
#if !CONFIG_CDEF_SINGLEPASS
          copy_rect(src, OD_FILT_BSTRIDE, colbuf[pli], OD_FILT_HBORDER,
                    rend + OD_FILT_VBORDER, OD_FILT_HBORDER);
#else
          copy_rect(src, CDEF_BSTRIDE, colbuf[pli], CDEF_HBORDER,
                    rend + CDEF_VBORDER, CDEF_HBORDER);
#endif
        }
/* Saving pixels in case we need to dering the superblock on the
    right. */
#if !CONFIG_CDEF_SINGLEPASS
        copy_rect(colbuf[pli], OD_FILT_HBORDER, src + hsize, OD_FILT_BSTRIDE,
                  rend + OD_FILT_VBORDER, OD_FILT_HBORDER);
#else
        copy_rect(colbuf[pli], CDEF_HBORDER, src + hsize, CDEF_BSTRIDE,
                  rend + CDEF_VBORDER, CDEF_HBORDER);
#endif
        copy_sb8_16(
            cm, &linebuf[pli][coffset], stride, xd->plane[pli].dst.buf,
#if !CONFIG_CDEF_SINGLEPASS
            (MI_SIZE_64X64 << mi_high_l2[pli]) * (sbr + 1) - OD_FILT_VBORDER,
            coffset, xd->plane[pli].dst.stride, OD_FILT_VBORDER, hsize);
#else
            (MI_SIZE_64X64 << mi_high_l2[pli]) * (fbr + 1) - CDEF_VBORDER,
            coffset, xd->plane[pli].dst.stride, CDEF_VBORDER, hsize);
#endif

        if (tile_top) {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(src, OD_FILT_BSTRIDE, OD_FILT_VBORDER,
                    hsize + 2 * OD_FILT_HBORDER, OD_DERING_VERY_LARGE);
#else
          fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, hsize + 2 * CDEF_HBORDER,
                    CDEF_VERY_LARGE);
#endif
        }
        if (tile_left) {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(src, OD_FILT_BSTRIDE, vsize + 2 * OD_FILT_VBORDER,
                    OD_FILT_HBORDER, OD_DERING_VERY_LARGE);
#else
          fill_rect(src, CDEF_BSTRIDE, vsize + 2 * CDEF_VBORDER, CDEF_HBORDER,
                    CDEF_VERY_LARGE);
#endif
        }
        if (tile_bottom) {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[(vsize + OD_FILT_VBORDER) * OD_FILT_BSTRIDE],
                    OD_FILT_BSTRIDE, OD_FILT_VBORDER,
                    hsize + 2 * OD_FILT_HBORDER, OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[(vsize + CDEF_VBORDER) * CDEF_BSTRIDE], CDEF_BSTRIDE,
                    CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
#endif
        }
        if (tile_right) {
#if !CONFIG_CDEF_SINGLEPASS
          fill_rect(&src[hsize + OD_FILT_HBORDER], OD_FILT_BSTRIDE,
                    vsize + 2 * OD_FILT_VBORDER, OD_FILT_HBORDER,
                    OD_DERING_VERY_LARGE);
#else
          fill_rect(&src[hsize + CDEF_HBORDER], CDEF_BSTRIDE,
                    vsize + 2 * CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
#endif
        }
#if CONFIG_HIGHBITDEPTH
        if (cm->use_highbitdepth) {
#if !CONFIG_CDEF_SINGLEPASS
          od_dering(
              (uint8_t *)&CONVERT_TO_SHORTPTR(
#else
          cdef_filter_fb(
              NULL,
              &CONVERT_TO_SHORTPTR(
#endif
                  xd->plane[pli]
                      .dst.buf)[xd->plane[pli].dst.stride *
#if !CONFIG_CDEF_SINGLEPASS
                                    (MI_SIZE_64X64 * sbr << mi_high_l2[pli]) +
                                (sbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
              xd->plane[pli].dst.stride, dst,
              &src[OD_FILT_VBORDER * OD_FILT_BSTRIDE + OD_FILT_HBORDER],
              xdec[pli], ydec[pli], dir, NULL, var, pli, dlist, dering_count,
              level, clpf_strength, clpf_damping, dering_damping, coeff_shift,
              0, 1);
#else
                                    (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                                (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
              xd->plane[pli].dst.stride,
              &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER], xdec[pli],
              ydec[pli], dir, NULL, var, pli, dlist, cdef_count, level,
              sec_strength, pri_damping, sec_damping, coeff_shift);
#endif
        } else {
#endif
#if !CONFIG_CDEF_SINGLEPASS
          od_dering(&xd->plane[pli]
                         .dst.buf[xd->plane[pli].dst.stride *
                                      (MI_SIZE_64X64 * sbr << mi_high_l2[pli]) +
                                  (sbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                    xd->plane[pli].dst.stride, dst,
                    &src[OD_FILT_VBORDER * OD_FILT_BSTRIDE + OD_FILT_HBORDER],
                    xdec[pli], ydec[pli], dir, NULL, var, pli, dlist,
                    dering_count, level, clpf_strength, clpf_damping,
                    dering_damping, coeff_shift, 0, 0);
#else
        cdef_filter_fb(
            &xd->plane[pli]
                 .dst.buf[xd->plane[pli].dst.stride *
                              (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                          (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
            NULL, xd->plane[pli].dst.stride,
            &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER], xdec[pli],
            ydec[pli], dir, NULL, var, pli, dlist, cdef_count, level,
            sec_strength, pri_damping, sec_damping, coeff_shift);
#endif

#if CONFIG_HIGHBITDEPTH
        }
#endif
      }
#if !CONFIG_CDEF_SINGLEPASS
      dering_left = 1;
#else
      cdef_left = 1;
#endif
    }
    {
      unsigned char *tmp;
#if !CONFIG_CDEF_SINGLEPASS
      tmp = prev_row_dering;
      prev_row_dering = curr_row_dering;
      curr_row_dering = tmp;
#else
      tmp = prev_row_cdef;
      prev_row_cdef = curr_row_cdef;
      curr_row_cdef = tmp;
#endif
    }
  }
#if !CONFIG_CDEF_SINGLEPASS
  aom_free(row_dering);
#else
  aom_free(row_cdef);
#endif
  for (pli = 0; pli < nplanes; pli++) {
    aom_free(linebuf[pli]);
    aom_free(colbuf[pli]);
  }
}
