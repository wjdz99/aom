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

#include <math.h>
#include <string.h>

#include "config/aom_dsp_rtcd.h"
#include "config/aom_scale_rtcd.h"

#include "aom/aom_integer.h"
#include "aom_ports/system_state.h"
#include "av1/common/av1_common_int.h"
#include "av1/common/reconinter.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/pickcdef.h"

#define REDUCED_PRI_STRENGTHS_LVL1 8
#define REDUCED_PRI_STRENGTHS_LVL2 5
#define REDUCED_SEC_STRENGTHS_LVL3 2

#define REDUCED_TOTAL_STRENGTHS_LVL1 \
  (REDUCED_PRI_STRENGTHS_LVL1 * CDEF_SEC_STRENGTHS)
#define REDUCED_TOTAL_STRENGTHS_LVL2 \
  (REDUCED_PRI_STRENGTHS_LVL2 * CDEF_SEC_STRENGTHS)
#define REDUCED_TOTAL_STRENGTHS_LVL3 \
  (REDUCED_PRI_STRENGTHS_LVL2 * REDUCED_SEC_STRENGTHS_LVL3)
#define TOTAL_STRENGTHS (CDEF_PRI_STRENGTHS * CDEF_SEC_STRENGTHS)

static const int priconv_lvl1[REDUCED_PRI_STRENGTHS_LVL1] = { 0, 1, 2,  3,
                                                              5, 7, 10, 13 };
static const int priconv_lvl2[REDUCED_PRI_STRENGTHS_LVL2] = { 0, 2, 4, 8, 14 };
static const int secconv_lvl3[REDUCED_SEC_STRENGTHS_LVL3] = { 0, 2 };
static const int nb_cdef_strengths[CDEF_PICK_METHODS] = {
  TOTAL_STRENGTHS, REDUCED_TOTAL_STRENGTHS_LVL1, REDUCED_TOTAL_STRENGTHS_LVL2,
  REDUCED_TOTAL_STRENGTHS_LVL3, TOTAL_STRENGTHS
};

#if CONFIG_CC_CDEF
static const int nb_cc_cdef_strengths[CDEF_PICK_METHODS] = {
  TOTAL_CCCDEF_STRENGTHS, REDUCED_TOTAL_STRENGTHS_LVL1,
  REDUCED_TOTAL_STRENGTHS_LVL2, REDUCED_TOTAL_STRENGTHS_LVL3,
  TOTAL_CCCDEF_STRENGTHS
};
#endif

// Get primary and secondary filter strength for the given strength index and
// search method
static INLINE void get_cdef_filter_strengths(CDEF_PICK_METHOD pick_method,
                                             int *pri_strength,
                                             int *sec_strength,
                                             int strength_idx) {
  const int tot_sec_filter = (pick_method == CDEF_FAST_SEARCH_LVL3)
                                 ? REDUCED_SEC_STRENGTHS_LVL3
                                 : CDEF_SEC_STRENGTHS;
  const int pri_idx = strength_idx / tot_sec_filter;
  const int sec_idx = strength_idx % tot_sec_filter;
  *pri_strength = pri_idx;
  *sec_strength = sec_idx;
  if (pick_method == CDEF_FULL_SEARCH) return;

  switch (pick_method) {
    case CDEF_FAST_SEARCH_LVL1: *pri_strength = priconv_lvl1[pri_idx]; break;
    case CDEF_FAST_SEARCH_LVL2: *pri_strength = priconv_lvl2[pri_idx]; break;
    case CDEF_FAST_SEARCH_LVL3:
      *pri_strength = priconv_lvl2[pri_idx];
      *sec_strength = secconv_lvl3[sec_idx];
      break;
    default: assert(0 && "Invalid CDEF search method");
  }
}

// Store CDEF filter strength calculated from strength index for given search
// method
#define STORE_CDEF_FILTER_STRENGTH(cdef_strength, pick_method, strength_idx) \
  get_cdef_filter_strengths((pick_method), &pri_strength, &sec_strength,     \
                            (strength_idx));                                 \
  cdef_strength = pri_strength * CDEF_SEC_STRENGTHS + sec_strength;

/* Search for the best strength to add as an option, knowing we
   already selected nb_strengths options. */
static uint64_t search_one(int *lev, int nb_strengths,
                           uint64_t mse[][TOTAL_STRENGTHS], int sb_count,
                           CDEF_PICK_METHOD pick_method) {
  uint64_t tot_mse[TOTAL_STRENGTHS];
  const int total_strengths = nb_cdef_strengths[pick_method];
  int i, j;
  uint64_t best_tot_mse = (uint64_t)1 << 63;
  int best_id = 0;
  memset(tot_mse, 0, sizeof(tot_mse));
  for (i = 0; i < sb_count; i++) {
    int gi;
    uint64_t best_mse = (uint64_t)1 << 63;
    /* Find best mse among already selected options. */
    for (gi = 0; gi < nb_strengths; gi++) {
      if (mse[i][lev[gi]] < best_mse) {
        best_mse = mse[i][lev[gi]];
      }
    }
    /* Find best mse when adding each possible new option. */
    for (j = 0; j < total_strengths; j++) {
      uint64_t best = best_mse;
      if (mse[i][j] < best) best = mse[i][j];
      tot_mse[j] += best;
    }
  }
  for (j = 0; j < total_strengths; j++) {
    if (tot_mse[j] < best_tot_mse) {
      best_tot_mse = tot_mse[j];
      best_id = j;
    }
  }
  lev[nb_strengths] = best_id;
  return best_tot_mse;
}

/* Search for the best luma+chroma strength to add as an option, knowing we
   already selected nb_strengths options. */
static uint64_t search_one_dual(int *lev0, int *lev1, int nb_strengths,
                                uint64_t (**mse)[TOTAL_STRENGTHS], int sb_count,
                                CDEF_PICK_METHOD pick_method) {
  uint64_t tot_mse[TOTAL_STRENGTHS][TOTAL_STRENGTHS];
  int i, j;
  uint64_t best_tot_mse = (uint64_t)1 << 63;
  int best_id0 = 0;
  int best_id1 = 0;
  const int total_strengths = nb_cdef_strengths[pick_method];
  memset(tot_mse, 0, sizeof(tot_mse));
  for (i = 0; i < sb_count; i++) {
    int gi;
    uint64_t best_mse = (uint64_t)1 << 63;
    /* Find best mse among already selected options. */
    for (gi = 0; gi < nb_strengths; gi++) {
      uint64_t curr = mse[0][i][lev0[gi]];
      curr += mse[1][i][lev1[gi]];

      if (curr < best_mse) {
        best_mse = curr;
      }
    }
    /* Find best mse when adding each possible new option. */
    for (j = 0; j < total_strengths; j++) {
      int k;

      for (k = 0; k < total_strengths; k++) {
        uint64_t best = best_mse;
        uint64_t curr = mse[0][i][j];
        curr += mse[1][i][k];
        if (curr < best) best = curr;
        tot_mse[j][k] += best;
      }
    }
  }
  for (j = 0; j < total_strengths; j++) {
    int k;

    for (k = 0; k < total_strengths; k++) {
      if (tot_mse[j][k] < best_tot_mse) {
        best_tot_mse = tot_mse[j][k];
        best_id0 = j;
        best_id1 = k;
      }
    }
  }
  lev0[nb_strengths] = best_id0;
  lev1[nb_strengths] = best_id1;

  return best_tot_mse;
}

/* Search for the set of strengths that minimizes mse. */
static uint64_t joint_strength_search(int *best_lev, int nb_strengths,
                                      uint64_t mse[][TOTAL_STRENGTHS],
                                      int sb_count,
                                      CDEF_PICK_METHOD pick_method) {
  uint64_t best_tot_mse;
  int fast = (pick_method >= CDEF_FAST_SEARCH_LVL1 &&
              pick_method <= CDEF_FAST_SEARCH_LVL3);
  int i;
  best_tot_mse = (uint64_t)1 << 63;
  /* Greedy search: add one strength options at a time. */
  for (i = 0; i < nb_strengths; i++) {
    best_tot_mse = search_one(best_lev, i, mse, sb_count, pick_method);
  }
  /* Trying to refine the greedy search by reconsidering each
     already-selected option. */
  if (!fast) {
    for (i = 0; i < 4 * nb_strengths; i++) {
      int j;
      for (j = 0; j < nb_strengths - 1; j++) best_lev[j] = best_lev[j + 1];
      best_tot_mse =
          search_one(best_lev, nb_strengths - 1, mse, sb_count, pick_method);
    }
  }
  return best_tot_mse;
}

/* Search for the set of luma+chroma strengths that minimizes mse. */
static uint64_t joint_strength_search_dual(int *best_lev0, int *best_lev1,
                                           int nb_strengths,
                                           uint64_t (**mse)[TOTAL_STRENGTHS],
                                           int sb_count,
                                           CDEF_PICK_METHOD pick_method) {
  uint64_t best_tot_mse;
  int i;
  best_tot_mse = (uint64_t)1 << 63;
  /* Greedy search: add one strength options at a time. */
  for (i = 0; i < nb_strengths; i++) {
    best_tot_mse =
        search_one_dual(best_lev0, best_lev1, i, mse, sb_count, pick_method);
  }
  /* Trying to refine the greedy search by reconsidering each
     already-selected option. */
  for (i = 0; i < 4 * nb_strengths; i++) {
    int j;
    for (j = 0; j < nb_strengths - 1; j++) {
      best_lev0[j] = best_lev0[j + 1];
      best_lev1[j] = best_lev1[j + 1];
    }
    best_tot_mse = search_one_dual(best_lev0, best_lev1, nb_strengths - 1, mse,
                                   sb_count, pick_method);
  }
  return best_tot_mse;
}

typedef void (*copy_fn_t)(uint16_t *dst, int dstride, const void *src,
                          int src_voffset, int src_hoffset, int sstride,
                          int vsize, int hsize);
typedef uint64_t (*compute_cdef_dist_t)(void *dst, int dstride, uint16_t *src,
                                        cdef_list *dlist, int cdef_count,
                                        BLOCK_SIZE bsize, int coeff_shift,
                                        int row, int col);

#if CONFIG_CC_CDEF
typedef uint64_t (*compute_cdef_dist_direction_t)(
    void *dst, int dstride, uint16_t *src, cdef_list *dlist, int cdef_count,
    BLOCK_SIZE bsize, int coeff_shift, int row, int col, int direction,
    int dir[CDEF_NBLOCKS][CDEF_NBLOCKS]);
#endif

static void copy_sb16_16_highbd(uint16_t *dst, int dstride, const void *src,
                                int src_voffset, int src_hoffset, int sstride,
                                int vsize, int hsize) {
  int r;
  const uint16_t *src16 = CONVERT_TO_SHORTPTR((uint8_t *)src);
  const uint16_t *base = &src16[src_voffset * sstride + src_hoffset];
  for (r = 0; r < vsize; r++)
    memcpy(dst + r * dstride, base + r * sstride, hsize * sizeof(*base));
}

static void copy_sb16_16(uint16_t *dst, int dstride, const void *src,
                         int src_voffset, int src_hoffset, int sstride,
                         int vsize, int hsize) {
  int r, c;
  const uint8_t *src8 = (uint8_t *)src;
  const uint8_t *base = &src8[src_voffset * sstride + src_hoffset];
  for (r = 0; r < vsize; r++)
    for (c = 0; c < hsize; c++)
      dst[r * dstride + c] = (uint16_t)base[r * sstride + c];
}

static INLINE void init_src_params(int *src_stride, int *width, int *height,
                                   int *width_log2, int *height_log2,
                                   BLOCK_SIZE bsize) {
  *src_stride = block_size_wide[bsize];
  *width = block_size_wide[bsize];
  *height = block_size_high[bsize];
  *width_log2 = MI_SIZE_LOG2 + mi_size_wide_log2[bsize];
  *height_log2 = MI_SIZE_LOG2 + mi_size_wide_log2[bsize];
}

/* Compute MSE only on the blocks we filtered. */
static uint64_t compute_cdef_dist_highbd(void *dst, int dstride, uint16_t *src,
                                         cdef_list *dlist, int cdef_count,
                                         BLOCK_SIZE bsize, int coeff_shift,
                                         int row, int col) {
  assert(bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
         bsize == BLOCK_8X8);
  uint64_t sum = 0;
  int bi, bx, by;
  uint16_t *dst16 = CONVERT_TO_SHORTPTR((uint8_t *)dst);
  uint16_t *dst_buff = &dst16[row * dstride + col];
  int src_stride, width, height, width_log2, height_log2;
  init_src_params(&src_stride, &width, &height, &width_log2, &height_log2,
                  bsize);
  for (bi = 0; bi < cdef_count; bi++) {
    by = dlist[bi].by;
    bx = dlist[bi].bx;
    sum += aom_mse_wxh_16bit_highbd(
        &dst_buff[(by << height_log2) * dstride + (bx << width_log2)], dstride,
        &src[bi << (height_log2 + width_log2)], src_stride, width, height);
  }
  return sum >> 2 * coeff_shift;
}

static uint64_t compute_cdef_dist(void *dst, int dstride, uint16_t *src,
                                  cdef_list *dlist, int cdef_count,
                                  BLOCK_SIZE bsize, int coeff_shift, int row,
                                  int col) {
  assert(bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
         bsize == BLOCK_8X8);
  uint64_t sum = 0;
  int bi, bx, by;
  uint8_t *dst8 = (uint8_t *)dst;
  uint8_t *dst_buff = &dst8[row * dstride + col];
  int src_stride, width, height, width_log2, height_log2;
  init_src_params(&src_stride, &width, &height, &width_log2, &height_log2,
                  bsize);
  for (bi = 0; bi < cdef_count; bi++) {
    by = dlist[bi].by;
    bx = dlist[bi].bx;
    sum += aom_mse_wxh_16bit(
        &dst_buff[(by << height_log2) * dstride + (bx << width_log2)], dstride,
        &src[bi << (height_log2 + width_log2)], src_stride, width, height);
  }
  return sum >> 2 * coeff_shift;
}

#if CONFIG_CC_CDEF
/* Compute MSE only on the blocks we filtered. */
static uint64_t compute_cdef_dist_highbd_direction(
    void *dst, int dstride, uint16_t *src, cdef_list *dlist, int cdef_count,
    BLOCK_SIZE bsize, int coeff_shift, int row, int col, int direction_target,
    int dir[CDEF_NBLOCKS][CDEF_NBLOCKS]) {
  assert(bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
         bsize == BLOCK_8X8);
  uint64_t sum = 0;
  int bi, bx, by;
  uint16_t *dst16 = CONVERT_TO_SHORTPTR((uint8_t *)dst);
  uint16_t *dst_buff = &dst16[row * dstride + col];
  int src_stride, width, height, width_log2, height_log2;
  init_src_params(&src_stride, &width, &height, &width_log2, &height_log2,
                  bsize);
  for (bi = 0; bi < cdef_count; bi++) {
    by = dlist[bi].by;
    bx = dlist[bi].bx;
    if (dir[by][bx] == direction_target) {
      sum += aom_mse_wxh_16bit_highbd(
          &dst_buff[(by << height_log2) * dstride + (bx << width_log2)],
          dstride, &src[bi << (height_log2 + width_log2)], src_stride, width,
          height);
    }
  }
  return sum >> 2 * coeff_shift;
}

static uint64_t compute_cdef_dist_direction(
    void *dst, int dstride, uint16_t *src, cdef_list *dlist, int cdef_count,
    BLOCK_SIZE bsize, int coeff_shift, int row, int col, int direction_target,
    int dir[CDEF_NBLOCKS][CDEF_NBLOCKS]) {
  assert(bsize == BLOCK_4X4 || bsize == BLOCK_4X8 || bsize == BLOCK_8X4 ||
         bsize == BLOCK_8X8);
  uint64_t sum = 0;
  int bi, bx, by;
  uint8_t *dst8 = (uint8_t *)dst;
  uint8_t *dst_buff = &dst8[row * dstride + col];
  int src_stride, width, height, width_log2, height_log2;
  init_src_params(&src_stride, &width, &height, &width_log2, &height_log2,
                  bsize);
  for (bi = 0; bi < cdef_count; bi++) {
    by = dlist[bi].by;
    bx = dlist[bi].bx;
    if (dir[by][bx] == direction_target) {
      sum += aom_mse_wxh_16bit(
          &dst_buff[(by << height_log2) * dstride + (bx << width_log2)],
          dstride, &src[bi << (height_log2 + width_log2)], src_stride, width,
          height);
    }
  }
  return sum >> 2 * coeff_shift;
}

#endif

static int sb_all_skip(const CommonModeInfoParams *const mi_params, int mi_row,
                       int mi_col) {
  const int maxr = AOMMIN(mi_params->mi_rows - mi_row, MI_SIZE_64X64);
  const int maxc = AOMMIN(mi_params->mi_cols - mi_col, MI_SIZE_64X64);
  const int stride = mi_params->mi_stride;
  MB_MODE_INFO **mbmi = mi_params->mi_grid_base + mi_row * stride + mi_col;
  for (int r = 0; r < maxr; ++r, mbmi += stride) {
    for (int c = 0; c < maxc; ++c) {
      if (!mbmi[c]->skip_txfm) return 0;
    }
  }
  return 1;
}

static void pick_cdef_from_qp(AV1_COMMON *const cm) {
  const int bd = cm->seq_params.bit_depth;
  const int q =
      av1_ac_quant_QTX(cm->quant_params.base_qindex, 0, bd) >> (bd - 8);
  CdefInfo *const cdef_info = &cm->cdef_info;
  cdef_info->cdef_bits = 0;
  cdef_info->nb_cdef_strengths = 1;
  cdef_info->cdef_damping = 3 + (cm->quant_params.base_qindex >> 6);

  int predicted_y_f1 = 0;
  int predicted_y_f2 = 0;
  int predicted_uv_f1 = 0;
  int predicted_uv_f2 = 0;
  aom_clear_system_state();
  if (!frame_is_intra_only(cm)) {
    predicted_y_f1 = clamp((int)roundf(q * q * -0.0000023593946f +
                                       q * 0.0068615186f + 0.02709886f),
                           0, 15);
    predicted_y_f2 = clamp((int)roundf(q * q * -0.00000057629734f +
                                       q * 0.0013993345f + 0.03831067f),
                           0, 3);
    predicted_uv_f1 = clamp((int)roundf(q * q * -0.0000007095069f +
                                        q * 0.0034628846f + 0.00887099f),
                            0, 15);
    predicted_uv_f2 = clamp((int)roundf(q * q * 0.00000023874085f +
                                        q * 0.00028223585f + 0.05576307f),
                            0, 3);
  } else {
    predicted_y_f1 = clamp(
        (int)roundf(q * q * 0.0000033731974f + q * 0.008070594f + 0.0187634f),
        0, 15);
    predicted_y_f2 = clamp(
        (int)roundf(q * q * 0.0000029167343f + q * 0.0027798624f + 0.0079405f),
        0, 3);
    predicted_uv_f1 = clamp(
        (int)roundf(q * q * -0.0000130790995f + q * 0.012892405f - 0.00748388f),
        0, 15);
    predicted_uv_f2 = clamp((int)roundf(q * q * 0.0000032651783f +
                                        q * 0.00035520183f + 0.00228092f),
                            0, 3);
  }
  cdef_info->cdef_strengths[0] =
      predicted_y_f1 * CDEF_SEC_STRENGTHS + predicted_y_f2;

  cdef_info->cdef_uv_strengths[0] =
      predicted_uv_f1 * CDEF_SEC_STRENGTHS + predicted_uv_f2;

  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int nvfb = (mi_params->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  const int nhfb = (mi_params->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  MB_MODE_INFO **mbmi = mi_params->mi_grid_base;
  for (int r = 0; r < nvfb; ++r) {
    for (int c = 0; c < nhfb; ++c) {
      mbmi[MI_SIZE_64X64 * c]->cdef_strength = 0;

#if CONFIG_CC_CDEF
      mbmi[MI_SIZE_64X64 * c]->cc_cdef_strength_index_fb[0] = 0;
      mbmi[MI_SIZE_64X64 * c]->cc_cdef_strength_index_fb[1] = 0;
#endif
    }
    mbmi += MI_SIZE_64X64 * mi_params->mi_stride;
  }
}

void av1_cdef_search(const YV12_BUFFER_CONFIG *frame,
                     const YV12_BUFFER_CONFIG *ref, AV1_COMMON *cm,
                     MACROBLOCKD *xd, CDEF_PICK_METHOD pick_method,
                     int rdmult) {
  if (pick_method == CDEF_PICK_FROM_Q) {
    pick_cdef_from_qp(cm);
    return;
  }

  cdef_list dlist[MI_SIZE_128X128 * MI_SIZE_128X128];

  int dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int var[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };

  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int nvfb = (mi_params->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  const int nhfb = (mi_params->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  int *sb_index = aom_malloc(nvfb * nhfb * sizeof(*sb_index));
  const int damping = 3 + (cm->quant_params.base_qindex >> 6);
  const int fast = (pick_method >= CDEF_FAST_SEARCH_LVL1 &&
                    pick_method <= CDEF_FAST_SEARCH_LVL3);
  const int total_strengths = nb_cdef_strengths[pick_method];
  DECLARE_ALIGNED(32, uint16_t, tmp_dst[1 << (MAX_SB_SIZE_LOG2 * 2)]);
  const int num_planes = av1_num_planes(cm);
  av1_setup_dst_planes(xd->plane, cm->seq_params.sb_size, frame, 0, 0, 0,
                       num_planes);

  uint64_t(*mse[2])[TOTAL_STRENGTHS];

  mse[0] = aom_malloc(sizeof(**mse) * nvfb * nhfb);
  mse[1] = aom_malloc(sizeof(**mse) * nvfb * nhfb);

  int bsize[3];
  int mi_wide_l2[3];
  int mi_high_l2[3];
  int xdec[3];
  int ydec[3];
  uint8_t *ref_buffer[3] = { ref->y_buffer, ref->u_buffer, ref->v_buffer };
  int ref_stride[3] = { ref->y_stride, ref->uv_stride, ref->uv_stride };

  for (int pli = 0; pli < num_planes; pli++) {
    xdec[pli] = xd->plane[pli].subsampling_x;
    ydec[pli] = xd->plane[pli].subsampling_y;
    bsize[pli] = ydec[pli] ? (xdec[pli] ? BLOCK_4X4 : BLOCK_8X4)
                           : (xdec[pli] ? BLOCK_4X8 : BLOCK_8X8);
    mi_wide_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_x;
    mi_high_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_y;
  }

  copy_fn_t copy_fn;
  compute_cdef_dist_t compute_cdef_dist_fn;

  if (cm->seq_params.use_highbitdepth) {
    copy_fn = copy_sb16_16_highbd;
    compute_cdef_dist_fn = compute_cdef_dist_highbd;
  } else {
    copy_fn = copy_sb16_16;
    compute_cdef_dist_fn = compute_cdef_dist;
  }

  DECLARE_ALIGNED(32, uint16_t, inbuf[CDEF_INBUF_SIZE]);
  uint16_t *const in = inbuf + CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER;
  const int coeff_shift = AOMMAX(cm->seq_params.bit_depth - 8, 0);
  int sb_count = 0;
  for (int fbr = 0; fbr < nvfb; ++fbr) {
    for (int fbc = 0; fbc < nhfb; ++fbc) {
      // No filtering if the entire filter block is skipped
      if (sb_all_skip(mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64))
        continue;

      const MB_MODE_INFO *const mbmi =
          mi_params->mi_grid_base[MI_SIZE_64X64 * fbr * mi_params->mi_stride +
                                  MI_SIZE_64X64 * fbc];
      if (((fbc & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64)) ||
          ((fbr & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_64X128)))
        continue;

      int nhb = AOMMIN(MI_SIZE_64X64, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
      int nvb = AOMMIN(MI_SIZE_64X64, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
      int hb_step = 1;
      int vb_step = 1;
      BLOCK_SIZE bs;
      if (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64 ||
          mbmi->sb_type == BLOCK_64X128) {
        bs = mbmi->sb_type;
        if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
          nhb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
          hb_step = 2;
        }
        if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
          nvb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
          vb_step = 2;
        }
      } else {
        bs = BLOCK_64X64;
      }

      const int cdef_count = av1_cdef_compute_sb_list(
          mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, bs);

      const int yoff = CDEF_VBORDER * (fbr != 0);
      const int xoff = CDEF_HBORDER * (fbc != 0);

      int dirinit = 0;

      for (int pli = 0; pli < num_planes; pli++) {
        for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf[i] = CDEF_VERY_LARGE;
        /* We avoid filtering the pixels for which some of the pixels to
           average are outside the frame. We could change the filter instead,
           but it would add special cases for any future vectorization. */
        const int ysize = (nvb << mi_high_l2[pli]) +
                          CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
        const int xsize = (nhb << mi_wide_l2[pli]) +
                          CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
        const int row = fbr * MI_SIZE_64X64 << mi_high_l2[pli];
        const int col = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
        for (int gi = 0; gi < total_strengths; gi++) {
          int pri_strength, sec_strength;
          get_cdef_filter_strengths(pick_method, &pri_strength, &sec_strength,
                                    gi);
          copy_fn(&in[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
                  xd->plane[pli].dst.buf, row - yoff, col - xoff,
                  xd->plane[pli].dst.stride, ysize, xsize);

          av1_cdef_filter_fb(
              NULL, tmp_dst, CDEF_BSTRIDE, in, xdec[pli], ydec[pli], dir,
              &dirinit, var, pli, dlist, cdef_count, pri_strength,
              sec_strength + (sec_strength == 3), damping, coeff_shift
#if CONFIG_CC_CDEF
              ,
              1
#endif
          );

          const uint64_t curr_mse = compute_cdef_dist_fn(
              ref_buffer[pli], ref_stride[pli], tmp_dst, dlist, cdef_count,
              bsize[pli], coeff_shift, row, col);

          if (pli < 2)
            mse[pli][sb_count][gi] = curr_mse;
          else
            mse[1][sb_count][gi] += curr_mse;
        }
      }
      sb_index[sb_count++] =
          MI_SIZE_64X64 * fbr * mi_params->mi_stride + MI_SIZE_64X64 * fbc;
    }
  }

  /* Search for different number of signalling bits. */
  int nb_strength_bits = 0;
  uint64_t best_rd = UINT64_MAX;
  CdefInfo *const cdef_info = &cm->cdef_info;
  for (int i = 0; i <= 3; i++) {
    int best_lev0[CDEF_MAX_STRENGTHS];
    int best_lev1[CDEF_MAX_STRENGTHS] = { 0 };
    const int nb_strengths = 1 << i;
    uint64_t tot_mse;
    if (num_planes > 1) {
      tot_mse = joint_strength_search_dual(best_lev0, best_lev1, nb_strengths,
                                           mse, sb_count, pick_method);

    } else {
      tot_mse = joint_strength_search(best_lev0, nb_strengths, mse[0], sb_count,
                                      pick_method);
    }

    const int total_bits = sb_count * i + nb_strengths * CDEF_STRENGTH_BITS *
                                              (num_planes > 1 ? 2 : 1);

    const int rate_cost = av1_cost_literal(total_bits);
    const uint64_t dist = tot_mse * 16;
    const uint64_t rd = RDCOST(rdmult, rate_cost, dist);
    if (rd < best_rd) {
      best_rd = rd;
      nb_strength_bits = i;
      memcpy(cdef_info->cdef_strengths, best_lev0,
             nb_strengths * sizeof(best_lev0[0]));
      if (num_planes > 1) {
        memcpy(cdef_info->cdef_uv_strengths, best_lev1,
               nb_strengths * sizeof(best_lev1[0]));
      }
    }
  }

  cdef_info->cdef_bits = nb_strength_bits;
  cdef_info->nb_cdef_strengths = 1 << nb_strength_bits;
  for (int i = 0; i < sb_count; i++) {
    uint64_t best_mse = UINT64_MAX;
    int best_gi = 0;
    for (int gi = 0; gi < cdef_info->nb_cdef_strengths; gi++) {
      uint64_t curr = mse[0][i][cdef_info->cdef_strengths[gi]];

      if (num_planes > 1) curr += mse[1][i][cdef_info->cdef_uv_strengths[gi]];

      if (curr < best_mse) {
        best_gi = gi;
        best_mse = curr;
      }
    }
    mi_params->mi_grid_base[sb_index[i]]->cdef_strength = best_gi;
  }

  if (fast) {
    for (int j = 0; j < cdef_info->nb_cdef_strengths; j++) {
      const int luma_strength = cdef_info->cdef_strengths[j];

      const int chroma_strength = cdef_info->cdef_uv_strengths[j];
      int pri_strength, sec_strength;

      STORE_CDEF_FILTER_STRENGTH(cdef_info->cdef_strengths[j], pick_method,
                                 luma_strength);
      STORE_CDEF_FILTER_STRENGTH(cdef_info->cdef_uv_strengths[j], pick_method,
                                 chroma_strength);
    }
  }

  cdef_info->cdef_damping = damping;

  aom_free(mse[0]);
  aom_free(mse[1]);
  aom_free(sb_index);
}

#if CONFIG_CC_CDEF
void av1_cc_cdef_search(const YV12_BUFFER_CONFIG *frame,
                        const YV12_BUFFER_CONFIG *ref, AV1_COMMON *cm,
                        MACROBLOCKD *xd, CDEF_PICK_METHOD pick_method,
                        int rdmult, int key_freq_max, int key_freq_min) {
  const int num_planes = av1_num_planes(cm);
  if (!cm->seq_params.enable_cc_cdef || num_planes == 1) {
    cm->cdef_info.cc_cdef_info.cccdef_frame_enable_flag[0] = 0;
    cm->cdef_info.cc_cdef_info.cccdef_frame_enable_flag[1] = 0;
    memset(cm->cdef_info.cc_cdef_info.filter_coeff, 0,
           MAX_NUMBER_OF_DIRECTIONS * MAX_NUMBER_OF_CCCDEF_FILTER_COEFF *
               sizeof(cm->cdef_info.cc_cdef_info.filter_coeff[0][0]));
    cm->cdef_info.cc_cdef_info.cccdef_bits[0] = 0;
    cm->cdef_info.cc_cdef_info.nb_cccdef_strengths[0] = 1;
    cm->cdef_info.cc_cdef_info.cccdef_bits[1] = 0;
    cm->cdef_info.cc_cdef_info.nb_cccdef_strengths[1] = 1;

    for (int plane = 0; plane < 2; plane++) {
      for (int dir = 0; dir < MAX_NUMBER_OF_DIRECTIONS; dir++) {
        cm->cdef_info.cc_cdef_info
            .cccdef_frame_new_filter_signal_flag[plane][dir] = 0;
        cm->cdef_info.cc_cdef_info.cccdef_frame_filter_idx_in_buf[plane][dir] =
            0;
      }
    }
    return;
  }

  cdef_list dlist[MI_SIZE_128X128 * MI_SIZE_128X128];
  int dir_chroma[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int var[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };
  int dir_luma[CDEF_NBLOCKS][CDEF_NBLOCKS] = { { 0 } };

  CdefInfo *const cdef_info = &cm->cdef_info;

  // set the CC-CDEF frame enable flag initialization
  cdef_info->cc_cdef_info.cccdef_frame_enable_flag[0] = 1;
  cdef_info->cc_cdef_info.cccdef_frame_enable_flag[1] = 1;

  const int num_filter_coeff = MAX_NUMBER_OF_CCCDEF_FILTER_COEFF;
  for (int plane = 0; plane < 2; plane++) {
    for (int dir = 0; dir < MAX_NUMBER_OF_DIRECTIONS; dir++) {
      cdef_info->cc_cdef_info.cccdef_frame_new_filter_signal_flag[plane][dir] =
          0;
      cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[plane][dir] = 0;
    }
  }

  if (cm->cur_frame->frame_type == KEY_FRAME ||
      cm->cur_frame->frame_type == INTRA_ONLY_FRAME) {
    av1_fill_cccdef_filter_coeff_buffer_with_default_filters(&cm->cdef_info);
  }

  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int nvfb = (mi_params->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  const int nhfb = (mi_params->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  int *sb_index = aom_malloc(nvfb * nhfb * sizeof(*sb_index));

  DECLARE_ALIGNED(32, uint16_t, tmp_rec[1 << (MAX_SB_SIZE_LOG2 * 2)]);
  DECLARE_ALIGNED(32, uint16_t, tmp_dst_ccdef[1 << (MAX_SB_SIZE_LOG2 * 2)]);

  static const int conv422[8] = { 7, 0, 2, 4, 5, 6, 6, 6 };
  static const int conv440[8] = { 1, 2, 2, 2, 3, 4, 6, 0 };

  av1_setup_dst_planes(xd->plane, cm->seq_params.sb_size, frame, 0, 0, 0,
                       num_planes);

  cccdefStats(*filter_block_stats[2])[MAX_NUMBER_OF_DIRECTIONS];
  filter_block_stats[0] =
      aom_malloc(sizeof(**filter_block_stats) * nvfb * nhfb);
  filter_block_stats[1] =
      aom_malloc(sizeof(**filter_block_stats) * nvfb * nhfb);

  memset(filter_block_stats[0], 0, sizeof(**filter_block_stats) * nvfb * nhfb);
  memset(filter_block_stats[1], 0, sizeof(**filter_block_stats) * nvfb * nhfb);

  int bsize[3];
  int mi_wide_l2[3];
  int mi_high_l2[3];
  int xdec[3];
  int ydec[3];
  uint8_t *ref_buffer[3] = { ref->y_buffer, ref->u_buffer, ref->v_buffer };
  int ref_stride[3] = { ref->y_stride, ref->uv_stride, ref->uv_stride };

  for (int pli = 0; pli < num_planes; pli++) {
    xdec[pli] = xd->plane[pli].subsampling_x;
    ydec[pli] = xd->plane[pli].subsampling_y;
    bsize[pli] = ydec[pli] ? (xdec[pli] ? BLOCK_4X4 : BLOCK_8X4)
                           : (xdec[pli] ? BLOCK_4X8 : BLOCK_8X8);
    mi_wide_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_x;
    mi_high_l2[pli] = MI_SIZE_LOG2 - xd->plane[pli].subsampling_y;
  }

  copy_fn_t copy_fn;
  compute_cdef_dist_t compute_cdef_dist_fn;
  compute_cdef_dist_direction_t compute_cdef_dist_fn_direction;

  if (cm->seq_params.use_highbitdepth) {
    copy_fn = copy_sb16_16_highbd;
    compute_cdef_dist_fn = compute_cdef_dist_highbd;
    compute_cdef_dist_fn_direction = compute_cdef_dist_highbd_direction;
  } else {
    copy_fn = copy_sb16_16;
    compute_cdef_dist_fn = compute_cdef_dist;
    compute_cdef_dist_fn_direction = compute_cdef_dist_direction;
  }

  DECLARE_ALIGNED(32, uint16_t, inbuf[CDEF_INBUF_SIZE]);
  uint16_t *const in_chroma =
      inbuf + CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER;

  DECLARE_ALIGNED(32, uint16_t, inbuf_luma[CDEF_INBUF_SIZE]);
  uint16_t *const in_luma =
      inbuf_luma + CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER;

  const int coeff_shift = AOMMAX(cm->seq_params.bit_depth - 8, 0);
  int fb_count = 0;
  int non_skip_fb_count = 0;
  for (int fbr = 0; fbr < nvfb; ++fbr) {
    for (int fbc = 0; fbc < nhfb; ++fbc) {
      // No filtering if the entire filter block is skipped
      if (sb_all_skip(mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64)) {
        fb_count++;
        continue;
      }

      const MB_MODE_INFO *const mbmi =
          mi_params->mi_grid_base[MI_SIZE_64X64 * fbr * mi_params->mi_stride +
                                  MI_SIZE_64X64 * fbc];
      if (((fbc & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64)) ||
          ((fbr & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_64X128)))
        continue;

      int nhb = AOMMIN(MI_SIZE_64X64, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
      int nvb = AOMMIN(MI_SIZE_64X64, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
      int hb_step = 1;
      int vb_step = 1;
      BLOCK_SIZE bs;
      if (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64 ||
          mbmi->sb_type == BLOCK_64X128) {
        bs = mbmi->sb_type;
        if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
          nhb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
          hb_step = 2;
        }
        if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
          nvb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
          vb_step = 2;
        }
      } else {
        bs = BLOCK_64X64;
      }

      const int cdef_count = av1_cdef_compute_sb_list(
          mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, bs);

      const int yoff = CDEF_VBORDER * (fbr != 0);
      const int xoff = CDEF_HBORDER * (fbc != 0);

      // Prepare luma input buffer
      // Directions based on the luma input buffer
      /* We avoid filtering the pixels for which some of the pixels to
              average are outside the frame. We could change the filter instead,
              but it would add special cases for any future vectorization. */
      for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf_luma[i] = CDEF_VERY_LARGE;

      const int ysize_luma =
          (nvb << mi_high_l2[0]) + CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
      const int xsize_luma =
          (nhb << mi_wide_l2[0]) + CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
      const int row_luma = fbr * MI_SIZE_64X64 << mi_high_l2[0];
      const int col_luma = fbc * MI_SIZE_64X64 << mi_wide_l2[0];

      copy_fn(&in_luma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
              xd->plane[0].dst.buf, row_luma - yoff, col_luma - xoff,
              xd->plane[0].dst.stride, ysize_luma, xsize_luma);

      // Get the directions for all 8x8 blocks within the 64x64 filter block
      // based on Luma
      for (int bi = 0; bi < cdef_count; bi++) {
        int by = dlist[bi].by;
        int bx = dlist[bi].bx;
        dir_luma[by][bx] =
            cdef_find_dir(&in_luma[8 * by * CDEF_BSTRIDE + 8 * bx],
                          CDEF_BSTRIDE, &var[by][bx], coeff_shift);

        // special care for 422
        if (xdec[0] != ydec[0]) {
          dir_chroma[by][bx] = (xdec[0] ? conv422 : conv440)[dir_luma[by][bx]];
        } else {
          dir_chroma[by][bx] = dir_luma[by][bx];
        }

        dir_luma[by][bx] = cdefdir_to_cccdef_dir[dir_luma[by][bx]];
      }

      for (int pli = 1; pli < num_planes; pli++) {
        for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf[i] = CDEF_VERY_LARGE;
        /* We avoid filtering the pixels for which some of the pixels to
           average are outside the frame. We could change the filter instead,
           but it would add special cases for any future vectorization. */
        const int ysize = (nvb << mi_high_l2[pli]) +
                          CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
        const int xsize = (nhb << mi_wide_l2[pli]) +
                          CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
        const int row = fbr * MI_SIZE_64X64 << mi_high_l2[pli];
        const int col = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];

        // Apply CDEF using best strength and damping
        int best_pri_strength, best_sec_strength;
        int mbmi_cdef_strength = mbmi->cdef_strength;
        int best_cdef_strength =
            cdef_info->cdef_uv_strengths[mbmi_cdef_strength];

        get_cdef_filter_strengths(pick_method, &best_pri_strength,
                                  &best_sec_strength, best_cdef_strength);
        copy_fn(&in_chroma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
                xd->plane[pli].dst.buf, row - yoff, col - xoff,
                xd->plane[pli].dst.stride, ysize, xsize);

        av1_cdef_filter_fb(NULL, tmp_rec, CDEF_BSTRIDE, in_chroma, xdec[pli],
                           ydec[pli], dir_chroma, NULL, var, pli, dlist,
                           cdef_count, best_pri_strength,
                           best_sec_strength + (best_sec_strength == 3),
                           cdef_info->cdef_damping, coeff_shift, 0);

        if (cm->seq_params.use_highbitdepth) {
          uint16_t *dst16 = CONVERT_TO_SHORTPTR((uint8_t *)ref_buffer[pli]);
          uint16_t *org = &dst16[row * ref_stride[pli] + col];

          av1_get_filter_statstics_64x64(
              tmp_rec, CDEF_BSTRIDE, in_luma, CDEF_BSTRIDE, NULL, org,
              ref_stride[pli], xdec[pli], ydec[pli], dir_luma, dlist,
              cdef_count, 0, filter_block_stats[pli - 1][fb_count]);
        } else {
          uint8_t *dst8 = (uint8_t *)ref_buffer[pli];
          uint8_t *org = &dst8[row * ref_stride[pli] + col];

          av1_get_filter_statstics_64x64(
              tmp_rec, CDEF_BSTRIDE, in_luma, CDEF_BSTRIDE, org, NULL,
              ref_stride[pli], xdec[pli], ydec[pli], dir_luma, dlist,
              cdef_count, 0, filter_block_stats[pli - 1][fb_count]);
        }

      }  // end for for pli

      fb_count++;
      sb_index[non_skip_fb_count++] =
          MI_SIZE_64X64 * fbr * mi_params->mi_stride + MI_SIZE_64X64 * fbc;
    }
  }

  // aggrate to frame stats
  cccdefStats frame_stats[2][MAX_NUMBER_OF_DIRECTIONS];  // plane, directions

  // process Cb
  get_frame_filter_stats(cm, filter_block_stats[0], frame_stats[0], 0);
  derive_filter_coeffs_from_stats(frame_stats[0],
                                  cdef_info->cc_cdef_info.filter_coeff[0], 0);

  // Process Cr
  get_frame_filter_stats(cm, filter_block_stats[1], frame_stats[1], 0);
  derive_filter_coeffs_from_stats(frame_stats[1],
                                  cdef_info->cc_cdef_info.filter_coeff[1], 0);

  // filter coefficients are selected
  // Now search filter strengths
  const int total_strengths = TOTAL_CCCDEF_STRENGTHS;
  const int number_of_directions = MAX_NUMBER_OF_DIRECTIONS;
  uint64_t(*mse[2])[TOTAL_CCCDEF_STRENGTHS];
  uint64_t frame_distortion[MAX_NUMBER_OF_DIRECTIONS]
                           [2 * (MAX_NUMBER_OF_TEMPORAL_FILTERS + 1)];
  int number_of_filter_sets[2][MAX_NUMBER_OF_DIRECTIONS];
  uint64_t unfilter_direction_distortion[2][MAX_NUMBER_OF_DIRECTIONS];

  for (int dir = 0; dir < number_of_directions; dir++) {
    CC_CdefFilterBuf *coeff_buffer = &cm->cdef_info.cccdef_filter_buf[dir];
    number_of_filter_sets[0][dir] =
        (1 + coeff_buffer->curr_number_of_filter_sets[0]);
    number_of_filter_sets[1][dir] =
        (1 + coeff_buffer->curr_number_of_filter_sets[1]);
    for (int filter_set = 0; filter_set < (number_of_filter_sets[0][dir] +
                                           number_of_filter_sets[1][dir]);
         filter_set++) {
      frame_distortion[dir][filter_set] = 0;
    }
    unfilter_direction_distortion[0][dir] = 0;
    unfilter_direction_distortion[1][dir] = 0;
  }

  for (int uv = 0; uv < 2; uv++) {
    mse[uv] = aom_malloc(sizeof(**mse) * nvfb * nhfb);
  }

  uint64_t unfilter_frame_distortion[2] = { 0, 0 };
  uint64_t filtered_frame_distortion[2] = { 0, 0 };

  uint64_t unfilter_frame_bits[2] = { 0, 0 };
  uint64_t filtered_frame_bits[2] = { 0, 0 };

  fb_count = 0;
  non_skip_fb_count = 0;
  for (int fbr = 0; fbr < nvfb; ++fbr) {
    for (int fbc = 0; fbc < nhfb; ++fbc) {
      // No filtering if the entire filter block is skipped
      if (sb_all_skip(mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64)) {
        fb_count++;
        continue;
      }

      const MB_MODE_INFO *const mbmi =
          mi_params->mi_grid_base[MI_SIZE_64X64 * fbr * mi_params->mi_stride +
                                  MI_SIZE_64X64 * fbc];
      if (((fbc & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64)) ||
          ((fbr & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_64X128)))
        continue;

      int nhb = AOMMIN(MI_SIZE_64X64, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
      int nvb = AOMMIN(MI_SIZE_64X64, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
      int hb_step = 1;
      int vb_step = 1;
      BLOCK_SIZE bs;
      if (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64 ||
          mbmi->sb_type == BLOCK_64X128) {
        bs = mbmi->sb_type;
        if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
          nhb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
          hb_step = 2;
        }
        if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
          nvb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
          vb_step = 2;
        }
      } else {
        bs = BLOCK_64X64;
      }

      const int cdef_count = av1_cdef_compute_sb_list(
          mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, bs);

      const int yoff = CDEF_VBORDER * (fbr != 0);
      const int xoff = CDEF_HBORDER * (fbc != 0);

      // Prepare luma input buffer
      // Directions based on the luma input buffer
      /* We avoid filtering the pixels for which some of the pixels to
              average are outside the frame. We could change the filter instead,
              but it would add special cases for any future vectorization. */
      int dirinit = 1;
      for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf_luma[i] = CDEF_VERY_LARGE;

      const int ysize_luma =
          (nvb << mi_high_l2[0]) + CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
      const int xsize_luma =
          (nhb << mi_wide_l2[0]) + CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
      const int row_luma = fbr * MI_SIZE_64X64 << mi_high_l2[0];
      const int col_luma = fbc * MI_SIZE_64X64 << mi_wide_l2[0];

      copy_fn(&in_luma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
              xd->plane[0].dst.buf, row_luma - yoff, col_luma - xoff,
              xd->plane[0].dst.stride, ysize_luma, xsize_luma);

      // Get the directions for all 8x8 blocks within the 64x64 filter block
      // based on Luma
      for (int bi = 0; bi < cdef_count; bi++) {
        int by = dlist[bi].by;
        int bx = dlist[bi].bx;
        dir_luma[by][bx] =
            cdef_find_dir(&in_luma[8 * by * CDEF_BSTRIDE + 8 * bx],
                          CDEF_BSTRIDE, &var[by][bx], coeff_shift);

        // special care for 422
        if (xdec[0] != ydec[0]) {
          dir_chroma[by][bx] = (xdec[0] ? conv422 : conv440)[dir_luma[by][bx]];
        } else {
          dir_chroma[by][bx] = dir_luma[by][bx];
        }

        dir_luma[by][bx] = cdefdir_to_cccdef_dir[dir_luma[by][bx]];
      }

      for (int pli = 1; pli < num_planes; pli++) {
        for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf[i] = CDEF_VERY_LARGE;
        /* We avoid filtering the pixels for which some of the pixels to
           average are outside the frame. We could change the filter instead,
           but it would add special cases for any future vectorization. */
        const int ysize = (nvb << mi_high_l2[pli]) +
                          CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
        const int xsize = (nhb << mi_wide_l2[pli]) +
                          CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
        const int row = fbr * MI_SIZE_64X64 << mi_high_l2[pli];
        const int col = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];

        // Apply CDEF using best strength and damping
        int best_pri_strength, best_sec_strength;
        int mbmi_cdef_strength = mbmi->cdef_strength;
        int best_cdef_strength =
            cdef_info->cdef_uv_strengths[mbmi_cdef_strength];

        get_cdef_filter_strengths(pick_method, &best_pri_strength,
                                  &best_sec_strength, best_cdef_strength);
        copy_fn(&in_chroma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
                xd->plane[pli].dst.buf, row - yoff, col - xoff,
                xd->plane[pli].dst.stride, ysize, xsize);

        av1_cdef_filter_fb(NULL, tmp_rec, CDEF_BSTRIDE, in_chroma, xdec[pli],
                           ydec[pli], dir_chroma, &dirinit, var, pli, dlist,
                           cdef_count, best_pri_strength,
                           best_sec_strength + (best_sec_strength == 3),
                           cdef_info->cdef_damping, coeff_shift, 0);

        const int bw_log2 = 3 - xdec[pli];
        const int bh_log2 = 3 - ydec[pli];
        size_t num_dst_pixels = (1 << bh_log2) * (1 << bw_log2) * cdef_count;

        for (int dir = 0; dir < number_of_directions; dir++) {
          CC_CdefFilterBuf *coeff_buffer =
              &cm->cdef_info.cccdef_filter_buf[dir];
          uint8_t direction_enable_flags[MAX_NUMBER_OF_DIRECTIONS] = { 0, 0, 0,
                                                                       0 };
          direction_enable_flags[dir] = 1;
          unfilter_direction_distortion[pli - 1][dir] +=
              compute_cdef_dist_fn_direction(
                  ref_buffer[pli], ref_stride[pli], tmp_rec, dlist, cdef_count,
                  bsize[pli], coeff_shift, row, col, dir, dir_luma);

          for (int filter_set = 0;
               filter_set < (number_of_filter_sets[pli - 1][dir]);
               filter_set++) {
            const short *const filter_coeff =
                (filter_set == (number_of_filter_sets[pli - 1][dir] - 1))
                    ? cdef_info->cc_cdef_info.filter_coeff[pli - 1]
                    : (&(coeff_buffer
                             ->buf[pli - 1]
                                  [filter_set *
                                   coeff_buffer
                                       ->number_of_coeff_ineach_set[pli - 1]]));

            const int single_direction_coeff =
                (filter_set == (number_of_filter_sets[pli - 1][dir] - 1)) ? 0
                                                                          : 1;

            const int uvoffset =
                (pli == 1) ? filter_set
                           : (number_of_filter_sets[0][dir] + filter_set);

            // Search strengths for CC-CDEF
            int gi = total_strengths - 1;

            // copy the output of the CDEF in the temporary buffer
            memcpy(tmp_dst_ccdef, tmp_rec, num_dst_pixels * sizeof(tmp_rec[0]));

            av1_cccdef_filter_fb(NULL, tmp_dst_ccdef, CDEF_BSTRIDE, in_luma,
                                 xdec[pli], ydec[pli], dir_luma, dlist,
                                 cdef_count, cm->seq_params.bit_depth, 1,
                                 filter_coeff, gi, 0, direction_enable_flags,
                                 single_direction_coeff);

            const uint64_t curr_mse = compute_cdef_dist_fn_direction(
                ref_buffer[pli], ref_stride[pli], tmp_dst_ccdef, dlist,
                cdef_count, bsize[pli], coeff_shift, row, col, dir, dir_luma);

            frame_distortion[dir][uvoffset] += curr_mse;

          }  // end of filter sets

        }  // end of direction

      }  // end for for pli

      fb_count++;
      sb_index[non_skip_fb_count++] =
          MI_SIZE_64X64 * fbr * mi_params->mi_stride + MI_SIZE_64X64 * fbc;
    }
  }

  /* Search filter coefficients for each direction */
  for (int uv = 0; uv < 2; uv++) {
    for (int dir = 0; dir < number_of_directions; dir++) {
      CC_CdefFilterBuf *coeff_buffer = &cm->cdef_info.cccdef_filter_buf[dir];

      int best_filter_set_idx = number_of_filter_sets[uv][dir] -
                                1;  // initially set the new filer as best
      uint64_t best_rd_coeff = UINT64_MAX;

      // filter_set = 0 ~ (number_of_filter_sets[uv][dir] - 2) Filters are in
      // history buffer
      // filter_set == (number_of_filter_sets[uv][dir] - 1)  New filter
      // filter_set == number_of_filter_sets[uv][dir] Filter OFF

      for (int filter_set = 0; filter_set <= (number_of_filter_sets[uv][dir]);
           filter_set++) {
        const int uvoff = (uv == 0)
                              ? filter_set
                              : (number_of_filter_sets[0][dir] + filter_set);

        int frame_bits = 0;
        frame_bits++;  // 1 bit to signal direction level ON/OFF

        if (filter_set != number_of_filter_sets[uv][dir]) {
          frame_bits++;  // 1 bit to signal signal new filter flag
          if (filter_set == (number_of_filter_sets[uv][dir] - 1)) {
            const short *const filter_coeff =
                &cdef_info->cc_cdef_info
                     .filter_coeff[uv][dir * num_filter_coeff];
            frame_bits += get_ccdef_coeff_total_bit_length(filter_coeff, 1);
          } else {
            frame_bits += CC_CDEF_FILTER_SET_IDX_BITS;
          }
        }

        uint64_t dist = (filter_set == number_of_filter_sets[uv][dir])
                            ? unfilter_direction_distortion[uv][dir]
                            : frame_distortion[dir][uvoff];

        dist = dist * 16;

        const int rate_cost = av1_cost_literal(frame_bits);
        const uint64_t rd = RDCOST(rdmult, rate_cost, dist);
        if (rd < best_rd_coeff) {
          best_rd_coeff = rd;
          best_filter_set_idx = filter_set;
        }
      }  // end of filter sets

      int force_new_filter = (coeff_buffer->number_of_new_filter_sets[uv] == 0);

      short *filter_coeff_dir =
          &cdef_info->cc_cdef_info.filter_coeff[uv][dir * num_filter_coeff];
      int at_least_one_coeff_nonzero = 0;
      for (int q = 0; q < num_filter_coeff; q++) {
        at_least_one_coeff_nonzero |= (filter_coeff_dir[q] != 0);
      }

      force_new_filter &= at_least_one_coeff_nonzero;

      force_new_filter &=
          (key_freq_min > 1 && key_freq_max > 1);  // Forcely signal new filter

      if (force_new_filter)
        best_filter_set_idx = (number_of_filter_sets[uv][dir] - 1);

      if (best_filter_set_idx == number_of_filter_sets[uv][dir]) {
        // Filter is off for this direction
        cdef_info->cc_cdef_info.cccdef_frame_direction_enable_flag[uv][dir] = 0;
        cdef_info->cc_cdef_info.cccdef_frame_new_filter_signal_flag[uv][dir] =
            0;
        cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[uv][dir] = 0;

        // Reset the filter coeff to 0
        for (int coeff = 0; coeff < (num_filter_coeff); coeff++)
          cdef_info->cc_cdef_info
              .filter_coeff[uv][dir * num_filter_coeff + coeff] = 0;
      } else if (best_filter_set_idx == (number_of_filter_sets[uv][dir] - 1)) {
        // Filter is ON for this direction
        // signal new filter coefficients
        cdef_info->cc_cdef_info.cccdef_frame_direction_enable_flag[uv][dir] = 1;
        cdef_info->cc_cdef_info.cccdef_frame_new_filter_signal_flag[uv][dir] =
            1;
        cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[uv][dir] = 0;
      } else {
        // Filter is ON for this direction
        // History filter is used
        cdef_info->cc_cdef_info.cccdef_frame_direction_enable_flag[uv][dir] = 1;
        cdef_info->cc_cdef_info.cccdef_frame_new_filter_signal_flag[uv][dir] =
            0;
        cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[uv][dir] =
            best_filter_set_idx;
        short *dst =
            &cdef_info->cc_cdef_info.filter_coeff[uv][dir * num_filter_coeff];

        copy_cccdef_filter_coeff_from_buffer(
            cm, cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[uv][dir],
            dst, uv, dir);  // copy the filter coeff
      }

    }  // for (int dir = 0; dir < number_of_directions; dir++)

  }  // end of uv

  // Filter coefficients are already derived
  // strength search start from here
  // search best strengths parameters

  unfilter_frame_distortion[0] = 0;
  unfilter_frame_distortion[1] = 0;

  fb_count = 0;
  non_skip_fb_count = 0;
  for (int fbr = 0; fbr < nvfb; ++fbr) {
    for (int fbc = 0; fbc < nhfb; ++fbc) {
      // No filtering if the entire filter block is skipped
      if (sb_all_skip(mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64)) {
        fb_count++;
        continue;
      }

      const MB_MODE_INFO *const mbmi =
          mi_params->mi_grid_base[MI_SIZE_64X64 * fbr * mi_params->mi_stride +
                                  MI_SIZE_64X64 * fbc];
      if (((fbc & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64)) ||
          ((fbr & 1) &&
           (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_64X128)))
        continue;

      int nhb = AOMMIN(MI_SIZE_64X64, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
      int nvb = AOMMIN(MI_SIZE_64X64, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
      int hb_step = 1;
      int vb_step = 1;
      BLOCK_SIZE bs;
      if (mbmi->sb_type == BLOCK_128X128 || mbmi->sb_type == BLOCK_128X64 ||
          mbmi->sb_type == BLOCK_64X128) {
        bs = mbmi->sb_type;
        if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
          nhb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_cols - MI_SIZE_64X64 * fbc);
          hb_step = 2;
        }
        if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
          nvb =
              AOMMIN(MI_SIZE_128X128, mi_params->mi_rows - MI_SIZE_64X64 * fbr);
          vb_step = 2;
        }
      } else {
        bs = BLOCK_64X64;
      }

      const int cdef_count = av1_cdef_compute_sb_list(
          mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, bs);

      const int yoff = CDEF_VBORDER * (fbr != 0);
      const int xoff = CDEF_HBORDER * (fbc != 0);

      // Prepare luma input buffer
      // Directions based on the luma input buffer
      /* We avoid filtering the pixels for which some of the pixels to
                      average are outside the frame. We could change the filter
         instead, but it would add special cases for any future vectorization.
       */
      int dirinit = 1;
      for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf_luma[i] = CDEF_VERY_LARGE;

      const int ysize_luma =
          (nvb << mi_high_l2[0]) + CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
      const int xsize_luma =
          (nhb << mi_wide_l2[0]) + CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
      const int row_luma = fbr * MI_SIZE_64X64 << mi_high_l2[0];
      const int col_luma = fbc * MI_SIZE_64X64 << mi_wide_l2[0];

      copy_fn(&in_luma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
              xd->plane[0].dst.buf, row_luma - yoff, col_luma - xoff,
              xd->plane[0].dst.stride, ysize_luma, xsize_luma);

      // Get the directions for all 8x8 blocks within the 64x64 filter block
      // based on Luma
      for (int bi = 0; bi < cdef_count; bi++) {
        int by = dlist[bi].by;
        int bx = dlist[bi].bx;
        dir_luma[by][bx] =
            cdef_find_dir(&in_luma[8 * by * CDEF_BSTRIDE + 8 * bx],
                          CDEF_BSTRIDE, &var[by][bx], coeff_shift);

        // special care for 422
        if (xdec[0] != ydec[0]) {
          dir_chroma[by][bx] = (xdec[0] ? conv422 : conv440)[dir_luma[by][bx]];
        } else {
          dir_chroma[by][bx] = dir_luma[by][bx];
        }

        dir_luma[by][bx] = cdefdir_to_cccdef_dir[dir_luma[by][bx]];
      }

      for (int pli = 1; pli < num_planes; pli++) {
        for (int i = 0; i < CDEF_INBUF_SIZE; i++) inbuf[i] = CDEF_VERY_LARGE;
        /* We avoid filtering the pixels for which some of the pixels to
         average are outside the frame. We could change the filter instead,
        but it would add special cases for any future vectorization. */
        const int ysize = (nvb << mi_high_l2[pli]) +
                          CDEF_VBORDER * (fbr + vb_step < nvfb) + yoff;
        const int xsize = (nhb << mi_wide_l2[pli]) +
                          CDEF_HBORDER * (fbc + hb_step < nhfb) + xoff;
        const int row = fbr * MI_SIZE_64X64 << mi_high_l2[pli];
        const int col = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];

        // Apply CDEF using best strength and damping
        int best_pri_strength, best_sec_strength;
        int mbmi_cdef_strength = mbmi->cdef_strength;
        int best_cdef_strength =
            cdef_info->cdef_uv_strengths[mbmi_cdef_strength];

        get_cdef_filter_strengths(pick_method, &best_pri_strength,
                                  &best_sec_strength, best_cdef_strength);
        copy_fn(&in_chroma[(-yoff * CDEF_BSTRIDE - xoff)], CDEF_BSTRIDE,
                xd->plane[pli].dst.buf, row - yoff, col - xoff,
                xd->plane[pli].dst.stride, ysize, xsize);

        av1_cdef_filter_fb(NULL, tmp_rec, CDEF_BSTRIDE, in_chroma, xdec[pli],
                           ydec[pli], dir_chroma, &dirinit, var, pli, dlist,
                           cdef_count, best_pri_strength,
                           best_sec_strength + (best_sec_strength == 3),
                           cdef_info->cdef_damping, coeff_shift, 0);

        const uint64_t curr_mse_after_cdef = compute_cdef_dist_fn(
            ref_buffer[pli], ref_stride[pli], tmp_rec, dlist, cdef_count,
            bsize[pli], coeff_shift, row, col);

        unfilter_frame_distortion[pli - 1] += curr_mse_after_cdef;

        const int bw_log2 = 3 - xdec[pli];
        const int bh_log2 = 3 - ydec[pli];
        size_t num_dst_pixels = (1 << bh_log2) * (1 << bw_log2) * cdef_count;

        const int direction_merge_mode = 0;
        const short *const filter_coeff =
            cdef_info->cc_cdef_info.filter_coeff[pli - 1];
        const uint8_t *const direction_enable_flags =
            &cdef_info->cc_cdef_info
                 .cccdef_frame_direction_enable_flag[pli - 1][0];

        const int uvoffset = (pli - 1);

        // Search strengths for CC-CDEF
        for (int gi = 0; gi < total_strengths; gi++) {
          // copy the output of the CDEF in the temporary buffer
          memcpy(tmp_dst_ccdef, tmp_rec, num_dst_pixels * sizeof(tmp_rec[0]));

          av1_cccdef_filter_fb(NULL, tmp_dst_ccdef, CDEF_BSTRIDE, in_luma,
                               xdec[pli], ydec[pli], dir_luma, dlist,
                               cdef_count, cm->seq_params.bit_depth, 1,
                               filter_coeff, gi, direction_merge_mode,
                               direction_enable_flags, 0);

          const uint64_t curr_mse = compute_cdef_dist_fn(
              ref_buffer[pli], ref_stride[pli], tmp_dst_ccdef, dlist,
              cdef_count, bsize[pli], coeff_shift, row, col);

          mse[uvoffset][non_skip_fb_count][gi] = curr_mse;
        }

      }  // end for for pli

      fb_count++;
      sb_index[non_skip_fb_count++] =
          MI_SIZE_64X64 * fbr * mi_params->mi_stride + MI_SIZE_64X64 * fbc;
    }
  }

  for (int uv = 0; uv < 2; uv++) {
    int nb_strength_bits = 0;
    uint64_t best_rd = UINT64_MAX;
    const int uvoffset = uv;
    int at_least_one_direction_on = 0;
    for (int dir = 0; dir < number_of_directions; dir++) {
      at_least_one_direction_on |=
          cdef_info->cc_cdef_info.cccdef_frame_direction_enable_flag[uv][dir];
    }

    if (at_least_one_direction_on) {
      // Search the strengths parameters
      for (int i = 0; i <= 2; i++) {
        int best_lev[TOTAL_CCCDEF_STRENGTHS];
        const int nb_strengths = 1 << i;
        uint64_t tot_mse;

        tot_mse = joint_strength_search_ccdef(best_lev, nb_strengths,
                                              mse[uvoffset], non_skip_fb_count);

        const int total_bits =
            non_skip_fb_count * i + nb_strengths * CC_CDEF_STRENGTH_VALUE_BITS;

        const int rate_cost = av1_cost_literal(total_bits);
        const uint64_t dist = tot_mse * 16;
        const uint64_t rd = RDCOST(rdmult, rate_cost, dist);
        if (rd < best_rd) {
          best_rd = rd;
          nb_strength_bits = i;

          memcpy(cdef_info->cc_cdef_info.cccdef_strengths_index_array_frame[uv],
                 best_lev, nb_strengths * sizeof(best_lev[0]));
        }
      }

      // May need to remove the duplicates entries in here
      cdef_info->cc_cdef_info.cccdef_bits[uv] = nb_strength_bits;
      cdef_info->cc_cdef_info.nb_cccdef_strengths[uv] =
          (1 << cdef_info->cc_cdef_info.cccdef_bits[uv]);

      filtered_frame_distortion[uv] = 0;

      // Find the index of the filter for each 64x64 filter block
      for (int fbindex = 0; fbindex < non_skip_fb_count; fbindex++) {
        uint64_t best_mse = UINT64_MAX;
        int best_strength_index = 0;
        for (int strength_index = 0;
             strength_index < cdef_info->cc_cdef_info.nb_cccdef_strengths[uv];
             strength_index++) {
          uint64_t curr =
              mse[uvoffset][fbindex]
                 [cdef_info->cc_cdef_info
                      .cccdef_strengths_index_array_frame[uv][strength_index]];

          if (curr < best_mse) {
            best_strength_index = strength_index;
            best_mse = curr;
          }
        }

#if FIX_FB_SIGNALING
        const int mbmi_cdef_strength =
            mi_params->mi_grid_base[sb_index[fbindex]]->cdef_strength;
        int level = cm->cdef_info.cdef_strengths[mbmi_cdef_strength] / 4;
        int sec_strength = cm->cdef_info.cdef_strengths[mbmi_cdef_strength] % 4;
        sec_strength += sec_strength == 3;

        int uv_level = cm->cdef_info.cdef_uv_strengths[mbmi_cdef_strength] / 4;
        int uv_sec_strength =
            cm->cdef_info.cdef_uv_strengths[mbmi_cdef_strength] % 4;
        uv_sec_strength += uv_sec_strength == 3;

        if (level == 0 && sec_strength == 0 && uv_level == 0 &&
            uv_sec_strength == 0) {
          best_strength_index = 0;
        }
#endif
        mi_params->mi_grid_base[sb_index[fbindex]]
            ->cc_cdef_strength_index_fb[uv] = best_strength_index;

        filtered_frame_distortion[uv] +=
            mse[uvoffset][fbindex]
               [cdef_info->cc_cdef_info.cccdef_strengths_index_array_frame
                    [uv][best_strength_index]];
      }

      // computing the number of bits for filter ON/OFF
      unfilter_frame_bits[uv]++;  // 1 bit to signal frame level ON/OFF
      filtered_frame_bits[uv]++;  // 1 bit to signal frame level ON/OFF

      for (int dir = 0; dir < number_of_directions; dir++) {
        filtered_frame_bits[uv]++;  // 1 bit to signal direction level ON/OFF

        if (cdef_info->cc_cdef_info
                .cccdef_frame_direction_enable_flag[uv][dir]) {
          filtered_frame_bits[uv]++;  // 1 bit to signal signal new filter flag

          if (cdef_info->cc_cdef_info
                  .cccdef_frame_new_filter_signal_flag[uv][dir]) {
            const short *const filter_coeff =
                &cdef_info->cc_cdef_info
                     .filter_coeff[uv][dir * num_filter_coeff];
            filtered_frame_bits[uv] +=
                get_ccdef_coeff_total_bit_length(filter_coeff, 1);
          } else {
            filtered_frame_bits[uv] += CC_CDEF_FILTER_SET_IDX_BITS;
          }
        }
      }

      // bits to signal the index of the strengths in FB
      filtered_frame_bits[uv] += 2;  // 2 bits to signal cccdef_bits
      int strength_bits = cm->cdef_info.cc_cdef_info.nb_cccdef_strengths[uv] *
                          CC_CDEF_STRENGTH_VALUE_BITS;
      filtered_frame_bits[uv] += strength_bits;
      filtered_frame_bits[uv] +=
          (non_skip_fb_count * cdef_info->cc_cdef_info.cccdef_bits[uv]);

      const uint64_t unfiltered_rd = RDCOST(rdmult, unfilter_frame_bits[uv],
                                            unfilter_frame_distortion[uv] * 16);
      const uint64_t filtered_rd = RDCOST(rdmult, filtered_frame_bits[uv],
                                          filtered_frame_distortion[uv] * 16);

      // printf(" plane = %5d, unfilter dist = %5d filter dist = %5d unfilter
      // bits = %5d filter bits = %5d unfiltered_rd = %5d filtered_rd = %5d \n",
      // uv, unfilter_frame_distortion[uv], filtered_frame_distortion[uv],
      // unfilter_frame_bits[uv], filtered_frame_bits[uv], unfiltered_rd,
      // filtered_rd);

      if (filtered_rd < unfiltered_rd)
        cdef_info->cc_cdef_info.cccdef_frame_enable_flag[uv] = 1;
      else
        cdef_info->cc_cdef_info.cccdef_frame_enable_flag[uv] = 0;

    } else {
      cdef_info->cc_cdef_info.cccdef_frame_enable_flag[uv] = 0;
    }

    // reset all if flag == 0
    if (!cdef_info->cc_cdef_info.cccdef_frame_enable_flag[uv]) {
      cdef_info->cc_cdef_info.cccdef_bits[uv] = 0;
      cdef_info->cc_cdef_info.nb_cccdef_strengths[uv] = 1;

      for (int i = 0; i < cdef_info->cc_cdef_info.nb_cccdef_strengths[uv];
           i++) {
        cdef_info->cc_cdef_info.cccdef_strengths_index_array_frame[uv][i] = 0;
      }

      for (int i = 0; i < non_skip_fb_count; i++) {
        mi_params->mi_grid_base[sb_index[i]]->cc_cdef_strength_index_fb[uv] = 0;
      }

      for (int dir = 0; dir < number_of_directions; dir++) {
        cdef_info->cc_cdef_info.cccdef_frame_direction_enable_flag[uv][dir] = 0;
        cdef_info->cc_cdef_info.cccdef_frame_new_filter_signal_flag[uv][dir] =
            0;
        cdef_info->cc_cdef_info.cccdef_frame_filter_idx_in_buf[uv][dir] = 0;
      }
    }

  }  // end of uv

  aom_free(filter_block_stats[0]);
  aom_free(filter_block_stats[1]);
  for (int uv = 0; uv < 2; uv++) aom_free(mse[uv]);

  aom_free(sb_index);
}

void av1_get_filter_statstics_64x64(
    uint16_t *dst16, int dstride, uint16_t *in, int instride, uint8_t *org8,
    uint16_t *org16, int orgstride, int xdec, int ydec,
    int dir[CDEF_NBLOCKS][CDEF_NBLOCKS], cdef_list *dlist, int cdef_count,
    bool is_rdo, cccdefStats filter_block_stats[MAX_NUMBER_OF_DIRECTIONS]) {
  int bi;
  int bx;
  int by;

  const int bw_log2 = 3 - xdec;
  const int bh_log2 = 3 - ydec;
  const int bw_log2_y = 3;
  const int bh_log2_y = 3;

  const int bsize =
      ydec ? (xdec ? BLOCK_4X4 : BLOCK_8X4) : (xdec ? BLOCK_4X8 : BLOCK_8X8);
  for (bi = 0; bi < cdef_count; bi++) {
    by = dlist[bi].by;
    bx = dlist[bi].bx;

    CHECK(dir[by][bx] >= MAX_NUMBER_OF_DIRECTIONS,
          " The direction exceed the maximum limit ");

    if (org8) {
      cccdef_filter_statistics_8x8_block(
          NULL,
          &dst16[is_rdo ? bi << (bw_log2 + bh_log2)
                        : (by << bh_log2) * dstride + (bx << bw_log2)],
          is_rdo ? 1 << bw_log2 : dstride,
          &in[(by * instride << bh_log2_y) + (bx << bw_log2_y)], instride,
          &org8[(by << bh_log2) * orgstride + (bx << bw_log2)], NULL, orgstride,
          dir[by][bx], bsize, xdec, ydec, filter_block_stats);
    } else {
      cccdef_filter_statistics_8x8_block(
          NULL,
          &dst16[is_rdo ? bi << (bw_log2 + bh_log2)
                        : (by << bh_log2) * dstride + (bx << bw_log2)],
          is_rdo ? 1 << bw_log2 : dstride,
          &in[(by * instride << bh_log2_y) + (bx << bw_log2_y)], instride, NULL,
          &org16[(by << bh_log2) * orgstride + (bx << bw_log2)], orgstride,
          dir[by][bx], bsize, xdec, ydec, filter_block_stats);
    }
  }
}

void cccdef_filter_statistics_8x8_block(
    uint8_t *dst8, uint16_t *dst16, int dstride, const uint16_t *in,
    const int instride, uint8_t *org8, uint16_t *org16, int orgstride, int dir,
    int bsize, int xdec, int ydec,
    cccdefStats filter_block_stats[MAX_NUMBER_OF_DIRECTIONS]) {
  int p, q, i, j;
  double diffs[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];
  int distortion;
  int number_coeff = MAX_NUMBER_OF_CCCDEF_FILTER_COEFF;
  for (p = 0; p < 4 << (bsize == BLOCK_8X8 || bsize == BLOCK_4X8); p++) {
    for (q = 0; q < 4 << (bsize == BLOCK_8X8 || bsize == BLOCK_8X4); q++) {
      int curr_unfiltered_chroma =
          (dst8) ? dst8[p * dstride + q] : dst16[p * dstride + q];
      int curr_original_chroma =
          (org8) ? org8[p * orgstride + q] : org16[p * orgstride + q];
      distortion = curr_original_chroma - curr_unfiltered_chroma;

      // i, j are the position of the corresponding luma position
      i = p << ydec;
      j = q << xdec;
      int coll_luma_pos = i * instride + j;

      int curr_luma = in[coll_luma_pos - 0];

      diffs[0] = (double)in[coll_luma_pos - 2 * instride];
      diffs[1] = (double)in[coll_luma_pos - instride];
      diffs[2] = (double)in[coll_luma_pos - 2];
      diffs[3] = (double)in[coll_luma_pos - 1];
      diffs[4] = (double)in[coll_luma_pos + 1];
      diffs[5] = (double)in[coll_luma_pos + 2];
      diffs[6] = (double)in[coll_luma_pos + instride];
      diffs[7] = (double)in[coll_luma_pos + 2 * instride];

      for (int k = 0; k < number_coeff; k++) {
        diffs[k] = (diffs[k] == CDEF_VERY_LARGE) ? curr_luma : diffs[k];
        diffs[k] -= curr_luma;
      }

      // cross-correlation
      for (int k = 0; k < number_coeff; k++) {
        filter_block_stats[dir].cross_correlation[k] +=
            diffs[k] * (double)distortion;
      }

      // auto-correlation
      for (int k = 0; k < number_coeff; k++) {
        for (int l = k; l < number_coeff; l++) {
          filter_block_stats[dir].auto_correlation[k][l] += diffs[k] * diffs[l];
        }
      }

      for (int k = 1; k < number_coeff; k++) {
        for (int l = 0; l < k; l++) {
          filter_block_stats[dir].auto_correlation[k][l] =
              filter_block_stats[dir].auto_correlation[l][k];
        }
      }
    }
  }
}

void get_frame_filter_stats(
    AV1_COMMON *cm, cccdefStats filter_block_stats[][MAX_NUMBER_OF_DIRECTIONS],
    cccdefStats frame_stats[MAX_NUMBER_OF_DIRECTIONS], int direction_merge) {
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int nvfb = (mi_params->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  const int nhfb = (mi_params->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
  const int numberDirections = MAX_NUMBER_OF_DIRECTIONS;
  const int number_of_coeff = MAX_NUMBER_OF_CCCDEF_FILTER_COEFF;

  // reset the frame stats
  for (int dir = 0; dir < numberDirections; dir++) {
    for (int k = 0; k < number_of_coeff; k++) {
      for (int l = 0; l < number_of_coeff; l++) {
        frame_stats[dir].auto_correlation[k][l] = 0;
      }
    }

    for (int k = 0; k < number_of_coeff; k++) {
      frame_stats[dir].cross_correlation[k] = 0;
    }
  }

  int fb_count = 0;
  for (int fbr = 0; fbr < nvfb; ++fbr) {
    for (int fbc = 0; fbc < nhfb; ++fbc) {
      // No filtering if the entire filter block is skipped
      if (sb_all_skip(mi_params, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64)) {
        fb_count++;
        continue;
      }
      for (int dir = 0; dir < numberDirections; dir++) {
        for (int k = 0; k < number_of_coeff; k++) {
          for (int l = 0; l < number_of_coeff; l++) {
            frame_stats[direction_merge ? 0 : dir].auto_correlation[k][l] +=
                filter_block_stats[fb_count][dir].auto_correlation[k][l];
          }
        }

        for (int k = 0; k < number_of_coeff; k++) {
          frame_stats[direction_merge ? 0 : dir].cross_correlation[k] +=
              filter_block_stats[fb_count][dir].cross_correlation[k];
        }
      }
      fb_count++;
    }
  }
}

void derive_filter_coeffs_from_stats(
    cccdefStats frame_stats[MAX_NUMBER_OF_DIRECTIONS],
    short f[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF * MAX_NUMBER_OF_DIRECTIONS],
    int direction_merge) {
  const int numberDirections = direction_merge ? 1 : MAX_NUMBER_OF_DIRECTIONS;
  const int number_of_coeff = MAX_NUMBER_OF_CCCDEF_FILTER_COEFF;
  const int cccdef_scale_shift = CC_CDEF_SCALE_SHIFT;
  const int scale = (1 << cccdef_scale_shift);
  const int max_coefficient_value = (1 << CC_CDEF_FILTER_COEFF_BITS) - 1;

  double filter_coefficients[MAX_NUMBER_OF_DIRECTIONS]
                            [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];

  for (int dir = 0; dir < numberDirections; dir++) {
    derive_ccdef_filter_coefficients(frame_stats[dir].auto_correlation,
                                     frame_stats[dir].cross_correlation,
                                     filter_coefficients[dir], number_of_coeff);
    // Derive integer filter coefficients
    for (int k = 0; k < number_of_coeff; k++) {
      int sign = filter_coefficients[dir][k] > 0 ? 1 : -1;
      int coeff =
          (int)(filter_coefficients[dir][k] * (double)sign * (double)scale +
                0.5);
      f[k + dir * number_of_coeff] =
          (coeff > max_coefficient_value) ? max_coefficient_value : coeff;
      f[k + dir * number_of_coeff] *= sign;
    }
  }
}

int matrix_decomposition(
    double auto_correlation[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                           [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
    double upper_matrix[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                       [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
    int number_of_coeff) {
  double q[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];
  double reg_squate = 0.0000001;
  for (int i = 0; i < number_of_coeff; i++) {
    for (int j = i; j < number_of_coeff; j++) {
      double sq = auto_correlation[i][j];

      if (i > 0) {
        for (int k = i - 1; k >= 0; k--) {
          sq -= upper_matrix[k][j] * upper_matrix[k][i];
        }
      }
      if (i == j) {
        if (sq <= reg_squate) {
          return 0;
        } else {
          q[i] = 1.0 / (upper_matrix[i][i] = sqrt(sq));
        }
      } else {
        upper_matrix[i][j] = sq * q[i];
        upper_matrix[j][i] = 0.0;
      }
    }
  }
  return 1;
}

int derive_ccdef_filter_coefficients(
    double auto_correlation[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                           [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
    double *cross_correlation, double *filter_coefficients,
    int number_of_coeff) {
  double vector_aux[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];
  double upper_matrix[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                     [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];

  int success = 1;

  if (matrix_decomposition(auto_correlation, upper_matrix, number_of_coeff)) {
    vector_aux[0] = cross_correlation[0] / upper_matrix[0][0];
    for (int i = 1; i < number_of_coeff; i++) {
      double sum = 0;
      for (int j = 0; j < i; j++) {
        sum += vector_aux[j] * upper_matrix[j][i];
      }
      vector_aux[i] = (cross_correlation[i] - sum) / upper_matrix[i][i];
    }

    number_of_coeff--;
    filter_coefficients[number_of_coeff] =
        vector_aux[number_of_coeff] /
        upper_matrix[number_of_coeff][number_of_coeff];

    for (int i = number_of_coeff - 1; i >= 0; i--) {
      double sum = 0;

      for (int j = i + 1; j <= number_of_coeff; j++) {
        sum += upper_matrix[i][j] * filter_coefficients[j];
      }
      filter_coefficients[i] = (vector_aux[i] - sum) / upper_matrix[i][i];
    }

  } else {
    success = 0;
    for (int i = 0; i < number_of_coeff; i++) {
      auto_correlation[i][i] += 0.0001;
    }

    success =
        matrix_decomposition(auto_correlation, upper_matrix, number_of_coeff);

    if (!success) {
      memset(filter_coefficients, 0, sizeof(double) * number_of_coeff);
      return 0;
    }

    vector_aux[0] = cross_correlation[0] / upper_matrix[0][0];
    for (int i = 1; i < number_of_coeff; i++) {
      double sum = 0;

      for (int j = 0; j < i; j++) {
        sum += vector_aux[j] * upper_matrix[j][i];
      }

      vector_aux[i] = (cross_correlation[i] - sum) / upper_matrix[i][i];
    }

    number_of_coeff--;
    filter_coefficients[number_of_coeff] =
        vector_aux[number_of_coeff] /
        upper_matrix[number_of_coeff][number_of_coeff];

    for (int i = number_of_coeff - 1; i >= 0; i--) {
      double sum = 0;

      for (int j = i + 1; j <= number_of_coeff; j++) {
        sum += upper_matrix[i][j] * filter_coefficients[j];
      }

      filter_coefficients[i] = (vector_aux[i] - sum) / upper_matrix[i][i];
    }
  }
  return success;
}

/* Search for the set of strengths that minimizes mse. */
uint64_t joint_strength_search_ccdef(int *best_lev, int nb_strengths,
                                     uint64_t mse[][TOTAL_CCCDEF_STRENGTHS],
                                     int sb_count) {
  uint64_t best_tot_mse;
  int i;
  best_tot_mse = (uint64_t)1 << 63;
  /* Greedy search: add one strength options at a time. */
  for (i = 0; i < nb_strengths; i++) {
    best_tot_mse = search_one_cccdef(best_lev, i, mse, sb_count);
  }
  /* Trying to refine the greedy search by reconsidering each
     already-selected option. */
  for (i = 0; i < 4 * nb_strengths; i++) {
    int j;
    for (j = 0; j < nb_strengths - 1; j++) best_lev[j] = best_lev[j + 1];
    best_tot_mse = search_one_cccdef(best_lev, nb_strengths - 1, mse, sb_count);
  }
  return best_tot_mse;
}

uint64_t search_one_cccdef(int *lev, int nb_strengths,
                           uint64_t mse[][TOTAL_CCCDEF_STRENGTHS],
                           int sb_count) {
  uint64_t tot_mse[TOTAL_CCCDEF_STRENGTHS];
  const int total_strengths = TOTAL_CCCDEF_STRENGTHS;
  int i, j;
  uint64_t best_tot_mse = (uint64_t)1 << 63;
  int best_id = 0;
  memset(tot_mse, 0, sizeof(tot_mse));
  for (i = 0; i < sb_count; i++) {
    int gi;
    uint64_t best_mse = (uint64_t)1 << 63;
    /* Find best mse among already selected options. */
    for (gi = 0; gi < nb_strengths; gi++) {
      if (mse[i][lev[gi]] < best_mse) {
        best_mse = mse[i][lev[gi]];
      }
    }
    /* Find best mse when adding each possible new option. */
    for (j = 0; j < total_strengths; j++) {
      uint64_t best = best_mse;
      if (mse[i][j] < best) best = mse[i][j];
      tot_mse[j] += best;
    }
  }
  for (j = 0; j < total_strengths; j++) {
    if (tot_mse[j] < best_tot_mse) {
      best_tot_mse = tot_mse[j];
      best_id = j;
    }
  }
  lev[nb_strengths] = best_id;
  return best_tot_mse;
}

int get_ccdef_coeff_golomb_length(short coeff) {
  int x = coeff + 1;
  int i = x;
  int length = 0;
  int total_bit_length = 0;

  while (i) {
    i >>= 1;
    ++length;
  }

  for (i = 0; i < length - 1; ++i) total_bit_length++;

  for (i = length - 1; i >= 0; --i) total_bit_length++;
  return total_bit_length;
}

int get_ccdef_coeff_total_bit_length(const short *const coeff,
                                     const int direction_merge_mode) {
  int bit_length = 0;

  int num_coeff =
      direction_merge_mode
          ? MAX_NUMBER_OF_CCCDEF_FILTER_COEFF
          : (MAX_NUMBER_OF_DIRECTIONS * MAX_NUMBER_OF_CCCDEF_FILTER_COEFF);

  for (int i = 0; i < num_coeff; i++) {
    int coeffabs = abs(coeff[i]);
    bit_length += get_ccdef_coeff_golomb_length(coeffabs);
    if (coeffabs) {
      bit_length++;  // sign bits
    }
  }

  return bit_length;
}
#endif
