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

#include <string.h>
#include <math.h>

#include "./aom_scale_rtcd.h"
#include "av1/common/dering.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/reconinter.h"
#include "av1/encoder/encoder.h"
#include "aom/aom_integer.h"

static double compute_dist(uint16_t *x, int xstride, uint16_t *y, int ystride,
                           int nhb, int nvb, int coeff_shift) {
  int i, j;
  double sum;
  sum = 0;
  for (i = 0; i < nvb << 3; i++) {
    for (j = 0; j < nhb << 3; j++) {
      double tmp;
      tmp = x[i * xstride + j] - y[i * ystride + j];
      sum += tmp * tmp;
    }
  }
  return sum / (double)(1 << 2 * coeff_shift);
}

int av1_dering_search(YV12_BUFFER_CONFIG *frame, const YV12_BUFFER_CONFIG *ref,
                      AV1_COMMON *cm, MACROBLOCKD *xd) {
  int r, c;
  int sbr, sbc;
  uint16_t *src;
  uint16_t *ref_coeff;
  dering_list dlist[MAX_MIB_SIZE * MAX_MIB_SIZE];
  int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS] = { { 0 } };
  int stride;
  int bsize[3];
  int dec[3];
  int pli;
  int level;
  int best_level;
  int dering_count;
  int coeff_shift = AOMMAX(cm->bit_depth - 8, 0);
  double best_tot_mse=0;
  int sb_count;
  int nvsb = (cm->mi_rows + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  int nhsb = (cm->mi_cols + MAX_MIB_SIZE - 1) / MAX_MIB_SIZE;
  int *sb_index = aom_malloc(nvsb * nhsb * sizeof(*sb_index));
  uint64_t(*mse)[DERING_STRENGTHS][CLPF_STRENGTHS] =
      aom_malloc(sizeof(*mse) * nvsb * nhsb);
  int clpf_damping = 3 + (cm->base_qindex >> 6);
  src = aom_memalign(32, sizeof(*src) * cm->mi_rows * cm->mi_cols * 64);
  ref_coeff =
      aom_memalign(32, sizeof(*ref_coeff) * cm->mi_rows * cm->mi_cols * 64);
  av1_setup_dst_planes(xd->plane, frame, 0, 0);
  for (pli = 0; pli < 3; pli++) {
    dec[pli] = xd->plane[pli].subsampling_x;
    bsize[pli] = OD_DERING_SIZE_LOG2 - dec[pli];
  }
  stride = cm->mi_cols << bsize[0];
  for (r = 0; r < cm->mi_rows << bsize[0]; ++r) {
    for (c = 0; c < cm->mi_cols << bsize[0]; ++c) {
#if CONFIG_AOM_HIGHBITDEPTH
      if (cm->use_highbitdepth) {
        src[r * stride + c] = CONVERT_TO_SHORTPTR(
            xd->plane[0].dst.buf)[r * xd->plane[0].dst.stride + c];
        ref_coeff[r * stride + c] =
            CONVERT_TO_SHORTPTR(ref->y_buffer)[r * ref->y_stride + c];
      } else {
#endif
        src[r * stride + c] =
            xd->plane[0].dst.buf[r * xd->plane[0].dst.stride + c];
        ref_coeff[r * stride + c] = ref->y_buffer[r * ref->y_stride + c];
#if CONFIG_AOM_HIGHBITDEPTH
      }
#endif
    }
  }
  /* Pick a base threshold based on the quantizer. The threshold will then be
     adjusted on a 64x64 basis. We use a threshold of the form T = a*Q^b,
     where a and b are derived empirically trying to optimize rate-distortion
     at different quantizer settings. */
  best_level = AOMMIN(
      MAX_DERING_LEVEL - 1,
      (int)floor(.5 +
                 .45 * pow(av1_ac_quant(cm->base_qindex, 0, cm->bit_depth) >>
                               (cm->bit_depth - 8),
                           0.6)));
  sb_count = 0;
  for (sbr = 0; sbr < nvsb; sbr++) {
    for (sbc = 0; sbc < nhsb; sbc++) {
      int nvb, nhb;
      int gi;
      DECLARE_ALIGNED(32, uint16_t, dst[MAX_MIB_SIZE * MAX_MIB_SIZE * 8 * 8]);
      DECLARE_ALIGNED(32, uint16_t,
                      tmp_dst[MAX_MIB_SIZE * MAX_MIB_SIZE * 8 * 8]);
      nhb = AOMMIN(MAX_MIB_SIZE, cm->mi_cols - MAX_MIB_SIZE * sbc);
      nvb = AOMMIN(MAX_MIB_SIZE, cm->mi_rows - MAX_MIB_SIZE * sbr);
      dering_count = sb_compute_dering_list(cm, sbr * MAX_MIB_SIZE,
                                            sbc * MAX_MIB_SIZE, dlist);
      if (dering_count == 0) continue;
      for (gi = 0; gi < DERING_STRENGTHS; gi++) {
        int threshold;
        DECLARE_ALIGNED(32, uint16_t, inbuf[OD_DERING_INBUF_SIZE]);
        uint16_t *in;
        int i, j;
        level = dering_level_table[gi];
        threshold = level << coeff_shift;
        for (r = 0; r < nvb << bsize[0]; r++) {
          for (c = 0; c < nhb << bsize[0]; c++) {
            dst[(r * MAX_MIB_SIZE << bsize[0]) + c] =
                src[((sbr * MAX_MIB_SIZE << bsize[0]) + r) * stride +
                    (sbc * MAX_MIB_SIZE << bsize[0]) + c];
          }
        }
        in = inbuf + OD_FILT_VBORDER * OD_FILT_BSTRIDE + OD_FILT_HBORDER;
        /* We avoid filtering the pixels for which some of the pixels to average
           are outside the frame. We could change the filter instead, but it
           would
           add special cases for any future vectorization. */
        for (i = 0; i < OD_DERING_INBUF_SIZE; i++)
          inbuf[i] = OD_DERING_VERY_LARGE;
        for (i = -OD_FILT_VBORDER * (sbr != 0);
             i < (nvb << bsize[0]) + OD_FILT_VBORDER * (sbr != nvsb - 1); i++) {
          for (j = -OD_FILT_HBORDER * (sbc != 0);
               j < (nhb << bsize[0]) + OD_FILT_HBORDER * (sbc != nhsb - 1);
               j++) {
            uint16_t *x;
            x = &src[(sbr * stride * MAX_MIB_SIZE << bsize[0]) +
                     (sbc * MAX_MIB_SIZE << bsize[0])];
            in[i * OD_FILT_BSTRIDE + j] = x[i * stride + j];
          }
        }
	for (i = 0; i < CLPF_STRENGTHS; i++) {
	  od_dering(tmp_dst, in, 0, dir, 0, dlist, dering_count, threshold,
		    i + (i == 3), clpf_damping, coeff_shift, 0);
	  copy_dering_16bit_to_16bit(dst, MAX_MIB_SIZE << bsize[0], tmp_dst,
				     dlist, dering_count, bsize[0]);
	  mse[sb_count][gi][i] = (int)compute_dist(
              dst, MAX_MIB_SIZE << bsize[0],
              &ref_coeff[(sbr * stride * MAX_MIB_SIZE << bsize[0]) +
                         (sbc * MAX_MIB_SIZE << bsize[0])],
              stride, nhb, nvb, coeff_shift);
          sb_index[sb_count] = MAX_MIB_SIZE * sbr * cm->mi_stride +
                                    MAX_MIB_SIZE * sbc;
	}
      }
      sb_count++;
    }
  }
  int i;
  int lev[DERING_REFINEMENT_LEVELS];
  int best_lev[DERING_REFINEMENT_LEVELS];
  int str[CLPF_REFINEMENT_LEVELS];
  int best_str[CLPF_REFINEMENT_LEVELS];
  best_tot_mse = (uint64_t)1 << 63;
  {
    int l0;
    for (l0=0;l0<DERING_STRENGTHS;l0++) {
      int l1;
      lev[0] = l0;
      for (l1=l0+1;l1<DERING_STRENGTHS;l1++) {
        int l2;
        lev[1] = l1;
        for (l2=l1+1;l2<DERING_STRENGTHS;l2++) {
          int l3;
          lev[2] = l2;
          for (l3=l2+1;l3<DERING_STRENGTHS;l3++) {
	    int cs0;
	    lev[3] = l3;
	    for (cs0=0; cs0 < CLPF_STRENGTHS; cs0++) {
	      int cs1;
	      str[0] = cs0;
	      for (cs1=cs0+1; cs1 < CLPF_STRENGTHS; cs1++) {
		double tot_mse = 0;
		str[1] = cs1;
		for (i=0;i<sb_count;i++) {
		  int gi;
		  int cs;
		  double best_mse = INT32_MAX;
		  for (gi = 0; gi < DERING_REFINEMENT_LEVELS; gi++) {
		    for (cs = 0; cs < CLPF_REFINEMENT_LEVELS; cs++) {
		      if (mse[i][lev[gi]][str[cs]] < best_mse) {
			best_mse = mse[i][lev[gi]][str[cs]];
		      }
		    }
		  }
		  tot_mse += best_mse;
		}
		if (tot_mse < best_tot_mse) {
		  for (i=0;i<DERING_REFINEMENT_LEVELS;i++) best_lev[i] = lev[i];
		  for (i=0;i<CLPF_REFINEMENT_LEVELS;i++) best_str[i] = str[i];
		  best_tot_mse = tot_mse;
		}
	      }
	    }
          }
        }
      }
    }
  }
  for (i=0;i<DERING_REFINEMENT_LEVELS;i++) lev[i] = best_lev[i];
  for (i=0;i<CLPF_REFINEMENT_LEVELS;i++) str[i] = best_str[i];
  for (i=0;i<sb_count;i++) {
    int gi, cs;
    int best_gi, best_clpf;
    double best_mse = INT32_MAX;
    best_gi = best_clpf = 0;
    for (gi = 0; gi < DERING_REFINEMENT_LEVELS; gi++) {
      for (cs = 0; cs < CLPF_REFINEMENT_LEVELS; cs++) {
	if (mse[i][lev[gi]][str[cs]] < best_mse) {
	  best_gi = gi;
	  best_clpf = cs;
	  best_mse = mse[i][lev[gi]][str[cs]];
	}
      }
    }
    cm->mi_grid_visible[sb_index[i]]
        ->mbmi.dering_gain = best_gi;
    cm->mi_grid_visible[sb_index[i]]
      ->mbmi.clpf_strength = best_clpf;
  }
  aom_free(src);
  aom_free(ref_coeff);
  aom_free(mse);
  aom_free(sb_index);
  return levels_to_id(best_lev, best_str);
}
