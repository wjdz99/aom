/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <stdio.h>
#include <limits.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"
#include "config/aom_scale_rtcd.h"

#include "aom/aom_integer.h"
#include "av1/common/av1_common_int.h"
#include "av1/common/blockd.h"
#include "av1/common/reconinter.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/rotation.h"


static void rotational_predictor(const uint8_t * src, int src_stride, uint8_t * dst,
                int dst_stride, const SubpelParams *subpel_params,
                int w, int h, ConvolveParams *conv_params,
                const InterpFilterParams *interp_filters[2],  int subsampling_x, int subsampling_y, int rot) {
  SubpelParams sp = *subpel_params;
  revert_scale_extra_bits(&sp);

  const InterpFilterParams *filter_params_x = interp_filters[0];
  const InterpFilterParams *filter_params_y = interp_filters[1];

  //   sp->subpel_x >>= SCALE_EXTRA_BITS;
  //sp->subpel_y >>= SCALE_EXTRA_BITS;

    for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
    const uint8_t *this_src  = src +  x;
    MV32 this_mv;
    this_mv.row = rotation_taps[(rot + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 0];
    this_mv.col = rotation_taps[(rot + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 1];
    this_mv.row *= (1 << (1 - subsampling_y));
    this_mv.col *= (1 << (1 - subsampling_x));


    sp.subpel_x += this_mv.col;
    sp.subpel_y += this_mv.row;
    this_src += (sp.subpel_y >> SUBPEL_BITS) * src_stride
                  + (sp.subpel_x >> SUBPEL_BITS);
    sp.subpel_x = sp.subpel_x & SUBPEL_MASK;
    sp.subpel_y = sp.subpel_y & SUBPEL_MASK;



  int16_t im_buf[MAX_FILTER_TAP] = { 0 };
  assert(w <= MAX_SB_SIZE && h <= MAX_SB_SIZE);

  const int fo_vert = filter_params_y->taps / 2 - 1;
  const int fo_horiz = filter_params_x->taps / 2 - 1;
  const int bd = 8;
  const int bits =
      FILTER_BITS * 2 - conv_params->round_0 - conv_params->round_1;

  // horizontal filter
  const uint8_t *src_horiz = this_src - fo_vert * src_stride;
  const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
      filter_params_x, sp.subpel_x & SUBPEL_MASK);
  for (int r = 0; r < MAX_FILTER_TAP; ++r) {
      int32_t sum = (1 << (bd + FILTER_BITS - 1));
      for (int k = 0; k < filter_params_x->taps; ++k) {
        sum += x_filter[k] * src_horiz[r * src_stride - fo_horiz + k];
      }
      assert(0 <= sum && sum < (1 << (bd + FILTER_BITS + 1)));
      im_buf[r] =
          (int16_t)ROUND_POWER_OF_TWO(sum, conv_params->round_0);
  }

  // vertical filter
  const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
      filter_params_y, sp.subpel_y & SUBPEL_MASK);
  const int offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;

      int32_t sum = 1 << offset_bits;
      for (int k = 0; k < filter_params_y->taps; ++k) {
        sum += y_filter[k] * im_buf[k];
      }
      assert(0 <= sum && sum < (1 << (offset_bits + 2)));
      int16_t res = ROUND_POWER_OF_TWO(sum, conv_params->round_1) -
                    ((1 << (offset_bits - conv_params->round_1)) +
                     (1 << (offset_bits - conv_params->round_1 - 1)));
      dst[y * dst_stride + x] = clip_pixel(ROUND_POWER_OF_TWO(res, bits));

  }
  src += src_stride;
  }

}

static void highbd_rotational_predictor(const uint8_t *src, int src_stride, uint8_t * dst,
                int dst_stride, const SubpelParams *subpel_params,
                int w, int h, ConvolveParams *conv_params,
                const InterpFilterParams *interp_filters[2], int bd, int subsampling_x, int subsampling_y, int rot){
(void) src;
(void) src_stride;
(void) dst;
(void) dst_stride;
(void) subpel_params;
(void) w;
(void) h;
(void) conv_params;
(void) interp_filters;
(void) bd;
(void) subsampling_x;
(void) subsampling_y;
(void) rot;

return;
}


void av1_enc_build_one_rotational_predictor(uint8_t *dst, int dst_stride,
                                       const MV *src_mv,  const int rot,
                                       InterPredParams *inter_pred_params) {
  SubpelParams subpel_params;
  uint8_t *src;
  int src_stride;
  enc_calc_subpel_params(src_mv, inter_pred_params, NULL, 0, 0, inter_pred_params->conv_params.do_average,
                          NULL, &src, &subpel_params, &src_stride);

  assert(IMPLIES(inter_pred_params->conv_params.is_compound,
                 inter_pred_params->conv_params.dst != NULL));


#if CONFIG_AV1_HIGHBITDEPTH
    if (inter_pred_params->use_hbd_buf) {
      highbd_rotational_predictor(src, src_stride, dst, dst_stride, &subpel_params,
                             inter_pred_params->block_width,
                             inter_pred_params->block_height,
                             &inter_pred_params->conv_params,
                             inter_pred_params->interp_filter_params,
                             inter_pred_params->bit_depth, inter_pred_params->subsampling_x, inter_pred_params->subsampling_y, rot);
    } else {
      rotational_predictor(src, src_stride, dst, dst_stride, &subpel_params,
                      inter_pred_params->block_width,
                      inter_pred_params->block_height,
                      &inter_pred_params->conv_params,
                      inter_pred_params->interp_filter_params, inter_pred_params->subsampling_x, inter_pred_params->subsampling_y, rot);
    }
#else
    rotational_predictor(src, src_stride, dst, dst_stride, &subpel_params,
                    inter_pred_params->block_width,
                    inter_pred_params->block_height,
                    &inter_pred_params->conv_params,
                    inter_pred_params->interp_filter_params, inter_pred_params->subsampling_x, inter_pred_params->subsampling_y, rot);
#endif

}

#if 0
static void rotational_predictor(uint8_t * src, int src_stride, uint8_t * dst, int dst_stride, SubpelParams *subpel_params,
                int w, int h, ConvolveParams *conv_params, const InterpFilterParams *interp_filters[2], int rot) {
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const InterpKernel *kernel = av1_filter_kernels[mi->mbmi.interp_filter];
  int ref = 0; // single ref


  // mv impact is already done in enc_calc_subpel_params. Here, only consider
  // rotational angle.


  int subpel_x, subpel_y;
  int x, y;

  // each pixel in the block
  for (y = 0; y < h; ++y) {
    for (x = 0; x < w; ++x) {
      uint8_t tmp[SUBPEL_TAPS];
      MV32 this_mv;
      const uint8_t *this_pre = src + x;


//      this_mv.row = mv.row + rotation_taps[rot + 16][y - h/2 + 32][x - w/2 + 32][0];
//      this_mv.col = mv.col + rotation_taps[rot + 16][y - h/2 + 32][x - w/2 + 32][1];
      this_mv.row = rotation_taps[(rot + 16) * 65 * 65 * 2 + (y - h/2 + 32) * 65 * 2 + (x - w/2 + 32) * 2 + 0];
      this_mv.col = rotation_taps[(rot + 16) * 65 * 65 * 2 + (y - h/2 + 32) * 65 * 2 + (x - w/2 + 32) * 2 + 1];

      this_mv.row *= (1 << (1 - pd->subsampling_y));
      this_mv.col *= (1 << (1 - pd->subsampling_x));

      subpel_x = this_mv.col & SUBPEL_MASK;
      subpel_y = this_mv.row & SUBPEL_MASK;
      this_pre += (this_mv.row >> SUBPEL_BITS) * src_stride
                  + (this_mv.col >> SUBPEL_BITS);

      // horiz
      {
        const int16_t *filter_x = kernel[subpel_x];
        int r;

        this_pre += - src_stride * 3 - 3;

        for (r = 0; r < SUBPEL_TAPS; ++r) {
          int k, sum = 0;
          for (k = 0; k < SUBPEL_TAPS; ++k) {
            sum += this_pre[k] * filter_x[k];
          }
          tmp[r] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
          this_pre += src_stride;
        }
      }

      // vert
      {
        const int16_t *filter_y = kernel[subpel_y];
        int k, sum = 0;
        for (k = 0; k < SUBPEL_TAPS; ++k) {
          sum += tmp[k] * filter_y[k];
          dst[y * dst_stride + x] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
        }
      }
    }
    src += src_stride;
  }
}
#endif


static const MV search_table[24] = {{0, -4}, {4, 0}, {0, 4}, {-4, 0},
                         {-4, -4}, {4, -4},{4, 4},{-4, 4},
                         {0, -2}, {2, 0}, {0, 2}, {-2, 0},
                         {-2, -2}, {2, -2},{2, 2},{-2, 2},
                         {0, -1}, {1, 0}, {0, 1}, {-1, 0},
                         {-1, -1}, {1, -1},{1, 1},{-1, 1}};

// this needs to use upsampled ref. so doesn't work now.
unsigned int av1_find_best_rotation(MACROBLOCKD *xd,
                           const SUBPEL_MOTION_SEARCH_PARAMS *ms_params,
                           unsigned int pre_besterr,
                           MV *bestmv,
                           int *best_rot) {
  const uint8_t *const src = ms_params->var_params.ms_buffers.src->buf;
  const int src_stride = ms_params->var_params.ms_buffers.src->stride;
  const uint8_t *const pre = ms_params->var_params.ms_buffers.ref->buf;
  const int pre_stride = ms_params->var_params.ms_buffers.ref->stride;
  const int w = ms_params->var_params.w;
  const int h = ms_params->var_params.h;
  const int subsampling_y = xd->plane[0].subsampling_y;
  const int subsampling_x = xd->plane[0].subsampling_x;
  // starting location
  int br = bestmv->row;
  int bc = bestmv->col;

  const SUBPEL_SEARCH_VAR_PARAMS *var_params = &ms_params->var_params;
  const SUBPEL_SEARCH_TYPE subpel_search_type = var_params->subpel_search_type;
  const InterpFilterParams *filter = av1_get_filter(subpel_search_type);
  const aom_variance_fn_ptr_t *vfp = var_params->vfp;

  uint8_t dst_buf[MAX_SB_SIZE * MAX_SB_SIZE];
  const int dst_stride = w;


  int subpel_x_q4, subpel_y_q4;
  int x, y;
  int t, n, m;


//  const int minc = VPXMAX(x->mv_col_min * 8,  - MV_MAX);
//  const int maxc = VPXMIN(x->mv_col_max * 8,  MV_MAX);
//  const int minr = VPXMAX(x->mv_row_min * 8, - MV_MAX);
//  const int maxr = VPXMIN(x->mv_row_max * 8,  + MV_MAX);
//  const int bwl = bw << 2;
//  const int bhl = bh << 2;
//  sin 12.8 < 0.25, to get 0.5(round to 1), then use bwl >> 1 to set limits
//  const int bd = VPXMAX(bwl, bhl) >> 2;

  unsigned int besterr = pre_besterr;
  *best_rot = 0;

//  {
//  // calculate this for each angle
//  if (bc < (minc + bd) || bc > (maxc - bd)  ||
//      br < (minr + bd) || br > (maxr - bd)) {
//    return INT_MAX;
//  }

  // calculate baseline error
  {
      t = 0;

      const uint8_t *pre0 = pre;
      // each pixel in the block
          for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
              uint8_t tmp[SUBPEL_TAPS];
              const uint8_t *this_pre = pre0 + x;

              int nr, nc;
              nr = br;
              nc = bc;
              nr *= (1 << (1 - subsampling_y));
              nc *= (1 << (1 - subsampling_x));

              subpel_x_q4 = nc & SUBPEL_MASK;
              subpel_y_q4 = nr & SUBPEL_MASK;
              this_pre += (nr >> SUBPEL_BITS) * pre_stride + (nc >> SUBPEL_BITS);

              // horiz
              {
                const int16_t *const filter_x =
                      av1_get_interp_filter_subpel_kernel(filter, subpel_x_q4);
                int r;

                this_pre += - pre_stride * (SUBPEL_TAPS / 2 - 1) - (SUBPEL_TAPS / 2 - 1);

                for (r = 0; r < SUBPEL_TAPS; ++r) {
                  int k, sum = 0;
                  for (k = 0; k < SUBPEL_TAPS; ++k) {
                    sum += this_pre[k] * filter_x[k];
                  }
                  tmp[r] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
                  this_pre += pre_stride;
                }
              }

              // vert
              {
                const int16_t *const filter_y =
                      av1_get_interp_filter_subpel_kernel(filter, subpel_y_q4);

                int k, sum = 0;
                for (k = 0; k < SUBPEL_TAPS; ++k) {
                  sum += tmp[k] * filter_y[k];
                }
                dst_buf[y * dst_stride + x] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
              }
            }
            pre0 += pre_stride;
          }

          unsigned int sse = 0;
          unsigned int var = vfp->vf(dst_buf, dst_stride, src, src_stride, &sse);
          besterr = var;
    }

  // rot search
  for (n = 1; n <= 16 ; n++) {
    int rounds = n ? 2 : 1;
    for (m = 0; m < rounds; m++) {

      if (m == 0) t = n;
      else t = -n;

     // printf("\n t = %d; (%d;%d) \n", t,br,bc);

      const uint8_t *pre0 = pre;
      // each pixel in the block
          for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
              uint8_t tmp[SUBPEL_TAPS];
              const uint8_t *this_pre = pre0 + x;

              const int rr = rotation_taps[(t + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 0];
              const int rc = rotation_taps[(t + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 1];

              //printf(" (%d,%d) ", rr, rc);

              int nr, nc;
              nr = br + rr;
              nc = bc + rc;
              nr *= (1 << (1 - subsampling_y));
              nc *= (1 << (1 - subsampling_x));

              subpel_x_q4 = nc & SUBPEL_MASK;
              subpel_y_q4 = nr & SUBPEL_MASK;
              this_pre += (nr >> SUBPEL_BITS) * pre_stride + (nc >> SUBPEL_BITS);

              // horiz
              {
                const int16_t *const filter_x =
                      av1_get_interp_filter_subpel_kernel(filter, subpel_x_q4);
                int r;

                this_pre += - pre_stride * (SUBPEL_TAPS / 2 - 1) - (SUBPEL_TAPS / 2 - 1);

                for (r = 0; r < SUBPEL_TAPS; ++r) {
                  int k, sum = 0;
                  for (k = 0; k < SUBPEL_TAPS; ++k) {
                    sum += this_pre[k] * filter_x[k];
                  }
                  tmp[r] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
                  this_pre += pre_stride;
                }
              }

              // vert
              {
                const int16_t *const filter_y =
                      av1_get_interp_filter_subpel_kernel(filter, subpel_y_q4);

                int k, sum = 0;
                for (k = 0; k < SUBPEL_TAPS; ++k) {
                  sum += tmp[k] * filter_y[k];
                }
                dst_buf[y * dst_stride + x] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
              }
            }
            pre0 += pre_stride;
          }

          unsigned int sse = 0;
          unsigned int var = vfp->vf(dst_buf, dst_stride, src, src_stride, &sse);

//         {
//            int i, j;
//            printf("\nsrc:\n");
//            for (i = 0; i < h; ++i) {
//              for (j = 0; j < w; ++j) {
//                printf("%d, ", src[i * src_stride +j]);
//              }
//              printf("\n");
//            }
//            printf("\ndst:\n");
//            for (i = 0; i < h; ++i) {
//              for (j = 0; j < w; ++j) {
//                printf("%d, ", abs(dst_buf[i * dst_stride +j] - src[i * src_stride +j]));
//              }
//              printf("\n");
//            }
//
//          }

      if (!t) {
        besterr = var;
      }
      if (var < besterr) {
        *best_rot = t;
        besterr = var;
      }
    }
  }
//  printf("\n  111111:  rot: %d; mv: %d,%d; besterr: %d; ", *best_rot, br, bc, besterr);


  //refine MV
  if (*best_rot)
    {
    int start = 8;

    // iterate 3 times
    for (int i = 0; i < 3; i++) {
      int best_idx = -1;

      for (m = 0; m < 8; m++) {
        int tr = br + search_table[start + m].row;
        int tc = bc + search_table[start + m].col;
        const uint8_t *pre0 = pre;

        for (y = 0; y < h; y++) {
          for (x = 0; x < w; x++) {
            uint8_t tmp[SUBPEL_TAPS];
            const uint8_t *this_pre = pre0 + x;
            int nr, nc;

            nr = tr + rotation_taps[(*best_rot + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 0];
            nc = tc + rotation_taps[(*best_rot + 16) * 32 * 32 * 2 + (y - h/2 + 16) * 32 * 2 + (x - w/2 + 16) * 2 + 1];
            nr *= (1 << (1 - subsampling_y));
            nc *= (1 << (1 - subsampling_x));

            subpel_x_q4 = nc & SUBPEL_MASK;
            subpel_y_q4 = nr & SUBPEL_MASK;
            this_pre += (nr >> SUBPEL_BITS) * pre_stride + (nc >> SUBPEL_BITS);

            // horiz
            {
              const int16_t *const filter_x =
                    av1_get_interp_filter_subpel_kernel(filter, subpel_x_q4);
              int r;

              this_pre += - pre_stride * (SUBPEL_TAPS / 2 - 1) - (SUBPEL_TAPS / 2 - 1);

              for (r = 0; r < SUBPEL_TAPS; ++r) {
                int k, sum = 0;
                for (k = 0; k < SUBPEL_TAPS; ++k) {
                  sum += this_pre[k] * filter_x[k];
                }
                tmp[r] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
                this_pre += pre_stride;
              }
            }

            // vert
            {
              const int16_t *const filter_y =
                    av1_get_interp_filter_subpel_kernel(filter, subpel_y_q4);

              int k, sum = 0;
              for (k = 0; k < SUBPEL_TAPS; ++k) {
                sum += tmp[k] * filter_y[k];
              }
              dst_buf[y * dst_stride + x] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
            }
          }
          pre0 += pre_stride;
        }

        unsigned int sse = 0;
        unsigned int var = vfp->vf(dst_buf, dst_stride, src, src_stride, &sse);

        if (var < besterr) {
          best_idx = m;
          besterr = var;
        }
      }

      if(best_idx == -1)
        break;

      if(best_idx >= 0) {
        br += search_table[start + best_idx].row;
        bc += search_table[start + best_idx].col;
      }
    }

    bestmv->row = br;
    bestmv->col = bc;
    }

//  printf("\n  222222:  rot: %d; mv: %d,%d; besterr: %d; ", *best_rot, br, bc, besterr);


  return besterr;
}
