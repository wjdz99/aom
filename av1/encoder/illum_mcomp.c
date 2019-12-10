/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"
#include "aom_ports/system_state.h"

#include "av1/common/common.h"
#include "av1/common/mvref_common.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/pred_common.h"
#include "av1/common/reconinter.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/partition_strategy.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/reconinter_enc.h"

static INLINE int mv_cost(const MV *mv, const MV *ref, const int *joint_cost,
                          int *const (*mvcost)[2],
                          MvSubpelPrecision precision) {
  (void)ref;
  return joint_cost[av1_get_mv_joint(mv)] + mvcost[precision][0][mv->row] +
         mvcost[precision][1][mv->col];
}

static int mvsad_err_cost(const MACROBLOCK *x, const MV *mv, const MV *ref,
                          MvSubpelPrecision precision, int sad_per_bit) {
  const MV diff = { (mv->row - ref->row) * 8, (mv->col - ref->col) * 8 };
  return ROUND_POWER_OF_TWO(
      (unsigned)mv_cost(&diff, ref, x->nmv_vec_cost, x->nmvcost, precision) *
          sad_per_bit,
      AV1_PROB_COST_SHIFT);
}

typedef int(mv_err_cost_fn)(
    const MV *mv, const MV *ref, MvSubpelPrecision max_precision,
    MvSubpelPrecision min_precision, const int *mvjcost,
    int *const (*mvcost)[2],
    const int (*dummy)[MV_SUBPEL_PRECISIONS - DISALLOW_ONE_DOWN_FLEX_MVRES],
    int error_per_bit);
#define PIXEL_TRANSFORM_ERROR_SCALE 4
static int mv_base_err_cost(
    const MV *mv, const MV *ref, MvSubpelPrecision max_precision,
    MvSubpelPrecision min_precision, const int *mvjcost,
    int *const (*mvcost)[2],
    const int (*dummy)[MV_SUBPEL_PRECISIONS - DISALLOW_ONE_DOWN_FLEX_MVRES],
    int error_per_bit) {
  (void)dummy;
  (void)min_precision;
  if (mvcost) {
    MV ref_ = *ref;
    const MV diff = { mv->row - ref_.row, mv->col - ref_.col };
    const int cost = mv_cost(&diff, &ref_, mvjcost, mvcost, max_precision);
    return (int)ROUND_POWER_OF_TWO_64((int64_t)cost * error_per_bit,
                                      RDDIV_BITS + AV1_PROB_COST_SHIFT -
                                          RD_EPB_SHIFT +
                                          PIXEL_TRANSFORM_ERROR_SCALE);
  }
  return 0;
}

static INLINE const uint8_t *get_buf_from_mv(const struct buf_2d *buf,
                                             const MV *mv) {
  return &buf->buf[mv->row * buf->stride + mv->col];
}

static INLINE int check_bounds(const MvLimits *mv_limits, int row, int col,
                               int range) {
  return ((row - range) >= mv_limits->row_min) &
         ((row + range) <= mv_limits->row_max) &
         ((col - range) >= mv_limits->col_min) &
         ((col + range) <= mv_limits->col_max);
}

static INLINE int is_mv_in(const MvLimits *mv_limits, const MV *mv) {
  return (mv->col >= mv_limits->col_min) && (mv->col <= mv_limits->col_max) &&
         (mv->row >= mv_limits->row_min) && (mv->row <= mv_limits->row_max);
}

#define CHECK_BETTER                                                          \
  {                                                                           \
    if (thissad < bestsad) {                                                  \
      if (use_mvcost)                                                         \
        thissad +=                                                            \
            mvsad_err_cost(x, &this_mv, &fcenter_mv, precision, sad_per_bit); \
      if (thissad < bestsad) {                                                \
        bestsad = thissad;                                                    \
        best_site = i;                                                        \
      }                                                                       \
    }                                                                         \
  }

#define MAX_PATTERN_SCALES 11
#define MAX_PATTERN_CANDIDATES 8  // max number of canddiates per scale
#define PATTERN_CANDIDATES_REF 3  // number of refinement candidates

// Calculate and return a sad+mvcost list around an integer best pel.
static INLINE void calc_int_cost_list(const AV1_COMMON *cm, const MACROBLOCK *x,
                                      const MV *const ref_mv, int sadpb,
                                      const aom_variance_fn_ptr_t *fn_ptr,
                                      const MV *best_mv, int *cost_list) {
  static const MV neighbors[4] = { { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } };
  const struct buf_2d *const what = &x->plane[0].src;
  const struct buf_2d *const in_what = &x->e_mbd.plane[0].pre[0];
  const MV fcenter_mv = { ref_mv->row >> 3, ref_mv->col >> 3 };
  const int br = best_mv->row;
  const int bc = best_mv->col;
  int i;
  unsigned int sse;
  const MV this_mv = { br, bc };
  const int(
      *flex_mv_costs)[MV_SUBPEL_PRECISIONS - DISALLOW_ONE_DOWN_FLEX_MVRES];
  const MACROBLOCKD *xd = &x->e_mbd;
  const MB_MODE_INFO *mbmi = xd->mi[0];
#if CONFIG_FLEX_MVRES
  const int use_flex_mv =
      is_flex_mv_precision_active(cm, mbmi->mode, mbmi->max_mv_precision);
  const int down_ctx = av1_get_mv_precision_down_context(cm, xd);
  flex_mv_costs = x->flex_mv_precision_costs[down_ctx];
  mv_err_cost_fn *mv_err_cost =
      use_flex_mv ? mv_flex_err_cost : mv_base_err_cost;
#else
  (void)cm;
  flex_mv_costs = NULL;
  mv_err_cost_fn *mv_err_cost = mv_base_err_cost;
#endif  // CONFIG_FLEX_MVRES

  cost_list[0] =
      fn_ptr->vf(what->buf, what->stride, get_buf_from_mv(in_what, &this_mv),
                 in_what->stride, &sse) +
      mvsad_err_cost(x, &this_mv, &fcenter_mv, mbmi->max_mv_precision, sadpb);
  if (check_bounds(&x->mv_limits, br, bc, 1)) {
    for (i = 0; i < 4; i++) {
      const MV neighbor_mv = { br + neighbors[i].row, bc + neighbors[i].col };
      cost_list[i + 1] =
          fn_ptr->vf(what->buf, what->stride,
                     get_buf_from_mv(in_what, &neighbor_mv), in_what->stride,
                     &sse) +
          mv_err_cost(&neighbor_mv, &fcenter_mv, mbmi->max_mv_precision,
                      MV_SUBPEL_NONE, x->nmv_vec_cost, x->nmvcost,
                      flex_mv_costs, x->errorperbit);
    }
  } else {
    for (i = 0; i < 4; i++) {
      const MV neighbor_mv = { br + neighbors[i].row, bc + neighbors[i].col };
      if (!is_mv_in(&x->mv_limits, &neighbor_mv))
        cost_list[i + 1] = INT_MAX;
      else
        cost_list[i + 1] =
            fn_ptr->vf(what->buf, what->stride,
                       get_buf_from_mv(in_what, &neighbor_mv), in_what->stride,
                       &sse) +
            mv_err_cost(&neighbor_mv, &fcenter_mv, mbmi->max_mv_precision,
                        MV_SUBPEL_NONE, x->nmv_vec_cost, x->nmvcost,
                        flex_mv_costs, x->errorperbit);
    }
  }
}

static INLINE void calc_int_sad_list(const AV1_COMMON *const cm,
                                     const MACROBLOCK *x,
                                     const MV *const ref_mv, int sadpb,
                                     const aom_variance_fn_ptr_t *fn_ptr,
                                     const MV *best_mv, int *cost_list,
                                     const int use_mvcost, const int bestsad) {
  (void)cm;
  static const MV neighbors[4] = { { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } };
  const struct buf_2d *const what = &x->plane[0].src;
  const struct buf_2d *const in_what = &x->e_mbd.plane[0].pre[0];
  const MV fcenter_mv = { ref_mv->row >> 3, ref_mv->col >> 3 };
  int i;
  const int br = best_mv->row;
  const int bc = best_mv->col;
  const MACROBLOCKD *xd = &x->e_mbd;
  const MB_MODE_INFO *mbmi = xd->mi[0];

  if (cost_list[0] == INT_MAX) {
    cost_list[0] = bestsad;
    if (check_bounds(&x->mv_limits, br, bc, 1)) {
      for (i = 0; i < 4; i++) {
        const MV this_mv = { br + neighbors[i].row, bc + neighbors[i].col };
        cost_list[i + 1] =
            fn_ptr->sdf(what->buf, what->stride,
                        get_buf_from_mv(in_what, &this_mv), in_what->stride);
      }
    } else {
      for (i = 0; i < 4; i++) {
        const MV this_mv = { br + neighbors[i].row, bc + neighbors[i].col };
        if (!is_mv_in(&x->mv_limits, &this_mv))
          cost_list[i + 1] = INT_MAX;
        else
          cost_list[i + 1] =
              fn_ptr->sdf(what->buf, what->stride,
                          get_buf_from_mv(in_what, &this_mv), in_what->stride);
      }
    }
  } else {
    if (use_mvcost) {
      for (i = 0; i < 4; i++) {
        const MV this_mv = { br + neighbors[i].row, bc + neighbors[i].col };
        if (cost_list[i + 1] != INT_MAX) {
          cost_list[i + 1] += mvsad_err_cost(x, &this_mv, &fcenter_mv,
                                             mbmi->max_mv_precision, sadpb);
        }
      }
    }
  }
}

// Generic pattern search function that searches over multiple scales.
// Each scale can have a different number of candidates and shape of
// candidates as indicated in the num_candidates and candidates arrays
// passed into this function
//
static int pattern_search(
    const AV1_COMP *cpi, MACROBLOCK *x, MV *start_mv, int search_param,
    int sad_per_bit, int do_init_search, int *cost_list,
    const aom_variance_fn_ptr_t *vfp, int use_mvcost, const MV *center_mv,
    const int num_candidates[MAX_PATTERN_SCALES],
    const MV candidates[MAX_PATTERN_SCALES][MAX_PATTERN_CANDIDATES]) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MB_MODE_INFO *mbmi = xd->mi[0];
  static const int search_param_to_steps[MAX_MVSEARCH_STEPS] = {
    10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
  };
  int i, s, t;
  const struct buf_2d *const what = &x->plane[0].src;
  const struct buf_2d *const in_what = &xd->plane[0].pre[0];
  const int last_is_4 = num_candidates[0] == 4;
  int br, bc;
  int bestsad = INT_MAX;
  int thissad;
  int k = -1;
  const MvSubpelPrecision precision = mbmi->max_mv_precision;
  const MV fcenter_mv = { center_mv->row >> 3, center_mv->col >> 3 };
  assert(search_param < MAX_MVSEARCH_STEPS);
  int best_init_s = search_param_to_steps[search_param];
  // adjust ref_mv to make sure it is within MV range
  clamp_mv(start_mv, x->mv_limits.col_min, x->mv_limits.col_max,
           x->mv_limits.row_min, x->mv_limits.row_max);
  br = start_mv->row;
  bc = start_mv->col;
  if (cost_list != NULL) {
    cost_list[0] = cost_list[1] = cost_list[2] = cost_list[3] = cost_list[4] =
        INT_MAX;
  }

  // Work out the start point for the search
  bestsad = vfp->sdf(what->buf, what->stride,
                     get_buf_from_mv(in_what, start_mv), in_what->stride) +
            mvsad_err_cost(x, start_mv, &fcenter_mv, mbmi->max_mv_precision,
                           sad_per_bit);

  // Search all possible scales upto the search param around the center point
  // pick the scale of the point that is best as the starting scale of
  // further steps around it.
  if (do_init_search) {
    s = best_init_s;
    best_init_s = -1;
    for (t = 0; t <= s; ++t) {
      int best_site = -1;
      if (check_bounds(&x->mv_limits, br, bc, 1 << t)) {
        for (i = 0; i < num_candidates[t]; i++) {
          const MV this_mv = { br + candidates[t][i].row,
                               bc + candidates[t][i].col };
          thissad =
              vfp->sdf(what->buf, what->stride,
                       get_buf_from_mv(in_what, &this_mv), in_what->stride);
          CHECK_BETTER
        }
      } else {
        for (i = 0; i < num_candidates[t]; i++) {
          const MV this_mv = { br + candidates[t][i].row,
                               bc + candidates[t][i].col };
          if (!is_mv_in(&x->mv_limits, &this_mv)) continue;
          thissad =
              vfp->sdf(what->buf, what->stride,
                       get_buf_from_mv(in_what, &this_mv), in_what->stride);
          CHECK_BETTER
        }
      }
      if (best_site == -1) {
        continue;
      } else {
        best_init_s = t;
        k = best_site;
      }
    }
    if (best_init_s != -1) {
      br += candidates[best_init_s][k].row;
      bc += candidates[best_init_s][k].col;
    }
  }

  // If the center point is still the best, just skip this and move to
  // the refinement step.
  if (best_init_s != -1) {
    const int last_s = (last_is_4 && cost_list != NULL);
    int best_site = -1;
    s = best_init_s;

    for (; s >= last_s; s--) {
      // No need to search all points the 1st time if initial search was used
      if (!do_init_search || s != best_init_s) {
        if (check_bounds(&x->mv_limits, br, bc, 1 << s)) {
          for (i = 0; i < num_candidates[s]; i++) {
            const MV this_mv = { br + candidates[s][i].row,
                                 bc + candidates[s][i].col };
            thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        } else {
          for (i = 0; i < num_candidates[s]; i++) {
            const MV this_mv = { br + candidates[s][i].row,
                                 bc + candidates[s][i].col };
            if (!is_mv_in(&x->mv_limits, &this_mv)) continue;
            thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        }

        if (best_site == -1) {
          continue;
        } else {
          br += candidates[s][best_site].row;
          bc += candidates[s][best_site].col;
          k = best_site;
        }
      }

      do {
        int next_chkpts_indices[PATTERN_CANDIDATES_REF];
        best_site = -1;
        next_chkpts_indices[0] = (k == 0) ? num_candidates[s] - 1 : k - 1;
        next_chkpts_indices[1] = k;
        next_chkpts_indices[2] = (k == num_candidates[s] - 1) ? 0 : k + 1;

        if (check_bounds(&x->mv_limits, br, bc, 1 << s)) {
          for (i = 0; i < PATTERN_CANDIDATES_REF; i++) {
            const MV this_mv = {
              br + candidates[s][next_chkpts_indices[i]].row,
              bc + candidates[s][next_chkpts_indices[i]].col
            };
            thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        } else {
          for (i = 0; i < PATTERN_CANDIDATES_REF; i++) {
            const MV this_mv = {
              br + candidates[s][next_chkpts_indices[i]].row,
              bc + candidates[s][next_chkpts_indices[i]].col
            };
            if (!is_mv_in(&x->mv_limits, &this_mv)) continue;
            thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        }

        if (best_site != -1) {
          k = next_chkpts_indices[best_site];
          br += candidates[s][k].row;
          bc += candidates[s][k].col;
        }
      } while (best_site != -1);
    }

    // Note: If we enter the if below, then cost_list must be non-NULL.
    if (s == 0) {
      cost_list[0] = bestsad;
      if (!do_init_search || s != best_init_s) {
        if (check_bounds(&x->mv_limits, br, bc, 1 << s)) {
          for (i = 0; i < num_candidates[s]; i++) {
            const MV this_mv = { br + candidates[s][i].row,
                                 bc + candidates[s][i].col };
            cost_list[i + 1] = thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        } else {
          for (i = 0; i < num_candidates[s]; i++) {
            const MV this_mv = { br + candidates[s][i].row,
                                 bc + candidates[s][i].col };
            if (!is_mv_in(&x->mv_limits, &this_mv)) continue;
            cost_list[i + 1] = thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        }

        if (best_site != -1) {
          br += candidates[s][best_site].row;
          bc += candidates[s][best_site].col;
          k = best_site;
        }
      }
      while (best_site != -1) {
        int next_chkpts_indices[PATTERN_CANDIDATES_REF];
        best_site = -1;
        next_chkpts_indices[0] = (k == 0) ? num_candidates[s] - 1 : k - 1;
        next_chkpts_indices[1] = k;
        next_chkpts_indices[2] = (k == num_candidates[s] - 1) ? 0 : k + 1;
        cost_list[1] = cost_list[2] = cost_list[3] = cost_list[4] = INT_MAX;
        cost_list[((k + 2) % 4) + 1] = cost_list[0];
        cost_list[0] = bestsad;

        if (check_bounds(&x->mv_limits, br, bc, 1 << s)) {
          for (i = 0; i < PATTERN_CANDIDATES_REF; i++) {
            const MV this_mv = {
              br + candidates[s][next_chkpts_indices[i]].row,
              bc + candidates[s][next_chkpts_indices[i]].col
            };
            cost_list[next_chkpts_indices[i] + 1] = thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        } else {
          for (i = 0; i < PATTERN_CANDIDATES_REF; i++) {
            const MV this_mv = {
              br + candidates[s][next_chkpts_indices[i]].row,
              bc + candidates[s][next_chkpts_indices[i]].col
            };
            if (!is_mv_in(&x->mv_limits, &this_mv)) {
              cost_list[next_chkpts_indices[i] + 1] = INT_MAX;
              continue;
            }
            cost_list[next_chkpts_indices[i] + 1] = thissad =
                vfp->sdf(what->buf, what->stride,
                         get_buf_from_mv(in_what, &this_mv), in_what->stride);
            CHECK_BETTER
          }
        }

        if (best_site != -1) {
          k = next_chkpts_indices[best_site];
          br += candidates[s][k].row;
          bc += candidates[s][k].col;
        }
      }
    }
  }

  // Returns the one-away integer pel cost/sad around the best as follows:
  // cost_list[0]: cost/sad at the best integer pel
  // cost_list[1]: cost/sad at delta {0, -1} (left)   from the best integer pel
  // cost_list[2]: cost/sad at delta { 1, 0} (bottom) from the best integer pel
  // cost_list[3]: cost/sad at delta { 0, 1} (right)  from the best integer pel
  // cost_list[4]: cost/sad at delta {-1, 0} (top)    from the best integer pel
  if (cost_list) {
    const MV best_int_mv = { br, bc };
    if (last_is_4) {
      calc_int_sad_list(&cpi->common, x, center_mv, sad_per_bit, vfp,
                        &best_int_mv, cost_list, use_mvcost, bestsad);
    } else {
      calc_int_cost_list(&cpi->common, x, center_mv, sad_per_bit, vfp,
                         &best_int_mv, cost_list);
    }
  }
  x->best_mv.as_mv.row = br;
  x->best_mv.as_mv.col = bc;
  return bestsad;
}

static int bigdia_search(const AV1_COMP *cpi, MACROBLOCK *x, MV *start_mv,
                         int search_param, int sad_per_bit, int do_init_search,
                         int *cost_list, const aom_variance_fn_ptr_t *vfp,
                         int use_mvcost, const MV *center_mv) {
  // First scale has 4-closest points, the rest have 8 points in diamond
  // shape at increasing scales
  static const int bigdia_num_candidates[MAX_PATTERN_SCALES] = {
    4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  };
  // Note that the largest candidate step at each scale is 2^scale
  /* clang-format off */
  static const MV
      bigdia_candidates[MAX_PATTERN_SCALES][MAX_PATTERN_CANDIDATES] = {
        { { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } },
        { { -1, -1 }, { 0, -2 }, { 1, -1 }, { 2, 0 }, { 1, 1 }, { 0, 2 },
          { -1, 1 }, { -2, 0 } },
        { { -2, -2 }, { 0, -4 }, { 2, -2 }, { 4, 0 }, { 2, 2 }, { 0, 4 },
          { -2, 2 }, { -4, 0 } },
        { { -4, -4 }, { 0, -8 }, { 4, -4 }, { 8, 0 }, { 4, 4 }, { 0, 8 },
          { -4, 4 }, { -8, 0 } },
        { { -8, -8 }, { 0, -16 }, { 8, -8 }, { 16, 0 }, { 8, 8 }, { 0, 16 },
          { -8, 8 }, { -16, 0 } },
        { { -16, -16 }, { 0, -32 }, { 16, -16 }, { 32, 0 }, { 16, 16 },
          { 0, 32 }, { -16, 16 }, { -32, 0 } },
        { { -32, -32 }, { 0, -64 }, { 32, -32 }, { 64, 0 }, { 32, 32 },
          { 0, 64 }, { -32, 32 }, { -64, 0 } },
        { { -64, -64 }, { 0, -128 }, { 64, -64 }, { 128, 0 }, { 64, 64 },
          { 0, 128 }, { -64, 64 }, { -128, 0 } },
        { { -128, -128 }, { 0, -256 }, { 128, -128 }, { 256, 0 }, { 128, 128 },
          { 0, 256 }, { -128, 128 }, { -256, 0 } },
        { { -256, -256 }, { 0, -512 }, { 256, -256 }, { 512, 0 }, { 256, 256 },
          { 0, 512 }, { -256, 256 }, { -512, 0 } },
        { { -512, -512 }, { 0, -1024 }, { 512, -512 }, { 1024, 0 },
          { 512, 512 }, { 0, 1024 }, { -512, 512 }, { -1024, 0 } },
      };
  /* clang-format on */
  return pattern_search(cpi, x, start_mv, search_param, sad_per_bit,
                        do_init_search, cost_list, vfp, use_mvcost, center_mv,
                        bigdia_num_candidates, bigdia_candidates);
}

static int square_search(const AV1_COMP *cpi, MACROBLOCK *x, MV *start_mv,
                         int search_param, int sad_per_bit, int do_init_search,
                         int *cost_list, const aom_variance_fn_ptr_t *vfp,
                         int use_mvcost, const MV *center_mv) {
  // All scales have 8 closest points in square shape
  static const int square_num_candidates[MAX_PATTERN_SCALES] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  };
  // Note that the largest candidate step at each scale is 2^scale
  /* clang-format off */
  static const MV
      square_candidates[MAX_PATTERN_SCALES][MAX_PATTERN_CANDIDATES] = {
        { { -1, -1 }, { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 }, { 0, 1 },
          { -1, 1 }, { -1, 0 } },
        { { -2, -2 }, { 0, -2 }, { 2, -2 }, { 2, 0 }, { 2, 2 }, { 0, 2 },
          { -2, 2 }, { -2, 0 } },
        { { -4, -4 }, { 0, -4 }, { 4, -4 }, { 4, 0 }, { 4, 4 }, { 0, 4 },
          { -4, 4 }, { -4, 0 } },
        { { -8, -8 }, { 0, -8 }, { 8, -8 }, { 8, 0 }, { 8, 8 }, { 0, 8 },
          { -8, 8 }, { -8, 0 } },
        { { -16, -16 }, { 0, -16 }, { 16, -16 }, { 16, 0 }, { 16, 16 },
          { 0, 16 }, { -16, 16 }, { -16, 0 } },
        { { -32, -32 }, { 0, -32 }, { 32, -32 }, { 32, 0 }, { 32, 32 },
          { 0, 32 }, { -32, 32 }, { -32, 0 } },
        { { -64, -64 }, { 0, -64 }, { 64, -64 }, { 64, 0 }, { 64, 64 },
          { 0, 64 }, { -64, 64 }, { -64, 0 } },
        { { -128, -128 }, { 0, -128 }, { 128, -128 }, { 128, 0 }, { 128, 128 },
          { 0, 128 }, { -128, 128 }, { -128, 0 } },
        { { -256, -256 }, { 0, -256 }, { 256, -256 }, { 256, 0 }, { 256, 256 },
          { 0, 256 }, { -256, 256 }, { -256, 0 } },
        { { -512, -512 }, { 0, -512 }, { 512, -512 }, { 512, 0 }, { 512, 512 },
          { 0, 512 }, { -512, 512 }, { -512, 0 } },
        { { -1024, -1024 }, { 0, -1024 }, { 1024, -1024 }, { 1024, 0 },
          { 1024, 1024 }, { 0, 1024 }, { -1024, 1024 }, { -1024, 0 } },
      };
  /* clang-format on */
  return pattern_search(cpi, x, start_mv, search_param, sad_per_bit,
                        do_init_search, cost_list, vfp, use_mvcost, center_mv,
                        square_num_candidates, square_candidates);
}

static int fast_hex_search(const AV1_COMP *cpi, MACROBLOCK *x, MV *ref_mv,
                           int search_param, int sad_per_bit,
                           int do_init_search,  // must be zero for fast_hex
                           int *cost_list, const aom_variance_fn_ptr_t *vfp,
                           int use_mvcost, const MV *center_mv) {
  return av1_hex_search(
      cpi, x, ref_mv, AOMMAX(MAX_MVSEARCH_STEPS - 2, search_param), sad_per_bit,
      do_init_search, cost_list, vfp, use_mvcost, center_mv);
}

static int fast_dia_search(const AV1_COMP *cpi, MACROBLOCK *x, MV *ref_mv,
                           int search_param, int sad_per_bit,
                           int do_init_search, int *cost_list,
                           const aom_variance_fn_ptr_t *vfp, int use_mvcost,
                           const MV *center_mv) {
  return bigdia_search(
      cpi, x, ref_mv, AOMMAX(MAX_MVSEARCH_STEPS - 2, search_param), sad_per_bit,
      do_init_search, cost_list, vfp, use_mvcost, center_mv);
}

#undef CHECK_BETTER

// Exhuastive motion search around a given centre position with a given
// step size.
static int exhuastive_mesh_search(const AV1_COMMON *cm, MACROBLOCK *x,
                                  MV *ref_mv, MV *best_mv, int range, int step,
                                  int sad_per_bit,
                                  const aom_variance_fn_ptr_t *fn_ptr,
                                  const MV *center_mv) {
  (void)cm;
  const MACROBLOCKD *const xd = &x->e_mbd;
  const MB_MODE_INFO *mbmi = xd->mi[0];
  const struct buf_2d *const what = &x->plane[0].src;
  const struct buf_2d *const in_what = &xd->plane[0].pre[0];
  MV fcenter_mv = { center_mv->row, center_mv->col };
  unsigned int best_sad = INT_MAX;
  int r, c, i;
  int start_col, end_col, start_row, end_row;
  int col_step = (step > 1) ? step : 4;

  assert(step >= 1);

  clamp_mv(&fcenter_mv, x->mv_limits.col_min, x->mv_limits.col_max,
           x->mv_limits.row_min, x->mv_limits.row_max);
  *best_mv = fcenter_mv;
  best_sad =
      fn_ptr->sdf(what->buf, what->stride,
                  get_buf_from_mv(in_what, &fcenter_mv), in_what->stride) +
      mvsad_err_cost(x, &fcenter_mv, ref_mv, mbmi->max_mv_precision,
                     sad_per_bit);
  start_row = AOMMAX(-range, x->mv_limits.row_min - fcenter_mv.row);
  start_col = AOMMAX(-range, x->mv_limits.col_min - fcenter_mv.col);
  end_row = AOMMIN(range, x->mv_limits.row_max - fcenter_mv.row);
  end_col = AOMMIN(range, x->mv_limits.col_max - fcenter_mv.col);

  for (r = start_row; r <= end_row; r += step) {
    for (c = start_col; c <= end_col; c += col_step) {
      // Step > 1 means we are not checking every location in this pass.
      if (step > 1) {
        const MV mv = { fcenter_mv.row + r, fcenter_mv.col + c };
        unsigned int sad =
            fn_ptr->sdf(what->buf, what->stride, get_buf_from_mv(in_what, &mv),
                        in_what->stride);
        if (sad < best_sad) {
          sad += mvsad_err_cost(x, &mv, ref_mv, mbmi->max_mv_precision,
                                sad_per_bit);
          if (sad < best_sad) {
            best_sad = sad;
            x->second_best_mv.as_mv = *best_mv;
            *best_mv = mv;
          }
        }
      } else {
        // 4 sads in a single call if we are checking every location
        if (c + 3 <= end_col) {
          unsigned int sads[4];
          const uint8_t *addrs[4];
          for (i = 0; i < 4; ++i) {
            const MV mv = { fcenter_mv.row + r, fcenter_mv.col + c + i };
            addrs[i] = get_buf_from_mv(in_what, &mv);
          }
          fn_ptr->sdx4df(what->buf, what->stride, addrs, in_what->stride, sads);

          for (i = 0; i < 4; ++i) {
            if (sads[i] < best_sad) {
              const MV mv = { fcenter_mv.row + r, fcenter_mv.col + c + i };
              const unsigned int sad =
                  sads[i] + mvsad_err_cost(x, &mv, ref_mv,
                                           mbmi->max_mv_precision, sad_per_bit);
              if (sad < best_sad) {
                best_sad = sad;
                x->second_best_mv.as_mv = *best_mv;
                *best_mv = mv;
              }
            }
          }
        } else {
          for (i = 0; i < end_col - c; ++i) {
            const MV mv = { fcenter_mv.row + r, fcenter_mv.col + c + i };
            unsigned int sad =
                fn_ptr->sdf(what->buf, what->stride,
                            get_buf_from_mv(in_what, &mv), in_what->stride);
            if (sad < best_sad) {
              sad += mvsad_err_cost(x, &mv, ref_mv, mbmi->max_mv_precision,
                                    sad_per_bit);
              if (sad < best_sad) {
                best_sad = sad;
                x->second_best_mv.as_mv = *best_mv;
                *best_mv = mv;
              }
            }
          }
        }
      }
    }
  }

  return best_sad;
}

/* do_refine: If last step (1-away) of n-step search doesn't pick the center
              point as the best match, we will do a final 1-away diamond
              refining search  */
static int full_pixel_diamond(const AV1_COMP *const cpi, MACROBLOCK *x,
                              MV *mvp_full, int step_param, int sadpb,
                              int further_steps, int do_refine, int *cost_list,
                              const aom_variance_fn_ptr_t *fn_ptr,
                              const MV *ref_mv, const search_site_config *cfg) {
  MV temp_mv;
  int thissme, n, num00 = 0;
  int bestsme =
      cpi->diamond_search_sad(&cpi->common, x, cfg, mvp_full, &temp_mv,
                              step_param, sadpb, &n, fn_ptr, ref_mv);
  if (bestsme < INT_MAX)
    bestsme = av1_get_mvpred_var(&cpi->common, x, &temp_mv, ref_mv, fn_ptr, 1);
  x->best_mv.as_mv = temp_mv;

  // If there won't be more n-step search, check to see if refining search is
  // needed.
  if (n > further_steps) do_refine = 0;

  while (n < further_steps) {
    ++n;

    if (num00) {
      num00--;
    } else {
      thissme = cpi->diamond_search_sad(&cpi->common, x, cfg, mvp_full,
                                        &temp_mv, step_param + n, sadpb, &num00,
                                        fn_ptr, ref_mv);
      if (thissme < INT_MAX)
        thissme =
            av1_get_mvpred_var(&cpi->common, x, &temp_mv, ref_mv, fn_ptr, 1);

      // check to see if refining search is needed.
      if (num00 > further_steps - n) do_refine = 0;

      if (thissme < bestsme) {
        bestsme = thissme;
        x->best_mv.as_mv = temp_mv;
      }
    }
  }

  // final 1-away diamond refining search
  if (do_refine) {
    const int search_range = 8;
    MV best_mv = x->best_mv.as_mv;
    thissme = av1_refining_search_sad(&cpi->common, x, &best_mv, sadpb,
                                      search_range, fn_ptr, ref_mv);
    if (thissme < INT_MAX)
      thissme =
          av1_get_mvpred_var(&cpi->common, x, &best_mv, ref_mv, fn_ptr, 1);
    if (thissme < bestsme) {
      bestsme = thissme;
      x->best_mv.as_mv = best_mv;
    }
  }

  // Return cost list.
  if (cost_list) {
    calc_int_cost_list(&cpi->common, x, ref_mv, sadpb, fn_ptr,
                       &x->best_mv.as_mv, cost_list);
  }
  return bestsme;
}

#define MIN_RANGE 7
#define MAX_RANGE 256
#define MIN_INTERVAL 1
// Runs an limited range exhaustive mesh search using a pattern set
// according to the encode speed profile.
static int full_pixel_exhaustive(const AV1_COMP *const cpi, MACROBLOCK *x,
                                 const MV *centre_mv_full, int sadpb,
                                 int *cost_list,
                                 const aom_variance_fn_ptr_t *fn_ptr,
                                 const MV *ref_mv, MV *dst_mv) {
  const SPEED_FEATURES *const sf = &cpi->sf;
  MV temp_mv = { centre_mv_full->row, centre_mv_full->col };
  MV f_ref_mv = { ref_mv->row >> 3, ref_mv->col >> 3 };
  int bestsme;
  int i;
  int interval = sf->mesh_patterns[0].interval;
  int range = sf->mesh_patterns[0].range;
  int baseline_interval_divisor;

  // Keep track of number of exhaustive calls (this frame in this thread).
  if (x->ex_search_count_ptr != NULL) ++(*x->ex_search_count_ptr);

  // Trap illegal values for interval and range for this function.
  if ((range < MIN_RANGE) || (range > MAX_RANGE) || (interval < MIN_INTERVAL) ||
      (interval > range))
    return INT_MAX;

  baseline_interval_divisor = range / interval;

  // Check size of proposed first range against magnitude of the centre
  // value used as a starting point.
  range = AOMMAX(range, (5 * AOMMAX(abs(temp_mv.row), abs(temp_mv.col))) / 4);
  range = AOMMIN(range, MAX_RANGE);
  interval = AOMMAX(interval, range / baseline_interval_divisor);

  // initial search
  bestsme = exhuastive_mesh_search(&cpi->common, x, &f_ref_mv, &temp_mv, range,
                                   interval, sadpb, fn_ptr, &temp_mv);

  if ((interval > MIN_INTERVAL) && (range > MIN_RANGE)) {
    // Progressive searches with range and step size decreasing each time
    // till we reach a step size of 1. Then break out.
    for (i = 1; i < MAX_MESH_STEP; ++i) {
      // First pass with coarser step and longer range
      bestsme = exhuastive_mesh_search(
          &cpi->common, x, &f_ref_mv, &temp_mv, sf->mesh_patterns[i].range,
          sf->mesh_patterns[i].interval, sadpb, fn_ptr, &temp_mv);

      if (sf->mesh_patterns[i].interval == 1) break;
    }
  }

  if (bestsme < INT_MAX)
    bestsme = av1_get_mvpred_var(&cpi->common, x, &temp_mv, ref_mv, fn_ptr, 1);
  *dst_mv = temp_mv;

  // Return cost list.
  if (cost_list) {
    calc_int_cost_list(&cpi->common, x, ref_mv, sadpb, fn_ptr, dst_mv,
                       cost_list);
  }
  return bestsme;
}

#define MIN_EX_SEARCH_LIMIT 128
static int is_exhaustive_allowed(const AV1_COMP *const cpi, MACROBLOCK *x) {
  const SPEED_FEATURES *const sf = &cpi->sf;
  int is_allowed = sf->allow_exhaustive_searches &&
                   (sf->exhaustive_searches_thresh < INT_MAX) &&
                   !cpi->rc.is_src_frame_alt_ref;
  if (x->m_search_count_ptr != NULL && x->ex_search_count_ptr != NULL) {
    const int max_ex =
        AOMMAX(MIN_EX_SEARCH_LIMIT,
               (*x->m_search_count_ptr * sf->max_exaustive_pct) / 100);
    is_allowed = *x->ex_search_count_ptr <= max_ex && is_allowed;
  }
  return is_allowed;
}

int av1_illum_full_pixel_search(const AV1_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE bsize, MV *mvp_full, int step_param,
                                int method, int run_mesh_search,
                                int error_per_bit, int *cost_list,
                                const MV *ref_mv, int var_max, int rd,
                                int x_pos, int y_pos, int intra,
                                const search_site_config *cfg) {
  const SPEED_FEATURES *const sf = &cpi->sf;
  const aom_variance_fn_ptr_t *fn_ptr = &cpi->fn_ptr[bsize];
  int var = 0;

  if (cost_list) {
    cost_list[0] = INT_MAX;
    cost_list[1] = INT_MAX;
    cost_list[2] = INT_MAX;
    cost_list[3] = INT_MAX;
    cost_list[4] = INT_MAX;
  }

  // Keep track of number of searches (this frame in this thread).
  if (x->m_search_count_ptr != NULL) ++(*x->m_search_count_ptr);

  switch (method) {
    case FAST_DIAMOND:
      var = fast_dia_search(cpi, x, mvp_full, step_param, error_per_bit, 0,
                            cost_list, fn_ptr, 1, ref_mv);
      break;
    case FAST_HEX:
      var = fast_hex_search(cpi, x, mvp_full, step_param, error_per_bit, 0,
                            cost_list, fn_ptr, 1, ref_mv);
      break;
    case HEX:
      var = av1_hex_search(cpi, x, mvp_full, step_param, error_per_bit, 1,
                           cost_list, fn_ptr, 1, ref_mv);
      break;
    case SQUARE:
      var = square_search(cpi, x, mvp_full, step_param, error_per_bit, 1,
                          cost_list, fn_ptr, 1, ref_mv);
      break;
    case BIGDIA:
      var = bigdia_search(cpi, x, mvp_full, step_param, error_per_bit, 1,
                          cost_list, fn_ptr, 1, ref_mv);
      break;
    case NSTEP:
      var = full_pixel_diamond(cpi, x, mvp_full, step_param, error_per_bit,
                               MAX_MVSEARCH_STEPS - 1 - step_param, 1,
                               cost_list, fn_ptr, ref_mv, cfg);

      // Should we allow a follow on exhaustive search?
      if (is_exhaustive_allowed(cpi, x)) {
        int exhuastive_thr = sf->exhaustive_searches_thresh;
        exhuastive_thr >>=
            10 - (mi_size_wide_log2[bsize] + mi_size_high_log2[bsize]);

        // Threshold variance for an exhaustive full search.
        if (var > exhuastive_thr) {
          int var_ex;
          MV tmp_mv_ex;
          var_ex =
              full_pixel_exhaustive(cpi, x, &x->best_mv.as_mv, error_per_bit,
                                    cost_list, fn_ptr, ref_mv, &tmp_mv_ex);

          if (var_ex < var) {
            var = var_ex;
            x->best_mv.as_mv = tmp_mv_ex;
          }
        }
      }
      break;
    default: assert(0 && "Invalid search method.");
  }

  // Should we allow a follow on exhaustive search?
  if (!run_mesh_search) {
    if (method == NSTEP) {
      if (is_exhaustive_allowed(cpi, x)) {
        int exhuastive_thr = sf->exhaustive_searches_thresh;
        exhuastive_thr >>=
            10 - (mi_size_wide_log2[bsize] + mi_size_high_log2[bsize]);
        // Threshold variance for an exhaustive full search.
        if (var > exhuastive_thr) run_mesh_search = 1;
      }
    }
  }

  if (run_mesh_search) {
    int var_ex;
    MV tmp_mv_ex;
    var_ex = full_pixel_exhaustive(cpi, x, &x->best_mv.as_mv, error_per_bit,
                                   cost_list, fn_ptr, ref_mv, &tmp_mv_ex);
    if (var_ex < var) {
      var = var_ex;
      x->best_mv.as_mv = tmp_mv_ex;
    }
  }

  if (method != NSTEP && rd && var < var_max)
    var = av1_get_mvpred_var(&cpi->common, x, &x->best_mv.as_mv, ref_mv, fn_ptr,
                             1);

  do {
    if (!intra || !av1_use_hash_me(&cpi->common)) break;

    // already single ME
    // get block size and original buffer of current block
    const int block_height = block_size_high[bsize];
    const int block_width = block_size_wide[bsize];
    if (block_height == block_width && x_pos >= 0 && y_pos >= 0) {
      if (block_width == 4 || block_width == 8 || block_width == 16 ||
          block_width == 32 || block_width == 64 || block_width == 128) {
        uint8_t *what = x->plane[0].src.buf;
        const int what_stride = x->plane[0].src.stride;
        uint32_t hash_value1, hash_value2;
        MV best_hash_mv;
        int best_hash_cost = INT_MAX;

        // for the hashMap
        hash_table *ref_frame_hash =
            intra ? &cpi->common.cur_frame->hash_table
                  : av1_get_ref_frame_hash_map(&cpi->common,
                                               x->e_mbd.mi[0]->ref_frame[0]);

        av1_get_block_hash_value(what, what_stride, block_width, &hash_value1,
                                 &hash_value2, is_cur_buf_hbd(&x->e_mbd), x);

        const int count = av1_hash_table_count(ref_frame_hash, hash_value1);
        // for intra, at lest one matching can be found, itself.
        if (count <= (intra ? 1 : 0)) {
          break;
        }

        Iterator iterator =
            av1_hash_get_first_iterator(ref_frame_hash, hash_value1);
        for (int i = 0; i < count; i++, aom_iterator_increment(&iterator)) {
          block_hash ref_block_hash =
              *(block_hash *)(aom_iterator_get(&iterator));
          if (hash_value2 == ref_block_hash.hash_value2) {
            // For intra, make sure the prediction is from valid area.
            if (intra) {
              const int mi_col = x_pos / MI_SIZE;
              const int mi_row = y_pos / MI_SIZE;
              const MV dv = { 8 * (ref_block_hash.y - y_pos),
                              8 * (ref_block_hash.x - x_pos) };
              if (!av1_is_dv_valid(dv, &cpi->common, &x->e_mbd, mi_row, mi_col,
                                   bsize, cpi->common.seq_params.mib_size_log2))
                continue;
            }
            MV hash_mv;
            hash_mv.col = ref_block_hash.x - x_pos;
            hash_mv.row = ref_block_hash.y - y_pos;
            if (!is_mv_in(&x->mv_limits, &hash_mv)) continue;
            const int refCost = av1_get_mvpred_var(&cpi->common, x, &hash_mv,
                                                   ref_mv, fn_ptr, 1);
            if (refCost < best_hash_cost) {
              best_hash_cost = refCost;
              best_hash_mv = hash_mv;
            }
          }
        }
        if (best_hash_cost < var) {
          x->second_best_mv = x->best_mv;
          x->best_mv.as_mv = best_hash_mv;
          var = best_hash_cost;
        }
      }
    }
  } while (0);

  return var;
}
