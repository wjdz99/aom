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
#define _POSIX_C_SOURCE 200112L  // rand_r()
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "av1/encoder/ransac.h"

#define MAX_MINPTS 16
#define MAX_DEGENERATE_ITER 10
#define MINPTS_MULTIPLIER 5

#define INLIER_THRESHOLD 1.0
#define MIN_TRIALS 20

// Number of points randomly chosen by ransac for model fitting in each
// iteration, for {translation, rotzoom, affine} and {hortrapezoid,
// vertrapezoid, homography} motions respectively.
#define NUM_INITIAL_POINTS_AFFINE_AND_SIMPLER 3
#define NUM_INITIAL_POINTS_HOMOGRAPHY_AND_TRAPEZOID 4

////////////////////////////////////////////////////////////////////////////////
// ransac
typedef int (*IsDegenerateFunc)(double *p);
typedef void (*NormalizeFunc)(double *p, int np, double *T);
typedef void (*DenormalizeFunc)(double *params, double *T1, double *T2);
typedef int (*FindTransformationFunc)(int points, double *points1,
                                      double *points2, double *params);
typedef void (*ProjectPointsDoubleFunc)(double *mat, double *points,
                                        double *proj, const int n,
                                        const int stride_points,
                                        const int stride_proj);

static void project_points_double_translation(double *mat, double *points,
                                              double *proj, const int n,
                                              const int stride_points,
                                              const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = x + mat[0];
    *(proj++) = y + mat[1];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_rotzoom(double *mat, double *points,
                                          double *proj, const int n,
                                          const int stride_points,
                                          const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = mat[2] * x + mat[3] * y + mat[0];
    *(proj++) = -mat[3] * x + mat[2] * y + mat[1];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_affine(double *mat, double *points,
                                         double *proj, const int n,
                                         const int stride_points,
                                         const int stride_proj) {
  int i;
  for (i = 0; i < n; ++i) {
    const double x = *(points++), y = *(points++);
    *(proj++) = mat[2] * x + mat[3] * y + mat[0];
    *(proj++) = mat[4] * x + mat[5] * y + mat[1];
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_hortrapezoid(double *mat, double *points,
                                               double *proj, const int n,
                                               const int stride_points,
                                               const int stride_proj) {
  int i;
  double x, y, Z, Z_inv;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    Z_inv = mat[7] * y + 1;
    assert(fabs(Z_inv) > 0.000001);
    Z = 1. / Z_inv;
    *(proj++) = (mat[2] * x + mat[3] * y + mat[0]) * Z;
    *(proj++) = (mat[5] * y + mat[1]) * Z;
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_vertrapezoid(double *mat, double *points,
                                               double *proj, const int n,
                                               const int stride_points,
                                               const int stride_proj) {
  int i;
  double x, y, Z, Z_inv;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    Z_inv = mat[6] * x + 1;
    assert(fabs(Z_inv) > 0.000001);
    Z = 1. / Z_inv;
    *(proj++) = (mat[2] * x + mat[0]) * Z;
    *(proj++) = (mat[4] * x + mat[5] * y + mat[1]) * Z;
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static void project_points_double_homography(double *mat, double *points,
                                             double *proj, const int n,
                                             const int stride_points,
                                             const int stride_proj) {
  int i;
  double x, y, Z, Z_inv;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    Z_inv = mat[6] * x + mat[7] * y + 1;
    assert(fabs(Z_inv) > 0.000001);
    Z = 1. / Z_inv;
    *(proj++) = (mat[2] * x + mat[3] * y + mat[0]) * Z;
    *(proj++) = (mat[4] * x + mat[5] * y + mat[1]) * Z;
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

static int get_rand_indices(int npoints, int minpts, int *indices,
                            unsigned int *seed) {
  int i, j;
  int ptr = rand_r(seed) % npoints;
  if (minpts > npoints) return 0;
  indices[0] = ptr;
  ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
  i = 1;
  while (i < minpts) {
    int index = rand_r(seed) % npoints;
    while (index) {
      ptr = (ptr == npoints - 1 ? 0 : ptr + 1);
      for (j = 0; j < i; ++j) {
        if (indices[j] == ptr) break;
      }
      if (j == i) index--;
    }
    indices[i++] = ptr;
  }
  return 1;
}

static int is_better_motion(int num_inliers, double variance,
                            int other_num_inliers, double other_variance) {
  if (num_inliers > other_num_inliers) return 1;
  if (num_inliers == other_num_inliers && variance < other_variance) return 1;

  return 0;
}

static void copy_points_at_indices(double *dest, const double *src,
                                   const int *indices, int num_points) {
  for (int i = 0; i < num_points; ++i) {
    const int index = indices[i];
    dest[i * 2] = src[index * 2];
    dest[i * 2 + 1] = src[index * 2 + 1];
  }
}

static int ransac(const int *matched_points, int npoints, int *number_of_inliers,
                  double *best_params, int num_desired_motions,
                  const int minpts, IsDegenerateFunc is_degenerate,
                  FindTransformationFunc find_transformation,
                  ProjectPointsDoubleFunc projectpoints) {
  static const double PROBABILITY_REQUIRED = 0.9;
  static const double EPS = 1e-12;
  static const double kInfiniteVariance = 1e12;

  int N = 10000, trial_count = 0;
  int ret_val = 0;
  //  static int rand_seed_offset = 0;
  unsigned int seed = (unsigned int)npoints;// + rand_seed_offset++;

  double params[MAX_PARAMDIM];
  WarpedMotionParams wm;
  int indices[MAX_MINPTS] = { 0 };
  int worst_kept_motion = 0;

  double *points1, *points2;
  double *corners1, *corners2;
  double *image1_coord;

  // Store the num inliers, variance, and inlier indices for the
  // num_desired_motions best transformations found so far.
  int *num_inliers_by_motion;
  double *variance_by_motion;
  int *inlier_indices_by_motion;

  // Stores the indices of the inlier points for the motion currently
  // under consideration.
  int *inlier_indices_this_motion;

  double *cnp1, *cnp2;

  *number_of_inliers = 0;
  if (npoints < minpts * MINPTS_MULTIPLIER || npoints == 0) {
    return 1;
  }

  memset(&wm, 0, sizeof(wm));

  points1 = (double *)aom_malloc(sizeof(*points1) * npoints * 2);
  points2 = (double *)aom_malloc(sizeof(*points2) * npoints * 2);
  corners1 = (double *)aom_malloc(sizeof(*corners1) * npoints * 2);
  corners2 = (double *)aom_malloc(sizeof(*corners2) * npoints * 2);
  image1_coord = (double *)aom_malloc(sizeof(*image1_coord) * npoints * 2);

  num_inliers_by_motion = (int *)aom_malloc(
      sizeof(*num_inliers_by_motion) * num_desired_motions);
  memset(num_inliers_by_motion, 0,
         sizeof(*num_inliers_by_motion) * num_desired_motions);

  variance_by_motion = (double *)aom_malloc(
      sizeof(*variance_by_motion) * num_desired_motions);
  for (int motion = 0; motion < num_desired_motions; ++motion) {
    variance_by_motion[motion] = kInfiniteVariance;
  }

  inlier_indices_by_motion = (int *)aom_malloc(
      sizeof(*inlier_indices_by_motion) * num_desired_motions * npoints);
  memset(inlier_indices_by_motion, 0,
         sizeof(*inlier_indices_by_motion) * num_desired_motions * npoints);

  inlier_indices_this_motion =
      (int *)aom_malloc(sizeof(*inlier_indices_by_motion) * npoints);
  memset(inlier_indices_by_motion, 0,
         sizeof(*inlier_indices_by_motion) * npoints);

  if (!(points1 && points2 && corners1 && corners2 && image1_coord &&
        num_inliers_by_motion && variance_by_motion &&
        inlier_indices_by_motion && inlier_indices_this_motion)) {
    ret_val = 1;
    goto finish_ransac;
  }

  cnp1 = corners1;
  cnp2 = corners2;
  for (int i = 0; i < npoints; ++i) {
    *(cnp1++) = *(matched_points++);
    *(cnp1++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
    *(cnp2++) = *(matched_points++);
  }
  matched_points -= 4 * npoints;

  while (N > trial_count) {
    int num_inliers = 0;
    double sum_distance = 0.0;
    double sum_distance_squared = 0.0;

    int degenerate = 1;
    int num_degenerate_iter = 0;
    while (degenerate) {
      num_degenerate_iter++;
      if (!get_rand_indices(npoints, minpts, indices, &seed)) {
        ret_val = 1;
        goto finish_ransac;
      }

      copy_points_at_indices(points1, corners1, indices, minpts);
      copy_points_at_indices(points2, corners2, indices, minpts);

      degenerate = is_degenerate(points1);
      if (num_degenerate_iter > MAX_DEGENERATE_ITER) {
        ret_val = 1;
        goto finish_ransac;
      }
    }

    if (find_transformation(minpts, points1, points2, params)) {
      trial_count++;
      continue;
    }

    projectpoints(params, corners1, image1_coord, npoints, 2, 2);

    for (int index = 0; index < npoints; ++index) {
      double dx = image1_coord[index * 2] - corners2[index * 2];
      double dy = image1_coord[index * 2 + 1] - corners2[index * 2 + 1];
      double distance = sqrt(dx * dx + dy * dy);

      if (distance < INLIER_THRESHOLD) {
        inlier_indices_this_motion[num_inliers++] = index;
        sum_distance += distance;
        sum_distance_squared += distance * distance;
      }
    }

    if (num_inliers >= num_inliers_by_motion[worst_kept_motion]) {
      int temp;
      double fracinliers, pNoOutliers, mean_distance, variance;
      mean_distance = sum_distance / ((double)num_inliers);
      variance = sum_distance_squared / ((double)num_inliers - 1.0) -
                 mean_distance * mean_distance * ((double)num_inliers) /
                     ((double)num_inliers - 1.0);
      if (is_better_motion(num_inliers, variance,
                           num_inliers_by_motion[worst_kept_motion],
                           variance_by_motion[worst_kept_motion])) {
        // This is one of the top num_desired_motions transformations
        // found so far. Remember the inlier points and the
        // variance. The parameters for each motion will be recomputed at
        // exit using only the inliers.
        num_inliers_by_motion[worst_kept_motion] = num_inliers;
        variance_by_motion[worst_kept_motion] = variance;
        memcpy(inlier_indices_by_motion + worst_kept_motion * npoints,
               inlier_indices_this_motion,
               sizeof(*inlier_indices_by_motion) * npoints);

        assert(npoints > 0);
        fracinliers = (double)num_inliers / (double)npoints;
        pNoOutliers = 1 - pow(fracinliers, minpts);
        pNoOutliers = fmax(EPS, pNoOutliers);
        pNoOutliers = fmin(1 - EPS, pNoOutliers);
        temp = (int)(log(1.0 - PROBABILITY_REQUIRED) / log(pNoOutliers));
        if (temp > 0 && temp < N) {
          N = AOMMAX(temp, MIN_TRIALS);
        }

        // Determine the new worst kept motion and its num_inliers and variance.
        for (int motion = 0; motion < num_desired_motions; ++motion) {
          if (is_better_motion(num_inliers_by_motion[worst_kept_motion],
                               variance_by_motion[worst_kept_motion],
                               num_inliers_by_motion[motion],
                               variance_by_motion[motion])) {
            worst_kept_motion = motion;
          }
        }
      }
    }
    trial_count++;
  }

  // Recompute the motions using only the inliers.
  for (int motion = 0; motion < num_desired_motions; ++motion) {
    const int num_inliers = num_inliers_by_motion[motion];
    const int *inlier_indices = inlier_indices_by_motion + motion * npoints;

    copy_points_at_indices(points1, corners1, inlier_indices, num_inliers);
    copy_points_at_indices(points2, corners2, inlier_indices, num_inliers);

    find_transformation(num_inliers, points1, points2,
                        best_params + (MAX_PARAMDIM - 1) * motion);
  }

finish_ransac:
  aom_free(points1);
  aom_free(points2);
  aom_free(corners1);
  aom_free(corners2);
  aom_free(image1_coord);
  aom_free(num_inliers_by_motion);
  aom_free(variance_by_motion);
  aom_free(inlier_indices_by_motion);
  aom_free(inlier_indices_this_motion);

  return ret_val;
}

static int is_collinear3(double *p1, double *p2, double *p3) {
  static const double collinear_eps = 1e-3;
  const double v =
      (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
  return fabs(v) < collinear_eps;
}

static int is_degenerate_translation(double *p) {
  return (p[0] - p[2]) * (p[0] - p[2]) + (p[1] - p[3]) * (p[1] - p[3]) <= 2;
}

static int is_degenerate_affine(double *p) {
  return is_collinear3(p, p + 2, p + 4);
}

static int is_degenerate_homography(double *p) {
  return is_collinear3(p, p + 2, p + 4) || is_collinear3(p, p + 2, p + 6) ||
         is_collinear3(p, p + 4, p + 6) || is_collinear3(p + 2, p + 4, p + 6);
}

int ransac_translation(int *matched_points, int npoints, int *number_of_inliers,
                       double *best_params, int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions, NUM_INITIAL_POINTS_AFFINE_AND_SIMPLER,
                is_degenerate_translation, find_translation,
                project_points_double_translation);
}

int ransac_rotzoom(int *matched_points, int npoints, int *number_of_inliers,
                   double *best_params, int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions, NUM_INITIAL_POINTS_AFFINE_AND_SIMPLER,
                is_degenerate_affine, find_rotzoom,
                project_points_double_rotzoom);
}

int ransac_affine(int *matched_points, int npoints, int *number_of_inliers,
                  double *best_params, int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions, NUM_INITIAL_POINTS_AFFINE_AND_SIMPLER,
                is_degenerate_affine, find_affine,
                project_points_double_affine);
}

int ransac_homography(int *matched_points, int npoints, int *number_of_inliers,
                      double *best_params, int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions,
                NUM_INITIAL_POINTS_HOMOGRAPHY_AND_TRAPEZOID,
                is_degenerate_homography, find_homography,
                project_points_double_homography);
}

int ransac_hortrapezoid(int *matched_points, int npoints,
                        int *number_of_inliers, double *best_params,
                        int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions,
                NUM_INITIAL_POINTS_HOMOGRAPHY_AND_TRAPEZOID,
                is_degenerate_homography, find_hortrapezoid,
                project_points_double_hortrapezoid);
}

int ransac_vertrapezoid(int *matched_points, int npoints,
                        int *number_of_inliers, double *best_params,
                        int num_desired_motions) {
  return ransac(matched_points, npoints, number_of_inliers, best_params,
                num_desired_motions,
                NUM_INITIAL_POINTS_HOMOGRAPHY_AND_TRAPEZOID,
                is_degenerate_homography, find_vertrapezoid,
                project_points_double_vertrapezoid);
}
