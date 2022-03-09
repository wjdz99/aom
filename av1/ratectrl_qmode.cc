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

#include "av1/ratectrl_qmode.h"

#include <vector>

namespace aom {

#define AOMMIN(x, y) (((x) < (y)) ? (x) : (y))
#define AOMMAX(x, y) (((x) > (y)) ? (x) : (y))

#define WINDOW_SIZE 7
#define HALF_WIN (WINDOW_SIZE / 2)

#define MIN_SHRINK_LEN 6  // the minimum length of gf if we are shrinking
#define SMOOTH_FILT_LEN 7
#define HALF_FILT_LEN (SMOOTH_FILT_LEN / 2)

// A 7-tap gaussian smooth filter
const double smooth_filt[SMOOTH_FILT_LEN] = { 0.006, 0.061, 0.242, 0.383,
                                              0.242, 0.061, 0.006 };

static int find_next_scenecut(const FIRSTPASS_STATS *const stats_start,
                              int first, int last) {
  // Identify unstable areas caused by scenecuts.
  // Find the max and 2nd max coded error, and the average of the rest frames.
  // If there is only one frame that yields a huge coded error, it is likely a
  // scenecut.
  double this_ratio, max_prev_ratio, max_next_ratio, max_prev_coded,
      max_next_coded;

  if (last - first == 0) return -1;

  for (int i = first; i <= last; i++) {
    if (stats_start[i].is_flash || (i > 0 && stats_start[i - 1].is_flash))
      continue;
    double temp_intra = AOMMAX(stats_start[i].intra_error, 0.01);
    this_ratio = stats_start[i].coded_error / temp_intra;
    // find the avg ratio in the preceding neighborhood
    max_prev_ratio = 0;
    max_prev_coded = 0;
    for (int j = AOMMAX(first, i - HALF_WIN); j < i; j++) {
      if (stats_start[j].is_flash || (j > 0 && stats_start[j - 1].is_flash))
        continue;
      temp_intra = AOMMAX(stats_start[j].intra_error, 0.01);
      double temp_ratio = stats_start[j].coded_error / temp_intra;
      if (temp_ratio > max_prev_ratio) {
        max_prev_ratio = temp_ratio;
      }
      if (stats_start[j].coded_error > max_prev_coded) {
        max_prev_coded = stats_start[j].coded_error;
      }
    }
    // find the avg ratio in the following neighborhood
    max_next_ratio = 0;
    max_next_coded = 0;
    for (int j = i + 1; j <= AOMMIN(i + HALF_WIN, last); j++) {
      if (stats_start[i].is_flash || (i > 0 && stats_start[i - 1].is_flash))
        continue;
      temp_intra = AOMMAX(stats_start[j].intra_error, 0.01);
      double temp_ratio = stats_start[j].coded_error / temp_intra;
      if (temp_ratio > max_next_ratio) {
        max_next_ratio = temp_ratio;
      }
      if (stats_start[j].coded_error > max_next_coded) {
        max_next_coded = stats_start[j].coded_error;
      }
    }

    if (max_prev_ratio < 0.001 && max_next_ratio < 0.001) {
      // the ratios are very small, only check a small fixed threshold
      if (this_ratio < 0.02) continue;
    } else {
      // check if this frame has a larger ratio than the neighborhood
      double max_sr = stats_start[i].sr_coded_error;
      if (i < last) max_sr = AOMMAX(max_sr, stats_start[i + 1].sr_coded_error);
      double max_sr_fr_ratio =
          max_sr / AOMMAX(stats_start[i].coded_error, 0.01);

      if (max_sr_fr_ratio > 1.2) continue;
      if (this_ratio < 2 * AOMMAX(max_prev_ratio, max_next_ratio) &&
          stats_start[i].coded_error <
              2 * AOMMAX(max_prev_coded, max_next_coded)) {
        continue;
      }
    }
    return i;
  }
  return -1;
}

// Get the average of stats inside a region.
static void analyze_region(const FIRSTPASS_STATS *stats, int k,
                           REGIONS *regions) {
  int i;
  regions[k].avg_cor_coeff = 0;
  regions[k].avg_sr_fr_ratio = 0;
  regions[k].avg_intra_err = 0;
  regions[k].avg_coded_err = 0;

  int check_first_sr = (k != 0);

  for (i = regions[k].start; i <= regions[k].last; i++) {
    if (i > regions[k].start || check_first_sr) {
      double num_frames =
          (double)(regions[k].last - regions[k].start + check_first_sr);
      double max_coded_error =
          AOMMAX(stats[i].coded_error, stats[i - 1].coded_error);
      double this_ratio =
          stats[i].sr_coded_error / AOMMAX(max_coded_error, 0.001);
      regions[k].avg_sr_fr_ratio += this_ratio / num_frames;
    }

    regions[k].avg_intra_err +=
        stats[i].intra_error / (double)(regions[k].last - regions[k].start + 1);
    regions[k].avg_coded_err +=
        stats[i].coded_error / (double)(regions[k].last - regions[k].start + 1);

    regions[k].avg_cor_coeff +=
        AOMMAX(stats[i].cor_coeff, 0.001) /
        (double)(regions[k].last - regions[k].start + 1);
    regions[k].avg_noise_var +=
        AOMMAX(stats[i].noise_var, 0.001) /
        (double)(regions[k].last - regions[k].start + 1);
  }
}

// Smooth filter intra_error and coded_error in firstpass stats.
// If stats[i].is_flash==1, the ith element should not be used in the filtering.
static void smooth_filter_stats(const FIRSTPASS_STATS *stats, int start_idx,
                                int last_idx, double *filt_intra_err,
                                double *filt_coded_err) {
  int i, j;
  for (i = start_idx; i <= last_idx; i++) {
    double total_wt = 0;
    for (j = -HALF_FILT_LEN; j <= HALF_FILT_LEN; j++) {
      int idx = AOMMIN(AOMMAX(i + j, start_idx), last_idx);
      if (stats[idx].is_flash) continue;

      filt_intra_err[i] +=
          smooth_filt[j + HALF_FILT_LEN] * stats[idx].intra_error;
      total_wt += smooth_filt[j + HALF_FILT_LEN];
    }
    if (total_wt > 0.01) {
      filt_intra_err[i] /= total_wt;
    } else {
      filt_intra_err[i] = stats[i].intra_error;
    }
  }
  for (i = start_idx; i <= last_idx; i++) {
    double total_wt = 0;
    for (j = -HALF_FILT_LEN; j <= HALF_FILT_LEN; j++) {
      int idx = AOMMIN(AOMMAX(i + j, start_idx), last_idx);
      // Coded error involves idx and idx - 1.
      if (stats[idx].is_flash || (idx > 0 && stats[idx - 1].is_flash)) continue;

      filt_coded_err[i] +=
          smooth_filt[j + HALF_FILT_LEN] * stats[idx].coded_error;
      total_wt += smooth_filt[j + HALF_FILT_LEN];
    }
    if (total_wt > 0.01) {
      filt_coded_err[i] /= total_wt;
    } else {
      filt_coded_err[i] = stats[i].coded_error;
    }
  }
}

// Calculate gradient
static void get_gradient(const double *values, int start, int last,
                         double *grad) {
  if (start == last) {
    grad[start] = 0;
    return;
  }
  for (int i = start; i <= last; i++) {
    int prev = AOMMAX(i - 1, start);
    int next = AOMMIN(i + 1, last);
    grad[i] = (values[next] - values[prev]) / (next - prev);
  }
}

// Calculate the regions stats of every region.
static void get_region_stats(const FIRSTPASS_STATS *stats, REGIONS *regions,
                             int num_regions) {
  for (int k = 0; k < num_regions; k++) {
    analyze_region(stats, k, regions);
  }
}

// Find tentative stable regions
static int find_stable_regions(const FIRSTPASS_STATS *stats,
                               const double *grad_coded, int this_start,
                               int this_last, REGIONS *regions) {
  int i, j, k = 0;
  regions[k].start = this_start;
  for (i = this_start; i <= this_last; i++) {
    // Check mean and variance of stats in a window
    double mean_intra = 0.001, var_intra = 0.001;
    double mean_coded = 0.001, var_coded = 0.001;
    int count = 0;
    for (j = -HALF_WIN; j <= HALF_WIN; j++) {
      int idx = AOMMIN(AOMMAX(i + j, this_start), this_last);
      if (stats[idx].is_flash || (idx > 0 && stats[idx - 1].is_flash)) continue;
      mean_intra += stats[idx].intra_error;
      var_intra += stats[idx].intra_error * stats[idx].intra_error;
      mean_coded += stats[idx].coded_error;
      var_coded += stats[idx].coded_error * stats[idx].coded_error;
      count++;
    }

    REGION_TYPES cur_type;
    if (count > 0) {
      mean_intra /= (double)count;
      var_intra /= (double)count;
      mean_coded /= (double)count;
      var_coded /= (double)count;
      int is_intra_stable = (var_intra / (mean_intra * mean_intra) < 1.03);
      int is_coded_stable = (var_coded / (mean_coded * mean_coded) < 1.04 &&
                             fabs(grad_coded[i]) / mean_coded < 0.05) ||
                            mean_coded / mean_intra < 0.05;
      int is_coded_small = mean_coded < 0.5 * mean_intra;
      cur_type = (is_intra_stable && is_coded_stable && is_coded_small)
                     ? STABLE_REGION
                     : HIGH_VAR_REGION;
    } else {
      cur_type = HIGH_VAR_REGION;
    }

    // mark a new region if type changes
    if (i == regions[k].start) {
      // first frame in the region
      regions[k].type = cur_type;
    } else if (cur_type != regions[k].type) {
      // Append a new region
      regions[k].last = i - 1;
      regions[k + 1].start = i;
      regions[k + 1].type = cur_type;
      k++;
    }
  }
  regions[k].last = this_last;
  return k + 1;
}

// Insert a region in the cur_region_idx. The start and last should both be in
// the current region. After insertion, the cur_region_idx will point to the
// last region that was splitted from the original region.
static void insert_region(int start, int last, REGION_TYPES type,
                          REGIONS *regions, int *num_regions,
                          int *cur_region_idx) {
  int k = *cur_region_idx;
  REGION_TYPES this_region_type = regions[k].type;
  int this_region_last = regions[k].last;
  int num_add = (start != regions[k].start) + (last != regions[k].last);
  // move the following regions further to the back
  for (int r = *num_regions - 1; r > k; r--) {
    regions[r + num_add] = regions[r];
  }
  *num_regions += num_add;
  if (start > regions[k].start) {
    regions[k].last = start - 1;
    k++;
    regions[k].start = start;
  }
  regions[k].type = type;
  if (last < this_region_last) {
    regions[k].last = last;
    k++;
    regions[k].start = last + 1;
    regions[k].last = this_region_last;
    regions[k].type = this_region_type;
  } else {
    regions[k].last = this_region_last;
  }
  *cur_region_idx = k;
}

// Remove the region with index next_region.
// parameter merge: 0: merge with previous; 1: merge with next; 2:
// merge with both, take type from previous if possible
// After removing, next_region will be the index of the next region.
static void remove_region(int merge, REGIONS *regions, int *num_regions,
                          int *next_region) {
  int k = *next_region;
  assert(k < *num_regions);
  if (*num_regions == 1) {
    *num_regions = 0;
    return;
  }
  if (k == 0) {
    merge = 1;
  } else if (k == *num_regions - 1) {
    merge = 0;
  }
  int num_merge = (merge == 2) ? 2 : 1;
  switch (merge) {
    case 0:
      regions[k - 1].last = regions[k].last;
      *next_region = k;
      break;
    case 1:
      regions[k + 1].start = regions[k].start;
      *next_region = k + 1;
      break;
    case 2:
      regions[k - 1].last = regions[k + 1].last;
      *next_region = k;
      break;
    default: assert(0);
  }
  *num_regions -= num_merge;
  for (k = *next_region - (merge == 1); k < *num_regions; k++) {
    regions[k] = regions[k + num_merge];
  }
}

// Clean up regions that should be removed or merged.
static void cleanup_regions(REGIONS *regions, int *num_regions) {
  int k = 0;
  while (k < *num_regions) {
    if ((k > 0 && regions[k - 1].type == regions[k].type &&
         regions[k].type != SCENECUT_REGION) ||
        regions[k].last < regions[k].start) {
      remove_region(0, regions, num_regions, &k);
    } else {
      k++;
    }
  }
}

// Remove regions that are of type and shorter than length.
// Merge it with its neighboring regions.
static void remove_short_regions(REGIONS *regions, int *num_regions,
                                 REGION_TYPES type, int length) {
  int k = 0;
  while (k < *num_regions && (*num_regions) > 1) {
    if ((regions[k].last - regions[k].start + 1 < length &&
         regions[k].type == type)) {
      // merge current region with the previous and next regions
      remove_region(2, regions, num_regions, &k);
    } else {
      k++;
    }
  }
  cleanup_regions(regions, num_regions);
}

static void adjust_unstable_region_bounds(const FIRSTPASS_STATS *stats,
                                          REGIONS *regions, int *num_regions) {
  int i, j, k;
  // Remove regions that are too short. Likely noise.
  remove_short_regions(regions, num_regions, STABLE_REGION, HALF_WIN);
  remove_short_regions(regions, num_regions, HIGH_VAR_REGION, HALF_WIN);

  get_region_stats(stats, regions, *num_regions);

  // Adjust region boundaries. The thresholds are empirically obtained, but
  // overall the performance is not very sensitive to small changes to them.
  for (k = 0; k < *num_regions; k++) {
    if (regions[k].type == STABLE_REGION) continue;
    if (k > 0) {
      // Adjust previous boundary.
      // First find the average intra/coded error in the previous
      // neighborhood.
      double avg_intra_err = 0;
      const int starti = AOMMAX(regions[k - 1].last - WINDOW_SIZE + 1,
                                regions[k - 1].start + 1);
      const int lasti = regions[k - 1].last;
      int counti = 0;
      for (i = starti; i <= lasti; i++) {
        avg_intra_err += stats[i].intra_error;
        counti++;
      }
      if (counti > 0) {
        avg_intra_err = AOMMAX(avg_intra_err / (double)counti, 0.001);
        int count_coded = 0, count_grad = 0;
        for (j = lasti + 1; j <= regions[k].last; j++) {
          const int intra_close =
              fabs(stats[j].intra_error - avg_intra_err) / avg_intra_err < 0.1;
          const int coded_small = stats[j].coded_error / avg_intra_err < 0.1;
          const int coeff_close = stats[j].cor_coeff > 0.995;
          if (!coeff_close || !coded_small) count_coded--;
          if (intra_close && count_coded >= 0 && count_grad >= 0) {
            // this frame probably belongs to the previous stable region
            regions[k - 1].last = j;
            regions[k].start = j + 1;
          } else {
            break;
          }
        }
      }
    }  // if k > 0
    if (k < *num_regions - 1) {
      // Adjust next boundary.
      // First find the average intra/coded error in the next neighborhood.
      double avg_intra_err = 0;
      const int starti = regions[k + 1].start;
      const int lasti = AOMMIN(regions[k + 1].last - 1,
                               regions[k + 1].start + WINDOW_SIZE - 1);
      int counti = 0;
      for (i = starti; i <= lasti; i++) {
        avg_intra_err += stats[i].intra_error;
        counti++;
      }
      if (counti > 0) {
        avg_intra_err = AOMMAX(avg_intra_err / (double)counti, 0.001);
        // At the boundary, coded error is large, but still the frame is stable
        int count_coded = 1, count_grad = 1;
        for (j = starti - 1; j >= regions[k].start; j--) {
          const int intra_close =
              fabs(stats[j].intra_error - avg_intra_err) / avg_intra_err < 0.1;
          const int coded_small =
              stats[j + 1].coded_error / avg_intra_err < 0.1;
          const int coeff_close = stats[j].cor_coeff > 0.995;
          if (!coeff_close || !coded_small) count_coded--;
          if (intra_close && count_coded >= 0 && count_grad >= 0) {
            // this frame probably belongs to the next stable region
            regions[k + 1].start = j;
            regions[k].last = j - 1;
          } else {
            break;
          }
        }
      }
    }  // if k < *num_regions - 1
  }    // end of loop over all regions

  cleanup_regions(regions, num_regions);
  remove_short_regions(regions, num_regions, HIGH_VAR_REGION, HALF_WIN);
  get_region_stats(stats, regions, *num_regions);

  // If a stable regions has higher error than neighboring high var regions,
  // or if the stable region has a lower average correlation,
  // then it should be merged with them
  k = 0;
  while (k < *num_regions && (*num_regions) > 1) {
    if (regions[k].type == STABLE_REGION &&
        (regions[k].last - regions[k].start + 1) < 2 * WINDOW_SIZE &&
        ((k > 0 &&  // previous regions
          (regions[k].avg_coded_err > regions[k - 1].avg_coded_err * 1.01 ||
           regions[k].avg_cor_coeff < regions[k - 1].avg_cor_coeff * 0.999)) &&
         (k < *num_regions - 1 &&  // next region
          (regions[k].avg_coded_err > regions[k + 1].avg_coded_err * 1.01 ||
           regions[k].avg_cor_coeff < regions[k + 1].avg_cor_coeff * 0.999)))) {
      // merge current region with the previous and next regions
      remove_region(2, regions, num_regions, &k);
      analyze_region(stats, k - 1, regions);
    } else if (regions[k].type == HIGH_VAR_REGION &&
               (regions[k].last - regions[k].start + 1) < 2 * WINDOW_SIZE &&
               ((k > 0 &&  // previous regions
                 (regions[k].avg_coded_err <
                      regions[k - 1].avg_coded_err * 0.99 ||
                  regions[k].avg_cor_coeff >
                      regions[k - 1].avg_cor_coeff * 1.001)) &&
                (k < *num_regions - 1 &&  // next region
                 (regions[k].avg_coded_err <
                      regions[k + 1].avg_coded_err * 0.99 ||
                  regions[k].avg_cor_coeff >
                      regions[k + 1].avg_cor_coeff * 1.001)))) {
      // merge current region with the previous and next regions
      remove_region(2, regions, num_regions, &k);
      analyze_region(stats, k - 1, regions);
    } else {
      k++;
    }
  }

  remove_short_regions(regions, num_regions, STABLE_REGION, WINDOW_SIZE);
  remove_short_regions(regions, num_regions, HIGH_VAR_REGION, HALF_WIN);
}

// Identify blending regions.
static void find_blending_regions(const FIRSTPASS_STATS *stats,
                                  REGIONS *regions, int *num_regions) {
  int i, k = 0;
  // Blending regions will have large content change, therefore will have a
  // large consistent change in intra error.
  int count_stable = 0;
  while (k < *num_regions) {
    if (regions[k].type == STABLE_REGION) {
      k++;
      count_stable++;
      continue;
    }
    int dir = 0;
    int start = 0, last;
    for (i = regions[k].start; i <= regions[k].last; i++) {
      // First mark the regions that has consistent large change of intra error.
      if (k == 0 && i == regions[k].start) continue;
      if (stats[i].is_flash || (i > 0 && stats[i - 1].is_flash)) continue;
      double grad = stats[i].intra_error - stats[i - 1].intra_error;
      int large_change = fabs(grad) / AOMMAX(stats[i].intra_error, 0.01) > 0.05;
      int this_dir = 0;
      if (large_change) {
        this_dir = (grad > 0) ? 1 : -1;
      }
      // the current trend continues
      if (dir == this_dir) continue;
      if (dir != 0) {
        // Mark the end of a new large change group and add it
        last = i - 1;
        insert_region(start, last, BLENDING_REGION, regions, num_regions, &k);
      }
      dir = this_dir;
      if (k == 0 && i == regions[k].start + 1) {
        start = i - 1;
      } else {
        start = i;
      }
    }
    if (dir != 0) {
      last = regions[k].last;
      insert_region(start, last, BLENDING_REGION, regions, num_regions, &k);
    }
    k++;
  }

  // If the blending region has very low correlation, mark it as high variance
  // since we probably cannot benefit from it anyways.
  get_region_stats(stats, regions, *num_regions);
  for (k = 0; k < *num_regions; k++) {
    if (regions[k].type != BLENDING_REGION) continue;
    if (regions[k].last == regions[k].start || regions[k].avg_cor_coeff < 0.6 ||
        count_stable == 0)
      regions[k].type = HIGH_VAR_REGION;
  }
  get_region_stats(stats, regions, *num_regions);

  // It is possible for blending to result in a "dip" in intra error (first
  // decrease then increase). Therefore we need to find the dip and combine the
  // two regions.
  k = 1;
  while (k < *num_regions) {
    if (k < *num_regions - 1 && regions[k].type == HIGH_VAR_REGION) {
      // Check if this short high variance regions is actually in the middle of
      // a blending region.
      if (regions[k - 1].type == BLENDING_REGION &&
          regions[k + 1].type == BLENDING_REGION &&
          regions[k].last - regions[k].start < 3) {
        int prev_dir = (stats[regions[k - 1].last].intra_error -
                        stats[regions[k - 1].last - 1].intra_error) > 0
                           ? 1
                           : -1;
        int next_dir = (stats[regions[k + 1].last].intra_error -
                        stats[regions[k + 1].last - 1].intra_error) > 0
                           ? 1
                           : -1;
        if (prev_dir < 0 && next_dir > 0) {
          // This is possibly a mid region of blending. Check the ratios
          double ratio_thres = AOMMIN(regions[k - 1].avg_sr_fr_ratio,
                                      regions[k + 1].avg_sr_fr_ratio) *
                               0.95;
          if (regions[k].avg_sr_fr_ratio > ratio_thres) {
            regions[k].type = BLENDING_REGION;
            remove_region(2, regions, num_regions, &k);
            analyze_region(stats, k - 1, regions);
            continue;
          }
        }
      }
    }
    // Check if we have a pair of consecutive blending regions.
    if (regions[k - 1].type == BLENDING_REGION &&
        regions[k].type == BLENDING_REGION) {
      int prev_dir = (stats[regions[k - 1].last].intra_error -
                      stats[regions[k - 1].last - 1].intra_error) > 0
                         ? 1
                         : -1;
      int next_dir = (stats[regions[k].last].intra_error -
                      stats[regions[k].last - 1].intra_error) > 0
                         ? 1
                         : -1;

      // if both are too short, no need to check
      int total_length = regions[k].last - regions[k - 1].start + 1;
      if (total_length < 4) {
        regions[k - 1].type = HIGH_VAR_REGION;
        k++;
        continue;
      }

      int to_merge = 0;
      if (prev_dir < 0 && next_dir > 0) {
        // In this case we check the last frame in the previous region.
        double prev_length =
            (double)(regions[k - 1].last - regions[k - 1].start + 1);
        double last_ratio, ratio_thres;
        if (prev_length < 2.01) {
          // if the previous region is very short
          double max_coded_error =
              AOMMAX(stats[regions[k - 1].last].coded_error,
                     stats[regions[k - 1].last - 1].coded_error);
          last_ratio = stats[regions[k - 1].last].sr_coded_error /
                       AOMMAX(max_coded_error, 0.001);
          ratio_thres = regions[k].avg_sr_fr_ratio * 0.95;
        } else {
          double max_coded_error =
              AOMMAX(stats[regions[k - 1].last].coded_error,
                     stats[regions[k - 1].last - 1].coded_error);
          last_ratio = stats[regions[k - 1].last].sr_coded_error /
                       AOMMAX(max_coded_error, 0.001);
          double prev_ratio =
              (regions[k - 1].avg_sr_fr_ratio * prev_length - last_ratio) /
              (prev_length - 1.0);
          ratio_thres = AOMMIN(prev_ratio, regions[k].avg_sr_fr_ratio) * 0.95;
        }
        if (last_ratio > ratio_thres) {
          to_merge = 1;
        }
      }

      if (to_merge) {
        remove_region(0, regions, num_regions, &k);
        analyze_region(stats, k - 1, regions);
        continue;
      } else {
        // These are possibly two separate blending regions. Mark the boundary
        // frame as HIGH_VAR_REGION to separate the two.
        int prev_k = k - 1;
        insert_region(regions[prev_k].last, regions[prev_k].last,
                      HIGH_VAR_REGION, regions, num_regions, &prev_k);
        analyze_region(stats, prev_k, regions);
        k = prev_k + 1;
        analyze_region(stats, k, regions);
      }
    }
    k++;
  }
  cleanup_regions(regions, num_regions);
}

// Clean up decision for blendings. Remove blending regions that are too short.
// Also if a very short high var region is between a blending and a stable
// region, just merge it with one of them.
static void cleanup_blendings(REGIONS *regions, int *num_regions) {
  int k = 0;
  while (k<*num_regions && * num_regions> 1) {
    int is_short_blending = regions[k].type == BLENDING_REGION &&
                            regions[k].last - regions[k].start + 1 < 5;
    int is_short_hv = regions[k].type == HIGH_VAR_REGION &&
                      regions[k].last - regions[k].start + 1 < 5;
    int has_stable_neighbor =
        ((k > 0 && regions[k - 1].type == STABLE_REGION) ||
         (k < *num_regions - 1 && regions[k + 1].type == STABLE_REGION));
    int has_blend_neighbor =
        ((k > 0 && regions[k - 1].type == BLENDING_REGION) ||
         (k < *num_regions - 1 && regions[k + 1].type == BLENDING_REGION));
    int total_neighbors = (k > 0) + (k < *num_regions - 1);

    if (is_short_blending ||
        (is_short_hv &&
         has_stable_neighbor + has_blend_neighbor >= total_neighbors)) {
      // Remove this region.Try to determine whether to combine it with the
      // previous or next region.
      int merge;
      double prev_diff =
          (k > 0)
              ? fabs(regions[k].avg_cor_coeff - regions[k - 1].avg_cor_coeff)
              : 1;
      double next_diff =
          (k < *num_regions - 1)
              ? fabs(regions[k].avg_cor_coeff - regions[k + 1].avg_cor_coeff)
              : 1;
      // merge == 0 means to merge with previous, 1 means to merge with next
      merge = prev_diff > next_diff;
      remove_region(merge, regions, num_regions, &k);
    } else {
      k++;
    }
  }
  cleanup_regions(regions, num_regions);
}

// Identify stable and unstable regions from first pass stats.
// Stats_start points to the first frame to analyze.
// Offset is the offset from the current frame to the frame stats_start is
// pointing to.
static void identify_regions(const FIRSTPASS_STATS *const stats_start,
                             int total_frames, int offset, REGIONS *regions,
                             int *total_regions) {
  int k;
  if (total_frames <= 1) return;

  // store the initial decisions
  REGIONS temp_regions[MAX_FIRSTPASS_ANALYSIS_FRAMES];
  av1_zero_array(temp_regions, MAX_FIRSTPASS_ANALYSIS_FRAMES);
  // buffers for filtered stats
  double filt_intra_err[MAX_FIRSTPASS_ANALYSIS_FRAMES] = { 0 };
  double filt_coded_err[MAX_FIRSTPASS_ANALYSIS_FRAMES] = { 0 };
  double grad_coded[MAX_FIRSTPASS_ANALYSIS_FRAMES] = { 0 };

  int cur_region = 0, this_start = 0, this_last;

  int next_scenecut = -1;
  do {
    // first get the obvious scenecuts
    next_scenecut =
        find_next_scenecut(stats_start, this_start, total_frames - 1);
    this_last = (next_scenecut >= 0) ? (next_scenecut - 1) : total_frames - 1;

    // low-pass filter the needed stats
    smooth_filter_stats(stats_start, this_start, this_last, filt_intra_err,
                        filt_coded_err);
    get_gradient(filt_coded_err, this_start, this_last, grad_coded);

    // find tentative stable regions and unstable regions
    int num_regions = find_stable_regions(stats_start, grad_coded, this_start,
                                          this_last, temp_regions);

    adjust_unstable_region_bounds(stats_start, temp_regions, &num_regions);

    get_region_stats(stats_start, temp_regions, num_regions);

    // Try to identify blending regions in the unstable regions
    find_blending_regions(stats_start, temp_regions, &num_regions);
    cleanup_blendings(temp_regions, &num_regions);

    // The flash points should all be considered high variance points
    k = 0;
    while (k < num_regions) {
      if (temp_regions[k].type != STABLE_REGION) {
        k++;
        continue;
      }
      int start = temp_regions[k].start;
      int last = temp_regions[k].last;
      for (int i = start; i <= last; i++) {
        if (stats_start[i].is_flash) {
          insert_region(i, i, HIGH_VAR_REGION, temp_regions, &num_regions, &k);
        }
      }
      k++;
    }
    cleanup_regions(temp_regions, &num_regions);

    // copy the regions in the scenecut group
    for (k = 0; k < num_regions; k++) {
      if (temp_regions[k].last < temp_regions[k].start &&
          k == num_regions - 1) {
        num_regions--;
        break;
      }
      regions[k + cur_region] = temp_regions[k];
    }
    cur_region += num_regions;

    // add the scenecut region
    if (next_scenecut > -1) {
      // add the scenecut region, and find the next scenecut
      regions[cur_region].type = SCENECUT_REGION;
      regions[cur_region].start = next_scenecut;
      regions[cur_region].last = next_scenecut;
      cur_region++;
      this_start = next_scenecut + 1;
    }
  } while (next_scenecut >= 0);

  *total_regions = cur_region;
  get_region_stats(stats_start, regions, *total_regions);

  for (k = 0; k < *total_regions; k++) {
    // If scenecuts are very minor, mark them as high variance.
    if (regions[k].type != SCENECUT_REGION ||
        regions[k].avg_cor_coeff *
                (1 - stats_start[regions[k].start].noise_var /
                         regions[k].avg_intra_err) <
            0.8) {
      continue;
    }
    regions[k].type = HIGH_VAR_REGION;
  }
  cleanup_regions(regions, total_regions);
  get_region_stats(stats_start, regions, *total_regions);

  for (k = 0; k < *total_regions; k++) {
    regions[k].start += offset;
    regions[k].last += offset;
  }
}

void AV1RateControlQMode::SetRcParam(const RateControlParam &rc_param) {
  rc_param_ = rc_param;
}

GopChunkList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  GopChunkList gop_chuck_list;
  std::vector<REGIONS> regions_list;
  for (auto &stat : firstpass_info) {
    REGIONS region;
    int total_regions = 0;
    identify_regions(&stat, 16, 0, &region, &total_regions);
    regions_list.push_back(region);
  }

  return gop_chuck_list;
}

std::vector<FrameEncodeParameters> AV1RateControlQMode::GetGopEncodeInfo(
    const TplGopStats &tpl_stats_list) {
  std::vector<FrameEncodeParameters> frame_encoder_param;
  (void)tpl_stats_list;
  return frame_encoder_param;
}
}  // namespace aom