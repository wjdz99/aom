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

#include <algorithm>
#include <cassert>
#include <climits>
#include <vector>

#include "av1/encoder/pass2_strategy.h"

namespace aom {

static GopFrame gop_frame_basic(int coding_idx, int order_idx,
                                bool is_key_frame, bool is_arf_frame,
                                bool is_golden_frame, bool is_show_frame) {
  GopFrame gop_frame;
  gop_frame.coding_idx = coding_idx;
  gop_frame.order_idx = order_idx;
  gop_frame.is_key_frame = is_key_frame;
  gop_frame.is_arf_frame = is_arf_frame;
  gop_frame.is_golden_frame = is_golden_frame;
  gop_frame.is_show_frame = is_show_frame;
  gop_frame.encode_ref_mode = EncodeRefMode::kRegular;
  gop_frame.colocated_ref_idx = -1;
  gop_frame.update_ref_idx = -1;
  return gop_frame;
}

// This function create gop frames with indices of display order from
// order_start to order_end - 1. The function will recursively introduce
// intermediate ARF untill maximum depth is met or the number of regular frames
// in between two ARFs are less than 3. Than the regular frames will be added
// into the gop_struct.
void construct_gop_multi_layer(GopStruct *gop_struct,
                               RefFrameManager *ref_frame_manager,
                               int max_depth, int depth, int order_start,
                               int order_end) {
  int coding_idx = static_cast<int>(gop_struct->gop_frame_list.size());
  GopFrame gop_frame;
  int num_frames = order_end - order_start;
  // If there are less than 3 frames, stop introducing ARF
  if (depth < max_depth && num_frames < 3) {
    int order_mid = (order_start + order_end) / 2;
    // intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 1, 0, 0);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kForward,
                                   EncodeRefMode::kRegular);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_start, order_mid);
    // show existing intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 0, 0, 1);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kNone,
                                   EncodeRefMode::kShowExisting);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_mid + 1, order_end);
  } else {
    // regular frame
    for (int i = order_start; i < order_end; ++i) {
      coding_idx = static_cast<int>(gop_struct->gop_frame_list.size());
      gop_frame = gop_frame_basic(coding_idx, i, 0, 0, 0, 1);
      ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kLast,
                                     EncodeRefMode::kRegular);
      gop_struct->gop_frame_list.push_back(gop_frame);
    }
  }
}

GopStruct construct_gop(RefFrameManager *ref_frame_manager,
                        int show_frame_count, bool has_key_frame) {
  GopStruct gop_struct;
  gop_struct.show_frame_count = show_frame_count;
  int order_start = 0;
  int order_arf = show_frame_count - 1;
  int coding_idx;
  GopFrame gop_frame;
  if (has_key_frame) {
    ref_frame_manager->Reset();
    coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
    gop_frame = gop_frame_basic(coding_idx, order_start, 1, 0, 1, 1);
    ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kBackward,
                                   EncodeRefMode::kRegular);
    gop_struct.gop_frame_list.push_back(gop_frame);
    order_start++;
  }
  // ARF
  coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 1, 1, 0);
  ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kForward,
                                 EncodeRefMode::kRegular);
  gop_struct.gop_frame_list.push_back(gop_frame);
  construct_gop_multi_layer(&gop_struct, ref_frame_manager,
                            ref_frame_manager->ForwardMaxSize(), 1, order_start,
                            order_arf);
  // Overlay
  coding_idx = static_cast<int>(gop_struct.gop_frame_list.size());
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 0, 0, 1);
  ref_frame_manager->UpdateFrame(&gop_frame, RefUpdateType::kNone,
                                 EncodeRefMode::kOverlay);
  gop_struct.gop_frame_list.push_back(gop_frame);
  return gop_struct;
}

void AV1RateControlQMode::SetRcParam(const RateControlParam &rc_param) {
  rc_param_ = rc_param;
}

// initialize GF_GROUP_STATS
static void init_gf_stats(GF_GROUP_STATS *gf_stats) {
  gf_stats->gf_group_err = 0.0;
  gf_stats->gf_group_raw_error = 0.0;
  gf_stats->gf_group_skip_pct = 0.0;
  gf_stats->gf_group_inactive_zone_rows = 0.0;

  gf_stats->mv_ratio_accumulator = 0.0;
  gf_stats->decay_accumulator = 1.0;
  gf_stats->zero_motion_accumulator = 1.0;
  gf_stats->loop_decay_rate = 1.0;
  gf_stats->last_loop_decay_rate = 1.0;
  gf_stats->this_frame_mv_in_out = 0.0;
  gf_stats->mv_in_out_accumulator = 0.0;
  gf_stats->abs_mv_in_out_accumulator = 0.0;

  gf_stats->avg_sr_coded_error = 0.0;
  gf_stats->avg_pcnt_second_ref = 0.0;
  gf_stats->avg_new_mv_count = 0.0;
  gf_stats->avg_wavelet_energy = 0.0;
  gf_stats->avg_raw_err_stdev = 0.0;
  gf_stats->non_zero_stdev_count = 0;
}

static int input_stats(FIRSTPASS_STATS *stats_in_end, FIRSTPASS_STATS *stats_in,
                       FIRSTPASS_STATS *fps) {
  if (stats_in >= stats_in_end) return EOF;

  *fps = *stats_in;
  ++stats_in;
  return 1;
}

static int find_regions_index(const REGIONS *regions, int num_regions,
                              int frame_idx) {
  for (int k = 0; k < num_regions; k++) {
    if (regions[k].start <= frame_idx && regions[k].last >= frame_idx) {
      return k;
    }
  }
  return -1;
}

#define MIN_SHRINK_LEN 6

/*!\brief Determine the length of future GF groups.
 *
 * \ingroup gf_group_algo
 * This function decides the gf group length of future frames in batch
 *
 * \param[in]    cpi              Top-level encoder structure
 * \param[in]    max_gop_length   Maximum length of the GF group
 * \param[in]    max_intervals    Maximum number of intervals to decide
 *
 * \return Nothing is returned. Instead, cpi->ppi->rc.gf_intervals is
 * changed to store the decided GF group lengths.
 */
static std::vector<int> av1_calculate_gf_length(
    RateControlParam &rc_param, FIRSTPASS_STATS *stats_in,
    FIRSTPASS_STATS *stats_in_start, FIRSTPASS_STATS *stats_in_end,
    std::vector<REGIONS> &regions_list, int total_regions, int max_intervals,
    int min_gf_interval, int max_gf_interval, int frames_since_key,
    int frames_to_key) {
  // RATE_CONTROL *const rc = &cpi->rc;
  // PRIMARY_RATE_CONTROL *const p_rc = &cpi->ppi->p_rc;
  // TWO_PASS *const twopass = &cpi->ppi->twopass;
  FIRSTPASS_STATS next_frame;
  std::vector<int> gf_intervals(MAX_NUM_GF_INTERVALS);
  FIRSTPASS_STATS *const start_pos = stats_in;
  FIRSTPASS_STATS *const stats = start_pos;
  int i;
  av1_zero(next_frame);

  // TODO(urvang): Try logic to vary min and max interval based on q.
  const int active_min_gf_interval = min_gf_interval;
  const int active_max_gf_interval =
      std::min(max_gf_interval, rc_param.min_gop_length);
  const int min_shrink_int = std::max(MIN_SHRINK_LEN, active_min_gf_interval);

  i = frames_since_key == 0;
  int count_cuts = 1;
  // If cpi->gf_state.arf_gf_boost_lst is 0, we are starting with a KF or GF.
  int cur_start = 0, cur_last;
  int cut_pos[MAX_NUM_GF_INTERVALS + 1] = { -1 };
  int cut_here;
  GF_GROUP_STATS gf_stats;
  init_gf_stats(&gf_stats);
  while (count_cuts < max_intervals + 1) {
    FIRSTPASS_STATS *stats_in_loop = stats_in;
    // reaches next key frame, break here
    if (i >= frames_to_key) {
      cut_here = 2;
    } else if (i - cur_start >= max_gf_interval) {
      // reached maximum len, but nothing special yet (almost static)
      // let's look at the next interval
      cut_here = 1;
    } else if (EOF == input_stats(stats_in_end, stats_in_loop, &next_frame)) {
      // reaches last frame, break
      cut_here = 2;
    }

    if (cut_here) {
      cur_last = i - 1;  // the current last frame in the gf group
      int ori_last = cur_last;
      // The region frame idx does not start from the same frame as cur_start
      // and cur_last. Need to offset them.
      int offset = frames_since_key;
      REGIONS *regions = regions_list.data();
      int num_regions = total_regions;

      int scenecut_idx = -1;
      // only try shrinking if interval smaller than active_max_gf_interval
      if (cur_last - cur_start <= active_max_gf_interval &&
          cur_last > cur_start) {
        // find the region indices of where the first and last frame belong.
        int k_start =
            find_regions_index(regions, num_regions, cur_start + offset);
        int k_last =
            find_regions_index(regions, num_regions, cur_last + offset);
        if (cur_start + offset == 0) k_start = 0;

        // See if we have a scenecut in between
        for (int r = k_start + 1; r <= k_last; r++) {
          if (regions[r].type == SCENECUT_REGION &&
              regions[r].last - offset - cur_start > active_min_gf_interval) {
            scenecut_idx = r;
            break;
          }
        }

        // if the found scenecut is very close to the end, ignore it.
        if (regions[num_regions - 1].last - regions[scenecut_idx].last < 4) {
          scenecut_idx = -1;
        }

        if (scenecut_idx != -1) {
          // If we have a scenecut, then stop at it.
          // TODO(bohanli): add logic here to stop before the scenecut and for
          // the next gop start from the scenecut with GF
          int is_minor_sc =
              (regions[scenecut_idx].avg_cor_coeff *
                   (1 - stats[regions[scenecut_idx].start - offset].noise_var /
                            regions[scenecut_idx].avg_intra_err) >
               0.6);
          cur_last = regions[scenecut_idx].last - offset - !is_minor_sc;
        } else {
          int is_last_analysed = (k_last == num_regions - 1) &&
                                 (cur_last + offset == regions[k_last].last);
          int not_enough_regions =
              k_last - k_start <=
              1 + (regions[k_start].type == SCENECUT_REGION);
          // if we are very close to the end, then do not shrink since it may
          // introduce intervals that are too short
          if (!(is_last_analysed && not_enough_regions)) {
            const double arf_length_factor = 0.1;
            double best_score = 0;
            int best_j = -1;
            const int first_frame = regions[0].start - offset;
            const int last_frame = regions[num_regions - 1].last - offset;
            // score of how much the arf helps the whole GOP
            double base_score = 0.0;
            // Accumulate base_score in
            for (int j = cur_start + 1; j < cur_start + min_shrink_int; j++) {
              if (stats + j >= stats_in_end) break;
              base_score = (base_score + 1.0) * stats[j].cor_coeff;
            }
            int met_blending = 0;   // Whether we have met blending areas before
            int last_blending = 0;  // Whether the previous frame if blending
            for (int j = cur_start + min_shrink_int; j <= cur_last; j++) {
              if (stats + j >= stats_in_end) break;
              base_score = (base_score + 1.0) * stats[j].cor_coeff;
              int this_reg =
                  find_regions_index(regions, num_regions, j + offset);
              if (this_reg < 0) continue;
              // A GOP should include at most 1 blending region.
              if (regions[this_reg].type == BLENDING_REGION) {
                last_blending = 1;
                if (met_blending) {
                  break;
                } else {
                  base_score = 0;
                  continue;
                }
              } else {
                if (last_blending) met_blending = 1;
                last_blending = 0;
              }

              // Add the factor of how good the neighborhood is for this
              // candidate arf.
              double this_score = arf_length_factor * base_score;
              double temp_accu_coeff = 1.0;
              // following frames
              int count_f = 0;
              for (int n = j + 1; n <= j + 3 && n <= last_frame; n++) {
                if (stats + n >= stats_in_end) break;
                temp_accu_coeff *= stats[n].cor_coeff;
                this_score +=
                    temp_accu_coeff *
                    (1 - stats[n].noise_var /
                             AOMMAX(regions[this_reg].avg_intra_err, 0.001));
                count_f++;
              }
              // preceding frames
              temp_accu_coeff = 1.0;
              for (int n = j; n > j - 3 * 2 + count_f && n > first_frame; n--) {
                if (stats + n < stats_in_start) break;
                temp_accu_coeff *= stats[n].cor_coeff;
                this_score +=
                    temp_accu_coeff *
                    (1 - stats[n].noise_var /
                             AOMMAX(regions[this_reg].avg_intra_err, 0.001));
              }

              if (this_score > best_score) {
                best_score = this_score;
                best_j = j;
              }
            }

            // For blending areas, move one more frame in case we missed the
            // first blending frame.
            int best_reg =
                find_regions_index(regions, num_regions, best_j + offset);
            if (best_reg < num_regions - 1 && best_reg > 0) {
              if (regions[best_reg - 1].type == BLENDING_REGION &&
                  regions[best_reg + 1].type == BLENDING_REGION) {
                if (best_j + offset == regions[best_reg].start &&
                    best_j + offset < regions[best_reg].last) {
                  best_j += 1;
                } else if (best_j + offset == regions[best_reg].last &&
                           best_j + offset > regions[best_reg].start) {
                  best_j -= 1;
                }
              }
            }

            if (cur_last - best_j < 2) best_j = cur_last;
            if (best_j > 0 && best_score > 0.1) cur_last = best_j;
            // if cannot find anything, just cut at the original place.
          }
        }
      }
      cut_pos[count_cuts] = cur_last;
      count_cuts++;

      // reset pointers to the shrinked location
      stats_in_loop = start_pos + cur_last;
      cur_start = cur_last;
      int cur_region_idx =
          find_regions_index(regions, num_regions, cur_start + 1 + offset);
      if (cur_region_idx >= 0)
        if (regions[cur_region_idx].type == SCENECUT_REGION) cur_start++;

      i = cur_last;

      if (cut_here > 1 && cur_last == ori_last) break;

      // reset accumulators
      init_gf_stats(&gf_stats);
    }
    ++i;
  }

  // save intervals
  // rc->intervals_till_gf_calculate_due = count_cuts - 1;
  for (int n = 1; n < count_cuts; n++) {
    gf_intervals[n - 1] = cut_pos[n] - cut_pos[n - 1];
  }

  return gf_intervals;
}

GopStructList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  std::vector<REGIONS> regions_list(MAX_FIRSTPASS_ANALYSIS_FRAMES);
  int total_regions = 0;
  // TODO(jianj): firstpass_info.size() should eventually be replaced
  // by the number of frames to the next KF.
  assert(firstpass_info.size() <= INT_MAX);
  av1_identify_regions(firstpass_info.data(),
                       std::min(static_cast<int>(firstpass_info.size()),
                                MAX_FIRSTPASS_ANALYSIS_FRAMES),
                       0, regions_list.data(), &total_regions);

  // A temporary simple implementation
  const int max_gop_show_frame_count = 16;
  int remaining_show_frame_count = static_cast<int>(firstpass_info.size());
  GopStructList gop_list;

  RefFrameManager ref_frame_manager(rc_param_.max_ref_frames);

  while (remaining_show_frame_count > 0) {
    int show_frame_count =
        std::min(remaining_show_frame_count, max_gop_show_frame_count);
    // TODO(angiebird): determine gop show frame count based on first pass stats
    // here.
    bool has_key_frame = gop_list.size() == 0;
    GopStruct gop =
        construct_gop(&ref_frame_manager, show_frame_count, has_key_frame);
    gop_list.push_back(gop);
    remaining_show_frame_count -= show_frame_count;
  }
  return gop_list;
}

std::vector<FrameEncodeParameters> AV1RateControlQMode::GetGopEncodeInfo(
    const TplGopStats &tpl_stats_list) {
  std::vector<FrameEncodeParameters> frame_encoder_param;
  (void)tpl_stats_list;
  return frame_encoder_param;
}

}  // namespace aom
