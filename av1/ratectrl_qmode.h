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

#ifndef AOM_AV1_RATECTRL_QMODE_H_
#define AOM_AV1_RATECTRL_QMODE_H_

#include <deque>
#include <queue>
#include <vector>
#include "av1/encoder/firstpass.h"
#include "av1/ratectrl_qmode_interface.h"
#include "av1/reference_manager.h"

namespace aom {

constexpr int kLayerDepthOffset = 1;
constexpr int kMinIntervalToAddArf = 3;
constexpr int kMinArfInterval = (kMinIntervalToAddArf + 1) / 2;

struct TplUnitDepStats {
  double propagation_cost;
  double intra_cost;
  double inter_cost;
  std::array<MotionVector, kBlockRefCount> mv;
  std::array<int, kBlockRefCount> ref_frame_index;
};

struct TplFrameDepStats {
  int unit_size;  // equivalent to min_block_size
  std::vector<std::vector<TplUnitDepStats>> unit_stats;
};

struct TplGopDepStats {
  std::vector<TplFrameDepStats> frame_dep_stats_list;
};

GopFrame GopFrameInvalid();

// gop frame type used for facilitate setting up GopFrame
// TODO(angiebird): Define names for forward key frame and
// key frame with overlay
enum class GopFrameType {
  kRegularKey,  // High quality key frame without overlay
  kRegularArf,  // High quality arf with strong filtering followed by an overlay
                // later
  kIntermediateArf,  // Good quality arf with weak or no filtering followed by a
                     // show_existing later
  kRegularLeaf,      // Regular leaf frame
  kShowExisting,     // Show_existing frame
  kOverlay           // Overlay frame
};

// Set up is_key_frame, is_arf_frame, is_show_frame, is_golden_frame and
// encode_ref_mode in GopFrame based on gop_frame_type
void SetGopFrameByType(GopFrameType gop_frame_type, GopFrame *gop_frame);

GopFrame GopFrameBasic(int global_coding_idx_offset,
                       int global_order_idx_offset, int coding_idx,
                       int order_idx, int depth, GopFrameType gop_frame_type);

GopStruct ConstructGop(RefFrameManager *ref_frame_manager, int show_frame_count,
                       bool has_key_frame, int global_coding_idx_offset,
                       int global_order_idx_offset);

TplFrameDepStats CreateTplFrameDepStats(int frame_height, int frame_width,
                                        int min_block_size);

TplUnitDepStats TplBlockStatsToDepStats(const TplBlockStats &block_stats,
                                        int unit_count);

TplFrameDepStats CreateTplFrameDepStatsWithoutPropagation(
    const TplFrameStats &frame_stats);

std::vector<int> GetKeyFrameList(const FirstpassInfo &first_pass_info);

double TplFrameDepStatsAccumulateIntraCost(
    const TplFrameDepStats &frame_dep_stats);

double TplFrameDepStatsAccumulate(const TplFrameDepStats &frame_dep_stats);

void TplFrameDepStatsPropagate(int coding_idx,
                               const RefFrameTable &ref_frame_table,
                               TplGopDepStats *tpl_gop_dep_stats);

int GetBlockOverlapArea(int r0, int c0, int r1, int c1, int size);

TplGopDepStats ComputeTplGopDepStats(
    const TplGopStats &tpl_gop_stats,
    const std::vector<RefFrameTable> &ref_frame_table_list,
    const int frame_count);

class AV1RateControlQMode : public AV1RateControlQModeInterface {
 public:
  Status SetRcParam(const RateControlParam &rc_param) override;
  GopStructList DetermineGopInfo(const FirstpassInfo &firstpass_info) override;
  GopEncodeInfo GetGopEncodeInfo(
      const GopStruct &gop_struct, const TplGopStats &tpl_gop_stats,
      const RefFrameTable &ref_frame_table_snapshot) override;

  // Public for testing only.
  // Returns snapshots of the ref frame before and after each frame in
  // gop_struct. The returned list will have n+1 entries for n frames.
  // If this is first GOP, ref_frame_table is ignored and all refs are assumed
  // invalid; otherwise ref_frame_table is used as the initial state.
  std::vector<RefFrameTable> GetRefFrameTableList(
      const GopStruct &gop_struct, RefFrameTable ref_frame_table);

 private:
  RateControlParam rc_param_;
};
}  // namespace aom

#endif  // AOM_AV1_RATECTRL_QMODE_H_
