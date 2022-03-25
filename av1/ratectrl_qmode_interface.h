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

#ifndef AOM_AV1_RATECTRL_QMODE_INTERFACE_H_
#define AOM_AV1_RATECTRL_QMODE_INTERFACE_H_

#include <vector>
#include "av1/encoder/firstpass.h"

namespace aom {

constexpr int kBlockRefCount = 2;

struct MotionVector {
  int row;
  int col;
  int subpel_bits;  // see get_fullpel_value
};

struct RateControlParam {
  int max_gop_length;
  int min_gop_length;
  int max_ref_frames;
  int base_q_index;
};

struct TplBlockStats {
  int width;
  int height;
  int row;
  int col;
  int64_t intra_cost;
  int64_t inter_cost;
  MotionVector mv[kBlockRefCount];
  int ref_frame_index[kBlockRefCount];
};

enum class EncodeRefMode {
  kRegular,
  kOverlay,
  kShowExisting,
};

struct GopFrame {
  // basic info
  int order_idx;         // Index in display order
  int coding_idx;        // Index in coding order
  bool is_key_frame;     // If this is key frame, reset reference buffers are
                         // required
  bool is_arf_frame;     // Is this a forward frame, a frame with order_idx
                         // higher than the current display order
  bool is_show_frame;    // Is this frame a show frame after coding
  bool is_golden_frame;  // Is this a high quality frame

  // reference frame info
  EncodeRefMode encode_ref_mode;
  int colocated_ref_idx;  // colocated_ref_idx == -1 when encode_ref_mode ==
                          // EncodeRefMode::kRegular
  int update_ref_idx;  // The reference index that this frame should be updated
                       // to. update_ref_idx == -1 when this frame will not
                       // serve as a reference frame
};

struct GopStruct {
  int show_frame_count;
  std::vector<GopFrame> gop_frame_list;
};

using GopStructList = std::vector<GopStruct>;

struct FrameEncodeParameters {
  int q_index;
  int rdmult;
};

using FirstpassInfo = std::vector<FIRSTPASS_STATS>;
using RefFrameTable = std::vector<GopFrame>;

struct TplFrameStats {
  int min_block_size;
  int frame_width;
  int frame_height;
  std::vector<TplBlockStats> block_stats_list;
};

struct TplGopStats {
  std::vector<TplFrameStats> frame_stats_list;
};

struct TplFrameDepStats {
  int unit_size;  // equivalent to min_block_size
  std::vector<std::vector<double>> unit_stats;
};

struct TplGopDepStats {
  std::vector<TplFrameDepStats> frame_dep_stats_list;
};

class AV1RateControlQModeInterface {
 public:
  AV1RateControlQModeInterface();
  virtual ~AV1RateControlQModeInterface();

  virtual void SetRcParam(const RateControlParam &rc_param) = 0;
  virtual GopStructList DetermineGopInfo(
      const FirstpassInfo &firstpass_stats_list) = 0;
  // Accept firstpass and tpl info from the encoder and return q index and
  // rdmult. This needs to be called with consecutive GOPs as returned by
  // DetermineGopInfo.
  virtual std::vector<FrameEncodeParameters> GetGopEncodeInfo(
      const GopStruct &gop_struct, const TplGopStats &tpl_gop_stats,
      const RefFrameTable &ref_frame_table_snapshot) = 0;
};  // class AV1RateCtrlQMode
}  // namespace aom

#endif  // AOM_AV1_RATECTRL_QMODE_INTERFACE_H_
