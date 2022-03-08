
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

#include <vector>
#include <queue>
#include "av1/ratectrl_qmode.h"

namespace aom {
struct RefFrameInfo {
  int ref_frame_idx;
  int order_idx;
  int coding_idx;
};

struct RefFrameManager {
  int max_ref_frames;
  std::queue<int> free_buf_idx_list;
  int arf_max_size;
  std::vector<RefFrameInfo> arf_stack;
  int gld_max_size;
  std::queue<RefFrameInfo> gld_queue;
  int lst_max_size;
  std::queue<RefFrameInfo> lst_queue;
};

void ref_frame_manager_init(RefFrameManager *ref_frame_manager,
                            int max_ref_frames) {
  ref_frame_manager->max_ref_frames = max_ref_frames;
  for (int i = 0; i < max_ref_frames; ++i) {
    ref_frame_manager->free_buf_idx_list.push(i);
  }
  ref_frame_manager->arf_max_size = max_ref_frames - 2;
  ref_frame_manager->lst_max_size = max_ref_frames - 2;
  ref_frame_manager->gld_max_size = max_ref_frames - 2;
}

void ref_frame_manager_set_gop_frame_ref_info(
    const RefFrameManager &ref_frame_manager, GopFrame *gop_frame,
    EncodeRefMode encode_ref_mode) {
  gop_frame->encode_ref_mode = encode_ref_mode;
  if (encode_ref_mode == Regular) {
    gop_frame->colocated_ref_frame_idx = -1;
  } else {
    gop_frame->colocated_ref_frame_idx =
        ref_frame_manager.arf_stack.back().ref_frame_idx;
  }
}

GopFrame gop_frame_basic(int coding_idx, int order_idx, int is_key_frame,
                         int is_arf_frame, int is_show_frame) {
  GopFrame gop_frame;
  gop_frame.coding_idx = coding_idx;
  gop_frame.order_idx = order_idx;
  gop_frame.is_key_frame = is_key_frame;
  gop_frame.is_arf_frame = is_arf_frame;
  gop_frame.is_show_frame = is_show_frame;
  return gop_frame;
}

void construct_gop_multi_layer(GopStruct *gop_struct, int max_depth, int depth,
                               int order_start, int order_end) {
  int coding_idx = gop_struct->gop_frame_list.size();
  // ARF
  GopFrame gop_frame = gop_frame_basic(coding_idx, order_end, 0, 1, 0);
  gop_struct->gop_frame_list.push_back(gop_frame);
  if (depth < max_depth || order_start < order_end) {
    int order_mid = (order_start + order_end) / 2;
    construct_gop_multi_layer(gop_struct, max_depth, depth + 1, order_start,
                              order_mid);
    construct_gop_multi_layer(gop_struct, max_depth, depth + 1, order_mid + 1,
                              order_end);
  } else {
    // Regular Frame
    for (int i = order_start; i < order_end; ++i) {
      coding_idx = gop_struct->gop_frame_list.size();
      GopFrame gop_frame = gop_frame_basic(coding_idx, i, 0, 0, 1);
      gop_struct->gop_frame_list.push_back(gop_frame);
    }
    // Overlay Frame
    coding_idx = gop_struct->gop_frame_list.size();
    GopFrame gop_frame = gop_frame_basic(coding_idx, order_end, 0, 0, 1);
    gop_struct->gop_frame_list.push_back(gop_frame);
  }
};

void construct_gop(GopStruct *gop_struct, RefFrameManager *ref_frame_manager,
                   int show_frame_count) {
  gop_struct->show_frame_count = show_frame_count;
  construct_gop_multi_layer(gop_struct, ref_frame_manager->arf_max_size, 0, 0,
                            show_frame_count - 1);
}

class AV1RateControlQMode : public AV1RateControlQModeInterface {
 public:
  AV1RateControlQMode() = default;
  ~AV1RateControlQMode() = default;
  void SetRcParam(const RateControlParam &rc_param);
  GopStructList DetermineGopInfo(const FirstpassInfo &firstpass_info);

 private:
  RateControlParam rc_param_;
};

void AV1RateControlQMode::SetRcParam(const RateControlParam &rc_param) {
  rc_param_ = rc_param;
}

GopStructList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  const int max_gop_show_frame_count = 16;
  int remain_show_frame_count = firstpass_info.size();
  GopStructList gop_list;

  RefFrameManager ref_frame_manager;
  ref_frame_manager_init(&ref_frame_manager, rc_param_.max_ref_frames);

  while (remain_show_frame_count > 0) {
    int show_frame_count =
        std::min(remain_show_frame_count, max_gop_show_frame_count);
    // TODO(angiebird): determine gop show frame count here.
    GopStruct gop;
    construct_gop(&gop, &ref_frame_manager, show_frame_count);
    gop_list.push_back(gop);
    remain_show_frame_count -= show_frame_count;
  }
  return gop_list;
}
}  // namespace aom
