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

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "av1/reference_manager.h"
#include "av1/ratectrl_qmode.h"

namespace aom {

void RefFrameManager::Reset() {
  free_ref_idx_list_.clear();
  for (int i = 0; i < kRefFrameTableSize; ++i) {
    free_ref_idx_list_.push_back(i);
  }
  forward_stack_.clear();
  backward_queue_.clear();
  last_queue_.clear();
}

int RefFrameManager::AllocateRefIdx() {
  if (free_ref_idx_list_.empty()) {
    size_t backward_size = backward_queue_.size();
    size_t last_size = last_queue_.size();
    if (last_size >= backward_size) {
      int ref_idx = last_queue_.front();
      last_queue_.pop_front();
      free_ref_idx_list_.push_back(ref_idx);
    } else {
      int ref_idx = backward_queue_.front();
      backward_queue_.pop_front();
      free_ref_idx_list_.push_back(ref_idx);
    }
  }

  int ref_idx = free_ref_idx_list_.front();
  free_ref_idx_list_.pop_front();
  return ref_idx;
}

int RefFrameManager::GetRefFrameCountByType(
    RefUpdateType ref_update_type) const {
  size_t cnt = 0;
  switch (ref_update_type) {
    case RefUpdateType::kForward: cnt = forward_stack_.size(); break;
    case RefUpdateType::kBackward: cnt = backward_queue_.size(); break;
    case RefUpdateType::kLast: cnt = last_queue_.size(); break;
    case RefUpdateType::kNone: cnt = 0; break;
  }
  return static_cast<int>(cnt);
}

int RefFrameManager::GetRefFrameCount() const {
  return GetRefFrameCountByType(RefUpdateType::kForward) +
         GetRefFrameCountByType(RefUpdateType::kBackward) +
         GetRefFrameCountByType(RefUpdateType::kLast);
}

// TODO(angiebird): Add unit test.
// Find the ref_idx corresponding to a ref_update_type.
// Return -1 if no ref frame is found.
// The priority_idx indicate closeness between the current frame and
// the ref frame in display order.
// For example, ref_update_type == kForward and priority_idx == 0 means
// find the closest ref frame in forward_stack_.
int RefFrameManager::GetRefFrameIdxByPriority(RefUpdateType ref_update_type,
                                              int priority_idx) const {
  if (ref_update_type == RefUpdateType::kForward) {
    int size = static_cast<int>(forward_stack_.size());
    if (priority_idx < size) {
      return forward_stack_[size - priority_idx - 1];
    }
  } else if (ref_update_type == RefUpdateType::kBackward) {
    int size = static_cast<int>(backward_queue_.size());
    if (priority_idx < size) {
      return backward_queue_[size - priority_idx - 1];
    }
  } else if (ref_update_type == RefUpdateType::kLast) {
    int size = static_cast<int>(last_queue_.size());
    if (priority_idx < size) {
      return last_queue_[size - priority_idx - 1];
    }
  }
  return -1;
}

// The priority_idx indicate closeness between the current frame and
// the ref frame in display order.
// For example, ref_update_type == kForward and priority_idx == 0 means
// find the closest ref frame in forward_stack_.
GopFrame RefFrameManager::GetRefFrameByPriority(RefUpdateType ref_update_type,
                                                int priority_idx) const {
  int ref_idx = GetRefFrameIdxByPriority(ref_update_type, priority_idx);
  if (ref_idx == -1) {
    return gop_frame_invalid();
  }
  assert(ref_frame_table_[ref_idx].update_ref_idx == ref_idx);
  return ref_frame_table_[ref_idx];
}

GopFrame RefFrameManager::GetRefFrameByIndex(int ref_idx) const {
  return ref_frame_table_[ref_idx];
}

ReferenceName get_ref_name(RefUpdateType ref_update_type, int priority_idx,
                           const std::set<ReferenceName> &used_name_set) {
  // TODO(angiebird): Find the better way to assign name lists.
  // Maybe sort the names based on how frequent each name is being used in the
  // past?
  const std::vector<ReferenceName> forward_name_list{
    ReferenceName::kBwdrefFrame, ReferenceName::kAltref2Frame,
    ReferenceName::kAltrefFrame, ReferenceName::kGoldenFrame,
    ReferenceName::kLastFrame,   ReferenceName::kLast2Frame,
    ReferenceName::kLast3Frame
  };
  const std::vector<ReferenceName> backward_name_list{
    ReferenceName::kGoldenFrame, ReferenceName::kLastFrame,
    ReferenceName::kLast2Frame,  ReferenceName::kLast3Frame,
    ReferenceName::kBwdrefFrame, ReferenceName::kAltref2Frame,
    ReferenceName::kAltrefFrame
  };
  const std::vector<ReferenceName> last_name_list{
    ReferenceName::kLastFrame,   ReferenceName::kLast2Frame,
    ReferenceName::kLast3Frame,  ReferenceName::kGoldenFrame,
    ReferenceName::kBwdrefFrame, ReferenceName::kAltref2Frame,
    ReferenceName::kAltrefFrame
  };

  const std::vector<ReferenceName> *name_list = nullptr;
  switch (ref_update_type) {
    case RefUpdateType::kForward: name_list = &forward_name_list; break;
    case RefUpdateType::kBackward: name_list = &backward_name_list; break;
    case RefUpdateType::kLast: name_list = &last_name_list; break;
    case RefUpdateType::kNone: break;
  }

  if (name_list) {
    const int name_list_size = static_cast<int>(name_list->size());
    for (int idx = priority_idx; idx < name_list_size; ++idx) {
      ReferenceName ref_name = name_list->at(idx);
      bool not_used = used_name_set.find(ref_name) == used_name_set.end();
      if (not_used) return ref_name;
    }
  }
  return ReferenceName::kNoneFrame;
}

std::vector<ReferenceFrame> RefFrameManager::GetRefFrameList() const {
  constexpr int round_robin_size = 3;
  const std::vector<RefUpdateType> round_robin_list{ RefUpdateType::kForward,
                                                     RefUpdateType::kBackward,
                                                     RefUpdateType::kLast };
  std::vector<int> priority_idx_list(round_robin_size, 0);
  int available_ref_frames = GetRefFrameCount();
  std::vector<ReferenceFrame> ref_frame_list;
  int ref_frame_count = 0;
  int round_robin_idx = 0;
  std::set<ReferenceName> used_name_set;
  while (ref_frame_count < max_ref_frames_ && available_ref_frames > 0) {
    const RefUpdateType ref_update_type = round_robin_list[round_robin_idx];
    int priority_idx = priority_idx_list[round_robin_idx];
    int ref_idx = GetRefFrameIdxByPriority(ref_update_type, priority_idx);
    if (ref_idx != -1) {
      const ReferenceName name =
          get_ref_name(ref_update_type, priority_idx, used_name_set);
      assert(name != ReferenceName::kNoneFrame);
      used_name_set.insert(name);
      ReferenceFrame ref_frame = { ref_idx, name };
      ref_frame_list.push_back(ref_frame);
      --available_ref_frames;
      ++priority_idx_list[round_robin_idx];
    }
    ref_frame_count = static_cast<int>(ref_frame_list.size());
    round_robin_idx = (round_robin_idx + 1) % round_robin_size;
  }
  return ref_frame_list;
}

void RefFrameManager::UpdateOrder(int global_order_idx) {
  cur_global_order_idx_ = global_order_idx;
  if (forward_stack_.empty()) {
    return;
  }
  int ref_idx = forward_stack_.back();
  const GopFrame &gf_frame = ref_frame_table_[ref_idx];
  if (gf_frame.global_order_idx <= global_order_idx) {
    forward_stack_.pop_back();
    if (gf_frame.is_golden_frame) {
      // high quality frame
      backward_queue_.push_back(ref_idx);
    } else {
      last_queue_.push_back(ref_idx);
    }
  }
}

int RefFrameManager::ColocatedRefIdx(int global_order_idx) {
  if (forward_stack_.size() == 0) return -1;
  int ref_idx = forward_stack_.back();
  int arf_global_order_idx = ref_frame_table_[ref_idx].global_order_idx;
  if (arf_global_order_idx == global_order_idx) {
    return ref_idx;
  }
  return -1;
}

static RefUpdateType infer_ref_update_type(const GopFrame &gop_frame,
                                           int cur_global_order_idx) {
  if (gop_frame.global_order_idx > cur_global_order_idx) {
    return RefUpdateType::kForward;
  }
  if (gop_frame.is_golden_frame) {
    return RefUpdateType::kBackward;
  }
  if (gop_frame.encode_ref_mode == EncodeRefMode::kShowExisting ||
      gop_frame.encode_ref_mode == EncodeRefMode::kOverlay) {
    return RefUpdateType::kNone;
  }
  return RefUpdateType::kLast;
}

// Find bits such that (1<<bits) > value
static int get_bits(int value) {
  int bits = 0;
  while ((1 << bits) <= value) {
    ++bits;
  }
  return bits;
}

// Compute two gop_frame's distance based on layer_depth difference,
// frame types and order_idx. Two frames with shorter distance are
// more likely to have similar probability model.
int gop_frame_probability_model_distance(const GopFrame &frame_a,
                                         const GopFrame &frame_b,
                                         int max_order_dist) {
  int dist = abs(frame_a.layer_depth - frame_b.layer_depth);

  dist = (dist << 1) + (frame_a.is_key_frame != frame_b.is_key_frame);
  dist = (dist << 1) + (frame_a.is_golden_frame != frame_b.is_golden_frame);
  dist = (dist << 1) + (frame_a.is_arf_frame != frame_b.is_arf_frame);
  dist = (dist << 1) + (frame_a.is_show_frame != frame_b.is_show_frame);
  dist = (dist << 1) + (frame_a.encode_ref_mode != frame_b.encode_ref_mode);

  int order_dist = abs(frame_a.order_idx - frame_b.order_idx);
  int order_dist_bits = get_bits(max_order_dist);
  dist = (dist << order_dist_bits) + order_dist;
  return dist;
}

// Pick primary_ref_idx for probability model.
int RefFrameManager::GetProbabilityModelRefIdx(const GopFrame &gop_frame) {
  assert(gop_frame.is_valid);

  int max_order_dist = 0;
  for (const GopFrame &ref_frame : ref_frame_table_) {
    if (ref_frame.is_valid) {
      const int order_dist = abs(gop_frame.order_idx - ref_frame.order_idx);
      max_order_dist = std::max(max_order_dist, order_dist);
    }
  }

  std::vector<std::pair<int, int>> candidate_list;
  for (int ref_idx = 0; ref_idx < static_cast<int>(ref_frame_table_.size());
       ++ref_idx) {
    const GopFrame &ref_frame = ref_frame_table_[ref_idx];
    if (ref_frame.is_valid) {
      assert(ref_frame.update_ref_idx == ref_idx);
      int dist = gop_frame_probability_model_distance(gop_frame, ref_frame,
                                                      max_order_dist);
      std::pair<int, int> candidate = { dist, ref_idx };
      candidate_list.push_back(candidate);
    }
  }

  sort(candidate_list.begin(), candidate_list.end());

  if (candidate_list.size() > 0) {
    return candidate_list[0].second;  // ref_idx
  }
  return -1;
}

void RefFrameManager::UpdateRefFrameTable(GopFrame *gop_frame) {
  gop_frame->probability_model_primary_ref_idx =
      GetProbabilityModelRefIdx(*gop_frame);
  gop_frame->ref_frame_list = GetRefFrameList();
  gop_frame->colocated_ref_idx = ColocatedRefIdx(gop_frame->global_order_idx);

  if (gop_frame->is_show_frame) {
    UpdateOrder(gop_frame->global_order_idx);
  }
  // Call infer_ref_update_type() after UpdateOrder() so that
  // cur_global_order_idx_ is up-to-date
  RefUpdateType ref_update_type =
      infer_ref_update_type(*gop_frame, cur_global_order_idx_);
  if (ref_update_type == RefUpdateType::kNone) {
    gop_frame->update_ref_idx = -1;
  } else {
    const int ref_idx = AllocateRefIdx();
    gop_frame->update_ref_idx = ref_idx;
    switch (ref_update_type) {
      case RefUpdateType::kForward: forward_stack_.push_back(ref_idx); break;
      case RefUpdateType::kBackward: backward_queue_.push_back(ref_idx); break;
      case RefUpdateType::kLast: last_queue_.push_back(ref_idx); break;
      case RefUpdateType::kNone: break;
    }
    ref_frame_table_[ref_idx] = *gop_frame;
  }
}

}  // namespace aom
