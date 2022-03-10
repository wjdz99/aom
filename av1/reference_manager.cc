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

#include "av1/reference_manager.h"

namespace aom {

void RefFrameManager::Reset() {
  ref_frame_table_.resize(max_ref_frames_);
  for (int i = 0; i < max_ref_frames_; ++i) {
    free_ref_idx_list_.push_back(i);
  }
}

// void ref_frame_manager_show(const RefFrameManager *ref_frame_manager) {
//   printf("=\n");
//   printf("forward: ");
//   for (auto ref_idx : ref_frame_manager->forward_stack) {
//     printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
//   }
//   printf("\n");
//   printf("backward: ");
//   for (auto ref_idx : ref_frame_manager->backward_queue) {
//     printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
//   }
//   printf("\n");
//   printf("last: ");
//   for (auto ref_idx : ref_frame_manager->last_queue) {
//     printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
//   }
//   printf("\n");
// }

void RefFrameManager::FreeRefIdx() {
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
}

void RefFrameManager::UpdateOrder(int order_idx) {
  if (!forward_stack_.empty()) {
    int ref_idx = forward_stack_.back();
    const GopFrame &gf_frame = ref_frame_table_[ref_idx];
    if (gf_frame.order_idx <= order_idx) {
      forward_stack_.pop_back();
      if (gf_frame.is_golden_frame) {
        // high quality frame
        backward_queue_.push_back(ref_idx);
      } else {
        last_queue_.push_back(ref_idx);
      }
    }
  }
}

int RefFrameManager::ColocatedRefIdx(int order_idx) {
  if (forward_stack_.size() == 0) return -1;
  int ref_idx = forward_stack_.back();
  int arf_order_idx = ref_frame_table_[ref_idx].order_idx;
  if (arf_order_idx == order_idx) {
    return ref_idx;
  }
  return -1;
}

void RefFrameManager::UpdateFrame(GopFrame *gop_frame,
                                  RefUpdateType ref_update_type,
                                  EncodeRefMode encode_ref_mode) {
  gop_frame->colocated_ref_idx = ColocatedRefIdx(gop_frame->order_idx);
  if (gop_frame->is_show_frame) {
    UpdateOrder(gop_frame->order_idx);
  }
  if (ref_update_type == RefUpdateType::kNone) {
    gop_frame->update_ref_idx = -1;
  } else {
    FreeRefIdx();
    int ref_idx = free_ref_idx_list_.front();
    free_ref_idx_list_.pop_front();
    gop_frame->update_ref_idx = ref_idx;
    if (ref_update_type == RefUpdateType::kForward) {
      forward_stack_.push_back(ref_idx);
    } else if (ref_update_type == RefUpdateType::kBackward) {
      backward_queue_.push_back(ref_idx);
    } else if (ref_update_type == RefUpdateType::kLast) {
      last_queue_.push_back(ref_idx);
    }
    ref_frame_table_[ref_idx] = *gop_frame;
  }
  gop_frame->encode_ref_mode = encode_ref_mode;
}
}  // namespace aom