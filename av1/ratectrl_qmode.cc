
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
#include <vector>
#include "av1/ratectrl_qmode.h"

namespace aom {
void ref_frame_manager_reset(RefFrameManager *ref_frame_manager) {
  ref_frame_manager->ref_frame_table.resize(ref_frame_manager->max_ref_frames);
  for (int i = 0; i < ref_frame_manager->max_ref_frames; ++i) {
    ref_frame_manager->free_ref_idx_list.push_back(i);
  }
}
void ref_frame_manager_init(RefFrameManager *ref_frame_manager,
                            int max_ref_frames) {
  ref_frame_manager->max_ref_frames = max_ref_frames;
  ref_frame_manager->forward_max_size = max_ref_frames - 2;
  ref_frame_manager_reset(ref_frame_manager);
}

enum class RefUpdateType { kForward, kBackward, kLast, kNone };

void ref_frame_manager_show(const RefFrameManager *ref_frame_manager) {
  printf("=\n");
  printf("forward: ");
  for (auto ref_idx : ref_frame_manager->forward_stack) {
    printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
  }
  printf("\n");
  printf("backward: ");
  for (auto ref_idx : ref_frame_manager->backward_queue) {
    printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
  }
  printf("\n");
  printf("last: ");
  for (auto ref_idx : ref_frame_manager->last_queue) {
    printf("%d ", ref_frame_manager->ref_frame_table[ref_idx].order_idx);
  }
  printf("\n");
}

void ref_frame_manager_free_ref_idx(RefFrameManager *ref_frame_manager) {
  if (ref_frame_manager->free_ref_idx_list.size() == 0) {
    int backward_size = ref_frame_manager->backward_queue.size();
    int last_size = ref_frame_manager->last_queue.size();
    if (last_size >= backward_size) {
      int ref_idx = ref_frame_manager->last_queue.front();
      ref_frame_manager->last_queue.pop_front();
      ref_frame_manager->free_ref_idx_list.push_back(ref_idx);
    } else {
      int ref_idx = ref_frame_manager->backward_queue.front();
      ref_frame_manager->backward_queue.pop_front();
      ref_frame_manager->free_ref_idx_list.push_back(ref_idx);
    }
  }
}

void ref_frame_manager_update_order(RefFrameManager *ref_frame_manager,
                                    int order_idx) {
  if (ref_frame_manager->forward_stack.size() > 0) {
    int ref_idx = ref_frame_manager->forward_stack.back();
    const GopFrame &gf_frame = ref_frame_manager->ref_frame_table[ref_idx];
    if (gf_frame.order_idx <= order_idx) {
      ref_frame_manager->forward_stack.pop_back();
      if (gf_frame.is_golden_frame) {
        // high quality frame
        ref_frame_manager->backward_queue.push_back(ref_idx);
      } else {
        ref_frame_manager->last_queue.push_back(ref_idx);
      }
    }
  }
}

int ref_frame_manager_colocated_ref_idx(
    const RefFrameManager *ref_frame_manager, int order_idx) {
  if (ref_frame_manager->forward_stack.size() == 0) return -1;
  int ref_idx = ref_frame_manager->forward_stack.back();
  int arf_order_idx = ref_frame_manager->ref_frame_table[ref_idx].order_idx;
  if (arf_order_idx == order_idx) {
    return ref_idx;
  }
  return -1;
}

void ref_frame_manager_update_frame(RefFrameManager *ref_frame_manager,
                                    GopFrame *gop_frame,
                                    RefUpdateType ref_update_type,
                                    EncodeRefMode encode_ref_mode) {
  gop_frame->colocated_ref_idx = ref_frame_manager_colocated_ref_idx(
      ref_frame_manager, gop_frame->order_idx);
  if (gop_frame->is_show_frame) {
    ref_frame_manager_update_order(ref_frame_manager, gop_frame->order_idx);
  }
  if (ref_update_type == RefUpdateType::kNone) {
    gop_frame->update_ref_idx = -1;
  } else {
    ref_frame_manager_free_ref_idx(ref_frame_manager);
    int ref_idx = ref_frame_manager->free_ref_idx_list.front();
    ref_frame_manager->free_ref_idx_list.pop_front();
    gop_frame->update_ref_idx = ref_idx;
    if (ref_update_type == RefUpdateType::kForward) {
      ref_frame_manager->forward_stack.push_back(ref_idx);
    } else if (ref_update_type == RefUpdateType::kBackward) {
      ref_frame_manager->backward_queue.push_back(ref_idx);
    } else if (ref_update_type == RefUpdateType::kLast) {
      ref_frame_manager->last_queue.push_back(ref_idx);
    }
    ref_frame_manager->ref_frame_table[ref_idx] = *gop_frame;
  }
  gop_frame->encode_ref_mode = encode_ref_mode;
}

static GopFrame gop_frame_basic(int coding_idx, int order_idx, int is_key_frame,
                                int is_arf_frame, int is_golden_frame,
                                int is_show_frame) {
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
  int coding_idx = gop_struct->gop_frame_list.size();
  GopFrame gop_frame;
  int num_frames = order_end - order_start;
  // If there are less than 3 frames, stop introducing ARF
  if (depth < max_depth && num_frames < 3) {
    int order_mid = (order_start + order_end) / 2;
    // intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 1, 0, 0);
    ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                   RefUpdateType::kForward,
                                   EncodeRefMode::kRegular);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_start, order_mid);
    // show existing intermediate ARF
    gop_frame = gop_frame_basic(coding_idx, order_mid, 0, 0, 0, 1);
    ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                   RefUpdateType::kNone,
                                   EncodeRefMode::kShowExisting);
    gop_struct->gop_frame_list.push_back(gop_frame);
    construct_gop_multi_layer(gop_struct, ref_frame_manager, max_depth,
                              depth + 1, order_mid + 1, order_end);
  } else {
    // regular frame
    for (int i = order_start; i < order_end; ++i) {
      coding_idx = gop_struct->gop_frame_list.size();
      gop_frame = gop_frame_basic(coding_idx, i, 0, 0, 0, 1);
      ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                     RefUpdateType::kLast,
                                     EncodeRefMode::kRegular);
      gop_struct->gop_frame_list.push_back(gop_frame);
    }
  }
}

GopStruct construct_gop(RefFrameManager *ref_frame_manager,
                        int show_frame_count, int has_key_frame) {
  GopStruct gop_struct;
  gop_struct.show_frame_count = show_frame_count;
  int order_start = 0;
  int order_arf = show_frame_count - 1;
  int coding_idx;
  GopFrame gop_frame;
  if (has_key_frame) {
    ref_frame_manager_reset(ref_frame_manager);
    coding_idx = gop_struct.gop_frame_list.size();
    gop_frame = gop_frame_basic(coding_idx, order_start, 1, 0, 1, 1);
    ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                   RefUpdateType::kBackward,
                                   EncodeRefMode::kRegular);
    gop_struct.gop_frame_list.push_back(gop_frame);
    order_start++;
  }
  // ARF
  coding_idx = gop_struct.gop_frame_list.size();
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 1, 1, 0);
  ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                 RefUpdateType::kForward,
                                 EncodeRefMode::kRegular);
  gop_struct.gop_frame_list.push_back(gop_frame);
  construct_gop_multi_layer(&gop_struct, ref_frame_manager,
                            ref_frame_manager->forward_max_size, 1, order_start,
                            order_arf);
  // Overlay
  coding_idx = gop_struct.gop_frame_list.size();
  gop_frame = gop_frame_basic(coding_idx, order_arf, 0, 0, 0, 1);
  ref_frame_manager_update_frame(ref_frame_manager, &gop_frame,
                                 RefUpdateType::kNone, EncodeRefMode::kOverlay);
  gop_struct.gop_frame_list.push_back(gop_frame);
  return gop_struct;
}

void AV1RateControlQMode::SetRcParam(const RateControlParam &rc_param) {
  rc_param_ = rc_param;
}

GopStructList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  // A temporary simple implementation
  const int max_gop_show_frame_count = 16;
  int remaining_show_frame_count = firstpass_info.size();
  GopStructList gop_list;

  RefFrameManager ref_frame_manager;
  ref_frame_manager_init(&ref_frame_manager, rc_param_.max_ref_frames);

  while (remaining_show_frame_count > 0) {
    int show_frame_count =
        std::min(remaining_show_frame_count, max_gop_show_frame_count);
    // TODO(angiebird): determine gop show frame count based on first pass stats
    // here.
    int has_key_frame = gop_list.size() == 0;
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
