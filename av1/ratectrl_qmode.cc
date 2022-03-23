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
#include <numeric>
#include <vector>

#include "av1/encoder/pass2_strategy.h"
#include "av1/encoder/tpl_model.h"

namespace aom {

GopFrame gop_frame_invalid() {
  GopFrame gop_frame = {};
  gop_frame.is_valid = false;
  gop_frame.coding_idx = -1;
  gop_frame.order_idx = -1;
  return gop_frame;
}

GopFrame gop_frame_basic(int coding_idx, int order_idx, bool is_key_frame,
                         bool is_arf_frame, bool is_golden_frame,
                         bool is_show_frame) {
  GopFrame gop_frame;
  gop_frame.is_valid = true;
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

GopStructList AV1RateControlQMode::DetermineGopInfo(
    const FirstpassInfo &firstpass_info) {
  std::vector<REGIONS> regions_list(MAX_FIRSTPASS_ANALYSIS_FRAMES);
  int total_regions = 0;
  // TODO(jianj): firstpass_info.size() should eventually be replaced
  // by the number of frames to the next KF.
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

TplFrameDepStats create_tpl_frame_dep_stats_empty(int frame_height,
                                                  int frame_width,
                                                  int min_block_size) {
  const int unit_rows =
      frame_height / min_block_size + !!(frame_height % min_block_size);
  const int unit_cols =
      frame_width / min_block_size + !!(frame_width % min_block_size);
  TplFrameDepStats frame_dep_stats;
  frame_dep_stats.unit_size = min_block_size;
  frame_dep_stats.unit_stats = std::vector<std::vector<double>>(
      unit_rows, std::vector<double>(unit_cols, 0));
  return frame_dep_stats;
}

TplFrameDepStats create_tpl_frame_dep_stats_wo_propagation(
    const TplFrameStats &frame_stats) {
  const int min_block_size = frame_stats.min_block_size;
  TplFrameDepStats frame_dep_stats = create_tpl_frame_dep_stats_empty(
      frame_stats.frame_height, frame_stats.frame_width, min_block_size);
  for (const TplBlockStats &block_stats : frame_stats.block_stats_list) {
    const int block_unit_rows = block_stats.height / min_block_size;
    const int block_unit_cols = block_stats.width / min_block_size;
    const int unit_count = block_unit_rows * block_unit_cols;
    const int block_unit_row = block_stats.row / min_block_size;
    const int block_unit_col = block_stats.col / min_block_size;
    const double cost_diff =
        (block_stats.inter_cost - block_stats.intra_cost) * 1.0 / unit_count;
    for (int r = 0; r < block_unit_rows; r++) {
      for (int c = 0; c < block_unit_cols; c++) {
        frame_dep_stats.unit_stats[block_unit_row + r][block_unit_col + c] =
            cost_diff;
      }
    }
  }
  return frame_dep_stats;
}

int get_ref_coding_idx_list(const TplBlockStats &block_stats,
                            const RefFrameTable &ref_frame_table,
                            int *ref_coding_idx_list) {
  int ref_frame_count = 0;
  for (int i = 0; i < kBlockRefCount; ++i) {
    ref_coding_idx_list[i] = -1;
    int ref_frame_index = block_stats.ref_frame_index[i];
    if (ref_frame_index != -1) {
      ref_coding_idx_list[i] = ref_frame_table[ref_frame_index].coding_idx;
      ref_frame_count++;
    }
  }
  return ref_frame_count;
}

int get_block_overlap_area(int r0, int c0, int r1, int c1, int size) {
  const int r_low = std::max(r0, r1);
  const int r_high = std::min(r0 + size, r1 + size);
  const int c_low = std::max(c0, c1);
  const int c_high = std::min(c0 + size, c1 + size);
  if (r_high >= r_low && c_high >= c_low) {
    return (r_high - r_low) * (c_high - c_low);
  }
  return 0;
}

double tpl_frame_stats_accumulate(const TplFrameStats &frame_stats) {
  double ref_sum_cost_diff = 0;
  for (auto &block_stats : frame_stats.block_stats_list) {
    ref_sum_cost_diff += block_stats.inter_cost - block_stats.intra_cost;
  }
  return ref_sum_cost_diff;
}

double tpl_frame_dep_stats_accumulate(const TplFrameDepStats &frame_dep_stats) {
  double sum = 0;
  for (const auto &row : frame_dep_stats.unit_stats) {
    sum = std::accumulate(row.begin(), row.end(), sum);
  }
  return sum;
}

// This is a generalization of GET_MV_RAWPEL that allows for an arbitrary number
// of fractional bits.
// TODO(angiebird): Add unit test to this function
int get_fullpel_value(int subpel_value, int subpel_bits) {
  const int subpel_scale = (1 << subpel_bits);
  const int sign = subpel_value >= 0 ? 1 : -1;
  int fullpel_value = (abs(subpel_value) + subpel_scale / 2) >> subpel_bits;
  fullpel_value *= sign;
  return fullpel_value;
}

void tpl_frame_dep_stats_propagate(const TplFrameStats &frame_stats,
                                   const RefFrameTable &ref_frame_table,
                                   TplGopDepStats *tpl_gop_dep_stats) {
  const int min_block_size = frame_stats.min_block_size;
  const int frame_unit_rows =
      frame_stats.frame_height / frame_stats.min_block_size;
  const int frame_unit_cols =
      frame_stats.frame_width / frame_stats.min_block_size;
  for (const TplBlockStats &block_stats : frame_stats.block_stats_list) {
    int ref_coding_idx_list[kBlockRefCount] = { -1, -1 };
    int ref_frame_count = get_ref_coding_idx_list(block_stats, ref_frame_table,
                                                  ref_coding_idx_list);
    if (ref_frame_count > 0) {
      double propagation_ratio = 1.0 / ref_frame_count;
      for (int i = 0; i < kBlockRefCount; ++i) {
        if (ref_coding_idx_list[i] != -1) {
          auto &ref_frame_dep_stats =
              tpl_gop_dep_stats->frame_dep_stats_list[ref_coding_idx_list[i]];
          const auto &mv = block_stats.mv[i];
          const int mv_row = get_fullpel_value(mv.row, mv.subpel_bits);
          const int mv_col = get_fullpel_value(mv.col, mv.subpel_bits);
          const int block_unit_rows = block_stats.height / min_block_size;
          const int block_unit_cols = block_stats.width / min_block_size;
          const int unit_count = block_unit_rows * block_unit_cols;
          const double cost_diff =
              (block_stats.inter_cost - block_stats.intra_cost) * 1.0 /
              unit_count;
          for (int r = 0; r < block_unit_rows; r++) {
            for (int c = 0; c < block_unit_cols; c++) {
              const int ref_block_row =
                  block_stats.row + r * min_block_size + mv_row;
              const int ref_block_col =
                  block_stats.col + c * min_block_size + mv_col;
              const int ref_unit_row_low = ref_block_row / min_block_size;
              const int ref_unit_col_low = ref_block_col / min_block_size;
              for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                  const int unit_row = ref_unit_row_low + j;
                  const int unit_col = ref_unit_col_low + k;
                  if (unit_row >= 0 && unit_row < frame_unit_rows &&
                      unit_col >= 0 && unit_col < frame_unit_cols) {
                    const int overlap_area = get_block_overlap_area(
                        unit_row * min_block_size, unit_col * min_block_size,
                        ref_block_row, ref_block_col, min_block_size);
                    const double overlap_ratio =
                        overlap_area * 1.0 / (min_block_size * min_block_size);
                    ref_frame_dep_stats.unit_stats[unit_row][unit_col] +=
                        cost_diff * overlap_ratio * propagation_ratio;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

// TODO(angiebird): Add unit test for this function
std::vector<RefFrameTable> get_ref_frame_table_list(
    const GopStruct &gop_struct, RefFrameTable ref_frame_table) {
  const int frame_count = static_cast<int>(gop_struct.gop_frame_list.size());
  std::vector<RefFrameTable> ref_frame_table_list;
  ref_frame_table_list.push_back(ref_frame_table);
  for (int coding_idx = 0; coding_idx < frame_count; coding_idx++) {
    const auto &gop_frame = gop_struct.gop_frame_list[coding_idx];
    if (gop_frame.update_ref_idx != -1) {
      ref_frame_table[gop_frame.update_ref_idx] = gop_frame;
    }
    ref_frame_table_list.push_back(ref_frame_table);
  }
  return ref_frame_table_list;
}

TplGopDepStats compute_tpl_gop_dep_stats(
    const TplGopStats &tpl_gop_stats,
    const std::vector<RefFrameTable> &ref_frame_table_list) {
  const int frame_count = static_cast<int>(ref_frame_table_list.size());

  // Create the struct to store TPL dependency stats
  TplGopDepStats tpl_gop_dep_stats;
  for (int coding_idx = 0; coding_idx < frame_count; coding_idx++) {
    tpl_gop_dep_stats.frame_dep_stats_list.push_back(
        create_tpl_frame_dep_stats_wo_propagation(
            tpl_gop_stats.frame_stats_list[coding_idx]));
  }

  // Back propagation
  for (int coding_idx = frame_count - 1; coding_idx >= 0; coding_idx--) {
    auto &ref_frame_table = ref_frame_table_list[coding_idx];
    // TODO(angiebird): Handle/test the case where reference frame
    // is in the previous GOP
    tpl_frame_dep_stats_propagate(tpl_gop_stats.frame_stats_list[coding_idx],
                                  ref_frame_table, &tpl_gop_dep_stats);
  }
  return tpl_gop_dep_stats;
}

GopEncodeInfo AV1RateControlQMode::GetGopEncodeInfo(
    const GopStruct &gop_struct, const TplGopStats &tpl_gop_stats,
    const RefFrameTable &ref_frame_table_snapshot_init) {
  const std::vector<RefFrameTable> ref_frame_table_list =
      get_ref_frame_table_list(gop_struct, ref_frame_table_snapshot_init);

  GopEncodeInfo gop_encode_info;
  gop_encode_info.final_snapshot = ref_frame_table_list.back();
  TplGopDepStats gop_dep_stats =
      compute_tpl_gop_dep_stats(tpl_gop_stats, ref_frame_table_list);
  const int frame_count =
      static_cast<int>(tpl_gop_stats.frame_stats_list.size());
  for (int i = 0; i < frame_count; i++) {
    const TplFrameStats &frame_stats = tpl_gop_stats.frame_stats_list[i];
    const TplFrameDepStats &frame_dep_stats =
        gop_dep_stats.frame_dep_stats_list[i];
    const double cost_without_propagation =
        tpl_frame_stats_accumulate(frame_stats);
    const double cost_with_propagation =
        tpl_frame_dep_stats_accumulate(frame_dep_stats);
    // TODO(angiebird): This part is still a draft. Check whether this makes
    // sense mathmatically.
    const double frame_importance =
        cost_with_propagation / cost_without_propagation;
    // Imitate the behavior of av1_tpl_get_qstep_ratio()
    const double qstep_ratio = sqrt(1 / frame_importance);
    FrameEncodeParameters param;
    param.q_index = av1_get_q_index_from_qstep_ratio(rc_param_.base_q_index,
                                                     qstep_ratio, AOM_BITS_8);
    // TODO(angiebird): Determine rdmult based on q_index
    param.rdmult = 1;
    gop_encode_info.param_list.push_back(param);
  }
  return gop_encode_info;
}

}  // namespace aom
