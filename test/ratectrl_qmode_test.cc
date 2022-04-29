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

#include <array>
#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "av1/ratectrl_qmode.h"
#include "av1/reference_manager.h"
#include "test/mock_ratectrl_qmode.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace aom {

using ::testing::ElementsAre;
using ::testing::Field;
using ::testing::Return;

void test_gop_display_order(const GopStruct &gop_struct) {
  // Test whether show frames' order indices are sequential
  int expected_order_idx = 0;
  int expected_show_frame_count = 0;
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    if (gop_frame.is_show_frame) {
      EXPECT_EQ(gop_frame.order_idx, expected_order_idx);
      expected_order_idx++;
      expected_show_frame_count++;
    }
  }
  EXPECT_EQ(gop_struct.show_frame_count, expected_show_frame_count);
}

void test_gop_global_order_idx(const GopStruct &gop_struct,
                               int global_order_idx_offset) {
  // Test whether show frames' global order indices are sequential
  EXPECT_EQ(gop_struct.global_order_idx_offset, global_order_idx_offset);
  int expected_global_order_idx = global_order_idx_offset;
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    if (gop_frame.is_show_frame) {
      EXPECT_EQ(gop_frame.global_order_idx, expected_global_order_idx);
      expected_global_order_idx++;
    }
  }
}

void test_gop_global_coding_idx(const GopStruct &gop_struct,
                                int global_coding_idx_offset) {
  EXPECT_EQ(gop_struct.global_coding_idx_offset, global_coding_idx_offset);
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    EXPECT_EQ(gop_frame.global_coding_idx,
              global_coding_idx_offset + gop_frame.coding_idx);
  }
}

void test_colocated_show_frame(const GopStruct &gop_struct) {
  // Test whether each non show frame has a colocated show frame
  int gop_size = static_cast<int>(gop_struct.gop_frame_list.size());
  for (int gop_idx = 0; gop_idx < gop_size; ++gop_idx) {
    auto &gop_frame = gop_struct.gop_frame_list[gop_idx];
    if (gop_frame.is_show_frame == 0) {
      bool found_colocated_ref_frame = false;
      for (int i = gop_idx + 1; i < gop_size; ++i) {
        auto &next_gop_frame = gop_struct.gop_frame_list[i];
        if (gop_frame.order_idx == next_gop_frame.order_idx) {
          found_colocated_ref_frame = true;
          EXPECT_EQ(gop_frame.update_ref_idx, next_gop_frame.colocated_ref_idx);
          EXPECT_TRUE(next_gop_frame.is_show_frame);
        }
        if (gop_frame.update_ref_idx == next_gop_frame.update_ref_idx) {
          break;
        }
      }
      EXPECT_TRUE(found_colocated_ref_frame);
    }
  }
}

void test_layer_depth(const GopStruct &gop_struct, int max_layer_depth) {
  int gop_size = static_cast<int>(gop_struct.gop_frame_list.size());
  for (int gop_idx = 0; gop_idx < gop_size; ++gop_idx) {
    const auto &gop_frame = gop_struct.gop_frame_list[gop_idx];
    if (gop_frame.is_key_frame) {
      EXPECT_EQ(gop_frame.layer_depth, 0);
    }

    if (gop_frame.is_arf_frame) {
      EXPECT_LT(gop_frame.layer_depth, max_layer_depth);
    }

    if (!gop_frame.is_key_frame && !gop_frame.is_arf_frame) {
      EXPECT_EQ(gop_frame.layer_depth, max_layer_depth);
    }
  }
}

void test_arf_interval(const GopStruct &gop_struct) {
  std::vector<int> arf_order_idx_list;
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    if (gop_frame.is_arf_frame) {
      arf_order_idx_list.push_back(gop_frame.order_idx);
    }
  }
  std::sort(arf_order_idx_list.begin(), arf_order_idx_list.end());
  int arf_count = static_cast<int>(arf_order_idx_list.size());
  for (int i = 1; i < arf_count; ++i) {
    int arf_interval = arf_order_idx_list[i] - arf_order_idx_list[i - 1];
    EXPECT_GE(arf_interval, kMinArfInterval);
  }
}

TEST(RateControlQModeTest, ConstructGopARF) {
  int show_frame_count = 16;
  const bool has_key_frame = false;
  const int global_coding_idx_offset = 5;
  const int global_order_idx_offset = 20;
  RefFrameManager ref_frame_manager(kRefFrameTableSize);
  GopStruct gop_struct =
      construct_gop(&ref_frame_manager, show_frame_count, has_key_frame,
                    global_coding_idx_offset, global_order_idx_offset);
  EXPECT_EQ(gop_struct.show_frame_count, show_frame_count);
  test_gop_display_order(gop_struct);
  test_gop_global_order_idx(gop_struct, global_order_idx_offset);
  test_gop_global_coding_idx(gop_struct, global_coding_idx_offset);
  test_colocated_show_frame(gop_struct);
  const int max_layer_depth =
      ref_frame_manager.ForwardMaxSize() + kLayerDepthOffset;
  test_layer_depth(gop_struct, max_layer_depth);
  test_arf_interval(gop_struct);
}

TEST(RateControlQModeTest, ConstructGopKey) {
  const int show_frame_count = 16;
  const int has_key_frame = 1;
  const int global_coding_idx_offset = 10;
  const int global_order_idx_offset = 8;
  RefFrameManager ref_frame_manager(kRefFrameTableSize);
  GopStruct gop_struct =
      construct_gop(&ref_frame_manager, show_frame_count, has_key_frame,
                    global_coding_idx_offset, global_order_idx_offset);
  EXPECT_EQ(gop_struct.show_frame_count, show_frame_count);
  test_gop_display_order(gop_struct);
  test_gop_global_order_idx(gop_struct, global_order_idx_offset);
  test_gop_global_coding_idx(gop_struct, global_coding_idx_offset);
  test_colocated_show_frame(gop_struct);
  const int max_layer_depth =
      ref_frame_manager.ForwardMaxSize() + kLayerDepthOffset;
  test_layer_depth(gop_struct, max_layer_depth);
  test_arf_interval(gop_struct);
}

static TplBlockStats create_toy_tpl_block_stats(int h, int w, int r, int c,
                                                int cost_diff) {
  TplBlockStats tpl_block_stats = {};
  tpl_block_stats.height = h;
  tpl_block_stats.width = w;
  tpl_block_stats.row = r;
  tpl_block_stats.col = c;
  // A random trick that makes inter_cost - intra_cost = cost_diff;
  tpl_block_stats.intra_cost = cost_diff / 2;
  tpl_block_stats.inter_cost = cost_diff + cost_diff / 2;
  tpl_block_stats.ref_frame_index = { -1, -1 };
  return tpl_block_stats;
}

static TplFrameStats create_toy_tpl_frame_stats_with_diff_sizes(
    int min_block_size, int max_block_size) {
  TplFrameStats frame_stats;
  const int max_h = max_block_size;
  const int max_w = max_h;
  const int count = max_block_size / min_block_size;
  frame_stats.min_block_size = min_block_size;
  frame_stats.frame_height = max_h * count;
  frame_stats.frame_width = max_w * count;
  for (int i = 0; i < count; ++i) {
    for (int j = 0; j < count; ++j) {
      int h = max_h >> i;
      int w = max_w >> j;
      for (int u = 0; u * h < max_h; ++u) {
        for (int v = 0; v * w < max_w; ++v) {
          int r = max_h * i + h * u;
          int c = max_w * j + w * v;
          int cost_diff = std::rand() % 16;
          TplBlockStats block_stats =
              create_toy_tpl_block_stats(h, w, r, c, cost_diff);
          frame_stats.block_stats_list.push_back(block_stats);
        }
      }
    }
  }
  return frame_stats;
}

static void augment_tpl_frame_stats_with_ref_frames(
    TplFrameStats *tpl_frame_stats,
    const std::array<int, kBlockRefCount> &ref_frame_index) {
  for (auto &block_stats : tpl_frame_stats->block_stats_list) {
    block_stats.ref_frame_index = ref_frame_index;
  }
}
static void augment_tpl_frame_stats_with_motion_vector(
    TplFrameStats *tpl_frame_stats,
    const std::array<MotionVector, kBlockRefCount> &mv) {
  for (auto &block_stats : tpl_frame_stats->block_stats_list) {
    block_stats.mv = mv;
  }
}

static RefFrameTable create_toy_ref_frame_table(int frame_count) {
  RefFrameTable ref_frame_table;
  const int ref_frame_table_size = static_cast<int>(ref_frame_table.size());
  EXPECT_LE(frame_count, ref_frame_table_size);
  for (int i = 0; i < frame_count; ++i) {
    ref_frame_table[i] =
        gop_frame_basic(0, 0, i, i, 0, GopFrameType::kRegularLeaf);
  }
  for (int i = frame_count; i < ref_frame_table_size; ++i) {
    ref_frame_table[i] = gop_frame_invalid();
  }
  return ref_frame_table;
}

static MotionVector create_fullpel_mv(int row, int col) {
  return { row, col, 0 };
}

TEST(RateControlQModeTest, CreateTplFrameDepStats) {
  TplFrameStats frame_stats = create_toy_tpl_frame_stats_with_diff_sizes(8, 16);
  TplFrameDepStats frame_dep_stats =
      create_tpl_frame_dep_stats_wo_propagation(frame_stats);
  EXPECT_EQ(frame_stats.min_block_size, frame_dep_stats.unit_size);
  const int unit_rows = static_cast<int>(frame_dep_stats.unit_stats.size());
  const int unit_cols = static_cast<int>(frame_dep_stats.unit_stats[0].size());
  EXPECT_EQ(frame_stats.frame_height, unit_rows * frame_dep_stats.unit_size);
  EXPECT_EQ(frame_stats.frame_width, unit_cols * frame_dep_stats.unit_size);
  const double sum_cost_diff = tpl_frame_dep_stats_accumulate(frame_dep_stats);

  const double ref_sum_cost_diff = tpl_frame_stats_accumulate(frame_stats);
  EXPECT_NEAR(sum_cost_diff, ref_sum_cost_diff, 0.0000001);
}

TEST(RateControlQModeTest, GetBlockOverlapArea) {
  const int size = 8;
  const int r0 = 8;
  const int c0 = 9;
  std::vector<int> r1 = { 8, 10, 16, 10, 8, 100 };
  std::vector<int> c1 = { 9, 12, 17, 5, 100, 9 };
  std::vector<int> ref_overlap = { 64, 30, 0, 24, 0, 0 };
  for (int i = 0; i < static_cast<int>(r1.size()); ++i) {
    const int overlap0 = get_block_overlap_area(r0, c0, r1[i], c1[i], size);
    const int overlap1 = get_block_overlap_area(r1[i], c1[i], r0, c0, size);
    EXPECT_EQ(overlap0, ref_overlap[i]);
    EXPECT_EQ(overlap1, ref_overlap[i]);
  }
}

TEST(RateControlQModeTest, TplFrameDepStatsPropagateSingleZeroMotion) {
  // cur frame with coding_idx 1 use ref frame with coding_idx 0
  const std::array<int, kBlockRefCount> ref_frame_index = { 0, -1 };
  TplFrameStats frame_stats = create_toy_tpl_frame_stats_with_diff_sizes(8, 16);
  augment_tpl_frame_stats_with_ref_frames(&frame_stats, ref_frame_index);

  TplGopDepStats gop_dep_stats;
  const int frame_count = 2;
  // ref frame with coding_idx 0
  TplFrameDepStats frame_dep_stats0 = create_tpl_frame_dep_stats_empty(
      frame_stats.frame_height, frame_stats.frame_width,
      frame_stats.min_block_size);
  gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats0);

  // cur frame with coding_idx 1
  const TplFrameDepStats frame_dep_stats1 =
      create_tpl_frame_dep_stats_wo_propagation(frame_stats);
  gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats1);

  const RefFrameTable ref_frame_table = create_toy_ref_frame_table(frame_count);
  tpl_frame_dep_stats_propagate(frame_stats, ref_frame_table, &gop_dep_stats);

  // cur frame with coding_idx 1
  const double ref_sum_cost_diff = tpl_frame_stats_accumulate(frame_stats);

  // ref frame with coding_idx 0
  const double sum_cost_diff =
      tpl_frame_dep_stats_accumulate(gop_dep_stats.frame_dep_stats_list[0]);

  // The sum_cost_diff between coding_idx 0 and coding_idx 1 should be equal
  // because every block in cur frame has zero motion, use ref frame with
  // coding_idx 0 for prediction, and ref frame itself is empty.
  EXPECT_NEAR(sum_cost_diff, ref_sum_cost_diff, 0.0000001);
}

TEST(RateControlQModeTest, TplFrameDepStatsPropagateCompoundZeroMotion) {
  // cur frame with coding_idx 2 use two ref frames with coding_idx 0 and 1
  const std::array<int, kBlockRefCount> ref_frame_index = { 0, 1 };
  TplFrameStats frame_stats = create_toy_tpl_frame_stats_with_diff_sizes(8, 16);
  augment_tpl_frame_stats_with_ref_frames(&frame_stats, ref_frame_index);

  TplGopDepStats gop_dep_stats;
  const int frame_count = 3;
  // ref frame with coding_idx 0
  const TplFrameDepStats frame_dep_stats0 = create_tpl_frame_dep_stats_empty(
      frame_stats.frame_height, frame_stats.frame_width,
      frame_stats.min_block_size);
  gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats0);

  // ref frame with coding_idx 1
  const TplFrameDepStats frame_dep_stats1 = create_tpl_frame_dep_stats_empty(
      frame_stats.frame_height, frame_stats.frame_width,
      frame_stats.min_block_size);
  gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats1);

  // cur frame with coding_idx 2
  const TplFrameDepStats frame_dep_stats2 =
      create_tpl_frame_dep_stats_wo_propagation(frame_stats);
  gop_dep_stats.frame_dep_stats_list.push_back(frame_dep_stats2);

  const RefFrameTable ref_frame_table = create_toy_ref_frame_table(frame_count);
  tpl_frame_dep_stats_propagate(frame_stats, ref_frame_table, &gop_dep_stats);

  // cur frame with coding_idx 1
  const double ref_sum_cost_diff = tpl_frame_stats_accumulate(frame_stats);

  // ref frame with coding_idx 0
  const double sum_cost_diff0 =
      tpl_frame_dep_stats_accumulate(gop_dep_stats.frame_dep_stats_list[0]);
  EXPECT_NEAR(sum_cost_diff0, ref_sum_cost_diff * 0.5, 0.0000001);

  // ref frame with coding_idx 1
  const double sum_cost_diff1 =
      tpl_frame_dep_stats_accumulate(gop_dep_stats.frame_dep_stats_list[1]);
  EXPECT_NEAR(sum_cost_diff1, ref_sum_cost_diff * 0.5, 0.0000001);
}

TEST(RateControlQModeTest, TplFrameDepStatsPropagateSingleWithMotion) {
  // cur frame with coding_idx 1 use ref frame with coding_idx 0
  const std::array<int, kBlockRefCount> ref_frame_index = { 0, -1 };
  const int min_block_size = 8;
  TplFrameStats frame_stats = create_toy_tpl_frame_stats_with_diff_sizes(
      min_block_size, min_block_size * 2);
  augment_tpl_frame_stats_with_ref_frames(&frame_stats, ref_frame_index);

  const int mv_row = min_block_size / 2;
  const int mv_col = min_block_size / 4;
  const double r_ratio = 1.0 / 2;
  const double c_ratio = 1.0 / 4;
  std::array<MotionVector, kBlockRefCount> mv;
  mv[0] = create_fullpel_mv(mv_row, mv_col);
  mv[1] = create_fullpel_mv(0, 0);
  augment_tpl_frame_stats_with_motion_vector(&frame_stats, mv);

  TplGopDepStats gop_dep_stats;
  const int frame_count = 2;
  // ref frame with coding_idx 0
  gop_dep_stats.frame_dep_stats_list.push_back(create_tpl_frame_dep_stats_empty(
      frame_stats.frame_height, frame_stats.frame_width,
      frame_stats.min_block_size));

  // cur frame with coding_idx 1
  gop_dep_stats.frame_dep_stats_list.push_back(
      create_tpl_frame_dep_stats_wo_propagation(frame_stats));

  const RefFrameTable ref_frame_table = create_toy_ref_frame_table(frame_count);
  tpl_frame_dep_stats_propagate(frame_stats, ref_frame_table, &gop_dep_stats);

  const auto &dep_stats0 = gop_dep_stats.frame_dep_stats_list[0];
  const auto &dep_stats1 = gop_dep_stats.frame_dep_stats_list[1];
  const int unit_rows = static_cast<int>(dep_stats0.unit_stats.size());
  const int unit_cols = static_cast<int>(dep_stats0.unit_stats[0].size());
  for (int r = 0; r < unit_rows; ++r) {
    for (int c = 0; c < unit_cols; ++c) {
      double ref_value = 0;
      ref_value += (1 - r_ratio) * (1 - c_ratio) * dep_stats1.unit_stats[r][c];
      if (r - 1 >= 0) {
        ref_value += r_ratio * (1 - c_ratio) * dep_stats1.unit_stats[r - 1][c];
      }
      if (c - 1 >= 0) {
        ref_value += (1 - r_ratio) * c_ratio * dep_stats1.unit_stats[r][c - 1];
      }
      if (r - 1 >= 0 && c - 1 >= 0) {
        ref_value += r_ratio * c_ratio * dep_stats1.unit_stats[r - 1][c - 1];
      }
      EXPECT_NEAR(dep_stats0.unit_stats[r][c], ref_value, 0.0000001);
    }
  }
}

TEST(RateControlQModeTest, ComputeTplGopDepStats) {
  TplGopStats tpl_gop_stats;
  std::vector<RefFrameTable> ref_frame_table_list;
  for (int i = 0; i < 3; i++) {
    // Use the previous frame as reference
    const std::array<int, kBlockRefCount> ref_frame_index = { i - 1, -1 };
    int min_block_size = 8;
    TplFrameStats frame_stats = create_toy_tpl_frame_stats_with_diff_sizes(
        min_block_size, min_block_size * 2);
    augment_tpl_frame_stats_with_ref_frames(&frame_stats, ref_frame_index);
    tpl_gop_stats.frame_stats_list.push_back(frame_stats);

    ref_frame_table_list.push_back(create_toy_ref_frame_table(i));
  }
  const TplGopDepStats &gop_dep_stats =
      compute_tpl_gop_dep_stats(tpl_gop_stats, ref_frame_table_list);

  double ref_sum = 0;
  for (int i = 2; i >= 0; i--) {
    // Due to the linear propagation with zero motion, we can add the
    // frame_stats value and use it as reference sum for dependency stats
    ref_sum += tpl_frame_stats_accumulate(tpl_gop_stats.frame_stats_list[i]);
    const double sum =
        tpl_frame_dep_stats_accumulate(gop_dep_stats.frame_dep_stats_list[i]);
    EXPECT_NEAR(sum, ref_sum, 0.0000001);
    break;
  }
}

TEST(RefFrameManagerTest, GetRefFrameCount) {
  const std::vector<int> order_idx_list = { 0, 4, 2, 1, 2, 3, 4 };
  const std::vector<GopFrameType> type_list = {
    GopFrameType::kRegularKey,      GopFrameType::kRegularArf,
    GopFrameType::kIntermediateArf, GopFrameType::kRegularLeaf,
    GopFrameType::kShowExisting,    GopFrameType::kRegularLeaf,
    GopFrameType::kOverlay
  };
  RefFrameManager ref_manager(kRefFrameTableSize);
  int coding_idx = 0;
  const int first_leaf_idx = 3;
  EXPECT_EQ(type_list[first_leaf_idx], GopFrameType::kRegularLeaf);
  // update reference frame until we see the first kRegularLeaf frame
  for (; coding_idx <= first_leaf_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }
  EXPECT_EQ(ref_manager.GetRefFrameCount(), 4);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 2);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 1);
  EXPECT_EQ(ref_manager.CurGlobalOrderIdx(), 1);

  // update reference frame until we see the first kShowExisting frame
  const int first_show_existing_idx = 4;
  EXPECT_EQ(type_list[first_show_existing_idx], GopFrameType::kShowExisting);
  for (; coding_idx <= first_show_existing_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }
  EXPECT_EQ(ref_manager.GetRefFrameCount(), 4);
  EXPECT_EQ(ref_manager.CurGlobalOrderIdx(), 2);
  // After the first kShowExisting, the kIntermediateArf should be moved from
  // kForward to kLast due to the cur_global_order_idx_ update
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 2);

  const int second_leaf_idx = 5;
  EXPECT_EQ(type_list[second_leaf_idx], GopFrameType::kRegularLeaf);
  for (; coding_idx <= second_leaf_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }
  EXPECT_EQ(ref_manager.GetRefFrameCount(), 5);
  EXPECT_EQ(ref_manager.CurGlobalOrderIdx(), 3);
  // An additional kRegularLeaf frame is added into kLast
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 3);

  const int first_overlay_idx = 6;
  EXPECT_EQ(type_list[first_overlay_idx], GopFrameType::kOverlay);
  for (; coding_idx <= first_overlay_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }

  EXPECT_EQ(ref_manager.GetRefFrameCount(), 5);
  EXPECT_EQ(ref_manager.CurGlobalOrderIdx(), 4);
  // After the kOverlay, the kRegularArf should be moved from
  // kForward to kBackward due to the cur_global_order_idx_ update
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 0);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 2);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 3);
}

void test_ref_frame_manager_priority(const RefFrameManager &ref_manager,
                                     RefUpdateType type) {
  int ref_count = ref_manager.GetRefFrameCountByType(type);
  int prev_global_order_idx = ref_manager.CurGlobalOrderIdx();
  // The lower the priority is, the closer the gop_frame.global_order_idx should
  // be with cur_global_order_idx_
  for (int priority = 0; priority < ref_count; ++priority) {
    GopFrame gop_frame = ref_manager.GetRefFrameByPriority(type, priority);
    EXPECT_EQ(gop_frame.is_valid, true);
    if (type == RefUpdateType::kForward) {
      EXPECT_GE(gop_frame.global_order_idx, prev_global_order_idx);
    } else {
      EXPECT_LE(gop_frame.global_order_idx, prev_global_order_idx);
    }
    prev_global_order_idx = gop_frame.global_order_idx;
  }
  GopFrame gop_frame =
      ref_manager.GetRefFrameByPriority(RefUpdateType::kForward, ref_count);
  EXPECT_EQ(gop_frame.is_valid, false);
}

TEST(RefFrameManagerTest, GetRefFrameByPriority) {
  const std::vector<int> order_idx_list = { 0, 4, 2, 1, 2, 3, 4 };
  const std::vector<GopFrameType> type_list = {
    GopFrameType::kRegularKey,      GopFrameType::kRegularArf,
    GopFrameType::kIntermediateArf, GopFrameType::kRegularLeaf,
    GopFrameType::kShowExisting,    GopFrameType::kRegularLeaf,
    GopFrameType::kOverlay
  };
  RefFrameManager ref_manager(kRefFrameTableSize);
  int coding_idx = 0;
  const int first_leaf_idx = 3;
  EXPECT_EQ(type_list[first_leaf_idx], GopFrameType::kRegularLeaf);
  // update reference frame until we see the first kRegularLeaf frame
  for (; coding_idx <= first_leaf_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 2);
  test_ref_frame_manager_priority(ref_manager, RefUpdateType::kForward);

  const int first_overlay_idx = 6;
  EXPECT_EQ(type_list[first_overlay_idx], GopFrameType::kOverlay);
  for (; coding_idx <= first_overlay_idx; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }

  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 2);
  test_ref_frame_manager_priority(ref_manager, RefUpdateType::kBackward);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 3);
  test_ref_frame_manager_priority(ref_manager, RefUpdateType::kLast);
}

TEST(RefFrameManagerTest, GetRefFrameListByPriority) {
  const std::vector<int> order_idx_list = { 0, 4, 2, 1 };
  const int frame_count = static_cast<int>(order_idx_list.size());
  const std::vector<GopFrameType> type_list = { GopFrameType::kRegularKey,
                                                GopFrameType::kRegularArf,
                                                GopFrameType::kIntermediateArf,
                                                GopFrameType::kRegularLeaf };
  RefFrameManager ref_manager(kRefFrameTableSize);
  for (int coding_idx = 0; coding_idx < frame_count; ++coding_idx) {
    GopFrame gop_frame = gop_frame_basic(
        0, 0, coding_idx, order_idx_list[coding_idx], 0, type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }
  EXPECT_EQ(ref_manager.GetRefFrameCount(), frame_count);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kForward), 2);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kBackward), 1);
  EXPECT_EQ(ref_manager.GetRefFrameCountByType(RefUpdateType::kLast), 1);
  std::vector<ReferenceFrame> ref_frame_list =
      ref_manager.GetRefFrameListByPriority();
  EXPECT_EQ(ref_frame_list.size(), order_idx_list.size());
  std::vector<int> expected_global_order_idx = { 2, 0, 1, 4 };
  std::vector<ReferenceName> expected_names = { ReferenceName::kBwdrefFrame,
                                                ReferenceName::kGoldenFrame,
                                                ReferenceName::kLastFrame,
                                                ReferenceName::kAltref2Frame };
  for (size_t i = 0; i < ref_frame_list.size(); ++i) {
    ReferenceFrame &ref_frame = ref_frame_list[i];
    GopFrame gop_frame = ref_manager.GetRefFrameByIndex(ref_frame.index);
    EXPECT_EQ(gop_frame.global_order_idx, expected_global_order_idx[i]);
    EXPECT_EQ(ref_frame.name, expected_names[i]);
  }
}

TEST(RefFrameManagerTest, GetPrimaryRefFrame) {
  const std::vector<int> order_idx_list = { 0, 4, 2, 1 };
  const int frame_count = static_cast<int>(order_idx_list.size());
  const std::vector<GopFrameType> type_list = { GopFrameType::kRegularKey,
                                                GopFrameType::kRegularArf,
                                                GopFrameType::kIntermediateArf,
                                                GopFrameType::kRegularLeaf };
  const std::vector<int> layer_depth_list = { 0, 2, 4, 6 };
  RefFrameManager ref_manager(kRefFrameTableSize);
  for (int coding_idx = 0; coding_idx < frame_count; ++coding_idx) {
    GopFrame gop_frame =
        gop_frame_basic(0, 0, coding_idx, order_idx_list[coding_idx],
                        layer_depth_list[coding_idx], type_list[coding_idx]);
    ref_manager.UpdateRefFrameTable(&gop_frame);
  }

  for (int i = 0; i < frame_count; ++i) {
    // Test frame that share the same layer depth with a reference frame
    int layer_depth = layer_depth_list[i];
    // Set different frame type
    GopFrameType type = type_list[(i + 1) % frame_count];
    GopFrame gop_frame = gop_frame_basic(0, 0, 0, 0, layer_depth, type);
    ReferenceFrame ref_frame = ref_manager.GetPrimaryRefFrame(gop_frame);
    GopFrame primary_ref_frame =
        ref_manager.GetRefFrameByIndex(ref_frame.index);
    // The GetPrimaryRefFrame should find the ref_frame with matched layer depth
    // because it's our first priority
    EXPECT_EQ(primary_ref_frame.layer_depth, gop_frame.layer_depth);
  }

  const std::vector<int> mid_layer_depth_list = { 1, 3, 5 };
  for (int i = 0; i < 3; ++i) {
    // Test frame that share the same frame type with a reference frame
    GopFrameType type = type_list[i];
    // Let the frame layer_depth sit in the middle of two reference frames
    int layer_depth = mid_layer_depth_list[i];
    GopFrame gop_frame = gop_frame_basic(0, 0, 0, 0, layer_depth, type);
    ReferenceFrame ref_frame = ref_manager.GetPrimaryRefFrame(gop_frame);
    GopFrame primary_ref_frame =
        ref_manager.GetRefFrameByIndex(ref_frame.index);
    // The GetPrimaryRefFrame should find the ref_frame with matched frame type
    // Here we use coding_idx to confirm that.
    EXPECT_EQ(primary_ref_frame.coding_idx, i);
  }
}

// Reads a whitespace-delimited string from stream, and parses it as a double.
// Returns an empty string if the entire string was successfully parsed as a
// double, or an error messaage if not.
std::string ReadDouble(std::istream &stream, double *value) {
  std::string word;
  stream >> word;
  if (word.empty()) {
    return "Unexpectedly reached end of input";
  }
  char *end;
  *value = std::strtod(word.c_str(), &end);
  if (*end != '\0') {
    return "Unexpected characters found: " + word;
  }
  return "";
}

TEST(RateControlQModeTest, TestKeyframeDetection) {
  FirstpassInfo firstpass_info;
  // These golden files are generated by the following command line:
  // ./aomenc --width=352 --height=288 --fps=30/1 --limit=250 --codec=av1
  // --cpu-used=3 --end-usage=q --cq-level=36 --threads=0 --profile=0
  // --lag-in-frames=35 --min-q=0 --max-q=63 --auto-alt-ref=1 --passes=2
  // --kf-max-dist=160 --kf-min-dist=0 --drop-frame=0
  // --static-thresh=0 --minsection-pct=0 --maxsection-pct=2000
  // --arnr-maxframes=7
  // --arnr-strength=5 --sharpness=0 --undershoot-pct=100 --overshoot-pct=100
  // --frame-parallel=0
  // --tile-columns=0 -o output.webm hantro_collage_w352h288.yuv
  // First pass stats are written out in av1_get_second_pass_params right after
  // calculate_gf_length.
  std::ifstream firstpass_stats_file("firstpass_stats");
  firstpass_info.num_mbs_16x16 = (352 / 16 + 1) * (288 / 16 + 1);
  std::string newline;
  while (std::getline(firstpass_stats_file, newline)) {
    std::istringstream iss(newline);
    FIRSTPASS_STATS firstpass_stats_input = {};
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.frame), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.weight), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.intra_error), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.frame_avg_wavelet_energy),
              "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.coded_error), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.sr_coded_error), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.pcnt_inter), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.pcnt_motion), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.pcnt_second_ref), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.pcnt_neutral), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.intra_skip_pct), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.inactive_zone_rows), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.inactive_zone_cols), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.MVr), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.mvr_abs), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.MVc), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.mvc_abs), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.MVrv), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.MVcv), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.mv_in_out_count), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.new_mv_count), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.duration), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.count), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.raw_error_stdev), "");
    iss >> firstpass_stats_input.is_flash;
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.noise_var), "");
    ASSERT_EQ(ReadDouble(iss, &firstpass_stats_input.cor_coeff), "");
    ASSERT_TRUE(iss.eof()) << "Too many fields on line "
                           << firstpass_info.stats_list.size() + 1 << "\n"
                           << newline;
    firstpass_info.stats_list.push_back(firstpass_stats_input);
  }
  EXPECT_THAT(get_key_frame_list(firstpass_info),
              ElementsAre(0, 30, 60, 90, 120, 150, 180, 210, 240));
}

// MockRateControlQMode is provided for the use of clients of libaom, but it's
// not expected that it will be used in any real libaom tests.
// This simple "toy" test exists solely to verify the integration of gmock into
// the aom build.
TEST(RateControlQModeTest, TestMock) {
  MockRateControlQMode mock_rc;
  EXPECT_CALL(mock_rc,
              DetermineGopInfo(Field(&FirstpassInfo::num_mbs_16x16, 1000)))
      .WillOnce(Return(GopStructList{ { 6, 0, 0, {} }, { 4, 0, 0, {} } }));
  FirstpassInfo firstpass_info = {};
  firstpass_info.num_mbs_16x16 = 1000;
  EXPECT_THAT(mock_rc.DetermineGopInfo(firstpass_info),
              ElementsAre(Field(&GopStruct::show_frame_count, 6),
                          Field(&GopStruct::show_frame_count, 4)));
}

}  // namespace aom

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  std::srand(0);
  return RUN_ALL_TESTS();
}
