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

#include <memory>

#include "av1/ratectrl_qmode.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
namespace aom {
TEST(RateControlQModeTest, ConstructGop) {
  int show_frame_count = 16;
  int max_ref_frames = 7;
  bool has_key_frame = false;
  RefFrameManager ref_frame_manager(max_ref_frames);
  GopStruct gop_struct =
      construct_gop(&ref_frame_manager, show_frame_count, has_key_frame);
  int ref_order_idx = 0;
  // Test whether show frames' order indices are sequential
  for (const auto &gop_frame : gop_struct.gop_frame_list) {
    if (gop_frame.is_show_frame) {
      EXPECT_EQ(gop_frame.order_idx, ref_order_idx);
      ref_order_idx++;
    }
  }

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
          EXPECT_EQ(next_gop_frame.is_show_frame, 1);
        }
        if (gop_frame.update_ref_idx == next_gop_frame.update_ref_idx) {
          break;
        }
      }
      EXPECT_TRUE(found_colocated_ref_frame);
    }
  }
}

}  // namespace aom

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
