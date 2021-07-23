/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdio.h>
#include <stdlib.h>

#include "av1/common/common.h"
#include "av1/encoder/firstpass.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {
TEST(FirstpassTest, FirstpassInfoInitWithExtBuf) {
  FIRSTPASS_INFO firstpass_info;
  FIRSTPASS_STATS ext_stats_buf[10];
  const int ref_stats_size = 10;
  for (int i = 0; i < ref_stats_size; ++i) {
    av1_zero(ext_stats_buf[i]);
    ext_stats_buf[i].frame = i;
  }
  aom_codec_err_t ret =
      av1_firstpass_info_init(&firstpass_info, ext_stats_buf, 10);
  EXPECT_EQ(firstpass_info.stats_count, ref_stats_size);
  EXPECT_EQ(ret, AOM_CODEC_OK);
}

TEST(FirstpassTest, FirstpassInfoInitWithStaticBuf) {
  FIRSTPASS_INFO firstpass_info;
  aom_codec_err_t ret = av1_firstpass_info_init(&firstpass_info, NULL, 0);
  EXPECT_EQ(firstpass_info.stats_count, 0);
  EXPECT_EQ(ret, AOM_CODEC_OK);
}

TEST(FirstpassTest, FirstpassInfoPushPop) {
  FIRSTPASS_INFO firstpass_info;
  av1_firstpass_info_init(&firstpass_info, NULL, 0);
  EXPECT_EQ(firstpass_info.stats_buf_size, 48);
  for (int i = 0; i < 48; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    stats.frame = i;
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, 48);
  for (int i = 0; i < 20; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    aom_codec_err_t ret = av1_firstpass_info_pop(&firstpass_info, &stats);
    EXPECT_EQ(stats.frame, i);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, 28);

  for (int i = 0; i < 20; ++i) {
    FIRSTPASS_STATS stats;
    av1_zero(stats);
    aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
    EXPECT_EQ(ret, AOM_CODEC_OK);
  }
  EXPECT_EQ(firstpass_info.stats_count, 48);

  EXPECT_EQ(firstpass_info.stats_count, firstpass_info.stats_buf_size);
  // Pusht a stats when the queue is full.
  FIRSTPASS_STATS stats;
  av1_zero(stats);
  aom_codec_err_t ret = av1_firstpass_info_push(&firstpass_info, &stats);
  EXPECT_EQ(ret, AOM_CODEC_ERROR);
}
}  // namespace
