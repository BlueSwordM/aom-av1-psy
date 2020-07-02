/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom/aom_codec.h"
#include "aom_dsp/aom_dsp_common.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/y4m_video_source.h"
#include "test/util.h"

namespace {
typedef struct {
  // Superblock size
  const unsigned int sb_size;
  // log2(number of tile rows)
  const unsigned int tile_rows;
  // log2(number of tile columns)
  const unsigned int tile_cols;
} uniformTileConfigParam;

static const uniformTileConfigParam uniformTileConfigParams[] = {
  { 128, 0, 0 }, { 128, 0, 2 }, { 128, 2, 0 }, { 128, 1, 2 }, { 128, 2, 2 },
  { 128, 3, 2 }, { 64, 0, 0 },  { 64, 0, 2 },  { 64, 2, 0 },  { 64, 1, 2 },
  { 64, 2, 2 },  { 64, 3, 3 },  { 64, 4, 4 }
};

// Find smallest k>=0 such that (blk_size << k) >= target
static INLINE int tile_log2(int blk_size, int target) {
  int k;
  for (k = 0; (blk_size << k) < target; k++) {
  }
  return k;
}

// This class is used to validate tile configuration for uniform spacing.
class UniformTileConfigTestLarge
    : public ::libaom_test::CodecTestWith3Params<
          libaom_test::TestMode, uniformTileConfigParam, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  UniformTileConfigTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        tile_config_param_(GET_PARAM(2)), end_usage_check_(GET_PARAM(3)) {
    tile_config_violated_ = false;
    max_tile_cols_log2_ = tile_log2(1, AOM_MAX_TILE_COLS);
    max_tile_rows_log2_ = tile_log2(1, AOM_MAX_TILE_ROWS);
  }
  virtual ~UniformTileConfigTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = end_usage_check_;
    cfg_.g_threads = 1;
    cfg_.g_lag_in_frames = 19;
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AV1E_SET_TILE_COLUMNS, tile_config_param_.tile_cols);
      encoder->Control(AV1E_SET_TILE_ROWS, tile_config_param_.tile_rows);
      encoder->Control(AOME_SET_CPUUSED, 5);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
      encoder->Control(AV1E_SET_SUPERBLOCK_SIZE,
                       tile_config_param_.sb_size == 64
                           ? AOM_SUPERBLOCK_SIZE_64X64
                           : AOM_SUPERBLOCK_SIZE_128X128);
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      aom_tile_info tile_info;
      int config_tile_columns = AOMMIN(1 << (int)tile_config_param_.tile_cols,
                                       1 << max_tile_cols_log2_);
      int config_tile_rows = AOMMIN(1 << (int)tile_config_param_.tile_rows,
                                    1 << max_tile_rows_log2_);

      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_TILE_INFO, &tile_info);
      if (tile_info.tile_columns != config_tile_columns ||
          tile_info.tile_rows != config_tile_rows) {
        tile_config_violated_ = true;
      }
    }
    return AOM_CODEC_OK == res_dec;
  }

  ::libaom_test::TestMode encoding_mode_;
  const uniformTileConfigParam tile_config_param_;
  int max_tile_cols_log2_;
  int max_tile_rows_log2_;
  bool tile_config_violated_;
  aom_rc_mode end_usage_check_;
};

TEST_P(UniformTileConfigTestLarge, UniformTileConfigTest) {
  ::libaom_test::Y4mVideoSource video("niklas_1280_720_30.y4m", 0, 1);
  ASSERT_NO_FATAL_FAILURE(video.Begin());

  int max_tiles_cols = video.img()->w / (int)tile_config_param_.sb_size;
  int max_tiles_rows = video.img()->h / (int)tile_config_param_.sb_size;
  max_tile_cols_log2_ = tile_log2(1, AOMMIN(max_tiles_cols, AOM_MAX_TILE_COLS));
  max_tile_rows_log2_ = tile_log2(1, AOMMIN(max_tiles_rows, AOM_MAX_TILE_ROWS));

  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(tile_config_violated_, false);
}

TEST_P(UniformTileConfigTestLarge, UniformTileConfigTestLowRes) {
  ::libaom_test::Y4mVideoSource video("screendata.y4m", 0, 1);
  ASSERT_NO_FATAL_FAILURE(video.Begin());

  int max_tiles_cols = video.img()->w / (int)tile_config_param_.sb_size;
  int max_tiles_rows = video.img()->h / (int)tile_config_param_.sb_size;
  max_tile_cols_log2_ = tile_log2(1, AOMMIN(max_tiles_cols, AOM_MAX_TILE_COLS));
  max_tile_rows_log2_ = tile_log2(1, AOMMIN(max_tiles_rows, AOM_MAX_TILE_ROWS));

  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(tile_config_violated_, false);
}

AV1_INSTANTIATE_TEST_CASE(UniformTileConfigTestLarge,
                          ::testing::Values(::libaom_test::kOnePassGood,
                                            ::libaom_test::kTwoPassGood),
                          ::testing::ValuesIn(uniformTileConfigParams),
                          ::testing::Values(AOM_Q, AOM_VBR, AOM_CBR, AOM_CQ));
}  // namespace
