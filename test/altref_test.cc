/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"
namespace {

class AltRefForcedKeyTestLarge
    : public ::libaom_test::CodecTestWith2Params<libaom_test::TestMode, int>,
      public ::libaom_test::EncoderTest {
 protected:
  AltRefForcedKeyTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        cpu_used_(GET_PARAM(2)), forced_kf_frame_num_(1), frame_num_(0) {}
  virtual ~AltRefForcedKeyTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    cfg_.rc_end_usage = AOM_VBR;
    cfg_.g_threads = 0;
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, cpu_used_);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
#if CONFIG_AV1_ENCODER
      // override test default for tile columns if necessary.
      if (GET_PARAM(0) == &libaom_test::kAV1) {
        encoder->Control(AV1E_SET_TILE_COLUMNS, 6);
      }
#endif
    }
    frame_flags_ =
        (video->frame() == forced_kf_frame_num_) ? AOM_EFLAG_FORCE_KF : 0;
  }

  virtual void FramePktHook(const aom_codec_cx_pkt_t *pkt) {
    if (frame_num_ == forced_kf_frame_num_) {
      ASSERT_EQ(pkt->data.frame.flags & AOM_FRAME_IS_KEY,
                static_cast<aom_codec_frame_flags_t>(AOM_FRAME_IS_KEY))
          << "Frame #" << frame_num_ << " isn't a keyframe!";
    }
    ++frame_num_;
  }

  ::libaom_test::TestMode encoding_mode_;
  int cpu_used_;
  unsigned int forced_kf_frame_num_;
  unsigned int frame_num_;
};

TEST_P(AltRefForcedKeyTestLarge, Frame1IsKey) {
  const aom_rational timebase = { 1, 30 };
  const int lag_values[] = { 3, 15, 25, -1 };

  forced_kf_frame_num_ = 1;
  for (int i = 0; lag_values[i] != -1; ++i) {
    frame_num_ = 0;
    cfg_.g_lag_in_frames = lag_values[i];
    libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       timebase.den, timebase.num, 0, 30);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }
}

TEST_P(AltRefForcedKeyTestLarge, ForcedFrameIsKey) {
  const aom_rational timebase = { 1, 30 };
  const int lag_values[] = { 3, 15, 25, -1 };

  for (int i = 0; lag_values[i] != -1; ++i) {
    frame_num_ = 0;
    forced_kf_frame_num_ = lag_values[i] - 1;
    cfg_.g_lag_in_frames = lag_values[i];
    libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       timebase.den, timebase.num, 0, 30);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }
}

AV1_INSTANTIATE_TEST_SUITE(AltRefForcedKeyTestLarge,
                           ::testing::Values(::libaom_test::kOnePassGood),
                           ::testing::Values(2, 5));

typedef struct {
  const unsigned int min_kf_dist;
  const unsigned int max_kf_dist;
  const unsigned int min_gf_interval;
  const unsigned int max_gf_interval;
  const unsigned int lag_in_frames;
  libaom_test::TestMode encoding_mode;
} AltRefTestParams;

static const AltRefTestParams TestParams[] = {
  { 0, 10, 4, 8, 10, ::libaom_test::kOnePassGood },
  { 0, 30, 8, 12, 16, ::libaom_test::kOnePassGood },
  { 30, 30, 12, 16, 25, ::libaom_test::kOnePassGood },
  { 0, 60, 12, 20, 25, ::libaom_test::kOnePassGood },
  { 60, 60, 16, 28, 30, ::libaom_test::kOnePassGood },
  { 0, 100, 16, 32, 35, ::libaom_test::kOnePassGood },
  { 0, 10, 4, 8, 10, ::libaom_test::kTwoPassGood },
  { 0, 30, 8, 12, 16, ::libaom_test::kTwoPassGood },
  { 30, 30, 12, 16, 25, ::libaom_test::kTwoPassGood },
  { 0, 60, 16, 24, 25, ::libaom_test::kTwoPassGood },
  { 60, 60, 20, 28, 30, ::libaom_test::kTwoPassGood },
  { 0, 100, 24, 32, 35, ::libaom_test::kTwoPassGood },
};

std::ostream &operator<<(std::ostream &os, const AltRefTestParams &test_arg) {
  return os << "AltRefTestParams { min_kf_dist:" << test_arg.min_kf_dist
            << " max_kf_dist:" << test_arg.max_kf_dist
            << " min_gf_interval:" << test_arg.min_gf_interval
            << " max_gf_interval:" << test_arg.max_gf_interval
            << " lag_in_frames:" << test_arg.lag_in_frames
            << " encoding_mode:" << test_arg.encoding_mode << " }";
}

// This class is used to check the presence of altref frame.
class AltRefFramePresenceTestLarge
    : public ::libaom_test::CodecTestWith2Params<AltRefTestParams, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  AltRefFramePresenceTestLarge()
      : EncoderTest(GET_PARAM(0)), altref_test_params_(GET_PARAM(1)),
        rc_end_usage_(GET_PARAM(2)) {
    is_arf_frame_present_ = 0;
  }
  virtual ~AltRefFramePresenceTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(altref_test_params_.encoding_mode);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = rc_end_usage_;
    cfg_.g_threads = 1;
    cfg_.kf_min_dist = altref_test_params_.min_kf_dist;
    cfg_.kf_max_dist = altref_test_params_.max_kf_dist;
    cfg_.g_lag_in_frames = altref_test_params_.lag_in_frames;
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, 5);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
      encoder->Control(AV1E_SET_MIN_GF_INTERVAL,
                       altref_test_params_.min_gf_interval);
      encoder->Control(AV1E_SET_MAX_GF_INTERVAL,
                       altref_test_params_.max_gf_interval);
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (is_arf_frame_present_ != 1 && AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_ALTREF_PRESENT,
                                    &is_arf_frame_present_);
    }
    return AOM_CODEC_OK == res_dec;
  }

  const AltRefTestParams altref_test_params_;
  int is_arf_frame_present_;
  aom_rc_mode rc_end_usage_;
};

TEST_P(AltRefFramePresenceTestLarge, AltRefFrameEncodePresenceTest) {
  libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                     cfg_.g_timebase.den, cfg_.g_timebase.num,
                                     0, 100);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(is_arf_frame_present_, 1);
}

AV1_INSTANTIATE_TEST_SUITE(AltRefFramePresenceTestLarge,
                           ::testing::ValuesIn(TestParams),
                           ::testing::Values(AOM_Q, AOM_VBR, AOM_CBR, AOM_CQ));

}  // namespace
