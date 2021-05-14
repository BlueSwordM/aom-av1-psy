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

#include <new>

#include "aom/aom_codec.h"
#include "aom/aom_external_partition.h"
#include "av1/common/blockd.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/y4m_video_source.h"
#include "test/util.h"

#if CONFIG_AV1_ENCODER
#if !CONFIG_REALTIME_ONLY
namespace {

constexpr int kFrameNum = 8;
constexpr int kVersion = 1;

typedef struct TestData {
  int version = kVersion;
} TestData;

typedef struct ToyModel {
  TestData *data;
  aom_ext_part_config_t config;
  aom_ext_part_funcs_t funcs;
} ToyModel;

aom_ext_part_status_t ext_part_create_model(
    void *priv, const aom_ext_part_config_t *part_config,
    aom_ext_part_model_t *ext_part_model) {
  TestData *received_data = reinterpret_cast<TestData *>(priv);
  EXPECT_EQ(received_data->version, kVersion);
  ToyModel *toy_model = new (std::nothrow) ToyModel;
  EXPECT_NE(toy_model, nullptr);
  toy_model->data = received_data;
  *ext_part_model = toy_model;
  EXPECT_EQ(part_config->superblock_size, BLOCK_64X64);
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_send_features(
    aom_ext_part_model_t ext_part_model,
    const aom_partition_features_t *part_features) {
  (void)ext_part_model;
  (void)part_features;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_get_partition_decision(
    aom_ext_part_model_t ext_part_model,
    aom_partition_decision_t *ext_part_decision) {
  (void)ext_part_model;
  (void)ext_part_decision;
  return AOM_EXT_PART_ERROR;
}

aom_ext_part_status_t ext_part_send_partition_stats(
    aom_ext_part_model_t ext_part_model,
    const aom_partition_stats_t *ext_part_stats) {
  (void)ext_part_model;
  (void)ext_part_stats;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_delete_model(
    aom_ext_part_model_t ext_part_model) {
  ToyModel *toy_model = static_cast<ToyModel *>(ext_part_model);
  EXPECT_EQ(toy_model->data->version, kVersion);
  delete toy_model;
  return AOM_EXT_PART_OK;
}

class ExternalPartitionTest
    : public ::libaom_test::CodecTestWith2Params<libaom_test::TestMode, int>,
      public ::libaom_test::EncoderTest {
 protected:
  ExternalPartitionTest()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        cpu_used_(GET_PARAM(2)), psnr_(0.0), nframes_(0) {}
  virtual ~ExternalPartitionTest() {}

  virtual void SetUp() {
    InitializeConfig(encoding_mode_);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = AOM_VBR;
    cfg_.g_threads = 1;
    cfg_.g_lag_in_frames = 4;
    cfg_.rc_target_bitrate = 400;
    init_flags_ = AOM_CODEC_USE_PSNR;
  }

  virtual bool DoDecode() const { return false; }

  virtual void BeginPassHook(unsigned int) {
    psnr_ = 0.0;
    nframes_ = 0;
  }

  virtual void PSNRPktHook(const aom_codec_cx_pkt_t *pkt) {
    psnr_ += pkt->data.psnr.psnr[0];
    nframes_++;
  }

  double GetAveragePsnr() const {
    if (nframes_) return psnr_ / nframes_;
    return 0.0;
  }

  void SetExternalPartition(bool use_external_partition) {
    use_external_partition_ = use_external_partition;
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      aom_ext_part_funcs_t ext_part_funcs;
      ext_part_funcs.priv = reinterpret_cast<void *>(&test_data_);
      ext_part_funcs.create_model = ext_part_create_model;
      ext_part_funcs.send_features = ext_part_send_features;
      ext_part_funcs.get_partition_decision = ext_part_get_partition_decision;
      ext_part_funcs.send_partition_stats = ext_part_send_partition_stats;
      ext_part_funcs.delete_model = ext_part_delete_model;

      encoder->Control(AOME_SET_CPUUSED, cpu_used_);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
      if (use_external_partition_) {
        encoder->Control(AV1E_SET_EXTERNAL_PARTITION, &ext_part_funcs);
      }
    }
  }

 private:
  libaom_test::TestMode encoding_mode_;
  int cpu_used_;
  double psnr_;
  unsigned int nframes_;
  bool use_external_partition_;
  TestData test_data_;
};

// Encode twice and expect the same psnr value.
// The first run is the baseline without external partition.
// The second run is to get partition decisions from the toy model we defined.
// Here, we let the partition decision return true for all stages.
// In this case, the external partition doesn't alter the original encoder
// behavior. So we expect the same encoding results.
TEST_P(ExternalPartitionTest, EncodeMatch) {
  ::libaom_test::Y4mVideoSource video("paris_352_288_30.y4m", 0, kFrameNum);
  SetExternalPartition(false);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  const double psnr = GetAveragePsnr();

  SetExternalPartition(true);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  const double psnr2 = GetAveragePsnr();

  EXPECT_DOUBLE_EQ(psnr, psnr2);
}

AV1_INSTANTIATE_TEST_SUITE(ExternalPartitionTest,
                           ::testing::Values(::libaom_test::kTwoPassGood),
                           ::testing::Values(4));  // cpu_used

}  // namespace
#endif  // !CONFIG_REALTIME_ONLY
#endif  // CONFIG_AV1_ENCODER
