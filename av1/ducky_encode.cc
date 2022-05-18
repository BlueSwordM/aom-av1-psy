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
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include "config/aom_config.h"

#include "aom/aom_encoder.h"

#include "av1/av1_cx_iface.h"
#include "av1/av1_iface_common.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/ethread.h"
#include "av1/encoder/firstpass.h"
#include "av1/ducky_encode.h"

#include "common/tools_common.h"

namespace aom {
struct EncoderResource {
  STATS_BUFFER_CTX stats_buf_ctx;
  FIRSTPASS_STATS *frame_stats_buffer;
  aom_image_t img;
  AV1_PRIMARY *ppi;
};

class DuckyEncode::EncodeImpl {
 public:
  VideoInfo video_info;
  int g_usage;
  enum aom_rc_mode rc_end_usage;
  aom_rational64_t timestamp_ratio;
  std::vector<FIRSTPASS_STATS> stats_list;
  EncoderResource enc_resource;
};

DuckyEncode::DuckyEncode(const VideoInfo &video_info) {
  impl_ptr_ = std::unique_ptr<EncodeImpl>(new EncodeImpl());
  impl_ptr_->video_info = video_info;
  impl_ptr_->g_usage = GOOD;
  impl_ptr_->rc_end_usage = AOM_Q;
  // TODO(angiebird): Set timestamp_ratio properly
  // timestamp_ratio.den = cfg->g_timebase.den;
  // timestamp_ratio.num = (int64_t)cfg->g_timebase.num * TICKS_PER_SEC;
  impl_ptr_->timestamp_ratio = { 1, 1 };
  // TODO(angiebird): How to set ptsvol and duration?
}

DuckyEncode::~DuckyEncode() {}

static AV1EncoderConfig GetEncoderConfig(const VideoInfo &video_info,
                                         int g_usage, aom_enc_pass pass) {
  const aom_codec_iface *codec = aom_codec_av1_cx();
  aom_codec_enc_cfg_t cfg;
  aom_codec_enc_config_default(codec, &cfg, g_usage);
  cfg.g_w = video_info.frame_width;
  cfg.g_h = video_info.frame_height;
  cfg.g_pass = pass;
  // g_timebase is the inverse of frame_rate
  cfg.g_timebase.num = video_info.frame_rate.den;
  cfg.g_timebase.den = video_info.frame_rate.num;
  AV1EncoderConfig oxcf = av1_get_encoder_config(&cfg);
  // TODO(angiebird): Why didn't we init use_highbitdepth in
  // av1_get_encoder_config()?
  oxcf.use_highbitdepth = 0;
  return oxcf;
}

static EncoderResource InitEncoder(const VideoInfo &video_info, int g_usage,
                                   enum aom_rc_mode rc_end_usage,
                                   aom_enc_pass pass) {
  EncoderResource enc_resource;
  aom_img_alloc(&enc_resource.img, video_info.img_fmt, video_info.frame_width,
                video_info.frame_height, /*align=*/1);
  const AV1EncoderConfig oxcf = GetEncoderConfig(video_info, g_usage, pass);
  av1_initialize_enc(g_usage, rc_end_usage);
  AV1_PRIMARY *ppi =
      av1_create_primary_compressor(nullptr,
                                    /*num_lap_buffers=*/0, &oxcf);
  enc_resource.ppi = ppi;

  assert(ppi != nullptr);
  // Turn off ppi->b_calculate_psnr to avoid calling generate_psnr_packet() in
  // av1_post_encode_updates().
  // TODO(angiebird): Modify generate_psnr_packet() to handle the case that
  // cpi->ppi->output_pkt_list = nullptr.
  ppi->b_calculate_psnr = 0;
  aom_codec_err_t res = av1_create_stats_buffer(
      &enc_resource.frame_stats_buffer, &enc_resource.stats_buf_ctx,
      /*num_lap_buffers=*/0);
  (void)res;
  assert(res == AOM_CODEC_OK);
  ppi->twopass.stats_buf_ctx = &enc_resource.stats_buf_ctx;
  BufferPool *buffer_pool = nullptr;
  res = av1_create_context_and_bufferpool(ppi, &ppi->cpi, &buffer_pool, &oxcf,
                                          ENCODE_STAGE, -1);
  // TODO(angiebird): Why didn't we set initial_dimensions in
  // av1_create_compressor()?
  ppi->cpi->initial_dimensions.width = oxcf.frm_dim_cfg.width;
  ppi->cpi->initial_dimensions.height = oxcf.frm_dim_cfg.height;
  assert(res == AOM_CODEC_OK);
  assert(ppi->cpi != nullptr);
  assert(buffer_pool != nullptr);

  ppi->cpi->twopass_frame.stats_in = ppi->twopass.stats_buf_ctx->stats_in_start;
  const AV1_COMP *cpi = ppi->cpi;
  SequenceHeader *seq_params = ppi->cpi->common.seq_params;

  ppi->seq_params_locked = 1;
  assert(ppi->lookahead == nullptr);

  int lag_in_frames = cpi->oxcf.gf_cfg.lag_in_frames;
  ppi->lookahead = av1_lookahead_init(
      cpi->oxcf.frm_dim_cfg.width, cpi->oxcf.frm_dim_cfg.height,
      seq_params->subsampling_x, seq_params->subsampling_y,
      seq_params->use_highbitdepth, lag_in_frames, cpi->oxcf.border_in_pixels,
      cpi->common.features.byte_alignment,
      /*num_lap_buffers=*/0, /*is_all_intra=*/0,
      cpi->oxcf.tool_cfg.enable_global_motion);
  assert(ppi->lookahead != nullptr);
  return enc_resource;
}

static void FreeEncoder(EncoderResource *enc_resource) {
  aom_img_free(&enc_resource->img);
  av1_destroy_stats_buffer(&enc_resource->stats_buf_ctx,
                           enc_resource->frame_stats_buffer);
  enc_resource->frame_stats_buffer = nullptr;
  enc_resource->stats_buf_ctx.total_stats = nullptr;
  enc_resource->stats_buf_ctx.total_left_stats = nullptr;
  BufferPool *buffer_pool = enc_resource->ppi->cpi->common.buffer_pool;
  av1_destroy_context_and_bufferpool(enc_resource->ppi->cpi, &buffer_pool);
  av1_remove_primary_compressor(enc_resource->ppi);
  enc_resource->ppi = nullptr;
}

std::vector<FIRSTPASS_STATS> DuckyEncode::ComputeFirstPassStats() {
  std::vector<FIRSTPASS_STATS> stats_list;
  aom_enc_pass pass = AOM_RC_FIRST_PASS;
  EncoderResource enc_resource = InitEncoder(
      impl_ptr_->video_info, impl_ptr_->g_usage, impl_ptr_->rc_end_usage, pass);
  AV1_PRIMARY *ppi = enc_resource.ppi;
  struct lookahead_ctx *lookahead = ppi->lookahead;
  int i;
  int frame_count = impl_ptr_->video_info.frame_count;
  FILE *in_file = fopen(impl_ptr_->video_info.file_path.c_str(), "r");
  aom_rational64_t timestamp_ratio = impl_ptr_->timestamp_ratio;
  // TODO(angiebird): Ideally, ComputeFirstPassStats() doesn't output
  // bitstream. Do we need bitstream buffer here?
  std::vector<uint8_t> buf(1000);
  for (i = 0; i < frame_count; ++i) {
    assert(!av1_lookahead_full(lookahead));
    if (aom_img_read(&enc_resource.img, in_file)) {
      // TODO(angiebird): Set ts_start/ts_end properly
      int64_t ts_start = i;
      int64_t ts_end = i + 1;
      YV12_BUFFER_CONFIG sd;
      image2yuvconfig(&enc_resource.img, &sd);
      av1_lookahead_push(lookahead, &sd, ts_start, ts_end,
                         /*use_highbitdepth=*/0, /*flags=*/0);
      AV1_COMP_DATA cpi_data = {};
      cpi_data.cx_data = buf.data();
      cpi_data.cx_data_sz = buf.size();
      cpi_data.frame_size = 0;
      cpi_data.flush = 1;  // Makes av1_get_compressed_data process a frame
      cpi_data.ts_frame_start = ts_start;
      cpi_data.ts_frame_end = ts_end;
      cpi_data.pop_lookahead = 1;
      cpi_data.timestamp_ratio = &timestamp_ratio;
      // av1_get_compressed_data only generates first pass stats not
      // compresses data
      int res = av1_get_compressed_data(ppi->cpi, &cpi_data);
      (void)res;
      assert(res == static_cast<int>(AOM_CODEC_OK));
      stats_list.push_back(*(ppi->twopass.stats_buf_ctx->stats_in_end - 1));
    }
  }
  av1_end_first_pass(ppi->cpi);

  FreeEncoder(&enc_resource);
  fclose(in_file);
  return stats_list;
}

void DuckyEncode::StartEncode(const std::vector<FIRSTPASS_STATS> &stats_list) {
  aom_enc_pass pass = AOM_RC_SECOND_PASS;
  impl_ptr_->enc_resource = InitEncoder(
      impl_ptr_->video_info, impl_ptr_->g_usage, impl_ptr_->rc_end_usage, pass);
  impl_ptr_->stats_list = stats_list;
}

// TODO(angiebird): Add this part in a follow-up CL
EncodeFrameResult DuckyEncode::EncodeFrame() {
  EncodeFrameResult encode_frame_result = {};
  return encode_frame_result;
}

void DuckyEncode::EndEncode() { FreeEncoder(&impl_ptr_->enc_resource); }
}  // namespace aom
