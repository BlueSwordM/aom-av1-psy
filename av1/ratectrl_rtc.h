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

#ifndef AOM_AV1_RATECTRL_RTC_H_
#define AOM_AV1_RATECTRL_RTC_H_

#include <cstdint>
#include <memory>

#include "aom/aomcx.h"
#include "aom/aom_encoder.h"
#include "aom_mem/aom_mem.h"
#include "av1/encoder/encoder.h"

namespace aom {

struct AV1RateControlRtcConfig {
 public:
  AV1RateControlRtcConfig() {
    width = 1280;
    height = 720;
    max_quantizer = 63;
    min_quantizer = 2;
    target_bandwidth = 1000;
    buf_initial_sz = 600;
    buf_optimal_sz = 600;
    buf_sz = 1000;
    undershoot_pct = overshoot_pct = 50;
    max_intra_bitrate_pct = 50;
    max_inter_bitrate_pct = 0;
    framerate = 30.0;
    ts_number_layers = 1;
    aq_mode = 0;
    layer_target_bitrate[0] = static_cast<int>(target_bandwidth);
    ts_rate_decimator[0] = 1;
    av1_zero(max_quantizers);
    av1_zero(min_quantizers);
    av1_zero(scaling_factor_den);
    av1_zero(scaling_factor_num);
    av1_zero(layer_target_bitrate);
    av1_zero(ts_rate_decimator);
    scaling_factor_num[0] = 1;
    scaling_factor_den[0] = 1;
    max_quantizers[0] = max_quantizer;
    min_quantizers[0] = min_quantizer;
  }

  int width;
  int height;
  // 0-63
  int max_quantizer;
  int min_quantizer;
  int64_t target_bandwidth;
  int64_t buf_initial_sz;
  int64_t buf_optimal_sz;
  int64_t buf_sz;
  int undershoot_pct;
  int overshoot_pct;
  int max_intra_bitrate_pct;
  int max_inter_bitrate_pct;
  double framerate;
  int layer_target_bitrate[AOM_MAX_LAYERS];
  int ts_rate_decimator[AOM_MAX_TS_LAYERS];
  int aq_mode;
  // Number of spatial layers
  int ss_number_layers;
  // Number of temporal layers
  int ts_number_layers;
  int max_quantizers[AOM_MAX_LAYERS];
  int min_quantizers[AOM_MAX_LAYERS];
  int scaling_factor_num[AOM_MAX_SS_LAYERS];
  int scaling_factor_den[AOM_MAX_SS_LAYERS];
};

struct AV1FrameParamsRTC {
  FRAME_TYPE frame_type;
  int spatial_layer_id;
  int temporal_layer_id;
};

class AV1RateControlRTC {
 public:
  static std::unique_ptr<AV1RateControlRTC> Create(
      const AV1RateControlRtcConfig &cfg);
  ~AV1RateControlRTC() {
    if (cpi_) {
      if (cpi_->svc.number_spatial_layers > 1 ||
          cpi_->svc.number_temporal_layers > 1) {
        for (int sl = 0; sl < cpi_->svc.number_spatial_layers; sl++) {
          for (int tl = 0; tl < cpi_->svc.number_temporal_layers; tl++) {
            int layer =
                LAYER_IDS_TO_IDX(sl, tl, cpi_->svc.number_temporal_layers);
            LAYER_CONTEXT *const lc = &cpi_->svc.layer_context[layer];
            aom_free(lc->map);
          }
        }
      }
      if (cpi_->oxcf.q_cfg.aq_mode == CYCLIC_REFRESH_AQ) {
        aom_free(cpi_->enc_seg.map);
        cpi_->enc_seg.map = nullptr;
        av1_cyclic_refresh_free(cpi_->cyclic_refresh);
      }
      aom_free(cpi_);
    }
  }

  void UpdateRateControl(const AV1RateControlRtcConfig &rc_cfg);
  // GetQP() needs to be called after ComputeQP() to get the latest QP
  int GetQP() const;
  signed char *GetCyclicRefreshMap() const;
  int *GetDeltaQ() const;
  void ComputeQP(const AV1FrameParamsRTC &frame_params);
  // Feedback to rate control with the size of current encoded frame
  void PostEncodeUpdate(uint64_t encoded_frame_size);

 private:
  AV1RateControlRTC() = default;
  void InitRateControl(const AV1RateControlRtcConfig &cfg);
  AV1_COMP *cpi_;
};

}  // namespace aom

#endif  // AOM_AV1_RATECTRL_RTC_H_
