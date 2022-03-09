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

#ifndef AOM_AV1_RATECTRL_QMODE_H_
#define AOM_AV1_RATECTRL_QMODE_H_

#include <deque>
#include <queue>
#include <vector>
#include "av1/encoder/firstpass.h"
#include "av1/ratectrl_qmode_interface.h"
#include "av1/reference_manager.h"

namespace aom {

GopStruct construct_gop(RefFrameManager *ref_frame_manager,
                        int show_frame_count, int has_key_frame);

class AV1RateControlQMode : public AV1RateControlQModeInterface {
 public:
  void SetRcParam(const RateControlParam &rc_param) override;
  GopStructList DetermineGopInfo(const FirstpassInfo &firstpass_info) override;
  virtual std::vector<FrameEncodeParameters> GetGopEncodeInfo(
      const TplGopStats &tpl_stats_list) override;

 private:
  RateControlParam rc_param_;
};
}  // namespace aom

#endif
