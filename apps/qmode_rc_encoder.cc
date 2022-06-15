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

#include <fstream>
#include <string>
#include <vector>

#include "av1/ducky_encode.h"
#include "av1/ratectrl_qmode.h"

extern "C" void usage_exit(void) { exit(EXIT_FAILURE); }

int main(int argc, const char **argv_) {
  (void)argc;
  (void)argv_;
  std::string input_file = "/export/hda3/Videos/derf/bus_cif.y4m";
  aom_rational_t frame_rate = { 30, 1 };
  aom::VideoInfo input_video = { 352, 288,       frame_rate, AOM_IMG_FMT_I420,
                                 55,  input_file };
  aom::DuckyEncode ducky_encode(input_video);
  aom::AV1RateControlQMode qmode_rc;
  aom::RateControlParam rc_param = {};
  rc_param.frame_height = 288;
  rc_param.frame_width = 352;
  rc_param.max_gop_show_frame_count = 16;
  rc_param.min_gop_show_frame_count = 4;
  rc_param.ref_frame_table_size = 7;

  const aom::Status status = qmode_rc.SetRcParam(rc_param);

  if (status.code != AOM_CODEC_OK) return -1;

  std::vector<FIRSTPASS_STATS> frame_stats =
      ducky_encode.ComputeFirstPassStats();

  fprintf(stderr, "obtained first pass\n\n");

  aom::FirstpassInfo firstpass_info = { (352 / 16 + 1) * (288 / 16 + 1),
                                        frame_stats };
  aom::GopStructList gop_list = qmode_rc.DetermineGopInfo(firstpass_info);
  ducky_encode.StartEncode(frame_stats);

  fprintf(stderr, "gop struct determined \n");

  std::vector<aom::TplGopStats> tpl_gop_stats_list =
      ducky_encode.ComputeTplStats(gop_list);

  std::vector<aom::GopEncodeInfo> gop_encode_info_list;

  // TODO(jingning): modularize AV1RateControlQMode to make the frame
  // QP decisions on a sequence basis.
  aom::RefFrameTable ref_frame_table;
  for (size_t i = 0; i < gop_list.size(); ++i) {
    aom::GopStruct &gop_struct = gop_list[i];
    aom::TplGopStats &tpl_gop_stats = tpl_gop_stats_list[i];
    (void)tpl_gop_stats;
    (void)gop_struct;

    // TODO(jingning): Extract the tpl stats through ducky_encode and make
    // frame encoding decisions.
    aom::GopEncodeInfo gop_encode_info;
    // = qmode_rc.GetGopEncodeInfo(gop_struct, tpl_gop_stats, ref_frame_table);
    ref_frame_table = gop_encode_info.final_snapshot;
    gop_encode_info_list.push_back(gop_encode_info);
  }

  ducky_encode.EndEncode();

  fprintf(stderr, "tpl stats completed.\n");

  // TODO(jingning): Re-enable the next final encoding stage once the TPL stats
  // collection is done.
  return 0;

  // Full encoding of the video sequence.
  // Do binary search with rc_param.base_q_index around this block.
  ducky_encode.StartEncode(frame_stats);
  ducky_encode.EncodeVideo(gop_list, gop_encode_info_list);
  ducky_encode.EndEncode();

  return 0;
}
