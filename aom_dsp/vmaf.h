/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_VMAF_H_
#define AOM_AOM_DSP_VMAF_H_

#include "aom_scale/yv12config.h"

#if CONFIG_USE_VMAF_RC
typedef struct VmafContext VmafContext;
typedef struct VmafModel VmafModel;
#endif

typedef struct {
  // Stores the scaling factors for rdmult when tuning for VMAF.
  // rdmult_scaling_factors[row * num_cols + col] stores the scaling factors for
  // 64x64 block at (row, col).
  double *rdmult_scaling_factors;

  // Stores the luma sse of the last frame.
  double last_frame_ysse;

  // Stores the VMAF of the last frame.
  double last_frame_vmaf;

  // Stores the filter strength of the last frame.
  double last_frame_unsharp_amount;

  // Stores the origial qindex before scaling.
  int original_qindex;

#if CONFIG_USE_VMAF_RC
  // VMAF model used in VMAF caculations.
  VmafModel *vmaf_model;
#endif
} TuneVMAFInfo;

void aom_calc_vmaf(const char *model_path, const YV12_BUFFER_CONFIG *source,
                   const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                   double *vmaf);

void aom_calc_vmaf_multi_frame(
    void *user_data, const char *model_path,
    int (*read_frame)(float *ref_data, float *main_data, float *temp_data,
                      int stride_byte, void *user_data),
    int frame_width, int frame_height, int bit_depth, double *vmaf);

#if CONFIG_USE_VMAF_RC
void aom_init_vmaf_rc(VmafModel **vmaf_model, const char *model_path);

void aom_calc_vmaf_rc(VmafModel *vmaf_model, const YV12_BUFFER_CONFIG *source,
                      const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                      int cal_vmaf_neg, double *vmaf);

void aom_close_vmaf_rc(VmafModel *vmaf_model);
#endif  // CONFIG_USE_VMAF_RC

#endif  // AOM_AOM_DSP_VMAF_H_
